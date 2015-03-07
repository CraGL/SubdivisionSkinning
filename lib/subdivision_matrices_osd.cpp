// c++ -g -std=c++11 -I../../opensubdiv -I../../regression subdivision_matrices_osd.cpp -o subdivision_matrices_osd -DSUBDIVISION_MATRICES_OSD_MAIN ../../build/lib/libosdCPU.a

#include "subdivision_matrices_osd.h"

#include <vector>
#include <cassert>
#include <memory> // unique_ptr
#include <algorithm> // fill

#include <osd/cpuVertexBuffer.h>
#include <osd/cpuComputeContext.h>
#include <osd/cpuComputeController.h>
#include <far/mesh.h>

namespace
{

/*
Given the number of vertices in the mesh 'num_vertices'
and a vector of faces 'faces', each element of which is a vector of vertex indices,
returns a new mesh object suitable for passing to 
Free the resulting mesh with "delete" when finished with it.
*/
OpenSubdiv::HbrMesh<OpenSubdiv::OsdVertex>* newMeshFromNumVerticesAndFaces( int num_vertices, const subdivision_matrix::faces_t& faces )
{
    assert( num_vertices > 0 );
    
    // Following: http://graphics.pixar.com/opensubdiv/forum.html?place=msg%2Fopensubdiv%2FzKq9vG_azHQ%2FDb2K7L23BykJ
    
    static OpenSubdiv::HbrCatmarkSubdivision<OpenSubdiv::OsdVertex> CATMULL_CLARK;
    
    auto mesh = new OpenSubdiv::HbrMesh<OpenSubdiv::OsdVertex>( &CATMULL_CLARK );
    
    OpenSubdiv::OsdVertex default_vertex_value;
    for( int vi = 0; vi < num_vertices; ++vi ) mesh->NewVertex( vi, default_vertex_value );
    
    /// TODO Q: Is this necessary?
    // int ptex_index = 0;
    int face_offset = 0;
    for( const auto& num_face_vertices : faces.first )
    {
        OpenSubdiv::HbrFace<OpenSubdiv::OsdVertex>* face = mesh->NewFace( num_face_vertices, &faces.second[face_offset], 0 );
        
        /*
        face->SetPtexIndex( ptex_index );
        if( f.size() != 4 ) {
            ptex_index += f.size();
        } else {
            ptex_index += 1;
        }
        */
        
        face_offset += num_face_vertices;
    }
    
    mesh->SetInterpolateBoundaryMethod( OpenSubdiv::HbrMesh<OpenSubdiv::OsdVertex>::k_InterpolateBoundaryAlwaysSharp );
    
    mesh->Finish();
    
    return mesh;
}

}

namespace subdivision_matrix
{

/*
Given the number of vertices in the mesh 'num_vertices'
and a vector of faces 'faces', each element of which is a vector of vertex indices, and
a positive integer 'level' indicating the level of refinement,
fills 'positions_out' with sparse vectors such that the position
for the i-th uv value in the original call to precomputeStencils() can be obtained by:
    \sum_j control_points[ positions_out[i][j].first ] * positions_out[i][j].second

The optional parameter 'quad_faces_out', if specified, is a sequence of quad faces
obtained by subdivision, where each face is four indices into 'positions_out'.
*/
void compute_subdivision_coefficients_for_mesh(
    int num_vertices,
    const faces_t& faces,
    const int level,
    std::vector< sparse_vector_t >& positions_out,
    std::vector< index_t >* quad_faces_out
    )
{
    auto mesh = newMeshFromNumVerticesAndFaces( num_vertices, faces );
    compute_subdivision_coefficients_for_mesh( num_vertices, mesh, level, positions_out, quad_faces_out );
    delete mesh;
}

void compute_subdivision_coefficients_for_mesh(
    int num_coarse_vertices,
    OpenSubdiv::HbrMesh<OpenSubdiv::OsdVertex>* hmesh,
    const int level,
    std::vector< sparse_vector_t >& positions_out,
    std::vector< index_t >* quad_faces_out
    )
{
    assert( num_coarse_vertices == hmesh->GetNumVertices() );
    
    assert( level > 0 );
    positions_out.clear();
    if( quad_faces_out ) quad_faces_out->clear();
    
    std::unique_ptr< OpenSubdiv::FarMesh<OpenSubdiv::OsdVertex> > fmesh;
    // Adaptive would give us limit vertex positions, however:
    //   Adaptive is only supported for Catmull-Clark subdivision.
    //   Adaptive doesn't support getting quads.
    // So, we use non-adaptive aka uniform subdivision.
    {
        // The comment in "osdutil/refiner.cpp" says:
        // XXX:gelder
        // Had problems with patchArrayVector in patchTables not working
        // unless firstLevel is passed as 1 here.  Bug in refiner?
        OpenSubdiv::FarMeshFactory<OpenSubdiv::OsdVertex> factory( hmesh, level, false /* adaptive */, 1 /* firstLevel */ );
        fmesh.reset( factory.Create() );
    }
    
    
    // Following "osdutil/refiner.cpp"
    
    /// Refiner::Initialize()
    // Subdivision tables describe the addition steps with coefficients
    // needed to perform subdivision
    const OpenSubdiv::FarSubdivisionTables* ftable = fmesh->GetSubdivisionTables();
    
    
    const int kNumDimensions = std::min( 100, num_coarse_vertices );
    {
        /// UniformEvaluator::Initialize()
        std::unique_ptr< OpenSubdiv::OsdCpuComputeContext > _computeContext( OpenSubdiv::OsdCpuComputeContext::Create(fmesh->GetSubdivisionTables(), fmesh->GetVertexEditTables()) );
        // The vertex buffer must have the FarMesh's vertex count, which includes
        // the refined vertices, not the HbrMesh's vertex count, which only has
        // the coarse vertices.
        std::unique_ptr< OpenSubdiv::OsdCpuVertexBuffer > _vertexBuffer( OpenSubdiv::OsdCpuVertexBuffer::Create( kNumDimensions, fmesh->GetNumVertices() ) );
        
        for( int coarse_index = 0; coarse_index < num_coarse_vertices; coarse_index += kNumDimensions )
        {
            /// UniformEvaluator::SetCoarsePositions()
            // Each vertex stores its indicator function.
            // The identity matrix is the indicator function.
            // We are iterating vertex-by-vertex, so set the coarse_index column of the
            // identity matrix.
            {
                // vbuffer: v0d0 v0d1 v0d2 ... v0dN v1d0 v1d1 v1d2 ... v1dN v2d0 v2d1 ...
                float* vbuffer = _vertexBuffer->BindCpuBuffer();
                // TODO Q: Do I need to fill it with zeros, or can I just
                //         unset the last iteration's 1?
                //         In other words, does Refine() overwrite these values?
                std::fill( vbuffer, vbuffer + kNumDimensions * num_coarse_vertices, 0. );
                // A slice of 'kNumDimensions' columns of the identity matrix:
                for( int i = coarse_index; i < std::min( coarse_index + kNumDimensions, num_coarse_vertices ); ++i )
                {
                    vbuffer[ i-coarse_index + i*kNumDimensions ] = 1.;
                }
            }
            
            /// UniformEvaluator::Refine()
            {
                OpenSubdiv::OsdCpuComputeController cpuComputeController;
                cpuComputeController.Refine(
                    _computeContext.get(),
                    fmesh->GetKernelBatches(),
                    _vertexBuffer.get()
                    );
                // TODO Q: Is this Synchronize() call necessary?
                cpuComputeController.Synchronize();
            }
            
            
            /// UniformEvaluator::GetRefinedPositions()
            {
                const int numRefinedVerts = ftable->GetNumVertices(level);
                const int firstVertexOffset = ftable->GetFirstVertexOffset(level);
                const float* vbuffer = _vertexBuffer->BindCpuBuffer();
                
                if (numRefinedVerts == 0) {
                    std::cerr << "ERROR: No refinement occurred.\n";
                    return;
                }
                
                if( 0 == coarse_index )
                {
                    // We clear it upon entry to this function.
                    assert( positions_out.empty() );
                    positions_out.resize( numRefinedVerts );
                    for( int vi = 0; vi < numRefinedVerts; ++vi ) positions_out.at( vi ).reserve( 16 );
                }
                else
                {
                    assert( positions_out.size() == numRefinedVerts );
                }
                
                const float* off = vbuffer + firstVertexOffset*kNumDimensions;
                const float kEps = 1e-10;
                for( int vi = 0; vi < numRefinedVerts; ++vi )
                {
                    for( int i = coarse_index; i < std::min( coarse_index + kNumDimensions, num_coarse_vertices ); ++i )
                    {
                        const auto& val = off[ i - coarse_index ];
                        if( fabs( val ) > kEps )
                        {
                            positions_out.at( vi ).push_back( std::make_pair( i, val ) );
                        }
                    }
                    off += kNumDimensions;
                }
            }
        }
    }
    
    // Faces
    if( quad_faces_out )
    {
        // We clear it upon entry to this function.
        assert( quad_faces_out->empty() );
        
        // From "osdutil/refiner.cpp"'s Refiner::Initialize():
        // Find quads array at given level
        const OpenSubdiv::FarPatchTables* ptables = fmesh->GetPatchTables();
        const OpenSubdiv::FarPatchTables::PatchArrayVector& parrays = ptables->GetPatchArrayVector();
        // parrays doesn't contain base mesh, so it starts with level==1
        assert( level-1 < parrays.size() );
        const OpenSubdiv::FarPatchTables::PatchArray& parray = parrays[level-1];
        
        // Global index of the first point in this array
        const int firstVertexOffset = ftable->GetFirstVertexOffset(level);
        
        
        /// From "osdutil/refiner.cpp"'s Refiner::GetRefinedQuads():
        const int numUniformQuads = (int) parray.GetNumPatches();
        quad_faces_out->resize(numUniformQuads * 4);
        
        const unsigned int *quadIndices = ptables->GetFaceVertices( level );
        for( int i = 0; i < numUniformQuads*4; ++i )
        {
            quad_faces_out->at(i) = quadIndices[i] - firstVertexOffset;
        }
    }
}

} // ~namespace subdivision_matrix()

#ifdef SUBDIVISION_MATRICES_OSD_MAIN
void print_sparse_vectors( const std::vector< subdivision_matrix::sparse_vector_t >& sparse_vectors, const std::string& label );
void print_faces( const std::vector< subdivision_matrix::index_t >& faces, const std::string& label );
void print_OBJ( const std::vector<float>& original_positions, const std::vector< subdivision_matrix::sparse_vector_t >& sparse_vectors, const std::vector< subdivision_matrix::index_t >& faces, const std::string& label );

#define TORUS 32, { {4, 5, 1, 0}, {5, 6, 2, 1}, {6, 7, 3, 2}, {7, 4, 0, 3}, {8, 9, 5, 4}, {9, 10, 6, 5}, {10, 11, 7, 6}, {11, 8, 4, 7}, {12, 13, 9, 8}, {13, 14, 10, 9}, {14, 15, 11, 10}, {15, 12, 8, 11}, {16, 17, 13, 12}, {17, 18, 14, 13}, {18, 19, 15, 14}, {19, 16, 12, 15}, {20, 21, 17, 16}, {21, 22, 18, 17}, {22, 23, 19, 18}, {23, 20, 16, 19}, {24, 25, 21, 20}, {25, 26, 22, 21}, {26, 27, 23, 22}, {27, 24, 20, 23}, {28, 29, 25, 24}, {29, 30, 26, 25}, {30, 31, 27, 26}, {31, 28, 24, 27}, {0, 1, 29, 28}, {1, 2, 30, 29}, {2, 3, 31, 30}, {3, 0, 28, 31} }
#define CUBE 8, { {0, 1, 3, 2}, {2, 3, 5, 4}, {4, 5, 7, 6}, {6, 7, 1, 0}, {1, 7, 5, 3}, {6, 0, 2, 4} }

void test_compute_coefficients( int level, std::vector< subdivision_matrix::sparse_vector_t >& position_coeffs )
{
    using namespace subdivision_matrix;
    
    // std::vector< sparse_vector_t > position_coeffs;
    std::vector< index_t > faces;
    compute_subdivision_coefficients_for_mesh(
        // Cube
        // CUBE,
        // Torus
        TORUS,
        level,
        position_coeffs, &faces
        );
    
    print_sparse_vectors( position_coeffs, "position_coeffs: " );
    print_faces( faces, "faces: " );
}

// -I/path/to/regression
#include <common/shape_utils.h>
#include <shapes/catmark_torus.h>
#include <shapes/catmark_cube.h>
void test_shape_utils( int level, std::vector< subdivision_matrix::sparse_vector_t >& position_coeffs )
{
    using namespace subdivision_matrix;
    
    // Create a torus.
    std::vector<float> orgPositions;
    auto mesh = simpleHbr<OpenSubdiv::OsdVertex>( catmark_torus.c_str(), kCatmark, orgPositions, true);
    // auto mesh = simpleHbr<OpenSubdiv::OsdVertex>( catmark_cube.c_str(), kCatmark, orgPositions, true);
    
    // std::vector< sparse_vector_t > position_coeffs;
    std::vector< index_t > faces;
    compute_subdivision_coefficients_for_mesh(
        mesh->GetNumVertices(),
        mesh,
        level,
        position_coeffs, &faces
        );
    
    print_sparse_vectors( position_coeffs, "position_coeffs: " );
    print_faces( faces, "faces: " );
    print_OBJ( orgPositions, position_coeffs, faces, "reconstructed OBJ" );
}

#include <iostream>
#include <cstdlib>
void usage( const char* argv0 )
{
    std::cerr << "Usage: " << argv0 << " test_shape_utils|test_compute_coefficients level\n";
    exit( -1 );
}
int main( int argc, char* argv[] )
{
    using std::string;
    
    if( argc != 3 ) usage( argv[0] );
    
    const int level = atoi( argv[2] );
    if( level <= 0 ) usage( argv[0] );
    
    const string test_name = string(argv[1]);
    
    std::vector< subdivision_matrix::sparse_vector_t > position_coeffs;
    if( string("test_shape_utils") == test_name )
    {
        test_shape_utils( level, position_coeffs );
    }
    else if( string("test_compute_coefficients") == test_name )
    {
        test_compute_coefficients( level, position_coeffs );
    }
    else
    {
        usage( argv[0] );
    }
    
    // TODO: Do something with position_coeffs
    
    return 0;
}

void print_sparse_vectors( const std::vector< subdivision_matrix::sparse_vector_t >& sparse_vectors, const std::string& label = "" )
{
    std::cout << label << "{";
    for( const auto vec : sparse_vectors )
    {
        std::cout << "\n    { ";
        for( const auto p : vec )
        {
            std::cout << "{ " << p.first << ", " << p.second << " }, ";
        }
        std::cout << "},";
    }
    std::cout << "}\n";
}
void print_faces( const std::vector< subdivision_matrix::index_t >& faces, const std::string& label )
{
    std::cout << label << "{\n";
    for( int fi = 0; fi < faces.size(); ++fi )
    {
        if( fi > 0 && fi % 4 == 0 ) std::cout << " },\n";
        if( fi % 4 == 0 ) std::cout << "    {";
        std::cout << ' ' << faces.at(fi);
    }
    std::cout << "}\n";
}
#include <valarray>
void print_OBJ( const std::vector<float>& original_positions, const std::vector< subdivision_matrix::sparse_vector_t >& sparse_vectors, const std::vector< subdivision_matrix::index_t >& faces, const std::string& label )
{
    std::cout << "# OBJ: " << label << '\n';
    
    assert( original_positions.size() % 3 == 0 );
    std::vector< std::valarray< float > > ovs( original_positions.size() / 3, { 0., 0., 0. } );
    for( int i = 0; i < ovs.size(); ++i )
    {
        assert( ovs.at( i ).size() == 3 );
        ovs.at( i )[ 0 ] = original_positions.at( 3*i + 0 );
        ovs.at( i )[ 1 ] = original_positions.at( 3*i + 1 );
        ovs.at( i )[ 2 ] = original_positions.at( 3*i + 2 );
    }
    
    for( int i = 0; i < sparse_vectors.size(); ++i )
    {
        std::valarray< float > pos( { 0., 0., 0. } );
        for( const auto& i_w : sparse_vectors.at( i ) )
        {
            pos += i_w.second * ovs.at( i_w.first );
        }
        std::cout << "v " << pos[0] << ' ' << pos[1] << ' ' << pos[2] << '\n';
    }
    
    assert( faces.size() % 4 == 0 );
    for( int i = 0; i+3 < faces.size(); i += 4 )
    {
        std::cout << "f " << (faces[i+0]+1) << ' ' << (faces[i+1]+1) << ' ' << (faces[i+2]+1) << ' ' << (faces[i+3]+1) << '\n';
    }
}

#endif

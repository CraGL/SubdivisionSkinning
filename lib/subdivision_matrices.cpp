// c++ -g -std=c++11 -I/Users/yotam/Work/ext/OpenSubdiv/opensubdiv -I/Users/yotam/Work/ext/OpenSubdiv/regression subdivision_matrices.cpp -o subdivision_matrices -DSUBDIVISION_MATRICES_MAIN

#include "subdivision_matrices.h"

#include <vector>
#include <cassert>

namespace
{

/*
Given the number of vertices in the mesh 'num_vertices'
and a vector of faces 'faces', each element of which is a vector of vertex indices,
returns a new mesh object suitable for passing to 
Free the resulting mesh with "delete" when finished with it.
*/
OpenSubdiv::HbrMesh<OpenSubdiv::FarStencilFactoryVertex>* newMeshFromNumVerticesAndFaces( int num_vertices, const subdivision_matrix::faces_t& faces )
{
    assert( num_vertices > 0 );
    
    // Following: http://graphics.pixar.com/opensubdiv/forum.html?place=msg%2Fopensubdiv%2FzKq9vG_azHQ%2FDb2K7L23BykJ
    
    static OpenSubdiv::HbrCatmarkSubdivision<OpenSubdiv::FarStencilFactoryVertex> CATMULL_CLARK;
    
    auto mesh = new OpenSubdiv::HbrMesh<OpenSubdiv::FarStencilFactoryVertex>( &CATMULL_CLARK );
    
    OpenSubdiv::FarStencilFactoryVertex default_vertex_value;
    for( int vi = 0; vi < num_vertices; ++vi ) mesh->NewVertex( vi, default_vertex_value );
    
    /// TODO Q: Is this necessary?
    // int ptex_index = 0;
    int face_offset = 0;
    for( const auto& num_face_vertices : faces.first )
    {
        OpenSubdiv::HbrFace<OpenSubdiv::FarStencilFactoryVertex>* face = mesh->NewFace( num_face_vertices, &faces.second[face_offset], 0 );
        
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
    
    mesh->SetInterpolateBoundaryMethod( OpenSubdiv::HbrMesh<OpenSubdiv::FarStencilFactoryVertex>::k_InterpolateBoundaryAlwaysSharp );
    
    mesh->Finish();
    
    return mesh;
}

}

namespace subdivision_matrix
{

/*
Given a desired number of "u" and "v" parameters sampling the range [0,1]x[0,1],
returns in the output vectors 'us' and 'vs' coordinates sampling the range
in equal intervals.
*/
void createUVs( int NUM_U, int NUM_V, std::vector< real_t >& us, std::vector< real_t >& vs )
{
    assert( NUM_U > 0 );
    assert( NUM_V > 0 );
    
    us.clear();
    vs.clear();
    
    us.reserve( NUM_U * NUM_V );
    vs.reserve( NUM_U * NUM_V );
    
    /// XXX NOTE: Keep this in sync with subdivision_limit_mesh::createUVsWithZeroAndOne()
    for( int ui = 0; ui < NUM_U; ++ui )
    {
        // Sample from .5 to .95, so that each sample has an equal piece of area.
        const real_t u = real_t( ui + 0.5 )/NUM_U;
        // const real_t u = real_t( ui )/(NUM_U-1);
        for( int vi = 0; vi < NUM_V; ++vi )
        {
            // Sample from .5 to .95, so that each sample has an equal piece of area.
            const real_t v = real_t( vi + 0.5 )/NUM_V;
            // const real_t v = real_t( vi )/(NUM_V-1);
            
            us.push_back( u );
            vs.push_back( v );
        }
    }
}

/*
Given a mesh,
two same-length vectors 'us' and 'vs' where the ( us[i], vs[i] ) are the uv locations,
and an optional integer parameter 'reflevel' ("In most cases only approximate evaluation
is done by linearly interpolating between limit values after refLevel subdivision steps."),
returns a FarStencilTables object suitable for calling UpdateValues() and UpdateDerivs().
*/
OpenSubdiv::FarStencilTables precomputeStencils(
    OpenSubdiv::HbrMesh<OpenSubdiv::FarStencilFactoryVertex>* mesh,
    const std::vector< float >& us,
    const std::vector< float >& vs,
    int reflevel = 5
    )
{
    assert( mesh );
    assert( us.size() == vs.size() );
    assert( !us.empty() );
    assert( reflevel > 0 );
    
    OpenSubdiv::FarStencilTables controlStencils;
    
    OpenSubdiv::FarStencilTablesFactory<> factory(mesh);
    
    // std::cout << "Stencil num faces: " << mesh->GetNumFaces() << '\n';
    
    // Check if we are using catmull-clark subdivision.
    // const bool catmark = 0 != dynamic_cast< OpenSubdiv::HbrCatmarkSubdivision<OpenSubdiv::FarStencilFactoryVertex>* >( mesh->GetSubdivision() );
    
    // XXXX UDPATE: mesh->GetNumFaces() changes during calls to AppendStencils()!
    //              We end up with lots of duplicate stencils if we don't store the
    //              initial value of GetNumFaces(). Better yet, store GetNumCoarseFaces().
    //              I checked the HbrMesh code and the initial non-NULL faces are all
    //              coarse. I'm not sure how a NULL face could get into the structure,
    //              so let's assert that faces aren't NULL.
    const int nfaces = mesh->GetNumCoarseFaces();
    for( int i = 0; i < nfaces; ++i ) {
    
        const OpenSubdiv::HbrFace<OpenSubdiv::FarStencilFactoryVertex>* f = mesh->GetFace(i);
        assert( f );
        
        const int nv = f->GetNumVertices();
		
        // TODO Q: Why doesn't this die with loop subdivision?
        // if( catmark && nv!=4 ) {
        if( nv!=4 ) {
    
            // if the face is not a quad, we have to iterate over sub-quad(rants)
            for (int j=0; j<f->GetNumVertices(); ++j) {
                
                // std::cout << "Non-quad face with " << nv << " vertices.\n";
    
                factory.SetCurrentFace(i,j);
                
                const int appended = factory.AppendStencils( &controlStencils, us.size(), &us[0], &vs[0], reflevel );
                // std::cout << "Appended " << appended << " stencils.\n";
            }
        } else {
    
            factory.SetCurrentFace(i);
            
            const int appended = factory.AppendStencils( &controlStencils, us.size(), &us[0], &vs[0], reflevel );
            // std::cout << "Appended " << appended << " stencils.\n";
        }
    }
    	
    return controlStencils;
}

/*
Given OpenSubdiv::FarStencilTables as returned by precomputeStencils()
and the number of control points 'num_control_points'
fills 'positions_out' with sparse vectors such that the position
for the i-th uv value in the original call to precomputeStencils() can be obtained by:
    \sum_j control_points[ positions_out[i][j].first ] * positions_out[i][j].second

The optional parameters 'du_out' and 'dv_out', if specified, are similar to 'positions_out'
except that the above summation results in du and dv.
*/
void sparseVectorsForPositionsAndDerivativesFromFarStencilTables(
    const OpenSubdiv::FarStencilTables& controlStencils,
    int num_control_points,
    std::vector< sparse_vector_t >& positions_out,
    std::vector< sparse_vector_t >* du_out = nullptr,
    std::vector< sparse_vector_t >* dv_out = nullptr
    )
{
    /// 1 Create a vector of control points, where each control point is a sparse vector containing only its integer index.
    /// 2 Create a vector of output coefficients, whose length is the total number
    ///   of all positions on the mesh we are evaluating.
    /// 3 Create adapters for passing to OpenSubdiv::FarStencilTables.UpdateValues().
    /// 4 Call OpenSubdiv::FarStencilTables.UpdateValues().
    /// 5 Repeat with OpenSubdiv::FarStencilTables.UpdateDerivs().
    
    assert( num_control_points > 0 );
    // The derivative outputs, du_out and dv_out,
    // must be either both null or both given.
    assert( ( du_out == nullptr ) == ( dv_out == nullptr ) );
    // du_out and dv_out must be distinct unless they are both null.
    assert( du_out != dv_out || du_out == nullptr );
    
    /// 1
    // controlPoints are just indices!
    std::vector< sparse_vector_t > controlPoints( num_control_points );
    static const real_t kMagicValue{ 31337 };
    for( int i = 0; i < num_control_points; ++i ) controlPoints[i].push_back( std::make_pair( i, kMagicValue ) );
    
    
    /// 2
    positions_out.clear();
    positions_out.resize( controlStencils.GetNumStencils() );
    
    
    /// 3
    struct SparseVectorAdapter
    {
        sparse_vector_t* data;
        
        void Clear() { assert( data ); data->clear(); }
        void AddWithWeight( const SparseVectorAdapter& other, real_t value )
        {
            assert( data );
            assert( other.data );
            assert( other.data->size() == 1 );
            assert( (*other.data)[0].second == kMagicValue );
            
            const index_t index = (*other.data)[0].first;
            data->push_back( std::make_pair( index, value ) );
        }
    };
    
    std::vector< SparseVectorAdapter > controlPoints_adapters( controlPoints.size() );
    for( int i = 0; i < controlPoints_adapters.size(); ++i ) controlPoints_adapters[i].data = &controlPoints[i];
    
    std::vector< SparseVectorAdapter > positions_adapters( positions_out.size() );
    for( int i = 0; i < positions_adapters.size(); ++i ) positions_adapters[i].data = &positions_out[i];
    
    
    /// 4
    controlStencils.UpdateValues< SparseVectorAdapter >( &controlPoints_adapters[0], &positions_adapters[0] );
    
    
    if( du_out && dv_out )
    {
        du_out->clear();
        dv_out->clear();
        du_out->resize( controlStencils.GetNumStencils() );
        dv_out->resize( controlStencils.GetNumStencils() );
        
        std::vector< SparseVectorAdapter > du_adapters( du_out->size() );
        std::vector< SparseVectorAdapter > dv_adapters( dv_out->size() );
        for( int i = 0; i < du_adapters.size(); ++i ) du_adapters[i].data = &(*du_out)[i];
        for( int i = 0; i < du_adapters.size(); ++i ) dv_adapters[i].data = &(*dv_out)[i];
        
        controlStencils.UpdateDerivs< SparseVectorAdapter >( &controlPoints_adapters[0], &du_adapters[0], &dv_adapters[0] );
    }
}

/*
Given the number of vertices in the mesh 'num_vertices'
and a vector of faces 'faces', each element of which is a vector of vertex indices, and
two same-length vectors 'us' and 'vs' where the ( us[i], vs[i] ) are the uv locations,
fills 'positions_out' with sparse vectors such that the position
for the i-th uv value in the original call to precomputeStencils() can be obtained by:
    \sum_j control_points[ positions_out[i][j].first ] * positions_out[i][j].second

The optional parameters 'du_out' and 'dv_out', if specified, are similar to 'positions_out'
except that the above summation results in du and dv.
*/
void compute_subdivision_coefficients_for_mesh(
    int num_vertices,
    const faces_t& faces,
    const std::vector< real_t >& us,
    const std::vector< real_t >& vs,
    std::vector< sparse_vector_t >& positions_out,
    std::vector< sparse_vector_t >* du_out,
    std::vector< sparse_vector_t >* dv_out
    )
{
	auto mesh = newMeshFromNumVerticesAndFaces( num_vertices, faces );
    compute_subdivision_coefficients_for_mesh( num_vertices, mesh, us, vs, positions_out, du_out, dv_out );
    delete mesh;
}

void compute_subdivision_coefficients_for_mesh(
    int num_vertices,
    OpenSubdiv::HbrMesh<OpenSubdiv::FarStencilFactoryVertex>* mesh,
    const std::vector< real_t >& us,
    const std::vector< real_t >& vs,
    std::vector< sparse_vector_t >& positions_out,
    std::vector< sparse_vector_t >* du_out,
    std::vector< sparse_vector_t >* dv_out
    )
{
    std::vector< float > real_us, real_vs;
	for ( auto u : us )
    	real_us.push_back( u );
    for ( auto v : vs )
    	real_vs.push_back( v );
    
    auto precomputed = precomputeStencils( mesh, real_us, real_vs );
    sparseVectorsForPositionsAndDerivativesFromFarStencilTables( precomputed, num_vertices, positions_out, du_out, dv_out );
}

} // ~namespace subdivision_matrix()

#ifdef SUBDIVISION_MATRICES_MAIN
void print_sparse_vectors( const std::vector< subdivision_matrix::sparse_vector_t >& sparse_vectors, const std::string& label );

#define TORUS 32, { {4, 5, 1, 0}, {5, 6, 2, 1}, {6, 7, 3, 2}, {7, 4, 0, 3}, {8, 9, 5, 4}, {9, 10, 6, 5}, {10, 11, 7, 6}, {11, 8, 4, 7}, {12, 13, 9, 8}, {13, 14, 10, 9}, {14, 15, 11, 10}, {15, 12, 8, 11}, {16, 17, 13, 12}, {17, 18, 14, 13}, {18, 19, 15, 14}, {19, 16, 12, 15}, {20, 21, 17, 16}, {21, 22, 18, 17}, {22, 23, 19, 18}, {23, 20, 16, 19}, {24, 25, 21, 20}, {25, 26, 22, 21}, {26, 27, 23, 22}, {27, 24, 20, 23}, {28, 29, 25, 24}, {29, 30, 26, 25}, {30, 31, 27, 26}, {31, 28, 24, 27}, {0, 1, 29, 28}, {1, 2, 30, 29}, {2, 3, 31, 30}, {3, 0, 28, 31} }
#define CUBE 8, { {0, 1, 3, 2}, {2, 3, 5, 4}, {4, 5, 7, 6}, {6, 7, 1, 0}, {1, 7, 5, 3}, {6, 0, 2, 4} }

void test_newMesh_sparseVectors( int resolution, std::vector< subdivision_matrix::sparse_vector_t >& position_coeffs )
{
    using namespace subdivision_matrix;
    
    auto mesh = newMeshFromNumVerticesAndFaces(
        // Cube
        // CUBE
        // Torus
        TORUS
        );
    
    std::vector< float > us, vs;
    createUVs( resolution, resolution, us, vs );
    OpenSubdiv::FarStencilTables controlStencils = precomputeStencils( mesh, us, vs );
    
    // std::vector< sparse_vector_t > position_coeffs;
    std::vector< sparse_vector_t > du_coeffs;
    std::vector< sparse_vector_t > dv_coeffs;
    sparseVectorsForPositionsAndDerivativesFromFarStencilTables( controlStencils, mesh->GetNumVertices(), position_coeffs, &du_coeffs, &dv_coeffs );
    
    print_sparse_vectors( position_coeffs, "position_coeffs: " );
    print_sparse_vectors( du_coeffs, "du_coeffs: " );
    print_sparse_vectors( dv_coeffs, "dv_coeffs: " );
    
    delete mesh;
}

void test_compute_coefficients( int resolution, std::vector< subdivision_matrix::sparse_vector_t >& position_coeffs )
{
    using namespace subdivision_matrix;
    
    std::vector< real_t > us, vs;
    createUVs( resolution, resolution, us, vs );
    
    // std::vector< sparse_vector_t > position_coeffs;
    std::vector< sparse_vector_t > du_coeffs;
    std::vector< sparse_vector_t > dv_coeffs;
    compute_subdivision_coefficients_for_mesh(
        // Cube
        // CUBE,
        // Torus
        TORUS,
        us, vs,
        position_coeffs, &du_coeffs, &dv_coeffs
        );
    
    print_sparse_vectors( position_coeffs, "position_coeffs: " );
    print_sparse_vectors( du_coeffs, "du_coeffs: " );
    print_sparse_vectors( dv_coeffs, "dv_coeffs: " );
}

// -I/path/to/regression
#include <common/shape_utils.h>
#include <shapes/catmark_torus.h>
#include <shapes/catmark_cube.h>
void test_shape_utils( int resolution, std::vector< subdivision_matrix::sparse_vector_t >& position_coeffs )
{
    using namespace subdivision_matrix;
    
    // Create a torus.
    std::vector<float> orgPositions;
    auto mesh = simpleHbr<OpenSubdiv::FarStencilFactoryVertex>( catmark_torus.c_str(), kCatmark, orgPositions, true);
    // auto mesh = simpleHbr<OpenSubdiv::FarStencilFactoryVertex>( catmark_cube.c_str(), kCatmark, orgPositions, true);
    
    std::vector< float > us, vs;
    createUVs( resolution, resolution, us, vs );
    OpenSubdiv::FarStencilTables controlStencils = precomputeStencils( mesh, us, vs );
    

//     auto p_weights = controlStencils.GetWeights();
//     auto du_weights = controlStencils.GetDuWeights();
//     auto dv_weights = controlStencils.GetDvWeights();

    
    // std::vector< sparse_vector_t > position_coeffs;
    std::vector< sparse_vector_t > du_coeffs;
    std::vector< sparse_vector_t > dv_coeffs;
    sparseVectorsForPositionsAndDerivativesFromFarStencilTables( controlStencils, mesh->GetNumVertices(), position_coeffs, &du_coeffs, &dv_coeffs );
    
    print_sparse_vectors( position_coeffs, "position_coeffs: " );
    print_sparse_vectors( du_coeffs, "du_coeffs: " );
    print_sparse_vectors( dv_coeffs, "dv_coeffs: " );
}

#include <iostream>
#include <cstdlib>
void usage( const char* argv0 )
{
    std::cerr << "Usage: " << argv0 << " test_shape_utils|test_compute_coefficients|test_newMesh_sparseVectors resolution\n";
    exit( -1 );
}
int main( int argc, char* argv[] )
{
    using std::string;
    
    if( argc != 3 ) usage( argv[0] );
    
    const int resolution = atoi( argv[2] );
    if( resolution < 0 ) usage( argv[0] );
    
    const string test_name = string(argv[1]);
    
    std::vector< subdivision_matrix::sparse_vector_t > position_coeffs;
    if( string("test_shape_utils") == test_name )
    {
        test_shape_utils( resolution, position_coeffs );
    }
    else if( string("test_compute_coefficients") == test_name )
    {
        test_compute_coefficients( resolution, position_coeffs );
    }
    else if( string("test_newMesh_sparseVectors") == test_name )
    {
        test_newMesh_sparseVectors( resolution, position_coeffs );
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
        std::cout << " },";
    }
    std::cout << "}\n";
}

#endif

// c++ -g -std=c++11 -I/Users/yotam/Work/ext/OpenSubdiv/opensubdiv -I/Users/yotam/Work/ext/OpenSubdiv/regression subdivision_limit_mesh.cpp subdivision_matrices.cpp -o subdivision_limit_mesh -DSUBDIVISION_LIMIT_MESH_MAIN

#include "subdivision_limit_mesh.h"
#include "save_obj.h"

#include <vector>
#include <cassert>

namespace subdivision_limit_mesh
{

void deduplicate_vertices_and_indices( subdivision_matrix::MatrixX3_t duplicated_vertices, std::vector< std::vector< real_t > >& unique_vertices, std::vector< std::vector< int > >& limit_quad_faces_out );
void efficient_deduplicate_vertices_and_indices( int num_u, int num_v, subdivision_matrix::MatrixX3_t duplicated_vertices, std::vector< std::vector< real_t > >& unique_vertices, std::vector< std::vector< int > >& limit_quad_faces_out ) ;

// The same as subdivision_matrices::createUVs() but samples 0 and 1.
// This is a helper function for compute_quad_limit_mesh_for_mesh().
void createUVsWithZeroAndOne( int NUM_U, int NUM_V, std::vector< real_t >& us, std::vector< real_t >& vs )
{
    assert( NUM_U > 0 );
    assert( NUM_V > 0 );
    
    us.clear();
    vs.clear();
    
    us.reserve( NUM_U * NUM_V );
    vs.reserve( NUM_U * NUM_V );
    
    /// 1 Create the u and v samples, which are identical to the ones in createUVs()
    ///   but have 0 and 1 at the beginning and end.
    std::vector< real_t > u_samples, v_samples;
    u_samples.reserve( NUM_U + 2 );
    u_samples.push_back( 0. );
    for( int ui = 0; ui < NUM_U; ++ui )
    {
        // Sample from .5 to .95, so that each sample has an equal piece of area.
        const real_t u = real_t( ui + 0.5 )/NUM_U;
        u_samples.push_back( u );
    }
    u_samples.push_back( 1.0 );
    assert( u_samples.size() == NUM_U + 2 );
    
    v_samples.reserve( NUM_U + 2 );
    v_samples.push_back( 0. );
    for( int vi = 0; vi < NUM_V; ++vi )
    {
        // Sample from .5 to .95, so that each sample has an equal piece of area.
        const real_t v = real_t( vi + 0.5 )/NUM_V;
        v_samples.push_back( v );
    }
    v_samples.push_back( 1.0 );
    assert( v_samples.size() == NUM_V + 2 );
    
    /// 2 Put the Cartesian product in us and vs.
    for( int ui = 0; ui < u_samples.size(); ++ui )
    {
        const real_t u = u_samples.at(ui);
        for( int vi = 0; vi < v_samples.size(); ++vi )
        {
            const real_t v = v_samples.at(vi);
            
            us.push_back( u );
            vs.push_back( v );
        }
    }
}

void compute_quad_limit_mesh_for_mesh(
    const subdivision_matrix::MatrixX3_t& vertices,
    OpenSubdiv::HbrMesh<OpenSubdiv::FarStencilFactoryVertex>* mesh,
    int num_u, int num_v,
    std::vector< std::vector< real_t > >& limit_vertices_out,
    std::vector< std::vector< int > >& limit_quad_faces_out
    )
{
    /// 1 createUVs just like subdivision_matrices::createUVs() except with 0 and 1 in each row/column.
    /// 2 Extract the M matrix for computing limit vertex positions.
    /// 3 Calculate the positions.
    /// 4 Create a quad mesh for every "num_u+2 x num_v+2" vertices (they sample a control mesh face).
    /// 5 De-duplicate the vertices.
    
    /// 1
    std::vector< real_t > us, vs;
    createUVsWithZeroAndOne( num_u, num_v, us, vs );
    
    
    /// 2
    subdivision_matrix::SparseMatrix_t positions;
    subdivision_matrix::compute_subdivision_coefficients_for_mesh(
        vertices.rows(), mesh,
        us, vs,
        positions
        );
    
    
    /// 3
    subdivision_matrix::MatrixX3_t intermediate_vertices = positions * vertices;
//     std::cout << intermediate_vertices.rows() << " " << intermediate_vertices.size() << std::endl;
//     limit_vertices_out = positions * vertices;
    
    
    /// 4
    // Because of how createUVsWithZeroAndOne() works, we know that
    // every "num_u+2 x num_v+2" vertices sample a control mesh quad face in a grid.
    // Create faces for the sub-quads that make up this grid.
    // The order of elements in 'limit_vertices_out' is, for each control mesh quad face:
    //     (0,0), (0,.5/num_v), ..., (0,1), (.5/num_u,0), (.5/num_u,.5/num_v), ..., (.5/num_u,1), ...
    // In other words, every 'num_v+2' vertices forms a column (strip of constant u value).
    // The number of vertices should be an integer multiple of the number of control mesh quad faces
    // (n.b. the number of control mesh quad faces may be increased by extraordinary vertex sub-faces).
    assert( intermediate_vertices.rows() % ( (num_u+2)*(num_v+2) ) == 0 );
    for( int off = 0; off < intermediate_vertices.rows(); off += (num_u+2)*(num_v+2) )
    {
        // Iterate over each quad and add it to the faces.
        for( int u_off = 1; u_off < num_u + 2; ++u_off )
        {
			for( int v_off = 1; v_off < num_v + 2; ++v_off )
			{
				// TODO: Add the face:
				// (off + (u_off-1,v_off-1)), (off + (u_off-1,v_off)), (off + (u_off,v_off)), (off + (u_off,v_off-1)), clockwise
				// To convert from (off + (ui,vi)) to an actual vertex index,
				// it is something like: off + ui*(num_u+2) + vi.
				std::vector< int > quad{(off + (u_off-1)*(num_u+2) + v_off-1), 
										(off + (u_off-1)*(num_u+2) + v_off), 
										(off + u_off*(num_u+2)	   + v_off), 
										(off + u_off*(num_u+2)	   + v_off-1)};
				limit_quad_faces_out.push_back( quad );
			}
       		
        }
    }
    
    /// 5
    // TODO: De-duplicate vertices with the same position.
// 	deduplicate_vertices_and_indices( intermediate_vertices, limit_vertices_out, limit_quad_faces_out );
	efficient_deduplicate_vertices_and_indices( num_u+2, num_v+2, intermediate_vertices, limit_vertices_out, limit_quad_faces_out );
	/// 6 Save the vertices and faces to an obj file.
	save_obj::save_mesh( "model.obj", limit_quad_faces_out, limit_vertices_out );
	
	for( auto quad: limit_quad_faces_out ) 
    {
    	for( auto index : quad )
    		std::cout << index << " ";
    	std::cout << std::endl;
    }
    /* test
    for( int i=0; i<intermediate_vertices.rows(); i++ ) {
    	std::vector< real_t > v;
    	for( int j=0; j<intermediate_vertices.cols(); j++ ) {
    		v.push_back( intermediate_vertices(i,j) );
    	}
    	limit_vertices_out.push_back( v );
    } 

    for( auto quad: limit_quad_faces_out ) 
    {
    	for( auto index : quad )
    		std::cout << index << " ";
    	std::cout << std::endl;
    }
    for( auto vertex : limit_vertices_out  ) 
    {
    	for( auto coord : vertex )
    		std::cout << coord << " ";
    	std::cout << std::endl;
    }
    std::cout << limit_vertices_out.size() << " " << old_num << std::endl;
 	std::cout << "test: " << old_num << std::endl;
	for( int i=0; i<old_num; i++)
		std::cout << old2new[i] << std::endl;
	*/
}

void deduplicate_vertices_and_indices( subdivision_matrix::MatrixX3_t duplicated_vertices, std::vector< std::vector< real_t > >& unique_vertices, std::vector< std::vector< int > >& limit_quad_faces_out ) 
{
	int old_num = duplicated_vertices.rows();
    int old2new [ old_num ];
    std::vector< real_t > vertex { 0, 0, 0 };

    for( int i=0; i<old_num; ++i )
    {
		vertex[0] = duplicated_vertices(i,0);
		vertex[1] = duplicated_vertices(i,1);
		vertex[2] = duplicated_vertices(i,2);
		
		int p = 0;
		for( const auto& new_vertex : unique_vertices ) {
			if ( new_vertex[0] == vertex[0] and new_vertex[1] == vertex[1] and new_vertex[2] == vertex[2] )
				break;
			p++;
		}

        if( p != unique_vertices.size() ) 
            old2new[ i ] = p;
        else {
            unique_vertices.push_back( vertex );
            old2new[ i ] = unique_vertices.size()-1;
        }
    } 

    // now re-write the vertex indices in faces.
    
    for( int i=0; i<limit_quad_faces_out.size(); i++ ) 
    {
    	for( int j=0; j<limit_quad_faces_out[i].size(); j++ )
    		limit_quad_faces_out[i][j] = old2new[ limit_quad_faces_out[i][j] ];
    }
}

void efficient_deduplicate_vertices_and_indices( int num_u, int num_v, subdivision_matrix::MatrixX3_t duplicated_vertices, std::vector< std::vector< real_t > >& unique_vertices, std::vector< std::vector< int > >& limit_quad_faces_out ) 
{
	int old_num = duplicated_vertices.rows();
    int old2new [ old_num ];
    struct lookup {
    	int index;
    	std::vector< real_t > pos;
    };
    std::vector< lookup > lookup_list;
    
    for( int i=0; i<old_num; ++i )
    {
    	std::vector< real_t > vertex { duplicated_vertices(i,0), duplicated_vertices(i,1), duplicated_vertices(i,2) };
		
		int res = i % ( num_u*num_v );
		if( res < num_v or res >= (num_u-1)*num_v or res % num_v == 0 or ( (res+1) % num_v ) == 0 ) {			
			int p = 0;
			for( auto element : lookup_list ) {
				std::vector< real_t > pos = element.pos;
				if ( pos[0] == vertex[0] and pos[1] == vertex[1] and pos[2] == vertex[2] )
					break;
				p++;
			}
			if( p != lookup_list.size() ) 
            	old2new[ i ] = lookup_list[p].index;
            else {
            	unique_vertices.push_back( vertex );
            	old2new[ i ] = unique_vertices.size()-1;
            	
            	lookup new_lookup;
            	new_lookup.index = old2new[ i ];
            	new_lookup.pos = vertex;
            	lookup_list.push_back( new_lookup );
            }
            
		}
		else {
			unique_vertices.push_back( vertex );
			old2new[ i ] = unique_vertices.size()-1;
		}
    } 

    // now re-write the vertex indices in faces.
    
    for( int i=0; i<limit_quad_faces_out.size(); i++ ) 
    {
    	for( int j=0; j<limit_quad_faces_out[i].size(); j++ )
    		limit_quad_faces_out[i][j] = old2new[ limit_quad_faces_out[i][j] ];
    }
}

} // ~namespace subdivision_limit_mesh()

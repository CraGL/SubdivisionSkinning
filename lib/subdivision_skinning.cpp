#include "subdivision_skinning.h"
#include "subdivision_engine.h"
#include "subdivision_limit_mesh.h"

#include <cassert>
#include <iostream>
#include <cstdlib>
#include "ColorMap.h"

#define TORUS_VS    1.25052, 0.517982, 0.353553, 0.597239, 0.247384, 0.353553, 0.597239, 0.247384, -0.353553, 1.25052, 0.517982, -0.353553, 0.517982, 1.25052, 0.353553, 0.247384, 0.597239, 0.353553, 0.247384, 0.597239, -0.353553, 0.517982, 1.25052, -0.353553, -0.517982, 1.25052, 0.353553, -0.247384, 0.597239, 0.353553, -0.247384, 0.597239, -0.353553, -0.517982, 1.25052, -0.353553, -1.25052, 0.517982, 0.353553, -0.597239, 0.247384, 0.353553, -0.597239, 0.247384, -0.353553, -1.25052, 0.517982, -0.353553, -1.25052, -0.517982, 0.353553, -0.597239, -0.247384, 0.353553, -0.597239, -0.247384, -0.353553, -1.25052, -0.517982, -0.353553, -0.517982, -1.25052, 0.353553, -0.247384, -0.597239, 0.353553, -0.247384, -0.597239, -0.353553, -0.517982, -1.25052, -0.353553, 0.517982, -1.25052, 0.353553, 0.247384, -0.597239, 0.353553, 0.247384, -0.597239, -0.353553, 0.517982, -1.25052, -0.353553, 1.25052, -0.517982, 0.353553, 0.597239, -0.247384, 0.353553, 0.597239, -0.247384, -0.353553, 1.25052, -0.517982, -0.353553
#define TORUS_FACES { {4, 5, 1, 0}, {5, 6, 2, 1}, {6, 7, 3, 2}, {7, 4, 0, 3}, {8, 9, 5, 4}, {9, 10, 6, 5}, {10, 11, 7, 6}, {11, 8, 4, 7}, {12, 13, 9, 8}, {13, 14, 10, 9}, {14, 15, 11, 10}, {15, 12, 8, 11}, {16, 17, 13, 12}, {17, 18, 14, 13}, {18, 19, 15, 14}, {19, 16, 12, 15}, {20, 21, 17, 16}, {21, 22, 18, 17}, {22, 23, 19, 18}, {23, 20, 16, 19}, {24, 25, 21, 20}, {25, 26, 22, 21}, {26, 27, 23, 22}, {27, 24, 20, 23}, {28, 29, 25, 24}, {29, 30, 26, 25}, {30, 31, 27, 26}, {31, 28, 24, 27}, {0, 1, 29, 28}, {1, 2, 30, 29}, {2, 3, 31, 30}, {3, 0, 28, 31} }

namespace
{
    const int kResolution = 10;
    const int kLevel = 3;
}

namespace subdivision_skinning
{
    struct engine_t {
    	int num_of_controls;
        subdivision_matrix::subdivision_control_mesh mesh;
//         subdivision_matrix::MatrixX3_t handle_positions;
        subdivision_matrix::SparseMatrix_t M_matrices;
        std::vector< subdivision_matrix::MatrixX4_t > W_matrices;
       	std::vector< subdivision_matrix::MatrixX4_t > naive_W_matrices;
       	std::vector< subdivision_matrix::MatrixX4_t > target_W_matrices;
    };
    
    // delete the result when you are finished with it.
    engine_t* new_engine_for_control_mesh( const vertices_t& vertices, const faces_t& faces, const handles_t& handle_positions, const char* weight_function, OpenSubdiv::HbrMesh<OpenSubdiv::FarStencilFactoryVertex>* hmesh )
    {
    	Tick *tick1 = new Tick( "preparation" );
    	assert( vertices.size() % 3 == 0 );
    	int size = vertices.size() / 3;
    	
    	// build subdivision_matrix::mesh
		subdivision_matrix::subdivision_control_mesh mesh;
        mesh.faces = faces;
        mesh.vs.resize( size, 3 );    
        for( int i = 0; i < size; i++ )
        	mesh.vs.row( i ) << vertices[i*3], vertices[i*3+1], vertices[i*3+2];
		std::cout << "vertices #: " << size << " faces #: " << faces.first.size() << std::endl;
		 	
    	// build engine.
    	engine_t* engine = new engine_t;
        engine->num_of_controls = size;
        engine->mesh = mesh;
        
    	std::vector< subdivision_matrix::real_t > us, vs;
    	int resolution = kResolution;
		subdivision_matrix::createUVs( resolution, resolution, us, vs );
	
		subdivision_matrix::SparseMatrix_t M_matrices, Du_matrices, Dv_matrices;
		if( hmesh ) {
		    subdivision_matrix::compute_subdivision_coefficients_for_mesh(
		        mesh.vs.rows(),
                hmesh,
                us, vs,
                M_matrices, &Du_matrices, &Dv_matrices );
		} else {
            subdivision_matrix::compute_subdivision_coefficients_for_mesh(
                mesh.vs.rows(),
                mesh.faces, 
                us, vs, 
                M_matrices, &Du_matrices, &Dv_matrices );
		}
		
		engine->M_matrices = M_matrices;
		
		// std::vector< std::vector< subdivision_matrix::real_t > > limit_vertices_out;
		// std::vector< std::vector< int > > limit_quad_faces_out;
		// subdivision_limit_mesh::compute_quad_limit_mesh_for_mesh( mesh.vs, hmesh, resolution, resolution, limit_vertices_out, limit_quad_faces_out );
        
        // Convert the handle positions' type.
        assert( handle_positions.cols() == 3 );
       	subdivision_matrix::MatrixX3_t subdivision_handles( handle_positions.rows(), 3  );
       	
       	for( int i = 0; i < handle_positions.rows(); i++ )
        	subdivision_handles.row( i ) << handle_positions(i,0), 
        									handle_positions(i,1), 
        									handle_positions(i,2);
//         engine->handle_positions = subdivision_handles;
        // prepare the precompute matrices and save them in engine.
        delete tick1;
        
        Tick *tick2 = new Tick( "Pre-compute naive controls." );
        subdivision_matrix::naive_prepare( mesh.vs, subdivision_handles, engine->naive_W_matrices );
		delete tick2;
		
		Tick *tick3 = new Tick( "Pre-compute our controls." );
		subdivision_matrix::prepare( mesh.vs, M_matrices, Du_matrices, Dv_matrices, 
										subdivision_handles, weight_function, 
										engine->W_matrices );
		delete tick3;
		
		Tick *tick4 = new Tick( "Pre-compute target surface." );								
		subdivision_matrix::naive_prepare( M_matrices*mesh.vs, subdivision_handles, engine->target_W_matrices );
		delete tick4;
		
			 
        return engine;
    }
    // delete the result when you are finished with it.
    engine_t* new_engine_for_control_mesh( const vertices_t& vertices, const faces_t& faces, const handles_t& handle_positions, const char* weight_function, OpenSubdiv::HbrMesh<OpenSubdiv::OsdVertex>* hmesh, std::vector< subdivision_matrix::index_t >* quad_faces_out )
    {
    	Tick *tick1 = new Tick( "preparation" );
    	assert( vertices.size() % 3 == 0 );
    	int size = vertices.size() / 3;
    	
    	// build subdivision_matrix::mesh
		subdivision_matrix::subdivision_control_mesh mesh;
        mesh.faces = faces;
        mesh.vs.resize( size, 3 );    
        for( int i = 0; i < size; i++ )
        	mesh.vs.row( i ) << vertices[i*3], vertices[i*3+1], vertices[i*3+2];
		std::cout << "vertices #: " << size << " faces #: " << faces.first.size() << std::endl;
		 	
    	// build engine.
    	engine_t* engine = new engine_t;
        engine->num_of_controls = size;
        engine->mesh = mesh;
        
		subdivision_matrix::SparseMatrix_t M_matrices;
		if( hmesh ) {
		    subdivision_matrix::compute_subdivision_coefficients_for_mesh(
		        mesh.vs.rows(),
                hmesh,
                kLevel,
                M_matrices,
                quad_faces_out
                );
		} else {
            subdivision_matrix::compute_subdivision_coefficients_for_mesh(
                mesh.vs.rows(),
                mesh.faces, 
                kLevel,
                M_matrices,
                quad_faces_out
                );
		}
		
		engine->M_matrices = M_matrices;
		
		// std::vector< std::vector< subdivision_matrix::real_t > > limit_vertices_out;
		// std::vector< std::vector< int > > limit_quad_faces_out;
		// subdivision_limit_mesh::compute_quad_limit_mesh_for_mesh( mesh.vs, hmesh, resolution, resolution, limit_vertices_out, limit_quad_faces_out );
        
        // Convert the handle positions' type.
        assert( handle_positions.cols() == 3 );
       	subdivision_matrix::MatrixX3_t subdivision_handles( handle_positions.rows(), 3  );
       	
       	for( int i = 0; i < handle_positions.rows(); i++ )
        	subdivision_handles.row( i ) << handle_positions(i,0), 
        									handle_positions(i,1), 
        									handle_positions(i,2);
//         engine->handle_positions = subdivision_handles;
        // prepare the precompute matrices and save them in engine.
        delete tick1;
        
        Tick *tick2 = new Tick( "Pre-compute naive controls." );
        subdivision_matrix::naive_prepare( mesh.vs, subdivision_handles, engine->naive_W_matrices );
		delete tick2;
		
		Tick *tick3 = new Tick( "Pre-compute our controls." );
		subdivision_matrix::prepare( mesh.vs, M_matrices,
										subdivision_handles, weight_function, 
										engine->W_matrices );
		delete tick3;
		
		Tick *tick4 = new Tick( "Pre-compute target surface." );								
		subdivision_matrix::naive_prepare( M_matrices*mesh.vs, subdivision_handles, engine->target_W_matrices );
		delete tick4;
		
			 
        return engine;
    }
    // delete the result when you are finished with it.
    engine_t* new_engine_for_control_mesh( const vertices_t& vertices, const faces_t& faces, const subdivision_matrix::SparseMatrix_t& M_matrices, const subdivision_matrix::MatrixXX_t& weights )
    {
    	Tick* tick1 = new Tick( "preparation" );
        
        assert( vertices.size() % 3 == 0 );
        const int size = vertices.size() / 3;
        
        // build engine.
        engine_t* engine = new engine_t;
        engine->num_of_controls = size;
        
        // build subdivision_matrix::mesh
        subdivision_matrix::subdivision_control_mesh& mesh = engine->mesh;
        mesh.faces = faces;
        mesh.vs.resize( size, 3 );
        for( int i = 0; i < size; i++ )
        {
            mesh.vs.row( i ) << vertices[i*3], vertices[i*3+1], vertices[i*3+2];
        }
        std::cout << "vertices #: " << size << " faces #: " << faces.first.size() << std::endl;
        
        engine->M_matrices = M_matrices;
		delete tick1;
        
        
        Tick* tick3 = new Tick( "Pre-compute our controls." );
        subdivision_matrix::prepare( mesh.vs, M_matrices, weights, engine->W_matrices );
		delete tick3;
		
		return engine;
    }
    
    void set_new_handle_positions( engine_t* engine, const handles_t& handle_positions, const char* weight_function )
    {
        assert( engine );
        
        // Convert the handle positions' type.
        assert( handle_positions.cols() == 3 );
       	subdivision_matrix::MatrixX3_t subdivision_handles = handle_positions;
       	// engine->handle_positions = subdivision_handles;
        
        // prepare the precompute matrices and save them in engine.
        Tick *tick2 = new Tick( "Pre-compute naive controls." );
        subdivision_matrix::naive_prepare( engine->mesh.vs, subdivision_handles, engine->naive_W_matrices );
		delete tick2;
		
		Tick *tick3 = new Tick( "Pre-compute our controls." );
		subdivision_matrix::prepare( engine->mesh.vs, engine->M_matrices,
										subdivision_handles, weight_function, 
										engine->W_matrices );
		delete tick3;
		
		Tick *tick4 = new Tick( "Pre-compute target surface." );								
		subdivision_matrix::naive_prepare( engine->M_matrices*engine->mesh.vs, subdivision_handles, engine->target_W_matrices );
		delete tick4;
    }
    
    void update_control_vertices_given_handle_transforms( const engine_t* engine, const std::vector< transform_t >& transforms, vertices_t& vertices_out )
    {
    	assert( engine );
        assert( transforms.size() > 0 );
        assert( transforms.size() == engine->W_matrices.size() );
        
//         Tick *tick = new Tick( "update." );
        // Convert the transforms' type.
		std::vector< subdivision_matrix::transform_t > subdivision_transforms;
		for ( auto transform : transforms ) {
			assert( transform.cols() == 4 && transform.rows() == 4 );
			
			subdivision_matrix::transform_t converted_transform;
			for ( int i = 0; i < 4; i++ )
				converted_transform.row(i) << transform(i,0), transform(i,1), transform(i,2), transform(i,3);
				
			subdivision_transforms.push_back( converted_transform );		
		}
        
        const int nverts = engine->num_of_controls;
        
        subdivision_matrix::MatrixX3_t deformed_controls( nverts, 3 );
        subdivision_matrix::solve( engine->W_matrices, subdivision_transforms, deformed_controls );
        
        vertices_out.resize( nverts * 3 );
        for( int vi = 0; vi < nverts; ++vi )
        {
            vertices_out[ 3*vi + 0 ] = deformed_controls(vi, 0);
            vertices_out[ 3*vi + 1 ] = deformed_controls(vi, 1);
            vertices_out[ 3*vi + 2 ] = deformed_controls(vi, 2);
        }
//         delete tick;
        /*
        subdivision_matrix::MatrixX3_t b_min(1, 3), b_max(1, 3);
        b_min = engine->mesh.vs.colwise().minCoeff();
    	b_max = engine->mesh.vs.colwise().maxCoeff();	
    	real_t unit = real_t( (b_max - b_min).norm() );
    	
    	subdivision_matrix::MatrixX3_t target = engine->M_matrices * engine->mesh.vs;
      	subdivision_matrix::real_t our_hausdorff;
      	our_hausdorff = subdivision_matrix::hausdorff_distance( target, engine->M_matrices*deformed_controls );
      	std::cout << "our hausdorff percentage:\n" << our_hausdorff / unit <<"\n\n";
      	*/
      	
    }
    
    
    void update_control_vertices_given_handle_transforms( const engine_t* engine, const std::vector< transform_t >& transforms, vertices_t& vertices_out, const char* approach, vertices_t& g_targetPositions, vertices_t& g_targetColors )
    {
        assert( engine );
        assert( transforms.size() > 0 );
        assert( transforms.size() == engine->W_matrices.size() );
        
        // Convert the transforms' type.
		std::vector< subdivision_matrix::transform_t > subdivision_transforms;
		for ( auto transform : transforms ) {
			assert( transform.cols() == 4 && transform.rows() == 4 );
			
			subdivision_matrix::transform_t converted_transform;
			for ( int i = 0; i < 4; i++ )
				converted_transform.row(i) << transform(i,0), transform(i,1), transform(i,2), transform(i,3);
				
			subdivision_transforms.push_back( converted_transform );		
		}
        
        const int nverts = engine->num_of_controls;
        
        // compute the deformed control points' positions.
        subdivision_matrix::MatrixX3_t deformed_controls( nverts, 3 ), naive_result( nverts, 3 ), result( nverts, 3 );

        Tick *tick1 = new Tick( "update naive controls." );
		subdivision_matrix::solve( engine->naive_W_matrices, subdivision_transforms, naive_result );	
		delete tick1;
		
		
		Tick *tick2 = new Tick( "update our controls." );
		subdivision_matrix::solve( engine->W_matrices, subdivision_transforms, result );
        delete tick2;
                
        if ( strcmp( approach, "naive") == 0 )
        	deformed_controls = naive_result;
        else
        	deformed_controls = result;
	
		Tick *tick3 = new Tick( "update target surface." );
		subdivision_matrix::MatrixX3_t target_positions;
		subdivision_matrix::solve( engine->target_W_matrices, subdivision_transforms, target_positions );
		delete tick3;
		
/*		
		subdivision_matrix::MatrixX3_t our_positions = engine->M_matrices * result;
		subdivision_matrix::MatrixX3_t naive_positions = engine->M_matrices * naive_result;		
        subdivision_matrix::real_t our_diff 
        	= subdivision_matrix::norm_of_piecewise_distance( target_positions, our_positions );
        subdivision_matrix::real_t naive_diff 
        	= subdivision_matrix::norm_of_piecewise_distance( target_positions, naive_positions );
        std::cout << "naive vs ours:\n" << naive_diff << ' ' << our_diff << "\n\n";
*/       
		
		Tick *tick4 = new Tick( "compute piecewise distance." );
		subdivision_matrix::VectorX_t our_diffs 
        	= subdivision_matrix::piecewise_distance( target_positions, engine->M_matrices*result );
        subdivision_matrix::VectorX_t naive_diffs 
        	= subdivision_matrix::piecewise_distance( target_positions, engine->M_matrices*naive_result );
	 	
        subdivision_matrix::VectorX_t diff_values;
        if ( strcmp( approach, "naive") == 0 )
        	diff_values = naive_diffs;
        else
        	diff_values = our_diffs;
        delete tick4;
                	
#ifdef HAUSDORFF_COMPARISON      	
      	subdivision_matrix::MatrixX3_t target = engine->M_matrices * engine->mesh.vs;
      	subdivision_matrix::real_t naive_hausdorff, our_hausdorff;
      	naive_hausdorff = subdivision_matrix::hausdorff_distance( target, engine->M_matrices*naive_result );
      	our_hausdorff = subdivision_matrix::hausdorff_distance( target, engine->M_matrices*result );
      	std::cout << "naive hausdorff vs our hausdorff:\n" << naive_hausdorff << " " << our_hausdorff <<"\n\n";
#endif       	
        // Convert the new control points' positions from matrix to an one-row array.
    	vertices_out.resize( nverts * 3 );
    	g_targetPositions.resize( target_positions.size() );
    	g_targetColors.resize( target_positions.size() );
    	
        for( int vi = 0; vi < nverts; ++vi )
        {
            vertices_out[ 3*vi + 0 ] = deformed_controls(vi, 0);
            vertices_out[ 3*vi + 1 ] = deformed_controls(vi, 1);
            vertices_out[ 3*vi + 2 ] = deformed_controls(vi, 2);
        }
        
    	subdivision_matrix::MatrixX3_t b_min(1, 3), b_max(1, 3); 
    	b_min = engine->mesh.vs.colwise().minCoeff();
    	b_max = engine->mesh.vs.colwise().maxCoeff();
    	
    	real_t unit = real_t( (b_max - b_min).norm() );
        
        std::cout << unit << std::endl;
        ColorMap< real_t > colorMap(0, 0.1*unit);
        for( int vi = 0; vi < target_positions.rows(); ++vi ) 
        {
        	g_targetPositions[ 3*vi + 0 ] = target_positions(vi, 0);
        	g_targetPositions[ 3*vi + 1 ] = target_positions(vi, 1);
        	g_targetPositions[ 3*vi + 2 ] = target_positions(vi, 2);
        	colorMap.Value2Color( diff_values( vi ), &g_targetColors[ 3*vi ] );

        }
        // XXX END Debugging: Apply the first transform to all vertices.
    }
    void delete_engine_t( engine_t* engine )
    {
        delete engine;
    }
    
    void limit_vertices_for_initial_sample_locations( const engine_t* engine, vertices_t& limit_vertices_out )
    {
        assert( engine );
        
        subdivision_matrix::MatrixX3_t limit_vertices = engine->M_matrices * engine->mesh.vs;
        
        limit_vertices_out.clear();
        limit_vertices_out.resize( limit_vertices.rows()*3, -31337 );
        for( int i = 0; i < limit_vertices.rows(); ++i )
        {
            limit_vertices_out.at( 3*i + 0 ) = limit_vertices( i, 0 );
            limit_vertices_out.at( 3*i + 1 ) = limit_vertices( i, 1 );
            limit_vertices_out.at( 3*i + 2 ) = limit_vertices( i, 2 );
        }
    }
}

#ifdef SUBDIVISION_SKINNING_MAIN
void test_our_approach()
{
	float vs[] = { TORUS_VS };
	std::vector<float> vertices (vs, vs + sizeof(vs) / sizeof(float) ); 
	subdivision_skinning::faces_t faces = TORUS_FACES;
	
	subdivision_skinning::handles_t handle_positions;
	handle_positions.resize(2,3);
	handle_positions << 0.8, 0.8, 0.7,
						0.2, -0.3, -0.2;
	
	const char* weight_function = "shepard";
	
	subdivision_skinning::engine_t *engine = subdivision_skinning::new_engine_for_control_mesh( vertices, faces, handle_positions, weight_function);
	
	std::vector< subdivision_skinning::transform_t > transforms;
	subdivision_skinning::transform_t t1, t2;
   	t1.setIdentity(); 
   	t2.setIdentity();
   	transforms.push_back( t1 );
   	transforms.push_back( t2 );
	
	subdivision_skinning::vertices_t vertices_out;
	
	update_control_vertices_given_handle_transforms( engine, transforms, vertices_out );
	
	for ( auto v : vertices_out )
		std::cout << v << ' ';
	std::cout << std::endl;
}


int main( int argc, char* argv[] )
{
	if ( argc < 2 || argv[1].strcmp( "ours" ) == 0 )
		test_our_approach();
	else if ( argv[1].strcmp( "naive" ) == 0 )
		test_naive_approach();
		
	return 0;
}
#endif

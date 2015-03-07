/*
clang++ -std=c++11 -stdlib=libc++ -g -I../../opensubdiv -I../../regression subdivision_matrices_eigen.cpp subdivision_matrices.cpp subdivision_engine.cpp -o subdivision -DSUBDIVISION_ENGINE_MAIN -I/opt/local/include/eigen3
*/

#include "subdivision_engine.h"
#include <iostream>
#include <fstream>
#include <cstdlib>

const float EPS = 1e-7;
const int resolution = 10;
#define CUBE_VS    0.000000, -1.414214, 1.000000, 1.414214, 0.000000, 1.000000, -1.414214, 0.000000, 1.000000, 0.000000, 1.414214, 1.000000, -1.414214, 0.000000, -1.000000, 0.000000, 1.414214, -1.000000, 0.000000, -1.414214, -1.000000, 1.414214, 0.000000, -1.000000
#define CUBE_FACES { {0, 1, 3, 2}, {2, 3, 5, 4}, {4, 5, 7, 6}, {6, 7, 1, 0}, {1, 7, 5, 3}, {6, 0, 2, 4} }

#define TORUS_VS    1.25052, 0.517982, 0.353553, 0.597239, 0.247384, 0.353553, 0.597239, 0.247384, -0.353553, 1.25052, 0.517982, -0.353553, 0.517982, 1.25052, 0.353553, 0.247384, 0.597239, 0.353553, 0.247384, 0.597239, -0.353553, 0.517982, 1.25052, -0.353553, -0.517982, 1.25052, 0.353553, -0.247384, 0.597239, 0.353553, -0.247384, 0.597239, -0.353553, -0.517982, 1.25052, -0.353553, -1.25052, 0.517982, 0.353553, -0.597239, 0.247384, 0.353553, -0.597239, 0.247384, -0.353553, -1.25052, 0.517982, -0.353553, -1.25052, -0.517982, 0.353553, -0.597239, -0.247384, 0.353553, -0.597239, -0.247384, -0.353553, -1.25052, -0.517982, -0.353553, -0.517982, -1.25052, 0.353553, -0.247384, -0.597239, 0.353553, -0.247384, -0.597239, -0.353553, -0.517982, -1.25052, -0.353553, 0.517982, -1.25052, 0.353553, 0.247384, -0.597239, 0.353553, 0.247384, -0.597239, -0.353553, 0.517982, -1.25052, -0.353553, 1.25052, -0.517982, 0.353553, 0.597239, -0.247384, 0.353553, 0.597239, -0.247384, -0.353553, 1.25052, -0.517982, -0.353553
#define TORUS_FACES { {4, 5, 1, 0}, {5, 6, 2, 1}, {6, 7, 3, 2}, {7, 4, 0, 3}, {8, 9, 5, 4}, {9, 10, 6, 5}, {10, 11, 7, 6}, {11, 8, 4, 7}, {12, 13, 9, 8}, {13, 14, 10, 9}, {14, 15, 11, 10}, {15, 12, 8, 11}, {16, 17, 13, 12}, {17, 18, 14, 13}, {18, 19, 15, 14}, {19, 16, 12, 15}, {20, 21, 17, 16}, {21, 22, 18, 17}, {22, 23, 19, 18}, {23, 20, 16, 19}, {24, 25, 21, 20}, {25, 26, 22, 21}, {26, 27, 23, 22}, {27, 24, 20, 23}, {28, 29, 25, 24}, {29, 30, 26, 25}, {30, 31, 27, 26}, {31, 28, 24, 27}, {0, 1, 29, 28}, {1, 2, 30, 29}, {2, 3, 31, 30}, {3, 0, 28, 31} }

namespace subdivision_matrix
{

	MatrixXX_t shepard( const MatrixX3_t& positions, const MatrixX3_t& handle_positions )
	{	
		int row_num = positions.rows();
		int col_num = handle_positions.rows();
	
		MatrixXX_t weights( row_num, col_num );
	
		for ( int i = 0; i < row_num; i++ ) {
			int flag = -1;
			for ( int j = 0; j < col_num; j++ ) {
				vertex_t diff = positions.row( i ) - handle_positions.row( j );
				weights( i, j ) = diff.dot( diff );
				// if a position is on the handle, mark the handle's index
				if ( weights( i, j ) <= EPS )
					flag = j;
			}
			// if a position is on any handle, make its weight 1.0.
			if ( flag >= 0 ) {
				for ( int j = 0; j < col_num; j++ ) {
					weights( i, j ) = ( flag == j ? 1.0 : 0.0 );
				}
				continue;
			}
		
			for ( int j = 0; j < col_num; j++ ) {
				weights( i, j ) = 1.0 / weights( i, j );
			}
			weights.row( i ) /= weights.row( i ).sum();
		
		}
	
		return weights;
	}

//-------------------------------------- PREPARE METHODS --------------------------------------
	void prepare( const subdivision_control_mesh& mesh, const MatrixX3_t& handle_positions, const char* weight_function,  
			std::vector< MatrixX4_t > & W_Array )
	{
		/// 1 Get the sparse subdivision coefficients at many places around the mesh.
		/// 2 Compute the weights and areas for every point around the mesh.
		
		/// 1
		std::vector< real_t > us, vs;
		createUVs( resolution, resolution, us, vs );
	
		SparseMatrix_t M_matrices, Du_matrices, Dv_matrices;
		compute_subdivision_coefficients_for_mesh(
			mesh.vs.rows(),
			mesh.faces, 
			us, vs, 
			M_matrices, &Du_matrices, &Dv_matrices );	

	
		std::cout << M_matrices.rows() << ' ' << M_matrices.cols() << std::endl;
	
		prepare( mesh.vs, M_matrices, Du_matrices, Dv_matrices, handle_positions, weight_function, W_Array );
	}
    
    void prepare( const MatrixX3_t& vs, const SparseMatrix_t& M_matrices,
		const MatrixXX_t& weights, std::vector< MatrixX4_t > & W_Array )
    {
        ///  areas is a column vector whose number of rows equals to the number of all positions.
	    DiagonalMatrix_t areas( M_matrices.rows() );
	    // std::cout << "areas size: " << M_matrices.rows() << '\n';
	    areas.setIdentity();
	    prepare( vs, M_matrices, areas, weights, W_Array );
    }
    
	void prepare( const MatrixX3_t& vs, const SparseMatrix_t& M_matrices, const DiagonalMatrix_t& areas,
		const MatrixXX_t& weights, std::vector< MatrixX4_t > & W_Array )
	{
		/// Compute the weights and areas for every point around the mesh.
		///  areas is a column vector whose number of rows equals to the number of all positions.
		
		/// build system.
		const MatrixXX_t system = ( M_matrices.transpose() * areas * M_matrices ).eval();
	
		/// set up LLT solver and examine if the system can be inverted.
		Eigen::LLT< MatrixXX_t > llt;
		llt.compute( system );
		if ( Eigen::Success == llt.info () )
			; //std::cout << "succeed." << std::endl;
		else if ( Eigen::NumericalIssue == llt.info () )
			std::cout << "numerical issue." << std::endl;
	
	
		/// set up an identity matrix for compute the inverse of system.	
		MatrixX4_t controls( vs.rows(), 4 );
		controls.leftCols( 3 ) = vs;
		controls.col( 3 ).setOnes();
	    
	    // std::cout << "M_matrices rows: " << M_matrices.rows() << '\n';
	    // std::cout << "M_matrices cols: " << M_matrices.cols() << '\n';
	    // std::cout << "weights size: " << weights.col(0).asDiagonal().rows() << '\n';
	    
	    for( int k = 0; k < weights.cols(); k++ ) {
		    const MatrixX4_t extra = ( controls.transpose() * M_matrices.transpose() * areas * weights.col(k).asDiagonal() * M_matrices ).transpose();
			W_Array.push_back( llt.solve( extra ) );
		}
	
		return;
	}
	void prepare( const MatrixX3_t& vs, const SparseMatrix_t& M_matrices, const DiagonalMatrix_t& areas,
		const MatrixX3_t handle_positions, const char* weight_function,  
		std::vector< MatrixX4_t > & W_Array )
	{
		/// Compute the weights and areas for every point around the mesh.
		///  areas is a column vector whose number of rows equals to the number of all positions.
		
		MatrixX3_t positions( M_matrices.rows(), 3 );
		positions = M_matrices * vs;
	    
		/// weights are a dense matrix whose number of rows equales to the number of all positions,
		/// and the number of columns equals to the number of handles.
		MatrixXX_t weights;
		if( strcmp( weight_function, "harmonic" ) == 0) {
	//         weights = harmonic( mesh, positions, handle_positions );

		} else if ( strcmp( weight_function, "shepard" ) == 0 ) {
			 weights = shepard( positions, handle_positions );
		}
	    
	    prepare( vs, M_matrices, areas, weights, W_Array );
    }
	
    void prepare( const MatrixX3_t& vs, const SparseMatrix_t& M_matrices, const SparseMatrix_t& Du_matrices, const SparseMatrix_t& Dv_matrices,
		const MatrixX3_t& handle_positions, const char* weight_function,  
		std::vector< MatrixX4_t > & W_Array )
	{
	    ///  areas is a column vector whose number of rows equals to the number of all positions.
		DiagonalMatrix_t areas( M_matrices.rows() );
		MatrixX3_t du_list = Du_matrices * vs;
		MatrixX3_t dv_list = Dv_matrices * vs;
		Eigen::Matrix< real_t, 1, 3 > vec1, vec2;
		for( int i = 0; i < M_matrices.rows(); ++i )
		{
			vec1 = du_list.row( i );
			vec2 = dv_list.row( i );
			// Divide by the number of rows, which is like multiplying by dt.
			areas.diagonal()[i] = ( vec1.cross( vec2 ).norm() ) / M_matrices.rows();
		}
		
		prepare( vs, M_matrices, areas, handle_positions, weight_function, W_Array );
	}
	void prepare( const MatrixX3_t& vs, const SparseMatrix_t& M_matrices,
		const MatrixX3_t& handle_positions, const char* weight_function,
		std::vector< MatrixX4_t > & W_Array )
	{
	    ///  areas is a column vector whose number of rows equals to the number of all positions.
	    DiagonalMatrix_t areas( M_matrices.rows() );
	    areas.setIdentity();
	    prepare( vs, M_matrices, areas, handle_positions, weight_function, W_Array );
	}

	void naive_prepare( const MatrixX3_t& vs, const MatrixX3_t& handle_positions, 
		std::vector< MatrixX4_t > & W_Array )
	{	
		MatrixXX_t weights = shepard( vs, handle_positions );
		
		// build the homogeneous coordinates for control points.
		MatrixX4_t controls( vs.rows(), 4 );
		controls.leftCols( 3 ) = vs;
		controls.col( 3 ).setOnes();
		for( int k = 0; k < weights.cols(); k++ ) {
			MatrixXX_t wi( weights.rows(), 4 );
			for ( int i = 0; i < 4; i++ )
				wi.col(i) = weights.col(k);
					
			W_Array.push_back( controls.cwiseProduct( wi ) );
		}
// 		std::cout << "naive prepare finished\n";
	}	

//-------------------------------------- PREPARE METHODS END--------------------------------------

	void solve(  const std::vector< MatrixX4_t >& W_Array, 
		const std::vector< transform_t >& transforms,
		MatrixX3_t & result )
	{
	
		assert( W_Array.size() == transforms.size() );
	
		MatrixXX_t rhs( W_Array[0] );
		rhs.setZero();
		 
		for( int i = 0; i < transforms.size(); i++ )
		{	
			rhs += ( W_Array[i] * transforms[i].transpose() );
		}

		result = rhs.leftCols( 3 );
	
		return;
	}
	
	
	void compute_target_surface( const MatrixX3_t& vs, const SparseMatrix_t& M_matrices, 
		const MatrixX3_t& handle_positions, const std::vector< transform_t >& transforms, 
		MatrixX3_t& result )
	{		
		MatrixX3_t positions( M_matrices.rows(), 3 );
		positions = M_matrices * vs;
		
		std::vector< MatrixX4_t > W_Array;
		naive_prepare( positions, handle_positions, W_Array );
		
		solve( W_Array, transforms, result ); 
		
// 		std::cout << "target surface computed.\n";
	}
	
	VectorX_t piecewise_distance( const MatrixX3_t& target_surface,  const MatrixX3_t& deformed_surface )
	{
		MatrixX3_t diff_mat = deformed_surface - target_surface;
		VectorX_t result = diff_mat.rowwise().norm();
		
		return result;
	}
	
	real_t norm_of_piecewise_distance( const MatrixX3_t& target_surface,  const MatrixX3_t& deformed_surface )
	{
		MatrixX3_t diff_mat = deformed_surface - target_surface;
		return diff_mat.norm();
	}

	
	real_t hausdorff_distance( const MatrixX3_t& target_surface,  const MatrixX3_t& deformed_surface )
	{
		assert( target_surface.rows() == deformed_surface.rows() );
		int num = target_surface.rows();
		
		VectorX_t min_diffs( num, 1 );
		for( int i = 0; i < num; i++ ) {
			VectorX_t diffs = ( deformed_surface.rowwise() - target_surface.row(i) ).rowwise().norm();
			min_diffs(i) = diffs.minCoeff();
		}
		
		return min_diffs.maxCoeff();
	}
}

#ifdef SUBDIVISION_ENGINE_MAIN
void usage( const char* argv0 )
{
    std::cerr << "Usage: " << argv0 << " test_subdivision_engine\n";
    exit( -1 );
}

void print_mesh( const subdivision_matrix::subdivision_control_mesh& mesh );

void test( const char *type )
{
	subdivision_matrix::MatrixX3_t handle_positions;
	handle_positions.resize(2,3);
	handle_positions << 0.8, 0.8, 0.7,
						0.2, -0.3, -0.2;
	
	// make a mesh
	subdivision_matrix::subdivision_control_mesh mesh;
	
	// use the cube! has problem!
// 	mesh.vs.resize( 8, 3 );
// 	mesh.vs << CUBE_VS;
// 	mesh.faces = CUBE_FACES;

	mesh.vs.resize( 32, 3 );
	mesh.vs << TORUS_VS;
	mesh.faces = TORUS_FACES;
	
	
	std::vector< subdivision_matrix::MatrixX4_t > W_Array;
	if ( strcmp( type, "naive") == 0 )
		subdivision_matrix::naive_prepare( mesh.vs, handle_positions, W_Array );
	else
    	subdivision_matrix::prepare( mesh, handle_positions, "shepard", W_Array );
    
    for (auto W_i : W_Array )
    	std::cout << W_i << std::endl;
    
   	std::vector< subdivision_matrix::transform_t > transforms;
   	subdivision_matrix::transform_t t1, t2;
   	t1.setIdentity(); 
   	t2.setIdentity();
   	transforms.push_back( t1 );
   	transforms.push_back( t2 );
   	
	subdivision_matrix::MatrixX3_t new_controls;
	subdivision_matrix::solve( W_Array, transforms, new_controls );
	
	std::cout << "result\n" << new_controls << std::endl;
}

int main( int argc, char* argv[] )
{
    using std::string;
    
//     if( argc != 1 ) usage( argv[0] );
    
    if( argc == 1 )
		test( "ours" );
	else
		test( argv[1]);

    return 0;
}
#endif
void print_mesh( const subdivision_matrix::subdivision_control_mesh& mesh )
{
	std::cout << "vertices\n";
	std::cout << mesh.vs << std::endl;
	std::cout << "faces\n";
	int face_offset = 0;
    for( const auto& num_face_vertices : mesh.faces.first )
    {
        for( int fvi = 0; fvi < num_face_vertices; ++fvi )
        {
            std::cout << mesh.faces.second.at( face_offset + fvi ) << ' ';
        }
        face_offset += num_face_vertices;
        std::cout << '\n';
    }
}

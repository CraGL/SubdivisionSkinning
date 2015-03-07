#ifndef __subdivision_engine_h__
#define __subdivision_engine_h__

#include "subdivision_matrices_eigen.h"
#include "subdivision_matrices_osd_eigen.h"
#include <Eigen/Dense>
#include <Eigen/Cholesky>
#include "Timing.hpp"

namespace subdivision_matrix
{	
	// Declare the mesh structure.
	struct subdivision_control_mesh
	{
		// This should be an N-by-3 dense matrix.
		MatrixX3_t vs;
		// This is a vector of faces.
		faces_t faces;
	};
	
	typedef Eigen::DiagonalMatrix< real_t, Eigen::Dynamic > DiagonalMatrix_t;
	
	void prepare( const subdivision_control_mesh& mesh, const MatrixX3_t& handle_positions, const char * weight_function,  
			std::vector< MatrixX4_t > & W_Array );
	
	void prepare( const MatrixX3_t& vs, const SparseMatrix_t& M_matrices, const SparseMatrix_t& Du_matrices, const SparseMatrix_t& Dv_matrices,
		const MatrixX3_t& handle_positions, const char* weight_function,  
		std::vector< MatrixX4_t > & W_Array );
	
	void prepare( const MatrixX3_t& vs, const SparseMatrix_t& M_matrices,
		const MatrixX3_t& handle_positions, const char* weight_function,
		std::vector< MatrixX4_t > & W_Array );
	
	void prepare( const MatrixX3_t& vs, const SparseMatrix_t& M_matrices, const DiagonalMatrix_t& areas,
		const MatrixXX_t& weights, std::vector< MatrixX4_t > & W_Array );
	
	void prepare( const MatrixX3_t& vs, const SparseMatrix_t& M_matrices,
		const MatrixXX_t& weights, std::vector< MatrixX4_t > & W_Array );
	
	void solve(  const std::vector< MatrixX4_t >& W_Array, const std::vector< transform_t >& transforms, MatrixX3_t & result );
	
	// naive approach
	void naive_prepare( const MatrixX3_t& vs, const MatrixX3_t& handle_positions, 
		std::vector< MatrixX4_t > & W_Array );
		
	void compute_target_surface( const MatrixX3_t& vs, const SparseMatrix_t& M_matrices, 
		const MatrixX3_t& handle_positions, const std::vector< transform_t >& transforms, 
		MatrixX3_t& result );
		
	VectorX_t piecewise_distance( const MatrixX3_t& target_surface,  const MatrixX3_t& deformed_surface );

	real_t norm_of_piecewise_distance( const MatrixX3_t& target_surface,  const MatrixX3_t& deformed_surface );
	
	real_t hausdorff_distance( const MatrixX3_t& target_surface,  const MatrixX3_t& deformed_surface );
}

#endif

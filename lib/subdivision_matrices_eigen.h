#ifndef __subdivision_matrices_eigen_h__
#define __subdivision_matrices_eigen_h__

#include "subdivision_matrices.h"
#include <Eigen/Core>
#include <Eigen/Sparse>

namespace subdivision_matrix
{

/*
Given the number of vertices in the mesh 'num_vertices'
and a vector of faces 'faces', each element of which is a vector of vertex indices, and
two same-length vectors 'us' and 'vs' where the ( us[i], vs[i] ) are the uv locations
with which to sample every face,
fills 'positions_out' with a matrix such that
the matrix multiplication of 'positions_out' times a num_vertices-by-K matrix of control points
yields the positions.

The optional parameters 'du_out' and 'dv_out', if specified, are similar to 'positions_out'
except that the matrix multiplication results in du and dv vectors.
*/
void compute_subdivision_coefficients_for_mesh(
    int num_vertices,
    const faces_t& faces,
    const std::vector< real_t >& us,
    const std::vector< real_t >& vs,
    SparseMatrix_t& positions_out,
    SparseMatrix_t* du_out = nullptr,
    SparseMatrix_t* dv_out = nullptr
    );
void compute_subdivision_coefficients_for_mesh(
    int num_vertices,
    OpenSubdiv::HbrMesh<OpenSubdiv::FarStencilFactoryVertex>* mesh,
    const std::vector< real_t >& us,
    const std::vector< real_t >& vs,
    SparseMatrix_t& positions_out,
    SparseMatrix_t* du_out = nullptr,
    SparseMatrix_t* dv_out = nullptr
    );

}

#ifdef SUBDIVISION_MATRICES_HEADER_ONLY
#include "subdivision_matrices_eigen.cpp"
#endif

#endif /* __subdivision_matrices_eigen_h__ */

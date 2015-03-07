#ifndef __subdivision_matrices_osd_eigen_h__
#define __subdivision_matrices_osd_eigen_h__

#include "subdivision_matrices_osd.h"
#include <Eigen/Core>
#include <Eigen/Sparse>

namespace subdivision_matrix
{

/*
Given the number of vertices in the mesh 'num_vertices'
and a vector of faces 'faces', each element of which is a vector of vertex indices, and
a positive integer 'level' indicating the level of refinement,
fills 'positions_out' with a matrix such that
the matrix multiplication of 'positions_out' times a num_vertices-by-K matrix of control points
yields the positions.

The optional parameter 'faces_out', if specified, is a sequence of quad faces
obtained by subdivision, where each face is four indices into 'positions_out':
    face0_vertex0 face0_vertex1 face0_vertex2 face0_vertex3 face1_vertex0 face1_vertex1 face1_vertex2 face1_vertex3 ...
*/
void compute_subdivision_coefficients_for_mesh(
    int num_vertices,
    const faces_t& faces,
    const int level,
    SparseMatrix_t& positions_out,
    std::vector< index_t >* quad_faces_out = nullptr
    );
void compute_subdivision_coefficients_for_mesh(
    int num_vertices,
    OpenSubdiv::HbrMesh<OpenSubdiv::OsdVertex>* mesh,
    const int level,
    SparseMatrix_t& positions_out,
    std::vector< index_t >* quad_faces_out = nullptr
    );

}

#ifdef SUBDIVISION_MATRICES_HEADER_ONLY
#include "subdivision_matrices_osd_eigen.cpp"
#endif

#endif /* __subdivision_matrices_osd_eigen_h__ */

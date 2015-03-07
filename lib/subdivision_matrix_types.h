#ifndef __subdivision_matrix_types_h__
#define __subdivision_matrix_types_h__

#include <Eigen/Core>
#include <Eigen/Sparse>
#include <vector>
#include <utility> // pair

namespace subdivision_matrix
{

typedef float real_t;
// typedef double long_real_t;
typedef int index_t;
typedef std::vector< std::pair< index_t, real_t > > sparse_vector_t;

// typedef Eigen::Matrix< real_t, Eigen::Dynamic, Eigen::Dynamic > Matrix_t;
typedef Eigen::SparseMatrix< real_t, Eigen::RowMajor > SparseMatrix_t;

typedef Eigen::Matrix< real_t, 1, 3> vertex_t;

typedef Eigen::Matrix< real_t, Eigen::Dynamic, Eigen::Dynamic > MatrixXX_t;
typedef Eigen::Matrix< real_t, Eigen::Dynamic, 4 > MatrixX4_t;
typedef Eigen::Matrix< real_t, Eigen::Dynamic, 3 > MatrixX3_t;
typedef Eigen::Matrix< real_t, Eigen::Dynamic, 1 > VectorX_t;

typedef std::pair< std::vector< int >, std::vector< int > > faces_t;
typedef Eigen::Matrix< real_t, 4, 4 > transform_t;

}

#endif /* __subdivision_matrix_types_h__ */

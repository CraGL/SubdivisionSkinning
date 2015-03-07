#ifndef __subdivision_matrices_h__
#define __subdivision_matrices_h__

#include "subdivision_matrix_types.h"
#include <vector>

// -I/path/to/OpenSubDiv
#include <osd/mesh.h>
#include <far/stencilTablesFactory.h>

namespace subdivision_matrix
{

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
    std::vector< sparse_vector_t >* du_out = nullptr,
    std::vector< sparse_vector_t >* dv_out = nullptr
    );
// The same as above, except with an HbrMesh pointer
void compute_subdivision_coefficients_for_mesh(
    int num_vertices,
    OpenSubdiv::HbrMesh<OpenSubdiv::FarStencilFactoryVertex>* mesh,
    const std::vector< real_t >& us,
    const std::vector< real_t >& vs,
    std::vector< sparse_vector_t >& positions_out,
    std::vector< sparse_vector_t >* du_out = nullptr,
    std::vector< sparse_vector_t >* dv_out = nullptr
    );

/*
Given a desired number of "u" and "v" parameters sampling the range [0,1]x[0,1],
returns in the output vectors 'us' and 'vs' coordinates sampling the range
in equal intervals.
*/
void createUVs( int num_u, int num_v, std::vector< real_t >& us, std::vector< real_t >& vs );

}

#ifdef SUBDIVISION_MATRICES_HEADER_ONLY
#include "subdivision_matrices.cpp"
#endif

#endif /* __subdivision_matrices_h__ */

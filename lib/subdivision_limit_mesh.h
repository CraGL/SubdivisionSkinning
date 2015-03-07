#ifndef __subdivision_limit_mesh_h__
#define __subdivision_limit_mesh_h__

#include "subdivision_matrices.h"
#include "subdivision_matrices_eigen.h"

namespace subdivision_limit_mesh
{

typedef subdivision_matrix::real_t real_t;
/*
Given the vertices of the mesh 'vertices',
an HbrMesh pointer 'mesh',
and the desired number of "u" and "v" parameters sampling the range [0,1]x[0,1],
returns a quad mesh for the limit surface sampled at least according
to 'num_u' and 'num_v';
the resulting vertices are placed in 'vertices_out' and
the resulting faces (vectors of CCW faces 0-indexed into vertices_out) in 'faces_out'.

NOTE: If the same 'num_u' and 'num_v' are passed to createUVs() followed by compute_subdivision_coefficients*(),
      then the mesh returned by this function will include all of the same limit positions (and more).
*/
void compute_quad_limit_mesh_for_mesh(
    const subdivision_matrix::MatrixX3_t& vertices,
    OpenSubdiv::HbrMesh<OpenSubdiv::FarStencilFactoryVertex>* mesh,
    int num_u, int num_v,
    std::vector< std::vector< real_t > >& vertices_out,
    std::vector< std::vector< int > >& faces_out
    );
}

#ifdef SUBDIVISION_LIMIT_MESH_HEADER_ONLY
#include "subdivision_limit_mesh.cpp"
#endif

#endif /* __subdivision_limit_mesh_h__ */

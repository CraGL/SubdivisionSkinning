#ifndef __subdivision_skinning_h__
#define __subdivision_skinning_h__

#include <vector>
#include <Eigen/Core>

// -I/path/to/OpenSubDiv
#include <osd/mesh.h>
#include <far/stencilTablesFactory.h>

#include "subdivision_engine.h"

namespace subdivision_skinning
{
    typedef float real_t;
    // A sequence of vertex values: x0 y0 z0 x1 y1 z1 ...
    typedef std::vector< real_t > vertices_t;
    // Each face in a faces_t is a sequence of vertex indices.
    using subdivision_matrix::faces_t;
    // An H-by-3 matrix with an x,y,z position for each handle.
    typedef Eigen::MatrixXf handles_t;
    // A transform is a 4x4 matrix.
    typedef Eigen::Matrix4f transform_t;
    
    // Forward declare the engine structure.
    struct engine_t;
    
    // Precomputes and returns a new engine_t for the given subdivision surface.
    // The resulting pointer can be passed to update_control_vertices_given_handle_transforms().
    // To destroy the engine_t and free the memory associated with it, call delete_engine_t().
    // The optional HbrMesh* provides tagging information.
    engine_t* new_engine_for_control_mesh( const vertices_t& vertices, const faces_t& faces, const handles_t& handle_positions, const char* weight_function, OpenSubdiv::HbrMesh<OpenSubdiv::FarStencilFactoryVertex>* mesh = 0 );
    engine_t* new_engine_for_control_mesh( const vertices_t& vertices, const faces_t& faces, const handles_t& handle_positions, const char* weight_function, OpenSubdiv::HbrMesh<OpenSubdiv::OsdVertex>* mesh, std::vector< subdivision_matrix::index_t >* quad_faces_out = 0 );
    engine_t* new_engine_for_control_mesh( const vertices_t& vertices, const faces_t& faces, const subdivision_matrix::SparseMatrix_t& M_matrices, const subdivision_matrix::MatrixXX_t& weights );
    // Given an existing 'engine_t', modifies it to use the new 'handle_positions' and 'weight_function'.
    // Subsequent calls to update_control_vertices_given_handle_transforms() must pass
    // 'handle_positions.rows()' transforms.
    void set_new_handle_positions( engine_t* engine, const handles_t& handle_positions, const char* weight_function );
    // Given a set of transform_t's, one for each handle_position,
    // returns in 'vertices_out' updated positions for the control points of the subdivision surface.
    void update_control_vertices_given_handle_transforms( const engine_t* engine, const std::vector< transform_t >& transforms, vertices_t& vertices_out );
    // update with test and all kinds of comparisons.
    void update_control_vertices_given_handle_transforms( const engine_t* engine, const std::vector< transform_t >& transforms, vertices_t& vertices_out, const char* approach, vertices_t& g_targetPositions, vertices_t& g_targetColors );
    // Frees the memory associated with an engine_t.
    void delete_engine_t( engine_t* engine );
    // naive approach wrapper
    
    // Upon return, replaces 'limit_vertices_out' with limit vertex positions
    // for all of our sample locations, based on the initial control mesh.
    // This is useful for an experiment.
    void limit_vertices_for_initial_sample_locations( const engine_t* engine, vertices_t& limit_vertices_out );
}

#endif /* __subdivision_skinning_h__ */

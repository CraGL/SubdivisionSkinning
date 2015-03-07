#include "subdivision_skinning.h"
#include "subdivision_matrices_osd_eigen.h"

#include <algorithm> // std::copy

// Static definitions.
namespace
{
    const int kDefaultLevel = 4;
    
    struct subdivision_evaluator_t
    {
        subdivision_skinning::vertices_t vertices;
        subdivision_skinning::faces_t faces;
        subdivision_matrix::SparseMatrix_t M_matrices;
        std::vector< subdivision_matrix::index_t > quad_faces;
    };
    
    struct subdivision_skinning_engine_t
    {
        subdivision_matrix::SparseMatrix_t M_matrices;
        int num_transforms;
        subdivision_skinning::engine_t* engine;
        
        subdivision_skinning_engine_t() : engine( 0 ) {}
        ~subdivision_skinning_engine_t()
        {
            if( engine ) subdivision_skinning::delete_engine_t( engine );
            engine = 0;
        }
    };
}


extern "C" {

#include "subdivision_skinning_wrapper.h"

void* new_subdivision_evaluator( int num_vertices, const subdivision_evaluator_real_t* vertices, int num_faces, const int* faces, int level )
{
    assert( num_vertices > 0 );
    assert( num_faces > 0 );
    assert( vertices );
    assert( faces );
    assert( level > 0 );
    
    // Create the result object.
    subdivision_evaluator_t* result = new subdivision_evaluator_t;
    
    // Copy the vertices.
    result->vertices.resize( num_vertices*3 );
    std::copy( vertices, vertices + num_vertices*3, result->vertices.begin() );
    
    // Copy the faces.
    result->faces.first.reserve( num_faces );
    result->faces.second.reserve( num_faces*4 );
    for( int fi = 0; fi < num_faces; ++fi )
    {
        // Number of vertices in the face:
        result->faces.first.push_back( 0 );
        for( int vi = 0; vi < 4; ++vi )
        {
            const int vertex_index = faces[ 4*fi + vi ];
            if( -1 != vertex_index )
            {
                result->faces.second.push_back( vertex_index );
                result->faces.first.back() += 1;
            }
        }
        assert( result->faces.first.back() == 3 || result->faces.first.back() == 4 );
    }
    
    // Create M_matrices and quads.
    subdivision_matrix::compute_subdivision_coefficients_for_mesh(
        num_vertices,
        result->faces,
        level,
        result->M_matrices,
        &result->quad_faces
        );
    
    return result;
}
void delete_subdivision_evaluator( void* subdivision_evaluator )
{
    delete static_cast< subdivision_evaluator_t* >( subdivision_evaluator );
}

int num_refined_quad_faces_of_subdivision_evaluator( const void* subdivision_evaluator_in )
{
    const subdivision_evaluator_t* subdivision_evaluator = static_cast< const subdivision_evaluator_t* >( subdivision_evaluator_in );
    assert( subdivision_evaluator );
    
    assert( subdivision_evaluator->quad_faces.size() % 4 == 0 );
    return subdivision_evaluator->quad_faces.size() / 4;
}
void get_refined_quad_faces_of_subdivision_evaluator( const void* subdivision_evaluator_in, int* refined_quad_faces_out )
{
    const subdivision_evaluator_t* subdivision_evaluator = static_cast< const subdivision_evaluator_t* >( subdivision_evaluator_in );
    assert( subdivision_evaluator );
    
    assert( refined_quad_faces_out );
    std::copy( subdivision_evaluator->quad_faces.begin(), subdivision_evaluator->quad_faces.end(), refined_quad_faces_out );
}

int num_refined_vertices_of_subdivision_evaluator( const void* subdivision_evaluator_in )
{
    const subdivision_evaluator_t* subdivision_evaluator = static_cast< const subdivision_evaluator_t* >( subdivision_evaluator_in );
    assert( subdivision_evaluator );
    
    return subdivision_evaluator->M_matrices.rows();
}
void get_refined_vertices_of_subdivision_evaluator( const void* subdivision_evaluator_in, subdivision_evaluator_real_t* refined_vertices_out )
{
    const subdivision_evaluator_t* subdivision_evaluator = static_cast< const subdivision_evaluator_t* >( subdivision_evaluator_in );
    assert( subdivision_evaluator );
    assert( refined_vertices_out );
    
    assert( subdivision_evaluator->vertices.size() == subdivision_evaluator->M_matrices.cols()*3 );
    const Eigen::Map< const Eigen::Matrix< subdivision_skinning::vertices_t::value_type, Eigen::Dynamic, 3, Eigen::RowMajor > > vs( &subdivision_evaluator->vertices[0], subdivision_evaluator->M_matrices.cols(), 3 );
    Eigen::Map< Eigen::Matrix< subdivision_evaluator_real_t, Eigen::Dynamic, 3, Eigen::RowMajor > > out( refined_vertices_out, subdivision_evaluator->M_matrices.rows(), 3 );
    out = ( subdivision_evaluator->M_matrices * vs.cast< subdivision_matrix::real_t >() ).cast< subdivision_evaluator_real_t >();
}
void get_refined_vertices_of_subdivision_evaluator_with_control_vertices( const void* subdivision_evaluator_in, int vertex_dimension, const subdivision_evaluator_real_t* control_vertices, subdivision_evaluator_real_t* refined_vertices_out )
{
    const subdivision_evaluator_t* subdivision_evaluator = static_cast< const subdivision_evaluator_t* >( subdivision_evaluator_in );
    assert( subdivision_evaluator );
    assert( vertex_dimension > 0 );
    assert( control_vertices );
    assert( refined_vertices_out );
    
    const Eigen::Map< const Eigen::Matrix< subdivision_evaluator_real_t, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor > > vs( control_vertices, subdivision_evaluator->M_matrices.cols(), vertex_dimension );
    Eigen::Map< Eigen::Matrix< subdivision_evaluator_real_t, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor > > out( refined_vertices_out, subdivision_evaluator->M_matrices.rows(), vertex_dimension );
    out = ( subdivision_evaluator->M_matrices * vs.cast< subdivision_matrix::real_t >() ).cast< subdivision_evaluator_real_t >();
}

void* new_subdivision_skinning_engine( const void* subdivision_evaluator_in, int num_transforms, const subdivision_evaluator_real_t* weights_in )
{
    assert( num_transforms > 0 );
    assert( weights_in );
    assert( subdivision_evaluator_in );
    const subdivision_evaluator_t* subdivision_evaluator = static_cast< const subdivision_evaluator_t* >( subdivision_evaluator_in );
    
    subdivision_skinning_engine_t* result = new subdivision_skinning_engine_t;
    result->num_transforms = num_transforms;
    result->M_matrices = subdivision_evaluator->M_matrices;
    const subdivision_matrix::MatrixXX_t weights = Eigen::Map< const Eigen::Matrix< subdivision_evaluator_real_t, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor > >( weights_in, subdivision_evaluator->M_matrices.rows(), num_transforms ).cast< subdivision_matrix::real_t >();
    result->engine = subdivision_skinning::new_engine_for_control_mesh( subdivision_evaluator->vertices, subdivision_evaluator->faces, subdivision_evaluator->M_matrices, weights );
    
    return result;
}
void delete_subdivision_skinning_engine( void* subdivision_skinning_engine )
{
    delete static_cast< subdivision_skinning_engine_t* >( subdivision_skinning_engine );
}

void compute_control_mesh_vertices_given_transforms_for_subdivision_skinning_engine( const void* subdivision_skinning_engine_in, const subdivision_evaluator_real_t* transforms_in, subdivision_evaluator_real_t* control_mesh_vertices_out )
{
    const subdivision_skinning_engine_t* engine = static_cast< const subdivision_skinning_engine_t* >( subdivision_skinning_engine_in );
    assert( engine );
    assert( engine->engine );
    assert( transforms_in );
    assert( control_mesh_vertices_out );
    
    // Copy the transforms.
    std::vector< subdivision_skinning::transform_t > transforms( engine->num_transforms );
    for( int i = 0; i < engine->num_transforms; ++i )
    {
        transforms.at(i) = Eigen::Map< const Eigen::Matrix< subdivision_evaluator_real_t, 4, 4, Eigen::RowMajor > >( transforms_in + 16*i, 4, 4 );
    }
    
    // Solve for the new control mesh vertices.
    subdivision_skinning::vertices_t control_mesh_vertices;
    update_control_vertices_given_handle_transforms( engine->engine, transforms, control_mesh_vertices );
    
    // Copy to the output.
    assert( control_mesh_vertices.size() == engine->M_matrices.cols()*3 );
    // Copy (and convert floating point formats if needed).
    Eigen::Map< Eigen::Matrix< subdivision_evaluator_real_t, Eigen::Dynamic, 3, Eigen::RowMajor > >( control_mesh_vertices_out, engine->M_matrices.cols(), 3 ) =
        Eigen::Map< Eigen::Matrix< subdivision_skinning::vertices_t::value_type, Eigen::Dynamic, 3, Eigen::RowMajor > >( &control_mesh_vertices[0], engine->M_matrices.cols(), 3 );
}

int num_refined_vertices_of_subdivision_skinning_engine( const void* subdivision_skinning_engine )
{
    const subdivision_skinning_engine_t* engine = static_cast< const subdivision_skinning_engine_t* >( subdivision_skinning_engine );
    assert( engine );
    
    return engine->M_matrices.rows();
}
void get_refined_vertices_of_subdivision_skinning_engine_with_control_vertices( const void* subdivision_skinning_engine, int vertex_dimension, const subdivision_evaluator_real_t* control_vertices, subdivision_evaluator_real_t* refined_vertices_out )
{
    const subdivision_skinning_engine_t* engine = static_cast< const subdivision_skinning_engine_t* >( subdivision_skinning_engine );
    assert( engine );
    assert( vertex_dimension > 0 );
    assert( control_vertices );
    assert( refined_vertices_out );
    
    const Eigen::Map< const Eigen::Matrix< subdivision_evaluator_real_t, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor > > vs( control_vertices, engine->M_matrices.cols(), vertex_dimension );
    Eigen::Map< Eigen::Matrix< subdivision_evaluator_real_t, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor > > out( refined_vertices_out, engine->M_matrices.rows(), vertex_dimension );
    out = ( engine->M_matrices * vs.cast< subdivision_matrix::real_t >() ).cast< subdivision_evaluator_real_t >();
}

}

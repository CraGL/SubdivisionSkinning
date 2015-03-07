// c++ -g -std=c++11 -I../OpenSubdiv/opensubdiv -I../OpenSubdiv/regression subdivision_matrices_eigen.cpp subdivision_matrices.cpp -o subdivision_matrices_eigen -DSUBDIVISION_MATRICES_EIGEN_MAIN -I/usr/local/include/eigen3

#include "subdivision_matrices_eigen.h"

#include <cassert>

namespace 
{
/*
Given the number of vertices 'num_vertices',
a vector of sparse_vector_t 'sparse_vectors' containing a sequence of ( index, coefficient ) pairs,
and an output Eigen::SparseMatrix 'matrix',
fills 'matrix' such that the matrix multiplication of 'matrix' times a num_vertices-by-K matrix of control points
yields the positions specified by 'sparse_vectors'.
*/
void convert_vector_of_sparse_vectors_to_matrix( int num_vertices, const std::vector< subdivision_matrix::sparse_vector_t >& sparse_vectors, subdivision_matrix::SparseMatrix_t& matrix )
{
    assert( num_vertices > 0 );
    
    // Clear the output matrix.
    matrix.resize( sparse_vectors.size(), num_vertices );
    matrix.setZero();
    
    // We will fill the matrix with a vector of triplets.
    std::vector< Eigen::Triplet< subdivision_matrix::real_t > > triplets;
    
    // Count the number of triplets we will need.
    int nnz = 0;
    for( const auto& sp : sparse_vectors ) nnz += sp.size();
    triplets.reserve( nnz );
    
    // Convert 'sparse_vectors' to triplets.
    for( int row = 0; row < sparse_vectors.size(); ++row )
    {
        for( const auto p : sparse_vectors[row] )
        {
            triplets.push_back( { row, p.first, p.second } );
        }
    }
    
    matrix.setFromTriplets( triplets.begin(), triplets.end() );
}
}

namespace subdivision_matrix
{

void compute_subdivision_coefficients_for_mesh(
    int num_vertices,
    const faces_t& faces,
    const std::vector< real_t >& us,
    const std::vector< real_t >& vs,
    SparseMatrix_t& positions_out,
    SparseMatrix_t* du_out,
    SparseMatrix_t* dv_out
    )
{
    // Intermediary output vectors.
    std::vector< sparse_vector_t > position_coeffs;
    std::vector< sparse_vector_t > du_coeffs;
    std::vector< sparse_vector_t > dv_coeffs;
    
    // Call the function we're wrapping.
    compute_subdivision_coefficients_for_mesh(
        num_vertices,
        faces,
        us, vs,
        position_coeffs,
        du_out ? &du_coeffs : nullptr,
        dv_out ? &dv_coeffs : nullptr
        );
    
    convert_vector_of_sparse_vectors_to_matrix( num_vertices, position_coeffs, positions_out );
    if( du_out )
    {
        assert( dv_out );
        
        convert_vector_of_sparse_vectors_to_matrix( num_vertices, du_coeffs, *du_out );
        convert_vector_of_sparse_vectors_to_matrix( num_vertices, dv_coeffs, *dv_out );
    }
}

void compute_subdivision_coefficients_for_mesh(
    int num_vertices,
    OpenSubdiv::HbrMesh<OpenSubdiv::FarStencilFactoryVertex>* mesh,
    const std::vector< real_t >& us,
    const std::vector< real_t >& vs,
    SparseMatrix_t& positions_out,
    SparseMatrix_t* du_out,
    SparseMatrix_t* dv_out
    )
{
    // Intermediary output vectors.
    std::vector< sparse_vector_t > position_coeffs;
    std::vector< sparse_vector_t > du_coeffs;
    std::vector< sparse_vector_t > dv_coeffs;
    
    // Call the function we're wrapping.
    compute_subdivision_coefficients_for_mesh(
        num_vertices,
        mesh,
        us, vs,
        position_coeffs,
        du_out ? &du_coeffs : nullptr,
        dv_out ? &dv_coeffs : nullptr
        );
    
    // NOTE: We don't call mesh->GetNumVertices() because that may not be the number of coarse vertices.
    convert_vector_of_sparse_vectors_to_matrix( num_vertices, position_coeffs, positions_out );
    if( du_out )
    {
        assert( dv_out );
        
        convert_vector_of_sparse_vectors_to_matrix( num_vertices, du_coeffs, *du_out );
        convert_vector_of_sparse_vectors_to_matrix( num_vertices, dv_coeffs, *dv_out );
    }
}

}

#ifdef SUBDIVISION_MATRICES_EIGEN_MAIN
#include <iostream>

#define TORUS 32, { {4, 5, 1, 0}, {5, 6, 2, 1}, {6, 7, 3, 2}, {7, 4, 0, 3}, {8, 9, 5, 4}, {9, 10, 6, 5}, {10, 11, 7, 6}, {11, 8, 4, 7}, {12, 13, 9, 8}, {13, 14, 10, 9}, {14, 15, 11, 10}, {15, 12, 8, 11}, {16, 17, 13, 12}, {17, 18, 14, 13}, {18, 19, 15, 14}, {19, 16, 12, 15}, {20, 21, 17, 16}, {21, 22, 18, 17}, {22, 23, 19, 18}, {23, 20, 16, 19}, {24, 25, 21, 20}, {25, 26, 22, 21}, {26, 27, 23, 22}, {27, 24, 20, 23}, {28, 29, 25, 24}, {29, 30, 26, 25}, {30, 31, 27, 26}, {31, 28, 24, 27}, {0, 1, 29, 28}, {1, 2, 30, 29}, {2, 3, 31, 30}, {3, 0, 28, 31} }

template< typename T >
void print_dims( const std::string& prefix, const T& m )
{
    std::cerr << prefix << m.rows() << "x" << m.cols() << '\n';
}

#include <valarray>
void save_compute_coefficients( int resolution, std::ostream& out )
{
    using namespace subdivision_matrix;
    
    std::vector< real_t > us, vs;
    createUVs( resolution, resolution, us, vs );
    
    SparseMatrix_t position_coeffs;
    compute_subdivision_coefficients_for_mesh(
        // Cube
        // CUBE,
        // Torus
        TORUS,
        us, vs,
        position_coeffs
        );
    
    std::vector< std::vector< real_t > > positions_orig( { { 1.25052, 0.517982, 0.353553 }, { 0.597239, 0.247384, 0.353553 }, { 0.597239, 0.247384, -0.353553 }, { 1.25052, 0.517982, -0.353553 }, { 0.517982, 1.25052, 0.353553 }, { 0.247384, 0.597239, 0.353553 }, { 0.247384, 0.597239, -0.353553 }, { 0.517982, 1.25052, -0.353553 }, { -0.517982, 1.25052, 0.353553 }, { -0.247384, 0.597239, 0.353553 }, { -0.247384, 0.597239, -0.353553 }, { -0.517982, 1.25052, -0.353553 }, { -1.25052, 0.517982, 0.353553 }, { -0.597239, 0.247384, 0.353553 }, { -0.597239, 0.247384, -0.353553 }, { -1.25052, 0.517982, -0.353553 }, { -1.25052, -0.517982, 0.353553 }, { -0.597239, -0.247384, 0.353553 }, { -0.597239, -0.247384, -0.353553 }, { -1.25052, -0.517982, -0.353553 }, { -0.517982, -1.25052, 0.353553 }, { -0.247384, -0.597239, 0.353553 }, { -0.247384, -0.597239, -0.353553 }, { -0.517982, -1.25052, -0.353553 }, { 0.517982, -1.25052, 0.353553 }, { 0.247384, -0.597239, 0.353553 }, { 0.247384, -0.597239, -0.353553 }, { 0.517982, -1.25052, -0.353553 }, { 1.25052, -0.517982, 0.353553 }, { 0.597239, -0.247384, 0.353553 }, { 0.597239, -0.247384, -0.353553 }, { 1.25052, -0.517982, -0.353553 } } );
    Eigen::Matrix< real_t, 3, Eigen::Dynamic > positions;
    positions.resize( 3, positions_orig.size() );
    for( int row = 0; row < positions.rows(); ++row )
    for( int col = 0; col < positions.cols(); ++col )
    {
        positions( row, col ) = positions_orig.at(col).at(row);
    }
    
    print_dims( "position_coeffs is ", position_coeffs );
    print_dims( "positions is ", positions );
    
    Eigen::MatrixXf all_pos = positions * position_coeffs.transpose();
    print_dims( "all_pos is ", all_pos );
    // out << all_pos << '\n';
    
    for( int i = 0; i < position_coeffs.rows(); ++i )
    {
        // std::cerr << position_coeffs.row(i) << '\n';
        
        Eigen::VectorXf pos = positions * position_coeffs.row(i).transpose();
        out << "v " << pos(0) << ' ' << pos(1) << ' ' << pos(2) << '\n';
    }
}

#include <cstdlib>
#include <fstream>
void usage( const char* argv0 )
{
    std::cerr << "Usage: " << argv0 << " save_compute_coefficients resolution\n";
    exit( -1 );
}
int main( int argc, char* argv[] )
{
    using std::string;
    
    if( argc != 3 ) usage( argv[0] );
    
    const int resolution = atoi( argv[2] );
    if( resolution < 0 ) usage( argv[0] );
    
    const string test_name = string(argv[1]);
    
    if( string("save_compute_coefficients") == test_name )
    {
        std::ofstream out( "torus-eigen.obj" );
        save_compute_coefficients( resolution, out );
    }
    else
    {
        usage( argv[0] );
    }
    
    return 0;
}
#endif

// c++ -g subdivision_skinning_wrapper_test.cpp -o subdivision_skinning_wrapper_test ../../build/lib/libsubdivision_skinning.a ../../build/lib/libosdCPU.a
// c++ -g subdivision_skinning_wrapper_test.cpp -o subdivision_skinning_wrapper_test ../../build/Debug/lib/libsubdivision_skinning.a ../../build/Debug/lib/libosdCPU.a

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <cstdlib> // atoi
#include <cassert>
#include <cmath>

#include <algorithm> // std::fill

extern "C"
{
#include "subdivision_skinning_wrapper.h"
}

namespace
{

subdivision_evaluator_real_t TORUS_VS[] = { 1.25052, 0.517982, 0.353553, 0.597239, 0.247384, 0.353553, 0.597239, 0.247384, -0.353553, 1.25052, 0.517982, -0.353553, 0.517982, 1.25052, 0.353553, 0.247384, 0.597239, 0.353553, 0.247384, 0.597239, -0.353553, 0.517982, 1.25052, -0.353553, -0.517982, 1.25052, 0.353553, -0.247384, 0.597239, 0.353553, -0.247384, 0.597239, -0.353553, -0.517982, 1.25052, -0.353553, -1.25052, 0.517982, 0.353553, -0.597239, 0.247384, 0.353553, -0.597239, 0.247384, -0.353553, -1.25052, 0.517982, -0.353553, -1.25052, -0.517982, 0.353553, -0.597239, -0.247384, 0.353553, -0.597239, -0.247384, -0.353553, -1.25052, -0.517982, -0.353553, -0.517982, -1.25052, 0.353553, -0.247384, -0.597239, 0.353553, -0.247384, -0.597239, -0.353553, -0.517982, -1.25052, -0.353553, 0.517982, -1.25052, 0.353553, 0.247384, -0.597239, 0.353553, 0.247384, -0.597239, -0.353553, 0.517982, -1.25052, -0.353553, 1.25052, -0.517982, 0.353553, 0.597239, -0.247384, 0.353553, 0.597239, -0.247384, -0.353553, 1.25052, -0.517982, -0.353553 };
// #define TORUS_FACES { {4, 5, 1, 0}, {5, 6, 2, 1}, {6, 7, 3, 2}, {7, 4, 0, 3}, {8, 9, 5, 4}, {9, 10, 6, 5}, {10, 11, 7, 6}, {11, 8, 4, 7}, {12, 13, 9, 8}, {13, 14, 10, 9}, {14, 15, 11, 10}, {15, 12, 8, 11}, {16, 17, 13, 12}, {17, 18, 14, 13}, {18, 19, 15, 14}, {19, 16, 12, 15}, {20, 21, 17, 16}, {21, 22, 18, 17}, {22, 23, 19, 18}, {23, 20, 16, 19}, {24, 25, 21, 20}, {25, 26, 22, 21}, {26, 27, 23, 22}, {27, 24, 20, 23}, {28, 29, 25, 24}, {29, 30, 26, 25}, {30, 31, 27, 26}, {31, 28, 24, 27}, {0, 1, 29, 28}, {1, 2, 30, 29}, {2, 3, 31, 30}, {3, 0, 28, 31} }
int TORUS_FACES[] = { 4, 5, 1, 0, 5, 6, 2, 1, 6, 7, 3, 2, 7, 4, 0, 3, 8, 9, 5, 4, 9, 10, 6, 5, 10, 11, 7, 6, 11, 8, 4, 7, 12, 13, 9, 8, 13, 14, 10, 9, 14, 15, 11, 10, 15, 12, 8, 11, 16, 17, 13, 12, 17, 18, 14, 13, 18, 19, 15, 14, 19, 16, 12, 15, 20, 21, 17, 16, 21, 22, 18, 17, 22, 23, 19, 18, 23, 20, 16, 19, 24, 25, 21, 20, 25, 26, 22, 21, 26, 27, 23, 22, 27, 24, 20, 23, 28, 29, 25, 24, 29, 30, 26, 25, 30, 31, 27, 26, 31, 28, 24, 27, 0, 1, 29, 28, 1, 2, 30, 29, 2, 3, 31, 30, 3, 0, 28, 31 };

void print_OBJ( std::ostream& out, int num_vs, subdivision_evaluator_real_t* vs, int num_faces, int* quad_faces, const std::string& header_message = "" )
{
    // Save an optional header message.
    out << "# OBJ: " << header_message << '\n';
    
    // Save vertices.
    for( int i = 0; i < num_vs; ++i )
    {
        out << "v " << vs[3*i+0] << ' ' << vs[3*i+1] << ' ' << vs[3*i+2] << '\n';
    }
    
    out << '\n';
    
    for( int i = 0; i < num_faces; ++i )
    {
        // Vertices are 1-indexed.
        out << "f " << (1+quad_faces[ 4*i + 0 ]) << ' ' << (1+quad_faces[ 4*i + 1 ]) << ' ' << (1+quad_faces[ 4*i + 2 ]);
        if( -1 != quad_faces[ 4*i + 3 ] ) out << ' ' << (1+quad_faces[ 4*i + 3 ]) << '\n';
        else out << '\n';
    }
}
void save_OBJ( const std::string& out_path, int num_vs, subdivision_evaluator_real_t* vs, int num_faces, int* quad_faces, const std::string& header_message = "" )
{
    std::ofstream out( out_path );
    print_OBJ( out, num_vs, vs, num_faces, quad_faces, header_message );
    std::cout << "Saved a mesh to: " << out_path << '\n';
}

bool load_OBJ( std::istream& in, int* num_vs_out, subdivision_evaluator_real_t** vs_out, int* num_faces_out, int** quad_faces_out )
{
    std::string line;
    std::istringstream linestr;
    
    std::vector< subdivision_evaluator_real_t > vs;
    std::vector< int > faces;
    
    int line_number = 0;
    std::string type;
    while( !( in >> std::ws ).eof() )
    {
        line_number += 1;
        
        std::getline( in, line );
        if( line.empty() ) continue;
        
        linestr.clear();
        linestr.str( line );
        linestr >> type;
        
        if( type == std::string("v") )
        {
            subdivision_evaluator_real_t x(-31337), y(-31337), z(-31337);
            linestr >> x >> y >> z;
            vs.push_back( x );
            vs.push_back( y );
            vs.push_back( z );
            if( !linestr )
            {
                std::cerr << "Invalid vertex encountered on line " << line_number << ".\n";
                return false;
            }
        }
        else if( type == std::string("f") )
        {
            int num_face_vertices = 0;
            while( !( linestr >> std::ws ).eof() )
            {
                std::string vi;
                linestr >> vi;
                num_face_vertices += 1;
                
                // Now look for a slash, in case there's a bundle.
                int slash = vi.find( "/" );
                if( std::string::npos == slash )
                {
                    faces.push_back( atoi( vi.c_str() ) );
                }
                else
                {
                    faces.push_back( atoi( vi.substr( 0, slash ).c_str() ) );
                }
                // OBJ are 1-indexed, and a negative number indexes from the end.
                // Nothing should equal 0.
                if( faces.back() == 0 )
                {
                    std::cerr << "Invalid face vertex index of 0 encountered on line " << line_number << ".\n";
                    return false;
                }
            }
            
            if( !( num_face_vertices == 3 || num_face_vertices == 4 ) )
            {
                std::cerr << "Invalid face without three or four vertices encountered on line " << line_number << ".\n";
                return false;
            }
            // A triangle has a -1 as the fourth vertex index.
            // UPDATE: push_back a 0, because 
            if( 3 == num_face_vertices ) faces.push_back( 0 );
        }
    }
    
    assert( vs.size() % 3 == 0 );
    *num_vs_out = vs.size()/3;
    *vs_out = new subdivision_evaluator_real_t[ vs.size() ];
    std::copy( vs.begin(), vs.end(), *vs_out );
    
    assert( faces.size() % 4 == 0 );
    *num_faces_out = faces.size()/4;
    *quad_faces_out = new int[ faces.size() ];
    for( int i = 0; i < faces.size(); ++i )
    {
        (*quad_faces_out)[i] =
            faces[i] >= 0
            ? faces[i] - 1
            // Negative indices add to the back.
            : vs.size() + faces[i]
            ;
    }
    
    return true;
}
bool load_OBJ( const std::string& in_path, int* num_vs_out, subdivision_evaluator_real_t** vs_out, int* num_faces_out, int** quad_faces_out )
{
    std::ifstream in( in_path );
    return load_OBJ( in, num_vs_out, vs_out, num_faces_out, quad_faces_out );
}

}


int main( int argc, const char* argv[] )
{
    subdivision_evaluator_real_t* vs = &TORUS_VS[0];
    int* faces = &TORUS_FACES[0];
    int num_vs = sizeof( TORUS_VS )/sizeof( subdivision_evaluator_real_t )/3;
    int num_faces = sizeof( TORUS_FACES )/sizeof( int )/4;
    
    if( argc == 2 )
    {
        const bool success = load_OBJ( argv[1], &num_vs, &vs, &num_faces, &faces );
        if( !success ) return -1;
    }
    
    std::cout << "Number of vertices: " << num_vs << '\n';
    std::cout << "Number of faces: " << num_faces << '\n';
    
    save_OBJ( "original.obj", num_vs, vs, num_faces, faces );
    
    void* eval = new_subdivision_evaluator( num_vs, vs, num_faces, faces, 4 );
    
    const int num_refined_faces = num_refined_quad_faces_of_subdivision_evaluator( eval );
    int* refined_faces = new int[ num_refined_faces*4 ];
    get_refined_quad_faces_of_subdivision_evaluator( eval, refined_faces );
    
    const int num_refined_vs = num_refined_vertices_of_subdivision_evaluator( eval );
    subdivision_evaluator_real_t* refined_vs = new subdivision_evaluator_real_t[ num_refined_vs*3 ];
    get_refined_vertices_of_subdivision_evaluator( eval, refined_vs );
    
    save_OBJ( "original-refined.obj", num_refined_vs, refined_vs, num_refined_faces, refined_faces );
    
    // Test the non-3 dimensional path by duplicating submitting 6-dimensional
    // vertices xyzxyz:
    {
        std::vector< subdivision_evaluator_real_t > control_vs6( num_vs*6 );
        for( int i = 0; i < num_vs; ++i )
        for( int c = 0; c < 3; ++c )
        {
            control_vs6[ 6*i + 3 + c ] = control_vs6[ 6*i + c ] = vs[ 3*i + c ];
        }
        std::vector< subdivision_evaluator_real_t > refined_vs6( num_refined_vs*6 );
        get_refined_vertices_of_subdivision_evaluator_with_control_vertices( eval, 6, &control_vs6[0], &refined_vs6[0] );
        // Verify the output:
        subdivision_evaluator_real_t total_diff = 0.;
        for( int i = 0; i < num_refined_vs; ++i )
        for( int c = 0; c < 3; ++c )
        {
            total_diff += fabs( refined_vs6[ 6*i + c ] - refined_vs6[ 6*i + 3 + c ] );
            total_diff += fabs( refined_vs6[ 6*i + c ] - refined_vs[ 3*i + c ] );
        }
        std::cout << "Total difference of 6-dimensional versus 3-dimensional (should be 0): " << total_diff << '\n';
    }
    
    subdivision_evaluator_real_t* weights = new subdivision_evaluator_real_t[ num_refined_vs ];
    std::fill( weights, weights + num_refined_vs, 1. );
    void* engine = new_subdivision_skinning_engine( eval, 1, weights );
    
    subdivision_evaluator_real_t* vs_modified = new subdivision_evaluator_real_t[ num_vs*3 ];
    /// This produces identical output as input:
    // const subdivision_evaluator_real_t transforms[] = { 1,0,0,0, 0,1,0,0, 0,0,1,0, 0,0,0,1 };
    /// This produces output whose coordinates are scaled by 2:
    const subdivision_evaluator_real_t transforms[] = { 2,0,0,0, 0,2,0,0, 0,0,2,0, 0,0,0,2 };
    compute_control_mesh_vertices_given_transforms_for_subdivision_skinning_engine( engine, (subdivision_evaluator_real_t*)transforms, vs_modified );
    save_OBJ( "modified.obj", num_vs, vs_modified, num_faces, faces );
    
    const int num_refined_vs2 = num_refined_vertices_of_subdivision_skinning_engine( engine );
    std::cout << "num_refined_vs eval: " << num_refined_vs << '\n';
    std::cout << "num_refined_vs engine: " << num_refined_vs2 << '\n';
    
    get_refined_vertices_of_subdivision_skinning_engine_with_control_vertices( engine, 3, vs_modified, refined_vs );
    save_OBJ( "modified-refined.obj", num_refined_vs, refined_vs, num_refined_faces, refined_faces );
    
    // Test the non-3 dimensional path by duplicating submitting 6-dimensional
    // vertices xyzxyz:
    {
        std::vector< subdivision_evaluator_real_t > control_vs6( num_vs*6 );
        for( int i = 0; i < num_vs; ++i )
        for( int c = 0; c < 3; ++c )
        {
            control_vs6[ 6*i + 3 + c ] = control_vs6[ 6*i + c ] = vs_modified[ 3*i + c ];
        }
        std::vector< subdivision_evaluator_real_t > refined_vs6( num_refined_vs*6 );
        get_refined_vertices_of_subdivision_skinning_engine_with_control_vertices( engine, 6, &control_vs6[0], &refined_vs6[0] );
        // Verify the output:
        subdivision_evaluator_real_t total_diff = 0.;
        for( int i = 0; i < num_refined_vs; ++i )
        for( int c = 0; c < 3; ++c )
        {
            total_diff += fabs( refined_vs6[ 6*i + c ] - refined_vs6[ 6*i + 3 + c ] );
            total_diff += fabs( refined_vs6[ 6*i + c ] - refined_vs[ 3*i + c ] );
        }
        std::cout << "Total difference of 6-dimensional versus 3-dimensional (should be 0): " << total_diff << '\n';
    }
    
    delete_subdivision_skinning_engine( engine );
    delete_subdivision_evaluator( eval );
    
    delete [] vs_modified;
    delete [] refined_faces;
    delete [] refined_vs;
    delete [] weights;
    
    if( argc == 2 )
    {
        delete[] vs;
        delete[] faces;
    }
    
    return 0;
}

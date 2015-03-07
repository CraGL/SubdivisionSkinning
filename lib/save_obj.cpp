#include "save_obj.h"

#include <iostream>
#include <fstream>

namespace save_obj
{
void save_mesh( const std::string& out_path, const std::vector< std::vector< int > >& faces, const std::vector< std::vector< real_t > >& vertices, const std::string& header_message )
{
	std::ofstream out( out_path );
    save_mesh( out, faces, vertices, header_message );
    std::cout << "Saved a mesh to: " << out_path << '\n';
}
void save_mesh( std::ostream& out, const std::vector< std::vector< int > >& faces, const std::vector< std::vector< real_t > >& vertices, const std::string& header_message )
{
    // Save an optional header message.
    out << header_message;
    
    // Save vertices.
    for( const auto& vert: vertices )
    {
        out << "v";
        for( const auto& coord: vert )
        {
            out << ' ' << coord;
        }
        out << '\n';
    }
    
    out << '\n';
    
    for( const auto& f: faces )
    {
        out << 'f';
        for( const auto& vi : f )
        {
            // Vertices are 1-indexed.
            out << ' ' << (vi+1);
        }
        out << '\n';
    }
}
}

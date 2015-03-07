#ifndef __save_obj_h__
#define __save_obj_h__

#include <vector>
#include <string>
#include <iosfwd>

namespace save_obj
{
    typedef float real_t;
    void save_mesh( const std::string& outpath, const std::vector< std::vector< int > >& faces, const std::vector< std::vector< float > >& vertices, const std::string& header_message = "" );
    void save_mesh( std::ostream& out, const std::vector< std::vector< int > >& faces, const std::vector< std::vector< float > >& vertices, const std::string& header_message = "" );
}

#endif /* __save_obj_h__ */

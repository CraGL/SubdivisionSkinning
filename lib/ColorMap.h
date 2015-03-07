#ifndef __ColorMap_h__
#define __ColorMap_h__

template < typename value_t > struct ColorMap
{   
    typedef float color_component_t;
    struct Color
    {
        color_component_t r;
        color_component_t g;
        color_component_t b;
        
        Color() : r(0), g(0), b(0) {}
        Color( color_component_t r0, color_component_t g0, color_component_t b0 ) : r(r0), g(g0), b(b0) {}
        
        static Color lerp( const Color& a, const Color& b, value_t t )
        {
            return Color(
                a.r + t*(b.r-a.r),
                a.g + t*(b.g-a.g),
                a.b + t*(b.b-a.b)
                );
        }
    };
    
    // The minimum value.
    value_t min_value;
    // The maximum value.
    value_t max_value;
    // The 0-th entry corresponds to values <= min_value.
    // The last entry corresponds to values >= max_value.
    std::vector< Color > colormap;
    
	ColorMap( value_t min, value_t max )
        : min_value( min ), max_value( max )
    {
        SetWhiteToRed();
    }
    
   	ColorMap()
        : min_value( 0 ), max_value( 1 )
    {
        SetWhiteToRed();
    }
    
    void SetWhiteToRed()
    {
        colormap.resize(2);
        colormap[0] = Color( 1,1,1 );
        colormap[1] = Color( 1,0,0 );
    }
    
    // Returns the computed color.
    Color Value2Color( const value_t val )
    {
        assert( colormap.size() >= 2 );
        assert( max_value >= min_value );
        
        /// 1 Find 'val' as a fraction from min_value to max_value.
        /// 2 Find the percent's index into colormap.
        /// 3 Find the percent between the two surrounding colors.
        /// 4 Lerp the two surrounding colors in the color map.
        
        static const double eps = 1e-7;
        
        if( fabs( max_value - min_value ) < eps )
        {
            return Color::lerp( colormap.at(0), colormap.at(1), .5 );
        }
        
        /// 1
        const value_t fraction = std::max( value_t(0), std::min( value_t(1), ( val - min_value )/(max_value - min_value) ) );
        /// 2
        const value_t real_index = fraction*( colormap.size()-1 );
        const int lower_index = std::max( 0, std::min( int(colormap.size())-2, int( real_index ) ) );
        /// 3
        value_t remaining_fraction = real_index - lower_index;
        assert( remaining_fraction > -eps );
        assert( remaining_fraction < 1+eps );
        remaining_fraction = std::max( value_t(0), std::min( value_t(1), remaining_fraction ) );
        /// 4
        return Color::lerp( colormap.at( lower_index ), colormap.at( lower_index+1 ), remaining_fraction );
    }
    // Directly puts the output into a pointer to three color_component_t's.
    void Value2Color( const value_t val, color_component_t* rgb_out )
    {
        const Color result = Value2Color( val );
        rgb_out[0] = result.r;
        rgb_out[1] = result.g;
        rgb_out[2] = result.b;
    }
};

#endif /* __ColorMap_h__ */

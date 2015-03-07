#ifndef __Timing_hpp__
#define __Timing_hpp__

void tic( const char* label = 0 ) ;
void toc( const char* label = 0 ) ;

struct Tick
{
    Tick( const char* alabel = 0 ) : label( alabel ) { tic( label ); }
    ~Tick() { toc( label ); }
    
    const char* label;
};

#endif /* __Timing_hpp__ */

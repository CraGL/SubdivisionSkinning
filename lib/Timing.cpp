#include "Timing.hpp"

#include <iostream>
#include <sys/time.h> // gettimeofday(), getrusage()
#include <sys/resource.h> // getrusage()

#define TICTOC 1

static double timeofday()
{
	timeval tp ;
	gettimeofday( &tp, NULL ) ;
	double result = (double)tp.tv_sec + 1e-6 * (double)tp.tv_usec ;
	return result ;
}
static double cpuusage()
{
	struct rusage rus ;
	getrusage( RUSAGE_SELF, &rus ) ;
	double user = (double)rus.ru_utime.tv_sec + 1e-6 * (double)rus.ru_utime.tv_usec ;
	double sys = (double)rus.ru_stime.tv_sec + 1e-6 * (double)rus.ru_stime.tv_usec ;
	
	return user + sys ;
}

static const int TICTOC_STACK_SIZE = 25 ;
static double sTicStart[ TICTOC_STACK_SIZE ] = { 0.0 } ;
static int sTicStartIndex = -1 ;
void tic( const char* txt )
{
#if TICTOC
	sTicStartIndex += 1 ;
	if( sTicStartIndex >= TICTOC_STACK_SIZE )
	{
		sTicStartIndex = TICTOC_STACK_SIZE - 1 ;
		std::cerr << "tic()/toc() max depth reached.  Re-using top of stack.\n";
	}
	
	if( txt ) std::cout << '[' << txt << "] begins\n";
	
	// sTicStart[ sTicStartIndex ] = timeofday() ;
	sTicStart[ sTicStartIndex ] = cpuusage() ;
#endif
}
void toc( const char* txt )
{
#if TICTOC
	if( sTicStartIndex >= 0 )
	{
		// double curtime = timeofday() ;
		double curtime = cpuusage() ;
		if( txt ) std::cout << '[' << txt << "] ";
		std::cout << "elapsed cpu time: " << (curtime - sTicStart[ sTicStartIndex ]) << std::endl ;
		sTicStartIndex -= 1 ;
	}
	else
	{
		std::cerr << "Called toc() with no matching tic().\n";
	}
#endif
}

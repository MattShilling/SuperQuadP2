#define _USE_MATH_DEFINES

#include <cmath>
#include <cstdlib>
#include <ctime>
#include <random>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <omp.h>

#include "test_rig.h"

#define NUMNODES 256
#define NUMT 8
#define N 4

double CalcSpeedup(double test_two, double test_one) {
  return test_two / test_one;
}

float CalcFp(double speedup, double threads) {
  return (threads / (threads - 1.0)) * (1. - (1. / speedup));
}

#define XMIN     -1.
#define XMAX      1.
#define YMIN     -1.
#define YMAX      1.

float Height( int, int );

int main( int argc, char *argv[ ] )
{

  omp_set_num_threads(NUMT);
	// the area of a single full-sized tile:

	float fullTileArea = (  ( ( XMAX - XMIN )/(float)(NUMNODES-1) )  *
				( ( YMAX - YMIN )/(float)(NUMNODES-1) )  );

	// sum up the weighted heights into the variable "volume"
	// using an OpenMP for loop and a reduction:

  float total_height = 0;

#pragma omp parallel for collapse(2), default(none), reduction(+ : total_height)
  for( int iv = 0; iv < NUMNODES; iv++ )
  {
    for( int iu = 0; iu < NUMNODES; iu++ )
    {
      float factor = 0.0;
      const bool side_iv = (iv == NUMNODES - 1) || (iv == 0);
      const bool side_iu = (iu == NUMNODES - 1) || (iu == 0);
      // XOR
      if  (!side_iv != !side_iu) {
        // We are on a side so only add 0.5
        factor = 0.5;
      } else if (side_iv && side_iu) {
        // We are on a corner so only add 0.25
        factor = 0.25;
      } else {
        factor = 1.0;
      }

      total_height += factor * Height(iu, iv);

    }
  }

  float total_volume = total_height * fullTileArea * 2.0;

  std::cout << total_volume << std::endl;

}


// iu,iv = 0 .. NUMNODES-1
float Height( int iu, int iv )
{
	float x = -1.  +  2.*(float)iu /(float)(NUMNODES-1);	// -1. to +1.
	float y = -1.  +  2.*(float)iv /(float)(NUMNODES-1);	// -1. to +1.

	float xn = pow( fabs(x), (double)N );
	float yn = pow( fabs(y), (double)N );
	float r = 1. - xn - yn;
	if( r < 0. )
	        return 0.;
	float height = pow( 1. - xn - yn, 1./(float)N );
	return height;
}
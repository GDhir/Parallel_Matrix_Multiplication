#include "lib_mylib_atb.h"
#include <stdio.h>

void atb_par(const double *__restrict__ A, const double *__restrict__ B, double *__restrict__ C, int Ni, int Nj, int Nk)
{
  int i, j, k;

  if( Ni == 1 && Nj == 1 && Nk == 1 ) {
     C[0] = A[0]*B[0];
  }
  else if( Ni <= 64 && Nj <= 64 && Nk >= 2048 ) {

      atb_par_histogram( A, B, C, Ni, Nj, Nk );

  }
  else if( Ni <= 16 && Nj >= 64*Ni ) {

    if( Nk <= 2 ) {

        // j tiling without unrolling
        atb_par_jtiled( A, B, C, Ni, Nj, Nk );

    }
    else {

        // j tiling with unrolling
        atb_par_jtiled_kunroll2( A, B, C, Ni, Nj, Nk );

    }

  }
  else {

    if( Ni <= 4 ) {

      if( Nk <= 2 ) {
        // no tiling no unrolling
        atb_par_ref( A, B, C, Ni, Nj, Nk );
      }
      else {
        // no tiling, k unrolled by 4
        atb_par_kunroll4( A, B, C, Ni, Nj, Nk );

      }


    }
    else if( Nk < 256 ) {
      // i tiling

      if( Nk <= 2 ) {

        // i tiling without unrolling
        atb_par_itiled( A, B, C, Ni, Nj, Nk );

      }
      else {

        // i tiling with unrolling
        atb_par_itiled_kunroll2( A, B, C, Ni, Nj, Nk );

      }
    }
    else {
      // ik tiling 
        // ik tiling with unrolling
        // printf("abcd");
        atb_par_iktiled_kunroll2( A, B, C, Ni, Nj, Nk );

    }
  }
}


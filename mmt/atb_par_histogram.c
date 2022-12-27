#include "stdio.h"
#include "stdlib.h"

void atb_par_histogram(double *__restrict__ a, double *__restrict__ b, double *__restrict__ c, int Ni, int Nj, int Nk)
{
   int i, j, k;

#pragma omp parallel private(i, j, k)
   {
      // int remk = Nk%6;
      double cval[Ni*Nj];

      for(i = 0; i < Ni*Nj; i++ ) {
         cval[i] = 0;
      }

      #pragma omp for schedule(static)
      for(k = 0; k < Nk; k++ ) {
         for (i = 0; i < Ni; i++) {
            for( j = 0; j < Nj; j++ ) {
               cval[ i*Nj + j ] += a[ k*Ni + i ]*b[ k*Nj + j ];
            }
         }
      }

      // #pragma omp critical
      {
         for (i = 0; i < Ni; i++) {
            for( j = 0; j < Nj; j++ ) {
               
               // printf( "%f, \t %f \n", c[ i*Nj + j ], cval[ i*Nj + j ] );
               #pragma omp atomic
               c[ i*Nj + j ] += cval[ i*Nj + j ];
               // printf( "%f, \t %f \n", c[ i*Nj + j ], cval[ i*Nj + j ] );

            }
         }
      }
   }
}

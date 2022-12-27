#include "stdio.h"
#include "stdlib.h"

void mmtt_par(int Ni, int Nj, int Nk, double *__restrict__ a, double *__restrict__ b, double *__restrict__ c)
{
   int i, j, k;

#pragma omp parallel private(i, j, k)
   {
      double cval[Ni*Nj];

      for(i = 0; i < Ni*Nj; i++ ) {
         cval[i] = 0;
      }

      #pragma omp for 
      for(k = 0; k < Nk; k++ ) {
         for (i = 0; i < Ni; i++) {
            for( j = 0; j < Nj; j++ ) {
               cval[ i*Nj + j ] += a[ k*Ni + i ]*b[ j*Nk + k ];
            }
         }
      }

      for (i = 0; i < Ni; i++) {
         for( j = 0; j < Nj; j++ ) {
            
            #pragma omp atomic
            c[ i*Nj + j ] += cval[ i*Nj + j ];

         }
      }
   }
}

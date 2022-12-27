#include "utils.h"

void atb_par_jtiled_kunroll2( const double *__restrict__ a, const double *__restrict__ b, double *__restrict__ c, int ni, int nj, int nk)
{
   int i, j, k, ii, jj, kk;
   int bs = 4;
   int remk = nk%2;

#pragma omp parallel private(i, j, k, jj)
   {
    #pragma omp for
    for (jj = 0; jj < nj; jj += bs) {
            for (i = 0; i < ni; i++) {

                for( k = 0; k < remk; k++ ) {
                    for( j = jj; j < min(nj, jj + bs); j++ ) {
                        c[ i*nj + j ] += a[ k*ni + i ]*b[ k*nj + j ];
                    }
                }

                for(k = remk; k < nk; k += 2 ) {
                    for( j = jj; j < min(nj, jj + bs); j++ ) {
                        c[ i*nj + j ] += a[ k*ni + i ]*b[ k*nj + j ];
                        c[ i*nj + j ] += a[ (k + 1)*ni + i ]*b[ (k + 1)*nj + j ];
                    }
                }

            }
        }
    }
}

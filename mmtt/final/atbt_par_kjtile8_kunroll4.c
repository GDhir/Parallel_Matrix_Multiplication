#include "utils.h"

void atbt_par_kjtile8_kunroll4(const double *__restrict__ a, const double *__restrict__ b, double *__restrict__ c, int ni, int nj, int nk)
{
   int i, j, k, ii, jj, kk;
   int bs = 8;
   int remk;

#pragma omp parallel private(i, j, k)
   {
    #pragma omp for
    for( jj = 0; jj < nj; jj += bs ) {
        for(kk = 0; kk < nk; kk += bs ) {
            for (i = 0; i < ni; i++) {

                    remk = ( min( nk, kk + bs ) - kk )%4;
                
                    for( k = kk; k < kk + remk; k++ ) {
                        for( j = jj; j < min( nj, jj + bs ); j++ ) {
                            c[ i*nj + j ] += a[ k*ni + i ]*b[ j*nk + k ];
                        }
                    }

                    for( k = kk + remk; k < min( nk, kk + bs ); k += 4 ) {
                        for( j = jj; j < min( nj, jj + bs ); j++ ) {
                            c[ i*nj + j ] += a[ k*ni + i ]*b[ j*nk + k ];
                            c[ i*nj + j ] += a[ (k + 1)*ni + i ]*b[ j*nk + k + 1 ];
                            c[ i*nj + j ] += a[ (k + 2)*ni + i ]*b[ j*nk + k + 2 ];
                            c[ i*nj + j ] += a[ (k + 3)*ni + i ]*b[ j*nk + k + 3 ];
                        }
                    }
                }
            }
        }
    }
}

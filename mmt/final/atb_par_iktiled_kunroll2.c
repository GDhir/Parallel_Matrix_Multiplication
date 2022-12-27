#include "utils.h"

void atb_par_iktiled_kunroll2(const double *__restrict__ a, const double *__restrict__ b, double *__restrict__ c, int ni, int nj, int nk)
{
   int i, j, k, ii, jj, kk;
   int bs = 4;
   int remk;

#pragma omp parallel private(i, j, k, ii, kk, remk)
   {
    #pragma omp for
    for (ii = 0; ii < ni; ii += bs) {
        for(kk = 0; kk < nk; kk += bs ) {
            for (i = ii; i < min( ni, ii + bs ); i++) {

                remk = (min( nk, kk + bs ) - kk)%2;
                for( k = kk; k < kk + remk; k++ ) {
                    for( j = 0; j < nj; j++ ) {
                        c[ i*nj + j ] += a[ k*ni + i ]*b[ k*nj + j ];
                    }
                }

                for(k = kk + remk; k < min( nk, kk + bs ); k+= 2 ) {
                    for( j = 0; j < nj; j++ ) {
                        c[ i*nj + j ] += a[ k*ni + i ]*b[ k*nj + j ];
                        c[ i*nj + j ] += a[ (k + 1)*ni + i ]*b[ (k + 1)*nj + j ];
                    }
                }
            }
        }
    }
    }
}

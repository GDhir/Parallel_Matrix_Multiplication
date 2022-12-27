#include "utils.h"

void atbt_par_ikjtile16_kunroll8(const double *__restrict__ A, const double *__restrict__ B, double *__restrict__ C, int ni, int nj, int nk)
{
   int i, j, k, ii, jj, kk;
   int bs = 16;
   int remk;

#pragma omp parallel private(i, j, k, ii, jj, kk, remk)
   {
    #pragma omp for
    for (ii = 0; ii < ni; ii += bs) {
        for(kk = 0; kk < nk; kk += bs ) {
            for( jj = 0; jj < nj; jj += bs ) {
                for (i = ii; i < min( ni, ii + bs ); i++) {

                    remk = (min( nk, kk + bs ) - kk)%8;
                    for( k = kk; k < kk + remk; k++ ) {
                        for( j = jj; j < min( nj, jj + bs ); j++ ) {
                            C[ i*nj + j ] += A[ k*ni + i ]*B[ j*nk + k ];
                        }
                    }

                    for(k = kk + remk; k < min( nk, kk + bs ); k+= 8 ) {
                        for( j = jj; j < min( nj, jj + bs ); j++ ) {
                            C[ i*nj + j ] += A[ k*ni + i ]*B[ j*nk + k ];
                            C[ i*nj + j ] += A[ (k + 1)*ni + i ]*B[ j*nk + k + 1 ];
                            C[ i*nj + j ] += A[ (k + 2)*ni + i ]*B[ j*nk + k + 2 ];
                            C[ i*nj + j ] += A[ (k + 3)*ni + i ]*B[ j*nk + k + 3 ];
                            C[ i*nj + j ] += A[ (k + 4)*ni + i ]*B[ j*nk + k + 4 ];
                            C[ i*nj + j ] += A[ (k + 5)*ni + i ]*B[ j*nk + k + 5 ];
                            C[ i*nj + j ] += A[ (k + 6)*ni + i ]*B[ j*nk + k + 6 ];
                            C[ i*nj + j ] += A[ (k + 7)*ni + i ]*B[ j*nk + k + 7 ];
                        }
                    }
                }
            }
        }
    }
}
}

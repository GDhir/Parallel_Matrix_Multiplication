#include"utils.h"

void atbt_par_ijtile16_kunroll4_kij(const double *__restrict__ A, const double *__restrict__ B, double *__restrict__ C, int Ni, int Nj, int Nk)
{
   int i, j, k, ii, jj, kk;
   int bs = 16;
   int remk = Nk%4;

#pragma omp parallel private(i, j, k, ii, jj)
   {
    #pragma omp for
    for (ii = 0; ii < Ni; ii += bs) {
        for( jj = 0; jj < Nj; jj += bs ) {
            for( k = 0; k < remk; k++ ) {
                for (i = ii; i < min( Ni, ii + bs ); i++) {
                    for( j = jj; j < min( Nj, jj + bs ); j++ ) {
                            C[ i*Nj + j ] += A[ k*Ni + i ]*B[ j*Nk + k ];
                    } 
                }
            }


            for( k = remk; k < Nk; k += 4 ) {
                for (i = ii; i < min( Ni, ii + bs ); i++) {
                    for( j = jj; j < min( Nj, jj + bs ); j++ ) {
                            C[ i*Nj + j ] += A[ k*Ni + i ]*B[ j*Nk + k ];
                            C[ i*Nj + j ] += A[ ( k + 1 )*Ni + i ]*B[ j*Nk + k + 1 ];
                            C[ i*Nj + j ] += A[ ( k + 2 )*Ni + i ]*B[ j*Nk + k + 2 ];
                            C[ i*Nj + j ] += A[ ( k + 3 )*Ni + i ]*B[ j*Nk + k + 3 ];
                    } 
                }
            }
        }
    }
}
}

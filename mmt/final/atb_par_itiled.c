#include "utils.h"

void atb_par_itiled( const double *__restrict__ a, const double *__restrict__ b, double *__restrict__ c, int ni, int nj, int nk)
{
   int i, j, k, ii, jj, kk;
   int bs = 4;

#pragma omp parallel private(i, j, k, ii)
   {
    #pragma omp for
    for (ii = 0; ii < ni; ii += bs) {
        for (i = ii; i < min( ni, ii + bs ); i++) {
            for(k = 0; k < nk; k++ ) {
                for( j = 0; j < nj; j++ ) {
                    c[ i*nj + j ] += a[ k*ni + i ]*b[ k*nj + j ];
                }
            }
        }
    }
    }
}

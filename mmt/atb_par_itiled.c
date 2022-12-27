int min( int a, int b ) {
    if (a > b)
        return b;
    else 
        return a;
}

void atb_par_itiled( double *__restrict__ a, double *__restrict__ b, double *__restrict__ c, int ni, int nj, int nk)
{
   int i, j, k, ii, jj, kk;
   int bs = 4;

#pragma omp parallel private(i, j, k)
   {
    #pragma omp for
    for (ii = 0; ii < ni; ii += bs) {
        // for(kk = 0; kk < nk; kk += bs ) {
            // for( jj = 0; jj < nj; jj += bs ) {
                for (i = ii; i < min( ni, ii + bs ); i++) {

                    // remk = (min( nk, kk + bs ) - kk)%2;
                    // for( k = kk; k < kk + remk; k++ ) {
                    //     for( j = jj; j < min( nj, jj + bs ); j++ ) {
                    //         c[ i*nj + j ] += a[ k*ni + i ]*b[ k*nj + j ];
                    //     }
                    // }

                    // for(k = kk + remk; k < min( nk, kk + bs ); k+= 2 ) {
                    // //    for( j = jj; j < min( nj, jj + bs ); j++ ) {
                    //     for( j = 0; j < nj; j++ ) {
                    //         c[ i*nj + j ] += a[ k*ni + i ]*b[ k*nj + j ];
                    //         c[ i*nj + j ] += a[ (k + 1)*ni + i ]*b[ (k + 1)*nj + j ];
                    //     }
                    // }
                // }

                    for(k = 0; k < nk; k++ ) {
                        for( j = 0; j < nj; j++ ) {
                            c[ i*nj + j ] += a[ k*ni + i ]*b[ k*nj + j ];
                        }
                    }


            }
        // }
    }
    }
}

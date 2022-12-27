int min( int a, int b ) {
    if (a > b)
        return b;
    else 
        return a;
}

void mmtt_par(int ni, int nj, int nk, double *__restrict__ a, double *__restrict__ b, double *__restrict__ c)
{
   int i, j, k, ii, jj, kk;
   int bs = 16;
   int remk = nk%8;

#pragma omp parallel private(i, j, k)
   {
    #pragma omp for
    for (jj = 0; jj < nj; jj += bs) {
            for( k = 0; k < remk; k++ ) {
                for( j = jj; j < min( nj, jj + bs ); j++ ) {
                    for (i = 0; i < ni; i++) {
                        c[ i*nj + j ] += a[ k*ni + i ]*b[ j*nk + k ];
                    }
                }
            }

            for( k = remk; k < nk; k += 8 ) {
                for( j = jj; j < min( nj, jj + bs ); j++ ) {

                    for (i = 0; i < ni; i++) {
                        c[ i*nj + j ] += a[ k*ni + i ]*b[ j*nk + k ];
                        c[ i*nj + j ] += a[ ( k + 1 )*ni + i ]*b[ j*nk + k + 1 ];
                        c[ i*nj + j ] += a[ ( k + 2 )*ni + i ]*b[ j*nk + k + 2 ];
                        c[ i*nj + j ] += a[ ( k + 3 )*ni + i ]*b[ j*nk + k + 3 ];
                        c[ i*nj + j ] += a[ ( k + 4 )*ni + i ]*b[ j*nk + k + 4 ];
                        c[ i*nj + j ] += a[ ( k + 5 )*ni + i ]*b[ j*nk + k + 5 ];
                        c[ i*nj + j ] += a[ ( k + 6 )*ni + i ]*b[ j*nk + k + 6 ];
                        c[ i*nj + j ] += a[ ( k + 7 )*ni + i ]*b[ j*nk + k + 7 ];
                    }
                }
            }
        }
   }
}

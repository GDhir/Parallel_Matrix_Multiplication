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
   int remk = nk%4;

#pragma omp parallel private(i, j, k)
   {
    #pragma omp for
    for (i = 0; i < ni; i++) {
        for( k = 0; k < remk; k++ ) {
            for (jj = 0; jj < nj; jj += bs) {
                for( j = jj; j < min( nj, jj + bs ); j++ ) {
                    c[ i*nj + j ] += a[ k*ni + i ]*b[ j*nk + k ];
                }
            }
        }

        for( k = remk; k < nk; k += 4 ) {
            for (jj = 0; jj < nj; jj += bs) {
                for( j = jj; j < min( nj, jj + bs ); j++ ) {
                    c[ i*nj + j ] += a[ k*ni + i ]*b[ j*nk + k ];
                    c[ i*nj + j ] += a[ ( k + 1 )*ni + i ]*b[ j*nk + k + 1 ];
                    c[ i*nj + j ] += a[ ( k + 2 )*ni + i ]*b[ j*nk + k + 2 ];
                    c[ i*nj + j ] += a[ ( k + 3 )*ni + i ]*b[ j*nk + k + 3 ];
                }
            }
        }
   }
}
}

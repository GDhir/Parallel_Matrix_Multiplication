int min( int a, int b ) {
    if (a > b)
        return b;
    else 
        return a;
}

void mmtt_par(int ni, int nj, int nk, double *__restrict__ a, double *__restrict__ b, double *__restrict__ c)
{
   int i, j, k, ii, jj, kk;
   int bs = 8;
   int remk = nk%2;

#pragma omp parallel private(i, j, k, ii, jj)
   {
    #pragma omp for
    for (ii = 0; ii < ni; ii += bs) {
        for( jj = 0; jj < nj; jj += bs ) {
            for (i = ii; i < min( ni, ii + bs ); i++) {
                for( k = 0; k < nk; k++ ) {
                    for( j = jj; j < min( nj, jj + bs ); j++ ) {
                            c[ i*nj + j ] += a[ k*ni + i ]*b[ j*nk + k ];
                    } 
                }
            }
        }
    }
}
}

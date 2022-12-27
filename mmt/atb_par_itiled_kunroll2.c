int min( int a, int b ) {
    if (a > b)
        return b;
    else 
        return a;
}

void atb_par_itiled_kunroll2(int ni, int nj, int nk, double *__restrict__ a, double *__restrict__ b, double *__restrict__ c)
{
   int i, j, k, ii, jj, kk;
   int bs = 4;
   int remk = nk%2;

#pragma omp parallel private(i, j, k)
   {
    #pragma omp for
    for (ii = 0; ii < ni; ii += bs) {
        for (i = ii; i < min( ni, ii + bs ); i++) {
            for(k = 0; k < remk; k++ ) {
                for( j = 0; j < nj; j++ ) {
                    c[ i*nj + j ] += a[ k*ni + i ]*b[ k*nj + j ];
                }
            }

            for( k = remk; k < nk; k += 2 ) {
                for( j = 0; j < nj; j++ ) {
                    c[ i*nj + j ] += a[ k*ni + i ]*b[ k*nj + j ];
                    c[ i*nj + j ] += a[ (k + 1)*ni + i ]*b[ (k + 1)*nj + j ];
                }
            }
        }
    }
}
}

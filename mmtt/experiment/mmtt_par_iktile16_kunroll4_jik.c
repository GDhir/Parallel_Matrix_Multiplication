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
   int remk = nk%2;

#pragma omp parallel private(i, j, k)
   {
    #pragma omp for
    for (ii = 0; ii < ni; ii += bs) {
        for(kk = 0; kk < nk; kk += bs ) {
            for( j = 0; j < nj; j++ ) {
                for (i = ii; i < min( ni, ii + bs ); i++) {

                    remk = (min( nk, kk + bs ) - kk)%2;
                    for( k = kk; k < kk + remk; k++ ) {
                            c[ i*nj + j ] += a[ k*ni + i ]*b[ j*nk + k ];
                    }

                    for(k = kk + remk; k < min( nk, kk + bs ); k+= 2 ) {
                    //    for( j = jj; j < min( nj, jj + bs ); j++ ) {
                            c[ i*nj + j ] += a[ k*ni + i ]*b[ j*nk + k ];
                            c[ i*nj + j ] += a[ (k + 1)*ni + i ]*b[ j*nk + k + 1 ];
                    }
                }
            }
        }
    }
    }
}

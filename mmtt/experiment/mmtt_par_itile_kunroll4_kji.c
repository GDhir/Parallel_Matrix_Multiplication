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
    for (ii = 0; ii < ni; ii += bs) {

        for( k = 0; k < remk; k++ ) {
            for( j = 0; j < nj; j++ ) {
                for (i = ii; i < min( ni, ii + bs ); i++) {
                    c[ i*nj + j ] += a[ k*ni + i ]*b[ j*nk + k ];
                }
            }
        }

        for( k = remk; k < nk; k += 4 ) {
            for( j = 0; j < nj; j++ ) {
                for (i = ii; i < min( ni, ii + bs ); i++) {
                    c[ i*nj + j ] += a[ k*ni + i ]*b[ j*nk + k ];
                    c[ i*nj + j ] += a[ (k + 1)*ni + i ]*b[ j*nk + (k + 1) ];
                    c[ i*nj + j ] += a[ (k + 2)*ni + i ]*b[ j*nk + (k + 2) ];
                    c[ i*nj + j ] += a[ (k + 3)*ni + i ]*b[ j*nk + (k + 3) ];
                }
            }
        }
    }
}
}

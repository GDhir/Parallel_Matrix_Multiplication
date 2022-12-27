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
   int remk;

#pragma omp parallel private(i, j, k, kk, remk)
   {
    #pragma omp for
    for (i = 0; i < ni; i++) {
        for(kk = 0; kk < nk; kk += bs ) {
            for( j = 0; j < nj; j++ ) {

                remk = (min( nk, kk + bs ) - kk)%8;
                for( k = kk; k < kk + remk; k++ ) {
                        c[ i*nj + j ] += a[ k*ni + i ]*b[ j*nk + k ];
                }

                for(k = kk + remk; k < min( nk, kk + bs ); k+= 8 ) {
                        c[ i*nj + j ] += a[ k*ni + i ]*b[ j*nk + k ];
                        c[ i*nj + j ] += a[ (k + 1)*ni + i ]*b[ j*nk + k + 1 ];
                        c[ i*nj + j ] += a[ (k + 2)*ni + i ]*b[ j*nk + k + 2 ];
                        c[ i*nj + j ] += a[ (k + 3)*ni + i ]*b[ j*nk + k + 3 ];
                        c[ i*nj + j ] += a[ (k + 4)*ni + i ]*b[ j*nk + k + 4 ];
                        c[ i*nj + j ] += a[ (k + 5)*ni + i ]*b[ j*nk + k + 5 ];
                        c[ i*nj + j ] += a[ (k + 6)*ni + i ]*b[ j*nk + k + 6 ];
                        c[ i*nj + j ] += a[ (k + 7)*ni + i ]*b[ j*nk + k + 7 ];
                }
            }
        }
    }
    }
}

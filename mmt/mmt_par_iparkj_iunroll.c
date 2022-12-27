void mmt_par(int ni, int nj, int nk, double *__restrict__ a, double *__restrict__ b, double *__restrict__ c)
{
   int i, j, k;

#pragma omp parallel private(i, j, k)
   {
      int remi = ni%4;
      
      for (i = 0; i < remi; i++) {
         for(k = 0; k < nk; k++ ) {
            for( j = 0; j < nj; j++ ) {
               c[ i*nj + j ] += a[ k*ni + i ]*b[ k*nj + j ];
            }
         }
      }

      #pragma omp for
      for (i = remi; i < ni; i += 4) {
         for(k = 0; k < nk; k++ ) {
            for( j = 0; j < nj; j++ ) {
               c[ i*nj + j ] += a[ k*ni + i ]*b[ k*nj + j ];
               c[ ( i + 1 )*nj + j ] += a[ k*ni + i + 1 ]*b[ k*nj + j ];
               c[ ( i + 2 )*nj + j ] += a[ k*ni + i + 2 ]*b[ k*nj + j ];
               c[ ( i + 3 )*nj + j ] += a[ k*ni + i + 3 ]*b[ k*nj + j ];
            }
         }
      }
   }
}

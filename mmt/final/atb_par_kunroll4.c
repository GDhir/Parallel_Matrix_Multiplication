void atb_par_kunroll4(const double *__restrict__ a, const double *__restrict__ b, double *__restrict__ c, int ni, int nj, int nk)
{
   int i, j, k;

#pragma omp parallel private(i, j, k)
   {
      int remk = nk%4;
      #pragma omp for
      for (i = 0; i < ni; i++) {

         for(k = 0; k < remk; k++ ) {
            for( j = 0; j < nj; j++ ) {
               c[ i*nj + j ] += a[ k*ni + i ]*b[ k*nj + j ];
            }
         }

         for( k = remk; k < nk; k += 4 ) {
            for( j = 0; j < nj; j++ ) {
               c[ i*nj + j ] += a[ k*ni + i ]*b[ k*nj + j ];
               c[ i*nj + j ] += a[ (k + 1)*ni + i ]*b[ (k + 1)*nj + j ];
               c[ i*nj + j ] += a[ (k + 2)*ni + i ]*b[ (k + 2)*nj + j ];
               c[ i*nj + j ] += a[ (k + 3)*ni + i ]*b[ (k + 3)*nj + j ];
            //    c[ i*nj + j ] += a[ (k + 4)*ni + i ]*b[ (k + 4)*nj + j ];
            }
         }
      }
   }
}

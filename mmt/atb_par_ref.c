void atb_par_ref(double *__restrict__ a, double *__restrict__ b, double *__restrict__ c, int ni, int nj, int nk)
{
   int i, j, k;

#pragma omp parallel private(i, j, k)
   {
      #pragma omp for
      for (i = 0; i < ni; i++) {
         for(k = 0; k < nk; k++ ) {
            for( j = 0; j < nj; j++ ) {
               c[ i*nj + j ] += a[ k*ni + i ]*b[ k*nj + j ];
            }
         }
      }
   }
}

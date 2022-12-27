void atbt_par_ref( const double *__restrict__ a, const double *__restrict__ b, double *__restrict__ c, int ni, int nj, int nk )
{
   int i, j, k;

#pragma omp parallel private(i, j, k)
   {
    #pragma omp for
      for (i = 0; i < ni; i++) {
         for( k = 0; k < nk; k++ ) {
            for(j = 0; j < nj; j++) {
               c[i * nj + j] += a[k * ni + i] * b[j * nk + k];
            }
         }
      }
   }
}

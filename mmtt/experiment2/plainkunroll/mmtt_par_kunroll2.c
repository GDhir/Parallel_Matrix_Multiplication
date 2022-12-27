void mmtt_par(int ni, int nj, int nk, double *__restrict__ a, double *__restrict__ b, double *__restrict__ c)
{
   int i, j, k;

#pragma omp parallel private(i, j, k)
   {
    int remk = nk%2;
    #pragma omp for
      for (i = 0; i < ni; i++) {

         for( k = 0; k < remk; k++ ) {
            for(j = 0; j < nj; j++) {
               c[i * nj + j] += a[k * ni + i] * b[j * nk + k];
            }
         }

         for(k = remk; k < nk; k += 2) {
            for(j = 0; j < nj; j++) {
               c[i * nj + j] += a[k * ni + i] * b[j * nk + k];
               c[i * nj + j] += a[(k + 1) * ni + i] * b[j * nk + k + 1];
            }
         }
      }
   }
}

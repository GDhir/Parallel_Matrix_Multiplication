void mmtt_par(int ni, int nj, int nk, double *__restrict__ a, double *__restrict__ b, double *__restrict__ c)
{
   int i, j, k;

#pragma omp parallel private(i, j, k)
   {
    //   int remk = n%4;
    #pragma omp for
      for (j = 0; j < nj; j++) {
         for(k = 0; k < nk; k++) {
            for(i = 0; i < ni; i++) {
               c[i * nj + j] += a[k * ni + i] * b[j * nk + k];
            }
         }
      }
   }
}

void atbt_par_kunroll8(const double *__restrict__ a, const double *__restrict__ b, double *__restrict__ c, int ni, int nj, int nk)
{
   int i, j, k;

#pragma omp parallel private(i, j, k)
   {
    int remk = nk%8;
    #pragma omp for
      for (i = 0; i < ni; i++) {

         for( k = 0; k < remk; k++ ) {
            for(j = 0; j < nj; j++) {
               c[i * nj + j] += a[k * ni + i] * b[j * nk + k];
            }
         }

         for(k = remk; k < nk; k += 8) {
            for(j = 0; j < nj; j++) {
               c[i * nj + j] += a[k * ni + i] * b[j * nk + k];
               c[i * nj + j] += a[(k + 1) * ni + i] * b[j * nk + k + 1];
               c[i * nj + j] += a[(k + 2) * ni + i] * b[j * nk + k + 2];
               c[i * nj + j] += a[(k + 3) * ni + i] * b[j * nk + k + 3];
               c[i * nj + j] += a[(k + 4) * ni + i] * b[j * nk + k + 4];
               c[i * nj + j] += a[(k + 5) * ni + i] * b[j * nk + k + 5];
               c[i * nj + j] += a[(k + 6) * ni + i] * b[j * nk + k + 6];
               c[i * nj + j] += a[(k + 7) * ni + i] * b[j * nk + k + 7];
            }
         }
      }
   }
}

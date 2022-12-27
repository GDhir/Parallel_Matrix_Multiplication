void mmt_par(int ni, int nj, int nk, double *__restrict__ a, double *__restrict__ b, double *__restrict__ c)
{
   int i, j, k;

#pragma omp parallel private(i, j, k)
   {
    double sum = 0;
      
      for (i = 0; i < ni; i++) {
            for( j = 0; j < nj; j++ ) {

                sum = 0;

                #pragma omp for 
                for(k = 0; k < nk; k++ ) {
                    sum = sum + a[ k*ni + i ]*b[ k*nj + j ];
                }

                #pragma omp atomic
                c[ i*nj + j ] += sum;

            }
      }
   }
}

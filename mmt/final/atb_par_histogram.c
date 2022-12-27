void atb_par_histogram(const double *__restrict__ a, const double *__restrict__ b, double *__restrict__ c, int Ni, int Nj, int Nk)
{
   int i, j, k;

#pragma omp parallel private(i, j, k)
   {
      double cval[Ni*Nj];

      for(i = 0; i < Ni*Nj; i++ ) {
         cval[i] = 0;
      }

      #pragma omp for 
      for(k = 0; k < Nk; k++ ) {
         for (i = 0; i < Ni; i++) {
            for( j = 0; j < Nj; j++ ) {
               cval[ i*Nj + j ] += a[ k*Ni + i ]*b[ k*Nj + j ];
            }
         }
      }

      for (i = 0; i < Ni; i++) {
         for( j = 0; j < Nj; j++ ) {
            
            #pragma omp atomic
            c[ i*Nj + j ] += cval[ i*Nj + j ];

         }
      }
   }
}

void atbt_par_histogram(const double *__restrict__ A, const double *__restrict__ B, double *__restrict__ C, int Ni, int Nj, int Nk)
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
               cval[ i*Nj + j ] += A[ k*Ni + i ]*B[ j*Nk + k ];
            }
         }
      }

      for (i = 0; i < Ni; i++) {
         for( j = 0; j < Nj; j++ ) {
            
            #pragma omp atomic
            C[ i*Nj + j ] += cval[ i*Nj + j ];

         }
      }
   }
}

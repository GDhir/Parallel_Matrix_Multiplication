// matrix multiply kernel: C = A^T * B
__global__ void atb(const double *A, const double *B, double *C, int Ni, int Nj, int Nk) {
    // Initially empty; will clearly not pass correctness
    
        int i = 2*blockIdx.y;
        int j = 2*blockIdx.x*blockDim.x+threadIdx.x;
        
        if ( ( i + 1 < Ni ) && ( j + blockDim.x < Nj ) ) { 
            double sum0 = 0;
            double sum1 = 0;
            double sum2 = 0;
            double sum3 = 0;
    
            for (int k = 0; k < Nk; ++k) {
                 sum0 += A[ k*Ni + i ]*B[ k*Nj + j ]; 
                 sum1 += A[ k*Ni + i + 1 ]*B[ k*Nj + j ]; 
                 sum2 += A[ k*Ni + i ]*B[ k*Nj + j + blockDim.x ]; 
                 sum3 += A[ k*Ni + i + 1 ]*B[ k*Nj + j + blockDim.x ]; 
            }
                 
            C[ i * Nj + j ] = sum0;
            C[ (i + 1) * Nj + j ] = sum1;
            C[ i * Nj + j + blockDim.x ] = sum2;
            C[ (i + 1) * Nj + j + blockDim.x ] = sum3;
        }
        else if( i + 1 < Ni && j < Nj ) {

            double sum0 = 0;
            double sum1 = 0;
    
            for (int k = 0; k < Nk; ++k) {
                 sum0 += A[ k*Ni + i ]*B[ k*Nj + j ]; 
                 sum1 += A[ k*Ni + i + 1 ]*B[ k*Nj + j ]; 
            }
                 
            C[ i * Nj + j ] = sum0;
            C[ (i + 1) * Nj + j ] = sum1;

        }
        else if( ( i < Ni ) && ( j + blockDim.x < Nj ) ) {

            double sum0 = 0;
            double sum1 = 0;
    
            for (int k = 0; k < Nk; ++k) {
                 sum0 += A[ k*Ni + i ]*B[ k*Nj + j ]; 
                 sum1 += A[ k*Ni + i ]*B[ k*Nj + j + blockDim.x ]; 
            }
                 
            C[ i * Nj + j ] = sum0;
            C[ i * Nj + j + blockDim.x ] = sum1;

        }
        else if( i < Ni && j < Nj ) {

            double sum = 0;

            for (int k = 0; k < Nk; ++k) {
                sum += A[ k*Ni + i ]*B[ k*Nj + j ]; 
            }
                
            C[ i * Nj + j ] = sum;

        }
    
    }
    
    
// matrix multiply kernel: C = A^T * B
__global__ void atbt(const double *A, const double *B, double *C, int Ni, int Nj, int Nk) {
    // Initially empty; will clearly not pass correctness
    
        int j = blockIdx.y;
        int i = 2*blockIdx.x*blockDim.x+threadIdx.x;
        
        if ( ( i + blockDim.x < Ni ) && ( j < Nj ) ) { 
            double sum0 = 0;
            double sum1 = 0;
    
            for (int k = 0; k < Nk; ++k) {
                sum0 += A[ k*Ni + i ]*B[ j*Nk + k ]; 
                sum1 += A[ k*Ni + i + blockDim.x ]*B[ j*Nk + k ]; 
            }
                 
            C[ i * Nj + j ] = sum0;
            C[ (i + blockDim.x) * Nj + j ] = sum1;
        }
        else if( ( i < Ni ) && ( j < Nj ) ) {
            double sum0 = 0;
    
            for (int k = 0; k < Nk; ++k) {
                sum0 += A[ k*Ni + i ]*B[ j*Nk + k ]; 
            }
                 
            C[ i * Nj + j ] = sum0;
        }
    
    }
    
    
// matrix multiply kernel: C = A^T * B
__global__ void atb(const double *A, const double *B, double *C, int Ni, int Nj, int Nk) {
    // Initially empty; will clearly not pass correctness
    
        int i = blockIdx.y*blockDim.y+threadIdx.y;
        int j = 2*blockIdx.x*blockDim.x+threadIdx.x;
        
        if ( ( i < Ni ) && ( j < Nj ) ) { 
            double sum0 = 0;
            double sum1 = 0;
    
            if( j + blockDim.x < Nj ) {
                for (int k = 0; k < Nk; ++k) {
                    sum0 += A[ k*Ni + i ]*B[ k*Nj + j ]; 
                    sum1 += A[ k*Ni + i ]*B[ k*Nj + j + blockDim.x ]; 
                }
            }
            else {

                for (int k = 0; k < Nk; ++k) {
                    sum0 += A[ k*Ni + i ]*B[ k*Nj + j ]; 
                }

            }

            C[ i * Nj + j ] = sum0;

            if( j + blockDim.x < Nj ) {
                C[ i * Nj + j + blockDim.x ] = sum1;
            }
        }
    
    }
    
    
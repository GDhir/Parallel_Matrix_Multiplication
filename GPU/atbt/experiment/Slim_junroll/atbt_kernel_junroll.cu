// matrix multiply kernel: C = A^T * B
__global__ void atbt(const double *A, const double *B, double *C, int Ni, int Nj, int Nk) {
    // Initially empty; will clearly not pass correctness
    
        int j = 2*blockIdx.y;
        int i = blockIdx.x*blockDim.x+threadIdx.x;
        
        if ( ( i < Ni ) && ( j + 1 < Nj ) ) { 
            double sum0 = 0;
            double sum1 = 0;
    
            for (int k = 0; k < Nk; ++k) {
                sum0 += A[ k*Ni + i ]*B[ j*Nk + k ]; 
                sum1 += A[ k*Ni + i ]*B[ (j + 1)*Nk + k ]; 
            }
                 
            C[ i * Nj + j ] = sum0;
            C[ i * Nj + j + 1 ] = sum1;
        }
        else if( i < Ni && j < Nj ) {

            double sum0 = 0;
    
            for (int k = 0; k < Nk; ++k) {
                sum0 += A[ k*Ni + i ]*B[ j*Nk + k ]; 
            }
                    
            C[ i * Nj + j ] = sum0;

        }
    
    }
    
    
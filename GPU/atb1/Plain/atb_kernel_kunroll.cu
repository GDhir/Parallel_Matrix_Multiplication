// matrix multiply kernel: C = A^T * B
__global__ void atb(const double *A, const double *B, double *C, int Ni, int Nj, int Nk) {
    // Initially empty; will clearly not pass correctness
    
        int i = blockIdx.y*blockDim.y+threadIdx.y;
        int j = blockIdx.x*blockDim.x+threadIdx.x;
        
        if ( ( i < Ni ) && ( j < Nj ) ) { 
            double sum = 0;
    
            int remk = Nk%2;
            
            for( int k = 0; k < remk; k++ ) {
                sum += A[ k*Ni + i ]*B[ k*Nj + j ]; 
            }

            for (int k = remk; k < Nk; k += 2) {
                 sum += A[ k*Ni + i ]*B[ k*Nj + j ]; 
                 sum += A[ (k + 1)*Ni + i ]*B[ (k + 1)*Nj + j ]; 
            }
                 
            C[ i * Nj + j ] = sum;
        }
    
    }
    
    
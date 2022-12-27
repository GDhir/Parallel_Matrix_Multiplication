// matrix multiply kernel: C = A^T * B
__global__ void atb(const double *A, const double *B, double *C, int Ni, int Nj, int Nk) {
    // Initially empty; will clearly not pass correctness
    
        int i = blockIdx.y;
        int j = blockIdx.x*blockDim.x+threadIdx.x;
        
        int remk = Nk%4;

        if ( ( i < Ni ) && ( j < Nj ) ) { 
            double sum = 0;
    
            for( int k = 0; k < remk; k++ ) {
                sum += A[ k*Ni + i ]*B[ k*Nj + j ]; 
            }

            for (int k = remk; k < Nk; k += 4) {
                 sum += A[ k*Ni + i ]*B[ k*Nj + j ]; 
                 sum += A[ (k + 1)*Ni + i ]*B[ (k + 1)*Nj + j ]; 
                 sum += A[ (k + 2)*Ni + i ]*B[ (k + 2)*Nj + j ]; 
                 sum += A[ (k + 3)*Ni + i ]*B[ (k + 3)*Nj + j ]; 


            }
                 
            C[ i * Nj + j ] = sum;
        }
    
    }
    
    
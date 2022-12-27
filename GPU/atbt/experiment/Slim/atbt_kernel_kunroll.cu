// matrix multiply kernel: C = A^T * B
__global__ void atbt(const double *A, const double *B, double *C, int Ni, int Nj, int Nk) {
    // Initially empty; will clearly not pass correctness
    
        int j = blockIdx.y;
        int i = blockIdx.x*blockDim.x+threadIdx.x;
        
        int remk = Nk%4;

        if ( ( i < Ni ) && ( j < Nj ) ) { 
            double sum = 0;
    
            for( int k = 0; k < remk; k++ ) {
                sum += A[ k*Ni + i ]*B[ j*Nk + k ]; 
            }

            for (int k = remk; k < Nk; k += 4) {
                sum += A[ k*Ni + i ]*B[ j*Nk + k ]; 
                sum += A[ (k + 1)*Ni + i ]*B[ j*Nk + k + 1 ]; 
                sum += A[ (k + 2)*Ni + i ]*B[ j*Nk + k + 2 ]; 
                sum += A[ (k + 3)*Ni + i ]*B[ j*Nk + k + 3 ]; 

            }
                 
            C[ i * Nj + j ] = sum;
        }
    
    }
    
    
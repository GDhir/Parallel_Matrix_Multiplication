// matrix multiply kernel: C = A^T * B
__global__ void atb(const double *A, const double *B, double *C, int Ni, int Nj, int Nk) {
// Initially empty; will clearly not pass correctness

    int i = blockIdx.y;
    int j = blockIdx.x*blockDim.x+threadIdx.x;
    
    if ( ( i < Ni ) && ( j < Nj ) ) { 
        double sum = 0;

        for (int k = 0; k < Nk; ++k) {
             sum += A[ k*Ni + i ]*B[ k*Nj + j ]; 
        }
             
        C[ i * Nj + j ] = sum;
    }

}


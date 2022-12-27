// matrix multiply kernel: C = A^T * B
__global__ void atbt_kernel_slim(const double *A, const double *B, double *C, int Ni, int Nj, int Nk) {
// Initially empty; will clearly not pass correctness

    int j = blockIdx.y;
    int i = blockIdx.x*blockDim.x+threadIdx.x;
    
    if ( ( i < Ni ) && ( j < Nj ) ) { 
        double sum = 0;

        for (int k = 0; k < Nk; ++k) {
             sum += A[ k*Ni + i ]*B[ j*Nk + k ]; 
        }
             
        C[ i * Nj + j ] = sum;
    }

}


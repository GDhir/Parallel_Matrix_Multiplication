 // matrix multiply kernel: C = A^T * B^T
__global__ void atbt(const double *A, const double *B, double *C, int Ni, int Nj, int Nk) {
    // Initially empty; clearly Will not pass correctness test
    
    int i = blockIdx.x*blockDim.x + threadIdx.x;
    int j = blockIdx.y*blockDim.y + threadIdx.y;

    if( i < Ni && j < Nj ) {

        double csub = 0;
    
        for( int k = 0; k < Nk; k++ ) {
            csub += A[ k*Ni + i ]*B[ j*Nk + k ];
        }

        C[ i*Nj + j ] = csub;
    }
        
    
}
    
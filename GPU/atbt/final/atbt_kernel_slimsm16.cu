#define BSIZE 16

// matrix multiply kernel: C = A^T * B^T
__global__ void atbt_kernel_slimsm16(const double *A, const double *B, double *C, int Ni, int Nj, int Nk) {
    // Initially empty; will clearly not pass correctness

        int tx = threadIdx.x;
        int j = blockIdx.y;
        int i = blockIdx.x*blockDim.x+tx;
        int idx = j*Nk + tx;
        
        __shared__ double Bs[BSIZE];

        int kend = Nk*(j + 1);

        double sum = 0;

        for( int ks = 0; ks < Nk; ks += BSIZE ) {

            if( j < Nj && idx + ks < kend )
                Bs[tx] = B[ idx + ks ];
            else
                Bs[tx] = 0;

            __syncthreads();

            if( i < Ni ) {
                for (int k = ks; k < ks + BSIZE; ++k) {
                    sum += A[k*Ni + i]*Bs[ k - ks ]; 
                }
            }

            __syncthreads();

        }

        if( i < Ni && j < Nj )
            C[ i * Nj + j ] = sum;
    
    }
    
    
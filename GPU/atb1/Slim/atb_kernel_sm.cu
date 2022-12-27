#define BSIZE 32

// matrix multiply kernel: C = A^T * B
__global__ void atb(const double *A, const double *B, double *C, int Ni, int Nj, int Nk) {
    // Initially empty; will clearly not pass correctness

        int tx = threadIdx.x;
        int i = blockIdx.y;
        int j = blockIdx.x*blockDim.x+tx;
        int idx = tx*Ni + i;

        int aEnd = Nk*Ni;
        
        __shared__ double As[BSIZE];

        // if ( ( i < Ni ) && ( j < Nj ) ) { 
            double sum = 0;
    
            for( int ks = 0; ks < Nk; ks += BSIZE ) {

                if( i < Ni && idx + ks*Ni < aEnd )
                    As[tx] = A[ ks*Ni + idx ];
                else
                    As[tx] = 0;

                __syncthreads();

                if( j < Nj ) {
                    for (int k = ks; k < min(Nk, ks + BSIZE); ++k) {
                        sum += As[k - ks]*B[ k*Nj + j ]; 
                    }
                }

                __syncthreads();

            }

            if( i < Ni && j < Nj )
                C[ i * Nj + j ] = sum;
        // }
    
    }
    
    
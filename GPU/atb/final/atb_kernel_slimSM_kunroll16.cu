#define BSIZE 16

// matrix multiply kernel: C = A^T * B
__global__ void atb_kernel_slimSM_kunroll16(const double *A, const double *B, double *C, int Ni, int Nj, int Nk) {
    // Initially empty; will clearly not pass correctness

        int tx = threadIdx.x;
        int i = blockIdx.y;
        int j = blockIdx.x*blockDim.x+tx;
        int idx = tx*Ni + i;
        int remk;

        int aEnd = Nk*Ni;
        
        __shared__ double As[BSIZE];

        double sum = 0;

        for( int ks = 0; ks < Nk; ks += BSIZE ) {

            if( idx + ks*Ni < aEnd )
                As[tx] = A[ ks*Ni + idx ];
            else
                As[tx] = 0;

            __syncthreads();


            if( j < Nj ) {
                for(int k = ks; k < ks + remk; k++ ) {
                    sum += As[k - ks]*B[ k*Nj + j ];
                }

                for (int k = ks + remk; k < ks + remk + BSIZE; k += 2) {
                    sum += As[k - ks]*B[ k*Nj + j ];
                    sum += As[k + 1 - ks]*B[ (k + 1)*Nj + j ]; 
                }
            }

            __syncthreads();

        }
            
        if( i < Ni && j < Nj ) {
            C[ i * Nj + j ] = sum;
        }
    
    }
    
    
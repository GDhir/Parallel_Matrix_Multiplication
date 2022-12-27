__global__ void atb(const double *A, const double *B, double *C, int Ni, int Nj, int Nk) {

    int bx = blockIdx.x;    int by = blockIdx.y;
    int tx = threadIdx.x;  int ty = threadIdx.y;
    // Starting index of A and B for the thread block
    int aBegin = blockDim.y*by;
    int bBegin = blockDim.y*bx;
    double Csub = 0;
    // Declaration of shared memory buffers 
    __shared__ double As[2][2];
    __shared__ double Bs[2][2];

    // Starting index of A and B for the thread
    // int aInd = aBegin + Ni * tx + ty;
    // int bInd = bBegin + Nj * ty + tx;

    int acol = blockDim.y*by + ty;
    int arow = tx;
    int bcol = blockDim.y*bx + tx;
    int brow = ty;

    // int cInd = bx*blockDim.x + by*blockDim.y*Nj + ty*Nj + tx;

    int col = bx*blockDim.x + tx;
    int row = by*blockDim.y + ty;

    for (int kt = 0; kt < Nk; kt += blockDim.y) {

        if( acol < Ni && arow < Nk ) {
            As[tx][ty] = A[ arow*Ni + acol ];
            
        }
        else {
            As[tx][ty] = 0;

        }

        if( bcol < Nj && brow < Nk ) {
            Bs[ty][tx] = B[ brow*Nj + bcol ];
        }
        else {
            Bs[ty][tx] = 0;
        }

        __syncthreads();

        if( col < Nj && row < Ni ) {
            for (int k = kt; k < kt + blockDim.y; ++k) { 
            
                Csub += As[k - kt][ty] * Bs[k - kt][tx];
            
            }
        }
        
        __syncthreads();
        arow += blockDim.y;
        brow += blockDim.y;
    }
    if( col < Nj && row < Ni ) {
        C[ row*Nj + col ] = Csub;
    }
}
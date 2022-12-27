__global__ void atb(const double *A, const double *B, double *C, int Ni, int Nj, int Nk) {

    int bx = blockIdx.x;    int by = blockIdx.y;
    int tx = threadIdx.x;  int ty = threadIdx.y;
    // Starting index of A and B for the thread block
    int aBegin = blockDim.y*by;
    int bBegin = blockDim.y*bx;
    double Csub = 0;
    // Declaration of shared memory buffers 
    __shared__ double As[32][32];
    __shared__ double Bs[32][32];

    // Starting index of A and B for the thread
    int aInd = aBegin + Ni * tx + ty;
    int bInd = bBegin + Nj * ty + tx;

    int cBegin = bx*blockDim.x + by*blockDim.y*Nj;
    int cInd = cBegin + ty*Nj + tx;


    int aEnd = Ni*Nk;
    int bEnd = Nj*Nk;
    int cEnd = Ni*Nj;

    for (int kt = 0; kt < Nk; kt += blockDim.y) {

        if( aBegin + ty < Ni && aInd < aEnd ) {
            As[tx][ty] = A[ aInd ];
            
        }
        else {
            As[tx][ty] = 0;

        }

        if( bBegin + tx < Nj && bInd < bEnd ) {
            Bs[ty][tx] = B[ bInd ];
        }
        else {
            Bs[ty][tx] = 0;
        }

        __syncthreads();

        if( bx*blockDim.x + tx < Nj && cInd < cEnd ) {
            for (int k = kt; k < kt + blockDim.y; ++k) { 
            
                Csub += As[k - kt][ty] * Bs[k - kt][tx];
            
            }
        }
        
        __syncthreads();
        aInd += Ni*blockDim.y;
        bInd += Nj*blockDim.y;
    }

    if( bx*blockDim.x + tx < Nj && cInd < cEnd ) {
        C[ cInd ] = Csub;
    }
}
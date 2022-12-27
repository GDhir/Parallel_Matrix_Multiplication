__global__ void atb(const double *A, const double *B, double *C, int Ni, int Nj, int Nk, int blk) {

    int bx = blockIdx.x;    int by = blockIdx.y;
    int tx = threadIdx.x;  int ty = threadIdx.y;
    // Starting index of A and B for the thread block
    int aBegin = blockDim.y*by;
    int bBegin = blockDim.y*bx;
    double Csub = 0;
    // Declaration of shared memory buffers 
    extern __shared__ double As[];
    double *Bs = (double *)As + (blk*blk);

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
            As[tx*blk + ty] = A[ aInd ];
            
        }
        else {
            As[tx*blk + ty] = 0;

        }

        if( bBegin + tx < Nj && bInd < bEnd ) {
            Bs[ty*blk + tx] = B[ bInd ];
        }
        else {
            Bs[ty*blk + tx] = 0;
        }

        __syncthreads();

        if( bx*blockDim.x + tx < Nj && cInd < cEnd ) {
            for (int k = kt; k < kt + blockDim.y; ++k) { 
            
                Csub += As[(k - kt)*blk + ty] * Bs[(k - kt)*blk + tx];
            
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
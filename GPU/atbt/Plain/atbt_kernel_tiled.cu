#define BSIZE 32

__global__ void atbt(const double *A, const double *B, double *C, int Ni, int Nj, int Nk) {

    int bx = blockIdx.x;    int by = blockIdx.y;
    int tx = threadIdx.x;  int ty = threadIdx.y;
    // Starting index of A and B for the thread block
    int aBegin = blockDim.y*bx;
    int bBegin = blockDim.y*by*Nk;
    double Csub = 0;
    // Declaration of shared memory buffers 
    __shared__ double As[BSIZE][BSIZE];
    __shared__ double Bs[BSIZE][BSIZE];

    // Starting index of A and B for the thread
    int aInd = aBegin + Ni * ty + tx;
    int bInd = bBegin + Nk * ty + tx;

    int cBegin = bx*blockDim.x*Nj + by*blockDim.y;
    int cInd = cBegin + tx*Nj + ty;

    int aEnd = Ni*Nk;
    int bEnd = Nj*Nk;
    int cEnd = Ni*Nj;

    int rowend = ( blockDim.y*by + ty + 1 )*Nk;

    for (int kt = 0; kt < Nk; kt += BSIZE) {

        if( aBegin + tx < Ni && aInd < aEnd ) {
            As[ty][tx] = A[ aInd ];
            
        }
        else {
            As[ty][tx] = 0;

        }

        if( bInd < rowend && bInd < bEnd ) {
            Bs[ty][tx] = B[ bInd ];
        }
        else {
            Bs[ty][tx] = 0;
        }

        __syncthreads();

        for (int k = kt; k < kt + blockDim.y; ++k) { 
        
            Csub += As[k - kt][tx] * Bs[ty][k - kt];
        
        }
        
        __syncthreads();
        aInd += Ni*blockDim.y;
        bInd += BSIZE;
    }

    if( by*blockDim.y + ty < Nj && cInd < cEnd ) {
        C[ cInd ] = Csub;
    }
}
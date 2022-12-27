#define BSIZE 32

__global__ void atb(const double *A, const double *B, double *C, int Ni, int Nj, int Nk) {

    int bx = blockIdx.x;    int by = 2*blockIdx.y;
    int tx = threadIdx.x;  int ty = 2*threadIdx.y;
    // Starting index of A and B for the thread block
    int aBegin = blockDim.y*by;
    int bBegin = blockDim.y*bx;
    double Csub0 = 0;
    double Csub1 = 0;
    // Declaration of shared memory buffers 
    __shared__ double As[BSIZE][2*BSIZE];
    __shared__ double Bs[BSIZE][BSIZE];

    // Starting index of A and B for the thread
    int aInd = aBegin + Ni * tx + ty;
    int bInd = bBegin + Nj * threadIdx.y + tx;

    int cBegin = bx*blockDim.x + by*blockDim.y*Nj;
    int cInd = cBegin + ty*Nj + tx;

    int aEnd = Ni*Nk;
    int bEnd = Nj*Nk;
    int cEnd = Ni*Nj;

    for (int kt = 0; kt < Nk; kt += BSIZE) {

        if( aBegin + ty + 1 < Ni && aInd < aEnd ) {
            As[tx][ty] = A[ aInd ];
            As[tx][ty + 1] = A[ aInd + 1 ];
            
        }
        else if( aBegin + ty < Ni && aInd < aEnd ) {
            As[tx][ty] = A[ aInd ];
            As[tx][ty + 1] = 0;
        }
        else {
            As[tx][ty] = 0;
            As[tx][ty + 1] = 0;
        }

        if( bBegin + tx < Nj && bInd < bEnd ) {
            Bs[threadIdx.y][tx] = B[ bInd ];
        }
        else {
            Bs[threadIdx.y][tx] = 0;
        }

        __syncthreads();

        // if( bx*blockDim.x + tx < Nj && cInd < cEnd ) {

        for (int k = kt; k < kt + BSIZE; k++) { 
        
            Csub0 += As[k - kt][ty]*Bs[k - kt][tx];
            Csub1 += As[k - kt][ty + 1] * Bs[k - kt][tx];
        
        }
        // }
        // else if( bx*blockDim.x + tx < Nj && cInd < cEnd ) {

        //     for (int k = kt; k < kt + BSIZE; k++) { 
            
        //         Csub0 += As[k - kt][ty] * Bs[k - kt][tx];
            
        //     }

        // }
        
        __syncthreads();
        aInd += Ni*blockDim.y;
        bInd += Nj*blockDim.y;
    }

    if( bx*blockDim.x + tx < Nj && cInd + Nj < cEnd ) {
        C[ cInd ] = Csub0;
        C[ cInd + Nj ] = Csub1;
    }
    else if( bx*blockDim.x + tx < Nj && cInd < cEnd ) {
        C[ cInd ] = Csub0;
    }
}
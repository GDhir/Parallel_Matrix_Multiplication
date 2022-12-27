#define BSIZE 32

__global__ void atbt_kernel_tiled_junroll(const double *A, const double *B, double *C, int Ni, int Nj, int Nk) {

    int bx = blockIdx.x;    int by = 2*blockIdx.y;
    int tx = threadIdx.x;  int ty = 2*threadIdx.y;
    int tyreal = threadIdx.y;
    // Starting index of A and B for the thread block
    int aBegin = blockDim.y*bx;
    int bBegin = blockDim.y*by*Nk;
    double Csub0 = 0;
    double Csub1 = 0;
    // Declaration of shared memory buffers 
    __shared__ double As[BSIZE][BSIZE];
    __shared__ double Bs[2*BSIZE][BSIZE];

    // Starting index of A and B for the thread
    int aInd = aBegin + Ni * tyreal + tx;
    int bInd = bBegin + Nk * ty + tx;

    int cBegin = bx*blockDim.x*Nj + by*blockDim.y;
    int cInd = cBegin + tx*Nj + ty;

    int aEnd = Ni*Nk;
    int bEnd = Nj*Nk;
    int cEnd = Ni*Nj;

    int rowend = ( blockDim.y*by + ty + 1 )*Nk;

    for (int kt = 0; kt < Nk; kt += BSIZE) {

        if( aBegin + tx < Ni && aInd < aEnd ) {
            As[tyreal][tx] = A[ aInd ];
            
        }
        else {
            As[tyreal][tx] = 0;

        }

        if( bInd < rowend && bInd + Nk < bEnd ) {
            Bs[ty + 1][tx] = B[ bInd + Nk ];
            Bs[ty][tx] = B[ bInd ];
        }
        else if( bInd < rowend && bInd < bEnd ) {
            Bs[ty + 1][tx] = 0;
            Bs[ty][tx] = B[bInd];
        }
        else {
            Bs[ty + 1][tx] = 0;
            Bs[ty][tx] = 0;
        }

        __syncthreads();

        for (int k = kt; k < kt + blockDim.y; ++k) { 
        
            Csub0 += As[k - kt][tx] * Bs[ty][k - kt];
            Csub1 += As[k - kt][tx] * Bs[ty + 1][k - kt];
        
        }
        
        __syncthreads();
        aInd += Ni*blockDim.y;
        bInd += BSIZE;
    }

    if( by*blockDim.y + ty + 1 < Nj && cInd < cEnd ) {
        C[ cInd ] = Csub0;
        C[ cInd + 1 ] = Csub1;
    }
    else if( by*blockDim.y + ty < Nj && cInd < cEnd ) {
        C[ cInd ] = Csub0;
    }
}
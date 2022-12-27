#define BSIZE 32

__global__ void atbt(const double *A, const double *B, double *C, int Ni, int Nj, int Nk) {

    int bx = 2*blockIdx.x;    int by = 2*blockIdx.y;
    int tx = threadIdx.x;  int ty = 2*threadIdx.y;
    int tyreal = threadIdx.y;
    // Starting index of A and B for the thread block
    int aBegin = blockDim.y*bx;
    int bBegin = blockDim.y*by*Nk;
    double Csub0 = 0;
    double Csub1 = 0;
    double Csub2 = 0;
    double Csub3 = 0;
    // Declaration of shared memory buffers 
    __shared__ double As[BSIZE][2*BSIZE];
    __shared__ double Bs[2*BSIZE][BSIZE];

    // Starting index of A and B for the thread
    int aInd = aBegin + Ni * tyreal + tx;
    int bInd = bBegin + Nk * ty + tx;

    int cBegin = bx*blockDim.x*Nj + by*blockDim.y;
    int cInd = cBegin + tx*Nj + ty;

    int aEnd = Ni*Nk;
    int bEnd = Nj*Nk;
    int cEnd = Ni*Nj;

    int remk = BSIZE%2;

    int rowend = ( blockDim.y*by + ty + 1 )*Nk;

    for (int kt = 0; kt < Nk; kt += BSIZE) {

        if( aBegin + tx + BSIZE < Ni && aInd < aEnd ) {
            As[tyreal][tx + BSIZE] = A[ aInd + BSIZE ];
            As[tyreal][tx] = A[ aInd ];
            
        }
        else if( aBegin + tx < Ni && aInd < aEnd ) {
            As[tyreal][tx] = A[ aInd ];
            As[tyreal][tx + BSIZE] = 0;
        }
        else {
            As[tyreal][tx] = 0;
            As[tyreal][tx + BSIZE] = 0;
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

        for( int k = kt; k < kt + remk; k++ ) {
            Csub0 += As[k - kt][tx] * Bs[ty][k - kt];
            Csub1 += As[k - kt][tx + BSIZE] * Bs[ty][k - kt];
            Csub2 += As[k - kt][tx] * Bs[ty + 1][k - kt];
            Csub3 += As[k - kt][tx + BSIZE] * Bs[ty + 1][k - kt];
        }

        for (int k = kt + remk; k < kt + BSIZE; k += 2) { 
        
            Csub0 += As[k - kt][tx] * Bs[ty][k - kt];
            Csub1 += As[k - kt][tx + BSIZE] * Bs[ty][k - kt];
            Csub2 += As[k - kt][tx] * Bs[ty + 1][k - kt];
            Csub3 += As[k - kt][tx + BSIZE] * Bs[ty + 1][k - kt];

            Csub0 += As[k + 1 - kt][tx] * Bs[ty][k + 1 - kt];
            Csub1 += As[k + 1 - kt][tx + BSIZE] * Bs[ty][k + 1 - kt];
            Csub2 += As[k + 1 - kt][tx] * Bs[ty + 1][k + 1 - kt];
            Csub3 += As[k + 1 - kt][tx + BSIZE] * Bs[ty + 1][k + 1 - kt];
        
        }
        
        __syncthreads();
        aInd += Ni*BSIZE;
        bInd += BSIZE;
    }

    if( by*blockDim.y + ty + 1 < Nj && cInd + Nj < cEnd ) {
        C[ cInd ] = Csub0;
        C[ cInd + Nj*BSIZE ] = Csub1;
        C[ cInd + 1 ] = Csub2;
        C[ cInd + Nj*BSIZE + 1 ] = Csub3;
    }
    else if( by*blockDim.y + ty + 1 < Nj && cInd < cEnd ) {
        C[ cInd ] = Csub0;
        C[ cInd + 1 ] = Csub2;
    }
    else if( by*blockDim.y + ty < Nj && cInd + Nj < cEnd ) {
        C[ cInd ] = Csub0;
        C[ cInd + Nj*BSIZE ] = Csub1;
    }
    else if( by*blockDim.y + ty < Nj && cInd < cEnd ) {
        C[ cInd ] = Csub0;
    }
}
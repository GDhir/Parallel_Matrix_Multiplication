#define BSIZE 32

__global__ void atb(const double *A, const double *B, double *C, int Ni, int Nj, int Nk) {

    int bx = 2*blockIdx.x;    int by = 2*blockIdx.y;
    int tx = threadIdx.x;  int ty = 2*threadIdx.y;
    int tyreal = threadIdx.y;
    // Starting index of A and B for the thread block
    int aBegin = blockDim.y*by;
    int bBegin = blockDim.y*bx;
    double Csub0 = 0;
    double Csub1 = 0;
    double Csub2 = 0;
    double Csub3 = 0;
    // Declaration of shared memory buffers 
    __shared__ double As[BSIZE][2*BSIZE];
    __shared__ double Bs[BSIZE][2*BSIZE];

    // Starting index of A and B for the thread
    int aInd = aBegin + Ni * tx + ty;
    int bInd = bBegin + Nj * tyreal + tx;

    int cBegin = bx*blockDim.x + by*blockDim.y*Nj;
    int cInd = cBegin + ty*Nj + tx;

    int aEnd = Ni*Nk;
    int bEnd = Nj*Nk;
    int cEnd = Ni*Nj;

    int remk;
    int k;

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

        if( bBegin + tx + BSIZE < Nj && bInd < bEnd ) {
            Bs[tyreal][tx] = B[ bInd ];
            Bs[tyreal][tx + BSIZE] = B[ bInd + BSIZE ];
        }
        else if( bBegin + tx < Nj && bInd < bEnd ) {
            Bs[tyreal][tx] = B[ bInd ];
            Bs[tyreal][tx + BSIZE] = 0;
        }
        else {
            Bs[tyreal][tx] = 0;
            Bs[tyreal][tx + BSIZE] = 0;
        }

        __syncthreads();

        // if( bx*blockDim.x + tx + BSIZE < Nj && cInd + Nj < cEnd ) {

            remk = BSIZE%2;

            for( k = kt; k < kt + remk; k++ ) {

                Csub0 += As[k - kt][ty] * Bs[k - kt][tx];
                Csub1 += As[k - kt][ty] * Bs[k - kt][tx + BSIZE];
                Csub2 += As[k - kt][ty + 1] * Bs[k - kt][tx];
                Csub3 += As[k - kt][ty + 1] * Bs[k - kt][tx + BSIZE];

            }


            for ( k = kt + remk; k < kt + BSIZE; k += 2) { 
            
                Csub0 += As[k - kt][ty] * Bs[k - kt][tx];
                Csub1 += As[k - kt][ty] * Bs[k - kt][tx + BSIZE];
                Csub2 += As[k - kt][ty + 1] * Bs[k - kt][tx];
                Csub3 += As[k - kt][ty + 1] * Bs[k - kt][tx + BSIZE];

                Csub0 += As[k + 1 - kt][ty] * Bs[k + 1 - kt][tx];
                Csub1 += As[k + 1 - kt][ty] * Bs[k + 1 - kt][tx + BSIZE];
                Csub2 += As[k + 1 - kt][ty + 1] * Bs[k + 1 - kt][tx];
                Csub3 += As[k + 1 - kt][ty + 1] * Bs[k + 1 - kt][tx + BSIZE];
            
            }
        // }
        // else if( bx*blockDim.x + tx + BSIZE < Nj && cInd < cEnd ) {

        //     for (int k = kt; k < kt + BSIZE; k++) { 
            
        //         Csub0 += As[k - kt][ty] * Bs[k - kt][tx];
        //         Csub1 += As[k - kt][ty] * Bs[k - kt][tx + BSIZE];
            
        //     }
        // }
        // else if( bx*blockDim.x + tx < Nj && cInd + Nj < cEnd ) {

        //     for (int k = kt; k < kt + BSIZE; k++) { 
            
        //         Csub0 += As[k - kt][ty] * Bs[k - kt][tx];
        //         Csub2 += As[k - kt][ty + 1] * Bs[k - kt][tx];
            
        //     }

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

    if( bx*blockDim.x + tx + BSIZE < Nj && cInd + Nj < cEnd ) {
        C[ cInd ] = Csub0;
        C[ cInd + BSIZE ] = Csub1;
        C[ cInd + Nj ] = Csub2;
        C[ cInd + BSIZE + Nj ] = Csub3;
    }
    else if( bx*blockDim.x + tx + BSIZE < Nj && cInd < cEnd ) {
        C[ cInd ] = Csub0;
        C[ cInd + BSIZE ] = Csub1;
    }
    else if( bx*blockDim.x + tx < Nj && cInd + Nj < cEnd ) {
        C[ cInd ] = Csub0;
        C[ cInd + Nj ] = Csub2;
    }
    else if( bx*blockDim.x + tx < Nj && cInd < cEnd ) {
        C[ cInd ] = Csub0;
    }
}
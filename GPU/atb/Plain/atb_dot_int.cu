#define BSIZE 64

__global__ void atb(const int *A, const int *B, int *D, int* C, int Ni, int Nj, int Nk, int i, int j) {

    int bx = blockIdx.x;
    int tx = threadIdx.x;

    __shared__ int cache[ BSIZE ];

    int kId = bx*blockDim.x + tx;

    if( kId < Nk ) {
        cache[ tx ] = A[ kId*Ni + i ]*B[ kId*Nj + j ];
    }
    else {
        cache[ tx ] = 0;
    }

    int sz = ceil( BSIZE/2.0 );

    __syncthreads();

    while( sz > 1 ) {

        if( tx < sz ) {

            cache[tx] += cache[tx + sz];
            cache[tx + sz] = 0;

        }

        sz = ceil( sz/2.0 );

        __syncthreads();

    }

    if( tx == 0 ) {

        D[bx] = cache[0] + cache[1];
        // atomicAdd1( &C[ i*Nj + j ], D[bx] );
    }

    __syncthreads();

    if( bx == 0 ) {

        while( kId + blockDim.x < gridDim.x ) {

            D[tx] += D[kId + blockDim.x];
            kId += blockDim.x;

        }


        if( gridDim.x >= blockDim.x )
            sz = ceil( blockDim.x/2.0 );
        else
            sz = ceil( gridDim.x/2.0 );

        __syncthreads();

        while( sz > 1 ) {

            if( tx < sz ) {

                D[tx] += D[tx + sz];
                D[tx + sz] = 0;

            }

            sz = ceil( sz/2.0 );

            __syncthreads();

        }


        if( tx == 0 ) {

            C[i*Nj + j] = D[0] + D[1];
            
        }

    }
}
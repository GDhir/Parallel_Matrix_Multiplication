#define BSIZE 512

__global__ void atb_kernel_dot(const double *A, const double *B, double *D, double* C, int Ni, int Nj, int Nk, int i, int j) {

    int bx = blockIdx.x;
    int tx = threadIdx.x;

    __shared__ double cache[ BSIZE ];

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

}
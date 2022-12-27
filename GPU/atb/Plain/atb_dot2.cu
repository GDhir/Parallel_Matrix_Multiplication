#define BSIZE 64

__device__ double atomicAdd1(double* address, double val)
{
    unsigned long long int* address_as_ull = (unsigned long long int*)address;
    unsigned long long int old = *address_as_ull, assumed;
    do {
        assumed = old;
        old = atomicCAS(address_as_ull, assumed,
                __double_as_longlong(val + __longlong_as_double(assumed)));
    } while (assumed != old);
    return __longlong_as_double(old);
}

__global__ void atb(const double *A, const double *B, double *D, double* C, int Ni, int Nj, int Nk, int i, int j) {

    int bx = 2*blockIdx.x;
    int tx = threadIdx.x;

    __shared__ double cache[ 2*BSIZE ];

    int kId = bx*blockDim.x + tx;

    if( kId + BSIZE < Nk ) {
        cache[ tx ] = A[ kId*Ni + i ]*B[ kId*Nj + j ];
        cache[ tx + BSIZE ] = A[ (kId + BSIZE)*Ni + i ]*B[ (kId + BSIZE)*Nj + j ];
    }
    else if( kId < Nk ) {
        cache[ tx ] = A[ kId*Ni + i ]*B[ kId*Nj + j ];
        cache[ tx + BSIZE ] = 0;
    }
    else {
        cache[ tx ] = 0;
        cache[ tx + BSIZE ] = 0;
    }

    int sz = BSIZE;

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

        __syncthreads();

    }
}
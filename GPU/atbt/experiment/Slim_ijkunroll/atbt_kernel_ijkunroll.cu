// matrix multiply kernel: C = A^T * B
__global__ void atbt(const double *A, const double *B, double *C, int Ni, int Nj, int Nk) {
    // Initially empty; will clearly not pass correctness
    
        int j = 2*blockIdx.y;
        int i = 2*blockIdx.x*blockDim.x+threadIdx.x;
        int remk = Nk%2;
        if ( ( i + blockDim.x < Ni ) && ( j + 1 < Nj ) ) { 
            double sum0 = 0;
            double sum1 = 0;
            double sum2 = 0;
            double sum3 = 0;

            for( int k = 0; k < remk; k++ ) {
                sum0 += A[ k*Ni + i ]*B[ j*Nk + k ]; 
                sum1 += A[ k*Ni + i + blockDim.x ]*B[ j*Nk + k ]; 
                sum2 += A[ k*Ni + i ]*B[ (j + 1)*Nk + k ]; 
                sum3 += A[ k*Ni + i + blockDim.x ]*B[ (j + 1)*Nk + k ];
            }
    
            for (int k = remk; k < Nk; k += 2) {
                sum0 += A[ k*Ni + i ]*B[ j*Nk + k ]; 
                sum1 += A[ k*Ni + i + blockDim.x ]*B[ j*Nk + k ]; 
                sum2 += A[ k*Ni + i ]*B[ (j + 1)*Nk + k ]; 
                sum3 += A[ k*Ni + i + blockDim.x ]*B[ (j + 1)*Nk + k ]; 

                sum0 += A[ (k + 1)*Ni + i ]*B[ j*Nk + k + 1 ]; 
                sum1 += A[ (k + 1)*Ni + i + blockDim.x ]*B[ j*Nk + k + 1 ]; 
                sum2 += A[ (k + 1)*Ni + i ]*B[ (j + 1)*Nk + k + 1 ]; 
                sum3 += A[ (k + 1)*Ni + i + blockDim.x ]*B[ (j + 1)*Nk + k + 1 ]; 
            }
                 
            C[ i * Nj + j ] = sum0;
            C[ (i + blockDim.x) * Nj + j ] = sum1;
            C[ i * Nj + j + 1 ] = sum2;
            C[ (i + blockDim.x) * Nj + j + 1 ] = sum3;
        }
        else if( i + blockDim.x < Ni && j < Nj ) {

            double sum0 = 0;
            double sum1 = 0;

            for( int k = 0; k < remk; k++ ) {
                sum0 += A[ k*Ni + i ]*B[ j*Nk + k ]; 
                sum1 += A[ k*Ni + i + blockDim.x ]*B[ j*Nk + k ]; 
            }
    
            for (int k = remk; k < Nk; k += 2) {
                 sum0 += A[ k*Ni + i ]*B[ j*Nk + k ]; 
                 sum1 += A[ k*Ni + i + blockDim.x ]*B[ j*Nk + k ]; 

                 sum0 += A[ (k + 1)*Ni + i ]*B[ j*Nk + (k + 1) ]; 
                 sum1 += A[ (k + 1)*Ni + i + blockDim.x ]*B[ j*Nk + (k + 1) ]; 
            }
                 
            C[ i * Nj + j ] = sum0;
            C[ (i + blockDim.x) * Nj + j ] = sum1;

        }
        else if( ( i < Ni ) && ( j + 1 < Nj ) ) {

            double sum0 = 0;
            double sum1 = 0;
    
            for( int k = 0; k < remk; k++ ) {
                sum0 += A[ k*Ni + i ]*B[ j*Nk + k ]; 
                sum1 += A[ k*Ni + i ]*B[ (j + 1)*Nk + k ]; 
            }

            for (int k = remk; k < Nk; k += 2) {
                 sum0 += A[ k*Ni + i ]*B[ j*Nk + k ]; 
                 sum1 += A[ k*Ni + i ]*B[ (j + 1)*Nk + k ]; 

                 sum0 += A[ (k + 1)*Ni + i ]*B[ j*Nk + k + 1 ]; 
                 sum1 += A[ (k + 1)*Ni + i ]*B[ (j + 1)*Nk + k + 1 ]; 
            }
                 
            C[ i * Nj + j ] = sum0;
            C[ i * Nj + j + 1 ] = sum1;

        }
        else if( i < Ni && j < Nj ) {

            double sum = 0;

            for( int k = 0; k < remk; k++ ) {
                sum += A[ k*Ni + i ]*B[ j*Nk + k ]; 
            }

            for (int k = remk; k < Nk; k += 2) {
                sum += A[ k*Ni + i ]*B[ j*Nk + k ]; 
                sum += A[ (k + 1)*Ni + i ]*B[ j*Nk + k + 1 ]; 
            }
                
            C[ i * Nj + j ] = sum;

        }
    
    }
    
    
// matrix multiply kernel: C = A^T * B
__global__ void atb(const double *A, const double *B, double *C, int Ni, int Nj, int Nk) {
    // Initially empty; will clearly not pass correctness
    
        int i = 2*( blockIdx.y*blockDim.y+threadIdx.y );
        int j = 2*blockIdx.x*blockDim.x+threadIdx.x;
        
        if ( ( i + 1 < Ni ) && ( j + blockDim.x < Nj ) ) { 
            double sum0 = 0;
            double sum1 = 0;
            double sum2 = 0;
            double sum3 = 0;

            int remk = Nk%2;

            for( int k = 0; k < remk; k++ ) {
                sum0 += A[ k*Ni + i ]*B[ k*Nj + j ]; 
                sum1 += A[ k*Ni + i + 1 ]*B[ k*Nj + j ]; 
                sum2 += A[ k*Ni + i ]*B[ k*Nj + j + blockDim.x ]; 
                sum3 += A[ k*Ni + i + 1 ]*B[ k*Nj + j + blockDim.x ]; 
            }
    
            for (int k = remk; k < Nk; k += 2) {
                sum0 += A[ k*Ni + i ]*B[ k*Nj + j ]; 
                sum1 += A[ k*Ni + i + 1 ]*B[ k*Nj + j ]; 
                sum2 += A[ k*Ni + i ]*B[ k*Nj + j + blockDim.x ]; 
                sum3 += A[ k*Ni + i + 1 ]*B[ k*Nj + j + blockDim.x ]; 

                sum0 += A[ ( k + 1 )*Ni + i ]*B[ (k + 1)*Nj + j ]; 
                sum1 += A[ ( k + 1 )*Ni + i + 1 ]*B[ (k + 1)*Nj + j ]; 
                sum2 += A[ ( k + 1 )*Ni + i ]*B[ (k + 1)*Nj + j + blockDim.x ]; 
                sum3 += A[ ( k + 1 )*Ni + i + 1 ]*B[ (k + 1)*Nj + j + blockDim.x ]; 
            }
                 
            C[ i * Nj + j ] = sum0;
            C[ (i + 1) * Nj + j ] = sum1;
            C[ i * Nj + j + blockDim.x ] = sum2;
            C[ (i + 1) * Nj + j + blockDim.x ] = sum3;
        }
        else if( ( i + 1 < Ni ) && ( j < Nj ) ) {

            double sum0 = 0;
            double sum1 = 0;

            int remk = Nk%2;

            for( int k = 0; k < remk; k++ ) {
                sum0 += A[ k*Ni + i ]*B[ k*Nj + j ]; 
                sum1 += A[ k*Ni + i + 1 ]*B[ k*Nj + j ]; 
            }
    
            for (int k = remk; k < Nk; k += 2) {
                sum0 += A[ k*Ni + i ]*B[ k*Nj + j ]; 
                sum1 += A[ k*Ni + i + 1 ]*B[ k*Nj + j ]; 

                sum0 += A[ ( k + 1 )*Ni + i ]*B[ (k + 1)*Nj + j ]; 
                sum1 += A[ ( k + 1 )*Ni + i + 1 ]*B[ (k + 1)*Nj + j ]; 
            }
                 
            C[ i * Nj + j ] = sum0;
            C[ (i + 1) * Nj + j ] = sum1;
        }
        else if( ( i < Ni ) && ( j + blockDim.x < Nj ) ) {

            double sum0 = 0;
            double sum1 = 0;
            
            int remk = Nk%2;

            for( int k = 0; k < remk; k++ ) {
                sum0 += A[ k*Ni + i ]*B[ k*Nj + j ]; 
                sum1 += A[ k*Ni + i ]*B[ k*Nj + j + blockDim.x ]; 
            }
    
            for (int k = remk; k < Nk; k += 2) {
                sum0 += A[ k*Ni + i ]*B[ k*Nj + j ]; 
                sum1 += A[ k*Ni + i ]*B[ k*Nj + j + blockDim.x ]; 

                sum0 += A[ ( k + 1 )*Ni + i ]*B[ (k + 1)*Nj + j ]; 
                sum1 += A[ ( k + 1 )*Ni + i ]*B[ (k + 1)*Nj + j + blockDim.x ]; 
            }
                 
            C[ i * Nj + j ] = sum0;
            C[ i * Nj + j + blockDim.x ] = sum1;

        }
        else if( ( i < Ni ) && ( j < Nj ) ) {

            double sum0 = 0;
            
            int remk = Nk%2;

            for( int k = 0; k < remk; k++ ) {
                sum0 += A[ k*Ni + i ]*B[ k*Nj + j ]; 
            }
    
            for (int k = remk; k < Nk; k += 2) {
                sum0 += A[ k*Ni + i ]*B[ k*Nj + j ]; 
                sum0 += A[ ( k + 1 )*Ni + i ]*B[ (k + 1)*Nj + j ]; 
            }
                 
            C[ i * Nj + j ] = sum0;

        }
    
    }
    
    
#include <stdio.h>
#include <time.h>
#define threshold 0.00001

#define FIXME1 1
#define FIXME2 2
#define FIXME3 3
#define FIXME4 4


void checkCUDAError(const char *msg);

cudaEvent_t start, stop;
float tstart, elapsedTime;

__global__ void atb(const int *A, const int *B, int *D, int*C, int Ni, int Nj, int Nk, int i, int j);

int main( int argc, char* argv[] ){

  int *h_A, *h_B, *h_C, *h_D, *h_Cref, *d_A, *d_B, *d_C, *d_D;
  int i,j,k,Ni,Nj,Nk;
  
  char* filename = argv[0];
  char* str2;
  str2 = ".txt";
  char * str3 = (char *) malloc(1 + strlen(filename)+ strlen(str2) );
  strcpy(str3, filename);
  strcat(str3, str2);

  FILE* fptr;
  fptr = fopen( str3, "w" );

  Ni = 4;
  Nj = 4;
  Nk = 8192;

  printf("Specify Matrix dimension Ni, Nj, Nk: %d %d %d", Ni, Nj, Nk);
  fprintf( fptr, "Specify Matrix dimension Ni, Nj, Nk: %d %d %d", Ni, Nj, Nk);
  
  // scanf("%d %d %d", &Ni,&Nj,&Nk);
  h_A = (int *) malloc(sizeof(int)*Ni*Nk);
  h_B = (int *) malloc(sizeof(int)*Nk*Nj);
  h_C = (int *) malloc(sizeof(int)*Ni*Nj);
  h_Cref = (int *) malloc(sizeof(int)*Ni*Nj);
  for (i=0; i<Ni; i++)
   for (k=0; k<Nk; k++)
    h_A[k*Ni+i] = 2;//i*Ni+i-1;
  for (k=0; k<Nk; k++)
   for (j=0; j<Nj; j++)
    h_B[k*Nj+j] = 2890;//j*Nj+j+1;
  for (i=0; i<Ni; i++)
   for (j=0; j<Nj; j++) {
    h_C[i*Nj+j] = 0;
    h_Cref[i*Nj+j] = 0;}

  for (i=0;i<Ni;i++)
   for (k=0;k<Nk;k++)
    for (j=0;j<Nj;j++)
  // h_Cref[i][j] += h_A[k][i]*h_B[k][j];
     h_Cref[i*Nj+j] += h_A[i+Ni*k]*h_B[k*Nj+j];
  // Allocate device memory and copy input data over to GPU
  cudaMalloc(&d_A, Nk*Ni*sizeof(int));
  cudaMalloc(&d_B, Nj*Nk*sizeof(int));
  cudaMalloc(&d_C, Ni*Nj*sizeof(int));
  checkCUDAError("cudaMalloc failure");
  cudaMemcpy(d_A, h_A, Nk*Ni*sizeof(int), cudaMemcpyHostToDevice);
  cudaMemcpy(d_B, h_B, Nj*Nk*sizeof(int), cudaMemcpyHostToDevice);
  checkCUDAError("cudaMemcpy H2D failure");

  int BSIZE = 64;
  int bsx = BSIZE;
  int bsy = 1;
  int gridsz = ceil( Nk/( (double)bsx ) );

  dim3 block(bsx, bsy);
  dim3 grid( gridsz, 1 );

    h_D = (int *) malloc(sizeof(int)*gridsz);
    for( int i = 0; i < gridsz; i++ ) {
        h_D[i] = 0;
    }
  cudaMalloc(&d_D, gridsz*sizeof(int));
  cudaMemcpy(d_D, h_D, gridsz*sizeof(int), cudaMemcpyHostToDevice);

  cudaEventCreate(&start);
  cudaEventCreate(&stop);
  for(int trial=0;trial<5;trial++)
  {

   cudaEventRecord(start);
   // Launch kernel

    for( int i = 0; i < Ni; i++ ) {
        for( int j = 0; j < Nj; j++ ) {

            

            atb<<<grid, block>>>(d_A, d_B, d_D, d_C, Ni,Nj,Nk, i, j);
            // cudaMemcpy(h_D, d_D, gridsz*sizeof(double), cudaMemcpyDeviceToHost);

            // double sum = 0;
            // for( int k = 0; k < gridsz; k++ ) {
            //     sum += h_D[k];
            // }

            // h_C[i*Nj + j] = sum;
        }
    }


   cudaEventRecord(stop);
   checkCUDAError("kernel launch");
   cudaEventSynchronize(stop);
   cudaEventElapsedTime(&elapsedTime, start,stop);
//   cudaDeviceSynchronize();
   // Copy results back to host
   cudaMemcpy(h_C, d_C, Ni*Nj*sizeof(int), cudaMemcpyDeviceToHost);
   checkCUDAError("cudaMemcpy D2H");
   for (int l = 0; l < Ni*Nj; l++) if (fabs((h_C[l] - h_Cref[l])/h_Cref[l])>threshold) {
    printf("Error: mismatch at linearized index %d, was: %f, should be: %f\n", l, h_C[l], h_Cref[l]); 
    fprintf(fptr, "Error: mismatch at linearized index %d, was: %f, should be: %f\n", l, h_C[l], h_Cref[l]);
    // return -1;
  }
   printf("<Ni=%d,Nj=%d,Nk=%d>: Trial %d: GFLOPS: %.2f\n",Ni,Nj,Nk,trial,2.0e-6*Ni*Nj*Nk/elapsedTime);
  }
}

void checkCUDAError(const char *msg)
{
    cudaError_t err = cudaGetLastError();
    if( cudaSuccess != err)
    {
        fprintf(stderr, "Cuda error: %s: %s.\n", msg, cudaGetErrorString( err) );
        exit(EXIT_FAILURE);
    }
}


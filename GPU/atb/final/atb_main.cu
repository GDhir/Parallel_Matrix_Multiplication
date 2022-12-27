#include <stdio.h>
#include <time.h>
#define threshold 0.0000001

#define FIXME1 1
#define FIXME2 2
#define FIXME3 3
#define FIXME4 4


void checkCUDAError(const char *msg);

cudaEvent_t start, stop;
float tstart, elapsedTime;

__global__ void atb(const double *A, const double *B, double *C, int Ni, int Nj, int Nk);
__global__ void atb_kernel_dot(const double *A, const double *B, double *D, double* C, int Ni, int Nj, int Nk, int i, int j);
__global__ void atb_kernel_tiled32(const double *A, const double *B, double *C, int Ni, int Nj, int Nk);
__global__ void atb_kernel_tiled_ijkunroll32(const double *A, const double *B, double *C, int Ni, int Nj, int Nk);
__global__ void atb_kernel_tiled_kunroll32(const double *A, const double *B, double *C, int Ni, int Nj, int Nk);
__global__ void atb_kernel_tiled32(const double *A, const double *B, double *C, int Ni, int Nj, int Nk);
__global__ void atb_kernel_slim(const double *A, const double *B, double *C, int Ni, int Nj, int Nk);
__global__ void atb_kernel_slimSM_kunroll8(const double *A, const double *B, double *C, int Ni, int Nj, int Nk);
__global__ void atb_kernel_slimSM_kunroll16(const double *A, const double *B, double *C, int Ni, int Nj, int Nk);
__global__ void atb_kernel_slimSM_kunroll32(const double *A, const double *B, double *C, int Ni, int Nj, int Nk);
__global__ void atb_kernel_tiled_ijkunroll16(const double *A, const double *B, double *C, int Ni, int Nj, int Nk);

void atb_seq(const double *__restrict__ A, const double *__restrict__ B, double *__restrict__ C, int Ni, int Nj, int Nk)
{
  int i, j, k;

  for (i = 0; i < Ni; i++)
   for (j = 0; j < Nj; j++)
    for (k = 0; k < Nk; k++)
// C[i][j] = C[i][j] + A[k][i]*B[k][j];
     C[i*Nj+j]=C[i*Nj+j]+A[k*Ni+i]*B[k*Nj+j];
}

void launcher( const double *d_A, const double *d_B, double *d_C, int Ni, int Nj, int Nk ) {

  int BSIZE = 32;
  int bsx;
  int bsy;

  if( Ni >= 180 && Nj >= 180 ) {

    if( Nk >= 32 ) {

      printf( "t1" );
      bsx = BSIZE;
      bsy = BSIZE;

      dim3 block(bsx, bsy);
      dim3 grid( ceil( Nj/( (double)2*bsx ) ), ceil( Ni/( (double)2*bsy ) ) );
      atb_kernel_tiled_ijkunroll32<<<grid, block>>>(d_A, d_B, d_C, Ni, Nj, Nk);
    
    }
    else if( Nk >= 16 ) {
      printf( "t2" );
      bsx = BSIZE/2;
      bsy = BSIZE/2;

      dim3 block(bsx, bsy);
      dim3 grid( ceil( Nj/( (double)2*bsx ) ), ceil( Ni/( (double)2*bsy ) ) );
      atb_kernel_tiled_ijkunroll16<<<grid, block>>>(d_A, d_B, d_C, Ni,Nj,Nk);

    }
    else if( Nk >= 8 ) {
      printf( "t3" );
      bsx = BSIZE/4;
      bsy = 1;

      dim3 block(bsx, bsy);
      dim3 grid( ceil( Nj/( (double)bsx ) ), ceil( Ni/( (double)bsy ) ) );
      atb_kernel_slimSM_kunroll8<<<grid, block>>>(d_A, d_B, d_C, Ni,Nj,Nk);

    }
    else {

      bsx = 4;
      bsy = 1;

      dim3 block(bsx, bsy);
      dim3 grid( ceil( Nj/( (double)bsx ) ), ceil( Ni/( (double)bsy ) ) );
      atb_kernel_slim<<<grid, block>>>(d_A, d_B, d_C, Ni,Nj,Nk);

    }

  }
  else if( Ni > 128 && Nk > 128 ) {

    if( Nj >= 512 ) {
      printf( "t4" );
      bsx = BSIZE;
      bsy = BSIZE;

      dim3 block(bsx, bsy);
      dim3 grid( ceil( Nj/( (double)2*bsx ) ), ceil( Ni/( (double)2*bsy ) ) );
      atb_kernel_tiled_ijkunroll32<<<grid, block>>>(d_A, d_B, d_C, Ni, Nj, Nk);

    }
    else {
      printf( "t5" );
      bsx = BSIZE;
      bsy = 1;

      dim3 block(bsx, bsy);
      dim3 grid( ceil( Nj/( (double)bsx ) ), ceil( Ni/( (double)bsy ) ) );
      atb_kernel_slimSM_kunroll32<<<grid, block>>>(d_A, d_B, d_C, Ni,Nj,Nk);

    }

  }
  else if( Ni > 64 && Nj > 256 && Nk < 1024 ) {
    printf( "t6" );
    bsx = BSIZE;
    bsy = BSIZE;

    dim3 block(bsx, bsy);
    dim3 grid( ceil( Nj/( (double)bsx ) ), ceil( Ni/( (double)bsy ) ) );
    atb_kernel_tiled_kunroll32<<<grid, block>>>(d_A, d_B, d_C, Ni, Nj, Nk);

  }
  else if( Ni > 64 && Nj > 1024 ) {

    if( Nk < 2048 ) {
      printf( "t7" );

      if( Nk > 16 ) {

        bsx = BSIZE;
        bsy = BSIZE;

        dim3 block(bsx, bsy);
        dim3 grid( ceil( Nj/( (double)bsx ) ), ceil( Ni/( (double)bsy ) ) );
        atb_kernel_tiled_kunroll32<<<grid, block>>>(d_A, d_B, d_C, Ni, Nj, Nk);
      }
      else {

        bsx = BSIZE/4;
        bsy = 1;

        dim3 block(bsx, bsy);
        dim3 grid( ceil( Nj/( (double)bsx ) ), ceil( Ni/( (double)bsy ) ) );
        atb_kernel_slimSM_kunroll8<<<grid, block>>>(d_A, d_B, d_C, Ni,Nj,Nk);

      }
    }
    else {
      printf( "t8" );
      bsx = BSIZE;
      bsy = BSIZE;

      dim3 block(bsx, bsy);
      dim3 grid( ceil( Nj/( (double)2*bsx ) ), ceil( Ni/( (double)2*bsy ) ) );
      atb_kernel_tiled_ijkunroll32<<<grid, block>>>(d_A, d_B, d_C, Ni, Nj, Nk);

    }

  }
  else if( Ni > 32 && Nj > 1024 && Nk > 1024 ) {
    printf( "t10" );
      bsx = BSIZE;
      bsy = BSIZE;

      dim3 block(bsx, bsy);
      dim3 grid( ceil( Nj/( (double)2*bsx ) ), ceil( Ni/( (double)2*bsy ) ) );
      atb_kernel_tiled_ijkunroll32<<<grid, block>>>(d_A, d_B, d_C, Ni, Nj, Nk);

  }
  else if( Ni > 32 && Nj > 4096 && Nk > 256 ) {
    printf( "t11" );
      bsx = BSIZE;
      bsy = BSIZE;

      dim3 block(bsx, bsy);
      dim3 grid( ceil( Nj/( (double)bsx ) ), ceil( Ni/( (double)bsy ) ) );
      atb_kernel_tiled_kunroll32<<<grid, block>>>(d_A, d_B, d_C, Ni, Nj, Nk);

  }
  else {

    if( Ni >= 64 && Nj >= 256 && Nk >= 32 ) {
      bsx = BSIZE;
      bsy = 1;

      dim3 block(bsx, bsy);
      dim3 grid( ceil( Nj/( (double)bsx ) ), ceil( Ni/( (double)bsy ) ) );
      atb_kernel_slimSM_kunroll32<<<grid, block>>>(d_A, d_B, d_C, Ni,Nj,Nk);
    }
    else if( Ni < 4 || Nj < 4  ) {
      bsx = 2;
      bsy = 1;

      dim3 block(bsx, bsy);
      dim3 grid( ceil( Nj/( (double)bsx ) ), ceil( Ni/( (double)bsy ) ) );
      atb_kernel_slim<<<grid, block>>>(d_A, d_B, d_C, Ni,Nj,Nk);
    }
    else if( Ni < 8 || Nj < 8 ) {

      bsx = 4;
      bsy = 1;

      dim3 block(bsx, bsy);
      dim3 grid( ceil( Nj/( (double)bsx ) ), ceil( Ni/( (double)bsy ) ) );
      atb_kernel_slim<<<grid, block>>>(d_A, d_B, d_C, Ni,Nj,Nk);

    }
    else if( Ni < 16 || Nj < 16  ) {

      bsx = 8;
      bsy = 1;

      dim3 block(bsx, bsy);
      dim3 grid( ceil( Nj/( (double)bsx ) ), ceil( Ni/( (double)bsy ) ) );
      atb_kernel_slimSM_kunroll8<<<grid, block>>>(d_A, d_B, d_C, Ni,Nj,Nk);

    }
    else if( Ni < 32 || Nj < 32  ) {

      bsx = 16;
      bsy = 1;

      dim3 block(bsx, bsy);
      dim3 grid( ceil( Nj/( (double)bsx ) ), ceil( Ni/( (double)bsy ) ) );
      atb_kernel_slimSM_kunroll16<<<grid, block>>>(d_A, d_B, d_C, Ni,Nj,Nk);

    }
    else {

      bsx = 32;
      bsy = 1;

      dim3 block(bsx, bsy);
      dim3 grid( ceil( Nj/( (double)bsx ) ), ceil( Ni/( (double)bsy ) ) );
      atb_kernel_slimSM_kunroll32<<<grid, block>>>(d_A, d_B, d_C, Ni,Nj,Nk);

    }

}


};

int main(){

  double *h_A, *h_B, *h_C, *h_Cref, *d_A, *d_B, *d_C;
  int i,j,k,Ni,Nj,Nk;

  printf("Specify Matrix dimension Ni, Nj, Nk: ");
  scanf("%d %d %d", &Ni,&Nj,&Nk);
  h_A = (double *) malloc(sizeof(double)*Ni*Nk);
  h_B = (double *) malloc(sizeof(double)*Nk*Nj);
  h_C = (double *) malloc(sizeof(double)*Ni*Nj);
  h_Cref = (double *) malloc(sizeof(double)*Ni*Nj);
  for (i=0; i<Ni; i++)
   for (k=0; k<Nk; k++)
    h_A[k*Ni+i] = k*Ni+i-1;
  for (k=0; k<Nk; k++)
   for (j=0; j<Nj; j++)
    h_B[k*Nj+j] = k*Nj+j+1;
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
  cudaMalloc(&d_A, Nk*Ni*sizeof(double));
  cudaMalloc(&d_B, Nj*Nk*sizeof(double));
  cudaMalloc(&d_C, Ni*Nj*sizeof(double));
  checkCUDAError("cudaMalloc failure");
  cudaMemcpy(d_A, h_A, Nk*Ni*sizeof(double), cudaMemcpyHostToDevice);
  cudaMemcpy(d_B, h_B, Nj*Nk*sizeof(double), cudaMemcpyHostToDevice);
  checkCUDAError("cudaMemcpy H2D failure");

  // dim3 block(FIXME1,FIXME2);
  // dim3 grid(FIXME3,FIXME4);

  cudaEventCreate(&start);
  cudaEventCreate(&stop);
  for(int trial=0;trial<5;trial++)
  {
   cudaEventRecord(start);
   // Launch kernel
   launcher(d_A, d_B, d_C, Ni, Nj, Nk);
  //  atb<<<grid, block>>>(d_A, d_B, d_C, Ni,Nj,Nk);
   cudaEventRecord(stop);
   checkCUDAError("kernel launch");
   cudaEventSynchronize(stop);
   cudaEventElapsedTime(&elapsedTime, start,stop);
//   cudaDeviceSynchronize();
   // Copy results back to host
   cudaMemcpy(h_C, d_C, Ni*Nj*sizeof(double), cudaMemcpyDeviceToHost);
   checkCUDAError("cudaMemcpy D2H");
   for (int l = 0; l < Ni*Nj; l++) if (fabs((h_C[l] - h_Cref[l])/h_Cref[l])>threshold) {printf("Error: mismatch at linearized index %d, was: %f, should be: %f\n", l, h_C[l], h_Cref[l]); return -1;}
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



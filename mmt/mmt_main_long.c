// Use "gcc -O3 -fopenmp" to compile

#include <math.h>
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <unistd.h>
#ifndef ni
#define ni (16)
#endif
#ifndef nj
#define nj (16)
#endif
#ifndef nk
#define nk (1024*1024)
#endif
#define NTrials (5)
#define threshold (0.0000001)

void compare(FILE* fptr, int Ni, int Nj, long double wref[Ni][Nj], long double w[Ni][Nj], int numt);
void mmt_seq(int Ni, int Nj, int Nk, long double *__restrict__ a, long double *__restrict__ b,
                long double *__restrict__ c);
void mmt_par(int Ni, int Nj, int Nk, long double *__restrict__ a, long double *__restrict__ b,
                long double *__restrict__ c);

long double c[ni][nj], b[nk][nj], a[nk][ni], cc[ni][nj];

int main(int argc, char *argv[]){
  double tstart,telapsed;
  FILE* fptr;
  fptr = fopen( "mmtperf.txt", "w" );

  int i,j,k,nt,trial,max_threads,num_cases;
  int nthr_32[9] = {1,2,4,8,10,12,14,15,31};
  int nthr_40[9] = {1,2,4,8,10,12,14,19,39};
  int nthr_48[9] = {1,2,4,8,10,15,20,23,47};
  int nthr_56[9] = {1,2,4,8,10,15,20,27,55};
  int nthr_64[9] = {1,2,4,8,10,15,20,31,63};
  int nthreads[9];
  double mint_par[9],maxt_par[9];
  double mint_seq,maxt_seq;
  
  printf("Matrix Size = %d x %d; NTrials=%d\n", ni, nj, NTrials);
  fprintf(fptr, "Matrix Size = %d x %d; NTrials=%d\n", ni, nj,NTrials);
  
  for(k=0;k<nk;k++)
   for(i=0;i<ni;i++)
   { a[k][i] = 1.1*(2*k+i);
  //  a[k][i] = (2*k+i);
   }

  for(k=0;k<nk;k++)
   for(j=0;j<nj;j++)
   { 
      b[k][j] = 1.2*(k+2*j);
      // b[k][j] = (k+2*j);
   }

  printf("Reference sequential code performance in GFLOPS");
  fprintf(fptr, "Reference sequential code performance in GFLOPS");
  mint_seq = 1e9; maxt_seq = 0;
  for(trial=0;trial<NTrials;trial++)
  {
   for(i=0;i<ni;i++) for(j=0;j<nj;j++) c[i][j] = 0;
   tstart = omp_get_wtime();
   mmt_seq(ni, nj, nk, &a[0][0], &b[0][0], &c[0][0]);
   telapsed = omp_get_wtime()-tstart;
   if (telapsed < mint_seq) mint_seq=telapsed;
   if (telapsed > maxt_seq) maxt_seq=telapsed;
  }
  printf(" Min: %.2f; Max: %.2f\n",2.0e-9*ni*nj*nk/maxt_seq,2.0e-9*ni*nj*nk/mint_seq);
  fprintf(fptr, " Min: %.2f; Max: %.2f\n",2.0e-9*ni*nj*nk/maxt_seq,2.0e-9*ni*nj*nk/mint_seq);

  max_threads = omp_get_max_threads();
  printf("Max Threads (from omp_get_max_threads) = %d\n",max_threads);
  fprintf(fptr, "Max Threads (from omp_get_max_threads) = %d\n",max_threads);
  switch (max_threads)
  {
	  case 32: for(i=0;i<9;i++) nthreads[i] = nthr_32[i]; num_cases=9; break;
	  case 40: for(i=0;i<9;i++) nthreads[i] = nthr_40[i]; num_cases=9; break;
	  case 48: for(i=0;i<9;i++) nthreads[i] = nthr_48[i]; num_cases=9; break;
	  case 56: for(i=0;i<9;i++) nthreads[i] = nthr_56[i]; num_cases=9; break;
	  case 64: for(i=0;i<9;i++) nthreads[i] = nthr_64[i]; num_cases=9; break;
	  default: {
                    nt = 1;i=0;
                    while (nt <= max_threads) {nthreads[i]=nt; i++; nt *=2;}
                    if (nthreads[i-1] < max_threads) {nthreads[i] = max_threads; i++;}
                    num_cases = i;
                    nthreads[num_cases-1]--;
                    nthreads[num_cases-2]--;
		   }
  }

  for (nt=0;nt<num_cases;nt ++)
  {
   omp_set_num_threads(nthreads[nt]);
   mint_par[nt] = 1e9; maxt_par[nt] = 0;
   for (trial=0;trial<NTrials;trial++)
   {
    for(i=0;i<ni;i++) for(j=0;j<nj;j++) cc[i][j] = 0;
    tstart = omp_get_wtime();
    mmt_par(ni,nj,nk, &a[0][0], &b[0][0], &cc[0][0]);
    telapsed = omp_get_wtime()-tstart;
    if (telapsed < mint_par[nt]) mint_par[nt]=telapsed;
    if (telapsed > maxt_par[nt]) maxt_par[nt]=telapsed;
   }
   compare(fptr, ni, nj, c,cc,nthreads[nt]);
  }
  printf("Performance (Best & Worst) of parallelized version: GFLOPS on ");
  fprintf(fptr, "Performance (Best & Worst) of parallelized version: GFLOPS on ");
  for (nt=0;nt<num_cases-1;nt++) {printf("%d/",nthreads[nt]); fprintf(fptr, "%d/",nthreads[nt]); }
  printf("%d threads\n",nthreads[num_cases-1]);
  printf("Best Performance (GFLOPS): ");
  fprintf(fptr, "%d threads\n",nthreads[num_cases-1]);
  fprintf(fptr, "Best Performance (GFLOPS): ");
  for (nt=0;nt<num_cases;nt++) {printf("%.2f ",2.0e-9*ni*nj*nk/mint_par[nt]);fprintf(fptr, "%.2f ",2.0e-9*ni*nj*nk/mint_par[nt]);}
  printf("\n");
  printf("Worst Performance (GFLOPS): ");
  fprintf(fptr, "\n");
  fprintf(fptr, "Worst Performance (GFLOPS): ");
  for (nt=0;nt<num_cases;nt++) {printf("%.2f ",2.0e-9*ni*nj*nk/maxt_par[nt]);fprintf(fptr, "%.2f ",2.0e-9*ni*nj*nk/maxt_par[nt]);}
  printf("\n");
  fprintf(fptr, "\n");
  fclose(fptr);
}

void compare(FILE* fptr, int Ni, int Nj, long double wref[Ni][Nj], long double w[Ni][Nj], int numt)
{
  long double maxdiff,this_diff;
  int numdiffs;
  int i,j;
  numdiffs = 0;
  maxdiff = 0;
  for (i=0;i<ni;i++)
    for (j=0;j<nj;j++)
      {
        this_diff = wref[i][j]-w[i][j];
        if (this_diff < 0) this_diff = -1.0*this_diff;
        if (this_diff>threshold)
          { numdiffs++;
            if (this_diff > maxdiff) maxdiff=this_diff;
          }

          // printf( "%f \t %f \n", w[i][j], wref[i][j] );
      }
  if (numdiffs > 0)
  { printf("Error when executing on %d threads; %d Differences found over threshold %f; Max Diff = %f\n",
           numt,numdiffs,threshold,maxdiff);
    printf("Exiting\n");
    fprintf(fptr, "Error when executing on %d threads; %d Differences found over threshold %f; Max Diff = %f\n",
           numt,numdiffs,threshold,maxdiff);
    fprintf(fptr, "Exiting\n"); exit;
  
  }
}

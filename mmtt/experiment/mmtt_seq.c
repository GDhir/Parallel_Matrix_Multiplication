void mmtt_seq(int ni, int nj, int nk, double *__restrict__ a, double *__restrict__ b,
                 double *__restrict__ c) {
int i, j, k;

for(i=0;i<ni;i++)
 for(j=0;j<nj;j++)
  for(k=0;k<nk;k++)
//    c[i][j] = c[i][j] + a[k][i]*b[j][k];
   c[i * nj + j] = c[i * nj + j] + a[k * ni + i] * b[j * nk + k];
}

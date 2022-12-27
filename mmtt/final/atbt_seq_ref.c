void atbt_seq_ref(const double *__restrict__ a, const double *__restrict__ b, double *__restrict__ c, int ni, int nj, int nk) {
int i, j, k;

for(i=0;i<ni;i++)
 for(j=0;j<nj;j++)
  for(k=0;k<nk;k++)
//    c[i][j] = c[i][j] + a[k][i]*b[j][k];
   c[i * nj + j] = c[i * nj + j] + a[k * ni + i] * b[j * nk + k];
}

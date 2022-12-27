void mmt_seq(int ni, int nj, int nk, long double *__restrict__ a, long double *__restrict__ b,
                 long double *__restrict__ c) {
int i, j, k;
long double y, t, temp;

long double cval[ni*nj];

for( int i = 0; i < ni*nj; i++ ) {
   cval[i] = 0;
}

for(i=0;i<ni;i++) {
 for(j=0;j<nj;j++) {
  for(k=0;k<nk;k++) {
//    c[i][j] = c[i][j] + a[k][j]*b[k][i];
   temp = a[k*ni+i]*b[k*nj+j];
   y = temp - cval[i*nj + j];
   t = c[i*nj + j] + y;

   cval[i*nj + j] = ( t - c[i * nj + j] ) - y;
   // c[i*nj+j] += a[k*ni+i]*b[k*nj+j];
   c[i*nj + j] = t;
  }
 }
}
}

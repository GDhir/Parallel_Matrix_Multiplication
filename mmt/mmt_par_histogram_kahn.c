#include "stdio.h"

void mmt_par(int ni, int nj, int nk, long double *__restrict__ a, long double *__restrict__ b, long double *__restrict__ c)
{
    int i, j, k;
    long double cshared[ni*nj];

    for (i = 0; i < ni * nj; i++)
    {
        cshared[i] = 0;
    }

#pragma omp parallel private(i, j, k)
    {
        long double cprivate[ni * nj];
        long double cval[ni * nj];
        long double y, t, temp;

        for (i = 0; i < ni * nj; i++)
        {
            cprivate[i] = 0;
            cval[i] = 0;
        }

        #pragma omp for schedule(static)
        for (k = 0; k < nk; k++)
        {
            for (i = 0; i < ni; i++)
            {
                for (j = 0; j < nj; j++)
                {

                    temp = a[k * ni + i] * b[k * nj + j];
                    y = temp - cval[i * nj + j];
                    t = cprivate[i * nj + j] + y;
                    cval[i * nj + j] = (t - cprivate[i * nj + j]) - y;
                    cprivate[i * nj + j] = t;
                }
            }
        }

        
        #pragma omp critical
        {
        for (i = 0; i < ni; i++)
        {
            for (j = 0; j < nj; j++)
            {

                    y = cprivate[i * nj + j] - cshared[i * nj + j];
                    t = c[i * nj + j] + y;
                    cshared[i * nj + j] = ( t - c[i * nj + j] ) - y;
                    c[i * nj + j] = t;
                
            }
        }
        }
    }
}

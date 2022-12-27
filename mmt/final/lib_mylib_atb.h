#ifndef libval
#define libval 10

void atb_par_histogram(const double *__restrict__ a,const double *__restrict__ b, double *__restrict__ c, int Ni, int Nj, int Nk);
void atb_par_iktiled(const double *__restrict__ a, const double *__restrict__ b, double *__restrict__ c, int ni, int nj, int nk);
void atb_par_iktiled_kunroll2(const double *__restrict__ a, const double *__restrict__ b, double *__restrict__ c, int ni, int nj, int nk);
void atb_par_itiled_kunroll2(const double *__restrict__ a, const double *__restrict__ b, double *__restrict__ c, int ni, int nj, int nk);
void atb_par_itiled( const double *__restrict__ a, const double *__restrict__ b, double *__restrict__ c, int ni, int nj, int nk);
void atb_par_jtiled( const double *__restrict__ a, const double *__restrict__ b, double *__restrict__ c, int ni, int nj, int nk);
void atb_par_kunroll4(const double *__restrict__ a, const double *__restrict__ b, double *__restrict__ c, int ni, int nj, int nk);
void atb_par_ref(const double *__restrict__ a, const double *__restrict__ b, double *__restrict__ c, int ni, int nj, int nk);
void atb_par_jtiled_kunroll2( const double *__restrict__ a, const double *__restrict__ b, double *__restrict__ c, int ni, int nj, int nk);

#endif
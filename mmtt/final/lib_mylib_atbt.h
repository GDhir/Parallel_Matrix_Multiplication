#ifndef libval
#define libval 30

void atbt_par_kunroll8(const double *__restrict__ a, const double *__restrict__ b, double *__restrict__ c, int ni, int nj, int nk);
void atbt_par_kunroll4(const double *__restrict__ a, const double *__restrict__ b, double *__restrict__ c, int ni, int nj, int nk);
void atbt_par_kunroll2(const double *__restrict__ a, const double *__restrict__ b, double *__restrict__ c, int ni, int nj, int nk);
void atbt_par_ktile_kunroll8(const double *__restrict__ a, const double *__restrict__ b, double *__restrict__ c, int ni, int nj, int nk);
void atbt_par_kjtile8_kunroll4(const double *__restrict__ a, const double *__restrict__ b, double *__restrict__ c, int ni, int nj, int nk);
void atbt_par_jtile16_kunroll8_ikj(const double *__restrict__ a, const double *__restrict__ b, double *__restrict__ c, int ni, int nj, int nk);
void atbt_par_ikjtile16_kunroll8(const double *__restrict__ A, const double *__restrict__ B, double *__restrict__ C, int ni, int nj, int nk);
void atbt_par_ikjtile16_kunroll4(const double *__restrict__ A, const double *__restrict__ B, double *__restrict__ C, int ni, int nj, int nk);
void atbt_par_ijtile16_kunroll4_kij(const double *__restrict__ A, const double *__restrict__ B, double *__restrict__ C, int Ni, int Nj, int Nk);
void atbt_par_histogram(const double *__restrict__ A, const double *__restrict__ B, double *__restrict__ C, int Ni, int Nj, int Nk);
void atbt_seq_ref(const double *__restrict__ a, const double *__restrict__ b, double *__restrict__ c, int ni, int nj, int nk);
void atbt_par_ref( const double *__restrict__ a, const double *__restrict__ b, double *__restrict__ c, int ni, int nj, int nk );

#endif
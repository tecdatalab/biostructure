#ifndef __SITUS_LIB_DDF
#define __SITUS_LIB_DDF

#ifdef __cplusplus
extern "C" {
#endif

/* header file for lib_ddf.c */
int compfunc (const void *a, const void *b);
void indexy (unsigned long n, double *arr, unsigned long *indx);
void hist_eq (double *p, unsigned long M);
void hist_eq_2 (double *p, double *q, unsigned long M);
void rotmat (double, double [3], double [3][3]);
void matvec (double [3][3], double [3], double [3], double [3]);
unsigned long bisearch (double *a, unsigned long n, double x);
int linsolver (double *A, int n, int *perm, double *par, double *b);
double J3(double *Qs, int m, int M, double b, double c, double k);

#ifdef __cplusplus
}
#endif

#endif

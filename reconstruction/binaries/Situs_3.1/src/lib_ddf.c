/*********************************************************************
*                           L I B _ D D F                            *
**********************************************************************
* Library is part of the Situs package URL: situs.biomachina.org     *
* (c) Julio Kovacs and Willy Wriggers, 2018                          *
**********************************************************************
*                                                                    *
* Misc routines currently used only by DDforge                       *
*                                                                    *
**********************************************************************
* See legal statement for terms of distribution                      *
*********************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "lib_ddf.h"

typedef struct {
  unsigned long index;
  double value;
} intdou;


/*===================================================================================*/
int compfunc (const void *a, const void *b) {
  /* comparison function to use with qsort, in order to obtain an index table
     by means of indexy  */
  if ( ((intdou*)a)->value <  ((intdou*)b)->value ) return -1;
  if ( ((intdou*)a)->value == ((intdou*)b)->value ) return  0;
  if ( ((intdou*)a)->value >  ((intdou*)b)->value ) return  1;
}


/*===================================================================================*/
void indexy(unsigned long n, double *arr, unsigned long *indx) {
  /* makes an index table for the array arr, using the standard library function qsort */

  unsigned long i;
  intdou *wksp;

  wksp  = (intdou *) malloc(n*sizeof(intdou));

  for (i=0; i<n; i++) {
    wksp[i].index = i;
    wksp[i].value = arr[i];
  }

  qsort(wksp, n, sizeof(intdou), compfunc);

  for (i=0; i<n; i++)
    indx[i] = wksp[i].index;

  free(wksp);
}


/*===================================================================================*/
void hist_eq(double *signal, unsigned long M) {
  /* replaces 'signal' by its histogram-equalized version */

  double Md, *wksp;
  unsigned long j, m, *iwksp;

  wksp  = (double *) malloc((M+1)*sizeof(double));
  iwksp = (unsigned long *) malloc((M+1)*sizeof(long));

  Md = (double) M;

  /* construct index table (permutation) for signal */
  indexy(M, signal, iwksp);

  /* sort signal using the index */
  for (j=0; j<M; j++)
    wksp[j] = signal[iwksp[j]];

  /* fill signal with the hist-eq signal */
  /* this would be the array of ranks, were it not by the possible
     repetition of values, dealt with in the while loop */
  for (m=0,j=0; m<M; m++) {
    while (j<M && wksp[j]<=wksp[m]) j++;
    signal[iwksp[m]] = j/Md;
  }

  free(wksp); free(iwksp);
}


/*===================================================================================*/
void hist_eq_2(double *signal, double *hesorted, unsigned long M) {
  /* replaces 'signal' by its histogram-equalized version,
     and puts the sorted hist-eq'ed signal in 'hesorted' */

  double Md, *wksp;
  unsigned long j, m, *iwksp;

  wksp  = (double *) malloc((M+1)*sizeof(double));
  iwksp = (unsigned long *) malloc((M+1)*sizeof(long));

  Md = (double) M;

  /* construct index table (permutation) for signal */
  indexy(M, signal, iwksp);

  /* sort signal using the index */
  for (j=0; j<M; j++)
    wksp[j] = signal[iwksp[j]];

  /* fill 'signal' with the hist-eq signal, and 'hesorted'
     with the sorted version */
  /* this would be the array of ranks, were it not by the possible
     repetition of values, dealt with in the while loop */
  for (m=0,j=0; m<M; m++) {
    while (j<M && wksp[j]<=wksp[m]) j++;
    hesorted[m] = signal[iwksp[m]] = j/Md;
  }

  free(wksp); free(iwksp);
}


/*==================================================*/
void rotmat(double a, double v[3], double rot[3][3]) {
  /* computes the rotation matrix corresponding to a
   rotation angle a (in radians) around the vector v */

  double sa, ca, uca, t, v1, v2, v3;
  double v11, v12, v13, v22, v23, v33;
  int m;

  t=sqrt(pow(v[0],2)+pow(v[1],2)+pow(v[2],2));
  for(m=0;m<3;m++)     /* normalize v */
    v[m] /= t;

  v1=v[0]; v2=v[1]; v3=v[2];
  v11=v1*v1; v12=v1*v2; v13=v1*v3;
  v22=v2*v2; v23=v2*v3; v33=v3*v3;
  ca=cos(a); sa=sin(a); uca=1.0-ca;

  rot[0][0] = ca+v11*uca;
  rot[0][1] = v12*uca-v3*sa;
  rot[0][2] = v13*uca+v2*sa;
  rot[1][0] = v12*uca+v3*sa;
  rot[1][1] = ca+v22*uca;
  rot[1][2] = v23*uca-v1*sa;
  rot[2][0] = v13*uca-v2*sa;
  rot[2][1] = v23*uca+v1*sa;
  rot[2][2] = ca+v33*uca;
}


/*=====================================================*/
void matvec(double mat[3][3], double vec[3], double ori[3], double tvp[3]) {
  /* multiplies mat times vec w/r/t the origin ori, and returns the
   result in tvp */

  double tmp;
  int m, n;

  for(m=0;m<3;m++){
    tmp=0.0;
    for(n=0;n<3;n++)
      tmp += mat[m][n]*(vec[n]-ori[n]);
    tvp[m] = tmp+ori[m];
  }
}


/*=====================================================================================*/
unsigned long bisearch (double *a, unsigned long n, double x) {
  /* Given an array a[0,...,n-1], and a value x, returns a value j such that
     x is between a[j] and a[j+1], but if x is out of range, it returns -1 or n-1.
     'a' must be monotonically increasing.
  */

  unsigned long low = -1, high = n, mid;

  while (low < high-1) {
    mid = low + (high - low)/2;
    if (a[mid] > x) high = mid;
    else low = mid;
  }

  if (x==a[0])   return 0;
  if (x==a[n-1]) return n-2;
  return low;
}


/*=====================================================================================*/
int linsolver(double *A, int n, int *perm, double *par, double *b) {
  /*
   Solves the set of n linear equations A*x = b by the LU decomposition method.
   Inputs: ----
   'A' is given as 1D row of values, of size n^2:  A[i][j] -> A[i*n+j].
   'b' is the right-hand side vector b.
   'n' is the number of unknowns.
   On exit: ----
   'A' is replaced by the LU decomposition of a row-wise permutation of A.
   'perm' is a vector that records the row permutation effected by the partial pivoting.
   'par' is the parity of the permutation 'perm' (+1 if even, -1 if odd).
   'b' is replaced by the solution vector x.
   This function returns -1 if the matrix is singular; 0 otherwise. */

  int i, imax, imaxtn, iperm, itn, j, jtn, jtnj, k, ktn, ku;
  double rownorm, tmp, dmax, cmax;
  double *rowscale;

  rowscale = (double *) malloc(n*sizeof(double));
  *par = 1.0;

  for (i=0,itn=0; i<n; i++,itn+=n) {    // itn = i*n
    rownorm = 0.0;
    for (j=0; j<n; j++)
      if ((tmp=fabs(A[itn+j])) > rownorm) rownorm=tmp;
    if (rownorm == 0.0) {
      return -1;    // singular matrix
    }
    rowscale[i] = 1.0/rownorm;
  }

  for (j=0,jtn=0; j<n; j++,jtn+=n) {    // jtn = j*n

    for (i=0,itn=0; i<j; i++,itn+=n) {
      tmp = A[itn+j];
      for (k=0,ktn=0; k<i; k++,ktn+=n) tmp -= A[itn+k] * A[ktn+j];  // ktn = k*n
      A[itn+j]=tmp;
    }

    dmax = 0.0;
    for (i=j,itn=jtn; i<n; i++,itn+=n) {
      tmp = A[itn+j];
      for (k=0,ktn=0; k<j; k++,ktn+=n) tmp -= A[itn+k] * A[ktn+j];
      A[itn+j] = tmp;
      if ((cmax=rowscale[i]*fabs(tmp)) >= dmax) {
        dmax = cmax;
        imax = i;
      }
    }

    if (j != imax) {
      imaxtn = imax*n;
      for (k=0; k<n; k++) {
        tmp = A[imaxtn+k];
        A[imaxtn+k] = A[jtn+k];
        A[jtn+k] = tmp;
      }
      *par = -(*par);
      rowscale[imax] = rowscale[j];
    }

    perm[j] = imax;
    jtnj = jtn+j;   // jtnj = j*n+j
    if (A[jtnj] == 0.0) A[jtnj] = 1.0e-20;
    if (j != n-1) {
      tmp = 1.0/A[jtnj];
      for (i=j+1,itn=jtn+n; i<n; i++,itn+=n)
        A[itn+j] *= tmp;
    }
  }

  free(rowscale);

  /* forward substitution: */
  ku = -1;
  for (i=0,itn=0; i<n; i++,itn+=n) {
    iperm = perm[i];
    tmp = b[iperm];
    b[iperm] = b[i];
    if (ku>=0) 
	    for (j=ku; j<=i-1; j++) tmp -= A[itn+j]*b[j];
    else if (tmp) ku=i;
    b[i] = tmp;
  }

  /* back substitution: */
  for (i=n-1,itn=(n-1)*n; i>=0; i--,itn-=n) {
    tmp = b[i];
    for (j=i+1; j<n; j++) 
	    tmp -= A[itn+j]*b[j];
    b[i] = tmp/A[itn+i];
  }

  return 0;

}


/*=====================================================================================*/
double J3(double *Qs, int m, int M, double b, double c, double k) {
  /* computes the function J3, used to perform the exponential regression
   involved in the determination of the stopping time */

  int i;
  double T1, T2, T3, f1, f2;

  T1=0.0; T2=0.0; T3=0.0;
  for (i=m; i<=M; i++) {
    f1 = exp(-k*i);
    f2 = i*f1;
    T1 += f2;
    T2 += f2*f1;
    T3 += Qs[i]*f2;
  }
  return (T1*b-T2*c)/T3 - 1.0;
}

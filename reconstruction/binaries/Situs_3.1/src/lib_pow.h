#ifndef __SITUS_LIB_POW
#define __SITUS_LIB_POW

#include "lib_vec.h"

/* header file for lib_pow.c */

#define SWAPPING(_a,_b,_type) \
  {\
    _type _tmp;\
    \
    _tmp = (_a);\
    (_a) = (_b);\
    (_b) = _tmp;\
  }

#ifdef __cplusplus
#extern "C" {
#endif

void powell(int *, double *, double *, int,
            double(*)(double *),
            unsigned, FILE *, double, double *);
int random_bit();
double square(double);
void fit_parabola_3p(int *, double *, double *, double *,
                     double, double, double, double, double, double);
int number_good(double);
int numcmp(double, double);
void insert_new_bracket_point(double *, double *, double *, double *,
                              double *, double *, double, double);
void expand_width(double *, double, double, double, double, double, double);
int xnew_is_new(double, double, double, double);
void golden_section_probe(int *, double *, double, double, double);
void insert_new_trap_point(double *, double *, double *, double *,
                           double *, double *, double, double);
void trap_minimum(double *, double *, double *,
                  double *, double *, double *, double,
                  double(*)(double), int);
void bracket_minimum(int *, double *, double *, double *, double *,
                     double *, double *,
                     double(*)(double));
void order_args(double *, double *, double *, double *, double *, double *);
void minimize_1d_raw(int *, double *, double *, double *,
                     double *, double *, double *, double,
                     double(*)(double), int);
double powell_line_function(double);
void powell_line_minimize(double *, double *, int, double *, double (*)(double *));

typedef struct Powell_Scratch *POWSCRPTR;
typedef struct Powell_Scratch {
  int           save_len;
  double        *buffer_vect;
  double        *save_vect;
  double        *save_direction;
  double (*save_pfunction)(double *vect, POWSCRPTR scratch);
  double        curr_shift[3];
  unsigned long l;
  PDB           *pdb;
  double        *phi_hi;
  double        *phi_du;
  double        corr_hi_lo;
} POWSCR;

typedef struct Powell_Result_Vector *POWITERRESPTR;
typedef struct Powell_Result_Vector {
  POWITERRESPTR next;
  int           iter;
  double        corr;
  double        *res;
} POWITERRES;

typedef struct Powell_Results {
  POWITERRESPTR head;
  POWITERRESPTR last;
} POWRES;

void powell_r(int *, double *, double *, int,
              double(*)(double *, POWSCR *),
              unsigned, FILE *, double, double *, POWSCR *, POWRES *);
void trap_minimum_r(double *, double *, double *,
                    double *, double *, double *, double,
                    double(*)(double, POWSCR *), int, POWSCR *);
void bracket_minimum_r(int *, double *, double *, double *, double *,
                       double *, double *,
                       double(*)(double, POWSCR *), POWSCR *);
void minimize_1d_raw_r(int *, double *, double *, double *,
                       double *, double *, double *, double,
                       double(*)(double, POWSCR *), int, POWSCR *);
double powell_line_function_r(double, POWSCR *);
void powell_line_minimize_r(double *, double *, int, double *,
                            double (*)(double *, POWSCR *), POWSCR *);

#ifdef __cplusplus
}
#endif

#endif

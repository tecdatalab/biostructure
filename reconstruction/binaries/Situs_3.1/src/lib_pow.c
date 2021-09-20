/*********************************************************************
*                           L I B _ P O W                            *
**********************************************************************
* Library is part of the Situs package URL: situs.biomachina.org     *
* (c) Jochen Heyd, Pablo Chacon, and Willy Wriggers, 2001-2005       *
**********************************************************************
*                                                                    *
* Powell's conjugate gradient optimization method.                   *
*                                                                    *
**********************************************************************
* Modified from WNLIB code, (c) 1996, Will Naylor                    *
* The original code by Will Naylor comes with its own license        *
* statement below.                                                   *
* The modifications are distributed under the terms of the           *
* Situs legal statement.                                             *
*********************************************************************/

#include "situs.h"
#include "lib_pow.h"
#include "lib_rnd.h"
#include "lib_vec.h"

#define GOLDEN_RATIO 1.618034
#define GOLDEN_SECTION 0.3819660
#define MAX_MAGNIFICATION 1.0
#define MAX_ITER_BRACKET 100


/* global variables */
static int g_save_len;
static double *g_buffer_vect, *g_save_vect, *g_save_direction;
static double (*g_save_pfunction)(double vect[]);

/*====================================================================*/
int random_bit()
{
  if (genrand() > 0.5) return 1;
  else return 0;
}

/*====================================================================*/
double square(double x)
{
  return (x * x);
}

/*====================================================================*/
/* fit parabola to 3 points (x,y) y = a*(x-x0)^2 + b; Return a, x0, b. */
void fit_parabola_3p(int *pcode, double *pa, double *px0, double *pb,
                     double x1, double y1, double x2, double y2,
                     double x3, double y3)
{
  double x12, x23, dy12, dy23, ddy, diff;

  /* sort the x's */
  if (!(x1 < x3)) {
    SWAPPING(x1, x3, double);
    SWAPPING(y1, y3, double);
  }
  if (x2 < x1) {
    SWAPPING(x1, x2, double);
    SWAPPING(y1, y2, double);
  } else if (x3 < x2) {
    SWAPPING(x2, x3, double);
    SWAPPING(y2, y3, double);
  }

  if ((x1 == x2) || (x2 == x3) || (x1 == x3)) {
    *pcode = (-1);
    return;
  }

  x12 = 0.5 * (x1 + x2);
  x23 = 0.5 * (x2 + x3);

  dy12 = (y2 - y1) / (x2 - x1);
  dy23 = (y3 - y2) / (x3 - x2);

  ddy = dy23 - dy12;
  *pa = 0.5 * ddy / (x23 - x12);

  if (ddy == 0.0) {
    *px0 = 0.0;
    *pb = y2;
  } else {
    *px0 = (dy23 * x12 - dy12 * x23) / ddy;
    diff = x2 - (*px0);
    *pb = y2 - (*pa) * diff * diff;
  }

  *pcode = 0;
}

/*====================================================================*/
int number_good(double x)
{
  return (x < 1.0e+50);
}

/*====================================================================*/
int numcmp(double x1, double x2)
{
  if (x1 < x2) return (-1);
  else if (x1 > x2) return (1);
  else if (x1 == x2) return (0);
  else {
    int good1, good2;
    good1 = number_good(x1);
    good2 = number_good(x2);
    if (good1 && (!good2)) return (-1);
    else if ((!good1) && good2)return (1);
    else return (0);
  }
}

/*====================================================================*/
void insert_new_bracket_point(double *pf0, double *pf1, double *pf2,
                              double *px0, double *px1, double *px2,
                              double xnew, double fnew)
{
  if (xnew < *px0) {
    /* small bracket found */
    if (numcmp(*pf0, *pf1) < 0) {
      /* make it as small as possible */
      *px2 = *px1;
      *pf2 = *pf1;
      *px1 = *px0;
      *pf1 = *pf0;
      *px0 = xnew;
      *pf0 = fnew;
    } else {
      /* make range as big as possible */
      *px0 = xnew;
      *pf0 = fnew;
    }
  } else if (xnew > *px2) {
    /* small bracket found */
    if (numcmp(*pf2, *pf1) < 0) {
      /* make it as small as possible */
      *px0 = *px1;
      *pf0 = *pf1;
      *px1 = *px2;
      *pf1 = *pf2;
      *px2 = xnew;
      *pf2 = fnew;
    } else {
      /* make range as big as possible */
      *px2 = xnew;
      *pf2 = fnew;
    }
  }
}

/*====================================================================*/
void expand_width(double *pwidth, double f0, double f1, double f2,
                  double x0, double x1, double x2)
{
  int code;
  double a, cen, b, diff;

  (*pwidth) *= GOLDEN_RATIO;

  fit_parabola_3p(&code, &a, &cen, &b, x0, f0, x1, x1, x2, f2);

  if ((code != 0) || (a <= 0.0))(*pwidth) *= MAX_MAGNIFICATION / GOLDEN_RATIO;
  else {
    diff = cen - x1;
    diff = fabs(diff);
    diff *= 2.0;

    if (diff > *pwidth)  {
      if (diff > (*pwidth)*MAX_MAGNIFICATION / GOLDEN_RATIO)  {
        (*pwidth) *= MAX_MAGNIFICATION / GOLDEN_RATIO;
      }
    } else (*pwidth) = diff;
  }
}

/*====================================================================*/
int xnew_is_new(double xnew, double x0, double x1, double x2)
{
  return ((xnew > x0) && (xnew != x1) && (xnew < x2));
}

/*====================================================================*/
void golden_section_probe(int *psuccess, double *pxnew,
                          double x0, double x1, double x2)
{
  double diff01, diff12;

  diff01 = x1 - x0;
  diff12 = x2 - x1;

  if (diff01 > diff12) *pxnew = x1 - GOLDEN_SECTION * diff01;
  else *pxnew = x1 + GOLDEN_SECTION * diff12;
  *psuccess = xnew_is_new(*pxnew, x0, x1, x2);
}

/*====================================================================*/
void insert_new_trap_point(double *pf0, double *pf1, double *pf2,
                           double *px0, double *px1, double *px2,
                           double xnew, double fnew)
{
  if (xnew < *px1) {
    if (numcmp(fnew, *pf1) < 0) {
      *px2 = *px1;
      *pf2 = *pf1;
      *px1 = xnew;
      *pf1 = fnew;
    } else {
      *px0 = xnew;
      *pf0 = fnew;
    }
  } else if (xnew > *px1) {
    if (numcmp(fnew, *pf1) < 0) {
      *px0 = *px1;
      *pf0 = *pf1;
      *px1 = xnew;
      *pf1 = fnew;
    } else {
      *px2 = xnew;
      *pf2 = fnew;
    }
  }
}

/*====================================================================*/
/* assume *pf0,*pf1,*pf2 already computed */
void trap_minimum(double *pf0, double *pf1, double *pf2,
                  double *px0, double *px1, double *px2, double f_goal,
                  double (*pfunction)(double x), int max_iterations)
{
  double xnew, fnew;
  int success;
  int improvement_count;

  improvement_count = 0;

  if ((numcmp(*pf0, f_goal) <= 0) ||
      (numcmp(*pf2, f_goal) <= 0))++improvement_count;

  while ((numcmp(*pf0, *pf1) > 0) || (numcmp(*pf2, *pf1) > 0)) {
    if (improvement_count >= max_iterations) break;
    golden_section_probe(&success, &xnew, *px0, *px1, *px2);
    if (!success) break;
    fnew = (*pfunction)(xnew);
    insert_new_trap_point(pf0, pf1, pf2, px0, px1, px2, xnew, fnew);

    if ((numcmp(*pf1, f_goal) <= 0) &&
        (numcmp(fnew, f_goal) <= 0)) ++improvement_count;
  }
}

/*====================================================================*/
void bracket_minimum(int *pcode, double *pf0, double *pf1, double *pf2,
                     double *px0, double *px1, double *px2,
                     double (*pfunction)(double x))
{
  int iteration_count, cmp;
  double xnew, fnew, widthg;

  widthg = *px2 - *px0;

  iteration_count = 0;

  for (;;) {
    /* minimum successfully bracketed */
    if ((numcmp(*pf1, *pf0) < 0) && (numcmp(*pf1, *pf2) < 0)) {
      *pcode = 0;
      return;
    }
    ++iteration_count;
    if (iteration_count > MAX_ITER_BRACKET) {
      *pcode = -2;
      return;
    }

    expand_width(&widthg, *pf0, *pf1, *pf2, *px0, *px1, *px2);
    cmp = numcmp(*pf0, *pf2);
    if (cmp < 0)xnew = *px0 - widthg;
    else if (cmp > 0) xnew = *px2 + widthg;
    else if (random_bit()) xnew = *px0 - widthg;
    else xnew = *px2 + widthg;

    fnew = (*pfunction)(xnew);
    insert_new_bracket_point(pf0, pf1, pf2, px0, px1, px2, xnew, fnew);
  }
}

/*====================================================================*/
void order_args(double *pf0, double *pf1, double *pf2, double *px0,
                double *px1, double *px2)
{
  if (*px0 > *px2) {
    SWAPPING(*px0, *px2, double);
    SWAPPING(*pf0, *pf2, double);
  }
  if (*px0 > *px1) {
    SWAPPING(*px0, *px1, double);
    SWAPPING(*pf0, *pf1, double);
  } else {
    if (*px1 > *px2) {
      SWAPPING(*px1, *px2, double);
      SWAPPING(*pf1, *pf2, double);
    }
  }
}

/*====================================================================*/
/* assume *pf0,*pf1,*pf2 already computed */
void minimize_1d_raw(int *pcode, double *pf0, double *pf1, double *pf2,
                     double *px0, double *px1, double *px2, double f_goal,
                     double (*pfunction)(double x), int max_iterations)
{
  double original_f1, original_x1;

  order_args(pf0, pf1, pf2, px0, px1, px2);

  original_f1 = *pf1;
  original_x1 = *px1;

  bracket_minimum(pcode, pf0, pf1, pf2, px0, px1, px2, pfunction);

  if (*pcode == 0)
    trap_minimum(pf0, pf1, pf2, px0, px1, px2, f_goal, pfunction, max_iterations);
  else if (*pcode == -2) {
    if ((original_f1 == *pf0) && (original_f1 == *pf1) && (original_f1 == *pf2)) {
      *px0 = *px1 = *px2 = original_x1;
      *pcode = 0;
    }
  }
}

/*====================================================================*/
double powell_line_function(double x)
{
  cp_vect(&g_save_vect, &g_buffer_vect, g_save_len);
  add_scaled_vect(g_save_vect, g_save_direction, x, g_save_len);
  return ((*g_save_pfunction)(g_save_vect));
}

/*====================================================================*/
void powell_line_minimize(double *vect, double *direction,
                          int len, double *pval_min,
                          double (*pfunction)(double *vect))
{
  double ax, bx, cx, fa, fb, fc;
  int code;

  cp_vect(&g_buffer_vect, &vect, len);
  g_save_vect = vect;
  g_save_len = len;
  g_save_direction = direction;
  g_save_pfunction = pfunction;

  ax = -1.0;
  bx = 0.0;
  cx = 1.0;
  fa = powell_line_function(ax);
  fb = *pval_min;
  fc = powell_line_function(cx);

  minimize_1d_raw(&code, &fa, &fb, &fc, &ax, &bx, &cx, fb, &powell_line_function, 50);

  cp_vect(&vect, &g_buffer_vect, len);

  if (*pval_min == fb) return; /* do not move if no improvement */

  add_scaled_vect(vect, direction, bx, len);

  *pval_min = fb;
}

/*====================================================================*/
/* conjugate gradient method, called with external score function and output file */
void powell(int *pcode, double *pval_min, double *vect,
            int len, double (*pfunction)(double *vect),
            unsigned max_iter, FILE *pow_outfile,
            double tolerance, double *init_directions)
{
  int i, ibig = 0, j, iteration;
  double t, fptt, fp, del;
  double *pt, *ptt, *xit;
  double **xi;
  double last_val_min;

  iteration = 0;

  do_vect(&g_buffer_vect, len);
  do_mat(&xi, len, len);
  do_vect(&pt, len);
  do_vect(&ptt, len);
  do_vect(&xit, len);

  /*initialize_powell_directions*/
  zero_mat(xi, len, len);
  for (i = 0; i < len; ++i) {
    xi[i][i] = init_directions[i];
  }

  *pval_min = (*pfunction)(vect);
  cp_vect(&pt, &vect, len);

  for (iteration = 0;; ++iteration) {
    last_val_min = *pval_min;
    fp = *pval_min;
    del = 0.0;
    for (i = 0; i < len; ++i) {
      fptt = *pval_min;
      if (init_directions[i] != 0.0) {
        powell_line_minimize(vect, xi[i], len, pval_min, pfunction);
        if (fabs(fptt - (*pval_min)) > del) {
          del = fabs(fptt - (*pval_min));
          ibig = i;
        }
      }
    }

    if (*pval_min == last_val_min) {
      *pcode = 0;
      free_vect_and_zero_ptr(&xit);
      free_vect_and_zero_ptr(&ptt);
      free_vect_and_zero_ptr(&pt);
      free_vect_and_zero_ptr(&g_buffer_vect);
      free_mat_and_zero_ptr(&xi);
      return;
    }

    for (i = 0; i < len; ++i) fprintf(pow_outfile, " %7.3f", vect[i]);

    /* uncomment one of the following options */
    /* maximization */
    fprintf(pow_outfile, "    %10.7E %d\n", -(*pval_min), iteration + 1);

    /* minimization */
    /* fprintf(pow_outfile,"%10.7E %d\n",(*pval_min), iteration+1); */
    if (iteration >= max_iter - 1) {
      *pcode = 1;
      free_vect_and_zero_ptr(&xit);
      free_vect_and_zero_ptr(&ptt);
      free_vect_and_zero_ptr(&pt);
      free_vect_and_zero_ptr(&g_buffer_vect);
      free_mat_and_zero_ptr(&xi);
      return;
    }

    if (2.0 * fabs(fp - (*pval_min)) <= tolerance * (fabs(fp) + fabs(*pval_min))) {
      *pcode = 0;
      free_vect_and_zero_ptr(&xit);
      free_vect_and_zero_ptr(&ptt);
      free_vect_and_zero_ptr(&pt);
      free_vect_and_zero_ptr(&g_buffer_vect);
      free_mat_and_zero_ptr(&xi);
      return;
    }

    for (j = 0; j < len; ++j) {
      ptt[j] = 2.0 * vect[j] - pt[j];
      xit[j] = vect[j] - pt[j];
    }

    cp_vect(&pt, &vect, len);
    fptt = (*pfunction)(ptt);
    if (fptt < fp) {
      t = 2.0 * (fp - 2.0 * (*pval_min) + fptt) * square(fp - (*pval_min) - del) - del * square(fp - fptt);
      if (t < 0.0) {
        powell_line_minimize(vect, xit, len, pval_min, pfunction);
        cp_vect(&(xi[ibig]), &xit, len);
      }
    }
  }
}

/*====================================================================*/
/* Thread-safe versions of the above functions added by J.H.          */
/*====================================================================*/

/*====================================================================*/
/* assume *pf0,*pf1,*pf2 already computed */
void trap_minimum_r(double *pf0, double *pf1, double *pf2,
                    double *px0, double *px1, double *px2, double f_goal,
                    double (*pfunction)(double x, POWSCR *scratch),
                    int max_iterations, POWSCR *scratch)
{

  double xnew, fnew;
  int success;
  int improvement_count;

  improvement_count = 0;

  if ((numcmp(*pf0, f_goal) <= 0) ||
      (numcmp(*pf2, f_goal) <= 0))++improvement_count;

  while ((numcmp(*pf0, *pf1) > 0) || (numcmp(*pf2, *pf1) > 0)) {
    if (improvement_count >= max_iterations) break;
    golden_section_probe(&success, &xnew, *px0, *px1, *px2);
    if (!success) break;
    fnew = (*pfunction)(xnew, scratch);
    insert_new_trap_point(pf0, pf1, pf2, px0, px1, px2, xnew, fnew);

    if ((numcmp(*pf1, f_goal) <= 0) &&
        (numcmp(fnew, f_goal) <= 0)) ++improvement_count;
  }
}

/*====================================================================*/
void bracket_minimum_r(int *pcode, double *pf0, double *pf1, double *pf2,
                       double *px0, double *px1, double *px2,
                       double (*pfunction)(double x, POWSCR *scratch),
                       POWSCR *scratch)
{

  int iteration_count, cmp;
  double xnew, fnew, widthg;

  widthg = *px2 - *px0;

  iteration_count = 0;

  for (;;) {
    /* minimum successfully bracketed */
    if ((numcmp(*pf1, *pf0) < 0) && (numcmp(*pf1, *pf2) < 0)) {
      *pcode = 0;
      return;
    }
    ++iteration_count;
    if (iteration_count > MAX_ITER_BRACKET) {
      *pcode = -2;
      return;
    }

    expand_width(&widthg, *pf0, *pf1, *pf2, *px0, *px1, *px2);
    cmp = numcmp(*pf0, *pf2);
    if (cmp < 0)xnew = *px0 - widthg;
    else if (cmp > 0) xnew = *px2 + widthg;
    else if (random_bit()) xnew = *px0 - widthg;
    else xnew = *px2 + widthg;

    fnew = (*pfunction)(xnew, scratch);
    insert_new_bracket_point(pf0, pf1, pf2, px0, px1, px2, xnew, fnew);
  }
}

/*====================================================================*/
/* assume *pf0,*pf1,*pf2 already computed */
void minimize_1d_raw_r(int *pcode, double *pf0, double *pf1, double *pf2,
                       double *px0, double *px1, double *px2, double f_goal,
                       double (*pfunction)(double x, POWSCR *scratch),
                       int max_iterations, POWSCR *scratch)
{

  double original_f1, original_x1;

  order_args(pf0, pf1, pf2, px0, px1, px2);

  original_f1 = *pf1;
  original_x1 = *px1;

  bracket_minimum_r(pcode, pf0, pf1, pf2, px0, px1, px2, pfunction, scratch);

  if (*pcode == 0)
    trap_minimum_r(pf0, pf1, pf2, px0, px1, px2, f_goal, pfunction, max_iterations, scratch);
  else if (*pcode == -2) {
    if ((original_f1 == *pf0) && (original_f1 == *pf1) && (original_f1 == *pf2)) {
      *px0 = *px1 = *px2 = original_x1;
      *pcode = 0;
    }
  }
}

/*====================================================================*/
double powell_line_function_r(double x, POWSCR *scratch)
{

  cp_vect(&scratch->save_vect, &scratch->buffer_vect, scratch->save_len);
  add_scaled_vect(scratch->save_vect, scratch->save_direction, x, scratch->save_len);
  return ((*scratch->save_pfunction)(scratch->save_vect, scratch));
}

/*====================================================================*/
void powell_line_minimize_r(double *vect, double *direction,
                            int len, double *pval_min,
                            double (*pfunction)(double *vect, POWSCR *scratch), POWSCR *scratch)
{

  double ax, bx, cx, fa, fb, fc;
  int code;

  cp_vect(&scratch->buffer_vect, &vect, len);
  scratch->save_vect = vect;
  scratch->save_len = len;
  scratch->save_direction = direction;
  scratch->save_pfunction = pfunction;

  ax = -1.0;
  bx = 0.0;
  cx = 1.0;
  fa = powell_line_function_r(ax, scratch);
  fb = *pval_min;
  fc = powell_line_function_r(cx, scratch);

  minimize_1d_raw_r(&code, &fa, &fb, &fc, &ax, &bx, &cx, fb,
                    &powell_line_function_r, 50, scratch);

  cp_vect(&vect, &scratch->buffer_vect, len);

  if (*pval_min == fb) return; /* do not move if no improvement */

  add_scaled_vect(vect, direction, bx, len);

  *pval_min = fb;
}

/*====================================================================*/
/* conjugate gradient method, called with external score function and output file */
void powell_r(int *pcode, double *pval_min, double *vect,
              int len, double (*pfunction)(double *vect, POWSCR *scratch),
              unsigned max_iter, FILE *pow_outfile,
              double tolerance, double *init_directions,
              POWSCR *scratch, POWRES *results)
{
  int i, ibig = 0, j, iteration;
  double t, fptt, fp, del;
  double *pt, *ptt, *xit;
  double **xi;
  double last_val_min;
  POWITERRES *curr_result = NULL;

  iteration = 0;

  do_vect(&scratch->buffer_vect, len);
  do_mat(&xi, len, len);
  do_vect(&pt, len);
  do_vect(&ptt, len);
  do_vect(&xit, len);

  /*initialize_powell_directions*/

  zero_mat(xi, len, len);
  for (i = 0; i < len; ++i) {
    xi[i][i] = init_directions[i];
  }

  *pval_min = (*pfunction)(vect, scratch);
  cp_vect(&pt, &vect, len);

  for (iteration = 0;; ++iteration) {
    last_val_min = *pval_min;
    fp = *pval_min;
    del = 0.0;
    for (i = 0; i < len; ++i) {
      fptt = *pval_min;
      if (init_directions[i] != 0.0) {
        powell_line_minimize_r(vect, xi[i], len, pval_min, pfunction, scratch);
        if (fabs(fptt - (*pval_min)) > del) {
          del = fabs(fptt - (*pval_min));
          ibig = i;
        }
      }
    }

    if (*pval_min == last_val_min) {
      *pcode = 0;
      free_vect_and_zero_ptr(&xit);
      free_vect_and_zero_ptr(&ptt);
      free_vect_and_zero_ptr(&pt);
      free_vect_and_zero_ptr(&(scratch->buffer_vect));
      free_mat_and_zero_ptr(&xi);
      return;
    }

    /* make new result vector */

    curr_result = (POWITERRES *) alloc_vect(sizeof(POWITERRES), 1);
    curr_result->iter = iteration + 1;
    curr_result->corr = -(*pval_min);
    curr_result->res  = (double *) alloc_vect(len, sizeof(double));
    curr_result->next = NULL;
    for (i = 0; i < len; ++i) curr_result->res[i] = vect[i];

    /* append it to the result struct */

    if (results->head == NULL) {
      results->head = curr_result;
      results->last = curr_result;
    } else {
      results->last->next = curr_result;
      results->last       = curr_result;
    }

    /* stop criteria */

    if (iteration >= max_iter - 1) {
      *pcode = 1;
      free_vect_and_zero_ptr(&xit);
      free_vect_and_zero_ptr(&ptt);
      free_vect_and_zero_ptr(&pt);
      free_vect_and_zero_ptr(&(scratch->buffer_vect));
      free_mat_and_zero_ptr(&xi);
      return;
    }

    if (2.0 * fabs(fp - (*pval_min)) <= tolerance * (fabs(fp) + fabs(*pval_min))) {
      *pcode = 0;
      free_vect_and_zero_ptr(&xit);
      free_vect_and_zero_ptr(&ptt);
      free_vect_and_zero_ptr(&pt);
      free_vect_and_zero_ptr(&(scratch->buffer_vect));
      free_mat_and_zero_ptr(&xi);
      return;
    }

    for (j = 0; j < len; ++j) {
      ptt[j] = 2.0 * vect[j] - pt[j];
      xit[j] = vect[j] - pt[j];
    }

    cp_vect(&pt, &vect, len);
    fptt = (*pfunction)(ptt, scratch);
    if (fptt < fp) {
      t = 2.0 * (fp - 2.0 * (*pval_min) + fptt) * square(fp - (*pval_min) - del) - del * square(fp - fptt);
      if (t < 0.0) {
        powell_line_minimize_r(vect, xit, len, pval_min, pfunction, scratch);
        cp_vect(&(xi[ibig]), &xit, len);
      }
    }
  }
}

/* Original license statement by Will Naylor regarding the unmodified functions not ending in _r:

BY DOWNLOADING OR COMPILING THIS LIBRARY, YOU ACCEPT AND AGREE TO THE TERMS
AND CONDITIONS PRINTED BELOW.IF YOU DO NOT AGREE, DO NOT DOWNLOAD OR
COMPILE THIS LIBRARY.

The author provides this C code in the hope
that it will be helpful, however, we assume no responsibility
for the use of this code, nor any responsibility for its support.
The software is distributed on an "AS IS" basis, without warranty.
Neither the authors nor the software developers
make any representation, or warranty, either express or implied, with
respect to the software programs and subroutines, their quality, accuracy,
or fitness for a specific purpose.Therefore, neither the authors nor the
software developers shall have any liability to
you or any other person or entity with respect to any liability, loss,
or damage caused or alleged to have been caused directly or indirectly by
the programs and subroutines contained in this library.This includes, but
is not limited to, interruption of service, loss of data, loss of classroom
time, loss of consulting or anticipatory profits, or consequential damages
from the use of these programs.

COPYRIGHT NOTICE:

The source code in this file is provided free of charge to anybody
who wants it.It is in the public domain and therefore may be used by
anybody for any purpose.This copyright notice and the above legal notice
may not be removed.

AUTHOR:

Will Naylor
PO Box 700759
San Jose, CA
95170-0759

WNLIB RELEASE:

16 May 96 Release 6.0

WEB SITE:

http://www.willnaylor.com/wnlib.html

*/


/*********************************************************************
*                           C O L O R E S                            *
**********************************************************************
* Program is part of the Situs package URL: situs.biomachina.org     *
* (c) Pablo Chacon, Jochen Heyd, Valerio Mariani, John Heumann, and  *
* Willy Wriggers 2001-2019                                           *
**********************************************************************
*                                                                    *
* General purpose FFT-accelerated fitting tool.                      *
*                                                                    *
**********************************************************************
* See legal statement for terms of distribution                      *
*********************************************************************/

#include "situs.h"
#include "lib_err.h"
#include "lib_rnd.h"
#include "lib_pow.h"
#include "lib_vec.h"
#include "lib_vwk.h"
#include "lib_pwk.h"
#include "lib_eul.h"
#include "lib_tim.h"
#include "lib_vio.h"
#include "lib_pio.h"
#include "lib_smp.h"
#include "fftw3.h"

#ifdef _OPENMP
#include <omp.h> /* openmp used here only to find number of default SMP threads, since POSIX does not support a portable solution*/
#endif

#define DEF_STEP  1.0    /* deflection step in degrees for checking angular variability about axes */
#define DEF_RANGE 30.0   /* deflection range in degrees for checking angular variability about axes */

/* Macros for fftw2 style access */
#define c_re(c) ((c)[0])
#define c_im(c) ((c)[1])

/* map and kernel related global variables */

static unsigned g_extx;                      /* map extent */
static unsigned g_exty;                      /* map extent */
static unsigned g_extz;                      /* map extent */
static unsigned long g_nvox;                 /* number of voxels */
static unsigned g_extx_half;                 /* half map extent */
static unsigned g_exty_half;                 /* half map extent */
static unsigned g_extz_half;                 /* half map extent */
static double g_width;                       /* voxel size in Angstroms */
static double g_origx;                       /* map origin */
static double g_origy;                       /* map origin */
static double g_origz;                       /* map origin */
static double *g_phi_lo;                     /* low resolution map */
static double *g_phi_hi;                     /* high resolution map */
static double *g_phi_du;                     /* dummy map */
static double *g_phi_fi;                     /* filter kernel */
static double *g_phi_ga;                     /* low-pass (Gaussian) kernel */
static double *g_phi_fx;                     /* filtered low-pass kernel */
static double g_norm_hi, g_norm_lo;          /* normalization factors */
static double g_center_map[3];               /* center of the map in real space coordinates */
static unsigned g_ext_ga;                    /* low-pass kernel linear extent */
static unsigned g_ext_fi;                    /* filter kernel linear extent */
static unsigned g_ext_fx;                    /* filtered low-pass kernel linear extent */
static unsigned long g_nvox_ga;              /* low-pass kernel voxel count */
static unsigned long g_nvox_fi;              /* filter kernel voxel count */

/* FFT and FFTW related global variables */

static unsigned long g_fftw_nvox_r2c;        /* FFTW real to complex voxel count */
static unsigned long g_fftw_nvox_c2r;        /* FFTW complex to real voxel count */
static fftw_complex *g_fftw_grid_a;          /* A structure factors */
static fftw_complex *g_fftw_grid_b;          /* B structure factors */
static fftw_plan g_fftw_plan_fwd_lo;         /* low res forward FFT plan */
static fftw_plan g_fftw_plan_rev_lo;         /* low res reverse FFT plan */
static fftw_plan g_fftw_plan_fwd_hi;         /* high res forward FFT plan */
static fftw_plan g_fftw_plan_rev_hi;         /* high res reverse FFT plan */
static fftw_plan g_fftw_plan_fwd_du;         /* dummy array forward FFT plan */
static double g_fftw_scale;                  /* inverse g_nvox */
static unsigned g_ignored[3];                /* zero margin that can be safely ignored in fast kernel convolutions */

/* PDB related global variables */

static PDB *g_pdb_original;                  /* PDB coordinates */
static PDB *g_pdb_move;                      /* PDB coordinates */
static PDB *g_pdb_save;                      /* PDB coordinates */
static unsigned g_num_atoms;                 /* number of PDB atoms */

/* options related global variables */

static double g_delta_rot;                   /* rotational sampling step in degrees */
static double g_target_res;                  /* resolution in A */
static double g_target_ani;                  /* resolution anisotropy factor */
static unsigned g_num_explored;              /* number of explored best fits */
static char g_pow_mode;                      /* Powell option */
static unsigned g_pow_max_iter;              /* Powell max number of iterations */
static double g_pow_tolerance;               /* Powell tolerance */
static double g_pow_delta_pos;               /* Powell initialization position  */
static double g_pow_delta_ang;               /* Powell initialization orientation */
static int g_euler_mode;                     /* Euler angle option */
static char g_eu_in_file[69];                /* Euler angle input string */
static double g_eu_range[3][2];              /* Euler angle range */
static double g_low_cutoff;                  /* low res. map cutoff */
static int g_corr_mode;                      /* correlation option */
static double g_size_fac;                    /* grid size expansion factor for FFT padding */
static int g_peak_opt;                       /* peak search method */
static int g_peak_sharp;                     /* peak sharpness estimation option */
static int g_pow_alg;                        /* correlation method for Powell optimization */
static int g_sculptor;                       /* write special Sculptor output for interactive use */

/* Parallel processing related global variables */

static int g_p_nprocs;                       /* number of processors */

#ifdef _SMP_
static pthread_mutex_t g_p_fft_score_mutex;  /* mutex for FFT score update */
static pthread_mutex_t g_p_fft_log_mutex;    /* mutex for FFT log update */
static pthread_mutex_t g_p_fft_status_mutex; /* semaphore for FFT status */
static pthread_mutex_t g_p_fft_plan_mutex;   /* mutex for FFTW planner */
static pthread_cond_t  g_p_fft_status_cond;  /* semaphore for FFT status */
static int             g_p_fft_status_pred;  /* count of FFTs performed */
static int             g_p_fft_status_count; /* count of FFTs performed */
typedef struct {                             /* struct for FFT results */
  double psi;
  double theta;
  double phi;
  double x;
  double y;
  double z;
  double corr_orig;
  double corr_opt;
} FFT_RESULTS;
typedef struct {                             /* struct for passing arguments to FFT subroutine */
  int thread_id;
  int idx_start;
  int idx_count;
  FFT_RESULTS *results;
} FFT_THREAD_DATA;
FFT_THREAD_DATA thread_fft_data_array[MAX_NUM_THREADS];

workq_t *g_p_workq;

static pthread_mutex_t g_p_pow_print_mutex;  /* mutex for Powell print section */

static void *search6d_fft_par(void *thread_arg);
#endif

/* parameter storage and output related global variables */
typedef struct {
  double score;
  double pos[3];
  double euler[3];
} FIT;
typedef struct {
  unsigned eu;
  float score;
} SAV;
typedef struct {
  unsigned long ifft;
  unsigned long ireal;
  unsigned ix;
  unsigned iy;
  unsigned iz;
} POS;
static SAV *g_hash_sav;                      /* hash table for Euler angles and score associated with point on lattice */
static FIT *g_best_fit;                      /* saved parameters and coefficients of best fits */
static POS *g_inside_list;                   /* inside target positions */
static unsigned g_inside_num;                /* inside target number */
static POS *g_inside_list_flipped;           /* inside target positions */
static unsigned g_inside_num_flipped;        /* inside target number */
static float *g_eulers;                      /* Euler angle matrix */
static unsigned long g_eulers_count;         /* number of Euler angles */
static FILE *g_outfile;                      /* frequently used output file */
static char *g_program = "colores";




/* functions list, many of these are similar or identical to collage but use global variables so not placed in a library */

static void search6d_fft(unsigned);
static void powell_optimization(void *arg);
static double powell_correlation(double *, POWSCR *);
static void get_centered_structure_and_radius(char *, double *);
static void get_low_map(char *);
static void read_options(int, char **);
static void peak_extract(SAV **, FIT **);
static void flip_quadrants(unsigned, unsigned, unsigned, unsigned *, unsigned *, unsigned *);
static void create_inside_molecule_poslist_flipped(double *);
static void create_inside_molecule_poslist(double *);
static void draw_line();



int main(int argc, char **argv)
{

  the_time itime, etime;
  time_t seed_time;
  double sigma1d;
  double corr_hi_lo;
  double trial_time;
  double save_b_re;
  double struct_max_radius;
  unsigned num_degenerate;
  unsigned long indv;
  unsigned indx, indy, indz;
  unsigned i, j;
  unsigned long iter, p, q, m, n;
  double curr_shift[3];
  double zero_shift[3] = {0.0, 0.0, 0.0};
  double deflection_deg, deflection_rad;
  char out_string[20];
  char corr_string[70];
  double   *pow_init;
  POWARG   *pow_args;
  POWSCR   *pow_scratch;
  double   test_pow_corr, test_pow_corr_orig;
  double   test_pow_time, test_pow_time_orig;
  int      test_pow_method;
  unsigned zeropad[3];
  unsigned extx, exty, extz;
  double origx, origy, origz;
  unsigned ext_ga_save;
  unsigned long nvox_ga_save;
  double *phi_ga_save, *phi_fx_save;
  double sigma_factor = 0;
  FILE *out;
#ifdef _SMP_
  pthread_t  threads[MAX_NUM_THREADS];
  int        threads_workload[MAX_NUM_THREADS];
  int        thread;
  int        rc;
  int        p_fft_percent;
  char       p_fft_progressbar[51];
#endif


  /*========================= INITIALIZATION =================================*/

  sgenrand((unsigned long) time(&seed_time));

  draw_line();
  read_options(argc, &(*argv)); /* sets a variety of global variables */

  draw_line();
  get_low_map(argv[1]);
  draw_line();
  get_centered_structure_and_radius(argv[2], &struct_max_radius);
  draw_line();

  /* calculate sigma of kernel (1D) */
  sigma1d = g_target_res / (2.0 * g_width * sqrt(3.0));

  /* create Gaussian kernels, modify sigma factor arguments as necessary */
  create_gaussian(&g_phi_ga, &g_nvox_ga, &g_ext_ga, sigma1d, 3.0);
  create_gaussian(&phi_ga_save, &nvox_ga_save, &ext_ga_save, sigma1d, 5.0);

  /* create filter (_fi) kernel (e.g. Laplacian) and indicate (_fx) sigma factor */
  switch (g_corr_mode) {
    case 0:
      create_identity(&g_phi_fi, &g_nvox_fi, &g_ext_fi);
      strcpy(corr_string, "REMARK    Standard Linear Correlation");
      sigma_factor = 3.0;
      break;
    case 1:
      create_laplacian(&g_phi_fi, &g_nvox_fi, &g_ext_fi);
      strcpy(corr_string, "REMARK    Laplacian Correlation");
      sigma_factor = 4.0;
      break;
  }

  /* create convolved filter (_fx) kernel  (e.g. Mexican hat) */
  do_vect(&phi_fx_save, (ext_ga_save * ext_ga_save * ext_ga_save));
  convolve_kernel_inside(&phi_fx_save, phi_ga_save, ext_ga_save, ext_ga_save, ext_ga_save, g_phi_fi, g_ext_fi);
  do_vect(&g_phi_fx, (ext_ga_save * ext_ga_save * ext_ga_save));
  shrink_to_sigma_factor(&g_phi_fx, &g_ext_fx, phi_fx_save, ext_ga_save, sigma1d, sigma_factor);
  free_vect_and_zero_ptr(&phi_ga_save);
  free_vect_and_zero_ptr(&phi_fx_save);

  /* zeropad low-resolution map (consider kernel size as well as FFT requirements) */

  /* add padding for FFT */
  zeropad[0] = ceil(g_extx * g_size_fac);
  zeropad[1] = ceil(g_exty * g_size_fac);
  zeropad[2] = ceil(g_extz * g_size_fac);

  /* save g_ignored margin (i.e. the margin that will be ignored and simply set to zero in fast kernel convolutions) */
  for (m = 0; m < 3; m++) g_ignored[m] = zeropad[m];

  /* finally, we add half of the convolution kernel size to the padding, this is the kernel buffer region around original map */
  for (m = 0; m < 3; m++) zeropad[m] += (g_ext_fx - 1) / 2;

  /* now add extra padding to g_phi_lo, changing the map size parameters */
  create_padded_map(&g_phi_du, &g_extx, &g_exty, &g_extz, &g_origx, &g_origy, &g_origz, &g_nvox,
                    g_phi_lo, g_extx, g_exty, g_extz, g_origx, g_origy, g_origz,
                    g_width, g_width, g_width * g_target_ani, zeropad);
  free_vect_and_zero_ptr(&g_phi_lo);

  /* set up padding required for FFTW real xform... see FFTW manual) */
  g_fftw_nvox_r2c = g_extz * g_exty * (2 * (g_extx / 2 + 1));
  g_fftw_nvox_c2r = g_extz * g_exty * (g_extx / 2 + 1);

  /* re-allocate g_phi_lo on FFTW grid and copy content from padded map */
  cp_vect_destroy(&g_phi_lo, &g_phi_du, g_nvox);

  /* initialize some variables */
  g_fftw_scale = 1.0 / (double)g_nvox;
  g_extx_half = (g_extx - 1) / 2;
  g_exty_half = (g_exty - 1) / 2;
  g_extz_half = (g_extz - 1) / 2;
  g_center_map[0] = (g_extx / 2.0) * g_width + g_origx;
  g_center_map[1] = (g_exty / 2.0) * g_width + g_origy;
  g_center_map[2] = (g_extz / 2.0) * g_width * g_target_ani + g_origz;

  create_inside_molecule_poslist_flipped(g_phi_lo);
  create_inside_molecule_poslist(g_phi_lo);

  printf("colores> Memory allocation for FFT.\n");

  /* memory allocation */
  do_vect(&g_phi_du, g_fftw_nvox_r2c);
  do_vect(&g_phi_hi, g_fftw_nvox_r2c);

  g_fftw_grid_a = (fftw_complex *) alloc_vect(g_fftw_nvox_c2r, sizeof(fftw_complex));
  g_fftw_grid_b = (fftw_complex *) alloc_vect(g_fftw_nvox_c2r, sizeof(fftw_complex));

  printf("colores> FFT planning...\n");

  /*
   * Precalculate optimal FFTW plans. FFTW3 requires separate plans
   * for each input/output pair, and overwrites those arrays during
   * planning, so we reinitialize afterwards.
   */
  g_fftw_plan_fwd_du = fftw_plan_dft_r2c_3d(g_extz, g_exty, g_extx,
                       g_phi_du, (fftw_complex *) g_fftw_grid_b,
                       FFTW_ESTIMATE);
  cp_vect(&g_phi_du, &g_phi_lo, g_nvox);

  g_fftw_plan_fwd_lo = fftw_plan_dft_r2c_3d(g_extz, g_exty, g_extx,
                       g_phi_lo, (fftw_complex *) g_fftw_grid_a,
                       FFTW_ESTIMATE);
  g_fftw_plan_rev_lo = fftw_plan_dft_c2r_3d(g_extz, g_exty, g_extx,
                       (fftw_complex *) g_fftw_grid_a, g_phi_lo,
                       FFTW_ESTIMATE);
  cp_vect(&g_phi_lo, &g_phi_du, g_nvox);
  /* zero_vect(g_phi_du, g_nvox); Not needed: initialized before use */
  g_fftw_plan_fwd_hi = fftw_plan_dft_r2c_3d(g_extz, g_exty, g_extx,
                       g_phi_hi, (fftw_complex *) g_fftw_grid_b,
                       FFTW_ESTIMATE);
  g_fftw_plan_rev_hi = fftw_plan_dft_c2r_3d(g_extz, g_exty, g_extx,
                       (fftw_complex *) g_fftw_grid_b, g_phi_hi,
                       FFTW_ESTIMATE);
  /* zero_vect(g_phi_hi,  g_nvox); Not needed: initialized before use */
  draw_line();

  /* generate exhaustive search data for Sculptor interactive use */
  if (g_sculptor) {
    /* write out the target volume */
    write_vol("col_target.sit", g_width, g_origx, g_origy, g_origz,
              g_extx, g_exty, g_extz,
              g_phi_lo);
    /* write out the probe PDB */
    g_pdb_move = (PDB *) alloc_vect(g_num_atoms, sizeof(PDB));

    for (i = 0; i < g_num_atoms; ++i) *(g_pdb_move + i) = *(g_pdb_original + i);
    translate(g_pdb_original, g_pdb_move, g_num_atoms, g_center_map[0], g_center_map[1], g_center_map[2]);
    write_pdb("col_probe.pdb", g_num_atoms, g_pdb_move);
    free_vect_and_zero_ptr(&g_pdb_move);
  }


  /*============================== TESTING ===================================*/


  printf("colores> Testing the maps and correlations.\n");

  /* create the dummy structure g_pdb_move */
  g_pdb_move = (PDB *) alloc_vect(g_num_atoms, sizeof(PDB));
  for (i = 0; i < g_num_atoms; ++i)
    *(g_pdb_move + i) = *(g_pdb_original + i);

  /* rotate g_pdb_move (if desired), project to grid, and convolve with Gaussian */
  rot_euler(g_pdb_original, g_pdb_move, g_num_atoms, 0.0, 0.0, 0.0);
  printf("colores> Projecting probe structure to lattice...\n");
  project_mass(&g_phi_hi, g_nvox, g_width, g_width, g_width * g_target_ani, g_extx, g_exty, g_extz, g_pdb_move, g_num_atoms, zero_shift, g_ignored);
  printf("colores> Low-pass-filtering probe map...\n");
  convolve_kernel_inside(&g_phi_hi, g_phi_hi, g_extx, g_exty, g_extz, g_phi_ga, g_ext_ga);
  printf("colores> Target and probe maps:\n");
  print_map_info(g_phi_lo, g_nvox);
  print_map_info(g_phi_hi, g_nvox);

  /* redo high-resolution map */
  printf("colores> Projecting probe structure to lattice...\n");
  project_mass(&g_phi_hi, g_nvox, g_width, g_width, g_width * g_target_ani, g_extx, g_exty, g_extz, g_pdb_move, g_num_atoms, zero_shift, g_ignored);

  printf("colores> Applying filters to target and probe maps...\n");
  switch (g_corr_mode) {
    case 0:
      convolve_kernel_inside(&g_phi_lo, g_phi_lo, g_extx, g_exty, g_extz, g_phi_fi, g_ext_fi);
      break;
    case 1:
      relax_laplacian(&g_phi_lo, g_extx, g_exty, g_extz, g_ignored, 5.0);
      convolve_kernel_inside_erode(&g_phi_lo, g_phi_lo, g_extx, g_exty, g_extz, g_phi_fi, g_ext_fi);
      break;
  }
  convolve_kernel_inside(&g_phi_hi, g_phi_hi, g_extx, g_exty, g_extz, g_phi_fx, g_ext_fx);

  /* normalization */
  printf("colores> Normalizing target and probe maps...\n");
  g_norm_hi = calc_norm(g_phi_hi, g_nvox);
  g_norm_lo = calc_norm(g_phi_lo, g_nvox);
  normalize(g_phi_hi, g_nvox, g_norm_hi);
  normalize(g_phi_lo, g_nvox, g_norm_lo);

  printf("colores> Target and probe maps:\n");
  print_map_info(g_phi_lo, g_nvox);
  print_map_info(g_phi_hi, g_nvox);

  printf("colores> Writing target and probe maps for inspection or debugging...\n");

  /* write hi res debugging map */
  if (g_target_ani != 1.0) {
    free_vect_and_zero_ptr(&g_phi_du);
    interpolate_map(&g_phi_du, &extx, &exty, &extz, &origx, &origy, &origz,
                    g_width, g_width, g_width, g_phi_hi, g_extx, g_exty, g_extz, g_origx,
                    g_origy, g_origz, g_width, g_width, g_width * g_target_ani);
    write_vol("col_hi_fil.sit", g_width, origx, origy, origz, extx, exty, extz, g_phi_du);
    free_vect_and_zero_ptr(&g_phi_du);
    do_vect(&g_phi_du, g_fftw_nvox_r2c);
  } else write_vol("col_hi_fil.sit", g_width, g_origx, g_origy, g_origz, g_extx, g_exty, g_extz, g_phi_hi);

  /* write lo res debugging map */
  if (g_target_ani != 1.0) {
    free_vect_and_zero_ptr(&g_phi_du);
    interpolate_map(&g_phi_du, &extx, &exty, &extz, &origx, &origy, &origz,
                    g_width, g_width, g_width, g_phi_lo, g_extx, g_exty, g_extz, g_origx,
                    g_origy, g_origz, g_width, g_width, g_width * g_target_ani);
    write_vol("col_lo_fil.sit", g_width, origx, origy, origz, extx, exty, extz, g_phi_du);
    free_vect_and_zero_ptr(&g_phi_du);
    do_vect(&g_phi_du, g_fftw_nvox_r2c);
  } else write_vol("col_lo_fil.sit", g_width, g_origx, g_origy, g_origz, g_extx, g_exty, g_extz, g_phi_lo);

  /* compute correlation */
  printf("colores> Computing correlation between maps in direct space...\n");
  corr_hi_lo = 0;
  for (m = 0; m < g_nvox; m++) {
    corr_hi_lo += (*(g_phi_hi + m)) * (*(g_phi_lo + m));
  }
  corr_hi_lo *= 1.0 / (1.0 * g_nvox);
  printf("colores> Correlation with structure centered in density map: %15.7E\n", corr_hi_lo);

  /* FFT correlation test */
  printf("colores> Computing correlation in Fourier space...\n");
  normalize(g_phi_hi, g_nvox, (double) g_nvox);
  /* [FFT(A)] */
  fftw_execute(g_fftw_plan_fwd_lo);
  /* [FFT(B)] */
  fftw_execute(g_fftw_plan_fwd_hi);
  /* C [IFFT(BxA*)] */
  for (q = 0; q < g_fftw_nvox_c2r; q++) {
    save_b_re = c_re(g_fftw_grid_b[q]);
    c_re(g_fftw_grid_b[q]) = (c_re(g_fftw_grid_a[q]) * c_re(g_fftw_grid_b[q]) +
                              c_im(g_fftw_grid_a[q]) * c_im(g_fftw_grid_b[q])) * g_fftw_scale;
    c_im(g_fftw_grid_b[q]) = (c_re(g_fftw_grid_a[q]) * c_im(g_fftw_grid_b[q]) -
                              c_im(g_fftw_grid_a[q]) * save_b_re) * g_fftw_scale;
  }

  /* compute inverse FFT */
  zero_vect(g_phi_hi, g_fftw_nvox_r2c);
  fftw_execute(g_fftw_plan_rev_hi);
  printf("colores> FFT correlation with structure centered in density map: %15.7E\n", *(g_phi_hi + 0));


  /*============= EULER ANGLE GENERATION AND TIME ESTIMATE ===================*/

  draw_line();

  /* Euler angle management */

  printf("colores> Getting Euler angles.\n");
  g_eulers_count = 0;
  switch (g_euler_mode) {
    case 0:
      eu_proportional(g_eu_range, g_delta_rot, &g_eulers_count, &g_eulers);
      break;
    case 1:
      eu_sparsed(g_eu_range, g_delta_rot, &g_eulers_count, &g_eulers);
      break;
    case 2:
      eu_spiral(g_eu_range, g_delta_rot, &g_eulers_count, &g_eulers);
      break;
    case 3:
      read_eulers(g_eu_in_file, &g_eulers_count, &g_eulers);
      break;
  }
  write_eulers("col_eulers.dat", g_eulers_count, g_eulers, g_eu_range, g_delta_rot);
  printf("colores> Total number of orientations sampled: %ld\n", g_eulers_count);
  printf("colores> Euler angles saved in file col_eulers.dat.\n");

  /* initialize scoring hash table */
  g_hash_sav = (SAV *) alloc_vect(g_nvox, sizeof(SAV));
  for (m = 0; m < g_nvox; m++)  {
    (*(g_hash_sav + m)).eu = 0;
    (*(g_hash_sav + m)).score = 0.0f;
  }

  draw_line();

  /* timing of one FFT and computation of [FFT(A)]* */
  itime = get_the_time();
  fftw_execute(g_fftw_plan_fwd_lo);

  /* [FFT(A)]* */
  etime = get_the_time();
  printf("colores> Time of one FFT calculation: %s\n",
         smart_sprint_time(time_to_sec(time_diff(etime, itime))));

  /* estimating time for full 6D search */
  g_outfile = fopen("col_rotate.log", "w");
  if (g_outfile == NULL) {
    error_open_filename(80150, g_program, "col_rotate.log");
  }
  /* perform the 6D search 5 times selecting randomly the triplets of Euler angle */
  trial_time = 0;
  for (i = 0; i < 5; i++) {
    itime = get_the_time();
    j = (unsigned)(g_eulers_count * genrand());
    if (fmod(j + 1, 100) == 0.0) j--; /* avoid progress output */
    search6d_fft(j);
    etime = get_the_time();
    trial_time += (double)time_to_sec(time_diff(etime, itime));
  }
  fclose(g_outfile);
  trial_time /= 5.0;
  printf("colores> Average time spent on each rotation: %s\n", smart_sprint_time(trial_time));
  printf("colores> Estimated time for full 6D (on-lattice) search: %s\n",
         smart_sprint_time(trial_time * (g_eulers_count * 1.0) / (g_p_nprocs * 1.0)));
  if (g_pow_mode) printf("colores> Off-lattice Powell optimization will take significant extra time.\n");



  /*=================== 6D SEARCH WITH 3D FFT ON LATTICE =====================*/

  draw_line();
  printf("colores> Starting 6D on-lattice search with 3D FFT scan of Euler angles.\n");
  itime = get_the_time();

  g_outfile = fopen("col_rotate.log", "w");
  if (g_outfile == NULL) {
    error_open_filename(80210, g_program, "col_rotate.log");
  }
  fprintf(g_outfile, "# This file contains parameters of best translational fits for each rotation.\n");
  fprintf(g_outfile, "# May be useful for debugging. Note that correlation values are not yet normalized. \n");
  fprintf(g_outfile, "# Shows Euler angles in degrees, offset of correlation peak from \n");
  fprintf(g_outfile, "# reference center (%.3f,%.3f,%.3f) in A,\n", g_center_map[0], g_center_map[1], g_center_map[2]);
  fprintf(g_outfile, "# correlation at the center, correlation at the peak.\n");
  fprintf(g_outfile, "# \n");
  fprintf(g_outfile, "# Psi    Theta    Phi       X        Y        Z        Corr(0,0,0)       Corr(X,Y,Z)\n");

  /* Starting 6D on-lattice search with FFT */
#ifdef _SMP_

  printf("colores> Searching using %d processors\n", g_p_nprocs);

  /* Figure out which thread gets what */

  for (thread = 0; thread < g_p_nprocs; thread++) {
    threads_workload[thread] = g_eulers_count / g_p_nprocs;
  }
  for (thread = 0; thread < g_eulers_count - g_p_nprocs * (g_eulers_count / g_p_nprocs); thread++) {
    threads_workload[thread]++;
  }

  /* Set up data for individual threads */

  i = 0;
  for (thread = 0; thread < g_p_nprocs; thread++) {
    thread_fft_data_array[thread].thread_id = thread;
    thread_fft_data_array[thread].idx_start = i;
    i += threads_workload[thread];
    thread_fft_data_array[thread].idx_count = i - thread_fft_data_array[thread].idx_start;
    thread_fft_data_array[thread].results = (FFT_RESULTS *) alloc_vect(thread_fft_data_array[thread].idx_count, sizeof(FFT_RESULTS));
  }

  /* Initialize mutexes and semaphore */

  rc = pthread_mutex_init(&g_p_fft_score_mutex, NULL);
  if (rc != 0) {
    error("colores> Error: init failed for g_p_fft_score_mutex\n");
  }
  rc = pthread_mutex_init(&g_p_fft_log_mutex, NULL);
  if (rc != 0) {
    error("colores> Error: init failed for g_p_fft_log_mutex\n");
  }
  rc = pthread_mutex_init(&g_p_fft_status_mutex, NULL);
  if (rc != 0) {
    error("colores> Error: init failed for g_p_fft_status_mutex\n");
  }
  rc = pthread_mutex_init(&g_p_fft_plan_mutex, NULL);
  if (rc != 0) {
    error("colores> Error: init failed for g_p_fft_plan_mutex\n");
  }
  rc = pthread_cond_init(&g_p_fft_status_cond, NULL);
  if (rc != 0) {
    error("colores> Error: init failed for g_p_fft_status_cond\n");
  }
  g_p_fft_status_count = 0;
  g_p_fft_status_pred  = 0;

  /* Spawn threads */

  for (thread = 0; thread < g_p_nprocs; thread++) {
    rc = pthread_create(&threads[thread], NULL,
                        (void *)&search6d_fft_par, (void *) &thread_fft_data_array[thread]);
    if (rc != 0) {
      error("colores> Failed to spawn thread %d\n", thread);
    }
  }

  /* Get status updates and inform user */

  while (g_p_fft_status_count < g_eulers_count) {

    rc = pthread_mutex_lock(&g_p_fft_status_mutex);
    if (rc != 0) {
      error("colores> Error: lock failed for g_p_fft_status_mutex\n");
    }
    while (g_p_fft_status_pred == 0) {
      rc = pthread_cond_wait(&g_p_fft_status_cond, &g_p_fft_status_mutex);
      if (rc != 0) {
        error("colores> Error: wait failed for g_p_fft_status_cond\n");
      }
    }
    rc = pthread_mutex_unlock(&g_p_fft_status_mutex);
    if (rc != 0) {
      error("colores> Error: unlock failed for g_p_fft_status_mutex\n");
    }
    g_p_fft_status_pred--;

    p_fft_percent = ceil((g_p_fft_status_count * 100.0) / (g_eulers_count * 1.0));
    p_fft_progressbar[0] = 0;
    for (i = 0; i < 50; i++) {
      if (i < p_fft_percent / 2) {
        strcat(p_fft_progressbar, "#");
      } else {
        strcat(p_fft_progressbar, ".");
      }
    }

    printf("colores> |%s| %6d/%ld | %3d%% done\r", p_fft_progressbar, g_p_fft_status_count,
           g_eulers_count, p_fft_percent);
    fflush(NULL);

  }

  printf("\n");

  /* Wait for threads to finish */

  for (thread = 0; thread < g_p_nprocs; thread++) {
    rc = pthread_join(threads[thread], NULL);
    if (rc != 0) {
      error("colores> Error: pthread_join failed\n");
    }
  }

  /* Destroy mutexes */

  rc = pthread_mutex_destroy(&g_p_fft_score_mutex);
  if (rc != 0) {
    error("colores> Error: destroy failed for g_p_fft_score_mutex\n");
  }
  rc = pthread_mutex_destroy(&g_p_fft_log_mutex);
  if (rc != 0) {
    error("colores> Error: destroy failed for g_p_fft_log_mutex\n");
  }
  rc = pthread_mutex_destroy(&g_p_fft_status_mutex);
  if (rc != 0) {
    error("colores> Error: destroy failed for g_p_fft_status_mutex\n");
  }
  rc = pthread_mutex_destroy(&g_p_fft_plan_mutex);
  if (rc != 0) {
    error("colores> Error: destroy failed for g_p_fft_plan_mutex\n");
  }
  rc = pthread_cond_destroy(&g_p_fft_status_cond);
  if (rc != 0) {
    error("colores> Error: destroy failed for g_p_fft_status_cond\n");
  }

  /* Write output */

  for (thread = 0; thread < g_p_nprocs; thread++) {
    for (i = 0; i < thread_fft_data_array[thread].idx_count; i++) {
      fprintf(g_outfile, "%7.3f %7.3f %7.3f %8.3f %8.3f %8.3f   %15.7E   %15.7E\n",
              thread_fft_data_array[thread].results[i].psi,
              thread_fft_data_array[thread].results[i].theta,
              thread_fft_data_array[thread].results[i].phi,
              thread_fft_data_array[thread].results[i].x,
              thread_fft_data_array[thread].results[i].y,
              thread_fft_data_array[thread].results[i].z,
              thread_fft_data_array[thread].results[i].corr_orig,
              thread_fft_data_array[thread].results[i].corr_opt);
    }
  }

  /* Free memory */

  for (thread = 0; thread < g_p_nprocs; thread++) {
    free_vect_and_zero_ptr(&(thread_fft_data_array[thread].results));
  }

#else
  for (i = 0; i < g_eulers_count; i++) search6d_fft(i);
#endif

  etime = get_the_time();
  fclose(g_outfile);
  printf("colores> Actual time spent on 6D on-lattice search: %s\n",
         smart_sprint_time(time_to_sec(time_diff(etime, itime))));
  draw_line();

  /* free memory */
  free_vect_and_zero_ptr(&g_fftw_grid_a);
  free_vect_and_zero_ptr(&g_fftw_grid_b);
  fftw_destroy_plan(g_fftw_plan_fwd_lo);
  fftw_destroy_plan(g_fftw_plan_rev_lo);
  fftw_destroy_plan(g_fftw_plan_fwd_hi);
  fftw_destroy_plan(g_fftw_plan_rev_hi);
  fftw_destroy_plan(g_fftw_plan_fwd_du);

  /* Detect peaks */
  printf("colores> Translation function peak detection. \n");
  peak_extract(&g_hash_sav, &g_best_fit);
  draw_line();


  /*==================== POWELL OFF-LATTICE SEARCH ===========================*/


  if (g_pow_mode) {

    printf("colores> Off-lattice search (Powell's optimization method). \n");

    /* Powell initialization directions*/
    do_vect(&pow_init, 6);
    if (g_pow_delta_pos < 0 || g_pow_delta_ang < 0) {
      g_pow_delta_pos = g_width * 0.25;
      if (g_delta_rot > 40)
        g_pow_delta_ang = 10 * ROT_CONV; /* never larger than 10 degrees */
      else g_pow_delta_ang = g_delta_rot * ROT_CONV * 0.25;

    }
    for (i = 0; i < 3; ++i) {
      pow_init[i] = g_pow_delta_pos;
      pow_init[i + 3] = g_pow_delta_ang;
    }

    /* Determine best method to calculate correlation */

    if (g_pow_alg == 0) {

      test_pow_method = 1;

      /* Prepare scratch space */

      pow_scratch = (POWSCR *) alloc_vect(sizeof(POWSCR), 1);
      pow_scratch->pdb = (PDB *) alloc_vect(g_num_atoms, sizeof(PDB));
      for (q = 0; q < g_num_atoms; ++q) *(pow_scratch->pdb + q) = *(g_pdb_original + q);
      do_vect(&pow_scratch->phi_du, g_fftw_nvox_r2c);
      do_vect(&pow_scratch->phi_hi, g_fftw_nvox_r2c);

      printf("colores> Determining most efficient correlation algorithm based on convergence and time...\n");
      g_pow_alg = 1;
      itime = get_the_time();
      test_pow_corr_orig = powell_correlation(pow_init, pow_scratch);
      test_pow_corr_orig = powell_correlation(pow_init, pow_scratch);
      etime = get_the_time();
      test_pow_time_orig = (double)time_to_sec(time_diff(etime, itime)) / 2.0;
      printf("colores>    Original algorithm: Correlation = %10.8f  Time = %s\n",
             -test_pow_corr_orig,
             smart_sprint_time(test_pow_time_orig));
      g_pow_alg = 2;
      itime = get_the_time();
      test_pow_corr = powell_correlation(pow_init, pow_scratch);
      test_pow_corr = powell_correlation(pow_init, pow_scratch);
      etime = get_the_time();
      test_pow_time = (double)time_to_sec(time_diff(etime, itime)) / 2.0;
      printf("colores>    Masked algorithm:   Correlation = %10.8f  Time = %s\n",
             -test_pow_corr,
             smart_sprint_time(test_pow_time));
      if ((fabs(test_pow_corr - test_pow_corr_orig) < g_pow_tolerance) &&
          (test_pow_time < test_pow_time_orig)) {
        test_pow_time_orig = test_pow_time;
        test_pow_method    = 2;
      }
      g_pow_alg = 3;
      itime = get_the_time();
      test_pow_corr = powell_correlation(pow_init, pow_scratch);
      test_pow_corr = powell_correlation(pow_init, pow_scratch);
      etime = get_the_time();
      test_pow_time = (double)time_to_sec(time_diff(etime, itime)) / 2.0;
      printf("colores>    One-step algorithm: Correlation = %10.8f  Time = %s\n",
             -test_pow_corr,
             smart_sprint_time(test_pow_time));
      if ((fabs(test_pow_corr - test_pow_corr_orig) <= g_pow_tolerance) &&
          (test_pow_time < test_pow_time_orig)) {
        test_pow_time_orig = test_pow_time;
        test_pow_method    = 3;
      }

      g_pow_alg = test_pow_method;

      /* Clean up scratch space */

      free_vect_and_zero_ptr(&(pow_scratch->phi_hi));
      free_vect_and_zero_ptr(&(pow_scratch->phi_du));
      free_vect_and_zero_ptr(&(pow_scratch->pdb));
      free_vect_and_zero_ptr(&pow_scratch);

    }

    switch (g_pow_alg) {
      case 1:
        // Original three step code (scales based on g_nvox)
        printf("colores> Using original three-step correlation function.\n");
        break;
      case 2:
        // Masked three step (scales based on g_nvox)
        printf("colores> Using masked three-step correlation function.\n");
        break;
      case 3:
        // One step code (scales based on number of atoms in the high res structure)
        printf("colores> Using one-step correlation function.\n");
        break;
      default:
        fprintf(stderr, "colores> Error: Did not understand Powell correlation method, check -pwcorr\n");
        exit(1);
    }

    /* Powell initialization time */
    itime = get_the_time();

    g_outfile = fopen("col_powell.log", "w");
    if (g_outfile == NULL) {
      error_open_filename(80220, g_program, "col_powell.log");
    }
    fprintf(g_outfile, "# This file contains information about positions and ");
    fprintf(g_outfile, "Euler angles during the Powell search.\n");
    fprintf(g_outfile, "# The Euler angles are printed in degrees.\n");
    if (g_target_ani != 1.0) {
      fprintf(g_outfile, "# The coordinates in the z-direction are compressed by the anisotropy factor %.3f\n",
              g_target_ani);
    }
    if (g_p_nprocs > 1) {
      fprintf(g_outfile, "# Parallel Powell optimization: The order of maxima is not preserved in the log file.\n");
    }
    fprintf(g_outfile, "\n");

    /* Call Powell */

#ifdef _SMP_
    if (g_p_nprocs > 1) {
      /* Initialize SMP parallel work queue */
      g_p_workq = workq_init(g_p_nprocs);
      /* Initialize printing mutex */
      rc = pthread_mutex_init(&g_p_pow_print_mutex, NULL);
      if (rc != 0) {
        error("colores> Error: init failed for g_p_pow_print_mutex\n");
      }
    }
#endif

    if (g_p_nprocs > 1) {
      printf("colores> Using %d processors in SMP mode.\n", g_p_nprocs);
      printf("colores> Parallel Powell optimization: The order of maxima is not preserved in the output.\n");
    }
    printf("colores> Shown are: offset (in A) from reference center (%.3f,%.3f,%.3f),\n",
           g_center_map[0], g_center_map[1], g_center_map[2]);
    printf("colores> Euler angles (in degrees), and correlation value.\n");
    if (g_target_ani != 1.0) {
      printf("colores> The coordinates in the z-direction are compressed by the anisotropy factor %.3f\n",
             g_target_ani);
    }
    printf("colores> \n");
    printf("colores> Performing optimizations...\n");
    printf("colores> \n");

    for (iter = 0; iter < g_num_explored; iter++) {

      pow_args = (POWARG *) alloc_vect(sizeof(POWARG), 1);
      pow_args->iter     = iter;
      pow_args->pow_init = pow_init;
#ifdef _SMP_
      if (g_p_nprocs > 1) {
        pow_args->print_mutex = g_p_pow_print_mutex;
      }
#endif

      if (g_p_nprocs < 2) {
        powell_optimization((void *) pow_args);
      } else {
#ifdef _SMP_
        workq_add(g_p_workq, powell_optimization, (void *) pow_args);
#else
        error("SMP code not compiled in!\n");
#endif
      }

    }

#ifdef _SMP_
    /* Destroy SMP parallel work queue */
    if (g_p_nprocs > 1) {
      workq_destroy(g_p_workq);
    }
#endif


    fclose(g_outfile);

    /* sort initial best fits */
    for (p = 0; p < g_num_explored; p++)
      for (m = p + 1; m < g_num_explored; m++)
        if (g_best_fit[m].score > g_best_fit[p].score)
          SWAPPING(g_best_fit[p], g_best_fit[m], FIT);

    /* print Powell time */
    etime = get_the_time();
    printf("colores> Powell optimization time (%d runs): %s\n", g_num_explored,
           smart_sprint_time(time_to_sec(time_diff(etime, itime))));

    /* check degeneracy */
    num_degenerate = 0;
    for (iter = 0; iter < g_num_explored - num_degenerate; iter++) {
      for (p = iter + 1; p < g_num_explored - num_degenerate; p++)
        if ((similar_eulers(g_best_fit[p].euler[0], g_best_fit[p].euler[1],
                            g_best_fit[p].euler[2], g_best_fit[iter].euler[0],
                            g_best_fit[iter].euler[1], g_best_fit[iter].euler[2])) &&
            (((g_best_fit[p].pos[0] - g_best_fit[iter].pos[0]) *
              (g_best_fit[p].pos[0] - g_best_fit[iter].pos[0]) +
              (g_best_fit[p].pos[1] - g_best_fit[iter].pos[1]) *
              (g_best_fit[p].pos[1] - g_best_fit[iter].pos[1]) +
              (g_best_fit[p].pos[2] - g_best_fit[iter].pos[2]) *
              (g_best_fit[p].pos[2] - g_best_fit[iter].pos[2])) < (0.25 * g_width * g_width))) {
          num_degenerate++;
          SWAPPING(g_best_fit[p], g_best_fit[g_num_explored - num_degenerate], FIT);
          p--;
        }
    }

    g_num_explored -= num_degenerate;
    if (num_degenerate > 0) printf("colores> Removing %d redundant fits, keeping %d unique fits.\n", num_degenerate, g_num_explored);

    /* sort remaining best fits */
    for (p = 0; p < g_num_explored; p++)
      for (m = p + 1; m < g_num_explored; m++)
        if (g_best_fit[m].score > g_best_fit[p].score)
          SWAPPING(g_best_fit[p], g_best_fit[m], FIT);

    draw_line();

  } /* if(g_pow_mode) */


  /*============= RENORMALIZE AND WRITE TRANSLATION FUNCTION =================*/

  /* renormalize all scores */
  for (m = 0; m < g_nvox; m++)(*(g_hash_sav + m)).score /= g_best_fit[0].score;
  printf("colores> Renormalizing correlation values by highest score.\n");

  double *angle_map;
  angle_map = alloc_vect(g_extx * g_exty * g_extz, sizeof(double));

  /* write translation function and store pre-Powell lattice */
  /* note that translations need to be inverted, see search6d_fft() */
  for (m = 0; m < g_nvox; m++)  { /* invert translation function */
    indv = m;
    indz = indv / (g_extx * g_exty);
    indv -= indz * (g_extx * g_exty);
    indy = indv / g_extx;
    indx = indv - indy * g_extx;
    indx = g_extx - indx - 1;
    indy = g_exty - indy - 1;
    indz = g_extz - indz - 1;
    indv = indx + indy * g_extx + indz * g_extx * g_exty;
    *(g_phi_hi + indv) = (*(g_hash_sav + m)).score;
    *(angle_map + indv) = (double)(*(g_hash_sav + m)).eu;
  }

  /* interpolate translation function and export */
  printf("colores> Writing translation function lattice to density file in Situs format.\n");
  if (g_target_ani != 1.0) {
    free_vect_and_zero_ptr(&g_phi_du);
    interpolate_map(&g_phi_du, &extx, &exty, &extz, &origx, &origy, &origz,
                    g_width, g_width, g_width, g_phi_hi, g_extx, g_exty, g_extz, g_origx,
                    g_origy, g_origz, g_width, g_width, g_width * g_target_ani);
    write_vol("col_trans.sit", g_width, origx, origy, origz, extx, exty, extz, g_phi_du);
    free_vect_and_zero_ptr(&g_phi_du);
    do_vect(&g_phi_du, g_fftw_nvox_r2c);
  } else write_vol("col_trans.sit", g_width, g_origx, g_origy, g_origz, g_extx, g_exty, g_extz, g_phi_hi);

  /* export angle index volume for sculptor */
  if (g_sculptor) write_vol("col_rot.sit", g_width, g_origx, g_origy, g_origz, g_extx, g_exty, g_extz, angle_map);
  free_vect_and_zero_ptr(&angle_map);

  /* writing master .eli file for sculptor */
  if (g_sculptor) {
    if ((out = fopen("col_exh_search.eli", "w")) == NULL) {
      printf("colores> Could not open file col_exh_search.eli for writing\n");
    }
    fprintf(out, "col_probe.pdb\n");
    fprintf(out, "col_target.sit\n");
    fprintf(out, "col_eulers.dat\n");
    fprintf(out, "col_trans.sit\n");
    fprintf(out, "col_rot.sit\n");
    fprintf(out, "ROT       0.0000     0.0000     0.0000\n");
    fprintf(out, "TRANS     0.0000     0.0000     0.0000\n");
    fclose(out);
  }

  /* write info on optimal rotations and their scores as a function of */
  /* position; note that rotations are not defined for zero score! */
  printf("colores> Writing translation function lattice information to log file.\n");
  g_outfile = fopen("col_trans.log", "w");
  if (g_outfile == NULL) {
    error_open_filename(80230, g_program, "col_trans.log");
  }
  fprintf(g_outfile, "# This file contains correlation data as a function of translation.\n");
  fprintf(g_outfile, "# Note that correlation values are normalized by globally optimal fit. \n");
  fprintf(g_outfile, "# Shows Euler angles in degrees, offset of correlation peak from \n");
  fprintf(g_outfile, "# reference center (%.3f,%.3f,%.3f) in A,\n", g_center_map[0], g_center_map[1], g_center_map[2]);
  fprintf(g_outfile, "# correlation at offset value.\n");
  fprintf(g_outfile, "# Psi    Theta    Phi       X        Y        Z        Corr(X,Y,Z)\n");
  fprintf(g_outfile, "# \n");
  for (indz = 0; indz < g_extz; indz++)
    for (indy = 0; indy < g_exty; indy++)
      for (indx = 0; indx < g_extx; indx++) {
        q = g_extx * g_exty * indz + g_extx * indy + indx;
        if (g_hash_sav[q].score > 0.0f)
          fprintf(g_outfile, "%7.3f %7.3f %7.3f %8.3f %8.3f %8.3f   %15.7E\n",
                  *(g_eulers + g_hash_sav[q].eu * 3 + 0) / ROT_CONV,
                  *(g_eulers + g_hash_sav[q].eu * 3 + 1) / ROT_CONV,
                  *(g_eulers + g_hash_sav[q].eu * 3 + 2) / ROT_CONV,
                  /* note that translations need to be inverted, see search6d_fft() */
                  g_width * g_extx_half - g_width * indx,
                  g_width * g_exty_half - g_width * indy,
                  (g_width * g_extz_half - g_width * indz)*g_target_ani,
                  g_hash_sav[q].score);
      }
  fclose(g_outfile);
  free_vect_and_zero_ptr(&g_hash_sav);

  /*=========================== SAVE BEST RESULTS ============================*/

  draw_line();
  printf("colores> Saving the best results.\n");

  g_pdb_save = (PDB *) alloc_vect(g_num_atoms, sizeof(PDB));

  /* now write the PDB files of the best fits */
  for (m = 0; m < g_num_explored; m++) {
    for (p = 0; p < g_num_atoms; ++p) *(g_pdb_save + p) = *(g_pdb_original + p);

    sprintf(out_string, "col_best_%03lu.pdb", m + 1);

    if (g_peak_sharp == 1) {
      printf("colores> Estimating peak sharpness and writing best fit no. %3ld to file %s. \n", m + 1, out_string);
    } else {
      printf("colores> Writing best fit no. %3ld to file %s. \n", m + 1, out_string);
    }

    g_outfile = fopen(out_string, "w");
    if (g_outfile == NULL) {
      error_open_filename(80320, g_program, out_string);
    }
    fprintf(g_outfile, "REMARK\nREMARK    File name %s\n", out_string);
    fprintf(g_outfile, "REMARK    Low-resolution fit of structure %s into map %s\n", argv[2], argv[1]);
    fprintf(g_outfile, "REMARK    \n");
    fprintf(g_outfile, "REMARK    Resolution anisotropy factor (Z vs. XY): %f\n", g_target_ani);
    fprintf(g_outfile, "REMARK    Resolution: %f Angstrom; density cutoff: %f\n", g_target_res, g_low_cutoff);
    fprintf(g_outfile, "REMARK    Angular step %f degrees, %ld Euler angles\n", g_delta_rot, g_eulers_count);
    fprintf(g_outfile, "%s\n", corr_string);
    fprintf(g_outfile, "REMARK    Unnormalized correlation coefficient: %f\n", g_best_fit[m].score);
    fprintf(g_outfile, "REMARK    Normalized correlation coefficient: %f\n", g_best_fit[m].score / g_best_fit[0].score);
    fprintf(g_outfile, "REMARK    (normalized with respect to the global best fit in file col_best1.pdb)\n");
    fprintf(g_outfile, "REMARK    \n");
    fprintf(g_outfile, "REMARK    Euler angles (Psi, Theta, Phi):  %7.3f %7.3f %7.3f deg.\n",
            g_best_fit[m].euler[0] / ROT_CONV,
            g_best_fit[m].euler[1] / ROT_CONV,
            g_best_fit[m].euler[2] / ROT_CONV);
    fprintf(g_outfile, "REMARK    Center position (X,Y,Z):            %7.3f %7.3f %7.3f A\n",
            g_best_fit[m].pos[0] + g_center_map[0],
            g_best_fit[m].pos[1] + g_center_map[1],
            g_target_ani * g_best_fit[m].pos[2] + g_center_map[2]);
    fprintf(g_outfile, "REMARK    \n");

    rot_euler(g_pdb_original, g_pdb_save, g_num_atoms, g_best_fit[m].euler[0], g_best_fit[m].euler[1], g_best_fit[m].euler[2]);

    if (g_sculptor) {
      if ((out = fopen("col_exh_search.eli", "a")) == NULL) {
        printf("colores> Could not open file col_exh_search.eli for appending\n");
      }
      fprintf(out, "CANDIDATE %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f\n",
              g_best_fit[m].pos[0], g_best_fit[m].pos[1], g_best_fit[m].pos[2],
              g_best_fit[m].euler[0] / ROT_CONV,
              g_best_fit[m].euler[1] / ROT_CONV,
              g_best_fit[m].euler[2] / ROT_CONV);
      fclose(out);
    }

    if (g_peak_sharp == 1) {

      fprintf(g_outfile, "REMARK    The following table contains the angular variability of the \n");
      fprintf(g_outfile, "REMARK    correlation values about the maximum.\n");
      fprintf(g_outfile, "REMARK    Angle range: +-%5.3f deg., %4.2f deg steps.\n", DEF_RANGE, DEF_STEP);
      fprintf(g_outfile, "REMARK    Note: Correlation values in the table are normalized by the maximum\n");
      fprintf(g_outfile, "REMARK    value in the current fit, not the global best fit!\n");
      fprintf(g_outfile, "REMARK    \n");
      fprintf(g_outfile, "REMARK    Angle  C(rotX) C(rotY) C(rotZ)\n");

      /* calculate angular variability around the XYZ axes */

      /* X axis */
      for (deflection_deg = -DEF_RANGE; deflection_deg < DEF_RANGE; deflection_deg += DEF_STEP) {
        deflection_rad = deflection_deg * ROT_CONV;
        rot_axis(g_pdb_save, g_pdb_move, g_num_atoms, 'X', deflection_rad);
        curr_shift[0] = g_best_fit[m].pos[0];
        curr_shift[1] = g_best_fit[m].pos[1];
        curr_shift[2] = g_best_fit[m].pos[2] * g_target_ani;
        project_mass(&g_phi_hi, g_nvox, g_width, g_width, g_width * g_target_ani, g_extx, g_exty, g_extz, g_pdb_move, g_num_atoms, curr_shift, g_ignored);
        convolve_kernel_inside_fast(&g_phi_du, g_phi_hi, g_extx, g_exty, g_extz, g_phi_fx, g_ext_fx, g_norm_hi, g_ignored);
        cp_vect(&g_phi_hi, &g_phi_du, g_nvox);
        corr_hi_lo = 0;
        for (n = 0; n < g_nvox; n++) { /* compute correlation */
          corr_hi_lo += (*(g_phi_hi + n)) * (*(g_phi_lo + n));
        }
        corr_hi_lo /= (double)g_nvox;
        fprintf(g_outfile, "REMARK    %7.3f %1.5f", deflection_deg, corr_hi_lo / g_best_fit[m].score);

        /* Y axis */
        rot_axis(g_pdb_save, g_pdb_move, g_num_atoms, 'Y', deflection_rad);
        curr_shift[0] = g_best_fit[m].pos[0];
        curr_shift[1] = g_best_fit[m].pos[1];
        curr_shift[2] = g_best_fit[m].pos[2] * g_target_ani;
        project_mass(&g_phi_hi, g_nvox, g_width, g_width, g_width * g_target_ani, g_extx, g_exty, g_extz, g_pdb_move, g_num_atoms, curr_shift, g_ignored);
        convolve_kernel_inside_fast(&g_phi_du, g_phi_hi, g_extx, g_exty, g_extz, g_phi_fx, g_ext_fx, g_norm_hi, g_ignored);
        cp_vect(&g_phi_hi, &g_phi_du, g_nvox);
        corr_hi_lo = 0;
        for (n = 0; n < g_nvox; n++) { /* compute correlation */
          corr_hi_lo += (*(g_phi_hi + n)) * (*(g_phi_lo + n));
        }
        corr_hi_lo /= (double)g_nvox;
        fprintf(g_outfile, " %1.5f ",  corr_hi_lo / g_best_fit[m].score);

        /* Z axis */
        rot_axis(g_pdb_save, g_pdb_move, g_num_atoms, 'Z', deflection_rad);
        curr_shift[0] = g_best_fit[m].pos[0];
        curr_shift[1] = g_best_fit[m].pos[1];
        curr_shift[2] = g_best_fit[m].pos[2] * g_target_ani;
        project_mass(&g_phi_hi, g_nvox, g_width, g_width, g_width * g_target_ani, g_extx, g_exty, g_extz, g_pdb_move, g_num_atoms, curr_shift, g_ignored);
        convolve_kernel_inside_fast(&g_phi_du, g_phi_hi, g_extx, g_exty, g_extz, g_phi_fx, g_ext_fx, g_norm_hi, g_ignored);
        cp_vect(&g_phi_hi, &g_phi_du, g_nvox);
        corr_hi_lo = 0;
        for (n = 0; n < g_nvox; n++) { /* compute correlation */
          corr_hi_lo += (*(g_phi_hi + n)) * (*(g_phi_lo + n));
        }
        corr_hi_lo /= (double)g_nvox;
        fprintf(g_outfile, "%1.5f\n", corr_hi_lo / g_best_fit[m].score);
      }

    }

    fprintf(g_outfile, "REMARK    \n");
    fclose(g_outfile);

    for (p = 0; p < g_num_atoms; ++p) *(g_pdb_save + p) = *(g_pdb_original + p);
    rot_euler(g_pdb_original, g_pdb_save, g_num_atoms, g_best_fit[m].euler[0], g_best_fit[m].euler[1], g_best_fit[m].euler[2]);
    translate(g_pdb_save, g_pdb_move, g_num_atoms, g_best_fit[m].pos[0] + g_center_map[0], g_best_fit[m].pos[1] + g_center_map[1], g_target_ani * g_best_fit[m].pos[2] + g_center_map[2]);
    append_pdb(out_string, g_num_atoms, g_pdb_move);
  }

  /* free memory */

  free_vect_and_zero_ptr(&g_phi_hi);
  free_vect_and_zero_ptr(&g_phi_lo);
  free_vect_and_zero_ptr(&g_phi_du);
  free_vect_and_zero_ptr(&g_pdb_original);
  free_vect_and_zero_ptr(&g_pdb_move);
  free_vect_and_zero_ptr(&g_pdb_save);
  free_vect_and_zero_ptr(&g_eulers);
  free_vect_and_zero_ptr(&g_inside_list);
  free_vect_and_zero_ptr(&g_inside_list_flipped);

  draw_line();
  printf("colores> Output files:\n");
  printf("   col_best*.pdb      => Best docking results in PDB format with info in header\n");
  printf("   col_eulers.dat     => colores-readable list of Euler angles\n");
  printf("   col_rotate.log     => Rotation function (unnormalized) log file \n");
  printf("   col_trans.log      => Translation function (norm. by best fit) log file\n");
  printf("   col_trans.sit      => Translation function (norm. by best fit) in Situs format\n");
  printf("   col_lo_fil.sit     => Filtered target volume in Situs format, just prior to correlation calculation \n");
  printf("   col_hi_fil.sit     => Filtered (and centered) probe structure in Situs format, just prior to correlation calculation \n");
  if (g_pow_mode) printf("   col_powell.log     => Powell optimization log file\n");
  if (g_sculptor) {
    printf("   col_exh_search.eli => Master file for interactive exploration with Sculptor\n");
    printf("   col_rot.sit        => Angle index map for interactive exploration with Sculptor, in Situs format\n");
    printf("   col_target.sit     => Target volume for interactive exploration with Sculptor, in Situs format\n");
  }

  draw_line();

  printf("colores> All done!\n");
  return 0;
}













/*========================== FUNCTIONS ===============================*/






/*====================================================================*/
static void search6d_fft(unsigned i)
{
  /* carries out a FFT scan of all Euler angles; requires all global variables to be set */

  double curr_score;
  unsigned long q, m;
  double save_b_re;
  double max_score;
  double rx, ry, rz;
  double zero_shift[3] = {0.0, 0.0, 0.0};
  unsigned grid_pos[3] = {0.0, 0.0, 0.0};

  rot_euler(g_pdb_original, g_pdb_move, g_num_atoms, *(g_eulers + 3 * i + 0), *(g_eulers + 3 * i + 1), *(g_eulers + 3 * i + 2));
  project_mass(&g_phi_hi, g_nvox, g_width, g_width, g_width * g_target_ani, g_extx, g_exty, g_extz, g_pdb_move, g_num_atoms, zero_shift, g_ignored);

  convolve_kernel_inside_fast(&g_phi_du, g_phi_hi, g_extx, g_exty, g_extz, g_phi_fx, g_ext_fx, g_norm_hi * g_nvox, g_ignored);

  /* [FFT(B)] */
  fftw_execute(g_fftw_plan_fwd_du);

  /* C [IFFT(BxA*)] */
  for (q = 0; q < g_fftw_nvox_c2r; q++) {
    save_b_re = c_re(g_fftw_grid_b[q]);
    c_re(g_fftw_grid_b[q]) = (c_re(g_fftw_grid_a[q]) * c_re(g_fftw_grid_b[q]) +
                              c_im(g_fftw_grid_a[q]) * c_im(g_fftw_grid_b[q])) * g_fftw_scale;
    c_im(g_fftw_grid_b[q]) = (c_re(g_fftw_grid_a[q]) * c_im(g_fftw_grid_b[q + 1]) -
                              c_im(g_fftw_grid_a[q]) * save_b_re) * g_fftw_scale;
  }

  /* assign the translation function to g_phi_hi */
  zero_vect(g_phi_hi, g_fftw_nvox_r2c);
  fftw_execute(g_fftw_plan_rev_hi);

  /* update the highest score found so far for each grid translation, and save corresponding Euler angle */
  /* note that g_phi_hi has the quadrants flipped (FFT effect).                                          */

  max_score = -1e20;
  for (m = 0; m < g_inside_num_flipped; m++) {
    curr_score = (*(g_phi_hi + g_inside_list_flipped[m].ifft));
    q = g_inside_list_flipped[m].ireal;
    if (curr_score > g_hash_sav[q].score) {
      g_hash_sav[q].score = curr_score;
      g_hash_sav[q].eu = i;
    }
    if (curr_score > max_score) {
      max_score = curr_score;
      grid_pos[0] = g_inside_list_flipped[m].ix;
      grid_pos[1] = g_inside_list_flipped[m].iy;
      grid_pos[2] = g_inside_list_flipped[m].iz;
    }
  }

  /* write rotation function i.e. highest translational peak found for each rotation                     */
  /* Note that offset from center is negative (see definiton of C(T) in Wriggers & Chacon, 2001,         */
  /* Structure 9:779-788, Fig. 3.).                                                                      */
  rx = g_width * g_extx_half - g_width * grid_pos[0];
  ry = g_width * g_exty_half - g_width * grid_pos[1];
  rz = (g_width * g_extz_half - g_width * grid_pos[2]) * g_target_ani;

  fprintf(g_outfile, "%7.3f %7.3f %7.3f %8.3f %8.3f %8.3f   %15.7E   %15.7E\n"
          , *(g_eulers + 3 * i + 0) / ROT_CONV, *(g_eulers + 3 * i + 1) / ROT_CONV, *(g_eulers + 3 * i + 2) / ROT_CONV,
          rx, ry, rz, *(g_phi_hi + 0), max_score);

  /* keep user informed about progress */
  if (fmod(i + 1, 100) == 0.0)  printf("colores> %d of %ld Euler angles processed with 3D FFT.\n", i + 1, g_eulers_count);
}


/*====================================================================*/
#ifdef _SMP_
static void *search6d_fft_par(void *thread_arg)
{

  FFT_THREAD_DATA *my_data;

  /* Define local variables */

  double *p_phi_hi;                     /* high resolution map */
  double *p_phi_du;                     /* dummy map */
  fftw_complex *p_fftw_grid_b;          /* B structure factors */
  PDB *p_pdb_move;                      /* PDB coordinates */
  fftw_plan p_fftw_plan_forward;        /* forward FFT plan setup */
  fftw_plan p_fftw_plan_reverse;        /* reverse FFT plan setup */

  int    idx_start, idx_count;
  int    rc;
  FFT_RESULTS *results;

  int    i, idx_results;
  double curr_score;
  unsigned long q, m;
  double save_b_re;
  double max_score;
  double rx, ry, rz;
  double zero_shift[3] = {0.0, 0.0, 0.0};
  unsigned grid_pos[3] = {0, 0, 0};

  /* Startup */

  my_data = (FFT_THREAD_DATA *) thread_arg;
  idx_start = my_data->idx_start;
  idx_count = my_data->idx_count;
  results   = my_data->results;

  idx_results = 0;

  /* Allocate local variables */

  do_vect(&p_phi_du, g_fftw_nvox_r2c);
  do_vect(&p_phi_hi, g_fftw_nvox_r2c);

  p_pdb_move = (PDB *) alloc_vect(g_num_atoms, sizeof(PDB));
  for (q = 0; q < g_num_atoms; ++q)
    *(p_pdb_move + q) = *(g_pdb_original + q);

  p_fftw_grid_b = (fftw_complex *) alloc_vect(g_fftw_nvox_c2r, sizeof(fftw_complex));

  if (pthread_mutex_lock(&g_p_fft_plan_mutex)) {
    error("colores> Error: lock failed for g_p_fft_plan_mutex\n");
  }
  p_fftw_plan_forward = fftw_plan_dft_r2c_3d(g_extz, g_exty, g_extx,
                        p_phi_du, (fftw_complex *) p_fftw_grid_b,
                        FFTW_ESTIMATE);
  p_fftw_plan_reverse = fftw_plan_dft_c2r_3d(g_extz, g_exty, g_extx,
                        (fftw_complex *) p_fftw_grid_b, p_phi_hi,
                        FFTW_ESTIMATE);
  if (pthread_mutex_unlock(&g_p_fft_plan_mutex)) {
    error("colores> Error: unlock failed for g_p_fft_plan_mutex\n");
  }

  /* Do actual work */

  for (i = idx_start; i < idx_start + idx_count; i++) {

    rot_euler(g_pdb_original, p_pdb_move, g_num_atoms,
              *(g_eulers + 3 * i + 0), *(g_eulers + 3 * i + 1), *(g_eulers + 3 * i + 2));
    project_mass(&p_phi_hi, g_nvox, g_width, g_width, g_width * g_target_ani,
                 g_extx, g_exty, g_extz, p_pdb_move, g_num_atoms, zero_shift, g_ignored);
    convolve_kernel_inside_fast(&p_phi_du, p_phi_hi, g_extx, g_exty, g_extz,
                                g_phi_fx, g_ext_fx, g_norm_hi * g_nvox, g_ignored);

    /* [FFT(B)] */
    fftw_execute(p_fftw_plan_forward);

    /* C [IFFT(BxA*)] */
    for (q = 0; q < g_fftw_nvox_c2r; q++) {
      save_b_re = c_re(p_fftw_grid_b[q]);
      c_re(p_fftw_grid_b[q]) = (c_re(g_fftw_grid_a[q]) * c_re(p_fftw_grid_b[q]) +
                                c_im(g_fftw_grid_a[q]) * c_im(p_fftw_grid_b[q])) * g_fftw_scale;
      c_im(p_fftw_grid_b[q]) = (c_re(g_fftw_grid_a[q]) * c_im(p_fftw_grid_b[q]) -
                                c_im(g_fftw_grid_a[q]) * save_b_re) * g_fftw_scale;
    }

    /* assign the translation function to p_phi_hi */

    zero_vect(p_phi_hi, g_fftw_nvox_r2c);

    fftw_execute(p_fftw_plan_reverse);

    /* update the highest score found so far for each grid translation,
       and save corresponding Euler angle.
       Note that p_phi_hi has the quadrants flipped (FFT effect). */

    max_score = -1e20;
    for (m = 0; m < g_inside_num_flipped; m++) {
      curr_score = (*(p_phi_hi + g_inside_list_flipped[m].ifft));
      q = g_inside_list_flipped[m].ireal;
      if (curr_score > g_hash_sav[q].score) { /* Check first */
        rc = pthread_mutex_lock(&g_p_fft_score_mutex);
        if (rc != 0) {
          error("colores> Error: lock failed for g_p_fft_score_mutex\n");
        }
        if (curr_score > g_hash_sav[q].score) { /* Check again since it might have changed */
          g_hash_sav[q].score = curr_score;
          g_hash_sav[q].eu    = i;
        }
        rc = pthread_mutex_unlock(&g_p_fft_score_mutex);
        if (rc != 0) {
          error("colores> Error: unlock failed for g_p_fft_score_mutex\n");
        }
      }
      if (curr_score > max_score) {
        max_score = curr_score;
        grid_pos[0] = g_inside_list_flipped[m].ix;
        grid_pos[1] = g_inside_list_flipped[m].iy;
        grid_pos[2] = g_inside_list_flipped[m].iz;
      }
    }

    /* write rotation function i.e. highest translational peak found for each rotation */
    /* Note that offset from center is negative (see definiton of C(T) in Wriggers &   */
    /* Chacon, 2001, Structure 9:779-788, Fig. 3.).                                    */

    rx = g_width * g_extx_half - g_width * grid_pos[0];
    ry = g_width * g_exty_half - g_width * grid_pos[1];
    rz = (g_width * g_extz_half - g_width * grid_pos[2]) * g_target_ani;

    /* Write results back */

    results[idx_results].psi       = *(g_eulers + 3 * i + 0) / ROT_CONV;
    results[idx_results].theta     = *(g_eulers + 3 * i + 1) / ROT_CONV;
    results[idx_results].phi       = *(g_eulers + 3 * i + 2) / ROT_CONV;
    results[idx_results].x         = rx;
    results[idx_results].y         = ry;
    results[idx_results].z         = rz;
    results[idx_results].corr_orig = *(p_phi_hi + 0);
    results[idx_results].corr_opt  = max_score;

    idx_results++;

    /* keep user informed about progress */

    rc = pthread_mutex_lock(&g_p_fft_status_mutex);
    if (rc != 0) {
      error("colores> Error: lock failed for g_p_fft_status_mutex\n");
    }
    g_p_fft_status_pred++;
    g_p_fft_status_count++;
    rc = pthread_cond_signal(&g_p_fft_status_cond);
    if (rc != 0) {
      error("colores> Error: signal failed for g_p_fft_status_cond\n");
    }
    rc = pthread_mutex_unlock(&g_p_fft_status_mutex);
    if (rc != 0) {
      error("colores> Error: unlock failed for g_p_fft_status_mutex\n");
    }

  }

  /* Free memory */

  free_vect_and_zero_ptr(&p_phi_du);
  free_vect_and_zero_ptr(&p_phi_hi);
  free_vect_and_zero_ptr(&p_pdb_move);
  free_vect_and_zero_ptr(&p_fftw_grid_b);
  if (pthread_mutex_lock(&g_p_fft_plan_mutex)) {
    error("colores> Error: lock failed for g_p_fft_plan_mutex!\n");
  }
  fftw_destroy_plan(p_fftw_plan_forward);
  fftw_destroy_plan(p_fftw_plan_reverse);
  if (pthread_mutex_unlock(&g_p_fft_plan_mutex)) {
    error("colores> Error: unlock failed for g_p_fft_plan_mutex!\n");
  }

  pthread_exit(0);
}
#endif

/*====================================================================*/
/* Wrapper for Powell optimization (thread-safe)                      */
static void powell_optimization(void *args)
{

  POWARG        *arguments;
  unsigned long iter;
  double        *pow_init;
  int           pow_code;
  double        pow_vect6d[6];
  double        pow_score;
  POWSCR        *pow_scratch = NULL;
  POWRES        *pow_results = NULL;
  POWITERRES    *curr_iter = NULL;
  POWITERRES    *last_iter = NULL;
  int           p, q;
#ifdef _SMP_
  int           rc;
#endif

  /* Unpack arguments */

  arguments = (POWARG *) args;
  iter      = arguments->iter;
  pow_init  = arguments->pow_init;

  /* Assign 6D Powell vector */

  for (p = 0; p < 3; p++) {
    pow_vect6d[p] = g_best_fit[iter].pos[p];
    pow_vect6d[p + 3] = g_best_fit[iter].euler[p];
  }

  /* Prepare scratch space */

  pow_scratch = (POWSCR *) alloc_vect(sizeof(POWSCR), 1);
  pow_scratch->pdb = (PDB *) alloc_vect(g_num_atoms, sizeof(PDB));
  for (q = 0; q < g_num_atoms; ++q) *(pow_scratch->pdb + q) = *(g_pdb_original + q);
  do_vect(&pow_scratch->phi_du, g_fftw_nvox_r2c);
  do_vect(&pow_scratch->phi_hi, g_fftw_nvox_r2c);

  /* Prepare result struct */

  pow_results = (POWRES *) alloc_vect(sizeof(POWRES), 1);
  pow_results->head = NULL;
  pow_results->last = NULL;

  /* Call Powell maximization function */

  powell_r(&pow_code, &pow_score, pow_vect6d, 6,
           powell_correlation, g_pow_max_iter, g_outfile,
           g_pow_tolerance, pow_init,
           pow_scratch, pow_results);

  /* Clean up scratch space */

  free_vect_and_zero_ptr(&pow_scratch->phi_hi);
  free_vect_and_zero_ptr(&pow_scratch->phi_du);
  free_vect_and_zero_ptr(&pow_scratch->pdb);
  free_vect_and_zero_ptr(&pow_scratch);

  /* ===================== Serialized ====================== */

#ifdef _SMP_
  if (g_p_nprocs > 1) {
    rc = pthread_mutex_lock(&arguments->print_mutex);
    if (rc != 0) {
      error("colores> Error: lock failed for arguments->print_mutex failed\n");
    }
  }
#endif

  /* Output initial point to both screen and log file */

  printf("colores> Powell optimization for score maximum no. %2ld.\n", iter + 1);
  printf("colores>   X       Y       Z       Psi     Theta   Phi      Correlation \n");
  printf("colores> %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f    %10.7E Initial\n",
         g_best_fit[iter].pos[0],
         g_best_fit[iter].pos[1],
         g_best_fit[iter].pos[2]*g_target_ani,
         g_best_fit[iter].euler[0] / ROT_CONV,
         g_best_fit[iter].euler[1] / ROT_CONV,
         g_best_fit[iter].euler[2] / ROT_CONV,
         g_best_fit[iter].score);

  fprintf(g_outfile, "Exploring score maximum no. %2ld\n", iter + 1);
  fprintf(g_outfile, "   X       Y       Z       Psi     Theta   Phi      Correlation   Iteration\n");
  fprintf(g_outfile, " %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f    %10.7E 0\n",
          g_best_fit[iter].pos[0],
          g_best_fit[iter].pos[1],
          g_best_fit[iter].pos[2]*g_target_ani,
          g_best_fit[iter].euler[0] / ROT_CONV,
          g_best_fit[iter].euler[1] / ROT_CONV,
          g_best_fit[iter].euler[2] / ROT_CONV,
          g_best_fit[iter].score);

  /* Write detailed results of optimization steps to log file */

  curr_iter = pow_results->head;
  while (curr_iter != NULL) {
    printf("colores> %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f    %10.7E %d\n",
           curr_iter->res[0],
           curr_iter->res[1],
           curr_iter->res[2]*g_target_ani,
           curr_iter->res[3] / ROT_CONV,
           curr_iter->res[4] / ROT_CONV,
           curr_iter->res[5] / ROT_CONV,
           curr_iter->corr,
           curr_iter->iter);
    fprintf(g_outfile, " %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f    %10.7E %d\n",
            curr_iter->res[0],
            curr_iter->res[1],
            curr_iter->res[2]*g_target_ani,
            curr_iter->res[3] / ROT_CONV,
            curr_iter->res[4] / ROT_CONV,
            curr_iter->res[5] / ROT_CONV,
            curr_iter->corr,
            curr_iter->iter);
    curr_iter = curr_iter->next;
  }

  /* save and write final angles, positions, and score */
  for (p = 0; p < 3; p++) {
    g_best_fit[iter].pos[p] = pow_results->last->res[p];
    g_best_fit[iter].euler[p] = fmod(pow_results->last->res[p + 3] + 2 * PI, 2 * PI);
  }
  g_best_fit[iter].score = -pow_score;

  remap_eulers((g_best_fit[iter].euler + 0),
               (g_best_fit[iter].euler + 1),
               (g_best_fit[iter].euler + 2),
               g_best_fit[iter].euler[0],
               g_best_fit[iter].euler[1],
               g_best_fit[iter].euler[2],
               0.0, 0.0, 0.0);

  printf("colores> %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f    %10.7E Final\n",
         g_best_fit[iter].pos[0],
         g_best_fit[iter].pos[1],
         g_best_fit[iter].pos[2]*g_target_ani,
         g_best_fit[iter].euler[0] / ROT_CONV,
         g_best_fit[iter].euler[1] / ROT_CONV,
         g_best_fit[iter].euler[2] / ROT_CONV,
         g_best_fit[iter].score);
  printf("colores>\n");

  fprintf(g_outfile, " %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f    %10.7E Final\n",
          g_best_fit[iter].pos[0],
          g_best_fit[iter].pos[1],
          g_best_fit[iter].pos[2]*g_target_ani,
          g_best_fit[iter].euler[0] / ROT_CONV,
          g_best_fit[iter].euler[1] / ROT_CONV,
          g_best_fit[iter].euler[2] / ROT_CONV,
          g_best_fit[iter].score);

#ifdef _SMP_
  if (g_p_nprocs > 1) {
    rc = pthread_mutex_unlock(&arguments->print_mutex);
    if (rc != 0) {
      error("colores> Error: unlock failed for arguments->print_mutex failed\n");
    }
  }
#endif

  /* ===================== END ====================== */

  /* Clean up */

  fflush(NULL);

  curr_iter = pow_results->head;
  while (curr_iter != NULL) {
    last_iter = curr_iter;
    curr_iter = curr_iter->next;
    free_vect_and_zero_ptr(&(last_iter->res));
    free_vect_and_zero_ptr(&last_iter);
  }

  free_vect_and_zero_ptr(&args);

  return;
}

/*====================================================================*/
static double powell_correlation(double *vect6d, POWSCR *scratch)
{
  /* computes current correlation value for Powell based on 6D coordinates */
  /* Note: returns NEGATIVE value since we wish to maximize the score */

  rot_euler(g_pdb_original, scratch->pdb, g_num_atoms, vect6d[3], vect6d[4], vect6d[5]);
  scratch->curr_shift[0] = vect6d[0];
  scratch->curr_shift[1] = vect6d[1];
  scratch->curr_shift[2] = vect6d[2] * g_target_ani;
  scratch->corr_hi_lo = 0;

  switch (g_pow_alg) {
    case 1:
      // Original three step code (scales based on g_nvox)
      project_mass(&scratch->phi_hi, g_nvox, g_width, g_width, g_width * g_target_ani,
                   g_extx, g_exty, g_extz, scratch->pdb, g_num_atoms, scratch->curr_shift, g_ignored);
      convolve_kernel_inside_fast(&scratch->phi_du, scratch->phi_hi, g_extx, g_exty, g_extz,
                                  g_phi_fx, g_ext_fx, g_norm_hi, g_ignored);
      for (scratch->l = 0; scratch->l < g_nvox; scratch->l++) {
        scratch->corr_hi_lo += (*(scratch->phi_du + scratch->l)) * (*(g_phi_lo + scratch->l));
      }
      break;
    case 2:
      // Modified three step code (uses mask when computing the correlation)
      project_mass(&scratch->phi_hi, g_nvox, g_width, g_width, g_width * g_target_ani,
                   g_extx, g_exty, g_extz, scratch->pdb, g_num_atoms, scratch->curr_shift, g_ignored);
      convolve_kernel_inside_fast(&scratch->phi_du, scratch->phi_hi, g_extx, g_exty, g_extz,
                                  g_phi_fx, g_ext_fx, g_norm_hi, g_ignored);
      /* compute correlation with mask */
      for (scratch->l = 0; scratch->l < g_inside_num; scratch->l++) {
        scratch->corr_hi_lo += (*(scratch->phi_du + g_inside_list[scratch->l].ireal))
                               * (*(g_phi_lo + g_inside_list[scratch->l].ireal));
      }
      break;
    case 3:
      // New one step code (scales with number of atoms in the high res structure)
      project_mass_convolve_kernel_corr(g_width, g_width, g_width * g_target_ani,
                                        g_extx, g_exty, g_extz, scratch->pdb, g_num_atoms, scratch->curr_shift,
                                        g_phi_fx, g_ext_fx, g_norm_hi, g_ignored, g_phi_lo, &scratch->corr_hi_lo);
      break;
    default:
      fprintf(stderr, "colores> Error: Did not understand Powell correlation method, check -pwcorr\n");
      exit(1);
  }
  scratch->corr_hi_lo /= (double)g_nvox;
  return -scratch->corr_hi_lo;
}




/*====================================================================*/
static void get_low_map(char *file_name)
{
  /* reads low resolution map from file_name */

  double cut_width;
  double new_width;

  printf("colores> Processing low-resolution map.\n");

  read_vol(file_name, &g_width, &g_origx, &g_origy, &g_origz, &g_extx, &g_exty, &g_extz, &g_phi_lo);
  g_nvox = g_extx * g_exty * g_extz;

  /* if spacing is too wide adjust resolution */
  if (g_width > g_target_res * 0.7) {
    g_target_res = 2.0 * g_width;
    printf("colores> Warning: Insufficient spatial sampling (voxel spacing too wide) for initially assigned map resolution.\n");
    printf("colores> Target resolution adjusted to 2x voxel spacing, %f.\n", g_target_res);
  }

  /* if spacing is too narrow adjust spacing and prepare for interpolation */
  cut_width = g_target_res * 0.2;
  if (g_width < cut_width) {
    new_width = g_target_res * 0.25;
  } else new_width = g_width;

  /* interpolate if anisotropic map or due to spacing considerations above */
  if (g_target_ani != 1.0 || g_width < cut_width) {
    interpolate_map(&g_phi_du, &g_extx, &g_exty, &g_extz, &g_origx, &g_origy, &g_origz,
                    new_width, new_width, new_width * g_target_ani, g_phi_lo, g_extx, g_exty, g_extz, g_origx,
                    g_origy, g_origz, g_width, g_width, g_width);
    g_nvox = g_extx * g_exty * g_extz;
    cp_vect_destroy(&g_phi_lo, &g_phi_du, g_nvox);
    g_width = new_width;
  }

  /* set density values below g_low_cutoff to zero */
  threshold(g_phi_lo, g_nvox, g_low_cutoff);

  /* shrink map about non-zero density and resize to odd intervals */
  shrink_margin(&g_phi_du, &g_extx, &g_exty, &g_extz, &g_origx, &g_origy, &g_origz, &g_nvox,
                g_phi_lo, g_extx, g_exty, g_extz, g_origx, g_origy, g_origz,
                g_width, g_width, g_width * g_target_ani);
  cp_vect_destroy(&g_phi_lo, &g_phi_du, g_nvox);

  print_map_info(g_phi_lo, g_nvox);
}



/*====================================================================*/
static void get_centered_structure_and_radius(char *file_name, double *max_radius)
{
  /* reads atomic structure and centers it */

  double cen_x;
  double cen_y;
  double cen_z;

  printf("colores> Processing atomic structure.\n");

  /* read PDB file */
  read_pdb(file_name, &g_num_atoms, &g_pdb_original);

  calc_center(g_pdb_original, g_num_atoms, &cen_x, &cen_y, &cen_z);
  *max_radius = calc_sphere(g_pdb_original, g_num_atoms, cen_x, cen_y, cen_z);
  printf("colores> Geometric center: %6.3f %6.3f %6.3f, radius: %6.3f Angstrom\n", cen_x, cen_y, cen_z, *max_radius);

  /* center structure */
  translate(g_pdb_original, g_pdb_original, g_num_atoms, -cen_x, -cen_y, -cen_z);
}



/*====================================================================*/
static void read_options(int argc, char **argv)
{
  /* print usage info and read options from input arguments */

  char *pos;
  int i;
  char option_string[2048];

  if (argc < 2) { /* print usage and options info */

    printf("colores> USAGE:   colores <Density map> <PDB structure>  -<options>\n");
    printf("colores>\n");
    printf("colores> OPTIONS:\n");
    printf("colores>\n");
    printf("colores>\t-res <float>\t   Target resolution in A [default: -res 15]\n");
    printf("colores>\n");
    printf("colores>\t-ani <float>\t   Resolution anisotropy factor [default: -ani 1]\n");
    printf("colores>\n");
    printf("colores>\t-cutoff <float>\t   Density map cutoff value [default:-cutoff 0.0]\n");
    printf("colores>\n");
    printf("colores>\t-corr <int>\t   Correlation method:\n");
    printf("colores>\n");
    printf("colores>\t\t\t\t-corr 0\t-->  Standard cross correlation [default for res < 10A]\n");
    printf("colores>\t\t\t\t-corr 1\t-->  Laplacian filtered correlation [default for res >= 10A] \n");
    printf("colores>\n");
    printf("colores>\t-sizef <float>\t   Grid size expansion factor for FFT zero padding  \n");
    printf("colores>\t\t\t   [default:-sizef 0.1 for standard and 0.2 for Laplacian correlation]\n");
    printf("colores>\n");
    printf("colores>\t-euler <int>\t   Euler angle generation method [default: -euler 0]\n");
    printf("colores>\n");
    printf("colores>\t\t\t\t-euler 0\t\t--> Proportional method\n");
    printf("colores>\t\t\t\t-euler 1\t\t--> Lattman et al. method\n");
    printf("colores>\t\t\t\t-euler 2\t\t--> Spiral method\n");
    printf("colores>\t\t\t\t-euler 3 <filename>\t--> Input file\ncolores>\n");
    printf("colores>\t-erang <float>x6   Euler angle (psi,theta,phi) range in degrees\n");
    printf("colores>\t\t\t   [default: -erang 0 360 0 180 0 360]\n");
    printf("colores>\n");
    printf("colores>\t-deg <float>\t   Angular sampling step in degrees [default: -deg 30.0]\n");
    printf("colores>\n");
    printf("colores>\t-sculptor \t   Save outout files for interactive exploration with Sculptor [default: Off]\n");
    printf("colores>\n");
    printf("colores>\t-nopowell \t   Powell maximization Off [default On]\n");
    printf("colores>\n");
    printf("colores>\t-explor <int>\t   Number of local maxima explored \n");
    printf("colores>\t\t\t    [default: -explor 10]\n");
    printf("colores>\n");
    printf("colores>\t-pwti <float>x2\t   Powell tolerance & max iterations\n");
    printf("colores>\t\t\t    [default: -pwti 1e-6 25]\n");
    printf("colores>\n");
    printf("colores>\t-pwtr <float>x2\t   Trans & Rot initial step size\n");
    printf("colores>\t\t\t    [default: .25 voxel spacing .25 angular sampling]\n");
    printf("colores>\n");
    printf("colores>\t-pwcorr <int>\t   Powell correlation algorithm options\n");
    printf("colores>\t\t\t    [default: -pwcorr 0]\n");
    printf("colores>\n");
    printf("colores>\t\t\t    0: Determined at runtime\n");
    printf("colores>\t\t\t    1: Original three-step code\n");
    printf("colores>\t\t\t    2: Three-step code with mask applied\n");
    printf("colores>\t\t\t    3: One-step code for small probe structures\n");
    printf("colores>\n");
    printf("colores>\t-peak <int>\t   Peak search options [default: -peak 0]\n");
    printf("colores>\n");
    printf("colores>\t\t\t    0: Original peak search by sort and filter \n");
    printf("colores>\t\t\t    1: Peak search by filter only \n");
    printf("colores>\n");
    printf("colores>\t-nopeaksharp\t   Peak sharpness estimation Off [default: On]\n");
    printf("colores>\n");
    printf("colores>\t-nprocs <int>\t   Number of parallel threads [default: usually the number of cores on the CPU]\n");
    printf("colores>\n");
    draw_line();
    exit(0);
  }


  /* now read options from arguments and assign global variables */

  sprintf(option_string, "\n");
  printf("colores> Options read:\n");
  for (i = 1; i < argc; i++) sprintf(option_string, "%s %s", option_string, argv[i]);

  /* target resolution */
  g_target_res = 15.0;
  if ((pos = (char *)strstr(option_string, " -res")) != NULL)
    sscanf(pos + 5, "%lf", &g_target_res);
  if (g_target_res < 0) {
    error_resolution_range(80400, g_program);
  }
  printf("colores> Target resolution %3.3f\n", g_target_res);

  /* resolution anisotropy */
  g_target_ani = 1.0;
  if ((pos = (char *)strstr(option_string, " -ani")) != NULL)
    sscanf(pos + 5, "%lf", &g_target_ani);
  if (g_target_ani < 0.001 || g_target_ani > 1000) {
    error_anisotropy_range(80410, g_program);
  }
  printf("colores> Resolution anisotropy %3.3f\n", g_target_ani);


  /* non-negative cutoff of the low-resolution map */
  g_low_cutoff = 0;
  if ((pos = (char *)strstr(option_string, " -cutoff")) != NULL)
    sscanf(pos + 8, "%lf", &g_low_cutoff);
  if (g_low_cutoff < 0) {
    g_low_cutoff = 0;
    printf("colores> Low-resolution map cutoff must be non-negative, assigning %3.3f\n", g_low_cutoff);
  } else printf("colores> Low-resolution map cutoff %3.3f\n", g_low_cutoff);

  /* correlation mode */
  if (g_target_res >= 10) g_corr_mode = 1;
  else g_corr_mode = 0;
  if ((pos = (char *)strstr(option_string, " -corr")) != NULL)
    sscanf(pos + 6, "%d", &g_corr_mode);
  switch (g_corr_mode) {
    case 0:
      printf("colores> Standard cross correlation\n");
      break;
    case 1:
      printf("colores> Laplacian filtered correlation\n");
      break;
    default:
      error_option(80420, g_program);
  }

  /* grid size expansion for FFT zero padding */
  /* 0.25 is safe upper bound given in Numerical Recipes, but empirically found 0.1 / 0.2 are sufficient for std and Laplacian correlation */
  if (g_corr_mode == 0) g_size_fac = 0.1;
  else g_size_fac = 0.2;
  if ((pos = (char *)strstr(option_string, " -sizef")) != NULL)
    sscanf(pos + 7, "%lf", &g_size_fac);
  if (g_size_fac < 0) {
    printf("colores> Warning: grid size expansion factor must be positive, will be set to zero\n");
    g_size_fac = 0.0;
  }
  if (g_size_fac > 1) {
    printf("colores> Warning: grid size expansion seems too large, will be set to one.\n");
    g_size_fac = 1.0;
  }
  printf("colores> FFT grid size expansion factor %1.3f (thickness of additional zero layer as fraction of map dimensions)\n", g_size_fac);
  if (g_size_fac > 0.25) printf("colores> Note: FFT requires no more than -sizef 0.25 for numerical stability.\n");

  /* Euler angle generation */
  g_euler_mode = 0;
  if ((pos = (char *)strstr(option_string, " -euler")) != NULL)
    sscanf(pos + 7, "%d", &g_euler_mode);

  switch (g_euler_mode) {
    case 0:
      printf("colores> Euler angles generation using Proportional method\n");
      break;
    case 1:
      printf("colores> Euler angles generation using Lattman et al. method\n");
      break;
    case 2:
      printf("colores> Euler angles generation using Spiral method\n");
      break;
    case 3: {
        sscanf((char *)(strstr(option_string, "-euler") + 9), "%s", g_eu_in_file);
        printf("colores> Euler angles are read from file %s\n", g_eu_in_file);
        break;
      }
    default:
      error_option(80430, g_program);
  }

  /* angular sampling */
  g_delta_rot = 30;
  if ((pos = (char *)strstr(option_string, " -deg")) != NULL) sscanf(pos + 5, "%lf", &g_delta_rot);
  printf("colores> Angular sampling accuracy %3.3f\n", g_delta_rot);
  if (g_delta_rot <= 1) {
    error_euler_sampling(80440, g_program);
  }

  /* Euler angle range */
  g_eu_range[0][0] = 0.;
  g_eu_range[1][0] = 0.;
  g_eu_range[2][0] = 0.;
  g_eu_range[0][1] = 360.;
  g_eu_range[1][1] = 180.;
  g_eu_range[2][1] = 360.;
  if ((pos = (char *)strstr(option_string, " -erang")) != NULL)
    sscanf(pos + 7, "%lf %lf %lf %lf %lf %lf", &g_eu_range[0][0], &(g_eu_range[0][1]), &g_eu_range[1][0], &g_eu_range[1][1], &g_eu_range[2][0], &g_eu_range[2][1]);

  if (g_eu_range[0][0] > g_eu_range[0][1] || g_eu_range[1][0] > g_eu_range[1][1] || g_eu_range[2][0] > g_eu_range[2][1]) {
    error_euler_below_start(80450, g_program);
    //fprintf(stderr, "colores> Error: Euler angle range end value below start value [e.c. 80450]\n");
    //exit(80450);
  }
  if (g_eu_range[0][0] < -360. || g_eu_range[1][0] < -360. || g_eu_range[2][0] < -360.) {
    error_euler_below_neg_360(80460, g_program);
    //fprintf(stderr, "colores> Error: Euler angle range start value below -360 [e.c. 80460]\n");
    //exit(80460);
  }
  if (g_eu_range[0][0] > 360. || g_eu_range[1][0] > 360. || g_eu_range[2][0] > 360.) {
    error_euler_above_pos_360(80470, g_program);
    //fprintf(stderr, "colores> Error: Euler angle range start value above +360 [e.c. 80470]\n");
    //exit(80470);
  }
  if (g_eu_range[0][1] - g_eu_range[0][0] > 360.) {
    error_psi_euler_range_above_360(80480, g_program);
    //fprintf(stderr, "colores> Error: First (psi) Euler range exceeds 360 [e.c. 80480]\n");
    //exit(80480);
  }
  if (g_eu_range[1][1] - g_eu_range[1][0] > 180.) {
    error_theta_euler_range_above_180(80490, g_program);
    //fprintf(stderr, "colores> Error: Second (theta) Euler range exceeds 180 [e.c. 80490]\n");
    //exit(80490);
  }
  if (g_eu_range[2][1] - g_eu_range[2][0] > 360.) {
    error_phi_euler_range_above_360(80500, g_program);
    //fprintf(stderr, "colores> Error: Third (phi) Euler range exceeds 360 [e.c. 80500]\n");
    //exit(80500);
  }
  printf("colores> Euler angle range: ");

  for (i = 0; i < 3; i++)
    printf("[%5.3f:%5.3f] ", g_eu_range[i][0], g_eu_range[i][1]);
  printf("\n");

  /* Sculptor interactive mode */
  g_sculptor = 0;
  if ((pos = (char *)strstr(option_string, " -sculptor")) != NULL) {
    if (g_target_ani < 0.99999 || g_target_ani > 1.00001) printf("colores> Warning: Unable to turn on Sculptor mode in case of anisotropic maps\n");
    else g_sculptor = 1;
  }
  if (g_sculptor) printf("colores> Sculptor mode ON\n");
  else printf("colores> Sculptor mode OFF\n");

  /* number of best fits */
  g_num_explored = 10;
  if ((pos = (char *)strstr(option_string, " -explor")) != NULL)
    sscanf(pos + 8, "%d", &g_num_explored);
  printf("colores> Number of best fits explored %2d\n", g_num_explored);

  /* peak search method */
  g_peak_opt = 0;
  if ((pos = (char *)strstr(option_string, " -peak")) != NULL)
    sscanf(pos + 6, "%d", &g_peak_opt);

  switch (g_peak_opt) {
    case 0:
      printf("colores> Original peak search by sort and filter\n");
      break;
    case 1:
      printf("colores> New peak search by filter only\n");
      break;
    default:
      error_option(80420, g_program);
  }

  /* Powell maximization toggle */
  g_pow_mode = 1;
  if ((pos = (char *)strstr(option_string, " -nopowell")) != NULL) {
    g_pow_mode = 0;
    printf("colores> Powell maximization OFF\n");
  } else  printf("colores> Powell maximization ON\n");

  /* Powell parameters */
  g_pow_tolerance = 1.0e-6;
  g_pow_max_iter = 25;
  if ((pos = (char *)strstr(option_string, " -pwti")) != NULL)
    sscanf(pos + 6, "%lf %d", &g_pow_tolerance, &g_pow_max_iter);
  printf("colores> Powell tolerance %1.2E  Max iterations %d\n", g_pow_tolerance, g_pow_max_iter);
  g_pow_delta_pos = -9;
  g_pow_delta_ang = -9;
  if ((pos = (char *)strstr(option_string, " -pwtr")) != NULL)
    sscanf(pos + 6, "%lf %lf", &g_pow_delta_pos, &g_pow_delta_ang);
  if (g_pow_delta_pos < 0 || g_pow_delta_ang < 0)
    printf("colores> Powell trans & rot initial step sizes set to default values\n");
  else printf("colores> Powell trans step size %1.3f, rot step size %4.3f \n", g_pow_delta_pos, g_pow_delta_ang);

  /* Powell correlation algorithm */
  g_pow_alg = 0;
  if ((pos = (char *)strstr(option_string, " -pwcorr")) != NULL) sscanf(pos + 8, "%d", &g_pow_alg);
  if (g_pow_alg == 0) {
    printf("colores> Powell correlation algorithm determined automatically\n");
  } else {
    printf("colores> Powell correlation algorithm: %d\n", g_pow_alg);
  }
  if (g_pow_alg > 3) {
    error("colores> Error: Invalid Powell correlation algorithm\n");
  }

  /* Sensitivity toggle */
  g_peak_sharp = 1;
  if ((pos = (char *)strstr(option_string, " -nopeaksharp")) != NULL) {
    g_peak_sharp = 0;
    printf("colores> Peak sharpness estimation OFF\n");
  } else printf("colores> Peak sharpness estimation ON\n");

  /* Number of threads */
  g_p_nprocs = 1;
#ifdef _SMP_
  if ((pos = (char *)strstr(option_string, " -nprocs")) != NULL) {
    sscanf(pos + 8, "%d", &g_p_nprocs);
    printf("colores> Number of parallel threads requested by user: %d\n", g_p_nprocs);
  } else {
#ifdef _OPENMP
    int nthreads_omp, tid_omp;
    #pragma omp parallel shared(nthreads_omp) private(tid_omp)
    {
      tid_omp = omp_get_thread_num();
      if (tid_omp == 0) nthreads_omp = omp_get_num_threads();
    }
    g_p_nprocs = nthreads_omp;
	if (g_p_nprocs > 16) g_p_nprocs = 16;
#else
    printf("colores> Warning: Unable to determine number of cores on CPU. Use -nprocs to set the number of threads. \n");
#endif
    printf("colores> Number of parallel threads automatically assigned by program: %d\n", g_p_nprocs);
  }
  if (g_p_nprocs <= 0) {
    error("colores> Error: Number of threads not valid.\n");
  }
#else
  fprintf(stderr, "colores> Warning: Parallel (SMP) thread support not compiled in, running in serial mode.\n");
#endif

}



/*====================================================================*/
static void peak_extract(SAV **ccr, FIT **found_peak_final)
{
  /* Translational space peak filter routine.                   */
  /* Called only once.                                          */
  /* Note that ccr has translations inverted due to definition  */
  /* of Fourier correlation. So we invert them back to the      */
  /* actual displacements.                                      */

  double smooth_filter[3][3][3] = {
    {
      {0.02777777777778, 0.02777777777778, 0.02777777777778},
      {0.02777777777778, 0.05555555555556, 0.02777777777778},
      {0.02777777777778, 0.02777777777778, 0.02777777777778}
    },
    {
      {0.02777777777778, 0.05555555555556, 0.02777777777778},
      {0.05555555555556, 0.11111111111111, 0.05555555555556},
      {0.02777777777778, 0.05555555555556, 0.02777777777778}
    },
    {
      {0.02777777777778, 0.02777777777778, 0.02777777777778},
      {0.02777777777778, 0.05555555555556, 0.02777777777778},
      {0.02777777777778, 0.02777777777778, 0.02777777777778}
    }
  };

  double peak_filter[3][3][3] = {
    {{ -1, -1, -1}, { -1, -1, -1}, { -1, -1, -1}},
    {{ -1, -1, -1}, { -1, 26, -1}, { -1, -1, -1}},
    {{ -1, -1, -1}, { -1, -1, -1}, { -1, -1, -1}}
  };

  double *map1, *map2;
  float curr_cc, filter_sig, filter_max, curr_msd;
  float curr_diff, curr_max, search_cut;
  unsigned n, m, p, v;
  unsigned peak_count, peak_count_all;
  FIT *found_peak;
  FIT *max_cc_peak;
  int r, s, t;
  unsigned long q, map_index, euler_index;
  float pos_diff;

  /* save the g_num_explored highest scoring voxels so we don't filter them out by accident */
  /* allocate */
  max_cc_peak = (FIT *) alloc_vect(g_num_explored, sizeof(FIT));

  /* initialize */
  for (peak_count = 0; peak_count < g_num_explored; peak_count++) {
    max_cc_peak[peak_count].score = -99999.0;
    for (m = 0; m < 3; m++) max_cc_peak[peak_count].euler[m] = 0.0;
    for (m = 0; m < 3; m++) max_cc_peak[peak_count].pos[m] = 0.0;
  }
  /* search */
  for (v = 0; v < g_inside_num_flipped; v++) {
    m = g_inside_list_flipped[v].iz;
    n = g_inside_list_flipped[v].iy;
    p = g_inside_list_flipped[v].ix;
    q = p + (g_extx) * (n + g_exty * m);
    curr_cc = (*(*ccr + q)).score;
    if (curr_cc > max_cc_peak[0].score) {
      euler_index = (*(*ccr + q)).eu;
      max_cc_peak[0].score = curr_cc;
      max_cc_peak[0].euler[0] = *(g_eulers + 3 * euler_index + 0);
      max_cc_peak[0].euler[1] = *(g_eulers + 3 * euler_index + 1);
      max_cc_peak[0].euler[2] = *(g_eulers + 3 * euler_index + 2);
      max_cc_peak[0].pos[0] = g_width * g_extx_half - g_width * p;
      max_cc_peak[0].pos[1] = g_width * g_exty_half - g_width * n;
      max_cc_peak[0].pos[2] = g_width * g_extz_half - g_width * m;
      /* now bring lowest scoring peak to front */
      for (peak_count = 1; peak_count < g_num_explored; peak_count++)
        if (max_cc_peak[0].score > max_cc_peak[peak_count].score)
          SWAPPING(max_cc_peak[peak_count], max_cc_peak[0], FIT);
    }
  }

  /* now start filtering */

  do_vect(&map1, (g_extz + 2) * (g_exty + 2) * (g_extx + 2));

  /* Smooth the saved CC(r) with Gaussian-style kernel that   */
  /* takes into account only the adjacent voxels.             */
  /* Smoothing reduces noise and thereby the number of peaks. */
  /* The smoothed map, map1, will be shifted by one voxel.    */

  for (v = 0; v < g_inside_num_flipped; v++) {
    m = g_inside_list_flipped[v].iz;
    n = g_inside_list_flipped[v].iy;
    p = g_inside_list_flipped[v].ix;
    q = p + (g_extx) * (n + g_exty * m);
    curr_cc = (*(*ccr + q)).score;
    if (curr_cc != 0.0f)
      for (r = -1; r < 2; r++)
        for (s = -1; s < 2; s++)
          for (t = -1; t < 2; t++) {
            map_index = (t + p + 1) + (g_extx + 2) * ((s + n + 1) + (g_exty + 2) * (r + m + 1));
            *(map1 + map_index) += smooth_filter[t + 1][s + 1][r + 1] * curr_cc;
          }
  }


  do_vect(&map2, (g_extz + 2) * (g_exty + 2) * (g_extx + 2));

  /* Apply Laplacian-style peak detection filter to map1;     */
  /* map2 will also be shifted by one voxel relative to CC(r) */

  for (v = 0; v < g_inside_num_flipped; v++) {
    m = g_inside_list_flipped[v].iz;
    n = g_inside_list_flipped[v].iy;
    p = g_inside_list_flipped[v].ix;
    q = (p + 1) + (g_extx + 2) * ((n + 1) + (g_exty + 2) * (m + 1));
    for (r = -1; r < 2; r++)
      for (s = -1; s < 2; s++)
        for (t = -1; t < 2; t++) {
          map_index = (t + p + 1) + (g_extx + 2) * ((s + n + 1) + (g_exty + 2) * (r + m + 1));
          *(map2 + q) += peak_filter[t + 1][s + 1][r + 1] * *(map1 + map_index);
        }
  }
  free_vect_and_zero_ptr(&map1);

  /* In the following we define as filter contrast the density difference     */
  /* between a voxel in the Laplacian filtered map2 and its 26 neighbors,     */
  /* averaged over the neighbors.                                             */
  /* First we compute the maximum and sigma of the filter contrast.           */
  /* Note: The average filter contrast is zero.                               */

  curr_msd = 0.0f;
  filter_max = 0.0f;

  for (v = 0; v < g_inside_num_flipped; v++) {
    m = g_inside_list_flipped[v].iz;
    n = g_inside_list_flipped[v].iy;
    p = g_inside_list_flipped[v].ix;
    q = (p + 1) + (g_extx + 2) * ((n + 1) + (g_exty + 2) * (m + 1));

    curr_max = 0.0f;
    for (r = -1; r < 2; r++)
      for (s = -1; s < 2; s++)
        for (t = -1; t < 2; t++) {
          map_index = (t + p + 1) + (g_extx + 2) * ((s + n + 1) + (g_exty + 2) * (r + m + 1));
          curr_max += *(map2 + q) - *(map2 + map_index);
        }
    curr_max /= 26.0f;
    if (curr_max > filter_max) filter_max = curr_max;
    curr_msd += curr_max * curr_max;
  }

  filter_sig = sqrt(curr_msd / ((float)g_inside_num_flipped));
  printf("colores> Peak filter contrast: maximum %f, sigma %f\n", filter_max, filter_sig);

  /* Now we have an idea of the filter contrast distribution. We limit our peak search */
  /* to above the noise level (2*filter_sig capped at 0.25 max) and eliminate redundant peaks. */

  search_cut = 2.0f * filter_sig;
  if (search_cut > filter_max * 0.25f) search_cut = filter_max * 0.25f;


  /* Count number of possibly redundant peaks that meet the criterion */
  peak_count = 0;

  for (v = 0; v < g_inside_num_flipped; v++) {
    m = g_inside_list_flipped[v].iz;
    n = g_inside_list_flipped[v].iy;
    p = g_inside_list_flipped[v].ix;
    q = (p + 1) + (g_extx + 2) * ((n + 1) + (g_exty + 2) * (m + 1));
    curr_diff = 0.0f;
    for (r = -1; r < 2; r++)
      for (s = -1; s < 2; s++)
        for (t = -1; t < 2; t++) {
          map_index = (t + p + 1) + (g_extx + 2) * ((s + n + 1) + (g_exty + 2) * (r + m + 1));
          curr_diff += *(map2 + q) - *(map2 + map_index);
        }
    curr_diff /= 26.0f;
    if (curr_diff > search_cut) peak_count++;
  }


  printf("colores> Contrast threshold: %f, candidate peaks: %d\n", search_cut, peak_count + g_num_explored);

  /* Allocate storage for found peaks */

  found_peak = (FIT *) alloc_vect((peak_count + g_num_explored), sizeof(FIT));

  /* Now extract the peaks that were counted before */

  peak_count = 0;

  for (v = 0; v < g_inside_num_flipped; v++) {
    m = g_inside_list_flipped[v].iz;
    n = g_inside_list_flipped[v].iy;
    p = g_inside_list_flipped[v].ix;
    q = (p + 1) + (g_extx + 2) * ((n + 1) + (g_exty + 2) * (m + 1));
    curr_diff = 0.0f;
    for (r = -1; r < 2; r++)
      for (s = -1; s < 2; s++)
        for (t = -1; t < 2; t++) {
          map_index = (t + p + 1) + (g_extx + 2) * ((s + n + 1) + (g_exty + 2) * (r + m + 1));
          curr_diff += *(map2 + q) - *(map2 + map_index);
        }
    curr_diff /= 26.0f;
    if (curr_diff > search_cut) {
      q = p + (g_extx) * (n + g_exty * (m));
      euler_index = (*(*ccr + q)).eu;
      found_peak[peak_count].score = (*(*ccr + q)).score;
      found_peak[peak_count].euler[0] = *(g_eulers + 3 * euler_index + 0);
      found_peak[peak_count].euler[1] = *(g_eulers + 3 * euler_index + 1);
      found_peak[peak_count].euler[2] = *(g_eulers + 3 * euler_index + 2);
      found_peak[peak_count].pos[0] = g_width * g_extx_half - g_width * p;
      found_peak[peak_count].pos[1] = g_width * g_exty_half - g_width * n;
      found_peak[peak_count].pos[2] = g_width * g_extz_half - g_width * m;
      peak_count++;
    }
  }

  free_vect_and_zero_ptr(&map2);

  /* If needed add the g_num_explored saved maximum scoring peaks */
  if (g_peak_opt == 0) {
    for (n = 0; n < g_num_explored; n++) {
      found_peak[peak_count + n].score = max_cc_peak[n].score;
      for (m = 0; m < 3; m++) found_peak[peak_count + n].euler[m] = max_cc_peak[n].euler[m];
      for (m = 0; m < 3; m++) found_peak[peak_count + n].pos[m] = max_cc_peak[n].pos[m];
    }
    peak_count += g_num_explored;
  }
  free_vect_and_zero_ptr(&max_cc_peak);

  /* Sort the extracted peaks as a function of correlation coefficient. */

  for (p = 0; p < peak_count; p++)
    for (m = p + 1; m < peak_count; m++)
      if (found_peak[m].score > found_peak[p].score)
        SWAPPING(found_peak[p], found_peak[m], FIT);

  peak_count_all = peak_count;

  /* Eliminate redundant peaks.                                 */
  /* So far we haven't looked at the peak distribution in       */
  /* translational space. Adjacent voxels may exhibit different */
  /* saved rotations, even though the rotational space has not  */
  /* been saved. In this sense rotational alternatives are      */
  /* still available. Therefore, we remove only those peaks     */
  /* that have both similar positions AND similar angles.       */
  /* This preserves any rotational alternatives that might      */
  /* still exist between nearby voxels.                         */

  for (p = 0; p < peak_count_all; p++)
    if (found_peak[p].score != -99999.0)
      for (r = p + 1; r < peak_count_all; r++)
        if (found_peak[r].score != -99999.0) {
          pos_diff = sqrt(((found_peak[p].pos[0] - found_peak[r].pos[0]) * (found_peak[p].pos[0] - found_peak[r].pos[0]) +
                           (found_peak[p].pos[1] - found_peak[r].pos[1]) * (found_peak[p].pos[1] - found_peak[r].pos[1]) +
                           (found_peak[p].pos[2] - found_peak[r].pos[2]) * (found_peak[p].pos[2] - found_peak[r].pos[2])));
          /* sparsification by spatial resolution */
          if ((pos_diff < g_target_res) &&
              (similar_eulers(found_peak[p].euler[0], found_peak[p].euler[1],
                              found_peak[p].euler[2], found_peak[r].euler[0],
                              found_peak[r].euler[1], found_peak[r].euler[2]))
             ) found_peak[r].score = -99999.0;
        }

  /* remove marked redundant peaks */
  peak_count = 0;
  for (p = 0; p < peak_count_all; p++)
    if (found_peak[p].score != -99999.0) {
      found_peak[peak_count].score = found_peak[p].score;
      found_peak[peak_count].euler[0] = found_peak[p].euler[0];
      found_peak[peak_count].euler[1] = found_peak[p].euler[1];
      found_peak[peak_count].euler[2] = found_peak[p].euler[2];
      found_peak[peak_count].pos[2] = found_peak[p].pos[2];
      found_peak[peak_count].pos[1] = found_peak[p].pos[1];
      found_peak[peak_count].pos[0] = found_peak[p].pos[0];
      peak_count++;
    }

  printf("colores> Found %d non-redundant peaks.\n", peak_count);

  /* Adjust g_num_explored if necessary. */
  if (peak_count < g_num_explored) {
    printf("colores> Warning: Peak count below -explor %d option.\n", g_num_explored);
    printf("colores> Resetting number of explored peaks to %d. \n", peak_count);
    g_num_explored = peak_count;
  }

  *(found_peak_final) = (FIT *) alloc_vect(g_num_explored, sizeof(FIT));

  for (r = 0; r < g_num_explored; r++) {
    *(*found_peak_final + r) = found_peak[r];
  }
  free_vect_and_zero_ptr(&found_peak);
}


/*====================================================================*/
static void flip_quadrants(unsigned rix, unsigned riy, unsigned riz,
                           unsigned *rox, unsigned *roy, unsigned *roz)
{
  /* converts indexing scheme of FFT lattice to that of real space lattice */
  /* quadrants are "flipped" in FFT scheme                                 */
  /* called once by create_inside_molecule_poslist                         */
  /* and once by create_inside_molecule_poslist_flipped                    */
  /* input: rix, riy, riz                                                  */
  /* output: rox, roy, roz                                                 */

  int pos[8];
  unsigned quad_boundary_x, quad_boundary_y, quad_boundary_z;

  quad_boundary_x = g_extx_half + 1;
  quad_boundary_y = g_exty_half + 1;
  quad_boundary_z = g_extz_half + 1;

  pos[0] = g_extx_half + rix;
  pos[1] = rix - quad_boundary_x;
  pos[2] = g_exty_half + riy;
  pos[3] = riy - quad_boundary_y;
  pos[4] = g_extz_half + riz;
  pos[5] = riz - quad_boundary_z;

  if (rix < quad_boundary_x) {
    if (riy < quad_boundary_y) {
      if (riz < quad_boundary_z) {
        *rox = pos[0];
        *roy = pos[2];
        *roz = pos[4];
      } else {
        *rox = pos[0];
        *roy = pos[2];
        *roz = pos[5];
      }
    } else {
      if (riz < quad_boundary_z) {
        *rox = pos[0];
        *roy = pos[3];
        *roz = pos[4];
      } else {
        *rox = pos[0];
        *roy = pos[3];
        *roz = pos[5];
      }
    }
  } else {
    if (riy < quad_boundary_y) {
      if (riz < quad_boundary_z) {
        *rox = pos[1];
        *roy = pos[2];
        *roz = pos[4];
      } else {
        *rox = pos[1];
        *roy = pos[2];
        *roz = pos[5];
      }
    } else {
      if (riz < quad_boundary_z) {
        *rox = pos[1];
        *roy = pos[3];
        *roz = pos[4];
      } else {
        *rox = pos[1];
        *roy = pos[3];
        *roz = pos[5];
      }
    }
  }
}


/*====================================================================*/
static void create_inside_molecule_poslist(double *phi)
{
  /* calculates all positions in map that have positive density or are buried */
  /* called once in preprocessing stage of colores                            */
  /* writes list to g_inside_list number to g_inside_num                      */

  unsigned curr, i, p, m, n;
  unsigned long q, s;
  unsigned check, ix, iy, iz;
  int ix2, iy2, iz2, erosion_shell_width;
  char *mask_inside, *mask_inside2;

  mask_inside = (char *) alloc_vect(g_nvox, sizeof(char));

  g_inside_num = 0;
  printf("colores> Identifying inside or buried voxels...\n");
  /* mark and count number of voxels */
  for (iz = 0; iz < g_extz; iz++)
    for (iy = 0; iy < g_exty; iy++)
      for (ix = 0; ix < g_extx; ix++) {
        q = ix + (g_extx) * (iy + g_exty * iz); /* actual lattice index */
        if (phi[q] > 0.000) { /* inside or on contour */
          mask_inside[q] = 1;
          g_inside_num++;
        } else { /* scan in all directions to check if inside */
          check = 0;
          for (i = 0; i < ix; i++)
            if (phi[i + (g_extx) * (iy + g_exty * iz)] > 0.0000) {
              check++;
              i = ix;
            }
          for (i = ix + 1; i < g_extx; i++)
            if (phi[i + (g_extx) * (iy + g_exty * iz)] > 0.0000) {
              check++;
              i = g_extx;
            }
          for (i = 0; i < iy; i++)
            if (phi[ix + (g_extx) * (i + g_exty * iz)] > 0.0000) {
              check++;
              i = iy;
            }
          for (i = iy; i < g_exty; i++)
            if (phi[ix + (g_extx) * (i + g_exty * iz)] > 0.0000) {
              check++;
              i = g_exty;
            }
          for (i = 0; i < iz; i++)
            if (phi[ix + (g_extx) * (iy + g_exty * i)] > 0.0000) {
              check++;
              i = iz;
            }
          for (i = iz; i < g_extz; i++)
            if (phi[ix + (g_extx) * (iy + g_exty * i)] > 0.0000) {
              check++;
              i = g_extz;
            }
          if (check >= 4) { /* probably inside */
            mask_inside[q] = 1;
            g_inside_num++;
          }
        }
      }

  /* erode surface of inside subset, we want central part only */

  mask_inside2 = (char *) alloc_vect(g_nvox, sizeof(char));
  for (q = 0; q < g_nvox; q++)
    mask_inside2[q] = mask_inside[q];

  erosion_shell_width = 1;
  for (iz = 0; iz < g_extz; iz++)
    for (iy = 0; iy < g_exty; iy++)
      for (ix = 0; ix < g_extx; ix++) {
        q = ix + (g_extx) * (iy + g_exty * iz); /* actual lattice index */
        if (mask_inside[q] == 1) {
          check = 0;
          for (iz2 = -erosion_shell_width; check == 0 && iz2 <= erosion_shell_width; iz2++)
            for (iy2 = -erosion_shell_width; check == 0 && iy2 <= erosion_shell_width; iy2++)
              for (ix2 = -erosion_shell_width; check == 0 && ix2 <= erosion_shell_width; ix2++) {
                if (ix + ix2 < g_extx && iy + iy2 < g_exty && iz + iz2 < g_extz &&
                    ix + ix2 >= 0 && iy + iy2 >= 0 && iz + iz2 >= 0) {
                  s = (ix + ix2) + (g_extx) * ((iy + iy2) + g_exty * (iz + iz2));
                  if (*(mask_inside + s) == 0)  check = 1; /* surface point */
                } else check = 2; /* out of bounds */
              }
          if ((check > 0) && (mask_inside2[q] == 1)) {
            mask_inside2[q] = 0;
            g_inside_num--;
          }
        }
      }
  free_vect_and_zero_ptr(&mask_inside);

  printf("colores> Found %d inside or buried voxels (out of a total of %lu).\n", g_inside_num, g_nvox);

  /* allocate inside poslist */
  g_inside_list = (POS *) alloc_vect(g_inside_num, sizeof(POS));

  /* now fill inside poslist */
  curr = 0;
  for (m = 0; m < g_extz; m++)
    for (n = 0; n < g_exty; n++)
      for (p = 0; p < g_extx; p++) {
        s = p + (g_extx) * (n + g_exty * m); /* flipped FFT lattice index */
        flip_quadrants(p, n, m, &ix, &iy, &iz);
        q = ix + (g_extx) * (iy + g_exty * iz); /* actual lattice index */
        if (mask_inside2[q] == 1) { /* inside */
          g_inside_list[curr].ifft = s;
          g_inside_list[curr].ireal = q;
          g_inside_list[curr].ix = ix;
          g_inside_list[curr].iy = iy;
          g_inside_list[curr].iz = iz;
          curr++;
        }
      }
  free_vect_and_zero_ptr(&mask_inside2);
}


/*====================================================================*/
static void draw_line()
{
  printf("_____________________________________________________________________________\n");
}


/*====================================================================*/
static void create_inside_molecule_poslist_flipped(double *phi)
{
  /* calculates all positions in map that have positive density or are buried */
  /* saves the positions flipped                */
  /* called once in preprocessing stage of colores                            */
  /* writes list to g_inside_list_flipped, number to g_inside_num_flipped     */

  unsigned curr, i, p, m, n;
  unsigned long q, s;
  unsigned check, ix, iy, iz;
  int ix2, iy2, iz2, erosion_shell_width;
  char *mask_inside, *mask_inside2, *mask_inside3;

  mask_inside = (char *) alloc_vect(g_nvox, sizeof(char));

  g_inside_num_flipped = 0;
  printf("colores> Identifying inside or buried voxels and creating flipped mask...\n");
  /* mark and count number of voxels */
  for (iz = 0; iz < g_extz; iz++)
    for (iy = 0; iy < g_exty; iy++)
      for (ix = 0; ix < g_extx; ix++) {
        q = ix + (g_extx) * (iy + g_exty * iz); /* actual lattice index */
        if (phi[q] > 0.000) { /* inside or on contour */
          mask_inside[q] = 1;
          g_inside_num_flipped++;
        } else { /* scan in all directions to check if inside */
          check = 0;
          for (i = 0; i < ix; i++)
            if (phi[i + (g_extx) * (iy + g_exty * iz)] > 0.0000) {
              check++;
              i = ix;
            }
          for (i = ix + 1; i < g_extx; i++)
            if (phi[i + (g_extx) * (iy + g_exty * iz)] > 0.0000) {
              check++;
              i = g_extx;
            }
          for (i = 0; i < iy; i++)
            if (phi[ix + (g_extx) * (i + g_exty * iz)] > 0.0000) {
              check++;
              i = iy;
            }
          for (i = iy; i < g_exty; i++)
            if (phi[ix + (g_extx) * (i + g_exty * iz)] > 0.0000) {
              check++;
              i = g_exty;
            }
          for (i = 0; i < iz; i++)
            if (phi[ix + (g_extx) * (iy + g_exty * i)] > 0.0000) {
              check++;
              i = iz;
            }
          for (i = iz; i < g_extz; i++)
            if (phi[ix + (g_extx) * (iy + g_exty * i)] > 0.0000) {
              check++;
              i = g_extz;
            }
          if (check >= 4) { /* probably inside */
            mask_inside[q] = 1;
            g_inside_num_flipped++;
          }
        }
      }

  /* erode surface of inside subset, we want central part only */

  mask_inside2 = (char *) alloc_vect(g_nvox, sizeof(char));
  for (q = 0; q < g_nvox; q++)
    mask_inside2[q] = mask_inside[q];

  erosion_shell_width = 1;
  for (iz = 0; iz < g_extz; iz++)
    for (iy = 0; iy < g_exty; iy++)
      for (ix = 0; ix < g_extx; ix++) {
        q = ix + (g_extx) * (iy + g_exty * iz); /* actual lattice index */
        if (mask_inside[q] == 1) {
          check = 0;
          for (iz2 = -erosion_shell_width; check == 0 && iz2 <= erosion_shell_width; iz2++)
            for (iy2 = -erosion_shell_width; check == 0 && iy2 <= erosion_shell_width; iy2++)
              for (ix2 = -erosion_shell_width; check == 0 && ix2 <= erosion_shell_width; ix2++) {
                if (ix + ix2 < g_extx && iy + iy2 < g_exty && iz + iz2 < g_extz &&
                    ix + ix2 >= 0 && iy + iy2 >= 0 && iz + iz2 >= 0) {
                  s = (ix + ix2) + (g_extx) * ((iy + iy2) + g_exty * (iz + iz2));
                  if (*(mask_inside + s) == 0)  check = 1; /* surface point */
                } else check = 2; /* out of bounds */
              }
          if ((check > 0) && (mask_inside2[q] == 1)) {
            mask_inside2[q] = 0;
            g_inside_num_flipped--;
          }
        }
      }
  free_vect_and_zero_ptr(&mask_inside);

  // ADDED LINES BEGIN (PLUS THE MASK_INSIDE3 DECLARATION ABOVE)

  /* mask is flipped along all three axes */
  mask_inside3 = (char *) alloc_vect(g_nvox, sizeof(char));
  for (iz = 0; iz < g_extz; iz++)
    for (iy = 0; iy < g_exty; iy++)
      for (ix = 0; ix < g_extx; ix++) {
        q = ix + (g_extx) * (iy + g_exty * iz);
        s = (g_extx - 1 - ix) + (g_extx) * ((g_exty - 1 - iy) + g_exty * (g_extz - 1 - iz));
        mask_inside3[q] = mask_inside2[s];
      }
  free_vect_and_zero_ptr(&mask_inside2);

  // ADDED LINES FINISH

  printf("colores> Found %d inside or buried voxels (out of a total of %lu).\n", g_inside_num_flipped, g_nvox);

  /* allocate inside poslist */
  g_inside_list_flipped = (POS *) alloc_vect(g_inside_num_flipped, sizeof(POS));

  /* now fill inside poslist */
  curr = 0;
  for (m = 0; m < g_extz; m++)
    for (n = 0; n < g_exty; n++)
      for (p = 0; p < g_extx; p++) {
        s = p + (g_extx) * (n + g_exty * m); /* flipped FFT lattice index */
        flip_quadrants(p, n, m, &ix, &iy, &iz);
        q = ix + (g_extx) * (iy + g_exty * iz); /* actual lattice index */
        if (mask_inside3[q] == 1) { /* inside */ // THIS LINE HAS BEEN MODIFIED
          g_inside_list_flipped[curr].ifft = s;
          g_inside_list_flipped[curr].ireal = q;
          g_inside_list_flipped[curr].ix = ix;
          g_inside_list_flipped[curr].iy = iy;
          g_inside_list_flipped[curr].iz = iz;
          curr++;
        }
      }
  free_vect_and_zero_ptr(&mask_inside3); // THIS LINE HAS BEEN MODIFIED
}

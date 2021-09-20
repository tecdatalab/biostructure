/*********************************************************************
*                         C O L L A G E                              *
**********************************************************************
* Program is part of the Situs package URL: situs.biomachina.org     *
* (c) Willy Wriggers, Jochen Heyd, Valerio Mariani, 2005-2011        *
**********************************************************************
*                                                                    *
* Calculates cross correlation and performs single-pass of Powell    *
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



/* local functions, some are similar or identical to colores but largely use global variables so not placed in a library */
static void multi_rot_euler_trans(unsigned, unsigned, unsigned, PDB *, 
                                  PDB *, double, double, double, double, 
                                  double, double);
static double powell_correlation_multi(double *, POWSCR *);
static void get_structures_and_partition(char **);
static void get_low_map(char *);
static void read_options(int, char **);
static void draw_line();
static void powell_optimization_multi(void *, char **);
static void create_inside_molecule_poslist(double *);
static double max_projection_density(double);

/* global variables */
typedef struct {
  double score;
  double pos[3];
  double euler[3];
} FIT;
typedef struct {
  unsigned long ifft;
  unsigned long ireal;
  unsigned ix;
  unsigned iy;
  unsigned iz;
} POS;
static unsigned g_extx;                      /* map extent */
static unsigned g_exty;                      /* map extent */
static unsigned g_extz;                      /* map extent */
static unsigned g_extx_half;                 /* half map extent */
static unsigned g_exty_half;                 /* half map extent */
static unsigned g_extz_half;                 /* half map extent */
static unsigned long g_nvox;                 /* number of voxels */
static double g_width;                       /* voxel size in Angstroms */
static double g_origx;                       /* map origin */
static double g_origy;                       /* map origin */
static double g_origz;                       /* map origin */
static double *g_phi_lo;                     /* low resolution map */
static double *g_phi_hi;                     /* high resolution map */
static double g_norm_hi;                     /* corresponding normalization */
static double *g_phi_du;                     /* dummy map */
static double *g_phi_fi;                     /* filter kernel */
static double *g_phi_ga;                     /* low-pass (Gaussian) kernel */
static double *g_phi_fx;                     /* filtered low-pass kernel */
static unsigned g_ext_ga;                    /* low-pass kernel linear extent */
static unsigned g_ext_fi;                    /* filter kernel linear extent */
static unsigned g_ext_fx;                    /* filtered low-pass kernel linear extent */
static unsigned long g_nvox_ga;              /* low-pass kernel voxel count */
static unsigned long g_nvox_fi;              /* filter kernel voxel count */
static PDB *g_pdb_original;                  /* input coordinates in PDB format, fragile */
static PDB *g_pdb_save;                      /* coordinates in PDB format used for output */
static unsigned g_num_atoms;                 /* number of total PDB atoms */
static unsigned g_num_parts;                 /* number of parts for multi-body */
static unsigned *g_parts_num_atoms;          /* number of atoms in each multi-body part */
static double g_com_glob[3];                 /* global COM, PDB frame of reference */
static double g_cen_glob[3];                 /* global COM, colores style (map) frame of reference */
static double *g_com_x;                      /* local COM of each part, PDB frame of reference */
static double *g_com_y;          /* local COM of each part, PDB frame of reference */
static double *g_com_z;          /* local COM of each part, PDB frame of reference */
static double *g_cen_x;                      /* local COM of each part, colores style (map) frame of reference */
static double *g_cen_y;          /* local COM of each part, colores style (map) frame of reference */
static double *g_cen_z;          /* local COM of each part, colores style (map) frame of reference */
static double g_target_res;                  /* resolution in A */
static double g_target_ani;                  /* resolution anisotropy factor */
static int g_boost_option;                   /* boost method for steric clash reduction: scale or power */
static double g_boost_fact;                  /* boost threshold for steric clash reduction */
static double g_boost_limit;                 /* boost limit for steric clash reduction */
static double g_boost_par;                   /* boost factor or exponent for steric clash reduction */
static unsigned g_pow_max_iter;              /* Powell max number of iterations */
static char g_pow_mode;                      /* Powell option */
static double g_pow_tolerance;               /* Powell tolerance */
static double g_pow_delta_pos;               /* Powell initialization position  */
static double g_pow_delta_ang;               /* Powell initialization orientation */
static double g_low_cutoff;                  /* low res. map cutoff */
static int g_corr_mode;                      /* correlation option */
static double g_size_fac;                    /* grid size expansion factor for zero padding */
static unsigned g_ignored[3];                /* zero margin which can be safely ignored in fast kernel convolution */
static int g_pow_alg;                        /* correlation method for Powell optimization */
static FIT *g_best_fit;                      /* saved parameters and coefficients of best fits */
static FILE *g_outfile;                      /* frequently used output file */
static char *g_program = "collage";
static POS *g_inside_list;                   /* inside target positions */
static unsigned g_inside_num;                /* inside target number */


int main(int argc, char **argv)
{

  double test_fraction, test_mass_hi, test_mass_pdb;
  the_time itime, etime;
  time_t seed_time;
  double *pow_init;
  unsigned i, j, p, lindex, uindex;
  double sigma1d;
  double *phi_ga_save, *phi_fx_save;
  double sigma_factor = 0;
  unsigned long nvox_ga_save;
  char corr_string[70];
  char out_string[20];
  unsigned ext_ga_save;
  unsigned zeropad[3];
  unsigned long m;
  double corr_hi_lo;
  unsigned long q;
  double test_pow_corr, test_pow_corr_orig;
  double test_pow_time, test_pow_time_orig;
  int test_pow_method;
  POWARG *pow_args;
  POWSCR *pow_scratch;
  double norm_hi, norm_lo;

  sgenrand((unsigned long) time(&seed_time));

  draw_line();
  read_options(argc, &(*argv)); /* sets a variety of global variables */

  draw_line();
  get_low_map(argv[1]);
  draw_line();
  get_structures_and_partition(argv);
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

  /* zeropad low-resolution map as desired to accommodate the PDB */

  /* add padding */
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

  do_vect(&g_phi_lo, g_nvox);
  cp_vect_destroy(&g_phi_lo, &g_phi_du, g_nvox);

  /* memory allocation */
  do_vect(&g_phi_hi, g_nvox);
  do_vect(&g_phi_du, g_nvox);


  /* ************************************************************************************************ */
  /* the following global test calculation uses a single centered structure as in colacor version 2.5 */
  /* this is just a test calculation, no boost is applied                                             */
  /* ************************************************************************************************ */

  /* center atomic system globally */
  translate(g_pdb_original, g_pdb_original, g_num_atoms, -g_com_glob[0], -g_com_glob[1], -g_com_glob[2]);

  printf("collage> Projecting probe structure to lattice...\n");
  project_mass(&g_phi_hi, g_nvox, g_width, g_width, g_width * g_target_ani, g_extx, g_exty, g_extz, g_pdb_original, g_num_atoms, g_cen_glob, g_ignored);

  if (g_num_parts > 1 && g_boost_par != 1) printf("collage> Note: -boost arguments will be ignored in this initial test.\n");

  /* sum up mass in overlap region */
  printf("collage> Computing fraction of PDB contained within the map (above cutoff density) ...\n");
  test_mass_pdb = calc_mass(g_pdb_original, g_num_atoms);
  if (test_mass_pdb == 0) {
    printf("collage> Error: Input PDB has zero mass.\n");
    exit(1);
  }
  test_mass_hi = 0;
  for (m = 0; m < g_nvox; m++) {
    if (*(g_phi_lo + m) > 0) test_mass_hi += *(g_phi_hi + m);
  }
  test_fraction = test_mass_hi / test_mass_pdb; /* this is safe here based on earlier checks of structure */
  printf("collage> Overlap fraction: %15.7E\n", test_fraction);
  if (test_fraction < 0.5) printf("collage> Warning: Less than half of the input PDB is contained within the map!\n");

  printf("collage> Applying filters to target and probe maps...\n");
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


  /* colores style normalization and correlation coefficient */
  printf("collage> Normalizing target and probe maps...\n");
  norm_hi = calc_norm(g_phi_hi, g_nvox);
  norm_lo = calc_norm(g_phi_lo, g_nvox);
  normalize(g_phi_hi, g_nvox, norm_hi);
  normalize(g_phi_lo, g_nvox, norm_lo);
  g_norm_hi = norm_hi; /* save for special case of single-molecule one step algorithm */
  printf("collage> Target and probe maps:\n");
  print_map_info(g_phi_lo, g_nvox);
  print_map_info(g_phi_hi, g_nvox);
  corr_hi_lo = 0;
  for (m = 0; m < g_nvox; m++) corr_hi_lo += (*(g_phi_hi + m)) * (*(g_phi_lo + m));
  corr_hi_lo *= 1.0 / (1.0 * g_nvox);
  printf("collage> Initial correlation coefficient: %15.7E\n", corr_hi_lo);


  /* ********************************************************** */
  /* the following optimization requires locally centered parts */
  /* ********************************************************** */

  if (g_pow_mode) {

    /* undo global centering and apply local centering */
    lindex = 0;
    for (i = 0; i < g_num_parts; ++i) {
      uindex = lindex + g_parts_num_atoms[i];
      for (j = lindex; j < uindex; ++j) {
        g_pdb_original[j].x += (g_com_glob[0] - g_com_x[i]);
        g_pdb_original[j].y += (g_com_glob[1] - g_com_y[i]);
        g_pdb_original[j].z += (g_com_glob[2] - g_com_z[i]);
      }
      lindex = uindex;
    }

    /* initialize some parameters */
    draw_line();
    g_best_fit = (FIT *) alloc_vect((g_num_parts), sizeof(FIT));
    for (i = 0; i < g_num_parts; ++i) {
      g_best_fit[i].pos[0] = 0.0;
      g_best_fit[i].pos[1] = 0.0;
      g_best_fit[i].pos[2] = 0.0;
      g_best_fit[i].euler[0] = 0.0;
      g_best_fit[i].euler[1] = 0.0;
      g_best_fit[i].euler[2] = 0.0;
      g_best_fit[i].score = corr_hi_lo;
    }
    g_extx_half = (g_extx - 1) / 2;
    g_exty_half = (g_exty - 1) / 2;
    g_extz_half = (g_extz - 1) / 2;
    create_inside_molecule_poslist(g_phi_lo);
    printf("collage> Powell conjugate gradient maximization.\n");

    /* set Powell start directions*/
    do_vect(&pow_init, 6 * g_num_parts);
    if (g_pow_delta_pos < 0.0) {
      g_pow_delta_pos = g_width * 0.25;
    }
    if (g_pow_delta_ang < 0.0) {
      g_pow_delta_ang = 3.5 * ROT_CONV;
    } else if (g_pow_delta_ang > 10)  {
      g_pow_delta_ang = 10.0 * ROT_CONV;
    } else  {
      g_pow_delta_ang *= ROT_CONV;
    }
    for (i = 0; i < g_num_parts; ++i) for (j = 0; j < 3; ++j) {
        pow_init[6 * i + j] = g_pow_delta_pos;
        pow_init[6 * i + j + 3] = g_pow_delta_ang;
      }

    /* determine best method to calculate correlation */
    if (g_pow_alg == 0 || (g_pow_alg == 3 && g_num_parts > 1)) {
      test_pow_method = 1;
      /* prepare scratch space */
      pow_scratch = (POWSCR *) alloc_vect(sizeof(POWSCR), 1);
      pow_scratch->pdb = (PDB *) alloc_vect(g_num_atoms, sizeof(PDB));
      for (q = 0; q < g_num_atoms; ++q) *(pow_scratch->pdb + q) = *(g_pdb_original + q);
      do_vect(&pow_scratch->phi_du, g_nvox);
      do_vect(&pow_scratch->phi_hi, g_nvox);
      if (g_pow_alg == 3 && g_num_parts > 1) printf("collage> Warning: One-step algorithm not supported for multi-body docking, \n");
      printf("collage> Determining most efficient correlation algorithm based on convergence and time...\n");
      g_pow_alg = 1;
      itime = get_the_time();
      test_pow_corr_orig = powell_correlation_multi(pow_init, pow_scratch);
      test_pow_corr_orig = powell_correlation_multi(pow_init, pow_scratch);
      etime = get_the_time();
      test_pow_time_orig = (double)time_to_sec(time_diff(etime, itime)) / 2.0;
      printf("collage>   Original algorithm: Correlation = %10.8f  Time = %s\n",
             -test_pow_corr_orig,
             smart_sprint_time(test_pow_time_orig));
      g_pow_alg = 2;
      itime = get_the_time();
      test_pow_corr = powell_correlation_multi(pow_init, pow_scratch);
      test_pow_corr = powell_correlation_multi(pow_init, pow_scratch);
      etime = get_the_time();
      test_pow_time = (double)time_to_sec(time_diff(etime, itime)) / 2.0;
      printf("collage>   Masked algorithm:   Correlation = %10.8f  Time = %s\n",
             -test_pow_corr,
             smart_sprint_time(test_pow_time));
      if ((fabs(test_pow_corr - test_pow_corr_orig) < g_pow_tolerance) &&
          (test_pow_time < test_pow_time_orig)) {
        test_pow_time_orig = test_pow_time;
        test_pow_method = 2;
      }
      if (g_num_parts == 1) {
        g_pow_alg = 3;
        itime = get_the_time();
        test_pow_corr = powell_correlation_multi(pow_init, pow_scratch);
        test_pow_corr = powell_correlation_multi(pow_init, pow_scratch);
        etime = get_the_time();
        test_pow_time = (double)time_to_sec(time_diff(etime, itime)) / 2.0;
        printf("collage>   One-step algorithm: Correlation = %10.8f  Time = %s\n",
               -test_pow_corr,
               smart_sprint_time(test_pow_time));
        if ((fabs(test_pow_corr - test_pow_corr_orig) <= g_pow_tolerance) &&
            (test_pow_time < test_pow_time_orig)) {
          test_pow_time_orig = test_pow_time;
          test_pow_method = 3;
        }
      }
      g_pow_alg = test_pow_method;
      /* clean up scratch space */
      free_vect_and_zero_ptr(&(pow_scratch->phi_hi));
      free_vect_and_zero_ptr(&(pow_scratch->phi_du));
      free_vect_and_zero_ptr(&(pow_scratch->pdb));
      free_vect_and_zero_ptr(&pow_scratch);
    }

    /* check if steric exclusion boosting of density, and modify g_pow_alg to select case */
    if (g_boost_option == 0 && g_boost_par != 1 && g_num_parts > 1) g_pow_alg += 3;
    if (g_boost_option == 1 && g_boost_par != 1 && g_num_parts > 1) g_pow_alg += 5;

    /* print user info on selected correlation method */
    switch (g_pow_alg) {
      case 1:
        // Original three step code (scales based on g_nvox)
        printf("collage> Using original three-step correlation function.\n");
        break;
      case 2:
        // Masked three step (scales based on g_nvox)
        printf("collage> Using masked three-step correlation function.\n");
        break;
      case 3:
        // One step code (scales based on number of atoms in the high res structure)
        printf("collage> Using one-step correlation function.\n");
        break;
      case 4: /* with scale boost, internal use only */
        // Original three step code (scales based on g_nvox) with boost
        printf("collage> Using original three-step correlation function with steric exclusion (scale) boost.\n");
        break;
      case 5: /* with scale boost, internal use only */
        // Masked three step (scales based on g_nvox) with boost
        printf("collage> Using masked three-step correlation function with steric exclusion (scale) boost.\n");
        break;
      case 6: /* with power boost, internal use only */
        // Original three step code (scales based on g_nvox) with boost
        printf("collage> Using original three-step correlation function with steric exclusion (power) boost.\n");
        break;
      case 7: /* with power boost, internal use only */
        // Masked three step (scales based on g_nvox) with boost
        printf("collage> Using masked three-step correlation function with steric exclusion (power) boost.\n");
        break;
      default:
        fprintf(stderr, "collage> Error: Did not understand Powell correlation method, check -pwcorr\n");
        exit(1);
    }

    /* save Powell start time and initialize log file */
    itime = get_the_time();
    g_outfile = fopen("cge_powell.log", "w");
    if (g_outfile == NULL) {
      error_open_filename(80220, g_program, "cge_powell.log");
    }
    fprintf(g_outfile, "# This file contains information about positions and Euler angles during the Powell optimization of input PDBs.\n");
    fprintf(g_outfile, "# Deviations in Angstrom, Euler angles in degrees.\n");
    if (g_target_ani != 1.0) fprintf(g_outfile, "# The coordinates in the z-direction will be compressed by the anisotropy factor %.3f\n", g_target_ani);
    fprintf(g_outfile, "# More info (files and options) can be found in the output PDB file headers.\n");
    fprintf(g_outfile, "\n");

    printf("collage> \n");
    printf("collage> Performing Powell conjugate gradient maximization...\n");
    printf("collage> \n");

    /* now run Powell optimization */
    pow_args = (POWARG *) alloc_vect(sizeof(POWARG), 1);
    pow_args->iter = 0;
    pow_args = (POWARG *) alloc_vect(sizeof(POWARG), 1);
    pow_args->pow_init = pow_init;
    powell_optimization_multi((void *) pow_args, argv);
    if (g_target_ani != 1.0) printf("collage> The coordinates in the z-direction are compressed by the anisotropy factor %.3f\n", g_target_ani);
    etime = get_the_time();
    printf("collage> Total optimization time: %s \n", smart_sprint_time(time_to_sec(time_diff(etime, itime))));
    fclose(g_outfile);
    draw_line();

    /* write result to multiple PDB files */
    lindex = 0;
    for (i = 0; i < g_num_parts; ++i) {
      sprintf(out_string, "cge_%03i.pdb", i + 1);
      printf("collage> Writing body number %i to file %s. \n", i + 1, out_string);
      g_outfile = fopen(out_string, "w");
      if (g_outfile == NULL) {
        error_open_filename(80320, g_program, out_string);
      }
      /* PDB file header */
      fprintf(g_outfile, "REMARK\nREMARK    File name %s\n", out_string);
      fprintf(g_outfile, "REMARK    Low-resolution fit of structure %s into map %s\n", argv[i + 2], argv[1]);
      if (g_num_parts > 1) fprintf(g_outfile, "REMARK    Part of a multi-body fitting with %i components.\n", g_num_parts);
      else fprintf(g_outfile, "REMARK    Single molecule fitting.\n");
      fprintf(g_outfile, "REMARK    \n");
      fprintf(g_outfile, "REMARK    Resolution anisotropy factor (Z vs. XY): %f\n", g_target_ani);
      fprintf(g_outfile, "REMARK    Resolution: %f Angstrom; density cutoff: %f\n", g_target_res, g_low_cutoff);
      fprintf(g_outfile, "%s\n", corr_string);
      if (g_num_parts > 1) fprintf(g_outfile, "REMARK    Total correlation coefficient for entire multi-body system: %f\n", g_best_fit[i].score);
      else   fprintf(g_outfile, "REMARK    correlation coefficient: %f\n", g_best_fit[i].score);
      fprintf(g_outfile, "REMARK    \n");
      fprintf(g_outfile, "REMARK    Euler angles (Psi, Theta, Phi):  %7.3f %7.3f %7.3f deg.\n",
              g_best_fit[i].euler[0] / ROT_CONV,
              g_best_fit[i].euler[1] / ROT_CONV,
              g_best_fit[i].euler[2] / ROT_CONV);
      fprintf(g_outfile, "REMARK    Center position (X,Y,Z):            %7.3f %7.3f %7.3f A\n",
              g_best_fit[i].pos[0] + g_com_x[i],
              g_best_fit[i].pos[1] + g_com_y[i],
              g_target_ani * g_best_fit[i].pos[2] + g_com_z[i]);
      fprintf(g_outfile, "REMARK    \n");
      fclose(g_outfile);
      /* PDB file data */
      g_pdb_save = (PDB *) alloc_vect(g_parts_num_atoms[i], sizeof(PDB));
      for (p = 0; p < g_parts_num_atoms[i]; ++p) *(g_pdb_save + p) = *(g_pdb_original + p + lindex);
      uindex = lindex + g_parts_num_atoms[i];
      multi_rot_euler_trans(lindex, uindex, lindex, g_pdb_original, g_pdb_save, g_best_fit[i].euler[0], g_best_fit[i].euler[1], g_best_fit[i].euler[2], g_best_fit[i].pos[0] + g_com_x[i], g_best_fit[i].pos[1] + g_com_y[i], g_target_ani * g_best_fit[i].pos[2] + g_com_z[i]);
      append_pdb(out_string, g_parts_num_atoms[i], g_pdb_save);
      lindex = uindex;
      free_vect_and_zero_ptr(&g_pdb_save);
    }
  }

  printf("collage> Writing Powell log to file %s. \n", "cge_powell.log");
  draw_line();
  printf("collage> All done!!!\n");
  return 0;
}







/*====================================================================*/
static void draw_line()
{
  printf("_____________________________________________________________________________\n");
}


/*====================================================================*/
static void get_low_map(char *file_name)
{
  /* reads low resolution map from file_name */

  double cut_width;
  double new_width;

  printf("collage> Processing low-resolution map.\n");

  read_vol(file_name, &g_width, &g_origx, &g_origy, &g_origz, &g_extx, &g_exty, &g_extz, &g_phi_lo);
  g_nvox = g_extx * g_exty * g_extz;


  /* if spacing is too wide adjust resolution */
  if (g_width > g_target_res * 0.7) {
    g_target_res = 2.0 * g_width;
    printf("collage> Warning: Insufficient spatial sampling (voxel spacing too wide) for initially assigned map resolution.\n");
    printf("collage> Target resolution adjusted to 2x voxel spacing, %f.\n", g_target_res);
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

  /* shrink map about non-zero density (inadvertently also resizes to odd intervals, as needed by FFT in colores, which is technically not required here) */
  shrink_margin(&g_phi_du, &g_extx, &g_exty, &g_extz, &g_origx, &g_origy, &g_origz, &g_nvox,
                g_phi_lo, g_extx, g_exty, g_extz, g_origx, g_origy, g_origz,
                g_width, g_width, g_width * g_target_ani);

  cp_vect_destroy(&g_phi_lo, &g_phi_du, g_nvox);

  /* set boost limit above which density will be amplified */
  g_boost_limit = g_boost_fact * max_projection_density(g_width);
  print_map_info(g_phi_lo, g_nvox);
}


/*====================================================================*/
void multi_rot_euler_trans(unsigned lindex, unsigned uindex, unsigned oindex, PDB *pdb_original, PDB *pdb_moved, double psi, double theta, double phi, double x0, double y0, double z0)
{
  /* rotates multi-body part by Euler angles and applies translation, no sanity check for indices */

  unsigned id;
  double rot_matrix[3][3], currx, curry, currz;

  get_rot_matrix(rot_matrix, psi, theta, phi);

  for (id = lindex; id < uindex; ++id) {
    currx = pdb_original[id].x;
    curry = pdb_original[id].y;
    currz = pdb_original[id].z;
    pdb_moved[id - oindex].x = currx * rot_matrix[0][0] +
                               curry * rot_matrix[0][1] +
                               currz * rot_matrix[0][2] + x0;
    pdb_moved[id - oindex].y = currx * rot_matrix[1][0] +
                               curry * rot_matrix[1][1] +
                               currz * rot_matrix[1][2] + y0;
    pdb_moved[id - oindex].z = currx * rot_matrix[2][0] +
                               curry * rot_matrix[2][1] +
                               currz * rot_matrix[2][2] + z0;
  }
}


/*====================================================================*/
static void get_structures_and_partition(char **argv)
{
  /* reads PDB files and assigns partitioning, appends them to single structure which is locally centered */
  /* replaces colores get_centered_structure_and_radius for multi-body support */

  int curr_arg, slen;
  unsigned curr_num_atoms, i;
  PDB *curr_pdb;
  DIR *mydir; /* see dirent.h */
  struct dirent *entry = NULL; /* see dirent.h */

  /* determine number of PDB files and append them to g_pdb_original */
  printf("collage> Processing atomic structure(s). PDB input file(s) detected:\n");
  g_num_parts = 0;
  g_num_atoms = 0;
  g_com_x = NULL;
  g_com_y = NULL;
  g_com_z = NULL;
  g_parts_num_atoms = NULL;

  if ((mydir = opendir(argv[2])) != NULL) { /* have a directory of PDB files */
    while ((entry = readdir(mydir))) {
      slen = (int)strlen(entry->d_name);
      if (slen > 4 && (strstr(entry->d_name + slen - 4, ".pdb") != NULL || strstr(entry->d_name + slen - 4, ".PDB") != NULL)) {
        ++g_num_parts;
        read_pdb_silent(entry->d_name, &curr_num_atoms, &curr_pdb);
        printf("collage>   %i. %s (%i atoms)\n", g_num_parts, entry->d_name, curr_num_atoms);
        g_num_atoms += curr_num_atoms;
        g_pdb_original = (PDB *) realloc(g_pdb_original, (g_num_atoms) * sizeof(PDB));
        if (g_pdb_original == NULL) error_memory_allocation(80621, g_program);
        for (i = 0; i < curr_num_atoms; i++) g_pdb_original[g_num_atoms - curr_num_atoms + i] = curr_pdb[i];
        g_com_x = (double *) realloc(g_com_x, (g_num_parts) * sizeof(double));
        g_com_y = (double *) realloc(g_com_y, (g_num_parts) * sizeof(double));
        g_com_z = (double *) realloc(g_com_z, (g_num_parts) * sizeof(double));
        if (g_com_x == NULL || g_com_y == NULL || g_com_z == NULL) error_memory_allocation(80622, g_program);
        calc_center(curr_pdb, curr_num_atoms, g_com_x + g_num_parts - 1, g_com_y + g_num_parts - 1, g_com_z + g_num_parts - 1);
        free_vect_and_zero_ptr(&curr_pdb);
        g_parts_num_atoms = (unsigned *) realloc(g_parts_num_atoms, (g_num_parts) * sizeof(unsigned));
        if (g_parts_num_atoms == NULL) error_memory_allocation(80623, g_program);
        g_parts_num_atoms[g_num_parts - 1] = curr_num_atoms;
      }
    }
    closedir(mydir);
  } else { /* have a list of PDB file arguments */
    curr_arg = 2;
    while (fopen(argv[curr_arg], "r") != NULL) {
      ++g_num_parts;
      read_pdb_silent(argv[curr_arg], &curr_num_atoms, &curr_pdb);
      printf("collage>   %i. %s (%i atoms)\n", g_num_parts, argv[curr_arg], curr_num_atoms);
      g_num_atoms += curr_num_atoms;
      g_pdb_original = (PDB *) realloc(g_pdb_original, (g_num_atoms) * sizeof(PDB));
      if (g_pdb_original == NULL) error_memory_allocation(80624, g_program);
      for (i = 0; i < curr_num_atoms; i++) g_pdb_original[g_num_atoms - curr_num_atoms + i] = curr_pdb[i];
      g_com_x = (double *) realloc(g_com_x, (g_num_parts) * sizeof(double));
      g_com_y = (double *) realloc(g_com_y, (g_num_parts) * sizeof(double));
      g_com_z = (double *) realloc(g_com_z, (g_num_parts) * sizeof(double));
      if (g_com_x == NULL || g_com_y == NULL || g_com_z == NULL) error_memory_allocation(80625, g_program);
      calc_center(curr_pdb, curr_num_atoms, g_com_x + g_num_parts - 1, g_com_y + g_num_parts - 1, g_com_z + g_num_parts - 1);
      free_vect_and_zero_ptr(&curr_pdb);
      g_parts_num_atoms = (unsigned *) realloc(g_parts_num_atoms, (g_num_parts) * sizeof(unsigned));
      if (g_parts_num_atoms == NULL) error_memory_allocation(80626, g_program);
      g_parts_num_atoms[g_num_parts - 1] = curr_num_atoms;
      ++curr_arg;
    }
  }

  /* sanity checks */
  if (g_parts_num_atoms == NULL) {
    fprintf(stderr, "collage> Error: No input PDB structures detected.\n");
    exit(1);
  } else {
    if (g_num_parts == 1 && g_boost_par != 1) printf("collage> Warning: -boost arguments will be ignored for single input structure.\n");
  }

  /* set (PDB and colores style) global centroids: */
  calc_center(g_pdb_original, g_num_atoms, g_com_glob, (g_com_glob + 1), (g_com_glob + 2));
  g_cen_glob[0] = g_com_glob[0] - 0.5 * g_extx * g_width - g_origx;
  g_cen_glob[1] = g_com_glob[1] - 0.5 * g_exty * g_width - g_origy;
  g_cen_glob[2] = g_com_glob[2] - 0.5 * g_extz * g_width * g_target_ani - g_origz;

  /* set colores style local centroids: */
  do_vect(&g_cen_x, g_num_parts);
  do_vect(&g_cen_y, g_num_parts);
  do_vect(&g_cen_z, g_num_parts);
  for (i = 0; i < g_num_parts; ++i) {
    g_cen_x[i] = g_com_x[i] - 0.5 * g_extx * g_width - g_origx;
    g_cen_y[i] = g_com_y[i] - 0.5 * g_exty * g_width - g_origy;
    g_cen_z[i] = g_com_z[i] - 0.5 * g_extz * g_width * g_target_ani - g_origz;
  }
}


/*====================================================================*/
static void read_options(int argc, char **argv)
{
  /* print usage info and read options from input arguments */

  char *pos;
  int i;
  char option_string[2048];

  if (argc < 2) {
    /* print usage and options info */

    printf("collage> USAGE:   collage <Density map> <PDB files(s) or Directory>  -<options>\n");
    printf("collage>\n");
    printf("collage> OPTIONS:\n");
    printf("collage>\n");
    printf("collage>  -res <float>           Target resolution in A [default: -res 15]\n");
    printf("collage>\n");
    printf("collage>  -ani <float>           Resolution anisotropy factor [default: -ani 1]\n");
    printf("collage>\n");
    printf("collage>  -cutoff <float>        Density map cutoff value [default:-cutoff 0.0]\n");
    printf("collage>\n");
    printf("collage>  -corr <int>            Correlation method:\n");
    printf("collage>\n");
    printf("collage>          -corr 0 -->    Standard cross correlation [default]\n");
    printf("collage>          -corr 1 -->    Laplacian filtered correlation \n");
    printf("collage>\n");
    printf("collage>  -nopowell              Powell maximization Off [default On]\n");
    printf("collage>\n");
    printf("collage>  -pwcorr <int>          Powell correlation algorithm options\n");
    printf("collage>                         [default: -pwcorr 0]\n");
    printf("collage>\n");
    printf("collage>                         0: Determined at runtime\n");
    printf("collage>                         1: Original three-step code\n");
    printf("collage>                         2: Three-step code with mask applied\n");
    printf("collage>                         3: One-step code (for single PDB only)\n");
    printf("collage>\n");
    printf("collage>  -boost <int> <float>x2 (Multi-fragment) steric exclusion option, threshold, boost parameter\n");
    printf("collage>                         [default: none] \n");
    printf("collage>\n");
    printf("collage>  -pwti <float> <int>    Powell tolerance & max iterations\n");
    printf("collage>                         [default: -pwti 1e-6 50]\n");
    printf("collage>\n");
    printf("collage>  -pwtr <float>x2        Trans & Rot initial step size\n");
    printf("collage>                         [default: .25 voxel spacing 3.75 angular sampling]\n");
    printf("collage>\n");
    draw_line();
    exit(0);
  }

  /* now read options from arguments and assign global variables */

  sprintf(option_string, " ");
  printf("collage> Options read:\n");
  for (i = 1; i < argc; i++)
    sprintf(option_string, "%s %s", option_string, argv[i]);

  /* target resolution */
  g_target_res = 15.0;
  if ((pos = (char *)strstr(option_string, " -res")) != NULL)
    sscanf(pos + 5, "%lf", &g_target_res);
  if (g_target_res < 0) {
    error_resolution_range(80400, g_program);
  }
  printf("collage> Target resolution %3.3f\n", g_target_res);

  /* resolution anisotropy */
  g_target_ani = 1.0;
  if ((pos = (char *)strstr(option_string, " -ani")) != NULL)
    sscanf(pos + 5, "%lf", &g_target_ani);
  if (g_target_ani < 0.001 || g_target_ani > 1000) {
    error_anisotropy_range(80410, g_program);
  }
  printf("collage> Resolution anisotropy %3.3f\n", g_target_ani);

  /* Powell correlation algorithm */
  g_pow_alg = 0;
  if ((pos = (char *)strstr(option_string, " -pwcorr")) != NULL) sscanf(pos + 8, "%d", &g_pow_alg);
  if (g_pow_alg == 0) {
    printf("collage> Powell correlation algorithm determined automatically\n");
  } else {
    printf("collage> Desired Powell correlation algorithm: %d\n", g_pow_alg);
  }
  if (g_pow_alg > 3) {
    fprintf(stderr, "collage> Error: Invalid Powell correlation algorithm [e.c. 80440]\n");
    exit(80440);
  }

  /* non-negative cutoff of the low-resolution map */
  g_low_cutoff = 0;
  if ((pos = (char *)strstr(option_string, " -cutoff")) != NULL)
    sscanf(pos + 8, "%lf", &g_low_cutoff);
  if (g_low_cutoff < 0) {
    g_low_cutoff = 0;
    printf("collage> Low-resolution map cutoff must be non-negative, assigning %3.3f\n", g_low_cutoff);
  } else printf("collage> Low-resolution map cutoff %3.3f\n", g_low_cutoff);

  /* Powell maximization toggle */
  g_pow_mode = 1;
  if ((pos = (char *)strstr(option_string, " -nopowell")) != NULL) {
    g_pow_mode = 0;
    printf("collage> Powell maximization OFF\n");
  } else  printf("collage> Powell maximization ON\n");

  /* colores style grid size expansion for extra zero padding (for debugging purposes) */
  g_size_fac = 0.0;
  if ((pos = (char *)strstr(option_string, " -sizef")) != NULL) {
    printf("collage> Warning: -sizef has been discontinued, use -sizef_debug if you want to test grid size expansion for debugging. \n");
    printf("collage> Grid size expansion factor %1.3f (thickness of additional zero layer as fraction of map dimensions)\n", g_size_fac);
  }
  if ((pos = (char *)strstr(option_string, " -sizef_debug")) != NULL) {
    sscanf(pos + 7, "%lf", &g_size_fac);
    if (g_size_fac < 0) {
      printf("collage> Warning, grid size expansion factor must be positive, will be set to zero\n");
      g_size_fac = 0.0;
    }
    printf("collage> Grid size expansion factor %1.3f (thickness of additional zero layer as fraction of map dimensions)\n", g_size_fac);
  }

  /* correlation mode */
  g_corr_mode = 0;
  if ((pos = (char *)strstr(option_string, " -corr")) != NULL)
    sscanf(pos + 6, "%d", &g_corr_mode);
  switch (g_corr_mode) {
    case 0:
      printf("collage> Standard cross correlation\n");
      break;
    case 1:
      printf("collage> Laplacian filtered correlation\n");
      break;
    default:
      error_option(80420, g_program);
  }

  /* steric exclusion boost parameters */
  g_boost_fact = 0.9;
  g_boost_par = 1.0;
  g_boost_option = 0;
  if ((pos = (char *)strstr(option_string, " -boost")) != NULL)
    sscanf(pos + 7, "%d %lf %lf", &g_boost_option, &g_boost_fact, &g_boost_par);
  if (g_boost_option < 0 || g_boost_option > 1 || g_boost_fact < 0) {
    printf("collage> Multi-fragment steric exclusion parameters out of range, set to default (no boost)\n");
    g_boost_fact = 0.9;
    g_boost_par = 1.0;
    g_boost_option = 0;
  }
  if (g_boost_par != 1) printf("collage> Multi-fragment steric exclusion boost option %d, threshold %6.3f parameter %6.3f\n", g_boost_option, g_boost_fact, g_boost_par);


  /* Powell parameters */
  g_pow_tolerance = 1.0e-6;
  g_pow_max_iter = 50;
  if ((pos = (char *)strstr(option_string, " -pwti")) != NULL)
    sscanf(pos + 6, "%lf %d", &g_pow_tolerance, &g_pow_max_iter);
  printf("collage> Powell tolerance %1.2E  Max iterations %d\n", g_pow_tolerance, g_pow_max_iter);
  g_pow_delta_pos = -9;
  g_pow_delta_ang = -9;
  if ((pos = (char *)strstr(option_string, " -pwtr")) != NULL)
    sscanf(pos + 6, "%lf %lf", &g_pow_delta_pos, &g_pow_delta_ang);
  if (g_pow_delta_pos < 0 || g_pow_delta_ang < 0)
    printf("collage> Powell trans & rot initial step sizes set to default values\n");
  else printf("collage> Powell trans step size %1.3f, rot step size %4.3f \n", g_pow_delta_pos, g_pow_delta_ang);
}


/*====================================================================*/
static double powell_correlation_multi(double *vect6n, POWSCR *scratch)
{
  /* computes current correlation value for Powell based on 6n dim coordinates */
  /* Note: returns NEGATIVE value since we wish to maximize the score */
  /* compare to colores 'powell_correlation' */

  unsigned id, uindex, lindex, id6;

  /* local (multi-body) rot and trans according to vect6n, with local shift correction */
  lindex = 0;
  for (id = 0; id < g_num_parts; ++id) {
    id6 = id * 6;
    uindex = lindex + g_parts_num_atoms[id];
    multi_rot_euler_trans(lindex, uindex, 0, g_pdb_original, scratch->pdb, vect6n[id6 + 3], vect6n[id6 + 4], vect6n[id6 + 5], vect6n[id6] + g_cen_x[id], vect6n[id6 + 1] + g_cen_y[id], vect6n[id6 + 2] + g_cen_z[id]);
    lindex = uindex;
  }

  /* scratch shift must not contain any local (multi-body specific) information so we set it to zero */
  scratch->curr_shift[0] = 0;
  scratch->curr_shift[1] = 0;
  scratch->curr_shift[2] = 0;

  scratch->corr_hi_lo = 0;
  switch (g_pow_alg) {
    case 1:
      /* original three step code (scales based on g_nvox) */
      project_mass(&scratch->phi_hi, g_nvox, g_width, g_width, g_width * g_target_ani,
                   g_extx, g_exty, g_extz, scratch->pdb, g_num_atoms, scratch->curr_shift, g_ignored);
      convolve_kernel_inside_fast(&scratch->phi_du, scratch->phi_hi, g_extx, g_exty, g_extz,
                                  g_phi_fx, g_ext_fx, 1.0, g_ignored);
      normalize(scratch->phi_du, g_nvox, calc_norm(scratch->phi_du, g_nvox));
      for (scratch->l = 0; scratch->l < g_nvox; scratch->l++) {
        scratch->corr_hi_lo += (*(scratch->phi_du + scratch->l)) * (*(g_phi_lo + scratch->l));
      }
      break;
    case 2:
      /* modified three step code (uses mask when computing the correlation) */
      project_mass(&scratch->phi_hi, g_nvox, g_width, g_width, g_width * g_target_ani,
                   g_extx, g_exty, g_extz, scratch->pdb, g_num_atoms, scratch->curr_shift, g_ignored);
      convolve_kernel_inside_fast(&scratch->phi_du, scratch->phi_hi, g_extx, g_exty, g_extz,
                                  g_phi_fx, g_ext_fx, 1.0, g_ignored);
      normalize(scratch->phi_du, g_nvox, calc_norm(scratch->phi_du, g_nvox));
      for (scratch->l = 0; scratch->l < g_inside_num; scratch->l++) {
        scratch->corr_hi_lo += (*(scratch->phi_du + g_inside_list[scratch->l].ireal)) * (*(g_phi_lo + g_inside_list[scratch->l].ireal));
      }
      break;
    case 3:
      /* one step code (scales based on number of atoms in the high res structure) */
      /* this only works for single body because pmckc cannot normalize projected map on the fly */
      project_mass_convolve_kernel_corr(g_width, g_width, g_width * g_target_ani,
                                        g_extx, g_exty, g_extz, scratch->pdb, g_num_atoms, scratch->curr_shift,
                                        g_phi_fx, g_ext_fx, g_norm_hi, g_ignored, g_phi_lo, &scratch->corr_hi_lo);
      break;
    case 4: /* with factor scaling boost, internal use only */
      /* original three step code (scales based on g_nvox) with boost */
      project_mass(&scratch->phi_hi, g_nvox, g_width, g_width, g_width * g_target_ani,
                   g_extx, g_exty, g_extz, scratch->pdb, g_num_atoms, scratch->curr_shift, g_ignored);
      boost_factor_high(scratch->phi_hi, g_nvox, g_boost_limit, g_boost_par);

      convolve_kernel_inside_fast(&scratch->phi_du, scratch->phi_hi, g_extx, g_exty, g_extz,
                                  g_phi_fx, g_ext_fx, 1.0, g_ignored);
      normalize(scratch->phi_du, g_nvox, calc_norm(scratch->phi_du, g_nvox));
      for (scratch->l = 0; scratch->l < g_nvox; scratch->l++) {
        scratch->corr_hi_lo += (*(scratch->phi_du + scratch->l)) * (*(g_phi_lo + scratch->l));
      }
      break;
    case 5: /* with factor scaling boost, internal use only */
      /* modified three step code (uses mask when computing the correlation) with boost */
      project_mass(&scratch->phi_hi, g_nvox, g_width, g_width, g_width * g_target_ani,
                   g_extx, g_exty, g_extz, scratch->pdb, g_num_atoms, scratch->curr_shift, g_ignored);
      boost_factor_high(scratch->phi_hi, g_nvox, g_boost_limit, g_boost_par);
      convolve_kernel_inside_fast(&scratch->phi_du, scratch->phi_hi, g_extx, g_exty, g_extz,
                                  g_phi_fx, g_ext_fx, 1.0, g_ignored);
      normalize(scratch->phi_du, g_nvox, calc_norm(scratch->phi_du, g_nvox));
      for (scratch->l = 0; scratch->l < g_inside_num; scratch->l++) {
        scratch->corr_hi_lo += (*(scratch->phi_du + g_inside_list[scratch->l].ireal)) * (*(g_phi_lo + g_inside_list[scratch->l].ireal));
      }
      break;
    case 6: /* with power boost, internal use only */
      /* original three step code (scales based on g_nvox) with boost */
      project_mass(&scratch->phi_hi, g_nvox, g_width, g_width, g_width * g_target_ani,
                   g_extx, g_exty, g_extz, scratch->pdb, g_num_atoms, scratch->curr_shift, g_ignored);
      boost_power_high(scratch->phi_hi, g_nvox, g_boost_limit, g_boost_par);

      convolve_kernel_inside_fast(&scratch->phi_du, scratch->phi_hi, g_extx, g_exty, g_extz,
                                  g_phi_fx, g_ext_fx, 1.0, g_ignored);
      normalize(scratch->phi_du, g_nvox, calc_norm(scratch->phi_du, g_nvox));
      for (scratch->l = 0; scratch->l < g_nvox; scratch->l++) {
        scratch->corr_hi_lo += (*(scratch->phi_du + scratch->l)) * (*(g_phi_lo + scratch->l));
      }
      break;
    case 7: /* with power boost, internal use only */
      /* modified three step code (uses mask when computing the correlation) with boost */
      project_mass(&scratch->phi_hi, g_nvox, g_width, g_width, g_width * g_target_ani,
                   g_extx, g_exty, g_extz, scratch->pdb, g_num_atoms, scratch->curr_shift, g_ignored);
      boost_power_high(scratch->phi_hi, g_nvox, g_boost_limit, g_boost_par);
      convolve_kernel_inside_fast(&scratch->phi_du, scratch->phi_hi, g_extx, g_exty, g_extz,
                                  g_phi_fx, g_ext_fx, 1.0, g_ignored);
      normalize(scratch->phi_du, g_nvox, calc_norm(scratch->phi_du, g_nvox));
      for (scratch->l = 0; scratch->l < g_inside_num; scratch->l++) {
        scratch->corr_hi_lo += (*(scratch->phi_du + g_inside_list[scratch->l].ireal)) * (*(g_phi_lo + g_inside_list[scratch->l].ireal));
      }
      break;
    default:
      fprintf(stderr, "collage> Error: Did not understand Powell correlation method, check -pwcorr\n");
      exit(1);
  }

  scratch->corr_hi_lo /= (double)g_nvox;
  return -scratch->corr_hi_lo;
}


/*====================================================================*/
static void powell_optimization_multi(void *args, char **argv)
{
  /* wrapper for multi-body Powell optimization, serialized but using thread-safe powell_r */

  POWARG        *arguments;
  double        *pow_init;
  int           pow_code;
  double        *pow_vect6n;
  double        pow_score;
  POWSCR        *pow_scratch = NULL;
  POWRES        *pow_results = NULL;
  POWITERRES    *curr_iter = NULL;
  POWITERRES    *last_iter = NULL;
  int           p, q, j;

  /* unpack arguments */
  arguments = (POWARG *) args;
  pow_init  = arguments->pow_init;

  /* create and assign 6n dim Powell vector to zero */
  do_vect(&pow_vect6n, 6 * g_num_parts);

  /* prepare scratch space */
  pow_scratch = (POWSCR *) alloc_vect(sizeof(POWSCR), 1);
  pow_scratch->pdb = (PDB *) alloc_vect(g_num_atoms, sizeof(PDB));
  for (q = 0; q < g_num_atoms; ++q) *(pow_scratch->pdb + q) = *(g_pdb_original + q);
  do_vect(&pow_scratch->phi_du, g_nvox);
  do_vect(&pow_scratch->phi_hi, g_nvox);

  /* prepare result struct */
  pow_results = (POWRES *) alloc_vect(sizeof(POWRES), 1);
  pow_results->head = NULL;
  pow_results->last = NULL;

  /* call Powell maximization function */
  powell_r(&pow_code, &pow_score, pow_vect6n, 6 * g_num_parts,
           powell_correlation_multi, g_pow_max_iter, g_outfile,
           g_pow_tolerance, pow_init,
           pow_scratch, pow_results);

  /* clean up scratch space */
  free_vect_and_zero_ptr(&(pow_scratch->phi_hi));
  free_vect_and_zero_ptr(&(pow_scratch->phi_du));
  free_vect_and_zero_ptr(&(pow_scratch->pdb));
  free_vect_and_zero_ptr(&pow_scratch);

  /* manage output */

  if (g_num_parts > 1) {
    printf("collage> Multi-body Powell\n");
    fprintf(g_outfile, "Multi-body Powell optimization\n");
  } else {
    printf("collage> Single-body Powell\n");
    fprintf(g_outfile, "Single-body Powell optimization\n");
  }

  /* output initial configuration info to both screen and log file */
  printf("collage>   X       Y       Z       Psi     Theta   Phi      Correlation   Iteration   PDB (Body) Nr.\n");
  fprintf(g_outfile, "   X       Y       Z       Psi     Theta   Phi      Correlation   Iteration   PDB\n");
  for (j = 0; j < g_num_parts; j++) {
    printf("collage> %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f    %10.7E Initial     %-7d\n",
           g_best_fit[j].pos[0],
           g_best_fit[j].pos[1],
           g_best_fit[j].pos[2]*g_target_ani,
           g_best_fit[j].euler[0] / ROT_CONV,
           g_best_fit[j].euler[1] / ROT_CONV,
           g_best_fit[j].euler[2] / ROT_CONV,
           g_best_fit[j].score, j + 1);
    fprintf(g_outfile, " %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f    %10.7E Initial     %s\n",
            g_best_fit[j].pos[0],
            g_best_fit[j].pos[1],
            g_best_fit[j].pos[2]*g_target_ani,
            g_best_fit[j].euler[0] / ROT_CONV,
            g_best_fit[j].euler[1] / ROT_CONV,
            g_best_fit[j].euler[2] / ROT_CONV,
            g_best_fit[j].score, argv[j + 2]);
  }

  /* write intermediate info to log file only */
  curr_iter = pow_results->head;
  while (curr_iter != NULL) {
    for (j = 0; j < g_num_parts; j++) {
      /*
      printf("collage> %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f    %10.7E %-7d     %-7d\n",
       curr_iter->res[6*j+0],
       curr_iter->res[6*j+1],
       curr_iter->res[6*j+2]*g_target_ani,
       curr_iter->res[6*j+3]/ROT_CONV,
       curr_iter->res[6*j+4]/ROT_CONV,
       curr_iter->res[6*j+5]/ROT_CONV,
       curr_iter->corr,
       curr_iter->iter,j+1);
      */
      fprintf(g_outfile, " %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f    %10.7E %-7d     %s\n",
              curr_iter->res[6 * j + 0],
              curr_iter->res[6 * j + 1],
              curr_iter->res[6 * j + 2]*g_target_ani,
              curr_iter->res[6 * j + 3] / ROT_CONV,
              curr_iter->res[6 * j + 4] / ROT_CONV,
              curr_iter->res[6 * j + 5] / ROT_CONV,
              curr_iter->corr,
              curr_iter->iter, argv[j + 2]);
    }
    curr_iter = curr_iter->next;
  }

  /* process final configuration info */
  for (j = 0; j < g_num_parts; j++) {
    for (p = 0; p < 3; p++) {
      g_best_fit[j].pos[p] = pow_results->last->res[6 * j + p];
      g_best_fit[j].euler[p] = fmod(pow_results->last->res[6 * j + p + 3] + 2 * PI, 2 * PI);
    }
    g_best_fit[j].score = pow_results->last->corr;
    remap_eulers((g_best_fit[j].euler + 0),
                 (g_best_fit[j].euler + 1),
                 (g_best_fit[j].euler + 2),
                 g_best_fit[j].euler[0],
                 g_best_fit[j].euler[1],
                 g_best_fit[j].euler[2],
                 0.0, 0.0, 0.0);
  }

  /* output final configuration info to both screen and log file */
  for (j = 0; j < g_num_parts; j++) {
    printf("collage> %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f    %10.7E Final       %-7d\n",
           g_best_fit[j].pos[0],
           g_best_fit[j].pos[1],
           g_best_fit[j].pos[2]*g_target_ani,
           g_best_fit[j].euler[0] / ROT_CONV,
           g_best_fit[j].euler[1] / ROT_CONV,
           g_best_fit[j].euler[2] / ROT_CONV,
           g_best_fit[j].score, j + 1);
    fprintf(g_outfile, " %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f    %10.7E Final       %s\n",
            g_best_fit[j].pos[0],
            g_best_fit[j].pos[1],
            g_best_fit[j].pos[2]*g_target_ani,
            g_best_fit[j].euler[0] / ROT_CONV,
            g_best_fit[j].euler[1] / ROT_CONV,
            g_best_fit[j].euler[2] / ROT_CONV,
            g_best_fit[j].score, argv[j + 2]);
  }


  /* clean up */
  printf("collage> Total number of Powell steps: %i \n", pow_results->last->iter);
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
static void create_inside_molecule_poslist(double *phi)
{
  /* calculates all positions in map that have positive density or are buried */
  /* called once in preprocessing stage of collage                            */
  /* writes list to g_inside_list, number to g_inside_num                     */

  unsigned curr, i;
  unsigned long q, s;
  unsigned check, ix, iy, iz;
  int ix2, iy2, iz2, erosion_shell_width;
  char *mask_inside, *mask_inside2;

  mask_inside = (char *) alloc_vect(g_nvox, sizeof(char));
  for (q = 0; q < g_nvox; q++) mask_inside[q] = 0;

  g_inside_num = 0;
  printf("collage> Identifying inside or buried voxels...\n");
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
  for (q = 0; q < g_nvox; q++) mask_inside2[q] = mask_inside[q];

  unsigned count = 0;

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
            count++;
            g_inside_num--;
          }
        }
      }
  free_vect_and_zero_ptr(&mask_inside);

  printf("collage> Found %d inside or buried voxels (out of a total of %lu).\n", g_inside_num, g_nvox);

  /* allocate inside poslist */
  g_inside_list = (POS *) alloc_vect(g_inside_num, sizeof(POS));

  /* now fill inside poslist */
  curr = 0;
  for (iz = 0; iz < g_extz; iz++)
    for (iy = 0; iy < g_exty; iy++)
      for (ix = 0; ix < g_extx; ix++) {
        q = ix + (g_extx) * (iy + g_exty * iz); /* actual lattice index */
        if (mask_inside2[q] == 1) { /* inside */
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
double max_projection_density(double width)
{
  /* empirical max density of structure projections to a lattice with voxel spacing 'width' */

  return (18901162820 * width * width * width * width * width -  413277197352 * width * width * width * width + 2599824151495 * width * width * width - 8584659866040 * width * width + 8276688322695 * width - 5841891340218) / (-394441476660.0);
}

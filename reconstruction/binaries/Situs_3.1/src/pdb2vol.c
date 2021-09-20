/*********************************************************************
 *                           P D B 2 V O L                            *
 **********************************************************************
 * Program is part of the Situs package (c) Willy Wriggers, 1998-2005 *
 * URL: situs.biomachina.org                                          *
 **********************************************************************
 *                                                                    *
 * Kernel convolution / Low-pass filtering                            *
 *                                                                    *
 **********************************************************************
 * See legal statement for terms of distribution                      *
 *********************************************************************/

#include "situs.h"
#include "lib_std.h"
#include "lib_pio.h"
#include "lib_pwk.h"
#include "lib_vio.h"
#include "lib_vwk.h"
#include "lib_err.h"

/* global variables */
static int g_extx,  g_exty;   /* kernel map */
static int g_extx2, g_exty2;  /* protein map */
static int g_extx3, g_exty3;  /* output map */

/* function declarations */
static int check_water(PDB *, int);
static int check_hydrogen(PDB *, int);
static int check_codebook(PDB *, int);
static int check_nondens(PDB *, int);
static unsigned long g1idz(int, int, int);
static unsigned long g2idz(int, int, int);
static unsigned long g3idz(int, int, int);

int main(int argc, char *argv[])
{
  double *phi = NULL, *phi2 = NULL, *phi3 = NULL;
  double mass_total;
  int indx, indy, indz, i, j, k;
  unsigned extz = 0, extz2, extz3, exth = 0;
  unsigned long nvox, nvox2, nvox3, count, indv, extxy2;
  double origx2, origy2, origz2, origx3, origy3, origz3, width, dval, dsqu;
  int x0, y0, z0, x1, y1, z1;
  double gx, gy, gz, a, b, c;
  double minx, miny, minz, maxx, maxy, maxz;
  double minfact, maxfact, currfact, kmsd, reso;
  int menu_mode, corr_mode = 0, margin;
  double cutoff, bvalue, cvalue;
  double kampl, sigmap, varmap, varp;
  double sig1, sig2, sig3, sig4, sig5, rh1, rh2, rh3, rh4, rh5, rc1, rc2, rc3, rc4, rc5;
  unsigned num1, num2;
  int num_nondens = 0;
  int num_hydrogen = 0;
  int num_water = 0;
  int num_codebook = 0;
  int water_mode, codebook_mode, mass_mode, bfact_mode, hydrogen_mode, nondens_mode;
  PDB *pdb1, *pdb2;

  if (argc != 3) {
    fprintf(stderr, "pdb2vol> Usage: pdb2vol inputfile (PDB format) outputfile (density map) \n");
    exit(1);
  }

  read_pdb(argv[1], &num1, &pdb1);

  /* scan PDB and check content */

  for (i = 0; i < num1; ++i) {
    num_nondens += check_nondens(pdb1, i);
    num_hydrogen += check_hydrogen(pdb1, i);
    num_codebook += check_codebook(pdb1, i);
    num_water += check_water(pdb1, i);
  }

  printf("pdb2vol> Found %d hydrogens, %d water atoms, %d codebook vectors, %d density atoms\n", num_hydrogen, num_water, num_codebook, num1 - num_nondens);

  /* logical assignment of actions based on counted numbers */
  /* care must be taken that all oprions are covered */

  hydrogen_mode = 0; /* ignore hydrogens */
  if (num_hydrogen > 0) printf("pdb2vol> Hydrogens will be ignored. \n");

  if (num_nondens < num1) { /* some density atoms found */

    nondens_mode = 0; /* non-density atoms ignored */
    if (num_nondens > 0) printf("pdb2vol> Non-density atoms will be ignored.\n");
    mass_mode = 1; /* mass_weighting on */
    printf("pdb2vol> Mass-weighting on.\n");
    bfact_mode = 0; /* b-factor thresholding off */
    printf("pdb2vol> B-factor thresholding off.\n");
    codebook_mode = 0; /* codebook vectors ignore - redundant, all nondens atoms will be ignored anyway */
    water_mode = 0; /* water atoms ignored - redundant, all nondens atoms will be ignored anyway */

  } else { /* only nondens atoms */

    nondens_mode = 1;
    if (num_water > 0) {
      printf("pdb2vol> %d water atoms found in file %s \n", num_water, argv[1]);
      mass_mode = 1; /* if water found automated mass_weighting */
      printf("pdb2vol> Mass-weighting on.\n");
      printf("pdb2vol> Do you want to exclude the water atoms?\n");
      printf("pdb2vol> \n");
      printf("pdb2vol>      1: No\n");
      printf("pdb2vol>      2: Yes\n");
      printf("pdb2vol> ");
      water_mode = 2 - readln_int(); /* 0 = water ignored */
      if (!water_mode) printf("pdb2vol> Water atoms will be ignored.\n");
    } else {
      water_mode = 0; /* water atoms ignored - redundant, no water found anyway */
      printf("pdb2vol> Do you want to mass-weight the atoms ?\n");
      printf("pdb2vol> \n");
      printf("pdb2vol>      1: No\n");
      printf("pdb2vol>      2: Yes\n");
      printf("pdb2vol> ");
      mass_mode = readln_int() - 1; /* 0 = mass weighting off */
    }
    printf("pdb2vol> Do you want to select atoms based on a B-factor threshold?\n");
    printf("pdb2vol> \n");
    printf("pdb2vol>      1: No\n");
    printf("pdb2vol>      2: Yes\n");
    printf("pdb2vol> ");
    bfact_mode = readln_int() - 1; /* 0 = B-factor thresholding off */
    if (num_codebook > 0) {
      printf("pdb2vol> %d codebook vectors found in file %s \n", num_codebook, argv[1]);
      printf("pdb2vol> Do you want to exclude the codebook vectors?\n");
      printf("pdb2vol> \n");
      printf("pdb2vol>      1: No\n");
      printf("pdb2vol>      2: Yes\n");
      printf("pdb2vol> ");
      codebook_mode = 2 - readln_int(); /* 0 = codebook ignored */
    } else codebook_mode = 0;
  }
  if (water_mode < 0 || water_mode > 1 || mass_mode < 0 || mass_mode > 1 || bfact_mode < 0 || bfact_mode > 1 ||
      codebook_mode < 0 || codebook_mode > 1 || hydrogen_mode < 0 || hydrogen_mode > 1 || nondens_mode < 0 || nondens_mode > 1) {
    error_option(60000, "pdb2vol");
  }


  /* optional B-factor threshold determination */

  cutoff = 1e20;
  if (bfact_mode) {
    minfact = 1e20;
    maxfact = -1e20;
    for (i = 0; i < num1; ++i) {
      if (!water_mode && check_water(pdb1, i)) continue;
      if (!hydrogen_mode && check_hydrogen(pdb1, i)) continue;
      if (!codebook_mode && check_codebook(pdb1, i)) continue;
      if (!nondens_mode && check_nondens(pdb1, i)) continue;
      currfact = pdb1[i].beta;
      if (currfact < minfact) minfact = currfact;
      if (currfact > maxfact) maxfact = currfact;
    }
    minfact = 0.01 * floor(0.5 + 100.0 * minfact);
    maxfact = 0.01 * floor(0.5 + 100.0 * maxfact);

    /* select cutoff and exit loop if < minfact */
    printf("pdb2vol> Range of crystallographic B-factors: %5.2f - %5.2f.\n", minfact, maxfact);
    printf("pdb2vol> Enter B-factor cutoff (only atoms below this value will be included): ");
    cutoff = readln_double();
    if (cutoff < minfact) {
      printf("pdb2vol> No atoms selected. Try again. Bye bye.\n");
      exit(1);
    }
  }

  /* now copy the useful atoms to pdb2 */

  pdb2 = (PDB *) alloc_vect(num1, sizeof(PDB));
  num2 = 0;
  for (i = 0; i < num1; ++i) {
    if (num2 >= num1) {
      error_atom_count(60010, "pdb2vol", num2, num1);
    }
    if (!water_mode && check_water(pdb1, i)) continue;
    if (!hydrogen_mode && check_hydrogen(pdb1, i)) continue;
    if (!codebook_mode && check_codebook(pdb1, i)) continue;
    if (!nondens_mode && check_nondens(pdb1, i)) continue;
    if (bfact_mode && (0.01 * floor(0.5 + 100.0 * pdb1[i].beta)) > (cutoff + 0.001)) continue;
    copy_atoms(pdb1, pdb2, i, num2, 1);
    if (!mass_mode) pdb2[num2].weight = 1;
    ++num2;
  }
  printf("pdb2vol> %d out of %d atoms selected for conversion.\n", num2, num1);
  printf("pdb2vol> \n");
  free_vect_and_zero_ptr(&pdb1);

  /* from now on we work with pdb2 */

  /* measure protein extent and read voxel spacing */

  minx = 1e20;
  miny = 1e20;
  minz = 1e20;
  maxx = -1e20;
  maxy = -1e20;
  maxz = -1e20;
  for (i = 0; i < num2; ++i) {
    if (minx > pdb2[i].x) minx = pdb2[i].x;
    if (maxx < pdb2[i].x) maxx = pdb2[i].x;
    if (miny > pdb2[i].y) miny = pdb2[i].y;
    if (maxy < pdb2[i].y) maxy = pdb2[i].y;
    if (minz > pdb2[i].z) minz = pdb2[i].z;
    if (maxz < pdb2[i].z) maxz = pdb2[i].z;
  }
  printf("pdb2vol> The input structure measures %6.3f x %6.3f x %6.3f Angstrom\n", maxx - minx, maxy - miny, maxz - minz);
  printf("pdb2vol> \n");
  printf("pdb2vol> Please enter the desired voxel spacing for the output map (in Angstrom): ");
  width = readln_double();
  if (width < 0.1) {
    printf("pdb2vol> Voxel spacing set to minimum value: 0.1 Angstrom\n");
    width = 0.1;
  }
  printf("pdb2vol> \n");

  /* bring lattice into register with origin */
  minx = width * floor(minx / width);
  maxx = width * ceil(maxx / width);
  miny = width * floor(miny / width);
  maxy = width * ceil(maxy / width);
  minz = width * floor(minz / width);
  maxz = width * ceil(maxz / width);

  /* allocate protein density map */
  margin = 2;
  g_extx2 = ceil((maxx - minx) / width) + 2 * margin + 1;
  g_exty2 = ceil((maxy - miny) / width) + 2 * margin + 1;
  extz2 = ceil((maxz - minz) / width) + 2 * margin + 1;
  origx2 = minx - margin * width;
  origy2 = miny - margin * width;
  origz2 = minz - margin * width;
  nvox2 = g_extx2 * g_exty2 * extz2;
  phi2 = (double *) alloc_vect(nvox2, sizeof(double));

  if (nondens_mode) { /* standard situation for biological atoms */
    printf("pdb2vol> Kernel width. Please enter (in Angstrom):\n");
    printf("pdb2vol>      (as pos. value) kernel half-max radius or \n");
    printf("pdb2vol>      (as neg. value) target resolution (2 sigma)\n");
    printf("pdb2vol> Now enter (signed) value: ");
    reso = readln_double();
  } else reso = width / 4.0; /* some fraction of spacing for dens atoms, so only interpolation is done below */

  if (reso < 0.0) {
    sig1 = reso / (-2.0);
    sig2 = reso / (-2.0);
    sig3 = reso / (-2.0);
    sig4 = reso / (-2.0);
    sig5 = reso / (-2.0);
    rh1 = sig1 * sqrt(log(2.0)) / sqrt(1.5);
    rh2 = sig2 / (exp((1.0 / 1.0) * log(2.0)) * sqrt(3.0 * (3.0 + 1.0) / (5.0 * (5.0 + 1.0))));
    rh3 = sig3 / (exp((1.0 / 1.5) * log(2.0)) * sqrt(3.0 * (3.0 + 1.5) / (5.0 * (5.0 + 1.5))));
    rh4 = sig4 / (exp((1.0 / 2.0) * log(2.0)) * sqrt(3.0 * (3.0 + 2.0) / (5.0 * (5.0 + 2.0))));
    rh5 = sig5 / (exp((1.0 / 60.0) * log(2.0)) * sqrt(3.0 * (3.0 + 60.0) / (5.0 * (5.0 + 60.0))));
  } else {
    rh1 = reso;
    rh2 = reso;
    rh3 = reso;
    rh4 = reso;
    rh5 = reso;
    sig1 = rh1 * sqrt(1.5) / sqrt(log(2.0));
    sig2 = rh2 * (exp((1.0 / 1.0) * log(2.0)) * sqrt(3.0 * (3.0 + 1.0) / (5.0 * (5.0 + 1.0))));
    sig3 = rh3 * (exp((1.0 / 1.5) * log(2.0)) * sqrt(3.0 * (3.0 + 1.5) / (5.0 * (5.0 + 1.5))));
    sig4 = rh4 * (exp((1.0 / 2.0) * log(2.0)) * sqrt(3.0 * (3.0 + 2.0) / (5.0 * (5.0 + 2.0))));
    sig5 = rh5 * (exp((1.0 / 60.0) * log(2.0)) * sqrt(3.0 * (3.0 + 60.0) / (5.0 * (5.0 + 60.0))));
  }
  rc1 = sqrt(3.0) * sig1;
  rc2 = (exp((1.0 / 1.0) * log(2.0))) * rh2;
  rc3 = (exp((1.0 / 1.5) * log(2.0))) * rh3;
  rc4 = (exp((1.0 / 2.0) * log(2.0))) * rh4;
  rc5 = (exp((1.0 / 60.0) * log(2.0))) * rh5;

  if (rh1 / width >= 1.0) {
    printf("pdb2vol> \n");
    printf("pdb2vol> Please select the type of smoothing kernel:\n");
    printf("pdb2vol> \n");
    printf("pdb2vol>      1: Gaussian, exp(-1.5 r^2 / sigma^2)\n");
    printf("pdb2vol>         sigma = %6.3fA, r-half = %6.3fA, r-cut = %6.3fA\n", sig1, rh1, rc1);
    printf("pdb2vol> \n");
    printf("pdb2vol>      2: Triangular, max(0, 1 - 0.5 |r| / r-half) \n");
    printf("pdb2vol>         sigma = %6.3fA, r-half = %6.3fA, r-cut = %6.3fA\n", sig2, rh2, rc2);
    printf("pdb2vol> \n");
    printf("pdb2vol>      3: Semi-Epanechnikov, max(0, 1 - 0.5 |r|^1.5 / r-half^1.5) \n");
    printf("pdb2vol>         sigma = %6.3fA, r-half = %6.3fA, r-cut = %6.3fA\n", sig3, rh3, rc3);
    printf("pdb2vol> \n");
    printf("pdb2vol>      4: Epanechnikov, max(0, 1 - 0.5 r^2 / r-half^2) \n");
    printf("pdb2vol>         sigma = %6.3fA, r-half = %6.3fA, r-cut = %6.3fA\n", sig4, rh4, rc4);
    printf("pdb2vol> \n");
    printf("pdb2vol>      5: Hard Sphere, max(0, 1 - 0.5 r^60 / r-half^60) \n");
    printf("pdb2vol>         sigma = %6.3fA, r-half = %6.3fA, r-cut = %6.3fA\n", sig5, rh5, rc5);
    printf("pdb2vol> ");
    menu_mode = readln_int();
    if (menu_mode < 1 || menu_mode > 5) {
      error_option(60120, "pdb2vol");
    }

    printf("pdb2vol> \n");
    printf("pdb2vol> Do you want to correct for lattice interpolation smoothing effects?\n");
    printf("pdb2vol> \n");
    printf("pdb2vol>      1: Yes (slightly lowers the kernel width to maintain target resolution) \n");
    printf("pdb2vol>      2: No \n");
    printf("pdb2vol> ");
    corr_mode = readln_int();
    if (corr_mode < 1 || corr_mode > 2) {
      error_option(60125, "pdb2vol");
    }

    printf("pdb2vol> \n");
    printf("pdb2vol> Finally, please enter the desired kernel amplitude (scaling factor): ");
    kampl = readln_double();
    printf("pdb2vol> \n");

  } else {
    menu_mode = 6;
    kampl = 1.0;
  }

  /* interpolate structure to protein map and keep track of variability */

  printf("pdb2vol> Projecting atoms to cubic lattice by trilinear interpolation... \n");
  for (count = 0; count < nvox2; count++) *(phi2 + count) = 0.0;
  varp = 0.0;
  mass_total = 0.0;
  for (i = 0; i < num2; ++i) {
    /* compute position within grid*/
    gx = margin + (pdb2[i].x - minx) / width;
    gy = margin + (pdb2[i].y - miny) / width;
    gz = margin + (pdb2[i].z - minz) / width;
    x0 = floor(gx);
    y0 = floor(gy);
    z0 = floor(gz);
    x1 = x0 + 1;
    y1 = y0 + 1;
    z1 = z0 + 1;
    /* interpolate */
    a = x1 - gx;
    b = y1 - gy;
    c = z1 - gz;

    *(phi2 + g2idz(z0, y0, x0)) += pdb2[i].weight * a * b * c;
    varp += pdb2[i].weight * a * b * c * ((1 - a) * (1 - a) + (1 - b) * (1 - b) + (1 - c) * (1 - c));
    *(phi2 + g2idz(z1, y0, x0)) += pdb2[i].weight * a * b * (1 - c);
    varp += pdb2[i].weight * a * b * (1 - c) * ((1 - a) * (1 - a) + (1 - b) * (1 - b) + c * c);
    *(phi2 + g2idz(z0, y1, x0)) += pdb2[i].weight * a * (1 - b) * c;
    varp += pdb2[i].weight * a * (1 - b) * c * ((1 - a) * (1 - a) + b * b + (1 - c) * (1 - c));
    *(phi2 + g2idz(z0, y0, x1)) += pdb2[i].weight * (1 - a) * b * c;
    varp += pdb2[i].weight * (1 - a) * b * c * (a * a + (1 - b) * (1 - b) + (1 - c) * (1 - c));
    *(phi2 + g2idz(z1, y1, x0)) += pdb2[i].weight * a * (1 - b) * (1 - c);
    varp += pdb2[i].weight * a * (1 - b) * (1 - c) * ((1 - a) * (1 - a) + b * b + c * c);
    *(phi2 + g2idz(z0, y1, x1)) += pdb2[i].weight * (1 - a) * (1 - b) * c;
    varp += pdb2[i].weight * (1 - a) * (1 - b) * c * (a * a + b * b + (1 - c) * (1 - c));
    *(phi2 + g2idz(z1, y0, x1)) += pdb2[i].weight * (1 - a) * b * (1 - c);
    varp += pdb2[i].weight * (1 - a) * b * (1 - c) * (a * a + (1 - b) * (1 - b) + c * c);
    *(phi2 + g2idz(z1, y1, x1)) += pdb2[i].weight * (1 - a) * (1 - b) * (1 - c);
    varp += pdb2[i].weight * (1 - a) * (1 - b) * (1 - c) * (a * a + b * b + c * c);

    mass_total += pdb2[i].weight;
  }
  varp /= mass_total;
  printf("pdb2vol> ... done. Lattice smoothing (sigma = atom rmsd): %6.3f Angstrom\n", width * sqrt(varp));

  /* compute lattice noise corrected kernel maps */

  switch (menu_mode) {
    case 1:
      printf("pdb2vol> \n");
      if (corr_mode == 2)printf("pdb2vol> Computing Gaussian kernel (no lattice correction) ... \n");
      else printf("pdb2vol> Computing Gaussian kernel (correcting sigma for lattice smoothing)... \n");
      kmsd = sig1 * sig1 / (width * width);

      reso = 2.0 * sig1;
      if (corr_mode == 2) reso = 2.0 * sqrt(sig1 * sig1 + varp * width * width); /* variances are additive for uncorrelated samples */
      varmap = kmsd;
      if (corr_mode == 1) varmap -= varp; /* variances are additive for uncorrelated samples */

      if (varmap < 0.0) {
        error_lattice_smoothing(60130, "pdb2vol");
      }
      sigmap = sqrt(varmap / 3.0); /* sigma-1D */
      exth = (int) ceil(3 * sigmap); /* truncate at 3 sigma-1D == sqrt(3) sigma-3D */
      g_extx = 2 * exth + 1;
      g_exty = g_extx;
      extz = g_extx;
      nvox = g_extx * g_exty * extz;
      phi = (double *) alloc_vect(nvox, sizeof(double));

      /* write Gaussian within 3 sigma-1D to map */
      bvalue = -1.0 / (2.0 * sigmap * sigmap);
      cvalue = 9.0 * sigmap * sigmap;
      for (count = 0; count < nvox; count++) *(phi + count) = 0.0;
      for (indz = 0; indz < extz; indz++)
        for (indy = 0; indy < g_exty; indy++)
          for (indx = 0; indx < g_extx; indx++) {
            dsqu = (indx - exth) * (indx - exth) + (indy - exth) * (indy - exth) + (indz - exth) * (indz - exth);
            if (dsqu < cvalue) *(phi + g1idz(indz, indy, indx)) = kampl * exp(dsqu * bvalue);
          }
      printf("pdb2vol> ... done. Kernel map extent %d x %d x %d voxels \n", g_extx, g_exty, extz);
      printf("pdb2vol> \n");
      break;
    case 2:
      printf("pdb2vol> \n");
      if (corr_mode == 2)printf("pdb2vol> Computing triangular kernel (no lattice correction) ... \n");
      else printf("pdb2vol> Computing triangular kernel (correcting sigma for lattice smoothing)... \n");
      kmsd = sig2 * sig2 / (width * width);

      reso = 2.0 * sig2;
      if (corr_mode == 2) reso = 2.0 * sqrt(sig2 * sig2 + varp * width * width); /* variances are additive for uncorrelated samples */
      varmap = kmsd;
      if (corr_mode == 1) varmap -= varp; /* variances are additive for uncorrelated samples */

      if (varmap < 0.0) {
        error_lattice_smoothing(60210, "pdb2vol");
      }
      exth = (int) ceil(rc2 / width);
      g_extx = 2 * exth + 1;
      g_exty = g_extx;
      extz = g_extx;
      nvox = g_extx * g_exty * extz;
      phi = (double *) alloc_vect(nvox, sizeof(double));

      /* write kernel to map */
      bvalue = 0.5 * exp(-1.0 * log(rh2)) * exp(1.0 * log(width));
      for (count = 0; count < nvox; count++) *(phi + count) = 0.0;
      for (indz = 0; indz < extz; indz++)
        for (indy = 0; indy < g_exty; indy++)
          for (indx = 0; indx < g_extx; indx++) {
            dsqu = (indx - exth) * (indx - exth) + (indy - exth) * (indy - exth) + (indz - exth) * (indz - exth);
            dsqu = exp((1.0 / 2.0) * log(dsqu));
            dval = kampl * (1.0 - dsqu * bvalue);
            if (dval < 0.0) dval = 0.0;
            *(phi + g1idz(indz, indy, indx)) = dval;
          }
      printf("pdb2vol> ... done. Kernel map extent %d x %d x %d voxels \n", g_extx, g_exty, extz);
      printf("pdb2vol> \n");
      break;
    case 3:
      printf("pdb2vol> \n");
      if (corr_mode == 2)printf("pdb2vol> Computing semi-Epanechnikov kernel (no lattice correction) ... \n");
      else printf("pdb2vol> Computing semi-Epanechnikov kernel (correcting sigma for lattice smoothing)... \n");
      kmsd = sig3 * sig3 / (width * width);

      reso = 2.0 * sig3;
      if (corr_mode == 2) reso = 2.0 * sqrt(sig3 * sig3 + varp * width * width); /* variances are additive for uncorrelated samples */
      varmap = kmsd;
      if (corr_mode == 1) varmap -= varp; /* variances are additive for uncorrelated samples */

      if (varmap < 0.0) {
        error_lattice_smoothing(60230, "pdb2vol");
      }
      exth = (int) ceil(rc3 / width);
      g_extx = 2 * exth + 1;
      g_exty = g_extx;
      extz = g_extx;
      nvox = g_extx * g_exty * extz;
      phi = (double *) alloc_vect(nvox, sizeof(double));

      /* write kernel to map */
      bvalue = 0.5 * exp(-1.5 * log(rh3)) * exp(1.5 * log(width));
      for (count = 0; count < nvox; count++) *(phi + count) = 0.0;
      for (indz = 0; indz < extz; indz++)
        for (indy = 0; indy < g_exty; indy++)
          for (indx = 0; indx < g_extx; indx++) {
            dsqu = (indx - exth) * (indx - exth) + (indy - exth) * (indy - exth) + (indz - exth) * (indz - exth);
            dsqu = exp((1.5 / 2.0) * log(dsqu));
            dval = kampl * (1.0 - dsqu * bvalue);
            if (dval < 0.0) dval = 0.0;
            *(phi + g1idz(indz, indy, indx)) = dval;
          }
      printf("pdb2vol> ... done. Kernel map extent %d x %d x %d voxels \n", g_extx, g_exty, extz);
      printf("pdb2vol> \n");
      break;
    case 4:
      printf("pdb2vol> \n");
      if (corr_mode == 2)printf("pdb2vol> Computing Epanechnikov kernel (no lattice correction) ... \n");
      else printf("pdb2vol> Computing Epanechnikov kernel (correcting sigma for lattice smoothing)... \n");
      kmsd = sig4 * sig4 / (width * width);

      reso = 2.0 * sig4;
      if (corr_mode == 2) reso = 2.0 * sqrt(sig4 * sig4 + varp * width * width); /* variances are additive for uncorrelated samples */
      varmap = kmsd;
      if (corr_mode == 1) varmap -= varp; /* variances are additive for uncorrelated samples */

      if (varmap < 0.0) {
        error_lattice_smoothing(60230, "pdb2vol");
      }
      exth = (int) ceil(rc4 / width);
      g_extx = 2 * exth + 1;
      g_exty = g_extx;
      extz = g_extx;
      nvox = g_extx * g_exty * extz;
      phi = (double *) alloc_vect(nvox, sizeof(double));

      /* write kernel to map */
      bvalue = 0.5 * exp(-2.0 * log(rh4)) * exp(2.0 * log(width));
      for (count = 0; count < nvox; count++) *(phi + count) = 0.0;
      for (indz = 0; indz < extz; indz++)
        for (indy = 0; indy < g_exty; indy++)
          for (indx = 0; indx < g_extx; indx++) {
            dsqu = (indx - exth) * (indx - exth) + (indy - exth) * (indy - exth) + (indz - exth) * (indz - exth);
            dsqu = exp((2.0 / 2.0) * log(dsqu));
            dval = kampl * (1.0 - dsqu * bvalue);
            if (dval < 0.0) dval = 0.0;
            *(phi + g1idz(indz, indy, indx)) = dval;
          }
      printf("pdb2vol> ... done. Kernel map extent %d x %d x %d voxels \n", g_extx, g_exty, extz);
      printf("pdb2vol> \n");
      break;
    case 5:
      printf("pdb2vol> \n");
      if (corr_mode == 2)printf("pdb2vol> Computing Hard Sphere kernel (no lattice correction) ... \n");
      else printf("pdb2vol> Computing Hard Sphere kernel (correcting sigma for lattice smoothing)... \n");
      kmsd = sig5 * sig5 / (width * width);
      reso = 2.0 * sig5;
      if (corr_mode == 2) reso = 2.0 * sqrt(sig5 * sig5 + varp * width * width); /* variances are additive for uncorrelated samples */
      varmap = kmsd;
      if (corr_mode == 1) varmap -= varp; /* variances are additive for uncorrelated samples */

      if (varmap < 0.0) {
        error_lattice_smoothing(60240, "pdb2vol");
      }
      exth = (int) ceil(rc5 / width);
      g_extx = 2 * exth + 1;
      g_exty = g_extx;
      extz = g_extx;
      nvox = g_extx * g_exty * extz;
      phi = (double *) alloc_vect(nvox, sizeof(double));

      /* write kernel to map */
      bvalue = 0.5 * exp(-60.0 * log(rh5)) * exp(60.0 * log(width));
      for (count = 0; count < nvox; count++) *(phi + count) = 0.0;
      for (indz = 0; indz < extz; indz++)
        for (indy = 0; indy < g_exty; indy++)
          for (indx = 0; indx < g_extx; indx++) {
            dsqu = (indx - exth) * (indx - exth) + (indy - exth) * (indy - exth) + (indz - exth) * (indz - exth);
            dsqu = exp((60.0 / 2.0) * log(dsqu));
            dval = kampl * (1.0 - dsqu * bvalue);
            if (dval < 0.0) dval = 0.0;
            *(phi + g1idz(indz, indy, indx)) = dval;
          }
      printf("pdb2vol> ... done. Kernel map extent %d x %d x %d voxels \n", g_extx, g_exty, extz);
      printf("pdb2vol> \n");
      break;
  }

  /* convolve and write output */

  switch (menu_mode) {
    case 1:
    case 2:
    case 3:
    case 4:
    case 5:
      printf("pdb2vol> Convolving lattice with kernel... \n");
      /* allocate output density map */
      g_extx3 = ceil((maxx - minx) / width) + 2 * exth + 2 * margin + 1;
      g_exty3 = ceil((maxy - miny) / width) + 2 * exth + 2 * margin + 1;
      extz3 = ceil((maxz - minz) / width) + 2 * exth + 2 * margin + 1;
      origx3 = minx - (exth + margin) * width;
      origy3 = miny - (exth + margin) * width;
      origz3 = minz - (exth + margin) * width;
      nvox3 = g_extx3 * g_exty3 * extz3;
      phi3 = (double *) alloc_vect(nvox3, sizeof(double));
      extxy2 = g_extx2 * g_exty2;
      for (count = 0; count < nvox2; count++) {
        dval = *(phi2 + count);
        if (dval != 0.0) {
          indv = count;
          k = indv / extxy2;
          indv -= k * extxy2;
          j = indv / g_extx2;
          i = indv - j * g_extx2;
          for (indz = 0; indz < extz; indz++) for (indy = 0; indy < g_exty; indy++) for (indx = 0; indx < g_extx; indx++)
                *(phi3 + g3idz(k + indz, j + indy, i + indx)) += *(phi + g1idz(indz, indy, indx)) * dval;
        }
      }
      printf("pdb2vol> ... done. Spatial resolution (2 sigma) of output map: %6.3fA \n", reso);
      if (corr_mode == 2)printf("pdb2vol> (slightly larger than target resolution due to uncorrected lattice smoothing)\n");
      printf("pdb2vol> \n");
      write_vol(argv[2], width, origx3, origy3, origz3, g_extx3, g_exty3, extz3, phi3);
      break;
    case 6:
      write_vol(argv[2], width, origx2, origy2, origz2, g_extx2, g_exty2, extz2, phi2);
      printf("pdb2vol> \n");
      if (nondens_mode) printf("pdb2vol> Warning: kernel width should be large relative to lattice voxel spacing! \n");
      if (nondens_mode) printf("pdb2vol> If kernel smoothing desired, increase kernel width or decrease lattice spacing. \n");
      printf("pdb2vol> Map projected to lattice only, and written to file %s.\n", argv[2]);
      printf("pdb2vol> Effective spatial resolution (2 sigma) of lattice projection: %6.3fA\n", width * sqrt(varp) * 2.0);

      break;
  }
  return 0;
}

/* checks if i-th atom in pdbx is a water atom */
static int check_water(PDB *pdbx, int i)
{
  if ((strcmp(pdbx[i].res, "TIP3") == 0) || (strcmp(pdbx[i].res, "HOH") == 0) || (strcmp(pdbx[i].res, "H2O") == 0)) return 1;
  else return 0;
}

/* checks if i-th atom in pdbx is NOT a vol2pdb density atom */
static int check_nondens(PDB *pdbx, int i)
{
  if (strcmp(pdbx[i].type, "DE") != 0 || strcmp(pdbx[i].loc, "NS") != 0) return 1;
  else return 0;
}

/* checks if i-th atom in pdbx is a codebook vector */
static int check_codebook(PDB *pdbx, int i)
{
  if ((strcmp(pdbx[i].type, "QV") == 0 && strcmp(pdbx[i].loc, "OL") == 0) || (strcmp(pdbx[i].type, "QP") == 0 && strcmp(pdbx[i].loc, "DB") == 0)) return 1;
  else return 0;
}

/* checks if i-th atom in pdbx is a hydrogen atom */
/* note: false hydrogen positives for Hg, Hf, Ho (these are very rare in proteins) */
static int check_hydrogen(PDB *pdbx, int i)
{
  if (pdbx[i].type[0] == 'H' || (pdbx[i].type[0] == ' ' && pdbx[i].type[1] == 'H')) return 1;
  else return 0;
}


static unsigned long g1idz(int k, int j, int i)
{
  return g_extx * g_exty * k + g_extx * j + i;
}

static unsigned long g2idz(int k, int j, int i)
{
  return g_extx2 * g_exty2 * k + g_extx2 * j + i;
}

static unsigned long g3idz(int k, int j, int i)
{
  return g_extx3 * g_exty3 * k + g_extx3 * j + i;
}




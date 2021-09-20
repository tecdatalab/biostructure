/*********************************************************************
*                          Q U A N P D B                             *
**********************************************************************
* Program is part of the Situs package (c) Willy Wriggers, 1998-2019 *
* URL: situs.biomachina.org                                          *
**********************************************************************
*                                                                    *
* Vector quantization (coarse-graining) of atomic structures.        *
*                                                                    *
**********************************************************************
* See legal statement for terms of distribution                      *
*********************************************************************/

#include "situs.h"
#include "lib_std.h"
#include "lib_pio.h"
#include "lib_rnd.h"
#include "lib_pwk.h"
#include "lib_err.h"

#define TMAX 100000     /* # of neural gas iteration steps */
#define NCYCLE 8        /* # of stat. independent clustering cycles */
#define NNMIN 2         /* minimum possible # of codebook vectors, hard limit */
#define NNREC 50        /* maximum recommended # of codebook vectors, soft limit triggers only warning (former NNMAX) */
#define SEED 256        /* seed for random number generator */
#define FLENGTH 1000    /* input file name length */

/* function declarations */

static int check_water(PDB *, int);
static int check_hydrogen(PDB *, int);
static int check_codebook(PDB *, int);
static int check_nondens(PDB *, int);
static int compar(int *, int *);
void alloc_codebook_storgage(int);
void free_codebook_storage();

/* type and global variable declarations */

typedef double Rseq3[3];
static Rseq3 vglob;
static Rseq3 *wglob = 0; /* Rseq3[nn] */
static Rseq3 *average = 0; /* Rseq3[nn] */
static Rseq3 *vecdata = 0; /* Rseq3[nndata] */
static int *tk = 0, *order = 0, *listcount = 0; /* int[nn] */
static int **list; /* int[nn][nndata] */
static unsigned char **cconn; /* unsigned char[nn][nn] */
static double *effrad = 0; /* double[nn] */
static double *mcl = 0, *rmsd = 0; /* double[nn] */
static double *dist; /* double[nn] */

/* the main-program */

int main(int argc, char *argv[])

{
  PDB *pdb1, *pdb2;
  double cutoff;
  int nn, nndata;
  double epsilon, lambda, randoff;
  double ei, ef, li, lf, wdiff;
  int i, j, k, cycle, count = 0, vormode;
  int selatom;
  double mind, currdist, minfact, maxfact, currfact;
  int currindex = 0;
  double xcum, ycum, zcum, rmscum;
  double xdev, ydev, zdev;
  double cumx, cumy, cumz, cum;
  int *windex;
  double currm, currnm, currdiff, rgyr, dist0;
  int skip, psfmode, conncount, opt, nextopt;
  int ch1, ch2, done;
  char psf_file [FLENGTH];
  char pdb_file [FLENGTH];
  Rseq3 com;
  FILE *fout;
  unsigned num1, num2;
  int num_nondens = 0;
  int num_hydrogen = 0;
  int num_water = 0;
  int num_codebook = 0;
  int water_mode, codebook_mode, mass_mode, bfact_mode, hydrogen_mode, nondens_mode;


  if (argc != 3) {
    fprintf(stderr, "quanpdb> Usage: quanpdb inputfile (PDB) outputfile (vectors)\n");
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

  printf("quanpdb> Found %d hydrogens, %d water atoms, %d codebook vectors, %d density atoms\n", num_hydrogen, num_water, num_codebook, num1 - num_nondens);

  /* logical assignment of actions based on counted numbers */
  /* care must be taken that all oprions are covered */

  hydrogen_mode = 0; /* ignore hydrogens */
  if (num_hydrogen > 0) printf("quanpdb> Hydrogens will be ignored. \n");

  if (num_nondens < num1) { /* some density atoms found */

    /* unfortunately the mass-weighting scheme requires chemical atoms and
       does not work with the occupancy encoded masses of density atoms, sorry */

    printf("quanpdb> Vector Quantization of density atoms not possible, use quanvol instead. Bye bye.\n");
    exit(1);

  } else { /* only nondens atoms */

    nondens_mode = 1;
    if (num_water > 0) {
      printf("quanpdb> %d water atoms found in file %s \n", num_water, argv[1]);
      mass_mode = 1; /* if water found automated mass_weighting */
      printf("quanpdb> Mass-weighting on.\n");
      printf("quanpdb> Do you want to exclude the water atoms?\n");
      printf("quanpdb> \n");
      printf("quanpdb>      1: No\n");
      printf("quanpdb>      2: Yes\n");
      printf("quanpdb> ");
      water_mode = 2 - readln_int(); /* 0 = water ignored */
      if (!water_mode) printf("quanpdb> Water atoms will be ignored.\n");
    } else {
      water_mode = 0; /* water atoms ignored - redundant, no water found anyway */
      printf("quanpdb> Do you want to mass-weight the atoms ?\n");
      printf("quanpdb> \n");
      printf("quanpdb>      1: No\n");
      printf("quanpdb>      2: Yes\n");
      printf("quanpdb> ");
      mass_mode = readln_int() - 1; /* 0 = mass weighting off */
    }
    printf("quanpdb> Do you want to select atoms based on a B-factor threshold?\n");
    printf("quanpdb> \n");
    printf("quanpdb>      1: No\n");
    printf("quanpdb>      2: Yes\n");
    printf("quanpdb> ");
    bfact_mode = readln_int() - 1; /* 0 = B-factor thresholding off */
    if (num_codebook > 0) {
      printf("quanpdb> %d codebook vectors found in file %s \n", num_codebook, argv[1]);
      printf("quanpdb> Do you want to exclude the codebook vectors?\n");
      printf("quanpdb> \n");
      printf("quanpdb>      1: No\n");
      printf("quanpdb>      2: Yes\n");
      printf("quanpdb> ");
      codebook_mode = 2 - readln_int(); /* 0 = codebook ignored */
    } else codebook_mode = 0;
  }
  if (water_mode < 0 || water_mode > 1 || mass_mode < 0 || mass_mode > 1 || bfact_mode < 0 || bfact_mode > 1 ||
      codebook_mode < 0 || codebook_mode > 1 || hydrogen_mode < 0 || hydrogen_mode > 1 || nondens_mode < 0 || nondens_mode > 1) {
    error_option(60000, "quanpdb");
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
    printf("quanpdb> Range of crystallographic B-factors: %5.2f - %5.2f.\n", minfact, maxfact);
    printf("quanpdb> Enter B-factor cutoff (only atoms below this value will be included): ");
    cutoff = readln_double();
    if (cutoff < minfact) {
      printf("quanpdb> No atoms selected. Try again. Bye bye.\n");
      exit(1);
    }
  }

  /* now copy the useful atoms to pdb2 */

  pdb2 = (PDB *) alloc_vect(1.3 * num1, sizeof(PDB));
  num2 = 0;
  for (i = 0; i < num1; ++i) {
    if (num2 >= 1.3 * num1) {
      error_atom_count(60010, "quanpdb", num2, num1);
    }
    if (!water_mode && check_water(pdb1, i)) continue;
    if (!hydrogen_mode && check_hydrogen(pdb1, i)) continue;
    if (!codebook_mode && check_codebook(pdb1, i)) continue;
    if (!nondens_mode && check_nondens(pdb1, i)) continue;
    if (bfact_mode && (0.01 * floor(0.5 + 100.0 * pdb1[i].beta)) > (cutoff + 0.001)) continue;
    if (!mass_mode) {
      copy_atoms(pdb1, pdb2, i, num2, 1);
      ++num2;
    } else {
      copy_atoms(pdb1, pdb2, i, num2, (int)floor(pdb1[i].weight / 14.0 + 0.5));
      num2 += (int)floor(pdb1[i].weight / 14.0 + 0.5);
    }
  }
  printf("quanpdb> %d equally weighted inputs out of originally %d atoms selected for conversion.\n", num2, num1);
  printf("quanpdb> \n");

  if (num2 < NNMIN) {
    fprintf(stderr, "quanpdb> Error: Internal number of atoms too small [e.c. 31040]\n");
    exit(31040);
  }

  /* we don't free the smaller pdb1, it is still needed needed below for Voronoi cell output */	
  /* although from now on we work mostly with the larger pdb2 */


  /* compute sphericity */
  cumx = 0;
  cumy = 0;
  cumz = 0;
  for (i = 0; i < num2; ++i) {
    cumx += pdb2[i].x;
    cumy += pdb2[i].y;
    cumz += pdb2[i].z;
  }
  cumx /= (double)num2;
  cumy /= (double)num2;
  cumz /= (double)num2;
  cum = 0;
  for (i = 0; i < num2; ++i) {
    cum += (pdb2[i].x - cumx) * (pdb2[i].x - cumx) + (pdb2[i].y - cumy) * (pdb2[i].y - cumy) + (pdb2[i].z - cumz) * (pdb2[i].z - cumz);
  }
  cum /= (double)num2;
  cum = sqrt(cum);
  printf("quanpdb> Sphericity of the atomic structure:%5.2f \n", num2 / (0.60 * cum * cum * cum));


  /* select nn and exit loop if nn == 0 */
  printf("quanpdb> Enter desired number of codebook vectors for data quantization: (0 to exit): ");
  nn = readln_int();
  if (nn == 0) {
    printf("quanpdb> Bye bye.\n");
    exit(1);
  }
  if (nn > NNREC) {
    printf("quanpdb> Note: number of codebook vectors (%d) exceeds suggested maximum of %d!\n", nn, NNREC);
  }
  if (nn < NNMIN) {
    nn = NNMIN;
    printf("quanpdb> Note: minimum number of codebook vectors (defined in situs.h): %d\n", nn);
  }
  alloc_codebook_storage(nn);

  ei = 1.0;
  ef = 0.01;
  li = 2.0 * nn;
  lf = 0.02;

  /* compute averages and fluctuations for selected nn */

  nndata = nn * NCYCLE;

  printf("quanpdb> Computing %d datasets, %d iterations each... \n", NCYCLE, TMAX);

  /* tessellate and store found codebook vectors in vecdata array */

  for (cycle = 0; cycle < NCYCLE; ++cycle) {
    printf("quanpdb> Now producing dataset %d\n", cycle + 1);

    for (i = 0; i < nn; ++i) {
      selatom = floor(num2 * genrand());
      randoff = 1.0 - 2.0 * genrand();
      wglob[i][0] = pdb2[selatom].x + randoff;
      wglob[i][1] = pdb2[selatom].y + randoff;
      wglob[i][2] = pdb2[selatom].z + randoff;
    }

    for (count = 0; count < TMAX; ++count) {
      selatom = floor(num2 * genrand());
      vglob[0] = pdb2[selatom].x;
      vglob[1] = pdb2[selatom].y;
      vglob[2] = pdb2[selatom].z;
      for (i = 0; i < nn; ++i) order[i] = i;
      qsort(order, nn, sizeof(int), (int (*)(const void *, const void *))compar);
      for (i = 0; i < nn; ++i) tk[order[i]] = i;
      epsilon = ei * exp((count / (double)TMAX) * log(ef / ei));
      lambda =  li * exp((count / (double)TMAX) * log(lf / li));
      for (i = 0; i < nn; ++i) for (j = 0; j < 3; ++j)
          wglob[i][j] = wglob[i][j] + epsilon * exp(-tk[i] / lambda) * (vglob[j] - wglob[i][j]);
    }
    for (i = 0; i < nn; ++i) for (j = 0; j < 3; ++j) vecdata[nn * cycle + i][j] = wglob[i][j];
  }

  /* cluster analysis of found codebook vectors */

  for (j = 0; j < 3; ++j) for (i = 0; i < nn; ++i)  wglob[i][j] = vecdata[i][j];

  do {   /* loop until LBG converged */
    /* determine clusters */
    for (i = 0; i < nn; ++i) listcount[i] = 0;
    for (j = 0; j < nndata; ++j) {
      mind = 1e20;
      for (i = 0; i < nn; ++i) {
        currdist = (wglob[i][0] - vecdata[j][0]) * (wglob[i][0] - vecdata[j][0]) +
                   (wglob[i][1] - vecdata[j][1]) * (wglob[i][1] - vecdata[j][1]) +
                   (wglob[i][2] - vecdata[j][2]) * (wglob[i][2] - vecdata[j][2]);
        if (currdist < mind) {
          mind = currdist;
          currindex = i;
        }
      }
      list[currindex][listcount[currindex]] = j;
      ++listcount[currindex];
    }

    /* compute cluster centers */
    for (currindex = 0; currindex < nn; ++currindex) {
      xcum = 0;
      ycum = 0;
      zcum = 0;
      for (i = 0; i < listcount[currindex]; ++i) {
        xcum += vecdata[list[currindex][i]][0];
        ycum += vecdata[list[currindex][i]][1];
        zcum += vecdata[list[currindex][i]][2];
      }
      if (listcount[currindex] > 0) {
        average[currindex][0] = xcum / (double)listcount[currindex];
        average[currindex][1] = ycum / (double)listcount[currindex];
        average[currindex][2] = zcum / (double)listcount[currindex];
      }
    }

    wdiff = 0;
    for (i = 0; i < nn; ++i) for (j = 0; j < 3; ++j) {
        wdiff += fabs(wglob[i][j] - average[i][j]);
        wglob[i][j] = average[i][j];
      }
  } while (wdiff > 1e-20);

  /* variability within each cluster */

  for (currindex = 0; currindex < nn; ++currindex) {
    rmscum = 0;
    for (i = 0; i < listcount[currindex]; ++i) {
      xdev = vecdata[list[currindex][i]][0] - wglob[currindex][0];
      ydev = vecdata[list[currindex][i]][1] - wglob[currindex][1];
      zdev = vecdata[list[currindex][i]][2] - wglob[currindex][2];
      rmscum += xdev * xdev + ydev * ydev + zdev * zdev;
    }
    if (listcount[currindex] > 0) {
      rmscum /= (double)listcount[currindex];
      rmsd[currindex] = sqrt(rmscum);
    } else {  /* rec. field contains one (or zero) */
      rmsd[currindex] = -1;
      printf("quanpdb> Warning: Voronoi cell %d contains no data points.\n", currindex + 1);
      printf("quanpdb> Vector variability in output file will be set to 99.99.\n");
    }
  }


  /* write output to pdb file */

  fout = fopen(argv[2], "w");
  if (fout == NULL) {
    fprintf(stderr, "quanpdb> Error: Can't open file! %s  [e.c. 31100]\n", argv[2]);
    exit(31100);
  }

  /* compute equivalent spherical radius */

  /* compute closest vector for each atom in new array windex*/
  windex = (int *)alloc_vect(num2, sizeof(int));	
  for (i = 0; i < num2; ++i) {
    k = 0;
    wdiff = 1e20;
    for (j = 0; j < nn; ++j) {
      currdiff = (wglob[j][0] - pdb2[i].x) * (wglob[j][0] - pdb2[i].x) +
                 (wglob[j][1] - pdb2[i].y) * (wglob[j][1] - pdb2[i].y) +
                 (wglob[j][2] - pdb2[i].z) * (wglob[j][2] - pdb2[i].z);
      if (wdiff > currdiff) {
        wdiff = currdiff;
        k = j;
      }
    }
    windex[i] = k;
  }

  /* compute effrad[] and mcl[], free windex*/
  for (j = 0; j < nn; ++j) mcl[j] = 0;
  for (j = 0; j < nn; ++j) effrad[j] = 0;
  for (i = 0; i < num2; ++i) {
    effrad[windex[i]] +=
      (wglob[windex[i]][0] - pdb2[i].x) * (wglob[windex[i]][0] - pdb2[i].x) +
      (wglob[windex[i]][1] - pdb2[i].y) * (wglob[windex[i]][1] - pdb2[i].y) +
      (wglob[windex[i]][2] - pdb2[i].z) * (wglob[windex[i]][2] - pdb2[i].z);
    mcl[windex[i]] += 1;
  }
  free_vect_and_zero_ptr(&windex);

  /* check cell size and compute radius */
  for (i = 0; i < nn; ++i) {
    if (mcl[i] == 0) {
      fprintf(stderr, "quanpdb> Error: Voronoi cell %d contains no density [e.c. 31200]\n", i + 1);
      exit(31200);
    } else effrad[i] = sqrt(5.0 * effrad[i] / (3.0 * mcl[i]));
  }

  skip = 0;
  for (i = 0; i < nn; ++i) {
    if (rmsd[i] < 0) {
      skip = 1;
      fprintf(fout, "ATOM  %5d QPDB QPDB%5d %11.*f%8.*f%8.*f %5.2f %5.2f      QPDB\n", i + 1, i + 1, coord_precision(wglob[i][0]), wglob[i][0], coord_precision(wglob[i][1]), wglob[i][1], coord_precision(wglob[i][2]), wglob[i][2], 99.99, effrad[i]);
    } else {
      fprintf(fout, "ATOM  %5d QPDB QPDB%5d %11.*f%8.*f%8.*f %5.2f %5.2f      QPDB\n", i + 1, i + 1, coord_precision(wglob[i][0]), wglob[i][0], coord_precision(wglob[i][1]), wglob[i][1], coord_precision(wglob[i][2]), wglob[i][2], rmsd[i], effrad[i]);
    }
  }
  fclose(fout);

  printf("quanpdb>\n");
  printf("quanpdb> Codebook vectors have been written to file %s\n", argv[2]);
  printf("quanpdb> The PDB B-factor field contains the equivalent spherical radii\n");
  printf("quanpdb> of the corresponding Voronoi cells (in Angstrom).\n");
  printf("quanpdb> Cluster analysis of the %d independent calculations:\n", NCYCLE);
  printf("quanpdb> The PDB occupancy field in %s contains the rms variabilities of the vectors.\n", argv[2]);

  if (skip == 0) {
    rmscum = 0;
    for (i = 0; i < nn; ++i) rmscum += rmsd[i];
    rmscum /= (double)nn;
    printf("quanpdb> Average rms fluctuation of the %d codebook vectors: %6.3f Angstrom\n", nn, rmscum);
  } else printf("quanpdb> Unable to compute average rms fluctuation, cluster analysis did not converge\n");

  /* write statistics of cluster size deviations*/
  k = 0;
  for (i = 0; i < nn; ++i) k += abs(listcount[i] - NCYCLE);
  if (k > 0) {
    printf("quanpdb> \n");
    printf("quanpdb> Warning: There are cluster size deviations:\n");
    printf("quanpdb> Expected cluster size: %d\n", NCYCLE);
    printf("quanpdb> Actual cluster sizes: ");
    for (i = 0; i < (nn - 1); ++i) printf("%d,", listcount[i]);
    printf("%d.\n", listcount[nn - 1]);
  }

  /* compute com and rgyr of the codebook vectors */
  for (j = 0; j < 3; ++j)
    com[j] = 0;
  for (i = 0; i < nn; ++i)
    for (j = 0; j < 3; ++j)
      com[j] += wglob[i][j];
  for (j = 0; j < 3; ++j)
    com[j] /= (double)nn;
  rgyr = 0;
  for (i = 0; i < nn; ++i)
    for (j = 0; j < 3; ++j)
      rgyr += (wglob[i][j] - com[j]) * (wglob[i][j] - com[j]);
  rgyr /= (double)nn;
  rgyr = sqrt(rgyr);
  printf("quanpdb> Radius of gyration of the %d codebook vectors: %6.3f Angstrom\n", nn, rgyr);

  /* compute or write connectivities if desired */

  printf("quanpdb> \n");
  printf("quanpdb> Do you want to learn nearest-neighbor connectivities?\n");
  printf("quanpdb> Choose one of the following options -\n");
  printf("quanpdb>      1: No. \n");
  printf("quanpdb>      2: Learn and save to a PSF file\n");
  printf("quanpdb>      3: Learn and save to a constraints file\n");
  printf("quanpdb>      4: Learn and save to both PSF and constraints files\n");
  printf("quanpdb> ");
  psfmode = readln_int();


  /* compute connectivities, if necessary */
  switch (psfmode) {

    case 2:
    case 3:
    case 4:

      for (j = 0; j < nn; ++j)
        for (i = 0; i < nn; ++i)
          cconn[i][j] = 0;
      for (j = 0; j < num2; ++j) {
        for (i = 0; i < nn; ++i)
          dist[i] = (pdb2[j].x - wglob[i][0]) * (pdb2[j].x - wglob[i][0]) +
                    (pdb2[j].y - wglob[i][1]) * (pdb2[j].y - wglob[i][1]) +
                    (pdb2[j].z - wglob[i][2]) * (pdb2[j].z - wglob[i][2]);
        currm = 1e20;
        currnm = 1e20;
        nextopt = 0;
        opt = 0;
        for (i = 0; i < nn; i++) {
          if (dist[i] <= currm) {
            currnm = currm;
            nextopt = opt;
            currm = dist[i];
            opt = i;
          } else if (dist[i] <= currnm) {
            currnm = dist[i];
            nextopt = i;
          }
        }
        cconn[opt][nextopt] = 1;
        cconn[nextopt][opt] = 1;
      }
      break;
  }

  /* write files, if necessary */
  switch (psfmode) {

    case 1:
    default:

      break;

    case 2:
    case 4:

      /* read filename and open */
      for (done = 0; done == 0;) {
        printf("quanpdb> Enter PSF filename: ");
        if (fgets(psf_file, FLENGTH, stdin) == NULL) {
          fprintf(stderr, "quanpdb> Error: Can't read filename [e.c. 31300]\n");
          exit(31300);
        }
        removespaces(psf_file, FLENGTH);
        fout = fopen(psf_file, "r");
        if (fout != NULL) {
          printf("quanpdb> Warning: File exists: %s \n", psf_file);
          printf("quanpdb> Do you want to overwrite? (yes/no) ");
          ch1 = getchar();
          for (;;)  {
            ch2 = getchar();
            if (ch2 == EOF) {
              fprintf(stderr, "quanpdb> Error: EOF while reading input [e.c. 31310]\n");
              exit(31310);
            }
            if (ch2 == '\n') break;
          }
          if (ch1 == 'y' || ch1 == 'Y') done = 1;
        } else {
          done = 1;
        }
      }

      fout = fopen(psf_file, "w");
      fprintf(fout, "PSF \n");
      fprintf(fout, " \n");
      fprintf(fout, "       1 !NTITLE\n");
      fprintf(fout, " REMARK Vector connectivity psf file, computed with Competitive Hebbian Rule \n");
      fprintf(fout, " \n");
      fprintf(fout, "%8d !NATOM\n", nn);
      for (i = 0; i < nn; ++i) fprintf(fout, "%8d QPDB %-4d QPDB QPDB QPDB   0.000000       0.00000           0\n", i + 1, i + 1);
      fprintf(fout, " \n");
      conncount = 0;
      for (j = 0; j < nn; ++j) for (i = 0; i < j; ++i) if (cconn[i][j] == 1) ++conncount;
      fprintf(fout, "%8d !NBOND: bonds\n", conncount);
      k = 0;
      for (j = 0; j < nn; ++j) for (i = 0; i < j; ++i) if (cconn[i][j] == 1) {
            if (k == 4) {
              k = 0;
              fprintf(fout, " \n");
            }
            ++k;
            fprintf(fout, "%8d%8d", i + 1, j + 1);
          }
      if (k > 0) fprintf(fout, "\n");
      fprintf(fout, "\n");
      fprintf(fout, "       0 !NTHETA: angles\n");
      fprintf(fout, "\n");
      fprintf(fout, "       0 !NPHI: dihedrals\n");
      fprintf(fout, "\n");
      fprintf(fout, "       0 !NIMPHI: impropers\n");
      fprintf(fout, "\n");
      fprintf(fout, "       0 !NDON: donors\n");
      fprintf(fout, "\n");
      fprintf(fout, "       0 !NACC: acceptors\n");
      fprintf(fout, "\n");
      fprintf(fout, "       0 !NNB\n");
      fprintf(fout, "\n");
      fprintf(fout, "       0       0 !NGRP\n");
      fprintf(fout, "\n");
      fclose(fout);
      printf("quanpdb> Connectivity data written to PSF file %s.\n", psf_file);
      if (psfmode == 2) break;

    case 3:

      /* read filename and open */
      for (done = 0; done == 0;) {
        printf("quanpdb> Enter connectivity filename: ");
        if (fgets(psf_file, FLENGTH, stdin) == NULL) {
          fprintf(stderr, "quanpdb> Error: Can't read filename [e.c. 31330]\n");
          exit(31330);
        }
        removespaces(psf_file, FLENGTH);
        fout = fopen(psf_file, "r");
        if (fout != NULL) {
          printf("quanpdb> Warning: File exists: %s \n", psf_file);
          printf("quanpdb> Do you want to overwrite? (yes/no) ");
          ch1 = getchar();
          for (;;)  {
            ch2 = getchar();
            if (ch2 == EOF) {
              fprintf(stderr, "quanpdb> Error: EOF while reading input [e.c. 31340]\n");
              exit(31340);
            }
            if (ch2 == '\n') break;
          }
          if (ch1 == 'y' || ch1 == 'Y') done = 1;
        } else {
          done = 1;
        }
      }

      fout = fopen(psf_file, "w");

      for (j = 0; j < nn; ++j) for (i = 0; i < j; ++i) if (cconn[i][j] == 1) {
            dist0 = (wglob[i][0] - wglob[j][0]) * (wglob[i][0] - wglob[j][0]) +
                    (wglob[i][1] - wglob[j][1]) * (wglob[i][1] - wglob[j][1]) +
                    (wglob[i][2] - wglob[j][2]) * (wglob[i][2] - wglob[j][2]);
            dist0 = sqrt(dist0);
            fprintf(fout, "%d %d %f\n", i + 1, j + 1, dist0);
          }
      fclose(fout);
      printf("quanpdb> Connectivity data / distances for all vectors written to file %s.\n", psf_file);
      break;
  }


  /* compute or write Voronoi cells if desired */

  printf("quanpdb> \n");
  printf("quanpdb> Do you want to save the Voronoi cells?\n");
  printf("quanpdb> Choose one of the following options -\n");
  printf("quanpdb>      1: No. I'm done \n");
  printf("quanpdb>      2: Yes. Save cells to a PDB file\n");
  printf("quanpdb> ");
  vormode = readln_int();

  /* compute Voronoi cells, if necessary */
  switch (vormode) {

    case 2:

      for (j = 0; j < num1; ++j) {

        /* note: do we really want all atoms here or should some be excluded? */
        for (i = 0; i < nn; ++i)
          dist[i] = (pdb1[j].x - wglob[i][0]) * (pdb1[j].x - wglob[i][0]) +
                    (pdb1[j].y - wglob[i][1]) * (pdb1[j].y - wglob[i][1]) +
                    (pdb1[j].z - wglob[i][2]) * (pdb1[j].z - wglob[i][2]);
        currm = 1e20;
        opt = 0;
        for (i = 0; i < nn; i++) {
          if (dist[i] <= currm) {
            currm = dist[i];
            opt = i;
          }
        }
        pdb1[j].beta = opt + 1;
      }
      break;
  }

  /* write file, if necessary */
  switch (vormode) {

    case 1:
    default:

      printf("quanpdb> Bye bye!\n");
      break;

    case 2:

      /* read filename and open */
      for (done = 0; done == 0;) {
        printf("quanpdb> Enter Voronoi cell PDB filename: ");
        if (fgets(pdb_file, FLENGTH, stdin) == NULL) {
          fprintf(stderr, "quanpdb> Error: Can't read filename [e.c. 31400]\n");
          exit(31400);
        }
        removespaces(pdb_file, FLENGTH);
        fout = fopen(pdb_file, "r");
        if (fout != NULL) {
          printf("quanpdb> Warning: File exists: %s \n", pdb_file);
          printf("quanpdb> Do you want to overwrite? (yes/no) ");
          ch1 = getchar();
          for (;;)  {
            ch2 = getchar();
            if (ch2 == EOF) {
              fprintf(stderr, "quanpdb> Error: EOF while reading input [e.c. 31410]\n");
              exit(31410);
            }
            if (ch2 == '\n') break;
          }
          if (ch1 == 'y' || ch1 == 'Y') done = 1;
        } else {
          done = 1;
        }
      }

      write_pdb(pdb_file, num1, pdb1);
      printf("quanpdb> Voronoi cell numbers written to B-factor field in file %s.\n", pdb_file);
      break;
  }
  free_codebook_storage(0);
  return 0;
}

/*====================================================================*/
static int compar(int *e1, int *e2)
/* distance comparison function for qsort */

{
  double ax, ay, az, d1, d2;

  ax = vglob[0] - wglob[*e1][0];
  ay = vglob[1] - wglob[*e1][1];
  az = vglob[2] - wglob[*e1][2];
  d1 = ax * ax + ay * ay + az * az;
  ax = vglob[0] - wglob[*e2][0];
  ay = vglob[1] - wglob[*e2][1];
  az = vglob[2] - wglob[*e2][2];
  d2 = ax * ax + ay * ay + az * az;
  if (d1 < d2) return -1;
  else return 1;
}

/*====================================================================*/
/* checks if i-th atom in pdbx is a water atom */
static int check_water(PDB *pdbx, int i)
{
  if ((strcmp(pdbx[i].res, "TIP3") == 0) || (strcmp(pdbx[i].res, "HOH") == 0) || (strcmp(pdbx[i].res, "H2O") == 0)) return 1;
  else return 0;
}

/*====================================================================*/
/* checks if i-th atom in pdbx is NOT a vol2pdb density atom */
static int check_nondens(PDB *pdbx, int i)
{
  if (strcmp(pdbx[i].type, "DE") != 0 || strcmp(pdbx[i].loc, "NS") != 0) return 1;
  else return 0;
}

/*====================================================================*/
/* checks if i-th atom in pdbx is a codebook vector */
static int check_codebook(PDB *pdbx, int i)
{
  if ((strcmp(pdbx[i].type, "QV") == 0 && strcmp(pdbx[i].loc, "OL") == 0) || (strcmp(pdbx[i].type, "QP") == 0 && strcmp(pdbx[i].loc, "DB") == 0)) return 1;
  else return 0;
}

/*====================================================================*/
/* checks if i-th atom in pdbx is a hydrogen atom */
/* note: false hydrogen positives for Hg, Hf, Ho (these are very rare in proteins) */
static int check_hydrogen(PDB *pdbx, int i)
{
  if (pdbx[i].type[0] == 'H' || (pdbx[i].type[0] == ' ' && pdbx[i].type[1] == 'H')) return 1;
  else return 0;
}

/*====================================================================*/
void alloc_codebook_storage(unsigned int nn)
{
  int nndata;
	
  nndata = nn * NCYCLE;
  
  /* Following are Rseq3[nn] */
  wglob = (Rseq3 *)alloc_vect(nn, sizeof(Rseq3));
  average = (Rseq3 *)alloc_vect(nn, sizeof(Rseq3));

  /* Following are Rseq3[nndata] */
  vecdata = (Rseq3 *)alloc_vect(nndata, sizeof(Rseq3));

  /* Following are int[nn] */
  tk = (int *)alloc_vect(nn, sizeof(int));
  order = (int *)alloc_vect(nn, sizeof(int));
  listcount = (int *)alloc_vect(nn, sizeof(int));

  /* Following are int[nn][nndata] */
  list = alloc_mat(nn, nndata, sizeof(int));

  /* Following are unsigned char[nn][nn] */
  cconn = alloc_mat(nn, nn, sizeof(unsigned char));

  /* Following are double[nn] */
  effrad = alloc_vect(nn, sizeof(double));
  mcl = alloc_vect(nn, sizeof(double));
  rmsd = alloc_vect(nn, sizeof(double));
  dist = alloc_vect(nn, sizeof(double));
}

/*====================================================================*/
void free_codebook_storage()
{
  free_vect_and_zero_ptr(&wglob);
  free_vect_and_zero_ptr(&average);

  free_vect_and_zero_ptr(&vecdata);

  free_vect_and_zero_ptr(&tk);
  free_vect_and_zero_ptr(&order);
  free_vect_and_zero_ptr(&listcount);

  free_mat_and_zero_ptr(&list);
  free_mat_and_zero_ptr(&cconn);
  
  free_vect_and_zero_ptr(&effrad);
  free_vect_and_zero_ptr(&mcl);
  free_vect_and_zero_ptr(&rmsd);
  free_vect_and_zero_ptr(&dist);
}

/*********************************************************************
 *                           P D B 2 S A X                            *
 **********************************************************************
 * Program is part of the Situs package (c) Willy Wriggers, 2005-12   *
 * URL: situs.biomachina.org                                          *
 **********************************************************************
 *                                                                    *
 * Creation of simulated SAXS bead models for SAXS validations        *
 *                                                                    *
 **********************************************************************
 * See legal statement for terms of distribution                      *
 *********************************************************************/

#include "situs.h"
#include "lib_std.h"
#include "lib_pio.h"
#include "lib_err.h"
#include "lib_vec.h"

#define AFRAC 0.1  /* bead dens. threshold: # atoms / 1A rad. bead  */
/* note: the average protein has AFRAC ~ 0.5     */

/* function declarations */
static void CopyAtoms(PDB *, PDB *, int, int, int);
static char *Element(int);
static int IsWater(PDB *, int);
static int IsIgnored(PDB *, int);


int main(int argc, char *argv[])
{
  FILE *fnout;
  double radius, sqrrad, cubrad, adist;
  double *box;
  unsigned long nlatt, ml, nl;
  double shift1, shift2, shift3;
  int i, j, k, l, count = 1, acount;
  double gx, gy, gz, a, b, c;
  double minx, miny, minz, maxx, maxy, maxz;
  double mindist, maxdist, currdist;
  int wmode, watmode, warnflag = 0;
  double cutoff;
  unsigned numA, numS, numU;
  PDB *pdbA, *pdbS;

  if (argc != 4) {
    fprintf(stderr, "pdb2sax> Usage: pdb2sax inputfile (PDB-atoms) outputfile (PDB-beads) radius (A)\n");
    exit(1);
  }

  read_pdb(argv[1], &numA, &pdbA);

  /* check for water molecules */
  j = 0;
  for (i = 0; i < numA; ++i) {
    if (IsWater(pdbA, i) && !IsIgnored(pdbA, i)) ++j;
  }
  if (j > 0) {
    printf("pdb2sax> %d water molecules found in file %s \n", j, argv[1]);
    printf("pdb2sax> Do you want to exclude the water molecules from the smoothing?\n");
    printf("pdb2sax> \n");
    printf("pdb2sax>      1: No\n");
    printf("pdb2sax>      2: Yes\n");
    printf("pdb2sax> ");
    watmode = readln_int();
    wmode = 2; /* if water found automatic mass-weighting */
  } else {
    watmode = 1;
    printf("pdb2sax> Do you want to mass-weight the atoms or apply a B-factor cutoff?\n");
    printf("pdb2sax> \n");
    printf("pdb2sax>      1: No\n");
    printf("pdb2sax>      2: Yes\n");
    printf("pdb2sax> ");
    wmode = readln_int();
  }
  if (watmode < 1 || watmode > 2) {
    error_option(60000, "pdb2sax");
  }

  switch (wmode) {
    case 1:
      read_pdb(argv[1], &numS, &pdbS);
      break;
    case 2:
      /* create the mass-weighted atom dataset */
      pdbS = (PDB *) alloc_vect((numA), sizeof(PDB));
      mindist = 1e20;
      for (i = 0; i < numA; ++i) {
        if ((watmode == 2) && IsWater(pdbA, i)) continue;
        currdist = pdbA[i].beta;
        if (currdist < mindist && !IsIgnored(pdbA, i)) {
          mindist = currdist;
        }
      }
      mindist = 0.01 * floor(0.5 + 100.0 * mindist);
      maxdist = -1e20;
      for (i = 0; i < numA; ++i) {
        if ((watmode == 2) && IsWater(pdbA, i)) continue;
        currdist = pdbA[i].beta;
        if (currdist > maxdist && !IsIgnored(pdbA, i)) {
          maxdist = currdist;
        }
      }
      maxdist = 0.01 * floor(0.5 + 100.0 * maxdist);
      /* select cutoff and exit loop if < mindist */
      printf("pdb2sax> Range of crystallographic B-factors: %5.2f - %5.2f.\n", mindist, maxdist);
      printf("pdb2sax> Enter B-factor cutoff (only atoms below this value will be included): ");
      cutoff = readln_double();
      if (cutoff < mindist) {
        printf("pdb2sax> Bye bye.\n");
        exit(1);
      }
      numS = 0;
      numU = 0;
      for (i = 0; i < numA; ++i) {
        if (numS >= numA) {
          error_atom_count(60010, "pdb2sax", numS, numA);
        }
        if ((watmode == 2) && IsWater(pdbA, i)) continue;
        if ((0.01 * floor(0.5 + 100.0 * pdbA[i].beta)) > (cutoff + 0.001) || IsIgnored(pdbA, i)) continue;
        count = 1;


        count = 1;

        for (j = 2; j < 30; ++j) if (strcmp(pdbA[i].type, Element(j)) == 0) {
            if (j < 3) count = 0;
            else if (j < 10) count = 1;
            else if (j < 16) count = 2;
            else if (j < 22) count = 3;
            else if (j < 28) count = 4;
            else count = 5;
            break;
          }
        if (j == 30) for (j = 2; j < 30; ++j) if (strncmp(pdbA[i].type, Element(j), 1) == 0) {
              if (j < 3) count = 0;
              else if (j < 10) count = 1;
              else if (j < 16) count = 2;
              else if (j < 22) count = 3;
              else if (j < 28) count = 4;
              else count = 5;
              break;
            }
        if (j == 30) {
          count = 1;
          warnflag = 1;
        }

        CopyAtoms(pdbA, pdbS, i, numS, count);
        numS += count;
        ++numU;
      }

      if (warnflag) printf("pdb2sax> Warning: There are unrecognized atom types (PDB field 13 and 14)\n");

      printf("pdb2sax> There are %d non-hydrogen atoms, represented by %d equally weighted input vectors\n", numU, numS);
      printf("pdb2sax> \n");
      break;
    default:
      error_option(60020, "pdb2sax");
  }

  free_vect_and_zero_ptr(&pdbA);

  /* measure protein extent and read radius */
  minx = 1e20;
  miny = 1e20;
  minz = 1e20;
  maxx = -1e20;
  maxy = -1e20;
  maxz = -1e20;
  for (i = 0; i < numS; ++i) {
    if (minx > pdbS[i].x) minx = pdbS[i].x;
    if (maxx < pdbS[i].x) maxx = pdbS[i].x;
    if (miny > pdbS[i].y) miny = pdbS[i].y;
    if (maxy < pdbS[i].y) maxy = pdbS[i].y;
    if (minz > pdbS[i].z) minz = pdbS[i].z;
    if (maxz < pdbS[i].z) maxz = pdbS[i].z;
  }
  printf("pdb2sax> The input structure measures %6.3f x %6.3f x %6.3f Angstrom\n", maxx - minx, maxy - miny, maxz - minz);
  printf("pdb2sax> \n");
  radius = atof(argv[3]);
  if (radius < 0.1) {
    printf("pdb2sax> Radius set to minimum value: 0.1 Angstrom\n");
    radius = 0.1;
  }
  sqrrad = radius * radius;
  cubrad = sqrrad * radius;


  /* center protein */
  for (i = 0; i < numS; ++i) {
    pdbS[i].x -= 0.5 * (maxx + minx);
    pdbS[i].y -= 0.5 * (maxy + miny);
    pdbS[i].z -= 0.5 * (maxz + minz);
  }

  /* compute grid size */
  a = ceil((maxx - minx + 2 * radius) / (2.0 * radius));
  b = ceil((maxy - miny + 2 * radius) / (radius * sqrt(3.0)));
  c = ceil((maxz - minz + 2 * radius) / (radius * sqrt(8.0 / 3.0)));
  if (fmod(a, 2) == 0.0) a++;
  if (fmod(b, 2) == 0.0) b++;
  if (fmod(c, 2) == 0.0) c++;
  nlatt = (unsigned long)((c * (b + 1) + b) * (a + 1) + a);

  /* allocate grid */
  box = (double *) alloc_vect(nlatt * 3, sizeof(double));

  /* initialize beads on grid that satisfy density criterion */

  ml = 0;
  for (k = 0; k <= c; k++) {
    if (fmod(k + 3, 3) == 0.0) {
      shift3 = 0;
      shift2 = 0;
    } else if (fmod(k + 2, 3) == 0.0) {
      shift2 = radius;
      shift3 = radius * sqrt(3.0) / 3.0;
    } else {
      shift2 = 0;
      shift3 = radius * 2.0 * sqrt(3.0) / 3.0;
    }
    for (j = 0; j <= b; j++)  {
      if (fmod(j + 2, 2) == 0.0) shift1 = 0;
      else shift1 = radius;
      for (i = 0; i <= a; i++) {
        gx = -a * radius + 2 * radius * i + shift1 + shift2;
        gy = -b * radius * sqrt(3) / 2.0 + radius * j * sqrt(3) + shift3;
        gz = -c * radius * sqrt(8.0 / 3.0) / 2.0 + radius * k * sqrt(8.0 / 3.0);
        acount = 0;
        for (l = 0; l < numS; ++l) {
          adist = (pdbS[l].x - gx) * (pdbS[l].x - gx) + (pdbS[l].y - gy) * (pdbS[l].y - gy) +
                  (pdbS[l].z - gz) * (pdbS[l].z - gz);
          if (adist < sqrrad) ++acount;
        }
        if (acount > AFRAC * cubrad) { /* save bead */
          *(box + 3 * ml)   = gx;
          *(box + 3 * ml + 1) = gy;
          *(box + 3 * ml + 2) = gz;
          ++ml;
        }
      }
    }
  }


  /* initialize rest with junk */
  for (nl = ml; nl < nlatt; ++nl) {
    *(box + 3 * nl)   = -1e20;
    *(box + 3 * nl + 1) = -1e20;
    *(box + 3 * nl + 2) = -1e20;
  }

  /* open output file and write out beads */
  fnout = fopen(argv[2], "w");
  if (fnout == NULL) {
    fprintf(stderr, "pdb2sax> Error: Can't open file! %s  [e.c. 81300]\n", argv[2]);
    exit(81300);
  }
  for (nl = 0; nl < ml; ++nl) {
    gx = *(box + 3 * nl);
    gy = *(box + 3 * nl + 1);
    gz = *(box + 3 * nl + 2);
    fprintf(fnout, "ATOM  %5ld  C   SPH %5ld %11.*f%8.*f%8.*f %5.2f  0.00\n", nl + 1, nl + 1, coord_precision(gx), gx, coord_precision(gy), gy, coord_precision(gz), gz, radius);
  }
  fclose(fnout);
  printf("pdb2sax> %ld beads written to file %s .\n", ml, argv[2]);


  return 0;
}

/* copies ktot duplicates of atom # pi of pdbP to pdbS starting at # si */
static void CopyAtoms(PDB *pdbP, PDB *pdbS, int ip, int is, int ktot)
{
  int j, k;

  for (k = 0; k < ktot; ++k) {
    pdbS[is + k].x = pdbP[ip].x;
    pdbS[is + k].y = pdbP[ip].y;
    pdbS[is + k].z = pdbP[ip].z;
    pdbS[is + k].weight = pdbP[ip].weight;
    for (j = 0; j < 4; ++j) pdbS[is + k].segid[j] = pdbP[ip].segid[j];
    pdbS[is + k].serial = is + k + 1;
    for (j = 0; j < 7; ++j) pdbS[is + k].recd[j] = pdbP[ip].recd[j];
    for (j = 0; j < 3; ++j) pdbS[is + k].type[j] = pdbP[ip].type[j];
    for (j = 0; j < 3; ++j) pdbS[is + k].loc[j] = pdbP[ip].loc[j];
    for (j = 0; j < 2; ++j) pdbS[is + k].alt[j] = pdbP[ip].alt[j];
    for (j = 0; j < 5; ++j) pdbS[is + k].res[j] = pdbP[ip].res[j];
    for (j = 0; j < 2; ++j) pdbS[is + k].chain[j] = pdbP[ip].chain[j];
    pdbS[is + k].seq = pdbP[ip].seq;
    for (j = 0; j < 2; ++j) pdbS[is + k].icode[j] = pdbP[ip].icode[j];
    pdbS[is + k].occupancy = pdbP[ip].occupancy;
    pdbS[is + k].beta = pdbP[ip].beta;
    pdbS[is + k].footnote = pdbP[ip].footnote;
    for (j = 0; j < 3; ++j) pdbS[is + k].element[j] = pdbP[ip].element[j];
    for (j = 0; j < 3; ++j) pdbS[is + k].charge[j] = pdbP[ip].charge[j];
  }
  return;
}

/* checks if i-th atom in pdbX is a water molecule */
static int IsWater(PDB *pdbX, int i)
{
  if ((strcmp(pdbX[i].res, "TIP3") == 0) || (strcmp(pdbX[i].res, "HOH") == 0) || (strcmp(pdbX[i].res, "H2O") == 0)) return 1;
  else return 0;
}

/* checks if i-th atom in pdbX is a hydrogen or a codebook vector */
/* Note: false hydrogen positives for Hg, Hf, Ho (these are very rare in proteins) */
static int IsIgnored(PDB *pdbX, int i)
{
  if (pdbX[i].type[0] == 'H' || pdbX[i].type[0] == 'h' || (strcmp(pdbX[i].type, "QV") == 0 && strcmp(pdbX[i].loc, "OL") == 0) || (strcmp(pdbX[i].type, "QP") == 0 && strcmp(pdbX[i].loc, "DB") == 0)) return 1;
  else return 0;
}


/* returns first 30 element symbols */
static char *Element(int n)
{
  static char *name[] = {
    "H", "HE", "LI", "BE", "B", "C", "N", "O", "F", "NE",
    "NA", "MG", "AL", "SI", "P", "S", "CL", "AR", "K", "CA",
    "SC", "TI", "V", "CR", "MN", "FE", "CO", "NI", "CU", "ZN"
  };

  if ((n < 0) || (n >= 30)) {
    error_out_of_index(31510, "pdb2sax");
    return 0;
  } else return name[n];
}

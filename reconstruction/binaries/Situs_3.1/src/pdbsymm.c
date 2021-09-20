/*********************************************************************
*                           P D B S Y M M                            *
**********************************************************************
* Program is part of the Situs package URL: situs.biomachina.org     *
* (c) Valerio Mariani and Willy Wriggers, 1998-2005                  *
**********************************************************************
*                                                                    *
* C, D, and H (helical) symmetry builder                             *
*                                                                    *
**********************************************************************
* See legal statement for terms of distribution                      *
*********************************************************************/

#include "situs.h"
#include "lib_std.h"
#include "lib_vio.h"
#include "lib_err.h"
#include "lib_pio.h"

int main(int argc, char *argv[])
{
  int ibeg, iend;
  double theta, zrise;
  int i, j, k, l, unit;
  unsigned numAtoms;
  PDB *pdbP, *pdbH;
  double vcos, vsin;
  char newsegid[5];
  char type;
  double xoff = 0, yoff = 0, zoff = 0;
  double *phi1;
  double origx, origy, origz, width;
  unsigned extx, exty, extz;
  int order;
  double cenx, ceny, cenz, cengx, cengy, cengz;

  if (argc < 3 || argc > 4) {
    fprintf(stderr, "pdbsymm> Usage: pdbsymm inputfile (PDB format) [optional: input density map] outputfile (PDB format) \n");
    exit(1);
  }

  read_pdb(argv[1], &numAtoms, &pdbP);

  printf("pdbsymm> Enter type of symmetry (C,D, or H): ");
  type = readln_char();
  switch (type) {
    case 'h':
    case 'H':
      if (argc == 4) {
        printf("pdbsymm> Reading density map to determine helical axis (assuming center of x,y plane in z-direction).\n");
        read_vol(argv[2], &width, &origx, &origy, &origz, &extx, &exty, &extz, &phi1);
        cengx = (extx - 1.0) / 2.0;
        cengy = (exty - 1.0) / 2.0;
        cenx = floor(cengx + 0.51);
        ceny = floor(cengy + 0.51);
        if (fabs(cengx - cenx) > 0.1) printf("pdbsymm> X-center will be rounded up to nearest voxel %i\n", 1 + (int)cenx);
        if (fabs(cengy - ceny) > 0.1) printf("pdbsymm> Y-center will be rounded up to nearest voxel %i\n", 1 + (int)ceny);
        xoff = origx + width * cenx;
        yoff = origy + width * ceny;
        printf("pdbsymm> Helical axis: x-coord: %f A, y-coord: %f A\n", xoff, yoff);
      }
      printf("pdbsymm> Enter helical rise per monomer (in Angstrom): ");
      zrise = readln_double();
      printf("pdbsymm> Enter angular twist per monomer (in degrees): ");
      theta = readln_double();
      printf("pdbsymm> Enter number of monomers before input monomer: ");
      ibeg = readln_int();
      ibeg *= -1;
      printf("pdbsymm> Enter number of monomers after input monomer: ");
      iend = readln_int();

      if (argc == 3) {
        printf("pdbsymm> Enter x-coord of helical axis in Angstrom: ");
        xoff = readln_double();
        printf("pdbsymm> Enter y-coord of helical axis in Angstrom: ");
        yoff = readln_double();
      }
      pdbH = (PDB *) alloc_vect(numAtoms * (iend - ibeg + 1), sizeof(PDB));

      k = 0;
      for (unit = ibeg; unit <= iend; unit++) {
        for (i = 0; i < numAtoms; i++) {
          /* transform coordinates */
          vcos = cos((double)unit * PI * theta / 180.0);
          vsin = sin((double)unit * PI * theta / 180.0);
          pdbH[i + numAtoms * k].x = xoff + vcos * (pdbP[i].x - xoff) - vsin * (pdbP[i].y - yoff);
          pdbH[i + numAtoms * k].y = yoff + vsin * (pdbP[i].x - xoff) + vcos * (pdbP[i].y - yoff);
          pdbH[i + numAtoms * k].z = pdbP[i].z + (double)unit * zrise;

          /* add number to segid, keeping two letters, if possible */
          sprintf(newsegid, "%d", k + 1);
          j = 0;
          while (newsegid[j] != '\0' && j < 4)  {
            pdbH[i + numAtoms * k].segid[j] = newsegid[j];
            ++j;
          }
          l = 0;
          while (pdbP[i].segid[l] != '\0' && j < 4 && l < 2)  {
            pdbH[i + numAtoms * k].segid[j] = pdbP[i].segid[l];
            ++j;
            ++l;
          }
          pdbH[i + numAtoms * k].segid[j] = '\0';

          /* update atom number */
          pdbH[i + numAtoms * k].serial = i + numAtoms * k;

          /* keep rest as is */
          for (j = 0; j < 7; ++j) pdbH[i + numAtoms * k].recd[j] = pdbP[i].recd[j];
          for (j = 0; j < 3; ++j) pdbH[i + numAtoms * k].type[j] = pdbP[i].type[j];
          for (j = 0; j < 3; ++j) pdbH[i + numAtoms * k].loc[j] = pdbP[i].loc[j];
          for (j = 0; j < 2; ++j) pdbH[i + numAtoms * k].alt[j] = pdbP[i].alt[j];
          for (j = 0; j < 5; ++j) pdbH[i + numAtoms * k].res[j] = pdbP[i].res[j];
          for (j = 0; j < 2; ++j) pdbH[i + numAtoms * k].chain[j] = pdbP[i].chain[j];
          pdbH[i + numAtoms * k].seq = pdbP[i].seq;
          for (j = 0; j < 2; ++j) pdbH[i + numAtoms * k].icode[j] = pdbP[i].icode[j];
          pdbH[i + numAtoms * k].occupancy = pdbP[i].occupancy;
          pdbH[i + numAtoms * k].beta = pdbP[i].beta;
          pdbH[i + numAtoms * k].footnote = pdbP[i].footnote;
          for (j = 0; j < 3; ++j) pdbH[i + numAtoms * k].element[j] = pdbP[i].element[j];
          for (j = 0; j < 3; ++j) pdbH[i + numAtoms * k].charge[j] = pdbP[i].charge[j];
        }
        ++k;
      }
      write_pdb(argv[argc - 1], numAtoms * k, pdbH);
      printf("pdbsymm> %d atoms written to file %s.\n", numAtoms * k, argv[argc - 1]);
      break;
    case 'c':
    case 'C':
      if (argc == 4) {
        printf("pdbsymm> Reading density map to determine C axis (assuming center of x,y plane in z-direction).\n");
        read_vol(argv[2], &width, &origx, &origy, &origz, &extx, &exty, &extz, &phi1);
        cengx = (extx - 1.0) / 2.0;
        cengy = (exty - 1.0) / 2.0;
        cenx = floor(cengx + 0.51);
        ceny = floor(cengy + 0.51);
        if (fabs(cengx - cenx) > 0.1) printf("pdbsymm> X-center will be rounded up to nearest voxel %i\n", 1 + (int)cenx);
        if (fabs(cengy - ceny) > 0.1) printf("pdbsymm> Y-center will be rounded up to nearest voxel %i\n", 1 + (int)ceny);
        xoff = origx + width * cenx;
        yoff = origy + width * ceny;
        printf("pdbsymm> Principal Axis: x-coord: %f A, y-coord: %f A\n", xoff, yoff);
      }
      printf("pdbsymm> Enter C symmetry order: ");
      order = readln_int();
      theta = (double)360.0 / (double)order;

      if (argc == 3) {
        printf("pdbsymm> Enter x-coord of principal axis in Angstrom: ");
        xoff = readln_double();
        printf("pdbsymm> Enter y-coord of principal axis in Angstrom: ");
        yoff = readln_double();
      }

      pdbH = (PDB *) alloc_vect(numAtoms * order, sizeof(PDB));
      k = 0;
      for (unit = 0; unit < order; unit++) {
        for (i = 0; i < numAtoms; i++) {
          /* transform coordinates */
          vcos = cos((double)unit * PI * theta / 180.0);
          vsin = sin((double)unit * PI * theta / 180.0);
          pdbH[i + numAtoms * k].x = xoff + vcos * (pdbP[i].x - xoff) - vsin * (pdbP[i].y - yoff);
          pdbH[i + numAtoms * k].y = yoff + vsin * (pdbP[i].x - xoff) + vcos * (pdbP[i].y - yoff);
          pdbH[i + numAtoms * k].z = pdbP[i].z;


          /* add number to segid, keeping two letters, if possible */
          sprintf(newsegid, "%d", k + 1);
          j = 0;
          while (newsegid[j] != '\0' && j < 4) {
            pdbH[i + numAtoms * k].segid[j] = newsegid[j];
            ++j;
          }
          l = 0;
          while (pdbP[i].segid[l] != '\0' && j < 4 && l < 2) {
            pdbH[i + numAtoms * k].segid[j] = pdbP[i].segid[l];
            ++j;
            ++l;
          }
          pdbH[i + numAtoms * k].segid[j] = '\0';

          /* update atom number */
          pdbH[i + numAtoms * k].serial = i + numAtoms * k;

          /* keep rest as is */
          for (j = 0; j < 7; ++j) pdbH[i + numAtoms * k].recd[j] = pdbP[i].recd[j];
          for (j = 0; j < 3; ++j) pdbH[i + numAtoms * k].type[j] = pdbP[i].type[j];
          for (j = 0; j < 3; ++j) pdbH[i + numAtoms * k].loc[j] = pdbP[i].loc[j];
          for (j = 0; j < 2; ++j) pdbH[i + numAtoms * k].alt[j] = pdbP[i].alt[j];
          for (j = 0; j < 5; ++j) pdbH[i + numAtoms * k].res[j] = pdbP[i].res[j];
          for (j = 0; j < 2; ++j) pdbH[i + numAtoms * k].chain[j] = pdbP[i].chain[j];
          pdbH[i + numAtoms * k].seq = pdbP[i].seq;
          for (j = 0; j < 2; ++j) pdbH[i + numAtoms * k].icode[j] = pdbP[i].icode[j];
          pdbH[i + numAtoms * k].occupancy = pdbP[i].occupancy;
          pdbH[i + numAtoms * k].beta = pdbP[i].beta;
          pdbH[i + numAtoms * k].footnote = pdbP[i].footnote;
          for (j = 0; j < 3; ++j) pdbH[i + numAtoms * k].element[j] = pdbP[i].element[j];
          for (j = 0; j < 3; ++j) pdbH[i + numAtoms * k].charge[j] = pdbP[i].charge[j];
        }
        ++k;
      }
      write_pdb(argv[argc - 1], numAtoms * k, pdbH);
      printf("pdbsymm> %d atoms written to file %s.\n", numAtoms * k, argv[argc - 1]);
      break;
    case 'd':
    case 'D':
      if (argc == 4) {
        printf("pdbsymm> Reading density map to determine principal and secondary axes.\n");
        read_vol(argv[2], &width, &origx, &origy, &origz, &extx, &exty, &extz, &phi1);
        cengx = (extx - 1.0) / 2.0;
        cengy = (exty - 1.0) / 2.0;
        cengz = (extz - 1.0) / 2.0;
        cenx = floor(cengx + 0.51);
        ceny = floor(cengy + 0.51);
        cenz = floor(cengz + 0.51);
        if (fabs(cengx - cenx) > 0.1) printf("pdbsymm> X-center will be rounded up to nearest voxel %i\n", 1 + (int)cenx);
        if (fabs(cengy - ceny) > 0.1) printf("pdbsymm> Y-center will be rounded up to nearest voxel %i\n", 1 + (int)ceny);
        if (fabs(cengz - cenz) > 0.1) printf("pdbsymm> Z-center will be rounded up to nearest voxel %i\n", 1 + (int)cenz);
        xoff = origx + width * cenx;
        yoff = origy + width * ceny;
        zoff = origz + width * cenz;
        printf("pdbsymm> Principal axis in z direction: x-coord: %f A, y-coord: %f A\n", xoff, yoff);
        printf("pdbsymm> Secondary order 2 axis: y-coord: %f A, z-coord: %f A\n", yoff, zoff);
      }

      printf("pdbsymm> Enter symmetry order for principal axis: ");
      order = readln_int();
      theta = (double)360.0 / (double)order;

      if (argc == 3) {
        printf("pdbsymm> Enter x-coord of principal axis in Angstrom: ");
        xoff = readln_double();
        printf("pdbsymm> Enter y-coord of principal and secondary axes in Angstrom: ");
        yoff = readln_double();
        printf("pdbsymm> Enter z-coord for secondary order 2 axis in Angstrom: ");
        zoff = readln_double();
      }

      pdbH = (PDB *) alloc_vect(numAtoms * order * 2, sizeof(PDB));
      k = 0;
      for (unit = 0; unit < (2 * order); unit++) {
        for (i = 0; i < numAtoms; i++) {
          if (unit < order) {
            /* transform coordinates for upper half */
            vcos = cos((double)unit * PI * theta / 180.0);
            vsin = sin((double)unit * PI * theta / 180.0);
            pdbH[i + numAtoms * k].x = xoff + vcos * (pdbP[i].x - xoff) - vsin * (pdbP[i].y - yoff);
            pdbH[i + numAtoms * k].y = yoff + vsin * (pdbP[i].x - xoff) + vcos * (pdbP[i].y - yoff);
            pdbH[i + numAtoms * k].z = pdbP[i].z;
          } else {
            vcos = cos((double)(unit - order) * PI * theta / 180.0);
            vsin = sin((double)(unit - order) * PI * theta / 180.0);
            pdbH[i + numAtoms * k].x = xoff + vcos * (pdbP[i].x - xoff) - vsin * (pdbP[i].y - yoff);
            pdbH[i + numAtoms * k].y = yoff - vsin * (pdbP[i].x - xoff) - vcos * (pdbP[i].y - yoff);
            pdbH[i + numAtoms * k].z = 2.0 * zoff - pdbP[i].z;
          }

          /* add number to segid, keeping two letters, if possible */
          sprintf(newsegid, "%d", k + 1);
          j = 0;
          while (newsegid[j] != '\0' && j < 4) {
            pdbH[i + numAtoms * k].segid[j] = newsegid[j];
            ++j;
          }
          l = 0;
          while (pdbP[i].segid[l] != '\0' && j < 4 && l < 2) {
            pdbH[i + numAtoms * k].segid[j] = pdbP[i].segid[l];
            ++j;
            ++l;
          }
          pdbH[i + numAtoms * k].segid[j] = '\0';

          /* update atom number */
          pdbH[i + numAtoms * k].serial = i + numAtoms * k;

          /* keep rest as is */
          for (j = 0; j < 7; ++j) pdbH[i + numAtoms * k].recd[j] = pdbP[i].recd[j];
          for (j = 0; j < 3; ++j) pdbH[i + numAtoms * k].type[j] = pdbP[i].type[j];
          for (j = 0; j < 3; ++j) pdbH[i + numAtoms * k].loc[j] = pdbP[i].loc[j];
          for (j = 0; j < 2; ++j) pdbH[i + numAtoms * k].alt[j] = pdbP[i].alt[j];
          for (j = 0; j < 5; ++j) pdbH[i + numAtoms * k].res[j] = pdbP[i].res[j];
          for (j = 0; j < 2; ++j) pdbH[i + numAtoms * k].chain[j] = pdbP[i].chain[j];
          pdbH[i + numAtoms * k].seq = pdbP[i].seq;
          for (j = 0; j < 2; ++j) pdbH[i + numAtoms * k].icode[j] = pdbP[i].icode[j];
          pdbH[i + numAtoms * k].occupancy = pdbP[i].occupancy;
          pdbH[i + numAtoms * k].beta = pdbP[i].beta;
          pdbH[i + numAtoms * k].footnote = pdbP[i].footnote;
          for (j = 0; j < 3; ++j) pdbH[i + numAtoms * k].element[j] = pdbP[i].element[j];
          for (j = 0; j < 3; ++j) pdbH[i + numAtoms * k].charge[j] = pdbP[i].charge[j];
        }
        ++k;
      }
      write_pdb(argv[argc - 1], numAtoms * k, pdbH);
      printf("pdbsymm> %d atoms written to file %s.\n", numAtoms * k, argv[argc - 1]);
      break;
    default:
      error_symmetry_option("pdbsymm");
      break;
  }
  return 0;
}



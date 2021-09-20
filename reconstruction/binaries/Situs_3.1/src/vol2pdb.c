/*********************************************************************
*                          V O L 2 P D B                             *
**********************************************************************
* Program is part of the Situs package (c) Willy Wriggers, 1998-2009 *
* URL: situs.biomachina.org                                          *
**********************************************************************
*                                                                    *
* Convert map to PDB file, write densities to occupancy              *
*                                                                    *
**********************************************************************
* See legal statement for terms of distribution                      *
*********************************************************************/

#include "situs.h"
#include "lib_vio.h"
#include "lib_pio.h"
#include "lib_vwk.h"
#include "lib_err.h"

int main(int argc, char *argv[])
{
  double *phi;
  unsigned extx, exty, extz;
  unsigned long nvox, count, extxy, indv, atom_count;
  unsigned indx, indy, indz;
  double origx, origy, origz, width;
  double mindens1, maxdens1, mindens2;
  PDB *pdbH;
  int j;
  char recd[7] = "ATOM  ";
  char de[3] = "DE";
  char ns[3] = "NS";
  char dens[5] = "DENS";
  char space1[2] = " ";
  char space2[3] = "  ";
  double cutdens;

  if (argc < 3 || argc > 4) {
    fprintf(stderr, "vol2pdb> Usage: vol2pdb input-map [optional: density cutoff] output-PDB\n");
    exit(1);
  }

  /* read map and set parameters */
  read_vol(argv[1], &width, &origx, &origy, &origz, &extx, &exty, &extz, &phi);
  nvox = extx * exty * extz;
  extxy = extx * exty;
  print_map_info(phi, nvox);
  maxdens1 = calc_max(phi, nvox);
  mindens1 = calc_min(phi, nvox);
  mindens2 = mindens1;

  /* assert useable max density values */
  if (maxdens1 <= 0.0) {
    fprintf(stderr, "vol2pdb> Error: Maximum map density is negative. Only positive maps are currently supported.\n");
    exit(1);
  }
  if (maxdens1 < 1.0 || maxdens1 > 99.99) {
    printf("vol2pdb> Warning: Density values will be rescaled to [0,99.99] interval to fit PDB occupancy field.\n");
    printf("vol2pdb> Scaling factor: %f .\n", 99.99 / maxdens1);
    normalize(phi, nvox, maxdens1 / 99.99);
    mindens2 = calc_min(phi, nvox);
  }

  /* set cutoff */
  if (argc == 3) {
    cutdens = 0.005;
  } else {
    sscanf(argv[2], "%lf", &cutdens);
    if (cutdens < 0.0) {
      fprintf(stderr, "vol2pdb> Warning: Negative cutoff densities are not supported. Cutoff value set to 0.005.\n");
      cutdens = 0.005;
    } else if (cutdens < 0.005) {
      fprintf(stderr, "vol2pdb> Cutoff value set to 0.005 (densities will be rounded to two decimals). \n");
      cutdens = 0.005;
    }
  }

  /* mindens warning */
  if (mindens2 < cutdens) {
    if (mindens1 == mindens2) printf("vol2pdb> Warning: Densities below cutoff %f will be ignored. \n", cutdens);
    else printf("vol2pdb> Warning: Rescaled densities below cutoff %f will be ignored. \n", cutdens);
  }

  /* allocate atoms above cutoff */
  atom_count = 0;
  for (count = 0; count < nvox; count++) if (*(phi + count) >= cutdens) ++atom_count;
  pdbH = (PDB *) alloc_vect(atom_count, sizeof(PDB));

  /* now assign PDB fields */
  atom_count = 0;
  for (count = 0; count < nvox; count++) {
    if (*(phi + count) >= cutdens) {
      pdbH[atom_count].serial = 1; /* will be renumbered anyway */
      for (j = 0; j < 7; j++) pdbH[atom_count].recd[j] = recd[j];
      for (j = 0; j < 3; j++) pdbH[atom_count].type[j] = de[j];
      for (j = 0; j < 3; j++) pdbH[atom_count].loc[j] = ns[j];
      for (j = 0; j < 2; j++) pdbH[atom_count].alt[j] = space1[j];
      for (j = 0; j < 5; j++) pdbH[atom_count].res[j] = dens[j];
      for (j = 0; j < 2; j++) pdbH[atom_count].chain[j] = space1[j];
      pdbH[atom_count].seq = 1;
      for (j = 0; j < 2; j++) pdbH[atom_count].icode[j] = space1[j];
      indv = count;
      indz = indv / extxy;
      indv -= indz * extxy;
      indy = indv / extx;
      indx = indv - indy * extx;
      pdbH[atom_count].x = indx * width + origx;
      pdbH[atom_count].y = indy * width + origy;
      pdbH[atom_count].z = indz * width + origz;
      pdbH[atom_count].occupancy = *(phi + count);
      pdbH[atom_count].beta = 0.0f;
      pdbH[atom_count].footnote = 0;
      for (j = 0; j < 5; j++) pdbH[atom_count].segid[j] = dens[j];
      for (j = 0; j < 3; j++) pdbH[atom_count].element[j] = space2[j];
      for (j = 0; j < 3; j++) pdbH[atom_count].charge[j] = space2[j];
      pdbH[atom_count].weight = 0.0f;
      ++atom_count;
    }
  }
  write_pdb(argv[argc - 1], atom_count, pdbH);
  printf("vol2pdb> The PDB occupancy field contains the densities above cutoff rounded to two decimals.\n");
  return 0;
}

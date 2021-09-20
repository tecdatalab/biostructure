/*********************************************************************
*                          V O L A V E R                             *
**********************************************************************
* Program is part of the Situs package (c) Willy Wriggers, 2011      *
* URL: situs.biomachina.org                                          *
**********************************************************************
*                                                                    *
* Multiple map averaging tool.                                       *
*                                                                    *
**********************************************************************
* See legal statement for terms of distribution                      *
*********************************************************************/

#include "situs.h"
#include "lib_vio.h"
#include "lib_vwk.h"
#include "lib_err.h"

int main(int argc, char *argv[])
{
  double *phi1, *phi2, *phi3, *pphi1;
  unsigned extx1, exty1, extz1;
  unsigned extx2, exty2, extz2;
  unsigned long nvox1, count;
  double origx1, origy1, origz1, width1;
  double origx2, origy2, origz2, width2;
  int i;
  int minx, miny, minz, maxx, maxy, maxz;

  if (argc < 4) {
    fprintf(stderr, "volaver> Usage: volaver infile1 infile2 [infile...] outfile \n");
    fprintf(stderr, "volaver>   infile1:     reference density map defining voxel spacing \n");
    fprintf(stderr, "volaver>   infile2:     at least one more input density map \n");
    fprintf(stderr, "volaver>   [infile...]: optional additional input density maps\n");
    fprintf(stderr, "volaver>   outfile:     averaged density map \n");
    exit(1);
  }

  /* read initial phi1 and set preliminary grid geometry */
  read_vol(argv[1], &width1, &origx1, &origy1, &origz1, &extx1, &exty1, &extz1, &phi1);
  nvox1 = extx1 * exty1 * extz1;

  for (i = 2; i < (argc - 1); ++i) {

    /* read next map, phi2 */
    read_vol(argv[i], &width2, &origx2, &origy2, &origz2, &extx2, &exty2, &extz2, &phi2);

    /* update phi1 grid geometry based on phi2 extent (see also get_pad_parameters in voledit.c) */
    if (origx2 < origx1)
      minx = 0 - (int) ceil((origx1 - origx2) / width1);
    else minx = 0;
    if (origy2 < origy1)
      miny = 0 - (int) ceil((origy1 - origy2) / width1);
    else miny = 0;
    if (origz2 < origz1)
      minz = 0 - (int) ceil((origz1 - origz2) / width1);
    else minz = 0;
    maxx = extx1 - 1 + (int) ceil(((origx2 + (extx2 - 1) * width2) - (origx1 + (extx1 - 1) * width1)) / width1);
    if (maxx < extx1) maxx = extx1 - 1;
    maxy = exty1 - 1 + (int) ceil(((origy2 + (exty2 - 1) * width2) - (origy1 + (exty1 - 1) * width1)) / width1);
    if (maxy < exty1) maxy = exty1 - 1;
    maxz = extz1 - 1 + (int) ceil(((origz2 + (extz2 - 1) * width2) - (origz1 + (extz1 - 1) * width1)) / width1);
    if (maxz < extz1) maxz = extz1 - 1;
    project_map_lattice(&pphi1, maxx - minx + 1, maxy - miny + 1, maxz - minz + 1,
                        origx1 + minx * width1, origy1 + miny * width1, origz1 + minz * width1, width1, width1, width1,
                        phi1, extx1, exty1, extz1, origx1, origy1, origz1, width1, width1, width1);
    free_vect_and_zero_ptr(&phi1);
    phi1 = pphi1;
    extx1 = maxx - minx + 1;
    exty1 = maxy - miny + 1;
    extz1 = maxz - minz + 1;
    origx1 += minx * width1;
    origy1 += miny * width1;
    origz1 += minz * width1;
    nvox1 = extx1 * exty1 * extz1;

    /* now project phi2 to updated phi1 and add */
    project_map_lattice(&phi3, extx1, exty1, extz1, origx1, origy1, origz1,
                        width1, width1, width1, phi2, extx2, exty2, extz2,
                        origx2, origy2, origz2, width2, width2, width2);
    free_vect_and_zero_ptr(&phi2);
    for (count = 0; count < nvox1; count++) *(phi1 + count) += *(phi3 + count);
    free_vect_and_zero_ptr(&phi3);
  }

  /* normalize and write output */
  for (count = 0; count < nvox1; count++)
    *(phi1 + count) /= (1.0 * (argc - 2));
  write_vol(argv[argc - 1], width1, origx1, origy1, origz1, extx1, exty1, extz1, phi1);
  return 0;
}


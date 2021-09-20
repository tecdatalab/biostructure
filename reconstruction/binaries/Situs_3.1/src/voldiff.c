/*********************************************************************
*                          V O L D I F F                             *
**********************************************************************
* Program is part of the Situs package (c) Willy Wriggers, 1998-2009 *
* URL: situs.biomachina.org                                          *
**********************************************************************
*                                                                    *
* Difference-mapping tool.                                           *
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
  double *phi1, *phi2, *phi3;
  unsigned extx1, exty1, extz1;
  unsigned extx2, exty2, extz2;
  unsigned long nvox1, count;
  double origx1, origy1, origz1, width1;
  double origx2, origy2, origz2, width2;

  if (argc != 4) {
    fprintf(stderr, "voldiff> Usage: voldiff density-map-1 density-map-2 output-density-map\n");
    exit(1);
  }

  read_vol(argv[1], &width1, &origx1, &origy1, &origz1, &extx1, &exty1, &extz1, &phi1);
  nvox1 = extx1 * exty1 * extz1;
  read_vol(argv[2], &width2, &origx2, &origy2, &origz2, &extx2, &exty2, &extz2, &phi2);

  project_map_lattice(&phi3, extx1, exty1, extz1, origx1, origy1, origz1,
                      width1, width1, width1, phi2, extx2, exty2, extz2,
                      origx2, origy2, origz2, width2, width2, width2);
  free_vect_and_zero_ptr(&phi2);
  for (count = 0; count < nvox1; count++)
    *(phi1 + count) -= *(phi3 + count);
  write_vol(argv[3], width1, origx1, origy1, origz1, extx1, exty1, extz1, phi1);
  return 0;
}


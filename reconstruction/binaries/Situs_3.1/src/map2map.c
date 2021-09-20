/*********************************************************************
*                         M A P 2 M A P                              *
**********************************************************************
* Program is part of the Situs package (c) Willy Wriggers, 1998-2011 *
* URL: situs.biomachina.org                                          *
**********************************************************************
*                                                                    *
* Map file format conversion.                                        *
*                                                                    *
**********************************************************************
* See legal statement for terms of distribution                      *
*********************************************************************/

#include "situs.h"
#include "lib_vio.h"
#include "lib_vwk.h"
#include "lib_vec.h"
#include "lib_std.h"
#include "lib_err.h"

int main(int argc, char *argv[])
{
  double porigx, porigy, porigz;
  double *phi;
  double *pphi;
  unsigned pextx, pexty, pextz;
  unsigned long nvox;
  int menumode, ordermode = 7, swapmode, cubic = 1, orom = 1;
  double pwidth, widthx, widthy, widthz;
  double alpha, beta, gamma;
  int nxstart = 0, nystart = 0, nzstart = 0;
  int ncstart = 0, nrstart = 0, nsstart = 0;
  unsigned nc, nr, ns;
  unsigned nx, ny, nz;
  int currext;
  double xorigin, yorigin, zorigin;
  char ac = 'X', ar = 'Y', as = 'Z';

  if (argc != 3) {
    printf("map2map> Usage: map2map inputfile outputfile\n");
    printf("map2map> Situs, MRC/CCP4,and SPIDER input file formats should be automatically recognized.\n");
    printf("map2map> Output file format will be either classic Situs (.sit or.situs suffix) or otherwise MRC/CCP4.\n");
    printf("map2map> See online user guide for details (http://situs.biomachina.org/fguide.html#map2map).\n");
    exit(1);
  }
  /************** process map format selection ******************/

  /* start with some sanity checks and modify menu accordingly */

  if (test_situs_header_and_suffix(argv[1]) == 0 && have_situs_suffix(argv[2])) { /* this makes little sense but we need to capture it */
    printf("map2map> Input and output map seem to have both Situs format. Are you sure? ... \n");
    menumode = 0;

  } else if (test_situs_header_and_suffix(argv[1]) == 0) { /* strict test, excludes fortuitous ASCII, can assume Situs input map */
    printf("map2map> Input map in classic Situs (ASCII text) format detected.\n");
    printf("map2map> Choose one of the following options (or 0 for classic menu):\n");
    printf("map2map> \n");
    printf("map2map>      1: Convert to MRC / CCP4 binary (auto)*\n");     /* can assume MRC/CCP4, above excludes Situs output option */
    printf("map2map>      2: Convert to MRC / CCP4 binary (manual)**\n");   /* can assume MRC/CCP4, above excludes Situs output option */
    printf("map2map>      3: Convert to SPIDER binary* \n");
    printf("map2map>      4: Convert to X-PLOR (ASCII editable text)* \n");
    printf("map2map> \n");
    printf("map2map>      *: automatic fill of header fields \n");
    printf("map2map>     **: manual assignment of header fields \n");
    printf("map2map> \n");
    printf("map2map> Enter selection: ");
    menumode = readln_int();
    switch (menumode) {
      case 0:
        break;
      case 1:
        menumode = 7;
        break;
      case 2:
        menumode = 8;
        break;
      case 3:
        menumode = 9;
        break;
      case 4:
        menumode = 10;
        break;
      default:
        menumode = 99;
    }

  } else if ((test_mrc(argv[1], 0) == 0 || test_mrc(argv[1], 1) == 0) && have_situs_suffix(argv[2])) { /* reasonable to assume MRC/CCP4 input, Situs output */
    printf("map2map> Input map in MRC or CCP4 binary format detected. Desired output format seems to be Situs based on file extension.\n");
    printf("map2map> Choose one of the following options (or 0 for classic menu):\n");
    printf("map2map> \n");
    printf("map2map>      1: Convert to classic Situs (auto)*\n");
    printf("map2map>      2: Convert to classic Situs (manual)**\n");
    printf("map2map> \n");
    printf("map2map>      *: automatic fill of header fields \n");
    printf("map2map>     **: manual assignment of header fields \n");
    printf("map2map> \n");
    printf("map2map> Enter selection: ");
    menumode = readln_int();
    switch (menumode) {
      case 0:
        break;
      case 1:
        menumode = 2;
        break;
      case 2:
        menumode = 3;
        break;
      default:
        menumode = 99;
    }

  } else if (test_mrc(argv[1], 0) == 0 || test_mrc(argv[1], 1) == 0) { /* reasonable to assume MRC/CCP4 input, non-Situs output */
    printf("map2map> Input map in MRC or CCP4 binary format detected.\n");
    printf("map2map> Choose one of the following options (or 0 for classic menu):\n");
    printf("map2map> \n");
    printf("map2map>      1: Wash self (auto)*\n");  /* can assume MRC/CCP4, above excludes Situs output option */
    printf("map2map>      2: Wash self or edit header (manual)**\n");  /* can assume MRC/CCP4, above excludes Situs output option */
    printf("map2map>      3: Convert to SPIDER binary* \n");
    printf("map2map>      4: Convert to X-PLOR (ASCII editable text)* \n");
    printf("map2map> \n");
    printf("map2map>      *: automatic fill of header fields \n");
    printf("map2map>     **: manual assignment of header fields \n");
    printf("map2map> \n");
    printf("map2map> Enter selection: ");
    menumode = readln_int();
    switch (menumode) {
      case 0:
        break;
      case 1:
        menumode = 7;
        break;
      case 2:
        menumode = 8;
        break;
      case 3:
        menumode = 9;
        break;
      case 4:
        menumode = 10;
        break;
      default:
        menumode = 99;
    }

  } else if (test_spider(argv[1], 0) == 0 || test_spider(argv[1], 1) == 0) { /* reasonable after above to assume SPIDER input map */
    printf("map2map> Input map in SPIDER binary format detected. The map will be manually converted to a Situs-\n");
    printf("map2map> compatible format based on the output file extension (Situs if .sit or .situs, else MRC/CCP4).\n");
    menumode = 4;

  } else { /* not enough info, other selection */
    printf("map2map> Input map format not yet recognized.\n");
    printf("map2map> Select one of the following input formats (or 0 for classic menu):\n");
    printf("map2map> \n");
    printf("map2map>      1: ASCII (editable text) file, sequential list of map densities** \n");
    printf("map2map>      2: X-PLOR map (ASCII editable text)* \n");
    printf("map2map>      3: Generic 32-bit binary (unknown map or header parameters) \n");
    printf("map2map> \n");
    printf("map2map>      *: automatic fill of header fields \n");
    printf("map2map>     **: manual assignment of header fields \n");
    printf("map2map> \n");
    printf("map2map> Enter selection: ");
    menumode = readln_int();
    switch (menumode) {
      case 0:
        break;
      case 1:
        menumode = 1;
        break;
      case 2:
        menumode = 5;
        break;
      case 3:
        menumode = 6;
        break;
      default:
        menumode = 99;
    }
  }

  if (menumode == 0) {
    printf("map2map> Entering classic map2map menu.\n");
    printf("map2map> Select one of the following options:\n");
    printf("map2map> \n");
    printf("map2map> Convert selected INPUT formats to Situs or MRC/CCP4*: \n");
    printf("map2map> \n");
    printf("map2map>      1: ASCII (editable text) file, sequential list of map densities** \n");
    printf("map2map>      2: MRC or CCP4 binary (auto)** \n");
    printf("map2map>      3: MRC or CCP4 binary (manual)*** \n");
    printf("map2map>      4: SPIDER binary*** \n");
    printf("map2map>      5: X-PLOR map (ASCII editable text)** \n");
    printf("map2map>      6: Generic 32-bit binary (unknown map or header parameters) \n");
    printf("map2map> \n");
    printf("map2map> Convert Situs or MRC/CCP4**** to one of these OUTPUT formats: \n");
    printf("map2map> \n");
    printf("map2map>      7: MRC / CCP4 binary (auto)**\n");
    printf("map2map>      8: MRC / CCP4 binary (manual)***\n");
    printf("map2map>      9: SPIDER binary** \n");
    printf("map2map>     10: X-PLOR (ASCII editable text)** \n");
    printf("map2map> \n");
    printf("map2map>      *: output format based on file extension (Situs if .sit or .situs, else MRC/CCP4)\n");
    printf("map2map>     **: automatic fill of header fields\n");
    printf("map2map>    ***: manual assignment of header fields\n");
    printf("map2map>   ****: input format will be automatically detected, either Situs or MRC/CCP4\n");
    printf("map2map> \n");
    printf("map2map> Enter selection: ");
    menumode = readln_int();
  }



  /************** process input map ******************/

  switch (menumode) {
    case 1: /* free format ASCII */
      /* get order and map size parameters */
      printf("map2map> \n");
      printf("map2map> Data order and axis permutation in file %s.\n", argv[1]);
      printf("map2map> Assign columns (C, fastest), rows (R), and sections (S, slowest) to X, Y, Z:\n");
      printf("map2map> \n");
      printf("map2map>         C  R  S = \n");
      printf("map2map> ------------------\n");
      printf("map2map>      1: X  Y  Z (no permutation)\n");
      printf("map2map>      2: X  Z  Y\n");
      printf("map2map>      3: Y  X  Z\n");
      printf("map2map>      4: Y  Z  X\n");
      printf("map2map>      5: Z  X  Y\n");
      printf("map2map>      6: Z  Y  X\n");
      printf("map2map> Enter selection number [1-6]: ");
      ordermode = readln_int();

      /* assign axis character */
      switch (ordermode) {
        case 1:
          ac = 'X';
          ar = 'Y';
          as = 'Z';
          break;
        case 2:
          ac = 'X';
          ar = 'Z';
          as = 'Y';
          break;
        case 3:
          ac = 'Y';
          ar = 'X';
          as = 'Z';
          break;
        case 4:
          ac = 'Y';
          ar = 'Z';
          as = 'X';
          break;
        case 5:
          ac = 'Z';
          ar = 'X';
          as = 'Y';
          break;
        case 6:
          ac = 'Z';
          ar = 'Y';
          as = 'X';
          break;
        default:
          error_option(70209, "map2map");
      }

      printf("map2map> Enter number of columns (%c fields): ", ac);
      currext = readln_int();
      if (currext < 1) {
        error_number_columns(70010, "map2map");
      } else nc = currext;
      printf("map2map> Enter number of rows (%c fields): ", ar);
      currext = readln_int();
      if (currext < 1) {
        error_number_rows(70020, "map2map");
      } else nr = currext;
      printf("map2map> Enter number of sections (%c fields): ", as);
      currext = readln_int();
      if (currext < 1) {
        error_number_sections(70030, "map2map");
      } else ns = currext;

      /* read map and permute */
      nvox = nc * nr * ns;
      read_ascii(argv[1], nvox, &phi);
      permute_map(ordermode, nc, nr, ns, &pextx, &pexty, &pextz, ncstart, nrstart, nsstart, &nxstart, &nystart, &nzstart, phi, &pphi);

      /* get lattice spacing and origin */
      printf("map2map> Enter desired cubic grid (lattice) spacing in Angstrom: ");
      pwidth = readln_double();
      if (pwidth <= 0) {
        error_number_spacing(70040, "map2map");
      }
      printf("map2map> \n");
      printf("map2map> Enter X-origin (coord of first voxel) in Angstrom (0 is recommended, a different value shifts the resulting map accordingly): ");
      porigx = readln_double();
      printf("map2map> \n");
      printf("map2map> Enter Y-origin (coord of first voxel) in Angstrom (0 is recommended, a different value shifts the resulting map accordingly): ");
      porigy = readln_double();
      printf("map2map> \n");
      printf("map2map> Enter Z-origin (coord of first voxel) in Angstrom (0 is recommended, a different value shifts the resulting map accordingly): ");
      porigz = readln_double();
      printf("map2map> \n");
      break;

    case 2: /* MRC or CCP4 format with automatically filled cubic grid parameters */

      read_mrc(argv[1], &orom, &cubic, &ordermode, &nc, &nr, &ns, &ncstart, &nrstart, &nsstart, &widthx, &widthy, &widthz, &xorigin, &yorigin, &zorigin, &alpha, &beta, &gamma, &phi);
      permute_map(ordermode, nc, nr, ns, &nx, &ny, &nz, ncstart, nrstart, nsstart, &nxstart, &nystart, &nzstart, phi, &pphi);
      permute_print_info(ordermode, nc, nr, ns, nx, ny, nz, ncstart, nrstart, nsstart, nxstart, nystart, nzstart);
      assert_cubic_map(orom, cubic, alpha, beta, gamma, widthx, widthy, widthz, nx, ny, nz, nxstart, nystart, nzstart, xorigin, yorigin, zorigin, &pwidth, &porigx, &porigy, &porigz, &pextx, &pexty, &pextz, &pphi);
      break;

    case 3: /* MRC or CCP4 format with manual override of cubic grid parameters */

      read_mrc(argv[1], &orom, &cubic, &ordermode, &nc, &nr, &ns, &ncstart, &nrstart, &nsstart, &widthx, &widthy, &widthz, &xorigin, &yorigin, &zorigin, &alpha, &beta, &gamma, &phi);
      nvox = nc * nr * ns;
      permute_map(ordermode, nc, nr, ns, &nx, &ny, &nz, ncstart, nrstart, nsstart, &nxstart, &nystart, &nzstart, phi, &pphi);
      permute_print_info(ordermode, nc, nr, ns, nx, ny, nz, ncstart, nrstart, nsstart, nxstart, nystart, nzstart);
      assert_cubic_map(orom, cubic, alpha, beta, gamma, widthx, widthy, widthz, nx, ny, nz, nxstart, nystart, nzstart, xorigin, yorigin, zorigin, &pwidth, &porigx, &porigy, &porigz, &pextx, &pexty, &pextz, &pphi);

      /* offer to override cubic grid parameters */
      printf("map2map> \n");
      printf("map2map> The currently assigned voxel spacing is %f Angstrom\n", pwidth);
      printf("map2map> To keep the field, enter this value again (different value rescales the map): ");
      pwidth = readln_double();
      if (pwidth <= 0) {
        error_number_spacing(70140, "map2map");
      }
      printf("map2map> \n");
      printf("map2map> The currently assigned X-origin (coord of first voxel) is %f Angstrom\n", porigx);
      printf("map2map> To keep the field, enter this value again (different value shifts the map): ");
      porigx = readln_double();
      printf("map2map> \n");
      printf("map2map> The currently assigned Y-origin (coord of first voxel) is %f Angstrom\n", porigy);
      printf("map2map> To keep the field, enter this value again (different value shifts the map): ");
      porigy = readln_double();
      printf("map2map> \n");
      printf("map2map> The currently assigned Z-origin (coord of first voxel) is %f Angstrom\n", porigz);
      printf("map2map> To keep the field, enter this value again (different value shifts the map): ");
      porigz = readln_double();
      printf("map2map> \n");
      printf("map2map> The currently assigned NX (# X fields) value is %d \n", pextx);
      printf("map2map> Enter the same or a new value: ");
      pextx = readln_int();
      printf("map2map> The currently assigned NY (# Y fields) value is %d \n", pexty);
      printf("map2map> Enter the same or a new value: ");
      pexty = readln_int();
      printf("map2map> The currently assigned NZ (# Z fields) value is %d \n", pextz);
      printf("map2map> Enter the same or a new value: ");
      pextz = readln_int();
      if (pextx * pexty * pextz != nvox) {
        fprintf(stderr, "map2map> Sorry, NX * NY * NZ must match the total number of voxels, %ld,  in the map.\n", nvox);
        fprintf(stderr, "map2map> Please try again, bye bye.\n");
        exit(1);
      }
      break;

    case 4: /* SPIDER format with manual assignment of cubic grid parameters */
      read_spider(argv[1], &pextx, &pexty, &pextz, &pphi);
      printf("map2map> \n");
      printf("map2map> SPIDER maps don't typically contain scale and origin information relative to a PDB coordinate system.\n");
      printf("map2map> Enter desired grid (lattice) spacing in Angstrom: ");
      pwidth = readln_double();
      if (pwidth <= 0) {
        error_number_spacing(70141, "map2map");
      }
      printf("map2map> \n");
      printf("map2map> Enter X-origin (coord of first voxel) in Angstrom (0 is recommended, a different value shifts the resulting map accordingly): ");
      porigx = readln_double();
      printf("map2map> \n");
      printf("map2map> Enter Y-origin (coord of first voxel) in Angstrom (0 is recommended, a different value shifts the resulting map accordingly): ");
      porigy = readln_double();
      printf("map2map> \n");
      printf("map2map> Enter Z-origin (coord of first voxel) in Angstrom (0 is recommended, a different value shifts the resulting map accordingly): ");
      porigz = readln_double();
      printf("map2map> \n");
      break;

    case 5: /* X-PLOR ASCII format */
      read_xplor(argv[1], &orom, &cubic, &widthx, &widthy, &widthz, &alpha, &beta, &gamma, &nxstart, &nystart, &nzstart, &nx, &ny, &nz, &pphi);
      assert_cubic_map(orom, cubic, alpha, beta, gamma, widthx, widthy, widthz, nx, ny, nz, nxstart, nystart, nzstart, 0.0, 0.0, 0.0, &pwidth, &porigx, &porigy, &porigz, &pextx, &pexty, &pextz, &pphi);
      break;

    case 6: /* generic binary */
      printf("map2map> \n");
      printf("map2map> Do you want to swap the byte order (endianism)?\n");
      printf("map2map> \n");
      printf("map2map>      1: No (order is correct for current machine architecture)\n");
      printf("map2map>      2: Yes (foreign machine binary)\n");
      printf("map2map> ");
      swapmode = readln_int() - 1;
      dump_binary_and_exit(argv[1], argv[2], swapmode);
      break;

    case 7:
    case 8:
    case 9:
    case 10: /* Situs or MRC/CCP4 input format determined by test functions */
      ordermode = 1;
      read_vol(argv[1], &pwidth, &porigx, &porigy, &porigz, &pextx, &pexty, &pextz, &pphi);
      break;

    default:
      fprintf(stderr, "map2map> Error: Unable to process option [e.c. 70209]\n");
      exit(70209);
  }

  /************** write output map ******************/

  switch (menumode) {
    case 1:
    case 2:
    case 3:
    case 4:
    case 5: /* Situs or MRC/CCP4 format output determined by file extension */
      write_vol(argv[2], pwidth, porigx, porigy, porigz, pextx, pexty, pextz, pphi);
      break;

    case 6: /* generic binary already dumped, this option should not be reached */
      fprintf(stderr, "map2map> Error: Unable to process option [e.c. 70205]\n");
      exit(70205);
      break;

    case 7: /* MRC/CCP4 format, auto */
      write_mrc(1, argv[2], pwidth, porigx, porigy, porigz, pextx, pexty, pextz, pphi);
      break;

    case 8: /* MRC/CCP4 format, manual */
      write_mrc(0, argv[2], pwidth, porigx, porigy, porigz, pextx, pexty, pextz, pphi);
      break;

    case 9: /* SPIDER format */
      write_spider(argv[2], pwidth, porigx, porigy, porigz, pextx, pexty, pextz, pphi);
      break;

    case 10: /* X-PLOR format */
      write_xplor(argv[2], pwidth, porigx, porigy, porigz, pextx, pexty, pextz, pphi);
      break;
  }
  printf("map2map> \n");
  printf("map2map> All done.\n");
  return 0;
}



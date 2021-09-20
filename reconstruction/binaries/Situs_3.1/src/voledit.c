/*********************************************************************
*                           V O L E D I T                            *
**********************************************************************
* Program is part of the Situs package URL: situs.biomachina.org     *
* (c) Paul Boyle, Jochen Heyd, Pablo Chacon, Willy Wriggers 1998-2011*
**********************************************************************
*                                                                    *
* Slice inspection and general editing tool for volumetric maps      *
*                                                                    *
**********************************************************************
* See legal statement for terms of distribution                      *
*********************************************************************/

#include "situs.h"
#include "lib_vio.h"
#include "lib_std.h"
#include "lib_vwk.h"
#include "lib_vec.h"
#include "lib_err.h"

#define NUM_VERTEX 100  /* number of polygon clipping vertices */
#define DEFAULT_ZERO 0.0 /* default background density */
#define FLENGTH 1000  /* file name length */


typedef struct listentry *LLISTPTR; /* list for non-recursive flood fill */
typedef struct listentry {
  int      x;
  int      y;
  int      z;
  LLISTPTR next;
} LLIST;


static unsigned long pindex(int, int, int, int, int, int, int);
static char check_inside_polygon(int [NUM_VERTEX][2], int , int , int);
static void zero_pixel_in_all_slices(int, int, int, int, int, int, int, 
                                     int, double *pphi);
static void get_crop_parameters(int *, int *, int *, int *, int *, int *, 
                                int, int, int);
static void get_pad_parameters(int *, int *, int *, int *, int *, int *, 
                               int, int, int);
static void update_parameters(int , double **, int *, int *, int *, 
                              unsigned *, unsigned *, unsigned *,
                              double *, double *, double *, int , int , int ,
                              int , int , int, double, int *, double *, int *, 
                              double *, unsigned long *);
static FILE *get_save_file_info(int, char *, char *);
static int neighbor(char *, int, int, int, int, int, int);
static void flood_fill_nonrecursive(double *, char *, double, int, int, int, 
                                    int, int, int);
static void flood_fill_nonrecursive_push(LLIST **tail, int x, int y, int z);

int main(int argc, char *argv[])
{
  int i, j, k, k_shift, k_old, l, done;
  int indx = 0, indy = 0, indz = 0;
  unsigned extx, exty, extz;
  unsigned extx2, exty2, extz2;
  unsigned long count = 0, nvox, ind1, ind2, nngbr, ncore;
  double origx, origy, origz, width;
  double origx2, origy2, origz2, width2;
  double *pphi, *pphi2 = NULL, *ppro = NULL;
  double cutoff = 0, cutoff_old = 0, maxdensity, mindensity;
  char cha = 0, chb = 0, chc = 0;
  char out_file[FLENGTH];
  char *slicefile = "output slice: ";
  char *mapfile = "modified map: ";
  int rangea = 0, rangeb = 0, rangec = 0;
  int minx = 0, maxx = 0, miny = 0, maxy = 0, minz = 0, maxz = 0;
  int slicemode, menumode, clipmode, slitextmode;
  int num_vertx;
  int stride;
  int shrinkflag;
  int xstart, ystart, zstart;
  double pmindist, pcurrdist;
  int xmindist, ymindist, zmindist;
  char *iphi; /* 0=unconnected, 1=contiguous voxel, 2=neighbor of contiguous voxel */
  int polygon[NUM_VERTEX][2];
  char check_polygon;
  FILE *fout;

  if (argc != 2) {
    fprintf(stderr, "voledit> Usage: voledit density-map \n");
    exit(1);
  }

  read_vol(argv[1], &width, &origx, &origy, &origz, &extx, &exty, &extz, &pphi);
  nvox = extx * exty * extz;
  mindensity = calc_min(pphi, nvox);
  maxdensity = calc_max(pphi, nvox);
  printf("voledit> Min. / max. density values: %f / %f\n", mindensity, maxdensity);
  cutoff = 0.5 * (maxdensity + mindensity);

  /* determine slicemode, initial value of k, and projection */

  printf("voledit> Choose one of the following three options -\n");
  printf("voledit>      1: (x,y)-cross section as function of z-value \n");
  printf("voledit>      2: (z,x)-cross section as function of y-value \n");
  printf("voledit>      3: (y,z)-cross section as function of x-value \n");
  printf("voledit> ");
  slicemode = readln_int();

  switch (slicemode) {
    case 1:
      rangea = extx;
      rangeb = exty;
      rangec = extz;
      cha = 'x';
      chb = 'y';
      chc = 'z';
      break;
    case 2:
      rangea = extz;
      rangeb = extx;
      rangec = exty;
      cha = 'z';
      chb = 'x';
      chc = 'y';
      break;
    case 3:
      rangea = exty;
      rangeb = extz;
      rangec = extx;
      cha = 'y';
      chb = 'z';
      chc = 'x';
      break;
    default:
      error_option(51010, "voledit");
  }
  k = rangec / 2;
  stride = ceil(rangea / 100.0);
  if (stride < ceil(rangeb / 100.0)) stride = ceil(rangeb / 100.0);
  if (stride < 1) stride = 1;
  do_vect(&ppro, (unsigned long)(rangea * rangeb));
  for (i = 0; i < rangea; i++) for (j = 0; j < rangeb; j++) {
      ind2 = i + j * rangea;
      for (l = 0; l < rangec; l++) {
        ind1 = pindex(slicemode, i, j, l, extx, exty, extz);
        *(ppro + ind2) += *(pphi + ind1);
      }
      *(ppro + ind2) /= ((double)rangec);
    }

  /* loop through the cross sections / projection */

  for (done = 0; done == 0;) {
    printf("voledit> Box coordinates (1,1,1)=(%f,%f,%f); (%d,%d,%d)=(%f,%f,%f)\n", origx, origy, origz,
           extx, exty, extz, origx + extx * width, origy + exty * width, origz + extz * width);

    k_old = k;
    if (k >= 0) { /* cross section */
      printf("voledit> Cross-section at %c = %d, display density level = %f, display voxel step = %d:\n", chc, k + 1, cutoff, stride);
      for (i = 0; i < rangea; i += stride) printf(" ");
      printf("  %c\n", chb);
      for (i = 0; i < rangea; i += stride) printf(" ");
      printf("  ^\n");
      for (j = (2 * stride * ((rangeb - (2 * stride)) / (2 * stride))); j >= 0; j -= (2 * stride)) {
        for (i = 0; i < rangea; i += stride) {
          ind1 = pindex(slicemode, i, j, k, extx, exty, extz);
          ind2 = pindex(slicemode, i, j + 1, k, extx, exty, extz);
          if ((*(pphi + ind1) > cutoff) && (*(pphi + ind2) > cutoff)) printf("O");
          if ((*(pphi + ind1) <= cutoff) && (*(pphi + ind2) > cutoff)) printf("^");
          if ((*(pphi + ind1) > cutoff) && (*(pphi + ind2) <= cutoff)) printf("u");
          if ((*(pphi + ind1) <= cutoff) && (*(pphi + ind2) <= cutoff)) {
            if ((j / stride) % 4 + (i / stride) % 4) printf(" ");
            else  printf(".");
          }
        }
        if ((j / stride) % 4) printf("\n");
        else printf("- %-3d\n", j + 1);
      }
    } else { /* projection */
      printf("voledit> Projection (%c-average) onto (%c,%c) plane, display density level = %f, display voxel step = %d:\n", chc, cha, chb, cutoff, stride);
      for (i = 0; i < rangea; i += stride) printf(" ");
      printf("  %c\n", chb);
      for (i = 0; i < rangea; i += stride) printf(" ");
      printf("  ^\n");
      for (j = (2 * stride * ((rangeb - (2 * stride)) / (2 * stride))); j >= 0; j -= (2 * stride)) {
        for (i = 0; i < rangea; i += stride) {
          ind1 = i + j * rangea;
          ind2 = i + (j + 1) * rangea;
          if ((*(ppro + ind1) > cutoff) && (*(ppro + ind2) > cutoff)) printf("O");
          if ((*(ppro + ind1) <= cutoff) && (*(ppro + ind2) > cutoff)) printf("^");
          if ((*(ppro + ind1) > cutoff) && (*(ppro + ind2) <= cutoff)) printf("u");
          if ((*(ppro + ind1) <= cutoff) && (*(ppro + ind2) <= cutoff)) {
            if ((j / stride) % 4 + (i / stride) % 4) printf(" ");
            else  printf(".");
          }
        }
        if ((j / stride) % 4) printf("\n");
        else printf("- %-3d\n", j + 1);
      }
    }
    for (i = 0; i < rangea; i += 4 * stride) {
      printf("|   ");
    }
    printf("\n");
    for (i = 0; (i < rangea && i < 1000); i += 4 * stride) if ((i / stride + 1) % 4) printf("%-3d ", i + 1);
    for (i = 1000; i < rangea; i += 4 * stride) if ((i / stride + 1) % 4) printf("%3.3d ", i + 1 - 1000 * ((i + 1) / 1000));
    printf(">%c\n", cha);


    /* ask for new input */

    if (k >= 0) printf("voledit> Browse slices (sections %c=1-%i, projection %c=0): pos/neg int for up/down or 0 for more options: ", chc, rangec, chc);
    else printf("voledit> Enter positive %c index to return to cross-sections or 0 for more options: ", chc);
    k_shift = readln_int();

    if (k_shift == 0) {
      printf("voledit> Choose one of the following options: \n");
      printf("voledit>      1:  Set cross-section %c index or select projection\n", chc);
      printf("voledit>      2:  Set display density level \n");
      printf("voledit>      3:  Set voxel step for display \n");
      printf("voledit>      4:  Cropping (cut away margin)\n");
      printf("voledit>      5:  Zero padding (extra margin)\n");
      printf("voledit>      6:  Interpolation (change voxel spacing)\n");
      printf("voledit>      7:  Polygon clippling\n");
      printf("voledit>      8:  Segmentation (floodfill) of contiguous volume\n");
      printf("voledit>      9:  Thresholding (set low densities to zero)\n");
      printf("voledit>     10:  Binary thresholding (set densities to zero or one)\n");
      printf("voledit>     11:  Save current slice\n");
      printf("voledit>     12:  Save entire map\n");
      printf("voledit>     13:  Quit voledit \n");
      printf("voledit> ");
      menumode = readln_int();

      switch (menumode) {

        case 1: /* select slice */

          printf("voledit> Enter cross-section %c index (from 1-%i, or 0 for projection): ", chc, rangec);
          k = readln_int() - 1;
          break;

        case 2: /* select display level */

          printf("voledit> Enter new display density level: ");
          cutoff = readln_double();
          break;

        case 3: /* select stride */

          printf("voledit> Enter new voxel step (stride) for map display (pos. integer): ");
          stride = readln_int();
          if (stride < 1) stride = 1;

          break;

        case 4: /* crop */

          get_crop_parameters(&minx, &miny, &minz, &maxx, &maxy, &maxz, 
                              extx, exty, extz);
          project_map_lattice(&pphi2, maxx - minx + 1, maxy - miny + 1, 
                              maxz - minz + 1, origx + minx * width, 
                              origy + miny * width, origz + minz * width, 
                              width, width, width, pphi, extx, exty, extz, 
                              origx, origy, origz, width, width, width);
          free_vect_and_zero_ptr(&pphi);
          pphi = pphi2;
          update_parameters(slicemode, &ppro, &rangea, &rangeb, &rangec, 
                            &extx, &exty, &extz, &origx, &origy, &origz,
                            maxx, minx, maxy, miny, maxz, minz, width, &k, 
                            &cutoff, &stride, pphi, &nvox);
          break;

        case 5: /* pad */

          get_pad_parameters(&minx, &miny, &minz, &maxx, &maxy, &maxz, 
                             extx, exty, extz);
          project_map_lattice(&pphi2, maxx - minx + 1, maxy - miny + 1, 
                              maxz - minz + 1, origx + minx * width, 
                              origy + miny * width, origz + minz * width, 
                              width, width, width, pphi, extx, exty, extz, 
                              origx, origy, origz, width, width, width);
          free_vect_and_zero_ptr(&pphi);
          pphi = pphi2;
          update_parameters(slicemode, &ppro, &rangea, &rangeb, &rangec, 
                            &extx, &exty, &extz, &origx, &origy, &origz,
                            maxx, minx, maxy, miny, maxz, minz, width, &k, 
                            &cutoff, &stride, pphi, &nvox);
          break;

        case 6: /* interpolate */

          printf("voledit> Existing voxel spacing: %f Angstrom. Enter desired new value: ", width);
          width2 = readln_double();
          if (width2 <= 0) {
            error_voxel_size(51011, "voledit");
          }
          interpolate_map(&pphi2, &extx2, &exty2, &extz2, &origx2, &origy2, &origz2,
                          width2, width2, width2, pphi, extx, exty, extz, origx,
                          origy, origz, width, width, width);
          free_vect_and_zero_ptr(&pphi);
          pphi = pphi2;
          width = width2;
          minx = 0;
          maxx = extx2 - 1;
          miny = 0;
          maxy = exty2 - 1;
          minz = 0;
          maxz = extz2 - 1;
          update_parameters(slicemode, &ppro, &rangea, &rangeb, &rangec, &extx, &exty, &extz, &origx, &origy, &origz,
                            maxx, minx, maxy, miny, maxz, minz, width, &k, &cutoff, &stride, pphi, &nvox);
          break;

        case 7: /* polygon clipping, see p. 242-5 of Computational Geometry in C (Joseph O'Rourke), 2nd ed., 1999 reprint */

          printf("voledit> Choose one of the following four options -\n");
          printf("voledit>      1: Clip inside, including edges\n");
          printf("voledit>      2: Clip inside, excluding edges\n");
          printf("voledit>      3: Clip outside, including edges\n");
          printf("voledit>      4: Clip outside, excluding edges\n");
          printf("voledit> ");
          clipmode = readln_int();
          printf("voledit> Input the number of vertices: ");
          num_vertx = readln_int();
          if (num_vertx >= NUM_VERTEX) {
            error_number_vertices(51430, "voledit", NUM_VERTEX);
          }
          printf("voledit> For %d vertices enter indices in two columns: %c [1-%d]  %c [1-%d]\n",
                 num_vertx, cha, rangea, chb, rangeb);
          for (j = 0; j < num_vertx; j++)
            scanf("%d %d", &polygon[j][0], &polygon[j][1]);
          for (j = 0; j < num_vertx; j++) {
            if (polygon[j][0] < 1 || polygon[j][0] > rangea)
              printf("voledit> Warning: Vertex outside map, %c = %d must be within [1-%d]\n", cha, polygon[j][0], rangea);
            if (polygon[j][1] < 1 || polygon[j][1] > rangeb)
              printf("voledit> Warning: Vertex outside map, %c = %d must be within [1-%d]\n", chb, polygon[j][1], rangeb);
          }
          for (j = 0; j < num_vertx; j++) { /* reset voxel indices to internal range */
            --polygon[j][0];
            --polygon[j][1];
          }
          for (j = 0; j < rangeb; ++j)
            for (i = 0; i < rangea; ++i) {
              check_polygon = check_inside_polygon(polygon, num_vertx, i, j);
              switch (clipmode) {
                case 1: /* clip inside, including edges */
                  if (check_polygon != 'o')
                    zero_pixel_in_all_slices(i, j, indz, slicemode, extx, 
                                             exty, extz, rangec, pphi);
                  break;
                case 2: /* clip inside, excluding edges */
                  if (check_polygon == 'i')
                    zero_pixel_in_all_slices(i, j, indz, slicemode, extx, 
                                             exty, extz, rangec, pphi);
                  break;
                case 3:  /* clip outside, including edges */
                  if (check_polygon != 'i')
                    zero_pixel_in_all_slices(i, j, indz, slicemode, extx, 
                                             exty, extz, rangec, pphi);
                  break;
                case 4: /* clip outside, excluding edges */
                  if (check_polygon == 'o')
                    zero_pixel_in_all_slices(i, j, indz, slicemode, extx, 
                                             exty, extz, rangec, pphi);
                  break;
              }
            }
          update_parameters(slicemode, &ppro, &rangea, &rangeb, &rangec, 
                            &extx, &exty, &extz, &origx, &origy, &origz,
                            extx - 1, 0, exty - 1, 0, extz - 1, 0, width, &k, 
                            &cutoff, &stride, pphi, &nvox);
          break;

        case 8: /* floodfill */

          printf("voledit> Min. / max. density values: %f / %f\n", mindensity, maxdensity);

          /* input cutoff level */
          printf("voledit> Enter desired threshold level for floodfill segmentation: ");
          cutoff = readln_double();

          if (cutoff < DEFAULT_ZERO) {
            cutoff = DEFAULT_ZERO;
            printf("voledit> Threshold set to default (minimum) value %f \n", cutoff);
          }
          if (cutoff < mindensity) printf("voledit> Warning: Threshold is below minimum density level \n");
          if (cutoff > maxdensity) printf("voledit> Warning: Threshold is above maximum density level \n");

          /* input startpoint */
          printf("voledit> Enter floodfill start index x (1-%d): ", extx);
          xstart = readln_int() - 1;
          if (xstart < 0) xstart = 0;
          if (xstart >= extx) xstart = extx - 1;
          printf("voledit> Enter floodfill start index y (1-%d): ", exty);
          ystart = readln_int() - 1;
          if (ystart < 0) ystart = 0;
          if (ystart >= exty) ystart = exty - 1;
          printf("voledit> Enter floodfill start index z (1-%d): ", extz);
          zstart = readln_int() - 1;
          if (zstart < 0) zstart = 0;
          if (zstart >= extz) zstart = extz - 1;
          printf("voledit> 3D coordinates of start grid point (%d,%d,%d): (%4.2f,%4.2f,%4.2f)\n",
                 xstart + 1, ystart + 1, zstart + 1, origx + width * xstart, origy + width * ystart, origz + width * zstart);

          printf("voledit> Shrink map about found volume? Choose one of the following options: \n");
          printf("voledit> \n");
          printf("voledit>      1: No (keep input map size) \n");
          printf("voledit>      2: Yes (shrink map) \n");
          shrinkflag = readln_int();
          if (shrinkflag != 1 && shrinkflag != 2) {
            error_option(44010, "voledit");
          }

          /* allocate memory and initialize integer array */
          iphi = (char *) alloc_vect(nvox, sizeof(char));
          /* find start value closest to start index */
          pmindist = 1e20;
          xmindist = -1;
          ymindist = -1;
          zmindist = -1;
          for (indz = 0; indz < extz; indz++)
            for (indy = 0; indy < exty; indy++)
              for (indx = 0; indx < extx; indx++)
                if (*(pphi + gidz_general(indz, indy, indx, exty, extx)) >= cutoff) {
                  pcurrdist = (indx - xstart) * (indx - xstart) + 
                              (indy - ystart) * (indy - ystart) + 
                              (indz - zstart) * (indz - zstart);
                  if (pcurrdist < pmindist) {
                    pmindist = pcurrdist;
                    xmindist = indx;
                    ymindist = indy;
                    zmindist = indz;
                  }
                }

          if (xmindist == -1) {
            printf("voledit> No volume above density threshold. Bye bye.\n");
            exit(1);
          }

          if (xmindist != xstart || ymindist != ystart || zmindist != zstart) {
            printf("voledit> Start voxel below threshold. Using nearest voxel above threshold: (%d,%d,%d).\n", xmindist + 1, ymindist + 1, zmindist + 1);
            xstart = xmindist;
            ystart = ymindist;
            zstart = zmindist;
          }

          /* floodfill */
          printf("voledit> Filling contiguous volume...\n");
          flood_fill_nonrecursive(pphi, iphi, cutoff, xstart, ystart, 
                                  zstart, extx, exty, extz);

          /* find neighbor elements of filled volume (these have values below cutoff) */
          for (indz = 0; indz < extz; indz++)
            for (indy = 0; indy < exty; indy++)
              for (indx = 0; indx < extx; indx++)
                if (neighbor(iphi, indx, indy, indz, extx, exty, extz)) 
                  *(iphi + gidz_general(indz, indy, indx, exty, extx)) = 2;

          /* set outliers and neighbors below DEFAULT_ZERO to DEFAULT_ZERO density */
          for (count = 0; count < nvox; count++) {
            if (*(iphi + count) == 0) *(pphi + count) = DEFAULT_ZERO;
            else if (*(iphi + count) == 2 && *(pphi + count) < DEFAULT_ZERO) 
              *(pphi + count) = DEFAULT_ZERO;
          }

          /* allocate memory and initialize integer array */
          iphi = (char *) alloc_vect(nvox, sizeof(char));
          /* find start value closest to start index */
          pmindist = 1e20;
          xmindist = -1;
          ymindist = -1;
          zmindist = -1;
          for (indz = 0; indz < extz; indz++)
            for (indy = 0; indy < exty; indy++)
              for (indx = 0; indx < extx; indx++)
                if (*(pphi + gidz_general(indz, indy, indx, exty, extx)) >= cutoff) {
                  pcurrdist = (indx - xstart) * (indx - xstart) + (indy - ystart) * (indy - ystart) + (indz - zstart) * (indz - zstart);
                  if (pcurrdist < pmindist) {
                    pmindist = pcurrdist;
                    xmindist = indx;
                    ymindist = indy;
                    zmindist = indz;
                  }
                }

          if (xmindist == -1) {
            printf("voledit> No volume above density threshold. Bye bye.\n");
            exit(1);
          }

          if (xmindist != xstart || ymindist != ystart || zmindist != zstart) {
            printf("voledit> Start voxel below threshold. Using nearest voxel above threshold: (%d,%d,%d).\n", xmindist + 1, ymindist + 1, zmindist + 1);
            xstart = xmindist;
            ystart = ymindist;
            zstart = zmindist;
          }

          /* floodfill */
          printf("voledit> Filling contiguous volume...\n");
          flood_fill_nonrecursive(pphi, iphi, cutoff, xstart, ystart, 
                                  zstart, extx, exty, extz);

          /* find neighbor elements of filled volume (these have values below cutoff) */
          for (indz = 0; indz < extz; indz++)
            for (indy = 0; indy < exty; indy++)
              for (indx = 0; indx < extx; indx++)
                if (neighbor(iphi, indx, indy, indz, extx, exty, extz)) 
                  *(iphi + gidz_general(indz, indy, indx, exty, extx)) = 2;

          /* set outliers and neighbors below DEFAULT_ZERO to DEFAULT_ZERO density */
          for (count = 0; count < nvox; count++) {
            if (*(iphi + count) == 0) *(pphi + count) = DEFAULT_ZERO;
            else if (*(iphi + count) == 2 && *(pphi + count) < DEFAULT_ZERO) 
              *(pphi + count) = DEFAULT_ZERO;
          }

          /* determine volumes of contiguous volume and remaining neighbors */
          ncore = 0;
          nngbr = 0;
          for (count = 0; count < nvox; count++) {
            if (*(iphi + count) == 1) ++ncore;
            if (*(iphi + count) == 2 && *(pphi + count) > DEFAULT_ZERO) ++nngbr;
          }
          printf("voledit> Contiguous volume, > %f: %lu voxels (%e Angstrom^3)\n", cutoff, ncore, (ncore * width * width * width));
          printf("voledit> Extracted volume, > %f: %lu voxels (%e Angstrom^3)\n", DEFAULT_ZERO, ncore + nngbr, ((ncore + nngbr) * width * width * width));

          /* determine boundaries of floodfill volume and neighbors */

          minx = extx - 1;
          maxx = 0;
          miny = exty - 1;
          maxy = 0;
          minz = extz - 1;
          maxz = 0;
          for (indz = 0; indz < extz; indz++)
            for (indy = 0; indy < exty; indy++)
              for (indx = 0; indx < extx; indx++)
                if (*(iphi + gidz_general(indz, indy, indx, exty, extx)) > 0) {
                  if (indx < minx) minx = indx;
                  if (indx > maxx) maxx = indx;
                  if (indy < miny) miny = indy;
                  if (indy > maxy) maxy = indy;
                  if (indz < minz) minz = indz;
                  if (indz > maxz) maxz = indz;
                }

          if ((maxx < minx) || (maxy < miny) || (maxz < minz)) {
            error_no_volume(44030, "voledit");
          }

          if (shrinkflag != 2) {
            minx = 0;
            maxx = extx - 1;
            miny = 0;
            maxy = exty - 1;
            minz = 0;
            maxz = extz - 1;
          }
          project_map_lattice(&pphi2, maxx - minx + 1, maxy - miny + 1, 
                              maxz - minz + 1, origx + minx * width, 
                              origy + miny * width, origz + minz * width, 
                              width, width, width, pphi, extx, exty, extz, 
                              origx, origy, origz, width, width, width);
          free_vect_and_zero_ptr(&pphi);
          pphi = pphi2;
          update_parameters(slicemode, &ppro, &rangea, &rangeb, &rangec, 
                            &extx, &exty, &extz, &origx, &origy, &origz,
                            maxx, minx, maxy, miny, maxz, minz, width, &k, 
                            &cutoff, &stride, pphi, &nvox);
          break;


        case 9: /* thresholding */

          printf("voledit> Enter threshold level (densities below this will be set to zero): ");
          cutoff = readln_double();
          threshold(pphi, nvox, cutoff);
          update_parameters(slicemode, &ppro, &rangea, &rangeb, &rangec, 
                            &extx, &exty, &extz, &origx, &origy, &origz,
                            extx - 1, 0, exty - 1, 0, extz - 1, 0, width, 
                            &k, &cutoff, &stride, pphi, &nvox);
          break;

        case 10: /* binary thresholding */

          printf("voledit> Enter binary threshold level (densities below will be set to zero, equal or above to one): ");
          cutoff = readln_double();
          step_threshold(pphi, nvox, cutoff);
          update_parameters(slicemode, &ppro, &rangea, &rangeb, &rangec, 
                            &extx, &exty, &extz, &origx, &origy, &origz,
                            extx - 1, 0, exty - 1, 0, extz - 1, 0, width, 
                            &k, &cutoff, &stride, pphi, &nvox);
          cutoff_old = cutoff;
          break;

        case 11: /* save slice */

          printf("voledit> Choose one of the following three slice output options -\n");
          printf("voledit>      1: three text columns: %c index, %c index, and densities\n", cha, chb);
          printf("voledit>      2: densities only, %c index changing fastest, newline after each element\n", cha);
          printf("voledit>      3: densities only, one row per %c index, newline after each row \n", chb);
          printf("voledit> ");
          slitextmode = readln_int();

          fout = get_save_file_info(done, out_file, slicefile);
          if (k >= 0) { /* save cross section */

            switch (slitextmode) {
              default:
                printf("voledit> Did not recognize option, using three text columns\n");
              case 1:
                for (j = 0; j < rangeb; ++j) for (i = 0; i < rangea; ++i) {
                    ind1 = pindex(slicemode, i, j, k, extx, exty, extz);
                    fprintf(fout, "%4d %4d %lf\n", i + 1, j + 1, *(pphi + ind1));
                  }
                printf("voledit> Saved cross section to file %s in three columns: %c, %c, density.\n", out_file, cha, chb);
                break;
              case 2:
                for (j = 0; j < rangeb; ++j) for (i = 0; i < rangea; ++i) {
                    ind1 = pindex(slicemode, i, j, k, extx, exty, extz);
                    fprintf(fout, "%lf\n", *(pphi + ind1));
                  }
                printf("voledit> Saved cross section to file %s in one density column, %c index changing fastest.\n", out_file, cha);
                break;
              case 3:
                for (j = 0; j < rangeb; ++j) {
                  for (i = 0; i < rangea; ++i) {
                    ind1 = pindex(slicemode, i, j, k, extx, exty, extz);
                    fprintf(fout, "%lf ", *(pphi + ind1));
                  }
                  fprintf(fout, "\n");
                }
                printf("voledit> Saved cross section to file %s in row format, one line per %c index.\n", out_file, chb);
                break;
            }

          } else { /* save projection */

            switch (slitextmode) {
              default:
                printf("voledit> Did not recognize option, using three text columns\n");
              case 1:
                for (j = 0; j < rangeb; ++j) for (i = 0; i < rangea; ++i) {
                    ind1 = i + j * rangea;
                    fprintf(fout, "%4d %4d %lf\n", i + 1, j + 1, *(ppro + ind1));
                  }
                printf("voledit> Saved projection to file %s in three columns: %c, %c, density.\n", out_file, cha, chb);
                break;
              case 2:
                for (j = 0; j < rangeb; ++j) for (i = 0; i < rangea; ++i) {
                    ind1 = i + j * rangea;
                    fprintf(fout, "%lf\n", *(ppro + ind1));
                  }
                printf("voledit> Saved projection to file %s in one density column, %c index changing fastest.\n", out_file, cha);
                break;
              case 3:
                for (j = 0; j < rangeb; ++j) {
                  for (i = 0; i < rangea; ++i) {
                    ind1 = i + j * rangea;
                    fprintf(fout, "%lf ", *(ppro + ind1));
                  }
                  fprintf(fout, "\n");
                }
                printf("voledit> Saved projection to file %s in row format, one line per %c index.\n", out_file, chb);
                break;
            }
          }
          fclose(fout);
          break;

        case 12: /* save map */

          fout = get_save_file_info(done, out_file, mapfile);
          write_vol(out_file, width, origx, origy, origz, extx, exty, extz, pphi);
          break;

        case 13: /* quit */

          printf("voledit> Bye bye!\n");
          done = 1;
          break;
      }
    }
    k += k_shift;
    if (k < -1) k = -1;
    if (k >= rangec) k = rangec - 1;
    if (k == -1 && k_old >= 0) { /* switching from cross-section to projection */
      cutoff_old = cutoff;
      mindensity = calc_min(ppro, (unsigned long)(rangea * rangeb));
      maxdensity = calc_max(ppro, (unsigned long)(rangea * rangeb));
      cutoff = 0.5 * (maxdensity + mindensity);
    }
    if (k >= 0 && k_old == -1) cutoff = cutoff_old; /* switching from projection back to cross-section */
  }
  return 0;
}


static unsigned long pindex(int slicemode, int i, int j, int k, int extx, int exty, int extz)
{

  /* returns index of cross section point */

  switch (slicemode) {
    case 1:
      return i + j * extx + k * extx * exty;
    case 2:
      return j + k * extx + i * extx * exty;
    case 3:
      return k + i * extx + j * extx * exty;
    default:
      error_option(51610, "voledit");
      return 0;
  }
}


static char check_inside_polygon(int poly[NUM_VERTEX][2], int n, int x, int y)
{

  /* Polygon clipping, based on page 244 of Computational Geometry in C */
  /* (Joseph O'Rourke), 2nd ed., 1999 reprint                           */

  int i, i1;
  int px, py, px1, py1;
  double ltersection; /* x intersection of e with ray */
  int r_cross = 0; /* number of right edge/ray crossings */
  int l_cross = 0; /* number of left edge/ray crossings */

  /* for each vertex */
  for (i = 0; i < n; i++) {
    /* center the polygon vertex on the point */
    px = poly[i][0] - x;
    py = poly[i][1] - y;

    /* check if the point is a vertex */
    if (px == 0 && py == 0) return 'v';
    i1 = (i + n - 1) % n;
    px1 = poly[i1][0] - x;
    py1 = poly[i1][1] - y;

    if ((py > 0) != (py1 > 0)) {
      ltersection = (px * py1 - px1 * py) / (double)(py1 - py);
      if (ltersection > 0) r_cross++;
    }

    if ((py  < 0) != (py1 < 0)) {
      ltersection = (px * py1 - px1 * py) / (double)(py1 - py);
      if (ltersection < 0) l_cross++;
    }
  }

  /* point on the edge if left and right cross counts are not the same parity. */
  if ((r_cross % 2) != (l_cross % 2))  return 'e';

  /* point inside if an odd number of crossings. */
  if ((r_cross % 2) == 1) return 'i';
  else  return 'o';
}

static void zero_pixel_in_all_slices(int i, int j, int indz, int slicemode, 
                                     int extx, int exty, int extz, 
                                     int rangec, double *pphi)
{

  /* sets same pixel in all slices to DEFAULT_ZERO */

  unsigned long ind1;

  for (indz = 0; indz < rangec; ++indz) {
    ind1 = pindex(slicemode, i, j, indz, extx, exty, extz);
    *(pphi + ind1) = DEFAULT_ZERO;
  }
}

static void update_parameters(int slicemode, double **ppro, int *rangea, 
                              int *rangeb, int *rangec, unsigned *extx,
                              unsigned *exty, unsigned *extz, double *origx, 
                              double *origy, double *origz,int maxx, int minx, 
                              int maxy, int miny, int maxz, int minz, 
                              double width, int *k, double *cutoff, 
                              int *stride, double *pphi, unsigned long *nvox)
{

  /* updates map and display parameters and the projection */


  double maxdensity, mindensity;
  int i, j, l;
  unsigned long ind1, ind2;

  *extx = maxx - minx + 1;
  *exty = maxy - miny + 1;
  *extz = maxz - minz + 1;
  *origx += minx * width;
  *origy += miny * width;
  *origz += minz * width;

  switch (slicemode) {
    case 1:
      *rangea = *extx;
      *rangeb = *exty;
      *rangec = *extz;
      break;
    case 2:
      *rangea = *extz;
      *rangeb = *extx;
      *rangec = *exty;
      break;
    case 3:
      *rangea = *exty;
      *rangeb = *extz;
      *rangec = *extx;
      break;
  }
  *k = (*rangec) / 2;
  *nvox = (*extx) * (*exty) * (*extz);
  mindensity = calc_min(pphi, *nvox);
  maxdensity = calc_max(pphi, *nvox);
  *cutoff = 0.5 * (maxdensity + mindensity);
  *stride = ceil(*rangea / 100.0);
  if (*stride < ceil(*rangeb / 100.0)) *stride = ceil(*rangeb / 100.0);
  if (*stride < 1) *stride = 1;
  free_vect_and_zero_ptr(ppro);
  do_vect(ppro, (unsigned long)(*rangea * *rangeb));
  for (i = 0; i < *rangea; i++) for (j = 0; j < *rangeb; j++) {
      ind2 = i + j * *rangea;
      for (l = 0; l < *rangec; l++) {
        ind1 = pindex(slicemode, i, j, l, *extx, *exty, *extz);
        *(*ppro + ind2) += *(pphi + ind1);
      }
      *(*ppro + ind2) /= ((double) * rangec);
    }
}

static FILE *get_save_file_info(int done, char *out_file, char *outputfile)
{

  /* get file name if output file */
  /* check if already exists */
  /* seems redundant with other library functions */

  FILE *fout;
  int ch1, ch2;

  for (done = 0; done == 0;) {
    printf("voledit> Enter filename for the %s", outputfile);
    if (fgets(out_file, FLENGTH, stdin) == NULL) {
      error_read_filename(51400, "voledit");
    }
    removespaces(out_file, FLENGTH);
    fout = fopen(out_file, "r");
    if (fout != NULL) {
      printf("voledit> Warning: File exists: %s \n", out_file);
      printf("voledit> Do you want to overwrite? (yes/no) ");
      ch1 = getchar();
      for (;;)  {
        ch2 = getchar();
        if (ch2 == EOF) {
          error_EOF(51410, "voledit");
        }
        if (ch2 == '\n') break;
      }
      if (ch1 == 'y' || ch1 == 'Y') done = 1;
    } else {
      done = 1;
    }
  }
  if ((fout = fopen(out_file, "w")) == NULL) {
    error_read_filename(51420, "voledit");
  }
  return fout;
}

static void get_pad_parameters(int *minx, int *miny, int *minz, int *maxx, 
                               int *maxy, int *maxz, int extx, int exty, 
                               int extz)
{

  /* gets parameters for zero padding */

  printf("voledit> Add zero padded margin to map. Enter integer voxel values:\n");
  printf("voledit> Lower x margin (>= 0): ");
  *minx = 0 - readln_int();
  if (*minx > 0) *minx = 0;
  printf("voledit> Upper x margin (>= 0): ");
  *maxx = extx - 1 + readln_int();
  if (*maxx < extx) *maxx = extx - 1;
  printf("voledit> Lower y margin (>= 0): ");
  *miny = 0 - readln_int();
  if (*miny > 0) *miny = 0;
  printf("voledit> Upper y margin (>= 0): ");
  *maxy = exty - 1 + readln_int();
  if (*maxy < exty) *maxy = exty - 1;
  printf("voledit> Lower z margin (>= 0): ");
  *minz = 0 - readln_int();
  if (*minz > 0) *minz = 0;
  printf("voledit> Upper z margin (>= 0): ");
  *maxz = extz - 1 + readln_int();
  if (*maxz < extz) *maxz = extz - 1;
}

static void get_crop_parameters(int *minx, int *miny, int *minz, int *maxx, 
                                int *maxy, int *maxz, int extx, int exty, 
                                int extz)
{

  /* gets parameters for cropping */

  printf("voledit> Voxel range of cropping. Enter new grid values:\n");
  printf("voledit> Minimum x value (1-%d): ", extx);
  *minx = readln_int() - 1;
  if (*minx < 0) *minx = 0;
  if (*minx >= extx) *minx = extx - 1;
  printf("voledit> Maximum x value (%d-%d): ", *minx + 1, extx);
  *maxx = readln_int() - 1;
  if (*maxx < *minx) *maxx = *minx;
  if (*maxx >= extx) *maxx = extx - 1;
  printf("voledit> Minimum y value (1-%d): ", exty);
  *miny = readln_int() - 1;
  if (*miny < 0) *miny = 0;
  if (*miny >= exty) *miny = exty - 1;
  printf("voledit> Maximum y value (%d-%d): ", *miny + 1, exty);
  *maxy = readln_int() - 1;
  if (*maxy < *miny) *maxy = *miny;
  if (*maxy >= exty) *maxy = exty - 1;
  printf("voledit> Minimum z value (1-%d): ", extz);
  *minz = readln_int() - 1;
  if (*minz < 0) *minz = 0;
  if (*minz >= extz) *minz = extz - 1;
  printf("voledit> Maximum z value (%d-%d): ", *minz + 1, extz);
  *maxz = readln_int() - 1;
  if (*maxz < *minz) *maxz = *minz;
  if (*maxz >= extz) *maxz = extz - 1;
}


static void flood_fill_nonrecursive(double *pphi, char *iphi, double cutoff,  
                                    int l, int m, int n, int extx, int exty, 
                                    int extz)
{

  /* non-recursive floodfill function */

  int idx, x, y, z;
  int count = 0;

  LLIST *list_head, *list_tail, *list_current;

  list_head       = (LLIST *) alloc_vect(1, sizeof(LLIST));
  list_head->x    = l;
  list_head->y    = m;
  list_head->z    = n;
  list_head->next = NULL;
  list_tail       = list_head;

  while (list_head != NULL) {

    count ++;
    x = list_head->x;
    y = list_head->y;
    z = list_head->z;

    idx = extx * exty * z + extx * y + x;

    if ((pphi[idx] >= cutoff) && (iphi[idx] == 0)) {

      iphi[idx] = 1;

      if (x > 0) {
        flood_fill_nonrecursive_push(&list_tail, x - 1, y, z);
        if (y > 0) {
          flood_fill_nonrecursive_push(&list_tail, x - 1, y - 1, z);
          if (z > 0) 
            flood_fill_nonrecursive_push(&list_tail, x - 1, y - 1, z - 1);
          if ((z + 1) < extz) 
            flood_fill_nonrecursive_push(&list_tail, x - 1, y - 1, z + 1);
        }
        if ((y + 1) < exty) {
          flood_fill_nonrecursive_push(&list_tail, x - 1, y + 1, z);
          if (z > 0) 
            flood_fill_nonrecursive_push(&list_tail, x - 1, y + 1, z - 1);
          if ((z + 1) < extz) 
            flood_fill_nonrecursive_push(&list_tail, x - 1, y + 1, z + 1);
        }
        if (z > 0) 
          flood_fill_nonrecursive_push(&list_tail, x - 1, y, z - 1);
        if ((z + 1) < extz) 
          flood_fill_nonrecursive_push(&list_tail, x - 1, y, z + 1);
      }
      if ((x + 1) < extx) {
        flood_fill_nonrecursive_push(&list_tail, x + 1, y, z);
        if (y > 0) {
          flood_fill_nonrecursive_push(&list_tail, x + 1, y - 1, z);
          if (z > 0) 
            flood_fill_nonrecursive_push(&list_tail, x + 1, y - 1, z - 1);
          if ((z + 1) < extz) 
            flood_fill_nonrecursive_push(&list_tail, x + 1, y - 1, z + 1);
        }
        if ((y + 1) < exty) {
          flood_fill_nonrecursive_push(&list_tail, x + 1, y + 1, z);
          if (z > 0) 
            flood_fill_nonrecursive_push(&list_tail, x + 1, y + 1, z - 1);
          if ((z + 1) < extz) 
            flood_fill_nonrecursive_push(&list_tail, x + 1, y + 1, z + 1);
        }
        if (z > 0) 
          flood_fill_nonrecursive_push(&list_tail, x + 1, y, z - 1);
        if ((z + 1) < extz) 
          flood_fill_nonrecursive_push(&list_tail, x + 1, y, z + 1);
      }
      if (y > 0) {
        flood_fill_nonrecursive_push(&list_tail, x, y - 1, z);
        if (z > 0) flood_fill_nonrecursive_push(&list_tail, x, y - 1, z - 1);
        if ((z + 1) < extz) 
          flood_fill_nonrecursive_push(&list_tail, x, y - 1, z + 1);
      }
      if ((y + 1) < exty) {
        flood_fill_nonrecursive_push(&list_tail, x, y + 1, z);
        if (z > 0) 
          flood_fill_nonrecursive_push(&list_tail, x, y + 1, z - 1);
        if ((z + 1) < extz) 
          flood_fill_nonrecursive_push(&list_tail, x, y + 1, z + 1);
      }
      if (z > 0) 
        flood_fill_nonrecursive_push(&list_tail, x, y, z - 1);
      if ((z + 1) < extz) 
        flood_fill_nonrecursive_push(&list_tail, x, y, z + 1);

    }

    list_current = list_head;
    list_head    = list_head->next;
    free_vect_and_zero_ptr(&list_current);

  }

  return;

}

static void flood_fill_nonrecursive_push(LLIST **tail, int x, int y, int z)
{

  /* adds new voxels which need to be visited to the linked list */

  LLIST *current;

  current = (LLIST *) alloc_vect(sizeof(LLIST), 1);

  current->next = NULL;
  current->x    = x;
  current->y    = y;
  current->z    = z;

  (*tail)->next = current;
  (*tail)       = current;

  return;

}

static int neighbor(char *iphi, int l, int m, int n, int extx, int exty, int extz)
{

  /* checks if index element is neighbor of floodfill element */

  if (*(iphi + gidz_general(n, m, l, exty, extx)) == 1) 
    return 0;

  if (l > 0) {
    if (*(iphi + gidz_general(n, m, l - 1, exty, extx)) == 1) 
      return 1;
    if (m > 0) {
      if (*(iphi + gidz_general(n, m - 1, l - 1, exty, extx)) == 1) 
        return 1;
      if (n > 0) 
        if (*(iphi + gidz_general(n - 1, m - 1, l - 1, exty, extx)) == 1) 
          return 1;
      if ((n + 1) < extz) 
        if (*(iphi + gidz_general(n + 1, m - 1, l - 1, exty, extx)) == 1) 
          return 1;
    }
    if ((m + 1) < exty) {
      if (*(iphi + gidz_general(n, m + 1, l - 1, exty, extx)) == 1) 
        return 1;
      if (n > 0) 
        if (*(iphi + gidz_general(n - 1, m + 1, l - 1, exty, extx)) == 1) 
          return 1;
      if ((n + 1) < extz) 
        if (*(iphi + gidz_general(n + 1, m + 1, l - 1, exty, extx)) == 1) 
          return 1;
    }
    if (n > 0) 
      if (*(iphi + gidz_general(n - 1, m, l - 1, exty, extx)) == 1) 
        return 1;
    if ((n + 1) < extz) 
      if (*(iphi + gidz_general(n + 1, m, l - 1, exty, extx)) == 1) 
        return 1;
  }
  if ((l + 1) < extx) {
    if (*(iphi + gidz_general(n, m, l + 1, exty, extx)) == 1) 
      return 1;
    if (m > 0) {
      if (*(iphi + gidz_general(n, m - 1, l + 1, exty, extx)) == 1) 
        return 1;
      if (n > 0) 
        if (*(iphi + gidz_general(n - 1, m - 1, l + 1, exty, extx)) == 1) 
          return 1;
      if ((n + 1) < extz) 
        if (*(iphi + gidz_general(n + 1, m - 1, l + 1, exty, extx)) == 1) 
          return 1;
    }
    if ((m + 1) < exty) {
      if (*(iphi + gidz_general(n, m + 1, l + 1, exty, extx)) == 1) 
        return 1;
      if (n > 0) 
        if (*(iphi + gidz_general(n - 1, m + 1, l + 1, exty, extx)) == 1) 
          return 1;
      if ((n + 1) < extz) 
        if (*(iphi + gidz_general(n + 1, m + 1, l + 1, exty, extx)) == 1) 
          return 1;
    }
    if (n > 0) 
      if (*(iphi + gidz_general(n - 1, m, l + 1, exty, extx)) == 1) 
        return 1;
    if ((n + 1) < extz) 
      if (*(iphi + gidz_general(n + 1, m, l + 1, exty, extx)) == 1) 
        return 1;
  }
  if (m > 0) {
    if (*(iphi + gidz_general(n, m - 1, l, exty, extx)) == 1) 
      return 1;
    if (n > 0) 
      if (*(iphi + gidz_general(n - 1, m - 1, l, exty, extx)) == 1) 
        return 1;
    if ((n + 1) < extz) 
      if (*(iphi + gidz_general(n + 1, m - 1, l, exty, extx)) == 1) 
        return 1;
  }
  if ((m + 1) < exty) {
    if (*(iphi + gidz_general(n, m + 1, l, exty, extx)) == 1) 
      return 1;
    if (n > 0) 
      if (*(iphi + gidz_general(n - 1, m + 1, l, exty, extx)) == 1) 
        return 1;
    if ((n + 1) < extz) 
      if (*(iphi + gidz_general(n + 1, m + 1, l, exty, extx)) == 1) 
        return 1;
  }
  if (n > 0) 
    if (*(iphi + gidz_general(n - 1, m, l, exty, extx)) == 1) 
      return 1;
  if ((n + 1) < extz) 
    if (*(iphi + gidz_general(n + 1, m, l, exty, extx)) == 1) 
      return 1;
  return 0;
}

/*********************************************************************
*                          V O L H I S T                             *
**********************************************************************
* Program is part of the Situs package, URL: situs.biomachina.org    *
* (c) Willy Wriggers and Paul Boyle 1998-2010                        *
**********************************************************************
*                                                                    *
* Map rescaling, shifting of background peak, histogram match.       *
*                                                                    *
**********************************************************************
* See legal statement for terms of distribution                      *
*********************************************************************/

#include "situs.h"
#include "lib_vio.h"
#include "lib_std.h"
#include "lib_vwk.h"
#include "lib_vec.h"

int main(int argc, char *argv[])
{

  double *phi1, *phi2, *phi2pre, *phi3;
  unsigned extx1, exty1, extz1;
  unsigned extx2, exty2, extz2;
  unsigned long nvox1, nvox2, count;
  double origx1, origy1, origz1, width1;
  double origx2, origy2, origz2, width2;
  double currdensity, gain, bias, bias_old, surf1, surf2, min1, min2, max1, max2;
  int nbins, done;

  if (argc < 2 || argc > 4) {
    fprintf(stderr, "volhist> Usage: volhist infile1 [[infile2] outfile]  (all files are density maps, [...] = optional)\n");
    fprintf(stderr, "volhist> Modes:\n");
    fprintf(stderr, "volhist>        volhist infile1                      (print histogram) \n");
    fprintf(stderr, "volhist>        volhist infile1 outfile              (rescale or shift densities) \n");
    fprintf(stderr, "volhist>        volhist infile1 infile2 outfile      (match histogram of infile2 to that of infile1) \n");
    exit(1);
  }

  /* mode 1, read infile1 and print histogram */
  if (argc == 2) {
    read_vol(argv[1], &width1, &origx1, &origy1, &origz1, &extx1, &exty1, &extz1, &phi1);
    nvox1 = extx1 * exty1 * extz1;
    nbins = print_histogram(&extx1, &exty1, &extz1, &phi1, -1);
  }

  /* mode 2, rescale and shift infile1 and write to outfile */
  if (argc == 3) {
    read_vol(argv[1], &width1, &origx1, &origy1, &origz1, &extx1, &exty1, &extz1, &phi1);
    nvox1 = extx1 * exty1 * extz1;
    nbins = print_histogram(&extx1, &exty1, &extz1, &phi1, -1);
    printf("volhist> Enter a scaling factor (gain) by which map densities will be multiplied: ");
    gain = readln_double();
    for (count = 0; count < nvox1; count++) *(phi1 + count) *= gain;
    printf("volhist> Enter offset density value (bias), this will be added after scaling by gain factor: ");
    bias = readln_double();
    for (count = 0; count < nvox1; count++) *(phi1 + count) += bias;
    printf("volhist> Calculating new voxel histogram\n");
    print_histogram(&extx1, &exty1, &extz1, &phi1, nbins);
    write_vol(argv[2], width1, origx1, origy1, origz1, extx1, exty1, extz1, phi1);
  }

  /* mode 3, match histogram of infile2 to that of infile1 and write to outfile */
  if (argc == 4) {

    printf("volhist> Exploring difference histogram of (isovalue thresholded) infile1 - (infile2 * gain + bias) \n");

    /* read infile1 and project infile2 to infile1 lattice */
    read_vol(argv[1], &width1, &origx1, &origy1, &origz1, &extx1, &exty1, &extz1, &phi1);
    nvox1 = extx1 * exty1 * extz1;
    nbins = print_histogram(&extx1, &exty1, &extz1, &phi1, -1);
    read_vol(argv[2], &width2, &origx2, &origy2, &origz2, &extx2, &exty2, &extz2, &phi2pre);
    nvox2 = extx2 * exty2 * extz2;
    print_histogram(&extx2, &exty2, &extz2, &phi2pre, nbins);
    project_map_lattice(&phi2, extx1, exty1, extz1, origx1, origy1, origz1,
                        width1, width1, width1, phi2pre, extx2, exty2, extz2,
                        origx2, origy2, origz2, width2, width2, width2);

    /* get surface isovalues, with sanity checks */
    printf("volhist> Enter surface isovalue (cutoff) > 0 for infile1 %s: ", argv[1]);
    surf1 = readln_double();
    if (surf1 <= 0) {
      fprintf(stderr, "volhist> Error: surface isovalue must be > 0 [e.c. 12001]\n");
      exit(12001);
    }
    min1 = calc_min(phi1, nvox1);
    max1 = calc_max(phi1, nvox1);
    if (surf1 < min1 || surf1 > max1) {
      fprintf(stderr, "volhist> Error: surface isovalue is outside density range of file %s [e.c. 12002]\n", argv[1]);
      exit(12002);
    }
    printf("volhist> Enter surface isovalue (cutoff) > 0 for infile2 %s: ", argv[2]);
    surf2 = readln_double();
    if (surf2 <= 0) {
      fprintf(stderr, "volhist> Error: surface isovalue must be > 0 [e.c. 12003]\n");
      exit(12003);
    }
    min2 = calc_min(phi2, nvox1);
    max2 = calc_max(phi2, nvox1);
    if (surf2 < min2 || surf2 > max2) {
      fprintf(stderr, "volhist> Error: surface isovalue is outside density range of file %s [e.c. 12004]\n", argv[2]);
      exit(12004);
    }

    /* threshold input maps */
    threshold(phi1, nvox1, surf1);
    threshold(phi2, nvox1, surf2);

    /* loop until satisfying bias is found */
    bias = 0;
    do_vect(&phi3, nvox1);
    for (done = 0; done == 0;) {
      gain = (surf1 - bias) / surf2;
      for (count = 0; count < nvox1; count++) {
        currdensity = *(phi2 + count) * gain + bias; /* affine projection of infile2 */
        if (currdensity < surf1) currdensity = 0; /* threshold, i.e. remove bias from projected infile2 */
        *(phi3 + count) =  *(phi1 + count) - currdensity; /* get difference map */
      }
      printf("volhist> Computing difference histogram of (isovalue thresholded) infile1 - (infile2 * gain + bias).\n");
      print_diff_histogram(&extx1, &exty1, &extz1, &phi3, nbins);
      printf("volhist> Gain, bias used for affine transformation that matches infile2 to infile1: %f , %f \n", gain, bias);
      printf("volhist> Enter new bias value < %f to try again, or any value >= %f to exit and save current transformation: ", surf1, surf1);
      bias_old = bias;
      bias = readln_double();
      if (bias >= surf1) {
        /* write output with limited zero thresholding only, user is responsible to properly segment (threshold) all maps with voledit / floodfill */
        if (bias_old < 0) printf("volhist> Writing non-negative densities to histogram matched outfile = infile2 * gain + bias using the following gain, bias: %f , %f \n", gain, bias_old);
        else printf("volhist> Writing densities to histogram matched outfile = infile2 * gain + bias using the following gain, bias: %f , %f \n", gain, bias_old);
        for (count = 0; count < nvox2; count++) *(phi2pre + count) = *(phi2pre + count) * gain + bias_old;
        threshold(phi2pre, nvox2, 0.0);
        write_vol(argv[3], width2, origx2, origy2, origz2, extx2, exty2, extz2, phi2pre);
        printf("volhist> New surface threshold for outfile %s: %f \n", argv[3], surf1);
        if (bias_old < 0) printf("volhist> Negative densities after affine transformation set to zero.\n");
        if (bias_old > 0) printf("volhist> Positive background level added by bias (keep in mind when subtracting or editing map).\n");
        printf("volhist> Bye bye!\n");
        done = 1;
      }
    }
  }
  return 0;
}

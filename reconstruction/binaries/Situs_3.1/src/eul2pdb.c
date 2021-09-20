/*********************************************************************
*                          E U L 2 P D B                             *
**********************************************************************
* Program is part of the Situs package (c) Valerio Mariani, 2005     *
* URL: situs.biomachina.org                                          *
**********************************************************************
*                                                                    *
* Generates PDB files for inspection of colores Euler angles.        *
*                                                                    *
**********************************************************************
* See legal statement for terms of distribution                      *
*********************************************************************/


#include "situs.h"
#include "lib_err.h"
#include "lib_pio.h"
#include "lib_pwk.h"
#include "lib_eul.h"
#include "lib_vec.h"
#include "lib_vwk.h"
#include "lib_std.h"

int main(int argc, char **argv)
{

  char *program = "eul2pdb";
  char *file1 = "eulers";
  char *file2 = "pdb";
  FILE   *filein;
  FILE   *fileout;
  PDB  *fake_pdb;
  PDB  *fake_pdb2;
  int u;
  char    record[41];
  char    field[41];
  double  psi = 0, theta = 0, phi = 0;

  if (argc != 3) {
    error_IO_files_2(program, file1, file2);
    exit(1);
  }

  filein = fopen(argv[1], "r");
  if (filein == NULL) {
    error_open_filename(99999, program, argv[1]);
  }

  fileout = fopen(argv[2], "w");
  if (fileout == NULL) {
    error_open_filename(99999, program, argv[2]);
  }

  fake_pdb = (PDB *) alloc_vect(1, sizeof(PDB));
  fake_pdb2 = (PDB *) alloc_vect(1, sizeof(PDB));

  u = 0;

  fake_pdb[0].x = 0.0;
  fake_pdb[0].y = 0.0;
  fake_pdb[0].z = 10.0;

  while (fgets(record, 41, filein) != NULL) {
    get_fld(record,  1, 10, field);
    psi   = (double) atof(field);
    get_fld(record, 11, 20, field);
    theta = (double) atof(field);
    get_fld(record, 21, 30, field);
    phi   = (double) atof(field);

    rot_euler(fake_pdb, fake_pdb2, 1, psi * ROT_CONV, theta * ROT_CONV, phi * ROT_CONV);

    fprintf(fileout, "ATOM %5d  EU  EUL E        %8.3f%8.3f%8.3f  1.00%6.2f           E\n", u,
            fake_pdb2[0].x, fake_pdb2[0].y, fake_pdb2[0].z, phi * ROT_CONV);
    u++;
  }

  printf("eul2pdb> %d triplets of Euler angles processed and written to the file: %s\n", u, argv[2]);

  free_vect_and_zero_ptr(&fake_pdb);
  free_vect_and_zero_ptr(&fake_pdb2);
  fclose(filein);
  fclose(fileout);
  return 0;
}

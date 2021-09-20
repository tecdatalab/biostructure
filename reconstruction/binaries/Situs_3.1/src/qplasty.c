/*********************************************************************
*                          Q P L A S T Y                             *
**********************************************************************
* Program is part of the Situs package, URL: situs.biomachina.org    *
* (c) Mirabela Rusu and Willy Wriggers 2001-2008                     *
**********************************************************************
*                                                                    *
* Interpolation of sparsely sampled displacements.                   *
*                                                                    *
**********************************************************************
* See legal statement for terms of distribution                      *
*********************************************************************/


#include "situs.h"
#include "lib_std.h"
#include "lib_pio.h"
#include "lib_vec.h"

#define FLENGTH 1000            /* file name length */

/* function declarations */
//compute matrix inverse
double determinator(double, double, double, double, double, double, double, double, double, double, double, double);
double detsarrus(double, double, double, double, double, double, double, double, double);
void gaussjordan(double *, int);

// write matrices to files.
void print_matrix(double *, int, int, char *);

//interpolations
void fill_matrix(double *, int, int, int, int, int, int, double , double , double, double , int, int,  double);
double cartesian_displacement(PDB, PDB, char);
double get_weight(double, int , double, double);
int get_calpha_index(int, int, PDB *);
void interpolate_splines(int, PDB *, int, PDB *, int, PDB *, int, PDB *, int, int, double, int);
void interpolate_IDW(int, PDB *, int, PDB *, int, PDB *, int, PDB *, int, double, int, int);


/* the main-program */
int main(int argc, char *argv[])
{
  unsigned int numA, numB, numV, numQ;
  PDB *pdbA, *pdbB, *pdbV, *pdbQ;
  int i, j, inmode, outmode, intermode;
  double vectorscale, vectorrmsd;
  int interpol_method, kernel_type, weighting_type, vectors_influence;
  double weighting_pow = -1.00;
  int rigid_residues = -1;
  double elastic_coef;

  /* interactivity mode */
  intermode = 1; // default or -byatom
  if (argc == 6) {
    if (strcmp(argv[5], "-byres") == 0) intermode = 2;
    if (strcmp(argv[5], "-interactive") == 0) intermode = 3;
  }

  if (argc < 5 || argc > 6) {
    fprintf(stderr, "qplasty> Usage: qplasty file1 file2 file3 file4 [options] \n");
    fprintf(stderr, "qplasty>        file1: InAtomFile    (coords)\n");
    fprintf(stderr, "qplasty>        file2: InVectorFile1 (coords)\n");
    if (intermode == 3) {
      fprintf(stderr, "qplasty>        file3: InVectorFile2 (displacements or coords)\n");
      fprintf(stderr, "qplasty>        file4: OutAtomFile   (displacements or coords)\n");
    } else {
      fprintf(stderr, "qplasty>        file3: InVectorFile2 (coords)\n");
      fprintf(stderr, "qplasty>        file4: OutAtomFile   (coords)\n");
    }
    fprintf(stderr, "qplasty>        [options]: optional flag for default parameters or full interactive mode: \n");
    fprintf(stderr, "qplasty>                   <default> or -byatom : global IDW by atom\n");
    fprintf(stderr, "qplasty>                   -byres : global IDW by residue, to reduce the number of broken bonds\n");
    fprintf(stderr, "qplasty>                   -interactive : free choice of parameters\n");
    exit(1);
  }

  read_pdb(argv[1], &numA, &pdbA);
  read_pdb(argv[1], &numB, &pdbB);
  read_pdb(argv[2], &numQ, &pdbQ);
  if (numQ < 4) {
    fprintf(stderr, "qplasty> Error: Must have at least 4 vectors in file %s [e.c. 52320]\n", argv[2]);
    exit(52320);
  }
  read_pdb(argv[3], &numV, &pdbV);
  if (numV != numQ) {
    fprintf(stderr, "qplasty> Error: Must have same number of vectors in %s and %s [e.c. 52330]\n", argv[2], argv[3]);
    exit(52330);
  }

  /* define coordinate input in interactive mode */
  if (intermode == 3) {
    fprintf(stderr, "qplasty> \n");
    fprintf(stderr, "qplasty> What is the content of file %s?\n", argv[3]);
    fprintf(stderr, "qplasty> Choose one of the following options -\n");
    fprintf(stderr, "qplasty>      1: Second set of vector coordinates\n");
    fprintf(stderr, "qplasty>      2: Vector displacements relative to %s\n", argv[2]);
    fprintf(stderr, "qplasty> ");
    inmode = readln_int();
  } else {
    inmode = 1;
  }
  switch (inmode) {
    case 1:
      for (j = 0; j < numV; ++j) { /* make pdbV contain displacements */
        pdbV[j].x -= pdbQ[j].x;
        pdbV[j].y -= pdbQ[j].y;
        pdbV[j].z -= pdbQ[j].z;
      }
      fprintf(stderr, "qplasty> %d vector displacements computed.\n", numV);
      break;
    case 2:
      break;
    default:
      fprintf(stderr, "qplasty> Error: can't identify option [e.c. 52340]\n");
      exit(52340);
  }

  /* rescaling of displacements in interactive mode */
  if (intermode == 3) {
    vectorrmsd = 0;
    for (i = 0; i < numV; ++i) {
      vectorrmsd += pdbV[i].x * pdbV[i].x + pdbV[i].y * pdbV[i].y + pdbV[i].z * pdbV[i].z;
    }
    vectorrmsd /= ((double) numV);
    vectorrmsd = sqrt(vectorrmsd);
    fprintf(stderr, "qplasty> The rms amplitude of the displacements is: %f \n", vectorrmsd);
    fprintf(stderr, "qplasty> Enter desired scaling factor that will be applied to the amplitudes: ");
    vectorscale = readln_double();
  } else {
    vectorscale = 1;
  }
  for (i = 0; i < numQ; ++i) {
    pdbV[i].x *= vectorscale;
    pdbV[i].y *= vectorscale;
    pdbV[i].z *= vectorscale;
  }

  /* choice of residue blocking in interactive mode */
  if (intermode == 3) {
    fprintf(stderr, "qplasty> Consider residues as rigid entities that follow the movement of the alpha charbon?\n");
    fprintf(stderr, "qplasty>      1: Yes (recommended to reduce the number of broken bonds)\n");
    fprintf(stderr, "qplasty>      2: No\n");
    fprintf(stderr, "qplasty> ");
    rigid_residues = readln_int();
  } else {
    if (intermode == 1) rigid_residues = 2;
    if (intermode == 2) rigid_residues = 1;
  }
  if (rigid_residues < 1 ||  rigid_residues > 2) {
    fprintf(stderr, "qplasty> Error: can't identify rigid residue option [e.c. 52350]\n");
    exit(52350);
  }

  /* choice of interpolation function in interactive mode */
  if (intermode == 3) {
    fprintf(stderr, "qplasty> \n");
    fprintf(stderr, "qplasty> What interpolation method to use?\n");
    fprintf(stderr, "qplasty>      1: Inverse Distance Weighting\n");
    fprintf(stderr, "qplasty>      2: Thin Plate Splines\n");
    fprintf(stderr, "qplasty>      3: Elastic Body Splines\n");
    fprintf(stderr, "qplasty> ");
    interpol_method = readln_int();
  } else {
    interpol_method = 1;
  }

  /* select interpolation parameters and kernels as needed */
  switch (interpol_method) {
    case 1: // IDW
      if (intermode == 3) {
        fprintf(stderr, "qplasty> \n");
        fprintf(stderr, "qplasty> What is the type of the weighting scheme?\n");
        fprintf(stderr, "qplasty>      1: Global \n");
        fprintf(stderr, "qplasty>      2: Local \n");
        fprintf(stderr, "qplasty> ");
        weighting_type = readln_int();
      } else {
        weighting_type = 1;
      }
      if (weighting_type < 1 ||  weighting_type > 2) {
        fprintf(stderr, "qplasty> Error: can't identify weighting type [e.c. 52360]\n");
        exit(52360);
      }
      if (weighting_type == 1) {
        if (intermode == 3) {
          fprintf(stderr, "qplasty> \n");
          fprintf(stderr, "qplasty> What is the exponent of the global weighting function (recommended value: 8)?\n");
          fprintf(stderr, "qplasty> ");
          weighting_pow = readln_double();
        } else {
          weighting_pow = 8;
        }
      }
      if (weighting_type == 2) {
        fprintf(stderr, "qplasty> \n");
        fprintf(stderr, "qplasty> What is the exponent of the local weighting function (recommended value: 4)?\n");
        fprintf(stderr, "qplasty> ");
        weighting_pow = readln_double();
        fprintf(stderr, "qplasty> \n");
        fprintf(stderr, "qplasty> How many (closest) feature vectors to use in the interpolation (min=3, max=%d, recommended=%d-%d)?\n", numV, (int)(numV * 0.5), (int)(numV * 0.9));
        fprintf(stderr, "qplasty> ");
        vectors_influence = readln_int();
      } else
        vectors_influence = numV;

      if (weighting_pow <= 1) {
        fprintf(stderr, "qplasty> Error: The weighting exponent should be larger that 1 [e.c. 52370]\n");
        exit(52370);
      }

      if (vectors_influence < 3 ||  vectors_influence > numV) {
        fprintf(stderr, "qplasty> Error: number of feature vectors out of bounds [e.c. 52380]\n");
        exit(52380);
      }

      interpolate_IDW(numA, pdbA, numB, pdbB, numQ, pdbQ, numV, pdbV, weighting_type, weighting_pow, vectors_influence, rigid_residues);
      break;
    case 2:// TPS
      kernel_type = 1;
      interpolate_splines(numA, pdbA, numB, pdbB, numQ, pdbQ, numV, pdbV, interpol_method, kernel_type, 0.0f, rigid_residues);
      break;
    case 3: // EBS
      fprintf(stderr, "qplasty> \n");
      fprintf(stderr, "qplasty> What type of kernel to use?\n");
      fprintf(stderr, "qplasty>      1: U(r) = ((11-12*nu)*|x|^2-3*x*x^T)*|x|\n");
      fprintf(stderr, "qplasty>      2: U(r) = (7-8*nu)*r(x)-x*x^T*1/r(x)  - recommended\n");
      fprintf(stderr, "qplasty> ");
      kernel_type = readln_int();
      if (kernel_type < 1 || kernel_type > 2) {
        fprintf(stderr, "qplasty> Error: can't identify kernel type [e.c. 52390]\n");
        exit(52390);
      }
      fprintf(stderr, "qplasty> \n");
      fprintf(stderr, "qplasty> What is the Poisson ratio nu (range: [0 (soft) - 0.5 (incompressible, recommended)])?\n");
      fprintf(stderr, "qplasty> ");
      elastic_coef = readln_double();
      if (elastic_coef < 0.0 || elastic_coef > 0.5) {
        fprintf(stderr, "qplasty> Error: nu must be between 0 and 0.5 [e.c. 52400]\n");
        exit(52400);
      }
      interpolate_splines(numA, pdbA, numB, pdbB, numQ, pdbQ, numV, pdbV, interpol_method, kernel_type, elastic_coef, rigid_residues);
      break;
    case 4:
      break;
    default:
      fprintf(stderr, "qplasty> Error: can't identify interpolation method [e.c. 52410]\n");
      exit(52410);
  }

  /* define coordinate output in interactive mode */
  if (intermode == 3) {
    fprintf(stderr, "qplasty> \n");
    fprintf(stderr, "qplasty> Do you want to add or export the interpolated displacements?\n");
    fprintf(stderr, "qplasty> Choose one of the following options -\n");
    fprintf(stderr, "qplasty>      1: Add displacements to atomic coordinates from file %s \n", argv[1]);
    fprintf(stderr, "qplasty>      2: Export displacements\n");
    fprintf(stderr, "qplasty> ");
    outmode = readln_int();
  } else {
    outmode = 1;
  }
  switch (outmode) {
    case 1:
      for (j = 0; j < numA; ++j) {
        pdbB[j].x += pdbA[j].x;
        pdbB[j].y += pdbA[j].y;
        pdbB[j].z += pdbA[j].z;
      }
      break;
    case 2:
      break;
    default:
      fprintf(stderr, "qplasty> Error: can't identify option [e.c. 52510]\n");
      exit(52510);
  }
  write_pdb(argv[4], numB, pdbB);
  fprintf(stderr, "qplasty> New PDB coordinates written to file %s.\n", argv[4]);
  return 0;
}




double determinator(double x1, double x2, double x3, double x4, double y1, double y2, double y3, double y4, double z1, double z2, double z3, double z4)
/* computes 4x4 determinant with last column (1,1,1,1) */
{
  return detsarrus(x1, y1, z1, x2, y2, z2, x3, y3, z3) - detsarrus(x1, y1, z1, x2, y2, z2, x4, y4, z4) + detsarrus(x1, y1, z1, x3, y3, z3, x4, y4, z4) - detsarrus(x2, y2, z2, x3, y3, z3, x4, y4, z4);
}

double detsarrus(double a11, double a12, double a13, double a21, double a22, double a23, double a31, double a32, double a33)
/* computes 3x3 determinant */
{
  return a11 * a22 * a33 + a12 * a23 * a31 + a13 * a21 * a32 - a13 * a22 * a31 - a11 * a23 * a32 - a12 * a21 * a33;
}

void gaussjordan(double *mat, int n)
{
  /* nxn real Matrix inversion by Gauss-Jordan elimination */
  int i, j, k, l, m, ic = 0, ir = 0;
  int *coli, *rowi, *pivi;
  double pivinv, currmax, dummy;
  double tmp;

  coli = (int *) alloc_vect(n, sizeof(int));
  rowi = (int *) alloc_vect(n, sizeof(int));
  pivi = (int *) alloc_vect(n, sizeof(int));

  for (i = 0; i < n; i++) {
    currmax = 0;
    for (j = 0; j < n; j++)
      if (pivi[j] != 1)
        for (k = 0; k < n; k++) {
          if (pivi[k] == 0) {
            if (fabs(mat[j * n + k]) >= currmax) {
              currmax = fabs(mat[j * n + k]);
              ir = j;
              ic = k;
            }
          } else if (pivi[k] > 1) {
            fprintf(stderr, "qplasty> Error: Matrix singular [e.c. 52940]\n");
            exit(52940);
          }
        }
    ++(pivi[ic]);
    if (ir != ic) for (l = 0; l < n; l++) {
        tmp = mat[ir * n + l];
        mat[ir * n + l] = mat[ic * n + l];
        mat[ic * n + l] = tmp;
      }
    rowi[i] = ir;
    coli[i] = ic;
    if (mat[ic * n + ic] == 0) {
      fprintf(stderr, "qplasty> Error: Matrix singular [e.c. 52950]\n");
      exit(52950);
    }
    pivinv = 1.0 / mat[ic * n + ic];
    mat[ic * n + ic] = 1;
    for (l = 0; l < n; l++) mat[ic * n + l] *= pivinv;
    for (m = 0; m < n; m++)
      if (m != ic) {
        dummy = mat[m * n + ic];
        mat[m * n + ic] = 0;
        for (l = 0; l < n; l++) mat[m * n + l] -= mat[ic * n + l] * dummy;
      }
  }
  for (l = n - 1; l >= 0; l--) if (rowi[l] != coli[l]) for (k = 0; k < n; k++) {
        tmp = mat[k * n + rowi[l]];
        mat[k * n + rowi[l]] = mat[k * n + coli[l]];
        mat[k * n + coli[l]] = tmp;
      }

  free_vect_and_zero_ptr(&coli);
  free_vect_and_zero_ptr(&rowi);
  free_vect_and_zero_ptr(&pivi);
}

/*
 *   Compute the cartesian displacement two atoms
 *   \param point1
 *   \param point2
 *   \param  axis = 'x','y','z'
 *   \return the value of the displacement on the corresponding axe
 */
double cartesian_displacement(PDB point1, PDB point2, char axis)
{
  switch (axis) {
    case 'x':
      return (point2.x - point1.x);
      break;
    case 'y':
      return (point2.y - point1.y);
      break;
    case 'z':
      return (point2.z - point1.z);
      break;
    default:
      fprintf(stderr, "qplasty> Error: Wrong axis! it has to be of the type 'x','y','z'[e.c. 52600]\n");
      exit(52600);
      break;
  }
}


double cartesian_distance(PDB point1, PDB point2)
{
  return sqrt((point1.x - point2.x) * (point1.x - point2.x) +
              (point1.y - point2.y) * (point1.y - point2.y) +
              (point1.z - point2.z) * (point1.z - point2.z));
}

double min(double a, double b)
{
  double minimum = a < b ? a : b;
  return minimum;
}

double max(double a, double b)
{
  double maximum = a > b ? a : b;
  return maximum;
}

/**
 * Fill the matrix from the position [starting_row,starting_col] till the position
 * [starting_row+rows,starting_col+cols] with the values coresponding to kernel type;
 * \param matrix the matrix to fill
 * \param matrix_rows number of rows of the matrix
 * \param matrix_cols number of columns of the matrix
 * \param starting_row,starting_col position from where the matrix is fill
 * \param rows,cols number of rows and columns that will be filled
 * \param dx,dy,dz the displacements on axes x, y, z
 * \param distance distance between points [starting_row/3, starting_col/3]
 * \param kernel_type = 1 - EBS; 2 - EBS second version; 1 TPS r; 0 add points [dx,dy, dz]*I
 * \param elastic_parameter - poisson ratio varies between 0.01 soft and 0.49 incompresible material
 **/
void fill_matrix(double *matrix, int matrix_rows, int matrix_cols, int starting_row, int starting_col, int rows, int cols, double dx , double dy , double dz, double distance,  int interpol_method, int kernel_type, double elastic_parameter)
{
  if (starting_row > matrix_rows || starting_row + rows > matrix_rows) {
    fprintf(stderr, "qplasty> Error: the starting row position and final row position have to be lower that the row number[e.c. 52610]\n");
    exit(52610);
  }
  if (starting_col > matrix_cols || starting_col + cols > matrix_cols) {
    fprintf(stderr, "qplasty> Error: the starting col position(%d) and final col position(%d+cols(%d)) have to be lower that the col number[e.c. 52620]\n", starting_col, matrix_cols, cols);
    exit(52620);
  }

  int i;
  double value;
  double distance_pow3;
  switch (interpol_method) {
    case 0: //displacement and point matrix
      if (cols == 1) {
        matrix[(starting_row + 0)*matrix_cols + starting_col + 0] = dx;
        matrix[(starting_row + 1)*matrix_cols + starting_col + 0] = dy;
        matrix[(starting_row + 2)*matrix_cols + starting_col + 0] = dz;
      } else {
        // used to fill the point matrix and in this case the displacements represent the values on axes x, y and z
        matrix[(starting_row + 0)*matrix_cols + starting_col + 0] = dx;
        matrix[(starting_row + 1)*matrix_cols + starting_col + 1] = dy;
        matrix[(starting_row + 2)*matrix_cols + starting_col + 2] = dz;
      }
      break;

    case 2://TPS
      value = distance;
      for (i = 0; i < rows; i++)
        matrix[(starting_row + i)*matrix_cols + starting_col + i] = value;
      break;
    case 3:
      switch (kernel_type) {
        case 1:
          distance_pow3 = distance * distance * distance;

          matrix[(starting_row + 0)*matrix_cols + starting_col + 0] = -3.0 * distance * dx * dx + (12.0 * (1.0 - elastic_parameter) - 1.0) * distance_pow3;
          matrix[(starting_row + 0)*matrix_cols + starting_col + 1] = -3.0 * distance * dx * dy;
          matrix[(starting_row + 0)*matrix_cols + starting_col + 2] = -3.0 * distance * dx * dz;

          matrix[(starting_row + 1)*matrix_cols + starting_col + 0] = -3.0 * distance * dy * dx;
          matrix[(starting_row + 1)*matrix_cols + starting_col + 1] = -3.0 * distance * dy * dy + (12.0 * (1.0 - elastic_parameter) - 1.0) * distance_pow3;
          matrix[(starting_row + 1)*matrix_cols + starting_col + 2] = -3.0 * distance * dy * dz;

          matrix[(starting_row + 2)*matrix_cols + starting_col + 0] = -3.0 * distance * dz * dx;
          matrix[(starting_row + 2)*matrix_cols + starting_col + 1] = -3.0 * distance * dz * dy;
          matrix[(starting_row + 2)*matrix_cols + starting_col + 2] = -3.0 * distance * dz * dz + (12.0 * (1.0 - elastic_parameter) - 1.0) * distance_pow3;

          break;
        case 2:

          matrix[(starting_row + 0)*matrix_cols + starting_col + 0] = -1 / distance * dx * dx + (8 * (1 - elastic_parameter) - 1) * distance;
          matrix[(starting_row + 0)*matrix_cols + starting_col + 1] = -1 / distance * dx * dy;
          matrix[(starting_row + 0)*matrix_cols + starting_col + 2] = -1 / distance * dx * dz;

          matrix[(starting_row + 1)*matrix_cols + starting_col + 0] = -1 / distance * dy * dx;
          matrix[(starting_row + 1)*matrix_cols + starting_col + 1] = -1 / distance * dy * dy + (8 * (1 - elastic_parameter) - 1) * distance;
          matrix[(starting_row + 1)*matrix_cols + starting_col + 2] = -1 / distance * dy * dz;

          matrix[(starting_row + 2)*matrix_cols + starting_col + 0] = -1 / distance * dz * dx;
          matrix[(starting_row + 2)*matrix_cols + starting_col + 1] = -1 / distance * dz * dy;
          matrix[(starting_row + 2)*matrix_cols + starting_col + 2] = -1 / distance * dz * dz + (8 * (1 - elastic_parameter) - 1) * distance;

          break;
      }
      break;
  }

}

void print_matrix(double *matrix, int matrix_rows, int matrix_cols, char *filename)
{
  FILE *file;
  file = fopen(filename, "w");
  if (file == NULL) {
    fprintf(stderr, "qplasty> Could not create file![55000]\n");
    exit(55000);
  }

  int i, j;
  for (i = 0; i < matrix_rows; i++) {
    for (j = 0; j < matrix_cols; j++) {
      fprintf(file, "%g\t", matrix[i * matrix_cols + j]);
    }
    fprintf(file, "\n");
  }
  fclose(file);

  fprintf(stderr, "qplasty> File \"%s\" written successfully\n", filename);

}

void compute_distance_matrix(int numQ, PDB *pdbQ, double *distance_matrix)
{
  int i, j;
  double dist_tmp;

  for (i = 0; i < numQ; i++) {
    for (j = i; j < numQ; j++) {
      if (i == j) {
        distance_matrix[i * numQ + j] = 0;
      } else {
        dist_tmp = cartesian_distance(pdbQ[i], pdbQ[j]);
        distance_matrix[i * numQ + j] = dist_tmp;
        distance_matrix[j * numQ + i] = dist_tmp;
      }
    }
  }

}

/**
 * For the atom index_atom determine and return its CA; If atom not part of an amino acid then just return it
 * \param index_atom the index of the interest atom
 * \param numA the count  of atoms in pdbA
 * \param pdbA the PDB structure containing the atom
 */
int get_calpha_index(int index_atom, int numA, PDB *pdbA)
{
  if (strcmp(pdbA[index_atom].type, "C") != 0 || strcmp(pdbA[index_atom].loc, "A") != 0) {
    // the atoms placed before index_atom
    int i = index_atom - 1;
    while (i >= 0 && pdbA[i].seq == pdbA[index_atom].seq) { //inside the same residue
      if (strcmp(pdbA[i].type, "C") == 0 && strcmp(pdbA[i].loc, "A") == 0)
        return i;
      i--;
    }

    //atoms after i
    i = index_atom + 1;
    while (i < numA && pdbA[i].seq == pdbA[index_atom].seq) { //inside the same residue
      if (strcmp(pdbA[i].type, "C") == 0 && strcmp(pdbA[i].loc, "A") == 0)
        return i;
      i++;
    }
  }
  return index_atom;
}

/**
 * interpolate
 * Method:
 *   Points matrix (P) defined as: [pix*I piy*I piz*I I]
 *   where pix - coordinate on x asis of point pi
 *        piy - coordinate on y asis of point pi
 *    piz - coordinate on z asis of point pi
 *    I   - identity matrix
 *  Kernel Matrix (K) is defined as a the matrix of kernel transformation on all pairewise distances
 *  The kernel transformation is a matrix 3*3 as defined in the elastic body spline.
 *  the (landscape) matrix = [K  P]
 *                           [Pt O]
 *  where Pt is the transposition of Points matrix
 *  The weight matrix (W) is the matrix of determined coeficients
 * \param numA atoms count in the structure A
 * \param pdbA structure in the original structure (to be flexed)
 * \param numB atoms count in the structure B
 * \param pdbB = pdbA; one the function is executed will contain the flexed coordinates
 * \param numQ atoms count in the structure Q
 * \param pdbQ codebook vectors of the original conformation(pdbA)
 * \param numV atoms count in the structure V
 * \param pdbV codebook vectors of the density map
 * \param interpol_method integer 1 TPS; 2 EBS  (Remark: this function can be used also for flexing by TPS)
 * \param elastic_coef float between 0.01 (soft rubber) - 0.49 (hard rubber)
 **/

void interpolate_splines(int numA, PDB *pdbA,  int numB, PDB *pdbB, int numQ, PDB *pdbQ, int numV, PDB *pdbV,
                         int interpol_method, int kernel_type, double elastic_coef, int rigid_residues)
{

  if (interpol_method == 2)
    fprintf(stderr, "qplasty>      Thin Plate Splines interpolation ");
  else
    fprintf(stderr, "qplasty>      Elastic Body Splines interpolation ");

  if (rigid_residues == 2)
    fprintf(stderr, "by atoms\n");
  else
    fprintf(stderr, "by residue\n");


  if (interpol_method == 2) {
    fprintf(stderr, "qplasty>      Kernel type: U(r)=r \n");
  } else {
    fprintf(stderr, "qplasty>      Kernel type: %d \n", kernel_type);
    fprintf(stderr, "qplasty>      Elasticity Coefficient: %g \n", elastic_coef);
  }

  int block_size               = 3;
  int blocks_number            = numV + 4;
  int matrix_size              = block_size * blocks_number;
  int weight_matrix_rows       = block_size * blocks_number;
  int weight_matrix_cols       = 1;
  int displacement_matrix_rows = block_size * blocks_number;
  int displacement_matrix_cols = 1;

  int    i, j;
  double  disp_x, disp_y, disp_z, dist;    //tmp var

  // allocate matrix memory
  double *matrix;
  matrix = (double *) alloc_vect(matrix_size * matrix_size, sizeof(double));

  double *displacement_matrix;
  displacement_matrix = (double *) alloc_vect(displacement_matrix_rows * displacement_matrix_cols, sizeof(double));

  double *weight_matrix;
  weight_matrix = (double *) alloc_vect(weight_matrix_rows * weight_matrix_cols, sizeof(double));

  double *distance_matrix;
  distance_matrix = (double *) alloc_vect(numV * numV, sizeof(double));

  compute_distance_matrix(numQ, pdbQ, distance_matrix);

  //initalize the matrices
  for (i = 0; i < numV; i++) {
    for (j = i + 1; j < numV; j++) {
      disp_x = cartesian_displacement(pdbQ[i], pdbQ[j], 'x');
      disp_y = cartesian_displacement(pdbQ[i], pdbQ[j], 'y');
      disp_z = cartesian_displacement(pdbQ[i], pdbQ[j], 'z');
      dist   = distance_matrix[i * numV + j];

      fill_matrix(matrix, matrix_size, matrix_size, i * block_size, j * block_size, block_size, block_size,
                  disp_x , disp_y, disp_z, dist, interpol_method, kernel_type, elastic_coef);
      fill_matrix(matrix, matrix_size, matrix_size, j * block_size, i * block_size, block_size, block_size,
                  disp_x , disp_y, disp_z, dist,  interpol_method, kernel_type, elastic_coef);
    }

    fill_matrix(matrix, matrix_size, matrix_size, (i)*block_size, numV * block_size, block_size, block_size,
                pdbQ[i].x , pdbQ[i].x, pdbQ[i].x, 0, 0, 0, elastic_coef);
    fill_matrix(matrix, matrix_size, matrix_size, (i)*block_size, (numV + 1)*block_size, block_size, block_size,
                pdbQ[i].y , pdbQ[i].y, pdbQ[i].y, 0, 0, 0, elastic_coef);
    fill_matrix(matrix, matrix_size, matrix_size, (i)*block_size, (numV + 2)*block_size, block_size, block_size,
                pdbQ[i].z , pdbQ[i].z, pdbQ[i].z, 0, 0, 0, elastic_coef);
    fill_matrix(matrix, matrix_size, matrix_size, (i)*block_size, (numV + 3)*block_size, block_size, block_size,
                1, 1, 1, 0, 0, 0, elastic_coef);

    fill_matrix(matrix, matrix_size, matrix_size, numV * block_size, i * block_size,  block_size, block_size,
                pdbQ[i].x , pdbQ[i].x, pdbQ[i].x, 0, 0, 0,  elastic_coef);
    fill_matrix(matrix, matrix_size, matrix_size, (numV + 1)*block_size, i * block_size,  block_size, block_size,
                pdbQ[i].y , pdbQ[i].y, pdbQ[i].y, 0, 0, 0, elastic_coef);
    fill_matrix(matrix, matrix_size, matrix_size, (numV + 2)*block_size, i * block_size,  block_size, block_size,
                pdbQ[i].z , pdbQ[i].z, pdbQ[i].z, 0, 0, 0, elastic_coef);
    fill_matrix(matrix, matrix_size, matrix_size, (numV + 3)*block_size, i * block_size,  block_size, block_size,
                1, 1, 1, 0, 0, 0, elastic_coef);

    fill_matrix(displacement_matrix, displacement_matrix_rows, displacement_matrix_cols, i * block_size, 0, block_size, 1,
                pdbV[i].x, pdbV[i].y, pdbV[i].z, 0, 0, 0, elastic_coef);
  }

  gaussjordan(matrix, matrix_size); //compute the inverse of the matrix

  //determine the coefficients
  for (i = 0; i < matrix_size; i++)
    for (j = 0; j < matrix_size; j++)
      weight_matrix[i] += matrix[i * matrix_size + j] * displacement_matrix[j];

  double *kernel_tmp_matrix;
  kernel_tmp_matrix = (double *) alloc_vect(block_size * block_size, sizeof(double));
  PDB atom;
  int index_ca;

  for (j = 0; j < numA; j++) {
    if (rigid_residues == 2) {
      atom = pdbA[j];
    } else { // the atoms move like the CA of the residue
      index_ca = get_calpha_index(j, numA, pdbA);
      atom = pdbA[index_ca];
    }

    //affine transformation
    pdbB[j].x =  weight_matrix[block_size * (numV + 3) + 0];
    pdbB[j].x += weight_matrix[block_size * (numV + 0) + 0] * atom.x +
                 weight_matrix[block_size * (numV + 1) + 0] * atom.y +
                 weight_matrix[block_size * (numV + 2) + 0] * atom.z;

    pdbB[j].y =  weight_matrix[block_size * (numV + 3) + 1];
    pdbB[j].y += weight_matrix[block_size * (numV + 0) + 1] * atom.x +
                 weight_matrix[block_size * (numV + 1) + 1] * atom.y +
                 weight_matrix[block_size * (numV + 2) + 1] * atom.z;

    pdbB[j].z =  weight_matrix[block_size * (numV + 3) + 2];
    pdbB[j].z += weight_matrix[block_size * (numV + 0) + 2] * atom.x +
                 weight_matrix[block_size * (numV + 1) + 2] * atom.y +
                 weight_matrix[block_size * (numV + 2) + 2] * atom.z;

    for (i = 0; i < numV; i++) {

      disp_x = cartesian_displacement(atom, pdbQ[i], 'x');
      disp_y = cartesian_displacement(atom, pdbQ[i], 'y');
      disp_z = cartesian_displacement(atom, pdbQ[i], 'z');
      dist = cartesian_distance(atom, pdbQ[i]);

      fill_matrix(kernel_tmp_matrix, block_size, block_size, 0, 0, block_size, block_size,
                  disp_x , disp_y, disp_z, dist, interpol_method, kernel_type, elastic_coef);

      pdbB[j].x +=
        kernel_tmp_matrix[0 + 0] * weight_matrix[i * block_size + 0] +
        kernel_tmp_matrix[0 + 1] * weight_matrix[i * block_size + 1] +
        kernel_tmp_matrix[0 + 2] * weight_matrix[i * block_size + 2];

      pdbB[j].y +=
        kernel_tmp_matrix[1 * 3 + 0] * weight_matrix[i * block_size + 0] +
        kernel_tmp_matrix[1 * 3 + 1] * weight_matrix[i * block_size + 1] +
        kernel_tmp_matrix[1 * 3 + 2] * weight_matrix[i * block_size + 2];

      pdbB[j].z +=
        kernel_tmp_matrix[2 * 3 + 0] * weight_matrix[i * block_size + 0] +
        kernel_tmp_matrix[2 * 3 + 1] * weight_matrix[i * block_size + 1] +
        kernel_tmp_matrix[2 * 3 + 2] * weight_matrix[i * block_size + 2];
    }
  }



  fprintf(stderr, "qplasty>      Interpolated %d atoms!\n", numA);
  fprintf(stderr, "qplasty> \n");

  free_vect_and_zero_ptr(&matrix);
  free_vect_and_zero_ptr(&displacement_matrix);
  free_vect_and_zero_ptr(&weight_matrix);
  free_vect_and_zero_ptr(&distance_matrix);
  free_vect_and_zero_ptr(&kernel_tmp_matrix);
}

/**
 *  Function computing the weight of a codebook vector based on the distance to the atom
 *  \param dist_atom_codebook Distance between the atom (generally CA) and codebook vector
 *  \param param indices the parameter of the used function (the power or coefficient)
 *  \param type Possible values: 1 -  1/(dist^param); 2 - 1-exp(-param/dist) more copmplicated
 **/
double get_weight(double dist_atom_codebook, int type, double weighting_pow, double dist_further_vect)
{
  //error check
  if (dist_atom_codebook == 0) {
    fprintf(stderr, "WARNING: The distance between codebook vectors should not be zero!\n");
  }

  switch (type) {
    case 1:
      return 1 / pow((double)dist_atom_codebook, (double)weighting_pow);
      break;
    case 2:
      if (dist_further_vect > dist_atom_codebook)
        return pow((double)(dist_further_vect - dist_atom_codebook) / (double)(dist_further_vect * dist_atom_codebook), (double)weighting_pow);
      else
        return 0;
      break;
    default:
      return 1 / pow((double)dist_atom_codebook, (double)weighting_pow);
  }
}

/**
 * Interpolate by Inverse Distance Weighting
 * \param numA atoms count in the structure A
 * \param pdbA structure in the original structure (to be flexed)
 * \param numB atoms count in the structure B
 * \param pdbB = pdbA; one the function is executed will contain the flexed coordinates
 * \param numQ atoms count in the structure Q
 * \param pdbQ codebook vectors of the original conformation(pdbA)
 * \param numV atoms count in the structure V
 * \param pdbV codebook vectors of the density map
 * \param weighting_type integer =1 for global inverse distance, =2 for local (requires parameter vector_influence)
 * \param weighting_pow float, recommended values: 8 for the global and 4 for the local weighting_pow
 * \param vectors_influence only used by the local weighting scheme, number of (closest) vectors used for the interpolation
 * \param rigid_residues 1 to for a rigid residue (moves as its CA) and 2 for all atoms positions are interpolated
 **/
void interpolate_IDW(int numA, PDB *pdbA,  int numB, PDB *pdbB, int numQ, PDB *pdbQ, int numV, PDB *pdbV,
                     int weighting_type, double weighting_pow, int vectors_influence, int rigid_residues)
{

  //print options
  fprintf(stderr, "qplasty>      Inverse Distance Weighting Interpolation ");

  if (rigid_residues == 2)
    fprintf(stderr, "by atoms\n");
  else
    fprintf(stderr, "by residue\n");

  if (weighting_type == 1)
    fprintf(stderr, "qplasty>      Weighting Scheme: Global \n");
  else {
    fprintf(stderr, "qplasty>      Weighting Scheme: Local \n");
    fprintf(stderr, "qplasty>      Vectors of Influence: %d \n", vectors_influence);
  }

  fprintf(stderr, "qplasty>      Weighting Exponent: %g \n", weighting_pow);

  double weights_sum, weights_sum_x, weights_sum_y, weights_sum_z, weight_tmp;
  double dist_threshold, dist;
  double dist_vec[numV];
  PDB   res_ca;
  int i, j, ii;
  int index_ca = 0;

  //printf("qplasty>      Starting the interpolation...\n");
  for (j = 0; j < numA; j++) {
    if (rigid_residues == 2) {
      res_ca = pdbA[j];
    } else {
      index_ca = get_calpha_index(j, numA, pdbA);
      res_ca = pdbA[index_ca];
    }

    weights_sum_x      = 0.0;
    weights_sum_y      = 0.0;
    weights_sum_z      = 0.0;
    weights_sum        = 0.0;

    if (weighting_type == 2) {
      //detemine distances to the nth feature vector
      for (i = 0; i < numV; i++)
        dist_vec[i] = cartesian_distance(res_ca, pdbQ[i]);

      //sort distances
      for (i = 0; i < numV - 1; i++) {
        for (ii = i + 1; ii < numV; ii++) {
          if (dist_vec[i] > dist_vec[ii]) {
            dist = dist_vec[i];
            dist_vec[i] = dist_vec[ii];
            dist_vec[ii] = dist;
          }
        }
      }

      dist_threshold = dist_vec[vectors_influence - 1];
    } else
      dist_threshold = 0;

    int atom_eq_cv = 0;
    for (i = 0; i < numV; i++) {
      dist = cartesian_distance(res_ca, pdbQ[i]);
      if (dist > 0) {
        weight_tmp = get_weight(dist, weighting_type, weighting_pow, dist_threshold);

        weights_sum_x      += pdbV[i].x * weight_tmp;
        weights_sum_y      += pdbV[i].y * weight_tmp;
        weights_sum_z      += pdbV[i].z * weight_tmp;

        weights_sum += weight_tmp;
      } else { // the atom j coincides with the codebook i thus the atom moves like the codebook
        pdbB[j].x = pdbV[i].x;
        pdbB[j].y = pdbV[i].y;
        pdbB[j].z = pdbV[i].z;
        atom_eq_cv = 1;
        i = numV;// exit
      }
    }

    if (atom_eq_cv == 0) {
      pdbB[j].x = weights_sum_x / weights_sum;
      pdbB[j].y = weights_sum_y / weights_sum;
      pdbB[j].z = weights_sum_z / weights_sum;
    }
  }
  fprintf(stderr, "qplasty>      All done!\nqplasty>\n");
}

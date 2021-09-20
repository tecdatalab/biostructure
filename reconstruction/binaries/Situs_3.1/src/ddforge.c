/*********************************************************************
*                           D D F O R G E                            *
**********************************************************************
* Program is part of the Situs package URL: situs.biomachina.org     *
* (c) Julio Kovacs and Willy Wriggers, 2018                           *
**********************************************************************
*                                                                    *
* Damped-Dynamics Flexible Fitting.                                  *
*                                                                    *
**********************************************************************
* See legal statement for terms of distribution                      *
*********************************************************************/

#include "situs.h"
#include "lib_ddf.h"
#include "lib_pio.h"
#include "lib_vio.h"
#include "lib_vec.h"

#define ptr_check(p,n){ \
    if ((p) == NULL){ \
        fprintf(stderr, "ddforge> Memory allocation failure. Exit %d.\n",n); \
        exit(n); \
    } \
}

#define max(A,B) ((A)>(B) ? (A):(B))
#define min(A,B) ((A)<(B) ? (A):(B))

#define mH 1.008
#define mC 12.011
#define mN 14.007
#define mO 15.999
#define mS 32.06

/* SCPP indicates the side-chain prediction program to use:
   1: SCATD;   2: SCWRL */
#define SCPP 2

typedef struct
{ int  rt;     /* residue type */
  int  nat;    /* # of pseudo-atoms */
  int  nan;    /* # of dihedral angles */
  int  epoc;   /* 1 if 1st res. of chain; 2 if last res. of chain; 0 otherwise */
  int  k1;     /* index of first atom */
} tri;

typedef struct
{ int     k;
  int     l;
  double  D;
} twid;

typedef struct
{ double  x;
  double  y;
  double  z;
} trd;

typedef struct
{ trd e;
  trd b;
} twovec;

/*
typedef struct
{ char    cht;
  int     rnt;
  double  x;
  double  y;
  double  z;
} Calpha;
*/

typedef struct
{ char    cht;
  int     rnt;
  double  xN;
  double  yN;
  double  zN;
  double  xCa;
  double  yCa;
  double  zCa;
  double  xC;
  double  yC;
  double  zC;
} NCaC;

typedef struct
{ char cho;
  int  rno;
  int  iko;
  char cht;
  int  rnt;
  int  ikt;
} respair;

typedef struct
{ int  code;
  char cho;
  int  rno;
  int  iko;
} fvars;

typedef struct
{ char cho;
  int  rno;
  int  iko;
} sres;

typedef struct
{ char cho;
  int ich;
  int ikf;
  int ikl;
  int ivf;
  int ivl;
} chdata;


/* functions defined in this file */
int rtype(PDB *, int);
int getqdots(PDB *pdb_model, int num_atoms, int num_res, double dc2, fvars *tfv,
             int size, int mnv, int nco1, int nco2, fvars *coords, int varlen,
             sres *statres, int srlen, respair *dconst, int dclen, int nchains,
             chdata *chinfo, tri *props, trd *wop, twovec *erb, trd *der,
             double *qdots, double damp_scale, double drag_scale);
int  calcL(int, int, int, double, PDB *, double, double [21]);
void calcV(int, int, long, double *[6], fvars *, double [21], double *);
int Wilson(int, int, int, PDB *, tri *, twovec *, fvars *, trd *);
void optpdb(int, PDB *, char *, int, char *);

static double elapsed_time(clock_t t0, clock_t t1){
  return (t1 - t0) / (double) CLOCKS_PER_SEC;
}

/* masses of pseudo-atoms 4 and 5 for each residue */
double m45[21][2];

/* Pi */
double Pi;


int main(int argc, char *argv[]){

  PDB  *pdb_model, *pdb_model_all, *pdb_target;
  PDB  *pdb_temp, *pdb0;
  tri  *props;
  twovec *erb;
  trd  *wop, *der;
  NCaC *tca;
  respair *rescor;  /* residue correspondence table */
  fvars *coords;    /* list of residues containing free variables */
  sres  *statres;   /* list of static residues */
  respair *dconst;  /* list of pairs of residues defining distance constraints */
  fvars *tfv;       /* table of free variables */
  chdata *chinfo;   /* list of chains and their data */

  double minarea, maxarea, temp, temp1, temp2, temp3, tr3;
  double vphi, vpsi, vchi;
  double cx, cy, cz, maxx, maxy, maxz, minx, miny, minz;
  double cutoff, da, dx, dy, dz, fx, fy, fz, relder;
  double correl, delang, dispcap, rmsdisp, maxdisp, mdave, maxvel;
  double damp_scale, drag_scale, force_scale, dadr_ratio;
  double dfs, dfc, mdf;                      /* scale factors */
  double Cel, Cvw, Chb, resol, frs, sigma, sigma_a, tsg12, deriv;
  double r_cut, diamstr, hi_atom, lo_atom, rq;
  double var1, var2, s0, s0_1, s0_2, alpha1, alpha2;
  double alc, blc, denom, slc, ts2, c2lc, rsd, flc;
  double gas_const, temperature, scale_factor;
  double sccx, sccy, sccz, scm, mass, scalefac;
  double ca[3], e[3], e1[3], e2[3], tv[3], tvp[3], rot[3][3];
  double p1[3], p2[3], *tempx, *tempy, *tempz, meangv, mingv;
  double *qdots, tiboc, dboc, ttime, deltat, curtime, smd, *Qs;
  double *mdarr, *spv, vv, Gs, rrc, *Gkernel, wrat;
  double *veloc;    /* velocity vector fields */
  double speedratio, rmsvel, rmsvelp, rmsvel0, spvv, cavv, rho, infty;
  double box, boy, boz, q0, q1, q2, thresint, thresval;
  double rc1[3], rc2[3], rc3[3], rc4[3], delx, dely, delz;
  double rcond, nclc, cosvw, norw, residual, nearends, dist_cut, dcfinal, dc2;
  double g_orix, g_oriy, g_oriz, xp, yp, zp, tz2, ty2;
  double sumf, sumg, sumg0, sumgold, sumd, sumdp, wpdbtar, Kforce;
  double cnorm, Fnorm, F0x, F0y, F0z, coef, alpha, alo2, a0, a1, a2, b1, ffdc, ffdc2;
  double g_width, *g_phi_lo, *g_phi_lo_sorted, *g_phi_hi, *g_phi_sorted;
  double g_lo_cutoff, g_hi_cutoff, isoval1, isoval2, sumgs, Lcapr;
  double intGk, summasses, essmaxlo, Mcurr, maxg, quant_lo;
  double intomax, discrep, discrep_a, incr, maxincr, epsilon;
  double deltadens, lodens, hidens, meand, mapforce2, coefb;
  double phiaaa, phiaab, phiaba, phiabb, phibaa, phibab, phibba, phibbb;
  double timing, xmin, xmax, ymin, ymax, zmin, zmax;
  double S0a, S1a, S2a, S3a, ber, cer, ker, sum1, sum2, deltak, f1;
  double bp, cp, kp, error, error1, ytime, rtime;
  int c, i, i1, i2, ii, iim, ip, j, jj, jjm, k, k1, k2, l, lp, m, n, p;
  int ks, kN, kC, kCa, kCb, kOx, ich, root_index, ffsl, argindex, mer, n2;
  int num_atoms, num_atoms_bb, num_atoms_all, num_atoms_target, Ngv, tna, natheor;
  int num_res, resnum, rno, R, nchains, mnv, code, flagwr, npt, nmdv;
  int flag, altmode, neig, istep, nsslr, nsslp, iconf, nconfs, vlt;
  int atc, direction, ifr, ifpa, num_atoms_loop, nco1, nco2, flagscopt, flagsc0, flagsc1;
  int *iwork, size, nu1, nu2, isnc, icol, nbins, sgnc, maxtype;
  int atin[6], fnz, nf1, info, iend, itemp, istat, imod, Nmod, maxiters;
  int nrows, nrhs, lda, ldb, nlvl, lwork, liwork, rank, heflag;
  int iseq, iseq1, iseq2, jrl, jrr, nsboc, mpt, coway, cosgn, iter;
  int nrestar, corlen, rnori, rntar, varlen, srlen, rdflag, dcflag, dclen, Nplus;
  int g_extx, g_exty, g_extz, ix, iy, iz, Gext, Gext2;
  int ixa, iya, iza, ixb, iyb, izb, ixg, iyg, izg;
  long index, indzy, exey, ind1, ind2, indza, indzb, indya, indyb, nipa;
  long VolC, numgs, fftype, fpower, fps, fpt, *indexes, *indexes2, indexc;
  long eeiz, eeizg, indzyg, indatom, nvoxlayer, ll;
  unsigned long m0, m1, m2, mq, g_nvox, L;
  char file_name[100], bname[100], line[101], chID[2], chori[2], chtar[2];
  char atom, resname[5], trans, rdfile[200], dcfile[200], chainid;
  char *exten;
  FILE *out, *mf, *fca;
  clock_t t0, t1, tb0, tb1;
  char usage_string[] = "ddforge>\n USAGE:\n   ddforge <pdb> <map> <min. disp. between outputs (A)> <map resolution (A)> <density threshold> <damp/drag ratio> <side-chain optim. (0|1|2|3)> <force-field cutoff distance (A)> <max. disp. per step (A)> [-rdfile <residue-data file name>] [-dcfile <distance-constraint file name>]\n --OR--\n   ddforge <pdb1> <pdb2> <min. disp. between outputs (A)> <damp/drag ratio> <side-chain optim. (0|1|2|3)> <force-field cutoff distance (A)> <max. disp. per step (A)> [-rdfile <residue-data file name>] [-dcfile <distance-constraint file name>]\n --OR--\n   ddforge <pdb1> <pdb2> <map> <min. disp. between outputs (A)> <map resolution (A)> <density threshold> <damp/drag ratio> <side-chain optim. (0|1|2|3)> <force-field cutoff distance (A)> <max. disp. per step (A)> [-rdfile <residue-data file name>] [-dcfile <distance-constraint file name>]\n";

  /* OLD usage_string:
   "ddforge>\n USAGE:\n   ddforge <pdb> <map> <residue data filename> <min. disp. between outputs (A)> <map resolution (A)> <density threshold> <damp/drag ratio> <side-chain optim. (0|1|2|3)> <force-field cutoff distance (A)> <max. disp. per step (A)> [-dcfile <distance-constraint file name>]\n --OR--\n   ddforge <pdb1> <pdb2> <residue data filename> <min. disp. between outputs (A)> <damp/drag ratio> <side-chain optim. (0|1|2|3)> <force-field cutoff distance (A)> <max. disp. per step (A)> [-dcfile <distance-constraint file name>]\n --OR--\n   ddforge <pdb1> <pdb2> <map> <residue data filename> <min. disp. between outputs (A)> <map resolution (A)> <density threshold> <damp/drag ratio> <side-chain optim. (0|1|2|3)> <force-field cutoff distance (A)> <max. disp. per step (A)> [-dcfile <distance-constraint file name>]\n"
  */

  if(argc < 8){
    fprintf(stderr, usage_string);
    fprintf(stderr,"ddforge>   side-chain optimization flag: 0-none, 1-initial only, 2-subsequent only, 3-both\n\n");
    exit(0);
  }

  mpt=0;

  exten = strrchr(argv[1],'.');
  if(exten==NULL){
    fprintf(stderr,"ddforge> No extension found for first argument --must be '.pdb'\n");
    exit(0);
  }
  exten++;
  if(strcmp(exten,"pdb")!=0){
    fprintf(stderr,"ddforge> Extension of first argument must be '.pdb'\n");
    exit(0);
  }

  exten = strrchr(argv[2],'.');
  if(exten==NULL){
    fprintf(stderr,"ddforge> No extension found for second argument --must be '.sit(us)' or '.mrc' or '.pdb'\n");
    exit(0);
  }
  exten++;
  if(have_situs_suffix(argv[2]) || strcmp(exten,"mrc")==0)     // if(strcmp(exten,"sit")==0)
    mpt = 1;
  else if(strcmp(exten,"pdb")==0)
    mpt = 2;
  else{
    fprintf(stderr,"ddforge> Extension of second argument not valid --must be '.sit(us)' or '.mrc' or '.pdb'\n");
    exit(0);
  }

  exten = strrchr(argv[3],'.');
  if (exten != NULL){
    exten++;
    if(have_situs_suffix(argv[3]) || strcmp(exten,"mrc")==0){   // if(strcmp(exten,"sit")==0)
      if(mpt==2)
        mpt = 3;
      else{
        fprintf(stderr,"ddforge> Two map files cannot be specified.\n");
        exit(0);
      }
    }
  }

  if( !(argc >= 10 && mpt==1) && !(argc >= 8 && mpt==2) && !(argc >= 11 && mpt==3) ){
    fprintf(stderr, usage_string);
    fprintf(stderr,"ddforge>   side-chain optimization flag: 0-none, 1-initial only, 2-subsequent only, 3-both\n\n");
    exit(0);
  }


  if(mpt==1){
    sscanf(argv[3],"%le", &dboc);          /* min displacement between output conformations (A) */
    sscanf(argv[4],"%le", &resol);         /* map resolution (A) */
    sscanf(argv[5],"%le", &g_lo_cutoff);   /* density values below this will be set to 0,
                                              and all others will be lowered by this amount */
    sscanf(argv[6],"%le", &dadr_ratio);    /* damping_scale/drag_scale */
    sscanf(argv[7],"%d",  &flagscopt);     /* flag for side-chain optimization (see below) */
    sscanf(argv[8],"%le", &ffdc);          /* force-field distance cutoff (A) */
    sscanf(argv[9],"%le", &maxdisp);       /* max displacement per time step (A) */
  }

  if(mpt==2){
    sscanf(argv[3],"%le", &dboc);          /* min displacement between output conformations (A) */
    sscanf(argv[4],"%le", &dadr_ratio);    /* damping_scale/drag_scale */
    sscanf(argv[5],"%d",  &flagscopt);     /* flag for side-chain optimization (see below) */
    sscanf(argv[6],"%le", &ffdc);          /* force-field distance cutoff (A) */
    sscanf(argv[7],"%le", &maxdisp);       /* max displacement per time step (A) */
  }

  if(mpt==3){
    sscanf(argv[4],"%le", &dboc);          /* min displacement between output conformations (A) */
    sscanf(argv[5],"%le", &resol);         /* map resolution (A) */
    sscanf(argv[6],"%le", &g_lo_cutoff);   /* density values below this will be set to 0,
                                              and all others will be lowered by this amount */
    sscanf(argv[7],"%le", &dadr_ratio);    /* damping_scale/drag_scale */
    sscanf(argv[8],"%d",  &flagscopt);     /* flag for side-chain optimization (see below) */
    sscanf(argv[9],"%le", &ffdc);          /* force-field distance cutoff (A) */
    sscanf(argv[10],"%le",&maxdisp);       /* max displacement per time step (A) */
  }

  /* flagscopt defines the side-chain optimization, according to:
       0: none
       1: initial only
       2: subsequent only
       3: both initial and subsequent
     This value is then separated into two flags: flagsc0 and flagsc1, each meaning:
       flagsc=1 activates the periodic side-chain prediction with SCATD or SCWRL
       flagsc=0 skips side-chain prediction */
  if(flagscopt==0){
    flagsc0 = 0;
    flagsc1 = 0;
  }
  else if(flagscopt==1){
    flagsc0 = 1;
    flagsc1 = 0;
  }
  else if(flagscopt==2){
    flagsc0 = 0;
    flagsc1 = 1;
  }
  else if(flagscopt==3){
    flagsc0 = 1;
    flagsc1 = 1;
  }
  else{
    fprintf(stderr,"ddforge> Side-chain optimization flag not valid. Must be 0, 1, 2 or 3.\n");
    exit(0);
  }

  printf("ddforge>\n");
  if(mpt==1){
    printf("ddforge> Atomic structure file: %s\n", argv[1]);
    printf("ddforge> Density map file:      %s\n", argv[2]);
    printf("ddforge> Min. displacement between output confs: %f\n", dboc);
    printf("ddforge> Map resolution: %f A\n", resol);
    printf("ddforge> Density threshold: %f\n", g_lo_cutoff);
    printf("ddforge> Damp/drag ratio: %f\n", dadr_ratio);
    printf("ddforge> Side-chain optimization flag: %d\n", flagscopt);
    printf("ddforge> Force-field distance cutoff: %f A\n", ffdc);
    printf("ddforge> Max. displacement per time step: %f A\n", maxdisp);
  }
  if(mpt==2){
    printf("ddforge> Starting structure file: %s\n", argv[1]);
    printf("ddforge> Target structure file:   %s\n", argv[2]);
    printf("ddforge> Min. displacement between output confs: %f\n", dboc);
    printf("ddforge> Damp/drag ratio: %f\n", dadr_ratio);
    printf("ddforge> Side-chain optimization flag: %d\n", flagscopt);
    printf("ddforge> Force-field distance cutoff: %f A\n", ffdc);
    printf("ddforge> Max. displacement per time step: %f A\n", maxdisp);
  }
  if(mpt==3){
    printf("ddforge> Starting structure file: %s\n", argv[1]);
    printf("ddforge> Target structure file:   %s\n", argv[2]);
    printf("ddforge> Density map file:        %s\n", argv[3]);
    printf("ddforge> Min. displacement between output confs: %f\n", dboc);
    printf("ddforge> Map resolution: %f A\n", resol);
    printf("ddforge> Density threshold: %f\n", g_lo_cutoff);
    printf("ddforge> Damp/drag ratio: %f\n", dadr_ratio);
    printf("ddforge> Side-chain optimization flag: %d\n", flagscopt);
    printf("ddforge> Force-field distance cutoff: %f A\n", ffdc);
    printf("ddforge> Max. displacement per time step: %f A\n", maxdisp);
  }

  Pi = 3.14159265;

  force_scale = 1.0;
  drag_scale  = 1.0;
  damp_scale  = drag_scale*dadr_ratio;

  wpdbtar=0.5;  /* between 0 and 1; used only if mpt=3 as the weight for the
                   pdb-produced force, and 1-wpdbtar will be the weight for the
                   map-produced force. This applies only for residues that have
                   a corresponding target residue; those that don't will be pulled
                   by the map independently of this weight. */

  Kforce = 1.0; /* desired ratio between the pdb force and map force, in the case mpt=3 */

  /* ffsl=1 defines x.A(x) to go as sqrt(x) for small x. This can cause instabilities
            especially when using a pdb as the target. */
  /* ffsl=2 defines x.A(x) to go as x for small x. This option is more stable when
            the target is an atomic structure. */
  ffsl = 2;

  /* parameters of the function A(x) in the force-field formula (OLD 1): */
  /* a0 = 1.5;                   desired position of the max of x.A(x) */
  /* alpha = 2.0;                exponent of x */
  /* alo2 = alpha/2.0; */
  /* a1 = 2.0*(alpha-1)*pow(a0,alpha-0.5);  coef of sqrt(x) */

  /* parameters of the function A(x) (OLD 2): */
  /* a0 = 1.5;                desired position of the max of x.A(x)
   if(ffsl==1)
   a1 = pow(a0,1.5);         coef of sqrt(x)
   if(ffsl==2)
   a1 = pow(a0,2);           coef of 1
   a2 = a0*a0; */

  /* parameters of the function A(x) (NEW): */
  a0 = 1.5;                   /* desired position of the max of x.A(x) */
  a1 = sqrt(a0);
  a2 = pow(a0,2);
  b1 = pow(a0,fpower);
  ffdc2 = pow(ffdc,2);        /* square of force field cutoff distance */
  /* fpower: exponent of (a0/r) to use in the force field calculation; must a positive integer: */
  fpower = 3;        /* 3 for 1/distance^2 force;  2 for 1/distance force */
  fps = fpower/2;    /* quotient */
  fpt = fpower % 2;  /* remainder */

  /* factor to scale down dist_cut (no longer used): */
  dfs = 0.5;

  /* if heflag=1, do histogram equalization of the maps: */
  heflag = 0;

  /* index of root to be applied to the density maps: */
  /* root_index = 1; */

  /* initial rms displacement at each time step (A): */ 
  rmsdisp = 0.5 * dboc;

  /* initial displacement cap: */
  dispcap = maxdisp;  /* instead of min(maxdisp, dboc); (2018/2/15) */

  /* npt+1 = number of points to compare directions of velocities: */
  npt = 2;

  /* nmdv+1 = number of points to average rmsdisp: */
  nmdv = 5;

  if(mpt==1 || mpt==3){
    /* threshold on the integral: */
    thresint = 0.001;

    /* force-field type:
       1: density attracts all atoms.
       2: density attracts atoms of similar values.
       3: density attracts atoms of similar values and
          similar gradient directions (not yet implemented).
    */
    fftype = 1;
    
    nbins = 100;    /* number of intervals to divide the density range or
                       the number of voxels in the map (for fftype=2 or 3) */
    
    /* ratio between cap on g_hi_cutoff and the max of g_phi_hi: */
    Lcapr = 1.0;

    /* way to find the threshold value for g_phi_hi: */
    /* 1: equalizing volumes at mean (or threshold) density (according to 'vlt');
       2: equalizing integrals and maxima (or quantiles, according to 'maxtype');
       3: simultaneous optimization of threshold and sigma,
          by equalizing volumes, integrals, and maxima (or quantiles)
     */
    coway = 3;

    /* choose where to measure the volume: */
    vlt = 1;

    /* max number of iterations to solve for b,c,k in the exponential regression
       for the stopping criterion: */
    maxiters = 5000;

    /* relative error tolerance in the iterative solution for b,c,k: */
    epsilon = 0.00000001;

    /* saturation values of the overlap function: */
    alpha1 = 0.9;
    alpha2 = 0.99;

    /* array for the overlap values along the trajectory, for use in the regression: */
    Qs = (double *) malloc(1*sizeof(double));
  }

  /* initialize m45 (masses of pseudo-atoms 4 and 5 for each residue type) */
  m45[ 0][0]= 14.00;   m45[ 0][1]= 46.0;   /* approx. mean of the 20 types */
  m45[ 1][0]= mC+3*mH; m45[ 1][1]= 0.0;
  m45[ 2][0]= mC+2*mH; m45[ 2][1]= mS+mH;
  m45[ 3][0]= mC+2*mH; m45[ 3][1]= mC+2*mO;
  m45[ 4][0]= mC+2*mH; m45[ 4][1]= 2*mC+2*mH+2*mO;
  m45[ 5][0]= mC+2*mH; m45[ 5][1]= 6*mC+5*mH;
  m45[ 6][0]= 0.0;     m45[ 6][1]= 0.0;
  m45[ 7][0]= mC+2*mH; m45[ 7][1]= 3*mC+2*mN+3*mH;
  m45[ 8][0]= mC+  mH; m45[ 8][1]= 3*mC+8*mH;
  m45[ 9][0]= mC+2*mH; m45[ 9][1]= 3*mC+mN+9*mH;
  m45[10][0]= mC+2*mH; m45[10][1]= 3*mC+7*mH;
  m45[11][0]= mC+2*mH; m45[11][1]= 2*mC+mS+5*mH;
  m45[12][0]= mC+2*mH; m45[12][1]= mC+mN+mO+2*mH;
  m45[13][0]= mC+2*mH; m45[13][1]= 2*mC+4*mH;
  m45[14][0]= mC+2*mH; m45[14][1]= 2*mC+mN+mO+4*mH;
  m45[15][0]= mC+2*mH; m45[15][1]= 3*mC+3*mN+9*mH;
  m45[16][0]= mC+2*mH; m45[16][1]= mO+mH;
  m45[17][0]= mC+  mH; m45[17][1]= mC+mO+4*mH;
  m45[18][0]= mC+  mH; m45[18][1]= 2*mC+6*mH;
  m45[19][0]= mC+2*mH; m45[19][1]= 8*mC+mN+6*mH;
  m45[20][0]= mC+2*mH; m45[20][1]= 6*mC+mO+5*mH;


  /* start measuring running time */
  tb0 = clock();


  /* read (starting) atomic structure */
  printf("ddforge>\nddforge> Reading %s\n", argv[1]);
  read_pdb(argv[1], &num_atoms_all, &pdb_model_all);

  strcpy(bname,argv[1]);
  exten = strrchr(bname,'.');
  *exten = '\0';  /* now bname is the name of the (starting) pdb file, without the extension ".pdb" */

  /* initial SCATD/SCWRL to make sure we use a complete atomic structure;
     the side-chain-optimized structure is written in a file with extension "deform.0.pdb" */
  if(flagsc0==1){
    printf("ddforge>\nddforge> Performing initial side-chain prediction\n");
    optpdb(num_atoms_all, pdb_model_all, bname, 0, file_name);
    free(pdb_model_all);
    read_pdb(file_name, &num_atoms_all, &pdb_model_all);
  }

  if(mpt==1 || mpt==3){
    /* read density map */
    if(mpt==1) argindex=2;
    if(mpt==3) argindex=3;
    printf("ddforge>\nddforge> Reading %s\n", argv[argindex]);
    read_vol(argv[argindex], &g_width, &g_orix, &g_oriy, &g_oriz,
            &g_extx, &g_exty, &g_extz, &g_phi_lo);
    g_nvox = g_extx * g_exty * g_extz;

    /* set density values below g_lo_cutoff to zero, and lower the rest of the map */
    printf("ddforge> Setting density values below %f to zero.\n", g_lo_cutoff);
    for (m=0; m<g_nvox; m++){
      if (g_phi_lo[m]<g_lo_cutoff)
        g_phi_lo[m] = 0.0;
      else
        g_phi_lo[m] -= g_lo_cutoff;
    }

/* there might be a better place for this hist.eq. - must be done consistently with
   that of g_phi_hi */
    if (heflag==1){
      /* replace g_phi_lo by its histogram-equalized version */
      hist_eq(g_phi_lo, g_nvox);
    }

    /* root of g_phi_lo */
    /* if(root_index != 1)
      for(m=0;m<g_nvox;m++)
        g_phi_lo[m] = pow(g_phi_lo[m], 1.0/root_index);
    */

    /* compute the integral of g_phi_lo */
    sumf = 0.0;
    for(m=0;m<g_nvox;m++)
      sumf += g_phi_lo[m];

    if(sumf==0.0){
      printf("After cut off, g_phi_lo is everywhere 0. Stopping.\n");
      fprintf(stderr,"After cut off, g_phi_lo is everywhere 0. Stopping.\n");
      exit(1);
    }

    /* normalize so that the integral becomes 1 */
    for(m=0;m<g_nvox;m++)
      g_phi_lo[m] /= sumf;

    sumf = 1.0;     /* new integral of g_phi_lo */

    /* allocate storage for the structure-induced map: */
    g_phi_hi = (double *) malloc(g_nvox * sizeof(double));
    ptr_check(g_phi_hi,2);
  
    /* allocate storage for sorted arrays: */
    g_phi_lo_sorted = (double *) malloc(g_nvox * sizeof(double));
    ptr_check(g_phi_lo_sorted,2);
  
    g_phi_sorted = (double *) malloc(g_nvox * sizeof(double));
    ptr_check(g_phi_sorted,2);
  
    /* allocate storage for 'indexes' and 'indexes2' arrays: */
    indexes = (unsigned long *) malloc(g_nvox * sizeof(long));
    ptr_check(indexes,2);
    indexes2 = (unsigned long *) malloc(g_nvox * sizeof(long));
    ptr_check(indexes2,2);

    /* index table for g_phi_lo: */
    indexy(g_nvox, g_phi_lo, indexes);

    /* apply index table (permutation) to sort g_phi_lo: */
    for (m=0; m<g_nvox; m++)
      g_phi_lo_sorted[m] = g_phi_lo[indexes[m]];

    /* compute 'essential max' of g_phi_lo: */
    sumgs = 0.0; numgs = 0;
    for(index=g_nvox-1; index>=0; index--){
      Mcurr = g_phi_lo_sorted[index];
      sumgs += Mcurr;
      numgs++;
      if(sumgs-Mcurr*numgs >= thresint) break;
    }
    essmaxlo = Mcurr;

    /* compute 'rq' quantile of values of g_phi_lo >0 */
    rq = 0.95;     /* 95% */
    for (m=0; m<g_nvox; m++)
      if (g_phi_lo_sorted[m]>0) break;
    m0 = m;         /* index of 1st positive value */
    L = g_nvox-m0;  /* number of positive values */
    /* index of element that provides the requested quantile: */
    mq = min(max(ceil(rq*L),1), L) + m0 - 1;
    quant_lo = g_phi_lo_sorted[mq];

    /* mean and min of values of g_phi_lo greater than 0: */
    meangv=0.0; Ngv=0; mingv=g_phi_lo_sorted[g_nvox-1];
    for(m=0; m<g_nvox; m++)
      if(g_phi_lo[m]>0.0){
        if(g_phi_lo[m]<mingv)
          mingv = g_phi_lo[m];
        meangv += g_phi_lo[m];
        Ngv++;
      }
    meangv /= Ngv;
  
    /* threshold on density values: */
    if(Ngv>0)
      thresval = thresint/Ngv;
    else
      thresval = 0.0;
  
    /* interval of density for each "attracting layer" (fftype=2 or 3): */
    deltadens = essmaxlo/nbins;

    /* number of voxels for each "attracting layer" (fftype=2 or 3): */
    nvoxlayer = Ngv/nbins;

    /* isocontour value where to measure the volume: */
    if (vlt==1) isoval1 = 0.0;     /* at threshold */
    if (vlt==2) isoval1 = meangv;  /* at mean density */

    /* volume inside the isosurface of value 'isoval1': */
    VolC=0;
    for(m=0; m<g_nvox; m++)
      if(g_phi_lo[m] > isoval1)
        VolC++;

    if(VolC==0){
      printf("Vol(C)=0. Stopping.\n");
      fprintf(stderr,"Vol(C)=0. Stopping.\n");
      exit(1);
    }
  }

  if(mpt==2 || mpt==3){
    /* read target atomic structure */
    printf("ddforge>\nddforge> Reading %s\n", argv[2]);
    read_pdb(argv[2], &num_atoms_target, &pdb_target);

    /* list of N, Calpha, C atoms of target PDB, and their coordinates */
    tca = (NCaC *) malloc(1*sizeof(NCaC));

    j=0;    /* residue index */
    for(i=0;i<num_atoms_target;i++){
      if( strcmp(pdb_target[i].type, "C")==0 && strcmp(pdb_target[i].loc, "A")==0 ){
        atc=1;
        tca = (NCaC *) realloc(tca, (j+1)*sizeof(NCaC));
        tca[j].cht = pdb_target[i].chain[0];
        tca[j].rnt = pdb_target[i].seq;
        tca[j].xCa = pdb_target[i].x;
        tca[j].yCa = pdb_target[i].y;
        tca[j].zCa = pdb_target[i].z;

        /* search for first atom of this residue: */
        for(k=i; k>=0 && pdb_target[k].seq == tca[j].rnt
                      && pdb_target[k].chain[0] == tca[j].cht; k--){;}
        k1=k+1;

        /* search for last atom of this residue: */
        for(k=i; k<num_atoms_target && pdb_target[k].seq==tca[j].rnt
                                    && pdb_target[k].chain[0]==tca[j].cht; k++){;}
        k2=k-1;

        /* search for N and C atoms of this residue, and store their coords: */
        for(k=k1; atc<3 && k<=k2; k++){
          if(strcmp(pdb_target[k].type, "N")==0 && strcmp(pdb_target[k].loc, "" )==0){
            atc++; 
            tca[j].xN = pdb_target[k].x;
            tca[j].yN = pdb_target[k].y;
            tca[j].zN = pdb_target[k].z;
            continue;
          }
          if(strcmp(pdb_target[k].type, "C")==0 && strcmp(pdb_target[k].loc, "" )==0){
            atc++;
            tca[j].xC = pdb_target[k].x;
            tca[j].yC = pdb_target[k].y;
            tca[j].zC = pdb_target[k].z;
            continue;
          }
        }
        
        j++; i=k2;
      }
    }
    nrestar=j;    /* number of residues of target pdb */
  }


  /* save original PDB structure */
  pdb0 =(PDB *) malloc(num_atoms_all*sizeof(PDB));
  ptr_check(pdb0,3);
  for(index=0;index<num_atoms_all;index++)
    pdb0[index] = pdb_model_all[index];

  pdb_model = (PDB *) malloc(1*sizeof(PDB));
  ptr_check(pdb_model,3);

  pdb_temp = (PDB *) malloc(5*sizeof(PDB));
  ptr_check(pdb_temp,4);

  props = (tri *) malloc(1*sizeof(tri));
  ptr_check(props,5);

  /* residue correspondence table: */
  rescor  = (respair *) malloc(1*sizeof(respair));
  ptr_check(rescor,66);

  coords  = (fvars *) malloc(1*sizeof(fvars));
  ptr_check(coords,6);

  /* static residues constraints */
  statres = (sres *)  malloc(1*sizeof(sres));
  ptr_check(statres,7);

  /* distance constraints */
  dconst  = (respair *) malloc(1*sizeof(respair));
  ptr_check(rescor,77);


  /* sum of masses of atomic structure: */
  summasses = 0.0;
  for(k=0; k<num_atoms_all; k++){
    /* mass of this atom */
    atom = pdb_model_all[k].type[strlen(pdb_model_all[k].type)-1];
    if(atom=='C')
      mass=mC;
    else if(atom=='H')
      continue;                  /* H atoms don't generate density */
    else if(atom=='O')
      mass=mO;
    else if(atom=='N')
      mass=mN;
    else if(atom=='S')
      mass=mS;
    else
      continue;                  /* unknown atom */
    summasses += mass;
  }


  /* test of indexy
  double aa[] = {-3.3, 4., 2., 9., 17., 43., 12., 0., -3., 4.4, 2.5, 9.5};
  // index table (permutation) for aa:
  indexy(12, aa, indexes2);
  // apply permutation to print a sorted version of aa:
  for (m=0; m<12; m++)
    printf("%f, ", aa[indexes2[m]]);
  printf("\n");
  exit(0);    */


  if(mpt==1 || mpt==3){

  /** determine optimal threshold and sigma for g_phi_hi **/

    /* choose between essmax or quantile: */
    maxtype = 2;

    /* ratio between integral and max of g_phi_lo: */
    if (maxtype==1)  /* use essmax */
      intomax = sumf/essmaxlo;
    if (maxtype==2)  /* use quantile */
      intomax = sumf/quant_lo;

    /* initial sigma: */
    sigma = 0.5/sqrt(3.0) * resol;

    iter = 0;   /* iteration  number */

    fprintf(stderr, "\nddforge> Optimization of sigma and threshold:\n");

    do{
      if (iter>0) free(Gkernel);

      /* parameters for Gaussian kernel */
      r_cut  = 3.3*sigma;              /* cutoff radius */
      Gs   = 0.5;                      /* grid spacing (in A) */
      Gext = 2*ceil(r_cut/Gs)+1;       /* extent of grid */
      Gext2 = Gext*Gext;
      rrc  = 0.5*(Gext-1)*Gs;          /* rounded r_cut */

      Gkernel = (double *) malloc(Gext2*Gext * sizeof(double));
      ptr_check(Gkernel, 1);

      tsg12 = 2.0*pow(sigma,2);

      /* constants for blurring the atomic structure: */
      wrat = g_width/Gs;
      box  = g_orix/Gs + 0.5*Gext;
      boy  = g_oriy/Gs + 0.5*Gext;
      boz  = g_oriz/Gs + 0.5*Gext;

      /* compute Gaussian kernel and its integral */
      intGk=0.0;
      for(iz=0;iz<Gext;iz++){
        zp = (iz-(Gext-1)/2)*Gs;
        eeiz = Gext2*iz;
        tz2 = pow(zp, 2);
        for(iy=0;iy<Gext;iy++){
          yp = (iy-(Gext-1)/2)*Gs;
          indzy = eeiz + Gext*iy;
          ty2 = pow(yp, 2);
          for(ix=0;ix<Gext;ix++){
            xp = (ix-(Gext-1)/2)*Gs;
            temp = pow(xp, 2) + ty2 + tz2;
            Gkernel[indzy+ix] = exp(-temp/tsg12);   /* kernel */
            intGk += Gkernel[indzy+ix];             /* integral */
          }
        }
      }

      /* correction to the integral of the kernel
       to account for the difference in grid spacings: */
      intGk *= pow(Gs/g_width,3);

      /* integral of structure-induced map (before thresholding) */
      sumg0 = summasses*intGk;

      /* compute g_phi_hi using current value of sigma */
      exey = g_extx * g_exty;
      for(m=0; m<g_nvox; m++)
        g_phi_hi[m] = 0.0;

      for(k=0; k<num_atoms_all; k++){
        atom = pdb_model_all[k].type[strlen(pdb_model_all[k].type)-1];
        if(atom=='C')
          mass=mC;
        else if(atom=='H')
          continue;                  /* H atoms don't generate density */
        else if(atom=='O')
          mass=mO;
        else if(atom=='N')
          mass=mN;
        else if(atom=='S')
          mass=mS;
        else
          continue;                  /* unknown atom */

        p1[0] = pdb_model_all[k].x;
        p1[1] = pdb_model_all[k].y;
        p1[2] = pdb_model_all[k].z;
        ixa = ceil((p1[0]-g_orix-rrc)/g_width);
        iya = ceil((p1[1]-g_oriy-rrc)/g_width);
        iza = ceil((p1[2]-g_oriz-rrc)/g_width);
        ixb = floor((p1[0]-g_orix+rrc)/g_width);
        iyb = floor((p1[1]-g_oriy+rrc)/g_width);
        izb = floor((p1[2]-g_oriz+rrc)/g_width);
        if(ixa>=g_extx || iya>=g_exty || iza>=g_extz) continue;
        if(ixb<0 || iyb<0 || izb<0) continue;
        if(ixa<0) ixa=0; if(iya<0) iya=0; if(iza<0) iza=0;
        if(ixb>=g_extx) ixb=g_extx-1;
        if(iyb>=g_exty) iyb=g_exty-1;
        if(izb>=g_extz) izb=g_extz-1;

        q0 = p1[0]/Gs; q1 = p1[1]/Gs; q2 = p1[2]/Gs;

        for(iz=iza;iz<=izb;iz++){
          eeiz = exey*iz;
          izg = floor(iz*wrat+boz-q2);
          eeizg = Gext2*izg;
          for(iy=iya;iy<=iyb;iy++){
            indzy = eeiz + g_extx*iy;
            iyg = floor(iy*wrat+boy-q1);
            indzyg = eeizg + Gext*iyg;
            for(ix=ixa;ix<=ixb;ix++){
              ixg = floor(ix*wrat+box-q0);
              g_phi_hi[indzy+ix] += mass * Gkernel[indzyg+ixg];
            }
          }
        }
      }

      /* for testing: 
      sprintf(file_name, "test_g_phi_hi_%d.sit", iter);
      write_vol(file_name, g_width, g_orix, g_oriy, g_oriz, g_extx, g_exty, g_extz, g_phi_hi);
      */

      /* old way to sort g_hi_hi:
       for(m=0; m<g_nvox; m++) g_phi_sorted[m] = g_phi_hi[m];
       quicksort(g_nvox, g_phi_sorted);
       */

      /* index table for g_phi_hi: */
      indexy(g_nvox, g_phi_hi, indexes2);

      /* apply index table (permutation) to sort g_phi_hi: */
      for (m=0; m<g_nvox; m++)
        g_phi_sorted[m] = g_phi_hi[indexes2[m]];

      /* isosurface value for g_phi_hi whose inside volume is equal to VolC: */
      isoval2 = g_phi_sorted[g_nvox-VolC];

      /* --- determine threshold using current sigma: --- */
      sumgs = 0.0; numgs = 0;

      if (coway==1 || coway==3){ /* volume-based or simultaneous volume&integral&max */
        if (vlt==1){    /* volume at threshold */
          /* sumgs,numgs are needed in the calc of discrepenacy below */
          for(index=g_nvox-1; index>=0; index--){
            Mcurr = g_phi_sorted[index];
            if(Mcurr > isoval2){
              sumgs += Mcurr; numgs++;
            }
            else break;
          }
          g_hi_cutoff = isoval2;
        }
        if (vlt==2){    /* volume at mean density */
          for(index=g_nvox-1; index>=0; index--){
            Mcurr = g_phi_sorted[index];
            if(Mcurr>0.0){
              sumgs += Mcurr; numgs++;
            }
            else break;
            if(sumgs/numgs <= isoval2) break;
          }
          g_hi_cutoff = g_phi_sorted[max(index,0)];
        }
        /* compute max (or quantile) of g_phi_hi: */
        if (maxtype==1)
          maxg = g_phi_sorted[g_nvox-1];
        if (maxtype==2){
          /* compute 'rq' quantile of values of g_phi_hi > g_hi_cutoff */
          for (m=0; m<g_nvox; m++)
            if (g_phi_sorted[m] > g_hi_cutoff) break;
          m0 = m;         /* index of 1st value >threshold */
          L = g_nvox-m0;  /* number of values >threshold */
          /* index of element that provides the requested quantile: */
          mq = min(max(ceil(rq*L),1), L) + m0 - 1;
          maxg = g_phi_sorted[mq];
        }
      }

      if (coway==2){  /* based on integral&max */
        if (maxtype==1)    /* max of g_phi_hi */
          maxg = g_phi_sorted[g_nvox-1];

        for (index=g_nvox-1; index>=0; index--){
          Mcurr = g_phi_sorted[index];
          if (Mcurr > 0.0){
            sumgs += Mcurr; numgs++;
          }
          else break;

          if (maxtype==2){  /* 'rq' quantile of values of g_phi_hi > Mcurr */
            for (m=index; m<g_nvox; m++)
              if (g_phi_sorted[m] > Mcurr) break;
            m0 = m;         /* index of 1st value >Mcurr */
            L = g_nvox-m0;  /* number of values >Mcurr */
            /* index of element that provides the requested quantile: */
            mq = min(max(ceil(rq*L),1), L) + m0 - 1;
            maxg = g_phi_sorted[mq];
          }

          if (Mcurr >= maxg) continue;
          if ((sumgs-Mcurr*numgs)/(maxg-Mcurr) >= intomax) break;
        }

        g_hi_cutoff = Mcurr;
      }

      /* test (for GroEL):
      g_hi_cutoff = 0.47*maxg; */

      fprintf(stderr, "ddforge> iter = %2d,  sigma = %5.2f,  threshold_g = %8.2e,  max_g = %8.2e,  thr/max(PDB) = %5.3f,  Gext = %d\n", iter, sigma, g_hi_cutoff, maxg, g_hi_cutoff/maxg, Gext);

      if (coway != 3) break;   /* g_hi_cutoff is done, so exit the 'do' loop... */

      /* ... otherwise, proceed with the computation of new sigma,
             and continue looping until convergence: */

      /* compute discrepancy wrt intomax (Eq 4 of the notes, p 15): */
      discrep = (sumgs-numgs*g_hi_cutoff)/(maxg-g_hi_cutoff) - intomax;

      /* max increment allowed: */
      maxincr = 0.2*sigma;

      /* compute increment of sigma by a secant step (or arbitrary when iter=0) */
      if (iter==0){
        if (discrep>0) incr = maxincr;
        else incr = -maxincr;
        discrep_a = discrep;
        sigma_a = sigma;
        sigma -= incr;
      }
      else {
        deriv = (discrep-discrep_a)/(sigma-sigma_a);
        incr = discrep/deriv;
        if (incr >  maxincr) incr =  maxincr;
        if (incr < -maxincr) incr = -maxincr;
        discrep_a = discrep;
        sigma_a = sigma;
        sigma -= incr;
      }
      
      iter++;

    } while (fabs(incr)>0.01 && iter<100);
  /** end of determination of threshold and sigma for g_phi_hi **/

  }

  fprintf(stderr, "\n");

  /* times initialization */
  istep=0; iconf=0; nsslp=0;
  curtime = 0.0;    /* current time */

  smd = 0.0;        /* current sum of rmsdisp's */


/* ===> */
newstep:

  chainid = pdb_model_all[0].chain[0];
  resnum  = pdb_model_all[0].seq;
  strcpy(resname, pdb_model_all[0].res);
  sccx=0.0; sccy=0.0; sccz=0.0; scm=0.0;
  atin[1]=atin[2]=atin[3]=atin[4]=atin[5]=0;
  j=0;   /* pseudo-atom index */
  k=0;   /* residue index */
  props[0].rt = rtype(pdb_model_all, 0);
  props[0].epoc = 1;   /* 1st residue of pdb */
  for(i=0;i<num_atoms_all;i++){
    if(   (pdb_model_all[i].seq     != resnum)
       || (pdb_model_all[i].chain[0] != chainid) ){
      if(atin[1] != 0){
        pdb_model = (PDB *) realloc(pdb_model, (j+1)*sizeof(PDB));
        ptr_check(pdb_model,311);
        pdb_model[j].x = pdb_temp[1].x;
        pdb_model[j].y = pdb_temp[1].y;
        pdb_model[j].z = pdb_temp[1].z;
        pdb_model[j].type[0] = 1;
        pdb_model[j].occupancy = mN+mH;  /* mass of pseudo-atom */
        pdb_model[j].chain[0] = chainid;
        pdb_model[j].chain[1] = '\0';
        pdb_model[j].seq = resnum;
        pdb_model[j].footnote = k;         /* residue index */
        strcpy(pdb_model[j].res, resname);
        j++;
      }
      if(atin[2] != 0){
        pdb_model = (PDB *) realloc(pdb_model, (j+1)*sizeof(PDB));
        ptr_check(pdb_model,312);
        pdb_model[j].x = pdb_temp[2].x;
        pdb_model[j].y = pdb_temp[2].y;
        pdb_model[j].z = pdb_temp[2].z;
        pdb_model[j].type[0] = 2;
        pdb_model[j].occupancy = mC+mH;  /* mass of pseudo-atom */
        if(props[k].rt==6) pdb_model[j].occupancy += mH;
        pdb_model[j].chain[0] = chainid;
        pdb_model[j].chain[1] = '\0';
        pdb_model[j].seq = resnum;
        pdb_model[j].footnote = k;         /* residue index */
        strcpy(pdb_model[j].res, resname);
        j++;
      }
      if(atin[3] != 0){
        pdb_model = (PDB *) realloc(pdb_model, (j+1)*sizeof(PDB));
        ptr_check(pdb_model,313);
        pdb_model[j].x = pdb_temp[3].x;
        pdb_model[j].y = pdb_temp[3].y;
        pdb_model[j].z = pdb_temp[3].z;
        pdb_model[j].type[0] = 3;
        pdb_model[j].occupancy = mC+mO;  /* mass of pseudo-atom */
        pdb_model[j].chain[0] = chainid;
        pdb_model[j].chain[1] = '\0';
        pdb_model[j].seq = resnum;
        pdb_model[j].footnote = k;         /* residue index */
        strcpy(pdb_model[j].res, resname);
        j++;
      }
      if(atin[4] != 0){
        pdb_model = (PDB *) realloc(pdb_model, (j+1)*sizeof(PDB));
        ptr_check(pdb_model,314);
        pdb_model[j].x = pdb_temp[4].x;
        pdb_model[j].y = pdb_temp[4].y;
        pdb_model[j].z = pdb_temp[4].z;
        pdb_model[j].type[0] = 4;
        pdb_model[j].occupancy = m45[props[k].rt][0];  /* mass of pseudo-atom */
        pdb_model[j].chain[0] = chainid;
        pdb_model[j].chain[1] = '\0';
        pdb_model[j].seq = resnum;
        pdb_model[j].footnote = k;         /* residue index */
        strcpy(pdb_model[j].res, resname);
        j++;
      }
      if(atin[5] != 0){
        sccx /= scm; sccy /= scm; sccz /= scm;
        pdb_model = (PDB *) realloc(pdb_model, (j+1)*sizeof(PDB));
        ptr_check(pdb_model,315);
        pdb_model[j].x = sccx;
        pdb_model[j].y = sccy;
        pdb_model[j].z = sccz;
        pdb_model[j].type[0] = 5;
        pdb_model[j].occupancy = m45[props[k].rt][1];  /* mass of pseudo-atom */
        pdb_model[j].chain[0] = chainid;
        pdb_model[j].chain[1] = '\0';
        pdb_model[j].seq = resnum;
        pdb_model[j].footnote = k;         /* residue index */
        strcpy(pdb_model[j].res, resname);
        j++;
      }

      k++;
      props = (tri *) realloc(props, (k+1)*sizeof(tri));
      ptr_check(props,51);
      props[k].rt = rtype(pdb_model_all, i);

      if(pdb_model_all[i].chain[0] != chainid){
        props[k-1].epoc = 2;    /* last residue of previous chain */
        props[k].epoc = 1;      /* first residue of new chain */
      }
      else
        props[k].epoc = 0;      /* new residue still on same chain */

      chainid = pdb_model_all[i].chain[0];
      resnum  = pdb_model_all[i].seq;
      strcpy(resname, pdb_model_all[i].res);
      sccx=0.0; sccy=0.0; sccz=0.0; scm=0.0;
      atin[1]=atin[2]=atin[3]=atin[4]=atin[5]=0;
    }

    if(strcmp(pdb_model_all[i].type, "N")==0 && strcmp(pdb_model_all[i].loc, "" )==0){
      pdb_temp[1].x = pdb_model_all[i].x;
      pdb_temp[1].y = pdb_model_all[i].y;
      pdb_temp[1].z = pdb_model_all[i].z;
      atin[1]=1;
    }
    else if(strcmp(pdb_model_all[i].type, "C")==0 && strcmp(pdb_model_all[i].loc, "A")==0){
      pdb_temp[2].x = pdb_model_all[i].x;
      pdb_temp[2].y = pdb_model_all[i].y;
      pdb_temp[2].z = pdb_model_all[i].z;
      atin[2]=1;
    }
    else if(strcmp(pdb_model_all[i].type, "C")==0 && strcmp(pdb_model_all[i].loc, "" )==0){
      pdb_temp[3].x = pdb_model_all[i].x;
      pdb_temp[3].y = pdb_model_all[i].y;
      pdb_temp[3].z = pdb_model_all[i].z;
      atin[3]=1;
    }
    else if(strcmp(pdb_model_all[i].type, "C")==0 && strcmp(pdb_model_all[i].loc, "B")==0){
      pdb_temp[4].x = pdb_model_all[i].x;
      pdb_temp[4].y = pdb_model_all[i].y;
      pdb_temp[4].z = pdb_model_all[i].z;
      atin[4]=1;
    }
    else{
      atom=pdb_model_all[i].type[strlen(pdb_model_all[i].type)-1];
      if(atom!='H' &&
         !(strcmp(pdb_model_all[i].type, "O")==0 &&
           (strcmp(pdb_model_all[i].loc, "" )==0 || strcmp(pdb_model_all[i].loc, "XT")==0))){
        if(atom=='C')
          mass=mC;
        else if(atom=='N')
          mass=mN;
        else if(atom=='O')
          mass=mO;
        else if(atom=='S')
          mass=mS;
        else{
          /* unknown atom type */
          /* printf("Warning: unknown atom type: '%s'\n", pdb_model_all[i].type); */
          /* mass = 14.5; */
          mass = mH;            /* assume H */
        }
        scm += mass;
        sccx += pdb_model_all[i].x * mass;
        sccy += pdb_model_all[i].y * mass;
        sccz += pdb_model_all[i].z * mass;
        atin[5]=1;
      }
    }
  }
  if(atin[1] != 0){
    pdb_model = (PDB *) realloc(pdb_model, (j+1)*sizeof(PDB));
    ptr_check(pdb_model,316);
    pdb_model[j].x = pdb_temp[1].x;
    pdb_model[j].y = pdb_temp[1].y;
    pdb_model[j].z = pdb_temp[1].z;
    pdb_model[j].type[0] = 1;
    pdb_model[j].occupancy = mN+mH;  /* mass of pseudo-atom */
    pdb_model[j].chain[0] = chainid;
    pdb_model[j].chain[1] = '\0';
    pdb_model[j].seq = resnum;
    pdb_model[j].footnote = k;         /* residue index */
    strcpy(pdb_model[j].res, resname);
    j++;
  }
  if(atin[2] != 0){
    pdb_model = (PDB *) realloc(pdb_model, (j+1)*sizeof(PDB));
    ptr_check(pdb_model,317);
    pdb_model[j].x = pdb_temp[2].x;
    pdb_model[j].y = pdb_temp[2].y;
    pdb_model[j].z = pdb_temp[2].z;
    pdb_model[j].type[0] = 2;
    pdb_model[j].occupancy = mC+mH;  /* mass of pseudo-atom */
    if(props[k].rt==6) pdb_model[j].occupancy += mH;
    pdb_model[j].chain[0] = chainid;
    pdb_model[j].chain[1] = '\0';
    pdb_model[j].seq = resnum;
    pdb_model[j].footnote = k;         /* residue index */
    strcpy(pdb_model[j].res, resname);
    j++;
  }
  if(atin[3] != 0){
    pdb_model = (PDB *) realloc(pdb_model, (j+1)*sizeof(PDB));
    ptr_check(pdb_model,318);
    pdb_model[j].x = pdb_temp[3].x;
    pdb_model[j].y = pdb_temp[3].y;
    pdb_model[j].z = pdb_temp[3].z;
    pdb_model[j].type[0] = 3;
    pdb_model[j].occupancy = mC+mO;  /* mass of pseudo-atom */
    pdb_model[j].chain[0] = chainid;
    pdb_model[j].chain[1] = '\0';
    pdb_model[j].seq = resnum;
    pdb_model[j].footnote = k;         /* residue index */
    strcpy(pdb_model[j].res, resname);
    j++;
  }
  if(atin[4] != 0){
    pdb_model = (PDB *) realloc(pdb_model, (j+1)*sizeof(PDB));
    ptr_check(pdb_model,319);
    pdb_model[j].x = pdb_temp[4].x;
    pdb_model[j].y = pdb_temp[4].y;
    pdb_model[j].z = pdb_temp[4].z;
    pdb_model[j].type[0] = 4;
    pdb_model[j].occupancy = m45[props[k].rt][0];  /* mass of pseudo-atom */
    pdb_model[j].chain[0] = chainid;
    pdb_model[j].chain[1] = '\0';
    pdb_model[j].seq = resnum;
    pdb_model[j].footnote = k;         /* residue index */
    strcpy(pdb_model[j].res, resname);
    j++;
  }
  if(atin[5] != 0){
    sccx /= scm; sccy /= scm; sccz /= scm;
    pdb_model = (PDB *) realloc(pdb_model, (j+1)*sizeof(PDB));
    ptr_check(pdb_model,320);
    pdb_model[j].x = sccx;
    pdb_model[j].y = sccy;
    pdb_model[j].z = sccz;
    pdb_model[j].type[0]=5;
    pdb_model[j].occupancy = m45[props[k].rt][1];  /* mass of pseudo-atom */
    pdb_model[j].chain[0] = chainid;
    pdb_model[j].chain[1] = '\0';
    pdb_model[j].seq = resnum;
    pdb_model[j].footnote = k;         /* residue index */
    strcpy(pdb_model[j].res, resname);
    j++;
  }

  props[k].epoc = 2;         /* last residue of pdb */

  num_res = k+1;             /* number of residues */
  num_atoms = j;             /* number of pseudo-atoms */
  tna = 3*num_atoms;
  num_atoms_bb = 3*num_res;  /* number of pseudo-atoms on backbone */

  natheor=0;                 /* number of pseudo-atoms that there should be */
  for(i=0;i<num_res;i++){
    if(props[i].rt == 1){         /* ala */
      props[i].nat=4;
      props[i].nan=2;
    }
    else if(props[i].rt == 6){    /* gly */
      props[i].nat=3;
      props[i].nan=2;
    }
    else if(props[i].rt == 13){   /* pro */
      props[i].nat=5;
      props[i].nan=1;
    }
    else{
      props[i].nat=5;
      props[i].nan=3;
    }
    natheor += props[i].nat;
    if((props[i].epoc==1 && props[i].rt != 13) || props[i].epoc==2)
      props[i].nan -= 1;
  }

  /* indices of first pseudo-atoms of each residue */
  props[0].k1=0;
  for(i=1;i<num_res;i++)
    props[i].k1 = props[i-1].k1+props[i-1].nat;


if(istep==0){

  printf("ddforge> Number of atoms: %d\n", num_atoms_all);
  printf("ddforge> Number of pseudo-atoms: %d\n", num_atoms);
  printf("ddforge> Number of residues: %d\n", num_res);

  if(natheor != num_atoms){
    printf("\nddforge> Atomic structure is missing heavy atoms: theor=%d, actual=%d. Exiting.\n\n",natheor,num_atoms);
    exit(1);
  }

  /* get diameter of atomic structure */
  /* OLD:
  diamstr = 0.0;
  for(k=0;k<num_atoms-1;k++)
    for(l=k+1;l<num_atoms;l++){
      temp = pow(pdb_model[k].x-pdb_model[l].x, 2) +
             pow(pdb_model[k].y-pdb_model[l].y, 2) +
             pow(pdb_model[k].z-pdb_model[l].z, 2);
      if(temp>diamstr) diamstr=temp;
    }
  diamstr=sqrt(diamstr); */
  /* NEW (not exact, but a reasonabe upper bound): */
  xmin = pdb_model[0].x; xmax = xmin;
  ymin = pdb_model[0].y; ymax = ymin;
  zmin = pdb_model[0].z; zmax = zmin;
  for(k=1; k<num_atoms; k++){
    if (pdb_model[k].x < xmin) xmin = pdb_model[k].x;
    if (pdb_model[k].x > xmax) xmax = pdb_model[k].x;
    if (pdb_model[k].y < ymin) ymin = pdb_model[k].y;
    if (pdb_model[k].y > ymax) ymax = pdb_model[k].y;
    if (pdb_model[k].z < zmin) zmin = pdb_model[k].z;
    if (pdb_model[k].z > zmax) zmax = pdb_model[k].z;
  }
  diamstr = sqrt(pow(xmax-xmin,2)+pow(ymax-ymin,2)+pow(zmax-zmin,2));

  /* initial cutoff distance for dampers: */
  dist_cut = max(4.0, 0.5 * diamstr);
  dc2 = pow(dist_cut,2);

  /* final value of dist_cut: */
  /* dcfinal = max(4.0, dist_cut/8.0); */
  if (mpt==1 || mpt==3){    /* resol = map resolution */
    /* dcfinal = max(4.0, resol);                too high */
    /* dcfinal = max(4.0, resol/sqrt(3.0));        good   */
    /* dcfinal = max(4.0, 0.5*resol/sqrt(3.0))   too low; */
    /* better version of the "good" above, and with 7 instead of 4,
       since the latter produces oscillations in the trajectory: */
    dcfinal = max(7.0, 2.0*sigma);
  }
  if (mpt==2)
    dcfinal = max(4.0, dist_cut/8.0);

  /* final speed to be attained as a condition to stop: */
  /* speed1  = 0.001 * force_scale/drag_scale/diamstr;
     speed2  = 0.1 * force_scale/drag_scale/diamstr; */
  speedratio = 0.0001;

  /* read residue-data file -- given as an option, after all the
   mandatory arguments: ddforge .............. [-rdfile xxx] [-dcfile yyy]  */
  /* in case of a pdb as target, the backbone atoms of each residue (of the origin pdb)
   in columns 2-3 will be pulled toward the corresponding atoms of the corresponding
   residue (of the target pdb) in columns 4-5, if these columns are non-empty */
  rdflag = 0;
  for (i=1; i<argc-1; i++){
    if (strcmp(argv[i],"-rdfile")==0){
      sprintf(rdfile, argv[i+1]);
      rdflag = 1; break;
    }
  }
  if (rdflag==0) {
    /* no residue-data file given => define all variables as free (i.e. everything is flexible): */
    fprintf(stderr, "ddforge> No residue-data file specified. Assuming all variables are free.\n");
    corlen = 0;
    srlen = 0;
    coords = (fvars *) realloc(coords, num_res*sizeof(fvars));
    ptr_check(coords,61);
    for (j=0; j<num_res; j++) {
      coords[j].code = 31;
      coords[j].cho  = pdb_model[props[j].k1].chain[0];
      coords[j].rno  = pdb_model[props[j].k1].seq;
    }
    varlen = num_res;
  } else {
    fca = fopen(rdfile, "r");
    if(fca == NULL){
      fprintf(stderr, "ddforge> Error: Can't read residue data file %s\n", rdfile);
      exit(1);
    }
    fgets(line, 100, fca);   /* reads headings */
    for(k=0,j=0,m=0; ;){
      if (fgets(line, 100, fca) == NULL)
        n = 0;
      else
        n = sscanf(line, "%d %s %d %s %d", &code, &chori, &rnori, &chtar, &rntar);
      if(n == 0){
        corlen = k;
        varlen = j;
        srlen  = m;
        break;
      }
      if(n==1 || n==2 || n==4) continue;
      if(n==5){
        rescor = (respair *) realloc(rescor, (k+1)*sizeof(respair));
        ptr_check(rescor,665);
        rescor[k].cho = chori[0];
        rescor[k].rno = rnori;
        rescor[k].cht = chtar[0];
        rescor[k].rnt = rntar;
        k++;
      }
      if(code & 31){    /* replace 31 by 63 if using the full-atom pdb
                           rather than the reduced representation */
        coords = (fvars *) realloc(coords, (j+1)*sizeof(fvars));
        ptr_check(coords,61);
        coords[j].code = code & 31;
        coords[j].cho  = chori[0];
        coords[j].rno  = rnori;
        j++;
      }
      if(code & 64){   /* static backbone atoms */
        statres = (sres *) realloc(statres, (m+1)*sizeof(sres));
        ptr_check(statres,71);
        statres[m].cho  = chori[0];
        statres[m].rno  = rnori;
        m++;
      }
    }
    fclose(fca);
    fprintf(stderr, "ddforge> Read residue-data file '%s'\n", rdfile);
  }
  if(mpt==2 && corlen==0){
    printf("Error: no target residues to pull structure. Exiting.");
    exit(1);
  }

  /* read distance-constraint file -- given as an option, after all the
   mandatory arguments: ddforge .............. [-rdfile xxx] [-dcfile yyy]  */
  dcflag = 0;
  for (i=1; i<argc-1; i++){
    if (strcmp(argv[i],"-dcfile")==0){
      sprintf(dcfile, argv[i+1]);
      dcflag = 1; break;
    }
  }
  if (dcflag==0){
    dclen = 0;
  } else {
    fca = fopen(dcfile, "r");
    if(fca == NULL){
      fprintf(stderr, "ddforge> Error: Can't read distance-constraint file %s\n", dcfile);
      exit(1);
    }
    fgets(line, 100, fca);   /* reads headings */
    for(c=0; ;){
      if (fgets(line, 100, fca) == NULL)
        n = 0;
      else
        n = sscanf(line, "%s %d %s %d", &chori, &rnori, &chtar, &rntar);
      if(n == 0){
        dclen = c;
        break;
      }
      if(n!=4) continue;
      dconst = (respair *) realloc(dconst, (c+1)*sizeof(respair));
      ptr_check(dconst,667);
      dconst[c].cho = chori[0];
      dconst[c].rno = rnori;
      dconst[c].cht = chtar[0];
      dconst[c].rnt = rntar;
      c++;
    }
    fclose(fca);
    printf("ddforge> Read distance-constraint file '%s'\n", dcfile);
    printf("ddforge> Number of distance constraints:  %d\n", dclen);
  }

  /* add residue indices to the rescor list, if non-empty */
  for(k=0;k<corlen;k++){
    /* iko = index of residue rescor[k].rno (origin pdb) */
    for(i=0;i<num_res;i++)
      if( (pdb_model[props[i].k1].seq == rescor[k].rno) &&
           pdb_model[props[i].k1].chain[0] == rescor[k].cho)
        break;
    rescor[k].iko = i;
    /* ikt = index of residue rescor[k].rnt (target pdb) */
    for(i=0;i<nrestar;i++)
      if( (tca[i].rnt == rescor[k].rnt) &&
           tca[i].cht == rescor[k].cht)
        break;
    rescor[k].ikt = i;
  }

  /* add residue indices to the coords list */
  for(k=0;k<varlen;k++){
    for(i=0;i<num_res;i++)
      if( (pdb_model[props[i].k1].seq == coords[k].rno) &&
           pdb_model[props[i].k1].chain[0] == coords[k].cho)
        break;
    coords[k].iko = i;
  }

  /* add residue indices to the statres list */
  for(k=0;k<srlen;k++){
    for(i=0;i<num_res;i++)
      if( (pdb_model[props[i].k1].seq == statres[k].rno) &&
           pdb_model[props[i].k1].chain[0] == statres[k].cho)
        break;
    statres[k].iko = i;
  }

  /* add residue indices to the dconst list */
  for(c=0; c<dclen; c++){
    for(i=0; i<num_res; i++)
      if( (pdb_model[props[i].k1].seq == dconst[c].rno) &&
         pdb_model[props[i].k1].chain[0] == dconst[c].cho)
        break;
    dconst[c].iko = i;
    for(i=0; i<num_res; i++)
      if( (pdb_model[props[i].k1].seq == dconst[c].rnt) &&
         pdb_model[props[i].k1].chain[0] == dconst[c].cht)
        break;
    dconst[c].ikt = i;
  }

  /* number of free variables, and correction of wrong codes */
  size=0;
  for(k=0;k<varlen;k++){
    if(coords[k].code &  1){             /* (x, y, z) free */
      if(props[coords[k].iko].epoc == 1)
        size += 3;
      else
        coords[k].code -= 1; /* correction */
    }
    if(coords[k].code &  2){             /* (lambda, theta) free */
      if(props[coords[k].iko].epoc == 1)
        size += 2;
      else
        coords[k].code -= 2; /* correction */
    }
    if(coords[k].code &  4){             /* phi free */
      if(props[coords[k].iko].rt != 13
         || props[coords[k].iko].epoc == 1)    /* added on 10/4/2016 */
        size += 1;
      else
        coords[k].code -= 4; /* correction */
    }
    if(coords[k].code &  8){             /* psi free */
      if(props[coords[k].iko].epoc != 2)
        size += 1;
      else
        coords[k].code -= 8; /* correction */
    }
    if(coords[k].code & 16){             /* chi1 free */
      if(props[coords[k].iko].rt !=  1 &&
         props[coords[k].iko].rt !=  6 &&
         props[coords[k].iko].rt != 13)
        size += 1;
      else
        coords[k].code -= 16; /* correction */
    }
  /* if(coords[k].code & 32) size +=  ;      all other chi's free (unused at this time) */
  }

  printf("ddforge> Number of free variables: %d\n", size);

  /* number of (scalar) constraints: */
  nco1 = srlen*6;   /* for static residues */
  nco2 = dclen;

  /* (the order of the main matrix will be size+nco1+nco2) */

  /* table of free variables: */
  tfv = (fvars *) malloc(size * sizeof(fvars));
  ptr_check(tfv,9);

  for(p=0,j=0;p<varlen;p++){
    chainid = coords[p].cho;
        rno = coords[p].rno;
          R = coords[p].iko;
    if(coords[p].code &  1){
      tfv[j].code = 10; tfv[j].cho = chainid; tfv[j].rno = rno; tfv[j].iko = R; j++;
      tfv[j].code = 11; tfv[j].cho = chainid; tfv[j].rno = rno; tfv[j].iko = R; j++;
      tfv[j].code = 12; tfv[j].cho = chainid; tfv[j].rno = rno; tfv[j].iko = R; j++;
    }
    if(coords[p].code &  2){
      tfv[j].code = 20; tfv[j].cho = chainid; tfv[j].rno = rno; tfv[j].iko = R; j++;
      tfv[j].code = 21; tfv[j].cho = chainid; tfv[j].rno = rno; tfv[j].iko = R; j++;
    }
    if(coords[p].code &  4){
      tfv[j].code = 30; tfv[j].cho = chainid; tfv[j].rno = rno; tfv[j].iko = R; j++;
    }
    if(coords[p].code &  8){
      tfv[j].code = 40; tfv[j].cho = chainid; tfv[j].rno = rno; tfv[j].iko = R; j++;
    }
    if(coords[p].code & 16){
      tfv[j].code = 50; tfv[j].cho = chainid; tfv[j].rno = rno; tfv[j].iko = R; j++;
    }
  }

  /* count # of chains */
  nchains=0;
  for(R=0; R<num_res; R++)
    if(props[R].epoc == 1) nchains++;
  printf("ddforge> Number of chains: %d\n\n",nchains);

  /* list of chains and their data */
  chinfo = (chdata *) malloc(nchains * sizeof(chdata));
  ptr_check(chinfo,10);

  /* search for first and last residues of each chain,
     and store the chain IDs */
  for(R=0,ich=0; R<num_res; R++){
    if(props[R].epoc == 1)
      chinfo[ich].ikf = R;
    if(props[R].epoc == 2){
      chinfo[ich].ikl = R;
      chinfo[ich].cho = pdb_model[props[R].k1].chain[0];
      ich++;
    }
  }

  /* test
  for(ich=0;ich<nchains;ich++){
    printf("chain[%d]: %c, %d, %d\n",ich,chinfo[ich].cho,chinfo[ich].ikf,chinfo[ich].ikl);
  }
  for(j=0; j<size; j++){
    printf("j=%d  tfv[j].cho='%c'\n", j, tfv[j].cho);
  } */

  /* search for first and last (free) vars of each chain */
  /* (-1 for a chain that has no free vars) */
  for(ich=0; ich<nchains; ich++){
    chinfo[ich].ivf = chinfo[ich].ivl = -1;
  }
  for(ich=0,j=0; ich<nchains && j<size; ich++){
    if(chinfo[ich].cho == tfv[j].cho){
      chinfo[ich].ivf = j;
      for(; j<size; j++)
        if(tfv[j].cho != chinfo[ich].cho) break;
      chinfo[ich].ivl = j-1;
    }
  }

  /* max # of sets (i.e. car and sph are each 1 set) of variables
     (all, free and fixed) in a chain:    */
  mnv=0;
  for(ich=0;ich<nchains;ich++)
    if(chinfo[ich].ikl-chinfo[ich].ikf > mnv)
      mnv = chinfo[ich].ikl-chinfo[ich].ikf;
  mnv = (mnv+1)*3+2;

  tempx = (double *) malloc(num_atoms * sizeof(double));
  tempy = (double *) malloc(num_atoms * sizeof(double));
  tempz = (double *) malloc(num_atoms * sizeof(double));

  /* velocity vector fields: */
  veloc = (double *) malloc((npt+1)*(3*num_atoms+1) * sizeof(double));

  /* array of scalar products of velocities at consecutive steps: */
  spv = (double *) malloc(npt * sizeof(double));

  /* array of rmsdisp values at successive steps: */
  mdarr = (double *) malloc((nmdv+1) * sizeof(double));

  
  /* vectors e_i and r_{beta(i)} used in the calculation of der */
  erb = (twovec *) malloc(size * sizeof(twovec));
  ptr_check(erb,18);

  /* Wilson's matrix (only one column is allocated) */
  /* der = (trd *) malloc(num_atoms*size * sizeof(trd)); */
  der = (trd *) malloc(num_atoms * sizeof(trd));
  ptr_check(der,11);

  /* force field */
  wop = (trd *) malloc(num_atoms * sizeof(trd));
  ptr_check(wop,8);

  /* array (\dot{q_1},...,\dot{q_M}) */
  qdots = (double *) malloc(size * sizeof(double));
  ptr_check(qdots,12);

  if(mpt==1)
    printf("ddforge> Conf#  Step#      Time       Speed       Disp.      Dist_cut  Overlap     Cos      Compute time (sec)\n");

  if(mpt==2)
    printf("ddforge> Conf#  Step#      Time       Speed       Disp.      Dist_cut    Rmsd      Cos      Compute time (sec)\n");

  if(mpt==3)
    printf("ddforge> Conf#  Step#      Time       Speed       Disp.      Dist_cut  Overlap    Rmsd      Cos      Compute time (sec)\n");

  fflush(stdout);

}  /* end of "if(istep==0)" */


  /* COMPUTE FORCE FIELD  *************************************************/

  if(mpt==1 || mpt==3){
    exey = g_extx * g_exty;

    /* reset all voxels to 0: */
    for(m=0;m<g_nvox;m++)
      g_phi_hi[m] = 0.0;

    /* compute g_phi_hi */
    for(k=0; k<num_atoms_all; k++){

      /* mass of this atom */
      atom = pdb_model_all[k].type[strlen(pdb_model_all[k].type)-1];
      if(atom=='C')
        mass=mC;
      else if(atom=='H')
        continue;                  /* H atoms don't generate density */
      else if(atom=='O')
        mass=mO;
      else if(atom=='N')
        mass=mN;
      else if(atom=='S')
        mass=mS;
      else
        continue;                  /* unknown atom */

      /* coords. of this atom */
      p1[0] = pdb_model_all[k].x;
      p1[1] = pdb_model_all[k].y;
      p1[2] = pdb_model_all[k].z;
  
      /* indices (on grid) of the lower limit of the region covered by
         the Gaussian mask centered at this atom: */
      ixa =  ceil((p1[0]-g_orix-rrc)/g_width);
      iya =  ceil((p1[1]-g_oriy-rrc)/g_width);
      iza =  ceil((p1[2]-g_oriz-rrc)/g_width);
  
      /* indices of the upper limit: */
      ixb = floor((p1[0]-g_orix+rrc)/g_width);
      iyb = floor((p1[1]-g_oriy+rrc)/g_width);
      izb = floor((p1[2]-g_oriz+rrc)/g_width);
  
      if(ixa>=g_extx || iya>=g_exty || iza>=g_extz) continue;
      if(ixb<0 || iyb<0 || izb<0) continue;
  
      if(ixa<0) ixa=0; if(iya<0) iya=0; if(iza<0) iza=0;
      if(ixb>=g_extx) ixb=g_extx-1;
      if(iyb>=g_exty) iyb=g_exty-1;
      if(izb>=g_extz) izb=g_extz-1;

      /* terms to be substracted below: */
      q0 = p1[0]/Gs; q1 = p1[1]/Gs; q2 = p1[2]/Gs;
  
      for(iz=iza;iz<=izb;iz++){
        eeiz = exey*iz;
        izg = floor(iz*wrat+boz-q2);
        eeizg = Gext2*izg;
        for(iy=iya;iy<=iyb;iy++){
          indzy = eeiz + g_extx*iy;
          iyg = floor(iy*wrat+boy-q1);
          indzyg = eeizg + Gext*iyg;
          for(ix=ixa;ix<=ixb;ix++){
            ixg = floor(ix*wrat+box-q0);
            g_phi_hi[indzy+ix] += mass * Gkernel[indzyg+ixg];
          }
        }
      }
    }

    /* above, Gkernel[indzyg+ixg] equals exp(-temp/tsg12) */

    /* root of g_phi_hi: */
    /* if(root_index != 1)
    for(m=0;m<g_nvox;m++)
      g_phi_hi[m] = pow(g_phi_hi[m], 1.0/root_index); */

    /* for testing: */
    /* write_vol("test_sit.sit", g_width, g_orix, g_oriy, g_oriz, g_extx, g_exty, g_extz, g_phi_hi);
    exit(1); */

  
    if(Lcapr>0.0){

      /* cap on g_hi_cutoff: */
      g_hi_cutoff = min(g_hi_cutoff, Lcapr*g_phi_sorted[g_nvox-1]);
    
      /* ratio between cutoff and max:
      rcmd = g_hi_cutoff/g_phi_sorted[g_nvox-1];  */
  
      /* cut off g_phi_hi: */
      /* set density values below g_hi_cutoff to zero, and lower the rest: */
      for(m=0; m<g_nvox; m++){
        if(g_phi_hi[m]<g_hi_cutoff)
          g_phi_hi[m] = 0.0;
        else
          g_phi_hi[m] -= g_hi_cutoff;
      }

      if (heflag==1){
        /* do histogram equalization of g_phi_hi */
        hist_eq(g_phi_hi, g_nvox);
      }

      /* compute the integral of g_phi_hi */
      for(m=0,sumg=0.0; m<g_nvox; m++)
        sumg += g_phi_hi[m];
    }
    else{
      sumg = sumg0;
      /* rcmd = 0.0; */
    }

    /* compute max (or quantile) of g_phi_hi: */
    maxg = 0.0;
    for(m=0; m<g_nvox; m++)
      if(g_phi_hi[m]>maxg)
        maxg = g_phi_hi[m];

    /* normalize g_phi_hi: */
    if(maxg != 0.0){
      cnorm = sumf/sumg;

      for(m=0; m<g_nvox; m++)
        g_phi_hi[m] *= cnorm;
      
      /* adjust maxg and sumg: */
      maxg *= cnorm;
      sumg *= cnorm;
    }

    /* compute Int(|f-g|)/(Int(f)+Int(g)) ("1-overlap"),
       to be printed out during the integration: */
    sumd = 0.0;
    for(m=0;m<g_nvox;m++)
      sumd += fabs(g_phi_lo[m]-g_phi_hi[m]);
    sumd /= (sumf+sumg);

    /* compute mean{|f-g| : f>0}/2:
    meand = 0.0;
    for(m=0; m<g_nvox; m++)
      if(g_phi_lo[m]>0.0)
        meand += fabs(g_phi_lo[m]-g_phi_hi[m]);
    meand /= (2*Ngv);  */

    /* define the force field proper: */
    for(k=0; k<num_atoms; k++){

      /* coords. of this atom */
      p1[0] = pdb_model[k].x;
      p1[1] = pdb_model[k].y;
      p1[2] = pdb_model[k].z;

      if(fftype==2 || fftype==3){
        /* indices of atom's position: */
        fx = (p1[0]-g_orix)/g_width;
        fy = (p1[1]-g_oriy)/g_width;
        fz = (p1[2]-g_oriz)/g_width;
        ixa = floor(fx);
        iya = floor(fy);
        iza = floor(fz);
       
         /* compute g(r_k), called below 'hi_atom': */
        if(ixa<0 || ixa>g_extx-2 || iya<0 || iya>g_exty-2 || iza<0 || iza>g_extz-2){
          hi_atom = 0.1 * maxg;   /* a low value, so that this atom will be attracted
                                     by periferic regions of the map */
        }
        else{
          dx = fx-ixa;
          dy = fy-iya;
          dz = fz-iza;
         
          indza = exey*iza;   indzb = indza+exey;
          indya = g_extx*iya; indyb = indya+g_extx;
          
          phiaaa = g_phi_hi[indza + indya + ixa];
          phiaab = g_phi_hi[indzb + indya + ixa];
          phiaba = g_phi_hi[indza + indyb + ixa];
          phiabb = g_phi_hi[indzb + indyb + ixa];
          phibaa = g_phi_hi[indza + indya + ixa+1];
          phibab = g_phi_hi[indzb + indya + ixa+1];
          phibba = g_phi_hi[indza + indyb + ixa+1];
          phibbb = g_phi_hi[indzb + indyb + ixa+1];

          /* interpolate: (if you can find a way using less than 7 products,
             I would like to know :-) */
          hi_atom =
            (((phibbb-phiabb-phibba+phiaba-phibab+phibaa+phiaab-phiaaa)*dz +
              phibba-phiaba-phibaa+phiaaa)*dy + phibaa-phiaaa)*dx +
            (phiaba-phiaaa)*dy +
            ((phibab-phibaa-phiaab+phiaaa)*dx + (phiabb-phiaba-phiaab+phiaaa)*dy +
             phiaab-phiaaa)*dz + phiaaa;
        }
      }

      F0x = 0.0; F0y = 0.0; F0z = 0.0;
      
      if(fftype==1){
        for(iz=0,eeiz=0,zp=g_oriz; iz<g_extz; iz++,eeiz+=exey,zp+=g_width){
          /* eeiz = exey*iz; zp = g_oriz + iz*g_width; (old) */
          tz2 = pow(zp-p1[2], 2);
          for(iy=0,indzy=eeiz,yp=g_oriy; iy<g_exty; iy++,indzy+=g_extx,yp+=g_width){
            /* indzy = eeiz + g_extx*iy; yp = g_oriy + iy*g_width; (old) */
            ty2 = pow(yp-p1[1], 2);
            for(ix=0,index=indzy,xp=g_orix; ix<g_extx; ix++,index++,xp+=g_width){
              /* index = indzy+ix; xp = g_orix + ix*g_width; (old) */

              if(g_phi_lo[index]==0.0) continue;

              coef = max(g_phi_lo[index]-g_phi_hi[index], 0.0);

           /* coef = g_phi_lo[index]-g_phi_hi[index]; */

           /* if(coef>0.0) sgnc = 1;
              else sgnc = -1;
              coef = sgnc * sqrt(fabs(coef)); */
              
           /* if(coef>0.0) coef = 1.0;
              if(coef<0.0) coef =-1.0; */

              /* if(fabs(coef)<thresval) continue; */
              /* OR: (I think the following is more accurate) */
              if(coef==0.0) continue;

              temp = pow(xp-p1[0], 2) + ty2 + tz2;
              if(temp==0.0) continue;

              /* coef /= (a1*pow(temp,0.25) + pow(temp,alo2)); (OLD) */

              if(temp>ffdc2) continue;       /* distance > cutoff */

              if (temp>a2){
                coef *= b1/pow(temp,fps);
                if(fpt>0) coef /= sqrt(temp);
              }
              else {
                if(ffsl==1)
                  coef *= a1/pow(temp,0.25);
                if(ffsl==2)
                  ;    /* multiply by 1 */
              }

              F0x += coef * (xp - p1[0]);
              F0y += coef * (yp - p1[1]);
              F0z += coef * (zp - p1[2]);
            }
          }
        }
      }

      if(fftype==2 || fftype==3){

        lodens = max(hi_atom-deltadens/2.0,mingv);
        /* lodens=g_phi_lo_sorted[0];  // old  */
        m1 = bisearch(g_phi_lo_sorted, g_nvox, lodens);
        m1 = min(max(0,m1),g_nvox-2);

        hidens = min(hi_atom+deltadens/2.0,essmaxlo);
        /* hidens=g_phi_lo_sorted[g_nvox-1];  //  old */
        m2 = bisearch(g_phi_lo_sorted, g_nvox, hidens);
        m2 = min(max(0,m2),g_nvox-2)+1;

        /* here one may want to have: m2-m1+1 >= nvoxlayer
           (by increasing m2 and decreasing m1 if the condition is not met),
           but I haven't tried this. */
        
        for(m=m1;m<=m2;m++){
          index = indexes[m];

          if(g_phi_lo[index]==0.0) continue;
          
          /* filter points (function 'gamma' in the notes): */
          /* if(g_phi_lo[index] > deltadens) continue; */

          coef = max(g_phi_lo[index]-g_phi_hi[index], 0.0);

          /* coef = g_phi_lo[index]-g_phi_hi[index]; */

          /* coef = fabs(g_phi_lo[index]-g_phi_hi[index]); */

          /* coef = 1.0; */

          /* if(fabs(coef)<thresval) continue; */
          /* OR: (I think the following is more accurate) */
          if(coef==0.0) continue;
      
          /* decompose index: */
          iz = index/exey;
          indexc = index-exey*iz;
          iy = indexc/g_extx;
          ix = indexc-g_extx*iy;
          
          /* get coordinates of point: */
          xp = g_orix + ix*g_width;
          yp = g_oriy + iy*g_width;
          zp = g_oriz + iz*g_width;

          temp = pow(xp-p1[0], 2) + pow(yp-p1[1], 2) + pow(zp-p1[2], 2);
          if(temp==0.0) continue;

          /* coef /= (a1*pow(temp,0.25) + pow(temp,alo2)); (OLD) */

          if(temp>ffdc2) continue;       /* distance > cutoff */

          if (temp>a2){
            coef *= b1/pow(temp,fps);
            if(fpt>0) coef /= sqrt(temp);
          }
          else {
            if(ffsl==1)
              coef *= a1/pow(temp,0.25);
            if(ffsl==2)
              ;    /* multiply by 1 */
          }

          F0x += coef * (xp - p1[0]);
          F0y += coef * (yp - p1[1]);
          F0z += coef * (zp - p1[2]);
        }
      }

      Fnorm = force_scale * pdb_model[k].occupancy;

      wop[k].x = Fnorm * F0x;
      wop[k].y = Fnorm * F0y;
      wop[k].z = Fnorm * F0z;
    }
  }    /* end of "if(mpt==1 || mpt==3)" */

  if(mpt==2 || mpt==3){
    /* if no map, set to 0 to cover those atoms not directly pulled
       by the force field: */
    if(mpt==2)
      for(k=0; k<num_atoms; k++){
        wop[k].x = 0.0;
        wop[k].y = 0.0;
        wop[k].z = 0.0;
      }

    /* backbone Rmsd of current structure with respect to target,
       to be printed out during the integration: */
    sumdp = 0.0;

    for(p=0;p<corlen;p++){
      /* index of the first atom of residue rescor[p].rno (origin pdb) */
      k1 = props[rescor[p].iko].k1;

      for(k=k1;k<k1+3;k++){
        /* coords. of atom k of origin pdb */
        p1[0] = pdb_model[k].x;
        p1[1] = pdb_model[k].y;
        p1[2] = pdb_model[k].z;

        /* coords. of corresp. atom of residue rescor[p].rnt (target pdb) */
        if(k==k1){
          p2[0] = tca[rescor[p].ikt].xN;
          p2[1] = tca[rescor[p].ikt].yN;
          p2[2] = tca[rescor[p].ikt].zN;
        }
        if(k==k1+1){
          p2[0] = tca[rescor[p].ikt].xCa;
          p2[1] = tca[rescor[p].ikt].yCa;
          p2[2] = tca[rescor[p].ikt].zCa;
        }
        if(k==k1+2){
          p2[0] = tca[rescor[p].ikt].xC;
          p2[1] = tca[rescor[p].ikt].yC;
          p2[2] = tca[rescor[p].ikt].zC;
        }

        /* square of distance: */
        temp = pow(p2[0]-p1[0],2) + pow(p2[1]-p1[1],2) + pow(p2[2]-p1[2],2);
        sumdp += temp;

        if(ffsl==1)
          coef = force_scale * pdb_model[k].occupancy * a1 / pow(temp,0.25);
        if(ffsl==2)
          coef = force_scale * pdb_model[k].occupancy;
        
        if(mpt==2){
          wop[k].x = coef * (p2[0]-p1[0]);
          wop[k].y = coef * (p2[1]-p1[1]);
          wop[k].z = coef * (p2[2]-p1[2]);
        }

        if(mpt==3){
          mapforce2 = pow(wop[k].x,2)+pow(wop[k].y,2)+pow(wop[k].z,2);
          coefb = Kforce*sqrt(mapforce2/temp);
          if(coefb>coef) coef=coefb;
          coef *= wpdbtar;
          wop[k].x *= (1-wpdbtar);
          wop[k].y *= (1-wpdbtar);
          wop[k].z *= (1-wpdbtar);
          wop[k].x += coef * (p2[0]-p1[0]);
          wop[k].y += coef * (p2[1]-p1[1]);
          wop[k].z += coef * (p2[2]-p1[2]);
        }
      }
    }
    sumdp = sqrt(sumdp/(3*corlen));
  }

  /**************** compute the vectors e_j and r_{beta(j)} ******************/
  for(p=0,j=0;p<varlen;p++){        /* p = index in the 'coords' array */
                                    /* j = index of the current free gen. variable */
    i  = coords[p].iko;             /* i = residue where q_j belongs */
    k1 = props[i].k1;               /* k1= index of 1st atom of residue i */

    /* backbone atoms of this residue: */
    rc1[0]=pdb_model[k1].x;   rc1[1]=pdb_model[k1].y;   rc1[2]=pdb_model[k1].z;
    rc2[0]=pdb_model[k1+1].x; rc2[1]=pdb_model[k1+1].y; rc2[2]=pdb_model[k1+1].z;
    rc3[0]=pdb_model[k1+2].x; rc3[1]=pdb_model[k1+2].y; rc3[2]=pdb_model[k1+2].z;

    if(coords[p].code &  1){               /* (q_j, q_{j+1}, q_{j+2}) = (x,y,z) of N_i=pdb_model[k1] */
      erb[j].e.x   = erb[j].e.y   = erb[j].e.z   = 0.0;
      erb[j+1].e.x = erb[j+1].e.y = erb[j+1].e.z = 0.0;
      erb[j+2].e.x = erb[j+2].e.y = erb[j+2].e.z = 0.0;
      erb[j].b.x   = 1.0; erb[j].b.y   = 0.0; erb[j].b.z   = 0.0;
      erb[j+1].b.x = 0.0; erb[j+1].b.y = 1.0; erb[j+1].b.z = 0.0;
      erb[j+2].b.x = 0.0; erb[j+2].b.y = 0.0; erb[j+2].b.z = 1.0;
      j += 3;
    }

    if(coords[p].code &  2){               /* (q_j, q_{j+1}) = (lambda,theta) */
      e[0] = pdb_model[k1].x; e[1] = pdb_model[k1].y; e[2] = pdb_model[k1].z;
      temp=sqrt(e[0]*e[0]+e[1]*e[1]+e[2]*e[2]);
      for(m=0;m<3;m++) e[m] /= temp;       /* e for lambda*/
      erb[j].e.x = e[0]; erb[j].e.y = e[1]; erb[j].e.z = e[2];
      erb[j].b.x = 0.0;  erb[j].b.y = 0.0;  erb[j].b.z = 0.0;
      j++;

      e[0] = rc1[1]*rc2[2]-rc1[2]*rc2[1];
      e[1] = rc1[2]*rc2[0]-rc1[0]*rc2[2];
      e[2] = rc1[0]*rc2[1]-rc1[1]*rc2[0];
      temp=sqrt(e[0]*e[0]+e[1]*e[1]+e[2]*e[2]);
      for(m=0;m<3;m++) e[m] /= temp;       /* e for theta*/
      erb[j].e.x = e[0];   erb[j].e.y = e[1];   erb[j].e.z = e[2];
      erb[j].b.x = rc1[0]; erb[j].b.y = rc1[1]; erb[j].b.z = rc1[2];
      j++;
    }

    if(coords[p].code &  4){               /* q_j = phi_i */
      e[0] = rc2[0]-rc1[0];
      e[1] = rc2[1]-rc1[1];
      e[2] = rc2[2]-rc1[2];
      temp=sqrt(e[0]*e[0]+e[1]*e[1]+e[2]*e[2]);
      for(m=0;m<3;m++) e[m] /= temp;       /* e for phi */
      erb[j].e.x = e[0];   erb[j].e.y = e[1];   erb[j].e.z = e[2];
      erb[j].b.x = rc1[0]; erb[j].b.y = rc1[1]; erb[j].b.z = rc1[2];
      j++;
    }

    if(coords[p].code &  8){               /* q_j  = psi_i */
      e[0] = rc3[0]-rc2[0];
      e[1] = rc3[1]-rc2[1];
      e[2] = rc3[2]-rc2[2];
      temp=sqrt(e[0]*e[0]+e[1]*e[1]+e[2]*e[2]);
      for(m=0;m<3;m++) e[m] /= temp;       /* e for psi */
      erb[j].e.x = e[0];   erb[j].e.y = e[1];   erb[j].e.z = e[2];
      erb[j].b.x = rc2[0]; erb[j].b.y = rc2[1]; erb[j].b.z = rc2[2];
      j++;
    }

    if(coords[p].code & 16){               /* q_j = chi_i */
      e[0] = pdb_model[k1+3].x - rc2[0];
      e[1] = pdb_model[k1+3].y - rc2[1];
      e[2] = pdb_model[k1+3].z - rc2[2];
      temp=sqrt(e[0]*e[0]+e[1]*e[1]+e[2]*e[2]);
      for(m=0;m<3;m++) e[m] /= temp;       /* e for chi */
      erb[j].e.x = e[0];   erb[j].e.y = e[1];   erb[j].e.z = e[2];
      erb[j].b.x = rc2[0]; erb[j].b.y = rc2[1]; erb[j].b.z = rc2[2];
      j++;
    }
  }

  /*************** SOLVE FOR THE QDOTS *****************/
  info = getqdots(pdb_model, num_atoms, num_res, dc2, tfv, size, mnv, nco1, nco2,
                  coords, varlen, statres, srlen, dconst, dclen, nchains, chinfo,
                  props, wop, erb, der, qdots, damp_scale, drag_scale);

  if(info != 0){
    printf("\nddforge> Error: singular matrix encountered. Stopping.\n");
    fprintf(stderr, "\nddforge> Error: singular matrix encountered. Stopping.\n");
    free(pdb_model);
    exit(2);
  }

  /* compute \dot{r_k} for all pseudo-atoms */
  for(k=0; k<num_atoms; k++){
    tempx[k]=0.0; tempy[k]=0.0; tempz[k]=0.0;
  }
  for(j=0;j<size;j++){
    Wilson(-1, j, num_atoms, pdb_model, props, erb, tfv, der);
    for(k=0; k<num_atoms; k++){
      tempx[k] += der[k].x * qdots[j];
      tempy[k] += der[k].y * qdots[j];
      tempz[k] += der[k].z * qdots[j];
    }
  }

  /* 3N-dim norm of the velocity vector field (normalized to be an rms value),
     and maximum velocity: */
  rmsvel=0.0; maxvel=0.0;
  for(k=0; k<num_atoms; k++){
    temp = pow(tempx[k],2)+pow(tempy[k],2)+pow(tempz[k],2);
    rmsvel += temp;
    if(temp>maxvel) maxvel=temp;
  }
  rmsvel = sqrt(rmsvel/num_atoms);
  maxvel = sqrt(maxvel);

  /* save initial velocity for ending criterion */
  if(istep==0)
    rmsvel0 = rmsvel;

  /* append velocity vectors to the veloc array,
     storing the 3N-dim norm in the first column: */
  if(istep<=npt){
    isnc = istep*(tna+1);
    veloc[isnc] = rmsvel;
    for(k=0; k<num_atoms; k++){
      veloc[isnc+(3*k+1)] = tempx[k];
      veloc[isnc+(3*k+2)] = tempy[k];
      veloc[isnc+(3*k+3)] = tempz[k];
    }
  }
  else{
    for(i=0;i<npt;i++){
      isnc = i*(tna+1);
      for(icol=0;icol<=tna;icol++)
        veloc[isnc+icol] = veloc[isnc+tna+1+icol];
    }
    isnc = npt*(tna+1);
    veloc[isnc] = rmsvel;
    for(k=0; k<num_atoms; k++){
      veloc[isnc+(3*k+1)] = tempx[k];
      veloc[isnc+(3*k+2)] = tempy[k];
      veloc[isnc+(3*k+3)] = tempz[k];
    }
  }

  /* get new rmsdisp: */
  if(istep>=npt){
    /* scalar product of the current and previous velocity vectors: */
    spvv = 0.0;
    isnc = (npt-1)*(tna+1);
    for(icol=1;icol<=tna;icol++)
      spvv += veloc[isnc+icol] * veloc[isnc+tna+1+icol];
    
    /* normalize, since the veloc[...+icol] above with icol>0 are not rms values */
    spvv /= num_atoms;

    /* previous rmsvel: */
    rmsvelp = veloc[isnc];

    /* cos of angle between current and previous velocity vectors: */
    if(rmsvel==0.0)
      cavv = 0.0;
    else
      cavv = spvv/(rmsvelp*rmsvel);

    infty=10.0;

    /* IDEA 1: */
/*  if(rmsvel==0.0)
      mdf = 0.0;
    else{
      rho = pow(veloc[isnc]/rmsvel,2) - 1.0;
      if(rho>0.0)
        mdf = (cavv+sqrt(rho+pow(cavv,2)))/rho;
      if(rho==0.0){
        if(cavv<0)
          mdf = -0.5/cavv;
        else
          mdf = infty;
      }
      if(rho<0.0){
        if(rho+pow(cavv,2)<0.0)
          mdf = infty;
        else{
          if(cavv<0.0)
            mdf = (cavv+sqrt(rho+pow(cavv,2)))/rho;
          else
            mdf = infty;
        }
      }
    }
*/

    /* IDEAS 2,2a,2b,2c,2d */
/*  if(rmsvel==0.0)
      mdf = 0.0;
    else{
      rho = veloc[isnc]/rmsvel;
      if(rho<=cavv)
        mdf = infty;
      else{
        mdf = fabs(cavv)/(rho-cavv);                   idea 2
        mdf = (1.0+fabs(cavv))/(2.0*(rho-cavv));       idea 2a
        mdf = (1.0+2.0*fabs(cavv))/(3.0*(rho-cavv));   idea 2b
        mdf = (2.0+fabs(cavv))/(3.0*(rho-cavv));       idea 2c 
        mdf = sqrt(fabs(cavv))/(rho-cavv);             idea 2d
      }
    }
*/

    /* IDEA 5 */
    if(rmsvel==0.0)
      mdf = 0.0;
    else{
      rho = rmsvelp/rmsvel;
      if(rho<=cavv)
        mdf = infty;
      else
        mdf = 1.0/(rho-cavv);
    }


    /* IDEA 6 */
/*  if(rmsvel==0.0)
      mdf = 0.0;
    else{
      rho = veloc[isnc]/rmsvel;
      if(cavv<0.0 || rho*cavv>1.0)
        mdf = cavv/(rho*cavv-1.0);
      else
        mdf = infty;
    }
*/

    /* IDEA 9 */
/*  rho = rmsvel/veloc[isnc];
    mdf = rho*(1.0 + rho*cavv/(1.0+rho));
*/


    mdf = min(1.2, mdf);                        /* cap on mdf */
    rmsdisp *= mdf;
    dispcap = min(maxdisp, mdf*dispcap);       /* NEW: 6/27/2017 */

    /* cap on rmsdisp: */
/*  rmsdisp = min(4.0, rmsdisp);    (OLD) */
    if(maxvel>0.0)
      rmsdisp = min(rmsdisp, dispcap*rmsvel/maxvel);
               /* (so that: maximum displacement <= dispcap) */
  }

  /* append rmsdisp to array and compute average: */
  if(istep<=nmdv)
    mdarr[istep]=rmsdisp;
  else{
    for(i=0;i<nmdv;i++)
      mdarr[i]=mdarr[i+1];
    mdarr[nmdv]=rmsdisp;
  }

  if(istep>=nmdv){
    mdave=0.0;
    for(i=0;i<=nmdv;i++)
      mdave += mdarr[i];
    mdave /= nmdv+1;
  }


  if(istep==0)
    flagwr=1;    /* current conformation is already in a file (original pdb) */
  else
    flagwr=0;    /* current conformation not yet written to file */

  tb1 = clock();
  timing = elapsed_time(tb0,tb1);

  if(smd >= dboc || istep==0){
    if(mpt==1)
      printf("ddforge> %5d  %5d  %10.3e  %10.3e  %10.3e  %10.3e   %5.2f   %8.5f      %10.3e\n",
             iconf, istep, curtime, rmsvel, rmsdisp, dist_cut, (1.0-sumd)*100, cavv, timing);
    if(mpt==2)
      printf("ddforge> %5d  %5d  %10.3e  %10.3e  %10.3e  %10.3e   %5.2f   %8.5f      %10.3e\n",
             iconf, istep, curtime, rmsvel, rmsdisp, dist_cut, sumdp, cavv, timing);
    if(mpt==3)
      printf("ddforge> %5d  %5d  %10.3e  %10.3e  %10.3e  %10.3e   %5.2f   %6.2f   %8.5f      %10.3e\n",
             iconf, istep, curtime, rmsvel, rmsdisp, dist_cut, (1.0-sumd)*100, sumdp, cavv, timing);
    fflush(stdout);
    if(flagwr==0){
      /* if(dist_cut>dcfinal && flagsc1==1){   **OLD**  */
      if(flagsc1==1){
        optpdb(num_atoms_all, pdb_model_all, bname, iconf, file_name);
        /* read in side-chain-optimized pdb: */
        free(pdb_model_all);
        read_pdb(file_name, &num_atoms_all, &pdb_model_all);
      }
      else{
        sprintf(file_name, "%s.deform.%d.pdb", bname, iconf);
        write_pdb(file_name, num_atoms_all, pdb_model_all);
      }
      flagwr = 1;         /* current conformation written to file */
    }
    iconf++;
    smd = 0.0;
    nsslp=0;       /* number of steps since latest printout */
  }
  else{
    if(mpt==1)
      printf("ddforge>        %5d  %10.3e  %10.3e  %10.3e  %10.3e   %5.2f   %8.5f      %10.3e\n",
                    istep, curtime, rmsvel, rmsdisp, dist_cut, (1.0-sumd)*100, cavv, timing);
    if(mpt==2)
      printf("ddforge>        %5d  %10.3e  %10.3e  %10.3e  %10.3e   %5.2f   %8.5f      %10.3e\n",
                    istep, curtime, rmsvel, rmsdisp, dist_cut, sumdp, cavv, timing);
    if(mpt==3)
      printf("ddforge>        %5d  %10.3e  %10.3e  %10.3e  %10.3e   %5.2f   %6.2f   %8.5f      %10.3e\n",
                    istep, curtime, rmsvel, rmsdisp, dist_cut, (1.0-sumd)*100, sumdp, cavv, timing);
    fflush(stdout);
  }

  tb1 = clock();
  timing = elapsed_time(tb0,tb1);

  if(mpt==1 || mpt==3){
    Qs = (double *) realloc(Qs, (istep+1)*sizeof(double));
    Qs[istep] = 1.0-sumd;
  }

  if (istep>=npt){

    /* check the stopping criterion based on the exponential regression ("er"): */
    mer = npt;   // mer could also be set to 0 or to the index where overlap is minimum
    if ((mpt==1 || mpt==3) && istep-mer>9){
      for (S0a=0.0,i=mer; i<=istep; i++)
        S0a += Qs[i];
      /* initial values for er: */
      ber = Qs[istep];
      cer = Qs[istep]-Qs[mer];
      if (cer<=0.0){ytime = -1.0; rtime = -1.0; goto coc;}
      sum1 = 0.0; sum2 = 0.0;
      for (i=mer; i<=istep; i++){
        sum1 += i*log(cer/(1.0-Qs[i]));
        sum2 += i*i;
      }
      ker = sum1/sum2;
      n2 = 0;
      while (n2<maxiters){
        /* make a single Newton step on J3 for k, using the previous ber, cer, ker as initial approx : */
        deltak = 0.1*ker;
        kp = ker - deltak/(J3(Qs, mer, istep, ber, cer, ker+deltak)/J3(Qs, mer, istep, ber, cer, ker) - 1.0);
        /* with this new k, compute new b, c : */
        S1a=0.0; S2a=0.0; S3a=0.0;
        for (i=mer; i<=istep; i++){
          f1 = exp(-kp*i);
          S1a += f1;
          S2a += pow(f1,2);
          S3a += Qs[i]*f1;
        }
        denom = pow(S1a,2) - (istep-mer+1)*S2a;
        if (denom==0.0){ytime = -1.0; rtime = -1.0; goto coc;}
        bp = (S1a*S3a - S0a*S2a)/denom;
        cp = (S3a*(istep-mer+1) - S0a*S1a)/denom;
        /* check convergence: */
        error = fabs((bp - ber)/ber);
        error1 = fabs((cp - cer)/cer);
        if (error1>error) error=error1;
        error1 = fabs((kp - ker)/ker);
        if (error1>error) error=error1;
        ber = bp; cer = cp; ker = kp;
        if(error < epsilon) break; else n2++;
      }
      if (n2<maxiters){
        ytime = log(1.0/(1.0-alpha1))/ker;
        rtime = log(1.0/(1.0-alpha2))/ker;
      } else {
        ytime = -1.0; rtime = -1.0;
      }
    } else {
      ytime = -1.0; rtime = -1.0;
    }
  coc:

    if (rmsvel==0.0 || iconf>100 || ((rmsvel<speedratio*rmsvel0*istep || nsslp>100) && dist_cut<dcfinal)
        || (rtime>0.0 && rtime<istep)){
      /* output current conformation if necessary, and end */
      if(flagwr==0){
        if(mpt==1)
          printf("ddforge> %5d  %5d  %10.3e  %10.3e  %10.3e  %10.3e   %5.2f   %8.5f      %10.3e\n",
                 iconf, istep, curtime, rmsvel, rmsdisp, dist_cut, (1.0-sumd)*100, cavv, timing);
        if(mpt==2)
          printf("ddforge> %5d  %5d  %10.3e  %10.3e  %10.3e  %10.3e   %5.2f   %8.5f      %10.3e\n",
                 iconf, istep, curtime, rmsvel, rmsdisp, dist_cut, sumdp, cavv, timing);
        if(mpt==3)
          printf("ddforge> %5d  %5d  %10.3e  %10.3e  %10.3e  %10.3e   %5.2f   %6.2f   %8.5f      %10.3e\n",
                 iconf, istep, curtime, rmsvel, rmsdisp, dist_cut, (1.0-sumd)*100, sumdp, cavv, timing);
        fflush(stdout);
        if(0)  /* to do scatd to the last conf, replace by: if(flagsc1==1) */
          optpdb(num_atoms_all, pdb_model_all, bname, iconf, file_name);
        else{
          sprintf(file_name, "%s.deform.%d.pdb", bname, iconf);
          write_pdb(file_name, num_atoms_all, pdb_model_all);
        }
      }
      if (rtime>0.0 && rtime<istep)
        printf("ddforge> \nddforge> Safe ending time step = %.0f\n", floor(ytime+0.5));
      else
        printf("ddforge> \nddforge> Safe ending time step could not be determined.\n");
      printf("ddforge> \nddforge> All done. Program exiting normally.\n\n");
    /*system("rm -f temp.pdb"); */
    /*free(pdb0); free(pdb_model); free(pdb_model_all); free(pdb_temp);
      free(props); free(erb); free(der); free(g_phi_hi); free(g_phi_sorted);
      free(coords); free(statres); free(wop); free(tfv); free(chinfo);
      free(tempx); free(tempy); free(tempz); free(qdots); */
      exit(0);
    }
    else {
      /* adjust dist_cut */
      /* way 1 (conservative, but slow during the first part of the trajectory):
      rho = min(rmsvel/rmsvelp, 1.0); */
      /* way 2 (using rmsvel0 as the reference, instead of rmsvelp): */
      rho = min(rmsvel/rmsvel0, 1.0);
      /* way 3 (weighing with istep -- less agressive than way 2):
      rho = min(istep*rmsvel/rmsvel0, 1.0); */

      dist_cut *= max(pow(rho, 0.4), 0.5);
      if(dist_cut<0.95*dcfinal) dist_cut = 0.95*dcfinal;
      dc2 = pow(dist_cut, 2);
      
      nsslr=0;     /* number of steps since latest relaxation */
    }
  }

  /* determine Delta t to get new atomic structure: */
  deltat = rmsdisp/rmsvel;


  /* apply Delta q's to the pdb structure -------------------------------------------- */

  j=0;     /* index of the current free variable */

  for(p=0;p<varlen;p++){        /* p = index in the 'coords' array */

    /* residue index: */
    i = coords[p].iko;

    /* chain ID of chain where residue i belongs: */
    chainid = coords[p].cho;

    /* search for 1st atom of residue i */
    for(k=0; k<num_atoms_all; k++)
      if(pdb_model_all[k].seq == coords[p].rno &&
         pdb_model_all[k].chain[0] == chainid){
        k1=k; break;
      }

    /* search for last atom of residue i */
    for(k=k1; k<num_atoms_all && pdb_model_all[k].seq==coords[p].rno
                              && pdb_model_all[k].chain[0]==chainid; k++){;}
    k2=k-1;


    /* determine indices of backbone atoms (N, Ca, C, O) of residue i: */
    /* first, initialize to -1 in case --for some reason-- an atom is absent: */
    kN=-1; kCa=-1; kC=-1; kOx=-1;
    for(atc=0,k=k1; atc<4 && k<=k2; k++){
      if(strcmp(pdb_model_all[k].type, "N")==0 && strcmp(pdb_model_all[k].loc, "" )==0){
        atc++; kN=k; continue;
      }
      if(strcmp(pdb_model_all[k].type, "C")==0 && strcmp(pdb_model_all[k].loc, "A")==0){
        atc++; kCa=k; continue;
      }
      if(strcmp(pdb_model_all[k].type, "C")==0 && strcmp(pdb_model_all[k].loc, "" )==0){
        atc++; kC=k; continue;
      }
      if(strcmp(pdb_model_all[k].type, "O")==0 &&
        (strcmp(pdb_model_all[k].loc , "" )==0 || strcmp(pdb_model_all[k].loc,"XT")==0)){
        atc++; kOx=k; continue;
      }
    }


    if(coords[p].code &  1){               /* (q_j, q_{j+1}, q_{j+2}) = (x,y,z) of N_i */
      delx = qdots[j+0]*deltat;
      dely = qdots[j+1]*deltat;
      delz = qdots[j+2]*deltat;
      for(k=k1;k<num_atoms_all;k++){
        if(pdb_model_all[k].chain[0] != chainid)
          break;                           /* finished current chain */
        pdb_model_all[k].x += delx;
        pdb_model_all[k].y += dely;
        pdb_model_all[k].z += delz;
      }
      j += 3;
    }


    if(coords[p].code &  2){               /* (q_j, q_{j+1}) = (lambda,theta) */

      rc1[0] = pdb_model_all[kN].x;  rc1[1] = pdb_model_all[kN].y;  rc1[2] = pdb_model_all[kN].z;
      e[0] = rc1[0]; e[1] = rc1[1]; e[2] = rc1[2];       /* e for lambda */
      rotmat(qdots[j]*deltat, e, rot);
      for(k=k1;k<num_atoms_all;k++){       /* q_j = lambda_i */
        if(pdb_model_all[k].chain[0] != chainid)
          break;                           /* finished current chain */
        if(k == kN) continue;
        tv[0]=pdb_model_all[k].x; tv[1]=pdb_model_all[k].y; tv[2]=pdb_model_all[k].z;
        matvec(rot, tv, rc1, tvp);
        pdb_model_all[k].x=tvp[0]; pdb_model_all[k].y=tvp[1]; pdb_model_all[k].z=tvp[2];
      }
      j++;

      rc2[0] = pdb_model_all[kCa].x; rc2[1] = pdb_model_all[kCa].y; rc2[2] = pdb_model_all[kCa].z;
      e[0] = rc1[1]*rc2[2]-rc1[2]*rc2[1];
      e[1] = rc1[2]*rc2[0]-rc1[0]*rc2[2];
      e[2] = rc1[0]*rc2[1]-rc1[1]*rc2[0];    /* e for theta*/
      rotmat(qdots[j]*deltat, e, rot);
      for(k=k1;k<num_atoms_all;k++){         /* q_j = theta_i */
        if(pdb_model_all[k].chain[0] != chainid)
          break;                             /* finished current chain */
        if(k == kN) continue;
        tv[0]=pdb_model_all[k].x; tv[1]=pdb_model_all[k].y; tv[2]=pdb_model_all[k].z;
        matvec(rot, tv, rc1, tvp);
        pdb_model_all[k].x=tvp[0]; pdb_model_all[k].y=tvp[1]; pdb_model_all[k].z=tvp[2];
      }
      j++;
    }


    rc1[0] = pdb_model_all[kN].x;  rc1[1] = pdb_model_all[kN].y;  rc1[2] = pdb_model_all[kN].z;
    rc2[0] = pdb_model_all[kCa].x; rc2[1] = pdb_model_all[kCa].y; rc2[2] = pdb_model_all[kCa].z;
    rc3[0] = pdb_model_all[kC].x;  rc3[1] = pdb_model_all[kC].y;  rc3[2] = pdb_model_all[kC].z;

    if(coords[p].code &  4){               /* q_j = phi_i */
      rc1[0] = pdb_model_all[kN].x;  rc1[1] = pdb_model_all[kN].y;  rc1[2] = pdb_model_all[kN].z;
      rc2[0] = pdb_model_all[kCa].x; rc2[1] = pdb_model_all[kCa].y; rc2[2] = pdb_model_all[kCa].z;
      e[0] = rc2[0]-rc1[0];
      e[1] = rc2[1]-rc1[1];
      e[2] = rc2[2]-rc1[2];                /* e for phi */
      rotmat(qdots[j]*deltat, e, rot);
      for(k=k1;k<num_atoms_all;k++){
        if(pdb_model_all[k].chain[0] != chainid)
          break;                           /* finished current chain */
        if(k==kN || k==kCa) continue;
        if(k<=k2 && pdb_model_all[k].type[strlen(pdb_model_all[k].type)-1]=='H'
                 && (strcmp(pdb_model_all[k].loc,  "N")==0 || strcmp(pdb_model_all[k].loc, "")==0)) continue;
        tv[0]=pdb_model_all[k].x; tv[1]=pdb_model_all[k].y; tv[2]=pdb_model_all[k].z;
        matvec(rot, tv, rc2, tvp);
        pdb_model_all[k].x=tvp[0]; pdb_model_all[k].y=tvp[1]; pdb_model_all[k].z=tvp[2];
      }
      j++;
    }

    if(coords[p].code &  8){               /* q_j  = psi_i */
      rc2[0] = pdb_model_all[kCa].x; rc2[1] = pdb_model_all[kCa].y; rc2[2] = pdb_model_all[kCa].z;
      rc3[0] = pdb_model_all[kC].x;  rc3[1] = pdb_model_all[kC].y;  rc3[2] = pdb_model_all[kC].z;
      e[0] = rc3[0]-rc2[0];
      e[1] = rc3[1]-rc2[1];
      e[2] = rc3[2]-rc2[2];                /* e for psi */
      rotmat(qdots[j]*deltat, e, rot);
      for(k=k1;k<num_atoms_all;k++){
        if(pdb_model_all[k].chain[0] != chainid)
          break;                           /* finished current chain */
        if(k>k2 || k==kOx){
          tv[0]=pdb_model_all[k].x; tv[1]=pdb_model_all[k].y; tv[2]=pdb_model_all[k].z;
          matvec(rot, tv, rc2, tvp);
          pdb_model_all[k].x=tvp[0]; pdb_model_all[k].y=tvp[1]; pdb_model_all[k].z=tvp[2];
        }
      }
      j++;
    }

    if(coords[p].code & 16){               /* q_j = chi_i */
      rc2[0] = pdb_model_all[kCa].x; rc2[1] = pdb_model_all[kCa].y; rc2[2] = pdb_model_all[kCa].z;

      /* index of Cb of residue i: */
      for(k=k1; k<=k2; k++){
        if(strcmp(pdb_model_all[k].type, "C")==0 && strcmp(pdb_model_all[k].loc, "B")==0){
          kCb=k; break;
        }
      }
      rc4[0] = pdb_model_all[kCb].x; rc4[1] = pdb_model_all[kCb].y; rc4[2] = pdb_model_all[kCb].z;

      e[0] = rc4[0] - rc2[0];
      e[1] = rc4[1] - rc2[1];
      e[2] = rc4[2] - rc2[2];              /* e for chi */
      rotmat(qdots[j]*deltat, e, rot);
      for(k=k1;k<=k2;k++){
        if(k==kN || k==kCa || k==kCb || k==kC || k==kOx) continue;
        if(pdb_model_all[k].type[strlen(pdb_model_all[k].type)-1]=='H' &&
           (strcmp(pdb_model_all[k].loc, "N")==0 || strcmp(pdb_model_all[k].loc, "")==0)) continue;
        if(pdb_model_all[k].type[strlen(pdb_model_all[k].type)-1]=='H' &&
           strcmp(pdb_model_all[k].loc,  "A")==0) continue;
        tv[0]=pdb_model_all[k].x; tv[1]=pdb_model_all[k].y; tv[2]=pdb_model_all[k].z;
        matvec(rot, tv, rc2, tvp);
        pdb_model_all[k].x=tvp[0]; pdb_model_all[k].y=tvp[1]; pdb_model_all[k].z=tvp[2];
      }
      j++;
    }
  }
  /*--------------------------------------------------------------------------------------*/

  istep++;
  nsslr++;
  nsslp++;
  curtime += deltat;    /* conformation just obtained corresponds to this time */
  smd += rmsdisp;

  goto newstep;

}


/*======================================================================*/
void optpdb(int num_atoms_all, PDB *pdb_model_all, char *bname, int iconf, char *file_name)
{
/* writes out a pdb model with the N, CA, C, O, CB atoms of every residue
   and calls SCATD or SCWRL to optimize side chains */
/* the SCRWL option was added on 9/29/2016 */

  int i, j, num_atoms_scpred;
  char scpred_cmd[1000];
  PDB *pdb_model_scpred;

#if SCPP == 1       /* use SCATD */
  char scpred_1[] = "~/Software/SideChain/SCATD -i temp.pdb -o ";
  char scpred_2[] = " -a ip -r radii/PARSE.siz > /dev/null";
#elif SCPP == 2     /* use SCWRL; option -h avoids hydrogens in the output */
  char scpred_1[] = "~/Software/SCWRL4/Scwrl4 -i temp.pdb -o ";
  char scpred_2[] = " -h > /dev/null";
#endif

  pdb_model_scpred = (PDB *) malloc(1*sizeof(PDB));
  ptr_check(pdb_model_scpred,34);

  j=0;   /* atom index in the 5-atom-per-residue model */

  /* build 5-atom-per-residue pdb model: */
  for(i=0;i<num_atoms_all;i++){
    if( (strcmp(pdb_model_all[i].type, "N")==0 && strcmp(pdb_model_all[i].loc, "" )==0) ||
        (strcmp(pdb_model_all[i].type, "C")==0 && strcmp(pdb_model_all[i].loc, "A")==0) ||
        (strcmp(pdb_model_all[i].type, "C")==0 && strcmp(pdb_model_all[i].loc, "" )==0) ||
        (strcmp(pdb_model_all[i].type, "O")==0 &&
                                              ( strcmp(pdb_model_all[i].loc, ""  )==0 ||
                                                strcmp(pdb_model_all[i].loc, "XT")==0   ) ) ||
        (strcmp(pdb_model_all[i].type, "C")==0 && strcmp(pdb_model_all[i].loc, "B")==0)   ){
      pdb_model_scpred = (PDB *) realloc(pdb_model_scpred, (j+1)*sizeof(PDB));
      ptr_check(pdb_model_scpred,341);
      pdb_model_scpred[j] = pdb_model_all[i];
      j++;
    }
  }
  num_atoms_scpred = j;

  write_pdb("temp.pdb", num_atoms_scpred, pdb_model_scpred);
  /* bname = name of input pdb file without the extension ".pdb" */
  sprintf(file_name, "%s.deform.%d.pdb", bname, iconf);
  sprintf(scpred_cmd, "%s%s%s", scpred_1, file_name, scpred_2);

  system(scpred_cmd);

  free(pdb_model_scpred);
}


/*======================================================================*/
int rtype(PDB *pdb_model, int i)
{
/* returns the residue type of pdb_model[i] */

  if(strcmp(pdb_model[i].res, "ALA")==0)
    return 1;
  else if(strcmp(pdb_model[i].res, "CYS")==0)
    return 2;
  else if(strcmp(pdb_model[i].res, "ASP")==0)
    return 3;
  else if(strcmp(pdb_model[i].res, "GLU")==0)
    return 4;
  else if(strcmp(pdb_model[i].res, "PHE")==0)
    return 5;
  else if(strcmp(pdb_model[i].res, "GLY")==0)
    return 6;
  else if(strcmp(pdb_model[i].res, "HIS")==0)
    return 7;
  else if(strcmp(pdb_model[i].res, "ILE")==0)
    return 8;
  else if(strcmp(pdb_model[i].res, "LYS")==0)
    return 9;
  else if(strcmp(pdb_model[i].res, "LEU")==0)
    return 10;
  else if(strcmp(pdb_model[i].res, "MET")==0)
    return 11;
  else if(strcmp(pdb_model[i].res, "ASN")==0)
    return 12;
  else if(strcmp(pdb_model[i].res, "PRO")==0)
    return 13;
  else if(strcmp(pdb_model[i].res, "GLN")==0)
    return 14;
  else if(strcmp(pdb_model[i].res, "ARG")==0)
    return 15;
  else if(strcmp(pdb_model[i].res, "SER")==0)
    return 16;
  else if(strcmp(pdb_model[i].res, "SEP")==0)
    return 16;
  else if(strcmp(pdb_model[i].res, "THR")==0)
    return 17;
  else if(strcmp(pdb_model[i].res, "TPO")==0)
    return 17;
  else if(strcmp(pdb_model[i].res, "VAL")==0)
    return 18;
  else if(strcmp(pdb_model[i].res, "TRP")==0)
    return 19;
  else if(strcmp(pdb_model[i].res, "TYR")==0)
    return 20;
  else{          /* unknown type */
    printf("Unknown residue type: '%s'. Stopping.\n",pdb_model[i].res);
    fprintf(stderr, "Unknown residue type. Stopping.\n");
    exit(5);
  }
}


/*====================================================================*/
int varind(fvars *tfv, int size, int resind, int varcode){
  /* returns the index of the (free) variable corresponding to residue index
     'resind' and the variable code 'varcode' by doing a lookup of the tfv array.
     -1 is returned if the specified resind/varcode is not a free var. */

  int j, ju, jm, jl, jo;

  jl=-1;
  ju=size;

  while(ju-jl > 1){
    jm = (ju+jl)>>2;
    if(resind >= tfv[jm].iko)
      jl=jm;
    else
      ju=jm;
  }

  for(j=ju-1;j>=0;j--)
    if(tfv[j].iko < resind) break;
  j++;

  for(jo=-1; j<size && tfv[j].iko==resind; j++)
    if(tfv[j].code==varcode){jo=j; break;}

  return jo;
}


/*====================================================================*/
/*********************************************************************
*                          GETQDOTS()                                *
**********************************************************************
*                                                                    *
*  Function to solve the linear system of equations for the qdot's   *
*                                                                    *
**********************************************************************
*                                                                    *
*  Needs Lapack 3.0 and Blast libraries.                             *
*                                                                    *
*********************************************************************/


/*OUTPUT **********************

qdots: time derivatives of the free generalized coords.

getqdots returns 0 if the calculation was successful,
         or a non-zero value otherwise.

*****************************/

int getqdots(PDB *pdb_model, int num_atoms, int num_res, double dc2, fvars *tfv,
             int size, int mnv, int nco1, int nco2, fvars *coords, int varlen,
             sres *statres, int srlen, respair *dconst, int dclen, int nchains,
             chdata *chinfo, tri *props, trd *wop, twovec *erb, trd *der,
             double *qdots, double damp_scale, double drag_scale){

  double currdist, r[3], rsqu, mfa, sum, temp, c1, c2, c3;
  double bk, mdx, mdy, mdz, prod, prod1, difer;
  int b, c, h, i, i1, i2, ian, j, k, kp, l, lp, k1, k2, k1p, kb, kc, m, n, sizex, p;
  int lepoc, nepoc, R, S, xi, eta, sigma, vv;
  int mu, nu, ich, ich1, ich2, ivar, nat2, nat3, nat4, siz2, siz3, sizx2, sizx3;
  int k1R, k1Rp, k1S, k1Sp;
  long index, ind2, ind3, ll, isi, jsi, ks, ls, kl, lk, nesi, nsix, msi, nsi;
  long ncnv, ch1nv, ch1ncnv, ch2nv, ch2ncnv, isize, k1Rna, k1Rpna;
  double mtot, mb1, mb2, mpr, mta, ri, rj, nr2;
  double e[3], y[3], v[3], w[3], wcv[3], vw[3], cosg, sk;
  double *vs[6];
  double *Aop, rc1[3], rc2[3], rc3[3], Bnorm, Fnorm, Gnorm;

  double *Y[10], *Z[10];     /* i.e., *(Y[m]), so Y is an array of pointers to double */
  
  double zcar[21], zsph[21], zphi[21], zpsi[21], zchi[21];

  double L[21];     /* not really needed in this version */
  double RR[21];    /* not really needed in this version */

  double *carall[21];     /* see description below, in the malloc statement */

  double *allvars[21][6];   /* *allvars[m][0]    not used */
                            /* *allvars[m][1] -> cartesian */
                            /* *allvars[m][2] -> spherical */
                            /* *allvars[m][3] -> phi */
                            /* *allvars[m][4] -> psi */
                            /* *allvars[m][5] -> chi */

  char chainid;

  FILE *fca;    /* for testing printouts */

  /* Lapack subroutine arguments */
  char uplo;
  int  info, lda, ldb, nrhs, *ipiv, lwork;
  double *qm, *work;


  for(m=0;m<6;m++){
    vs[m]  = (double *) malloc(size * sizeof(double));
    ptr_check(vs[m],14);
  }

  /* order of system matrix: */
  sizex=size+nco1+nco2;

  Aop  = (double *) malloc(sizex*sizex * sizeof(double));
  ptr_check(Aop,15);
  for(ll=0; ll<sizex*sizex; ll++)
    Aop[ll]=0.0;

  /* no used with the LU routine:
  ipiv = (int *) malloc(sizex * sizeof(int));
  ptr_check(ipiv,16);
  */

  qm = (double *) malloc(sizex * sizeof(double));
  ptr_check(qm,17);

  for(m=0;m<10;m++){
    Y[m]  = (double *) malloc(num_atoms * sizeof(double));
    ptr_check(Y[m],19);
  }

  for(m=0;m<10;m++){
    Z[m]  = (double *) malloc(size * sizeof(double));
    ptr_check(Z[m],20);
  }


/* array to temporarily store the R_ij for i= each of the 5 types of variables
   (given by ivar) of a chain and j= all vars of a different chain */
  for(ivar=1;ivar<=5;ivar++)
    for(m=0;m<21;m++){
      allvars[m][ivar] = (double *) malloc(mnv * sizeof(double));
      ptr_check(allvars[m][ivar],23);
    }

/* array to store the R_ij for pairs of vars where q_i = c(ich1) and
   q_j = all vars of chains ich2!=ich1  (free and fixed),
   for use in the first part of the scanning of pairs of vars belonging to the same chain.
   R_{c(ich1),q_j} -> carall[m][ich1*nchains*mnv + ich2*mnv + ivar]
   where ivar is given by:
                          0  for q_j = c(ich2)
                          1  for q_j = s(ich2)
   3*(R-chinfo[ich2].ikf)+2  for q_j = phi_R(ich2)
   3*(R-chinfo[ich2].ikf)+3  for q_j = psi_R(ich2)
   3*(R-chinfo[ich2].ikf)+4  for q_j = chi_R(ich2)
*/
  for(m=0;m<21;m++){
    carall[m] = (double *) malloc(nchains*nchains*mnv * sizeof(double));
    ptr_check(carall[m],24);
  }


  /* store the 6-dim factors (e_j, e_j x r_{beta(j)}): */
  for(j=0;j<size;j++){
    vs[0][j] = erb[j].e.x;
    vs[1][j] = erb[j].e.y;
    vs[2][j] = erb[j].e.z;
    vs[3][j] = erb[j].e.y * erb[j].b.z - erb[j].e.z * erb[j].b.y;
    vs[4][j] = erb[j].e.z * erb[j].b.x - erb[j].e.x * erb[j].b.z;
    vs[5][j] = erb[j].e.x * erb[j].b.y - erb[j].e.y * erb[j].b.x;
  }


  /***************************************************************************/
  /***** build matrix Bij (only upper triangular part)
             using recursive approach (adapted from Go's papers) *****/

  for(k=0;k<num_atoms;k++){
    bk = drag_scale * 0.1 * pdb_model[k].occupancy;
    Y[6][k] = bk*pdb_model[k].x;
    Y[7][k] = bk*pdb_model[k].y;
    Y[8][k] = bk*pdb_model[k].z;
    Y[9][k] = bk;
    Y[0][k] = Y[6][k]*pdb_model[k].x;
    Y[1][k] = Y[6][k]*pdb_model[k].y;
    Y[2][k] = Y[6][k]*pdb_model[k].z;
    Y[3][k] = Y[7][k]*pdb_model[k].y;
    Y[4][k] = Y[7][k]*pdb_model[k].z;
    Y[5][k] = Y[8][k]*pdb_model[k].z;
  }

  /* scan residues (backwards) to perform recursion: */
  /* j = index of next free variable */
  /* p = index in the 'coords' array */

  for(R=num_res-1,j=size-1,p=varlen-1;R>=0;R--){

    k1 = props[R].k1;          /* index of 1st pseudoatom of residue R */

    if(props[R].epoc != 2)
      k1p = props[R+1].k1;     /* index of 1st pseudoatom of residue R+1 */

    if(props[R].rt == 1 || props[R].rt == 6)    /* chi_R */
      for(m=0;m<10;m++)
        zchi[m] = 0.0;
    else
      for(m=0;m<10;m++)
        zchi[m] = Y[m][k1+4];

    if(props[R].epoc == 2)                      /* psi_R */
      for(m=0;m<10;m++)
        zpsi[m] = 0.0;
    else
      for(m=0;m<10;m++)
        zpsi[m] = zphi[m] + Y[m][k1p] + Y[m][k1p+1];

    for(m=0;m<10;m++)                           /* phi_R */
      zphi[m] = zpsi[m] + zchi[m] + Y[m][k1+2];
    if(props[R].rt != 6)
      for(m=0;m<10;m++)
        zphi[m] += Y[m][k1+3];

    if(props[R].epoc == 1){
      for(m=0;m<10;m++)                         /* sph */
        zsph[m] = zphi[m] + Y[m][k1+1];

      for(m=0;m<10;m++)                         /* car */
        zcar[m] = zsph[m] + Y[m][k1];
    }

    /* fill in the Z array for the free vars of this residue: */
    if(R==coords[p].iko){
      if(coords[p].code & 16){                  /* chi_R */
        for(m=0;m<10;m++) Z[m][j] = zchi[m];
        j--;
      }
      if(coords[p].code &  8){                  /* psi_R */
        for(m=0;m<10;m++) Z[m][j] = zpsi[m];
        j--;
      }
      if(coords[p].code &  4){                  /* phi_R */
        for(m=0;m<10;m++) Z[m][j] = zphi[m];
        j--;
      }
      if(coords[p].code &  2){                  /* sph */
        for(m=0;m<10;m++){
          Z[m][j]   = zsph[m];
          Z[m][j-1] = zsph[m];
        }
        j -= 2;
      }
      if(coords[p].code &  1){                  /* car */
        for(m=0;m<10;m++){
          Z[m][j]   = zcar[m];
          Z[m][j-1] = zcar[m];
          Z[m][j-2] = zcar[m];
        }
        j -= 3;
      }
      p--;
    }

  }


  /* scan pairs of free vars and compute Bij */
  for(i=0,isi=0; i<size; i++,isi+=sizex){
    for(j=i;j<size;j++){
      if( tfv[j].cho  != tfv[i].cho ||
         (tfv[j].code == 30 && tfv[i].code == 50 && tfv[j].iko >  tfv[i].iko) ||
         (tfv[i].code == 40 && tfv[j].code == 50 && tfv[i].iko == tfv[j].iko) ||
         (tfv[j].code == 40 && tfv[i].code == 50 && tfv[j].iko >  tfv[i].iko) ||
         (tfv[j].code == 50 && tfv[i].code == 50 && tfv[j].iko != tfv[i].iko) ){
        /* (in this case, M_i inter M_j = empty, hence B_{ij} = 0) */
        Aop[isi+j] = 0.0;
      }
      else{
        /* (in this case, M_j is a subset of M_i, hence Z_{ij} = Z_j) */
        if(tfv[i].code >= 20 && tfv[j].code >= 20){    /* q_i and q_j angular */
          /* first term: */
          c3 = erb[i].e.x * erb[j].e.x + erb[i].e.y * erb[j].e.y + erb[i].e.z * erb[j].e.z;
          sum  = c3 * (Z[0][j]+Z[3][j]+Z[5][j]);

          /* second term: */
          sum -= erb[i].e.x * erb[j].e.x * Z[0][j];
          sum -= erb[i].e.y * erb[j].e.y * Z[3][j];
          sum -= erb[i].e.z * erb[j].e.z * Z[5][j];
          sum -= (erb[i].e.x * erb[j].e.y + erb[i].e.y * erb[j].e.x)*Z[1][j];
          sum -= (erb[i].e.x * erb[j].e.z + erb[i].e.z * erb[j].e.x)*Z[2][j];
          sum -= (erb[i].e.y * erb[j].e.z + erb[i].e.z * erb[j].e.y)*Z[4][j];

          /* third term: */
          v[0] = vs[3][i]; w[0] = vs[3][j];
          v[1] = vs[4][i]; w[1] = vs[4][j];
          v[2] = vs[5][i]; w[2] = vs[5][j];

          sum += (v[0]*w[0]+v[1]*w[1]+v[2]*w[2]) * Z[9][j];

          /* fourth term: */
          c1 = erb[i].e.x * erb[j].b.x + erb[i].e.y * erb[j].b.y + erb[i].e.z * erb[j].b.z;
          c2 = erb[j].e.x * erb[i].b.x + erb[j].e.y * erb[i].b.y + erb[j].e.z * erb[i].b.z;
          v[0] = c1 * erb[j].e.x + c2 * erb[i].e.x - c3 * (erb[i].b.x + erb[j].b.x);
          v[1] = c1 * erb[j].e.y + c2 * erb[i].e.y - c3 * (erb[i].b.y + erb[j].b.y);
          v[2] = c1 * erb[j].e.z + c2 * erb[i].e.z - c3 * (erb[i].b.z + erb[j].b.z);
          sum += v[0] * Z[6][j] + v[1] * Z[7][j] + v[2] * Z[8][j];

          /* Bij: */
          Aop[isi+j] = sum;
        }
        if(tfv[i].code >= 20 && tfv[j].code < 20){     /* q_i angular, q_j linear */
        /* (this case should never occur) */
          xi = tfv[j].code - 10;   /* 0, 1, 2, for x, y, z, respectively */
          if(xi==0){
            v[0] = erb[i].e.y * erb[i].b.z - erb[i].e.z * erb[i].b.y;
            w[0] = erb[i].e.y * Z[8][j]    - erb[i].e.z * Z[7][j];
          }
          if(xi==1){
            v[1] = erb[i].e.z * erb[i].b.x - erb[i].e.x * erb[i].b.z;
            w[1] = erb[i].e.z * Z[6][j]    - erb[i].e.x * Z[8][j];
          }
          if(xi==2){
            v[2] = erb[i].e.x * erb[i].b.y - erb[i].e.y * erb[i].b.x;
            w[2] = erb[i].e.x * Z[7][j]    - erb[i].e.y * Z[6][j];
          }
          Aop[isi+j] = w[xi] - v[xi] * Z[9][j];
        }
        if(tfv[i].code < 20 && tfv[j].code >= 20){     /* q_i linear, q_j angular */
          xi = tfv[i].code - 10;   /* 0, 1, 2, for x, y, z, respectively */
          if(xi==0){
            v[0] = erb[j].e.y * erb[j].b.z - erb[j].e.z * erb[j].b.y;
            w[0] = erb[j].e.y * Z[8][j]    - erb[j].e.z * Z[7][j];
          }
          if(xi==1){
            v[1] = erb[j].e.z * erb[j].b.x - erb[j].e.x * erb[j].b.z;
            w[1] = erb[j].e.z * Z[6][j]    - erb[j].e.x * Z[8][j];
          }
          if(xi==2){
            v[2] = erb[j].e.x * erb[j].b.y - erb[j].e.y * erb[j].b.x;
            w[2] = erb[j].e.x * Z[7][j]    - erb[j].e.y * Z[6][j];
          }
          Aop[isi+j] = w[xi] - v[xi] * Z[9][j];
        }
        if(tfv[i].code < 20 && tfv[j].code < 20){      /* q_i and q_j linear */
          if(tfv[i].code == tfv[j].code)
            Aop[isi+j] = Z[9][j];
          else
            Aop[isi+j] = 0.0;
        }
      }
    }
  }


/* test:
Bnorm = 0.0;
  for(i=0,isi=0; i<size; i++,isi+=sizex){
    for(j=i;j<size;j++){
      Bnorm += pow(Aop[isi+j],2);
    }
  }
Bnorm = sqrt(Bnorm);
fprintf(stderr,"\nNorm of matrix B = %e    ",Bnorm);
*/


/***** add matrix Vij (only upper triangular part)
       using recursive approach *****/

  /* SCAN PAIRS OF VARIABLES TO PERFORM RECURSION */

  nat2 = num_atoms*2;
  nat3 = num_atoms*3;
  nat4 = num_atoms*4;

  siz2 = size*2;
  siz3 = size*3;
  
  sizx2 = sizex*2;
  sizx3 = sizex*3;

  ncnv= nchains*mnv;

  /* FIRST, scan pairs of DIFFERENT chains */

  for(ich1=0;ich1<nchains-1;ich1++){

    ch1nv = ich1*mnv;
    ch1ncnv = ch1nv*nchains;

    for(ich2=ich1+1;ich2<nchains;ich2++){

      ch2nv = ich2*mnv;
      ch2ncnv = ch2nv*nchains;

      i = chinfo[ich1].ivl;
      isize = i*size;
      isi = i*sizex;

      for(R=chinfo[ich1].ikl; R>=chinfo[ich1].ikf; R--){

        ivar = 3*(R-chinfo[ich1].ikf);

        k1R = props[R].k1;          /* index of 1st pseudoatom of residue R */
        if(props[R].epoc != 2)
          k1Rp = props[R+1].k1;     /* index of 1st pseudoatom of residue R+1 */

        k1Rna  = k1R*num_atoms;
        k1Rpna = k1Rp*num_atoms;


    /* chi_R (ich1) ******************************************************/

        vv=0;
        j=chinfo[ich2].ivl;

        for(S=chinfo[ich2].ikl; S>=chinfo[ich2].ikf; S--){

          k1S = props[S].k1;          /* index of 1st pseudoatom of residue S */
          if(props[S].epoc != 2)
            k1Sp = props[S+1].k1;     /* index of 1st pseudoatom of residue S+1 */

          if(props[R].rt == 1 || props[R].rt == 6 ||
             props[S].rt == 1 || props[S].rt == 6)    /* chi_S (ich2) */
            for(m=0;m<21;m++)
              zchi[m] = 0.0;
          else
            calcL(k1R+4, k1S+4, 0, dc2, pdb_model, damp_scale, zchi);
            /* for(m=0;m<21;m++) zchi[m] = L[m][k1Rna+nat4+k1S+4]; */
          for(m=0;m<21;m++) allvars[m][5][vv] = zchi[m];
          vv++;

          /* if current vars are free, compute V_ij and add it to Aop_ij */
          if(i>=0 && tfv[i].iko==R && tfv[i].code==50 &&
             j>=0 && tfv[j].iko==S && tfv[j].code==50){
            calcV(i, j, isi, vs, tfv, zchi, Aop);
            /* for(m=0;m<21;m++) RR[m][isize+j] = zchi[m]; */
            j--;
          }


          if(props[R].rt == 1 || props[R].rt == 6 ||
             props[S].epoc == 2)                         /* psi_S (ich2) */
            for(m=0;m<21;m++)
              zpsi[m] = 0.0;
          else{
            for(m=0;m<21;m++) zpsi[m] = zphi[m];
            calcL(k1R+4, k1Sp,   1, dc2, pdb_model, damp_scale, zpsi);
            calcL(k1R+4, k1Sp+1, 1, dc2, pdb_model, damp_scale, zpsi);
            /* for(m=0;m<21;m++) zpsi[m] = zphi[m] + L[m][k1Rna+nat4+k1Sp] + L[m][k1Rna+nat4+k1Sp+1]; */
          }

          for(m=0;m<21;m++) allvars[m][5][vv] = zpsi[m];
          vv++;

          /* if current vars are free, compute V_ij and add it to Aop_ij */
          if(i>=0 && tfv[i].iko==R && tfv[i].code==50 &&
             j>=0 && tfv[j].iko==S && tfv[j].code==40){
            calcV(i, j, isi, vs, tfv, zpsi, Aop);
            /* for(m=0;m<21;m++) RR[m][isize+j] = zpsi[m]; */
            j--;
          }


          if(props[R].rt == 1 || props[R].rt == 6)        /* phi_S (ich2) */
            for(m=0;m<21;m++)
              zphi[m] = 0.0;
          else{
            for(m=0;m<21;m++) zphi[m] = zpsi[m] + zchi[m];
            calcL(k1R+4, k1S+2, 1, dc2, pdb_model, damp_scale, zphi);
            /* for(m=0;m<21;m++) zphi[m] = zpsi[m] + zchi[m] + L[m][k1Rna+nat4+k1S+2]; */
            if(props[S].rt != 6)
              calcL(k1R+4, k1S+3, 1, dc2, pdb_model, damp_scale, zphi);
              /* for(m=0;m<21;m++) zphi[m] += L[m][k1Rna+nat4+k1S+3]; */
          }

          for(m=0;m<21;m++) allvars[m][5][vv] = zphi[m];
          vv++;

          /* if current vars are free, compute V_ij and add it to Aop_ij */
          if(i>=0 && tfv[i].iko==R && tfv[i].code==50 &&
             j>=0 && tfv[j].iko==S && tfv[j].code==30){
            calcV(i, j, isi, vs, tfv, zphi, Aop);
            /* for(m=0;m<21;m++) RR[m][isize+j] = zphi[m]; */
            j--;
          }


          if(props[S].epoc != 1) continue;    /* pass to next S */

          /* otherwise...                          sph_S(ich2)  &  car_S(ich2) */

          if(props[R].rt == 1 || props[R].rt == 6)
            for(m=0;m<21;m++){zsph[m] = 0.0; zcar[m] = 0.0;}
          else{
            for(m=0;m<21;m++) zsph[m] = zphi[m];
            calcL(k1R+4, k1S+1, 1, dc2, pdb_model, damp_scale, zsph);
            for(m=0;m<21;m++) zcar[m] = zsph[m];
            calcL(k1R+4, k1S, 1, dc2, pdb_model, damp_scale, zcar);
            /*for(m=0;m<21;m++){
              zsph[m] = zphi[m] + L[m][k1Rna+nat4+k1S+1];
              zcar[m] = zsph[m] + L[m][k1Rna+nat4+k1S];
            }*/
          }

          for(m=0;m<21;m++){
            allvars[m][5][vv]     = zsph[m];
            allvars[m][5][vv+1]   = zcar[m];
            carall[m][ch2ncnv+ch1nv+ivar+4] = zcar[m];
          }
          vv += 2;

          /* if current vars are free, compute V_ij and add it to Aop_ij */
          if(i>=0 && tfv[i].iko==R && tfv[i].code==50 &&
             j>=0 && tfv[j].iko==S && tfv[j].code/10==2){
            calcV(i, j,   isi, vs, tfv, zsph, Aop);
            calcV(i, j-1, isi, vs, tfv, zsph, Aop);
            /* for(m=0;m<21;m++){
                 RR[m][isize+j]   = zsph[m];
                 RR[m][isize+j-1] = zsph[m];
               } */
            j -= 2;
          }

          /* if current vars are free, compute V_ij and add it to Aop_ij */
          if(i>=0 && tfv[i].iko==R && tfv[i].code==50 &&
             j>=0 && tfv[j].iko==S && tfv[j].code/10==1){
            calcV(i, j,   isi, vs, tfv, zcar, Aop);
            calcV(i, j-1, isi, vs, tfv, zcar, Aop);
            calcV(i, j-2, isi, vs, tfv, zcar, Aop);
            /* for(m=0;m<21;m++){
                 RR[m][isize+j]   = zcar[m];
                 RR[m][isize+j-1] = zcar[m];
                 RR[m][isize+j-2] = zcar[m];
               } */
            j -= 3;
          }

        }      /* end of loop over S */

        if(i>=0 && tfv[i].iko==R && tfv[i].code==50){
          i--;
          isize -= size;
          isi -= sizex;
        }

   /* psi_R (ich1) ******************************************************/

        vv=0;
        j=chinfo[ich2].ivl;

        for(S=chinfo[ich2].ikl; S>=chinfo[ich2].ikf; S--){

          k1S = props[S].k1;          /* index of 1st pseudoatom of residue S */
          if(props[S].epoc != 2)
            k1Sp = props[S+1].k1;     /* index of 1st pseudoatom of residue S+1 */

          if(props[R].epoc == 2 ||
             props[S].rt == 1 || props[S].rt == 6)       /* chi_S (ich2) */
            for(m=0;m<21;m++)
              zchi[m] = 0.0;
          else{
            for(m=0;m<21;m++) zchi[m] = allvars[m][3][vv];
            calcL(k1Rp,   k1S+4, 1, dc2, pdb_model, damp_scale, zchi);
            calcL(k1Rp+1, k1S+4, 1, dc2, pdb_model, damp_scale, zchi);
            /* for(m=0;m<21;m++)
                 zchi[m] = allvars[m][3][vv] + L[m][k1Rpna+k1S+4]
                                             + L[m][k1Rpna+num_atoms+k1S+4]; */
          }

          for(m=0;m<21;m++) allvars[m][4][vv] = zchi[m];
          vv++;

          /* if current vars are free, compute V_ij and add it to Aop_ij */
          if(i>=0 && tfv[i].iko==R && tfv[i].code==40 &&
             j>=0 && tfv[j].iko==S && tfv[j].code==50){
            calcV(i, j, isi, vs, tfv, zchi, Aop);
            /* for(m=0;m<21;m++) RR[m][isize+j] = zchi[m]; */
            j--;
          }


          if(props[R].epoc == 2 || props[S].epoc == 2)     /* psi_S (ich2) */
            for(m=0;m<21;m++)
              zpsi[m] = 0.0;
          else{
            for(m=0;m<21;m++) zpsi[m] = zphi[m] + allvars[m][3][vv] - allvars[m][3][vv-2];
            calcL(k1Rp,   k1Sp,   1, dc2, pdb_model, damp_scale, zpsi);
            calcL(k1Rp,   k1Sp+1, 1, dc2, pdb_model, damp_scale, zpsi);
            calcL(k1Rp+1, k1Sp,   1, dc2, pdb_model, damp_scale, zpsi);
            calcL(k1Rp+1, k1Sp+1, 1, dc2, pdb_model, damp_scale, zpsi);
            /* for(m=0;m<21;m++)
                 zpsi[m] = zphi[m] + allvars[m][3][vv] - allvars[m][3][vv-2]
                                + L[m][k1Rpna+k1Sp] + L[m][k1Rpna+k1Sp+1]
                                + L[m][k1Rpna+num_atoms+k1Sp]
                                + L[m][k1Rpna+num_atoms+k1Sp+1]; */
          }

          for(m=0;m<21;m++) allvars[m][4][vv] = zpsi[m];
          vv++;

          /* if current vars are free, compute V_ij and add it to Aop_ij */
          if(i>=0 && tfv[i].iko==R && tfv[i].code==40 &&
             j>=0 && tfv[j].iko==S && tfv[j].code==40){
            calcV(i, j, isi, vs, tfv, zpsi, Aop);
            /* for(m=0;m<21;m++) RR[m][isize+j] = zpsi[m]; */
            j--;
          }

          if(props[R].epoc == 2)                            /* phi_S (ich2) */
            for(m=0;m<21;m++)
              zphi[m] = 0.0;
          else{
            for(m=0;m<21;m++) zphi[m] = zpsi[m] + zchi[m] + allvars[m][3][vv]
                                    - allvars[m][3][vv-1] - allvars[m][3][vv-2];
            calcL(k1Rp,   k1S+2, 1, dc2, pdb_model, damp_scale, zphi);
            calcL(k1Rp+1, k1S+2, 1, dc2, pdb_model, damp_scale, zphi);
            if(props[S].rt != 6){
              calcL(k1Rp,   k1S+3, 1, dc2, pdb_model, damp_scale, zphi);
              calcL(k1Rp+1, k1S+3, 1, dc2, pdb_model, damp_scale, zphi);
            }
            /* for(m=0;m<21;m++)
                 zphi[m] = zpsi[m] + zchi[m] + allvars[m][3][vv]
                           - allvars[m][3][vv-1] - allvars[m][3][vv-2]
                           + L[m][k1Rpna+k1S+2] + L[m][k1Rpna+num_atoms+k1S+2];
            if(props[S].rt != 6)
              for(m=0;m<21;m++)
                zphi[m] += L[m][k1Rpna+k1S+3] + L[m][k1Rpna+num_atoms+k1S+3]; */
          }

          for(m=0;m<21;m++) allvars[m][4][vv] = zphi[m];
          vv++;

          /* if current vars are free, compute V_ij and add it to Aop_ij */
          if(i>=0 && tfv[i].iko==R && tfv[i].code==40 &&
             j>=0 && tfv[j].iko==S && tfv[j].code==30){
            calcV(i, j, isi, vs, tfv, zphi, Aop);
            /* for(m=0;m<21;m++) RR[m][isize+j] = zphi[m]; */
            j--;
          }

          if(props[S].epoc != 1) continue;    /* pass to next S */

          if(props[R].epoc == 2)                          /* sph_S(ich2) */
            for(m=0;m<21;m++) zsph[m] = 0.0;
          else{
            for(m=0;m<21;m++) zsph[m] = zphi[m] + allvars[m][3][vv] - allvars[m][3][vv-1];
            calcL(k1Rp,   k1S+1, 1, dc2, pdb_model, damp_scale, zsph);
            calcL(k1Rp+1, k1S+1, 1, dc2, pdb_model, damp_scale, zsph);
            /* for(m=0;m<21;m++)
                 zsph[m] = zphi[m] + allvars[m][3][vv] - allvars[m][3][vv-1]
                         + L[m][k1Rpna+k1S+1] + L[m][k1Rpna+num_atoms+k1S+1]; */
          }

          for(m=0;m<21;m++) allvars[m][4][vv] = zsph[m];
          vv++;

          /* if current vars are free, compute V_ij and add it to Aop_ij */
          if(i>=0 && tfv[i].iko==R && tfv[i].code==40 &&
             j>=0 && tfv[j].iko==S && tfv[j].code/10==2){
            calcV(i, j,   isi, vs, tfv, zsph, Aop);
            calcV(i, j-1, isi, vs, tfv, zsph, Aop);
            /* for(m=0;m<21;m++){
                 RR[m][isize+j]   = zsph[m];
                 RR[m][isize+j-1] = zsph[m];
               } */
            j -= 2;
          }

          if(props[R].epoc == 2)                         /* car_S(ich2) */
            for(m=0;m<21;m++) zcar[m] = 0.0;
          else{
            for(m=0;m<21;m++) zcar[m] = zsph[m] + allvars[m][3][vv] - allvars[m][3][vv-1];
            calcL(k1Rp,   k1S, 1, dc2, pdb_model, damp_scale, zcar);
            calcL(k1Rp+1, k1S, 1, dc2, pdb_model, damp_scale, zcar);
            /* for(m=0;m<21;m++)
                 zcar[m] = zsph[m] + allvars[m][3][vv] - allvars[m][3][vv-1]
                          + L[m][k1Rpna+k1S] + L[m][k1Rpna+num_atoms+k1S]; */
          }

          for(m=0;m<21;m++){
            allvars[m][4][vv] = zcar[m];
             carall[m][ch2ncnv+ch1nv+ivar+3] = zcar[m];
          }
          vv++;

          /* if current vars are free, compute V_ij and add it to Aop_ij */
          if(i>=0 && tfv[i].iko==R && tfv[i].code==40 &&
             j>=0 && tfv[j].iko==S && tfv[j].code/10==1){
            calcV(i, j,   isi, vs, tfv, zcar, Aop);
            calcV(i, j-1, isi, vs, tfv, zcar, Aop);
            calcV(i, j-2, isi, vs, tfv, zcar, Aop);
            /* for(m=0;m<21;m++){
                 RR[m][isize+j]   = zcar[m];
                 RR[m][isize+j-1] = zcar[m];
                 RR[m][isize+j-2] = zcar[m];
               } */
            j -= 3;
          }

        }    /* end of loop over S */

        if(i>=0 && tfv[i].iko==R && tfv[i].code==40){
          i--;
          isize -= size;
          isi -= sizex;
        }

    /* phi_R (ich1) ******************************************************/

        vv=0;
        j=chinfo[ich2].ivl;

        for(S=chinfo[ich2].ikl; S>=chinfo[ich2].ikf; S--){

          k1S = props[S].k1;          /* index of 1st pseudoatom of residue S */
          if(props[S].epoc != 2)
            k1Sp = props[S+1].k1;     /* index of 1st pseudoatom of residue S+1 */

          if(props[S].rt == 1 || props[S].rt == 6)       /* chi_S (ich2) */
            for(m=0;m<21;m++)
              zchi[m] = 0.0;
          else{
            for(m=0;m<21;m++) zchi[m] = allvars[m][4][vv] + allvars[m][5][vv];
            calcL(k1R+2, k1S+4, 1, dc2, pdb_model, damp_scale, zchi);
            if(props[R].rt != 6)
              calcL(k1R+3, k1S+4, 1, dc2, pdb_model, damp_scale, zchi);
            /* for(m=0;m<21;m++)
                zchi[m] = allvars[m][4][vv] + allvars[m][5][vv]
                          + L[m][k1Rna+nat2+k1S+4];
            if(props[R].rt != 6)
              for(m=0;m<21;m++)
                zchi[m] += L[m][k1Rna+nat3+k1S+4]; */
          }

          for(m=0;m<21;m++) allvars[m][3][vv] = zchi[m];
          vv++;

          /* if current vars are free, compute V_ij and add it to Aop_ij */
          if(i>=0 && tfv[i].iko==R && tfv[i].code==30 &&
             j>=0 && tfv[j].iko==S && tfv[j].code==50){
            calcV(i, j, isi, vs, tfv, zchi, Aop);
            /* for(m=0;m<21;m++) RR[m][isize+j] = zchi[m]; */
            j--;
          }


          if(props[S].epoc == 2)                           /* psi_S (ich2) */
            for(m=0;m<21;m++)
              zpsi[m] = 0.0;
          else{
            for(m=0;m<21;m++) zpsi[m] = zphi[m] + allvars[m][4][vv]   + allvars[m][5][vv]
                                                - allvars[m][4][vv-2] - allvars[m][5][vv-2];
            calcL(k1R+2, k1Sp,   1, dc2, pdb_model, damp_scale, zpsi);
            calcL(k1R+2, k1Sp+1, 1, dc2, pdb_model, damp_scale, zpsi);
            if(props[R].rt != 6){
              calcL(k1R+3, k1Sp,   1, dc2, pdb_model, damp_scale, zpsi);
              calcL(k1R+3, k1Sp+1, 1, dc2, pdb_model, damp_scale, zpsi);
            }
            /* for(m=0;m<21;m++)
              zpsi[m] = zphi[m] + allvars[m][4][vv] + allvars[m][5][vv]
                        - allvars[m][4][vv-2] - allvars[m][5][vv-2]
                        + L[m][k1Rna+nat2+k1Sp] + L[m][k1Rna+nat2+k1Sp+1];
            if(props[R].rt != 6)
              for(m=0;m<21;m++)
                zpsi[m] += L[m][k1Rna+nat3+k1Sp] + L[m][k1Rna+nat3+k1Sp+1]; */
          }

          for(m=0;m<21;m++) allvars[m][3][vv] = zpsi[m];
          vv++;

          /* if current vars are free, compute V_ij and add it to Aop_ij */
          if(i>=0 && tfv[i].iko==R && tfv[i].code==30 &&
             j>=0 && tfv[j].iko==S && tfv[j].code==40){
            calcV(i, j, isi, vs, tfv, zpsi, Aop);
            /* for(m=0;m<21;m++) RR[m][isize+j] = zpsi[m]; */
            j--;
          }


          for(m=0;m<21;m++)                                  /* phi_S (ich2) */
            zphi[m] = zpsi[m] + zchi[m]
                      + allvars[m][4][vv] - allvars[m][4][vv-1] - allvars[m][4][vv-2]
                      + allvars[m][5][vv] - allvars[m][5][vv-1] - allvars[m][5][vv-2];
          calcL(k1R+2, k1S+2, 1, dc2, pdb_model, damp_scale, zphi);
          if(props[S].rt != 6)
            calcL(k1R+2, k1S+3, 1, dc2, pdb_model, damp_scale, zphi);
          if(props[R].rt != 6)
            calcL(k1R+3, k1S+2, 1, dc2, pdb_model, damp_scale, zphi);
          if(props[R].rt != 6 && props[S].rt != 6)
            calcL(k1R+3, k1S+3, 1, dc2, pdb_model, damp_scale, zphi);
          /* for(m=0;m<21;m++)
              zphi[m] = zpsi[m] + zchi[m]
                       + allvars[m][4][vv] - allvars[m][4][vv-1] - allvars[m][4][vv-2]
                       + allvars[m][5][vv] - allvars[m][5][vv-1] - allvars[m][5][vv-2]
                       + L[m][k1Rna+nat2+k1S+2];
          if(props[S].rt != 6)
            for(m=0;m<21;m++)
              zphi[m] += L[m][k1Rna+nat2+k1S+3];
          if(props[R].rt != 6)
            for(m=0;m<21;m++)
              zphi[m] += L[m][k1Rna+nat3+k1S+2];
          if(props[R].rt != 6 && props[S].rt != 6)
            for(m=0;m<21;m++)
              zphi[m] += L[m][k1Rna+nat3+k1S+3];
          if(props[S].rt != 1 && props[S].rt != 6){   [these terms are now included above]
            for(m=0;m<21;m++)
              zphi[m] += L[m][k1Rna+nat2+k1S+4];
            if(props[R].rt != 6)
              for(m=0;m<21;m++)
                zphi[m] += L[m][k1Rna+nat3+k1S+4];
          } */

          for(m=0;m<21;m++) allvars[m][3][vv] = zphi[m];
          vv++;

          /* if current vars are free, compute V_ij and add it to Aop_ij */
          if(i>=0 && tfv[i].iko==R && tfv[i].code==30 &&
             j>=0 && tfv[j].iko==S && tfv[j].code==30){
            calcV(i, j, isi, vs, tfv, zphi, Aop);
            /* for(m=0;m<21;m++) RR[m][isize+j] = zphi[m]; */
            j--;
          }


          if(props[S].epoc != 1) continue;    /* pass to next S */

          for(m=0;m<21;m++)                                  /* sph_S(ich2) */
            zsph[m] = zphi[m] + allvars[m][4][vv] - allvars[m][4][vv-1];
          calcL(k1R+2, k1S+1, 1, dc2, pdb_model, damp_scale, zsph);
          if(props[R].rt != 6){
            calcL(k1R+3, k1S+1, 1, dc2, pdb_model, damp_scale, zsph);
            if(props[R].rt != 1)
              calcL(k1R+4, k1S+1, 1, dc2, pdb_model, damp_scale, zsph);
          }
          /* for(m=0;m<21;m++)
              zsph[m] = zphi[m] + allvars[m][4][vv] - allvars[m][4][vv-1]
                       + L[m][k1Rna+nat2+k1S+1];
          if(props[R].rt != 6){
            for(m=0;m<21;m++)
              zsph[m] += L[m][k1Rna+nat3+k1S+1];
            if(props[R].rt != 1)
              for(m=0;m<21;m++)
                zsph[m] += L[m][k1Rna+nat4+k1S+1];
          } */

          for(m=0;m<21;m++) allvars[m][3][vv] = zsph[m];
          vv++;

          /* if current vars are free, compute V_ij and add it to Aop_ij */
          if(i>=0 && tfv[i].iko==R && tfv[i].code==30 &&
             j>=0 && tfv[j].iko==S && tfv[j].code/10==2){
            calcV(i, j,   isi, vs, tfv, zsph, Aop);
            calcV(i, j-1, isi, vs, tfv, zsph, Aop);
            /* for(m=0;m<21;m++){
                 RR[m][isize+j]   = zsph[m];
                 RR[m][isize+j-1] = zsph[m];
               } */
            j -= 2;
          }

          for(m=0;m<21;m++)                                /* car_S(ich2) */
            zcar[m] = zsph[m] + allvars[m][4][vv] - allvars[m][4][vv-1];
          calcL(k1R+2, k1S, 1, dc2, pdb_model, damp_scale, zcar);
          if(props[R].rt != 6){
            calcL(k1R+3, k1S, 1, dc2, pdb_model, damp_scale, zcar);
            if(props[R].rt != 1)
              calcL(k1R+4, k1S, 1, dc2, pdb_model, damp_scale, zcar);
          }
          /* for(m=0;m<21;m++)
              zcar[m] = zsph[m] + allvars[m][4][vv] - allvars[m][4][vv-1]
                       + L[m][k1Rna+nat2+k1S];
          if(props[R].rt != 6){
            for(m=0;m<21;m++)
              zcar[m] += L[m][k1Rna+nat3+k1S];
            if(props[R].rt != 1)
              for(m=0;m<21;m++)
                zcar[m] += L[m][k1Rna+nat4+k1S];
          } */

          for(m=0;m<21;m++){
            allvars[m][3][vv] = zcar[m];
             carall[m][ch2ncnv+ch1nv+ivar+2] = zcar[m];
          }
          vv++;

          /* if current vars are free, compute V_ij and add it to Aop_ij */
          if(i>=0 && tfv[i].iko==R && tfv[i].code==30 &&
             j>=0 && tfv[j].iko==S && tfv[j].code/10==1){
            calcV(i, j,   isi, vs, tfv, zcar, Aop);
            calcV(i, j-1, isi, vs, tfv, zcar, Aop);
            calcV(i, j-2, isi, vs, tfv, zcar, Aop);
            /* for(m=0;m<21;m++){ 
                 RR[m][isize+j]   = zcar[m];
                 RR[m][isize+j-1] = zcar[m];
                 RR[m][isize+j-2] = zcar[m];
               } */
            j -= 3;
          }

        }     /* end of loop over S */

        if(i>=0 && tfv[i].iko==R && tfv[i].code==30){
          i--;
          isize -= size;
          isi -= sizex;
        }

        if(props[R].epoc != 1) continue;    /* pass to next R */

    /* sph_R (ich1) ******************************************************/

        vv=0;
        j=chinfo[ich2].ivl;

        for(S=chinfo[ich2].ikl; S>=chinfo[ich2].ikf; S--){

          k1S = props[S].k1;          /* index of 1st pseudoatom of residue S */
          if(props[S].epoc != 2)
            k1Sp = props[S+1].k1;     /* index of 1st pseudoatom of residue S+1 */

          if(props[S].rt == 1 || props[S].rt == 6)       /* chi_S (ich2) */
            for(m=0;m<21;m++)
              zchi[m] = 0.0;
          else{
            for(m=0;m<21;m++) zchi[m] = allvars[m][3][vv];
            calcL(k1R+1, k1S+4, 1, dc2, pdb_model, damp_scale, zchi);
            /* for(m=0;m<21;m++) zchi[m] = allvars[m][3][vv] + L[m][k1Rna+num_atoms+k1S+4]; */
          }

          for(m=0;m<21;m++) allvars[m][2][vv] = zchi[m];
          vv++;

          /* if current vars are free, compute V_ij and add it to Aop_ij */
          if(i>=0 && tfv[i].iko==R && tfv[i].code/10==2 &&
             j>=0 && tfv[j].iko==S && tfv[j].code==50){
            calcV(i,   j, isi,       vs, tfv, zchi, Aop);
            calcV(i-1, j, isi-sizex, vs, tfv, zchi, Aop);
            /* for(m=0;m<21;m++){
                 RR[m][isize+j]      = zchi[m];
                 RR[m][isize-size+j] = zchi[m];
               } */
            j--;
          }


          if(props[S].epoc == 2)                           /* psi_S (ich2) */
            for(m=0;m<21;m++)
              zpsi[m] = 0.0;
          else{
            for(m=0;m<21;m++)
              zpsi[m] = zphi[m] + allvars[m][3][vv] - allvars[m][3][vv-2];
            calcL(k1R+1, k1Sp,   1, dc2, pdb_model, damp_scale, zpsi);
            calcL(k1R+1, k1Sp+1, 1, dc2, pdb_model, damp_scale, zpsi);
            /* for(m=0;m<21;m++)
                zpsi[m] = zphi[m] + allvars[m][3][vv] - allvars[m][3][vv-2]
                         + L[m][k1Rna+num_atoms+k1Sp] + L[m][k1Rna+num_atoms+k1Sp+1]; */
          }

          for(m=0;m<21;m++) allvars[m][2][vv] = zpsi[m];
          vv++;

          /* if current vars are free, compute V_ij and add it to Aop_ij */
          if(i>=0 && tfv[i].iko==R && tfv[i].code/10==2 &&
             j>=0 && tfv[j].iko==S && tfv[j].code==40){
            calcV(i,   j, isi,       vs, tfv, zpsi, Aop);
            calcV(i-1, j, isi-sizex, vs, tfv, zpsi, Aop);
            /* for(m=0;m<21;m++){
                 RR[m][isize+j]      = zpsi[m];
                 RR[m][isize-size+j] = zpsi[m];
               } */
            j--;
          }


          for(m=0;m<21;m++)                                  /* phi_S (ich2) */
            zphi[m]  =   zpsi[m] + zchi[m] + allvars[m][3][vv]
                       - allvars[m][3][vv-1] - allvars[m][3][vv-2];
          calcL(k1R+1, k1S+2, 1, dc2, pdb_model, damp_scale, zphi);
          if(props[S].rt != 6)
            calcL(k1R+1, k1S+3, 1, dc2, pdb_model, damp_scale, zphi);
          /* for(m=0;m<21;m++)
              zphi[m]  =   zpsi[m] + zchi[m] + allvars[m][3][vv]
                       - allvars[m][3][vv-1] - allvars[m][3][vv-2]
                       + L[m][k1Rna+num_atoms+k1S+2];
          if(props[S].rt != 6){
            for(m=0;m<21;m++)
              zphi[m] += L[m][k1Rna+num_atoms+k1S+3];
            if(props[S].rt != 1)     [this term is now included above]
              for(m=0;m<21;m++)
                zphi[m] += L[m][k1Rna+num_atoms+k1S+4];
          } */

          for(m=0;m<21;m++) allvars[m][2][vv] = zphi[m];
          vv++;

          /* if current vars are free, compute V_ij and add it to Aop_ij */
          if(i>=0 && tfv[i].iko==R && tfv[i].code/10==2 &&
             j>=0 && tfv[j].iko==S && tfv[j].code==30){
            calcV(i,   j, isi,       vs, tfv, zphi, Aop);
            calcV(i-1, j, isi-sizex, vs, tfv, zphi, Aop);
            /* for(m=0;m<21;m++){
                 RR[m][isize+j]      = zphi[m];
                 RR[m][isize-size+j] = zphi[m];
               } */
            j--;
          }


          if(props[S].epoc != 1) continue;    /* pass to next S */

          for(m=0;m<21;m++)                                  /* sph_S(ich2) */
            zsph[m] = zphi[m] + allvars[m][3][vv] - allvars[m][3][vv-1];
          calcL(k1R+1, k1S+1, 1, dc2, pdb_model, damp_scale, zsph);
          /* for(m=0;m<21;m++)
              zsph[m] = zphi[m] + allvars[m][3][vv] - allvars[m][3][vv-1]
                       + L[m][k1Rna+num_atoms+k1S+1]; */

          for(m=0;m<21;m++) allvars[m][2][vv] = zsph[m];
          vv++;

          /* if current vars are free, compute V_ij and add it to Aop_ij */
          if(i>=0 && tfv[i].iko==R && tfv[i].code/10==2 &&
             j>=0 && tfv[j].iko==S && tfv[j].code/10==2){
            calcV(i,   j,   isi,       vs, tfv, zsph, Aop);
            calcV(i,   j-1, isi,       vs, tfv, zsph, Aop);
            calcV(i-1, j,   isi-sizex, vs, tfv, zsph, Aop);
            calcV(i-1, j-1, isi-sizex, vs, tfv, zsph, Aop);
            /* for(m=0;m<21;m++){
                 RR[m][isize+j]        = zsph[m];
                 RR[m][isize+j-1]      = zsph[m];
                 RR[m][isize-size+j]   = zsph[m];
                 RR[m][isize-size+j-1] = zsph[m];
               } */
            j -= 2;
          }

          for(m=0;m<21;m++)                                /* car_S(ich2) */
            zcar[m] = zsph[m] + allvars[m][3][vv] - allvars[m][3][vv-1];
          calcL(k1R+1, k1S, 1, dc2, pdb_model, damp_scale, zcar);
          /* for(m=0;m<21;m++)
              zcar[m] = zsph[m] + allvars[m][3][vv] - allvars[m][3][vv-1]
                       + L[m][k1Rna+num_atoms+k1S]; */

          for(m=0;m<21;m++){
            allvars[m][2][vv] = zcar[m];
             carall[m][ch2ncnv+ch1nv+1] = zcar[m];
          }
          vv++;

          /* if current vars are free, compute V_ij and add it to Aop_ij */
          if(i>=0 && tfv[i].iko==R && tfv[i].code/10==2 &&
             j>=0 && tfv[j].iko==S && tfv[j].code/10==1){
            calcV(i,   j,   isi,       vs, tfv, zcar, Aop);
            calcV(i,   j-1, isi,       vs, tfv, zcar, Aop);
            calcV(i,   j-2, isi,       vs, tfv, zcar, Aop);
            calcV(i-1, j,   isi-sizex, vs, tfv, zcar, Aop);
            calcV(i-1, j-1, isi-sizex, vs, tfv, zcar, Aop);
            calcV(i-1, j-2, isi-sizex, vs, tfv, zcar, Aop);
            /* for(m=0;m<21;m++){
                 RR[m][isize+j]        = zcar[m];
                 RR[m][isize+j-1]      = zcar[m];
                 RR[m][isize+j-2]      = zcar[m];
                 RR[m][isize-size+j]   = zcar[m];
                 RR[m][isize-size+j-1] = zcar[m];
                 RR[m][isize-size+j-2] = zcar[m];
               } */
            j -= 3;
          }

        }     /* end of loop over S */

        if(i>=0 && tfv[i].iko==R && tfv[i].code/10==2){
          i -= 2;
          isize -= siz2;
          isi -= sizx2;
        }


    /* car_R (ich1) ******************************************************/

        vv=0;
        j=chinfo[ich2].ivl;

        for(S=chinfo[ich2].ikl; S>=chinfo[ich2].ikf; S--){

          ivar = S-chinfo[ich2].ikf;
          ivar = (ivar<<1) + ivar;

          k1S = props[S].k1;          /* index of 1st pseudoatom of residue S */
          if(props[S].epoc != 2)
            k1Sp = props[S+1].k1;     /* index of 1st pseudoatom of residue S+1 */

          if(props[S].rt == 1 || props[S].rt == 6)       /* chi_S (ich2) */
            for(m=0;m<21;m++)
              zchi[m] = 0.0;
          else{
            for(m=0;m<21;m++) zchi[m] = allvars[m][2][vv];
            calcL(k1R, k1S+4, 1, dc2, pdb_model, damp_scale, zchi);
            /*for(m=0;m<21;m++)
                zchi[m] = allvars[m][2][vv] + L[m][k1Rna+k1S+4]; */
          }

          for(m=0;m<21;m++){
            allvars[m][1][vv] = zchi[m];
             carall[m][ch1ncnv+ch2nv+ivar+4] = zchi[m];
          }
          vv++;

          /* if current vars are free, compute V_ij and add it to Aop_ij */
          if(i>=0 && tfv[i].iko==R && tfv[i].code/10==1 &&
             j>=0 && tfv[j].iko==S && tfv[j].code==50){
            calcV(i,   j, isi,       vs, tfv, zchi, Aop);
            calcV(i-1, j, isi-sizex, vs, tfv, zchi, Aop);
            calcV(i-2, j, isi-sizx2, vs, tfv, zchi, Aop);
            /* for(m=0;m<21;m++){
                 RR[m][isize+j]      = zchi[m];
                 RR[m][isize-size+j] = zchi[m];
                 RR[m][isize-siz2+j] = zchi[m];
               } */
            j--;
          }


          if(props[S].epoc == 2)                           /* psi_S (ich2) */
            for(m=0;m<21;m++)
              zpsi[m] = 0.0;
          else{
            for(m=0;m<21;m++)
              zpsi[m] = zphi[m] + allvars[m][2][vv] - allvars[m][2][vv-2];
            calcL(k1R, k1Sp,   1, dc2, pdb_model, damp_scale, zpsi);
            calcL(k1R, k1Sp+1, 1, dc2, pdb_model, damp_scale, zpsi);
            /* for(m=0;m<21;m++)
                zpsi[m] = zphi[m] + allvars[m][2][vv] - allvars[m][2][vv-2]
                         + L[m][k1Rna+k1Sp] + L[m][k1Rna+k1Sp+1]; */
          }

          for(m=0;m<21;m++){
            allvars[m][1][vv] = zpsi[m];
             carall[m][ch1ncnv+ch2nv+ivar+3] = zpsi[m];
          }
          vv++;

          /* if current vars are free, compute V_ij and add it to Aop_ij */
          if(i>=0 && tfv[i].iko==R && tfv[i].code/10==1 &&
             j>=0 && tfv[j].iko==S && tfv[j].code==40){
            calcV(i,   j, isi,       vs, tfv, zpsi, Aop);
            calcV(i-1, j, isi-sizex, vs, tfv, zpsi, Aop);
            calcV(i-2, j, isi-sizx2, vs, tfv, zpsi, Aop);
            /* for(m=0;m<21;m++){
                 RR[m][isize+j]      = zpsi[m];
                 RR[m][isize-size+j] = zpsi[m];
                 RR[m][isize-siz2+j] = zpsi[m];
               } */
            j--;
          }


          for(m=0;m<21;m++)                                  /* phi_S (ich2) */
            zphi[m]  =   zpsi[m] + zchi[m]   + allvars[m][2][vv] 
                       - allvars[m][2][vv-1] - allvars[m][2][vv-2];
          calcL(k1R, k1S+2, 1, dc2, pdb_model, damp_scale, zphi);
          if(props[S].rt != 6)
            calcL(k1R, k1S+3, 1, dc2, pdb_model, damp_scale, zphi);
          /* for(m=0;m<21;m++)
              zphi[m]  =   zpsi[m] + zchi[m] + allvars[m][2][vv] 
                         - allvars[m][2][vv-1] - allvars[m][2][vv-2]
                         + L[m][k1Rna+k1S+2];
          if(props[S].rt != 6){
            for(m=0;m<21;m++)
              zphi[m] += L[m][k1Rna+k1S+3];
            if(props[S].rt != 1)   [this term is now included above]
              for(m=0;m<21;m++)
                zphi[m] += L[m][k1Rna+k1S+4];
          } */

          for(m=0;m<21;m++){
            allvars[m][1][vv] = zphi[m];
             carall[m][ch1ncnv+ch2nv+ivar+2] = zphi[m];
          }
          vv++;

          /* if current vars are free, compute V_ij and add it to Aop_ij */
          if(i>=0 && tfv[i].iko==R && tfv[i].code/10==1 &&
             j>=0 && tfv[j].iko==S && tfv[j].code==30){
            calcV(i,   j, isi,       vs, tfv, zphi, Aop);
            calcV(i-1, j, isi-sizex, vs, tfv, zphi, Aop);
            calcV(i-2, j, isi-sizx2, vs, tfv, zphi, Aop);
            /* for(m=0;m<21;m++){
                 RR[m][isize+j]      = zphi[m];
                 RR[m][isize-size+j] = zphi[m];
                 RR[m][isize-siz2+j] = zphi[m];
               } */
            j--;
          }


          if(props[S].epoc != 1) continue;    /* pass to next S */

          for(m=0;m<21;m++)                                  /* sph_S(ich2) */
            zsph[m] = zphi[m] + allvars[m][2][vv] - allvars[m][2][vv-1];
          calcL(k1R, k1S+1, 1, dc2, pdb_model, damp_scale, zsph);
          /* for(m=0;m<21;m++)
              zsph[m] = zphi[m] + allvars[m][2][vv] - allvars[m][2][vv-1]
                       + L[m][k1Rna+k1S+1]; */

          for(m=0;m<21;m++){
            allvars[m][1][vv] = zsph[m];
             carall[m][ch1ncnv+ch2nv+1] = zsph[m];
          }
          vv++;

          /* if current vars are free, compute V_ij and add it to Aop_ij */
          if(i>=0 && tfv[i].iko==R && tfv[i].code/10==1 &&
             j>=0 && tfv[j].iko==S && tfv[j].code/10==2){
            calcV(i,   j,   isi,       vs, tfv, zsph, Aop);
            calcV(i,   j-1, isi,       vs, tfv, zsph, Aop);
            calcV(i-1, j,   isi-sizex, vs, tfv, zsph, Aop);
            calcV(i-1, j-1, isi-sizex, vs, tfv, zsph, Aop);
            calcV(i-2, j,   isi-sizx2, vs, tfv, zsph, Aop);
            calcV(i-2, j-1, isi-sizx2, vs, tfv, zsph, Aop);
            /* for(m=0;m<21;m++){
                 RR[m][isize+j]        = zsph[m];
                 RR[m][isize+j-1]      = zsph[m];
                 RR[m][isize-size+j]   = zsph[m];
                 RR[m][isize-size+j-1] = zsph[m];
                 RR[m][isize-siz2+j]   = zsph[m];
                 RR[m][isize-siz2+j-1] = zsph[m];
               } */
            j -= 2;
          }


          for(m=0;m<21;m++)                                /* car_S(ich2) */
            zcar[m] = zsph[m] + allvars[m][2][vv] - allvars[m][2][vv-1];
          calcL(k1R, k1S, 1, dc2, pdb_model, damp_scale, zcar);
          
          /* for(m=0;m<21;m++)
              zcar[m] = zsph[m] + allvars[m][2][vv] - allvars[m][2][vv-1]
                       + L[m][k1Rna+k1S]; */

          for(m=0;m<21;m++){
            allvars[m][1][vv] = zcar[m];
             carall[m][ch1ncnv+ch2nv] = zcar[m];
             carall[m][ch2ncnv+ch1nv] = zcar[m];
          }
          vv++;

          /* if current vars are free, compute V_ij and add it to Aop_ij */
          if(i>=0 && tfv[i].iko==R && tfv[i].code/10==1 &&
             j>=0 && tfv[j].iko==S && tfv[j].code/10==1){
            calcV(i,   j,   isi,       vs, tfv, zcar, Aop);
            calcV(i,   j-1, isi,       vs, tfv, zcar, Aop);
            calcV(i,   j-2, isi,       vs, tfv, zcar, Aop);
            calcV(i-1, j,   isi-sizex, vs, tfv, zcar, Aop);
            calcV(i-1, j-1, isi-sizex, vs, tfv, zcar, Aop);
            calcV(i-1, j-2, isi-sizex, vs, tfv, zcar, Aop);
            calcV(i-2, j,   isi-sizx2, vs, tfv, zcar, Aop);
            calcV(i-2, j-1, isi-sizx2, vs, tfv, zcar, Aop);
            calcV(i-2, j-2, isi-sizx2, vs, tfv, zcar, Aop);
            /* for(m=0;m<21;m++){
                 RR[m][isize+j]        = zcar[m];
                 RR[m][isize+j-1]      = zcar[m];
                 RR[m][isize+j-2]      = zcar[m];
                 RR[m][isize-size+j]   = zcar[m];
                 RR[m][isize-size+j-1] = zcar[m];
                 RR[m][isize-size+j-2] = zcar[m];
                 RR[m][isize-siz2+j]   = zcar[m];
                 RR[m][isize-siz2+j-1] = zcar[m];
                 RR[m][isize-siz2+j-2] = zcar[m];
               } */
            j -= 3;
          }

        }     /* end of loop over S */

        if(i>=0 && tfv[i].iko==R && tfv[i].code/10==1){
          i -= 3;
          isize -= siz3;
          isi -= sizx3;
        }

      }      /* end of loop over R    */
    }        /* end of loop over ich2 */
  }          /* end of loop over ich1 */



  /* SECONDLY, scan pairs of vars in the SAME chain */

  for(ich1=0;ich1<nchains;ich1++){

    i = chinfo[ich1].ivf;
    isize = i*size;
    isi = i*sizex;

    ch1nv = ich1*mnv;

    for(R=chinfo[ich1].ikf; R<=chinfo[ich1].ikl; R++){

      k1R = props[R].k1;          /* index of 1st pseudoatom of residue R */
      if(props[R].epoc != 2)
        k1Rp = props[R+1].k1;     /* index of 1st pseudoatom of residue R+1 */

      k1Rna  = k1R  * num_atoms;
      k1Rpna = k1Rp * num_atoms;

      if(props[R].epoc != 1) goto phiR;


  /* car_R (ich1) ******************************************************/

      vv=0;
      j=chinfo[ich1].ivl;

      for(S=chinfo[ich1].ikl; S>=R; S--){

        ivar = S-chinfo[ich1].ikf;
        ivar = (ivar<<1) + ivar;

        k1S = props[S].k1;          /* index of 1st pseudoatom of residue S */
        if(props[S].epoc != 2)
          k1Sp = props[S+1].k1;     /* index of 1st pseudoatom of residue S+1 */

        for(m=0;m<21;m++)                                 /* chi_S (ich1) */
          zchi[m] = 0.0;
        for(ich2=0,index=ch1nv+ivar+4; ich2<nchains; ich2++,index+=ncnv){
          if(ich2==ich1) continue;
          for(m=0;m<21;m++)
            zchi[m] += carall[m][index];
        }

        for(m=0;m<21;m++) allvars[m][1][vv] = zchi[m];
        vv++;

        /* if current vars are free, compute V_ij and add it to Aop_ij */
        if(i>=0 && tfv[i].iko==R && tfv[i].code/10==1 &&
           j>=0 && tfv[j].iko==S && tfv[j].code==50){
          calcV(i,   j, isi,       vs, tfv, zchi, Aop);
          calcV(i+1, j, isi+sizex, vs, tfv, zchi, Aop);
          calcV(i+2, j, isi+sizx2, vs, tfv, zchi, Aop);
          /* for(m=0;m<21;m++){
               RR[m][isize+j]      = zchi[m];
               RR[m][isize+size+j] = zchi[m];
               RR[m][isize+siz2+j] = zchi[m];
             } */
          j--;
        }


        for(m=0;m<21;m++)                                 /* psi_S (ich1) */
          zpsi[m] = 0.0;
        for(ich2=0,index=ch1nv+ivar+3; ich2<nchains; ich2++,index+=ncnv){
          if(ich2==ich1) continue;
          for(m=0;m<21;m++)
            zpsi[m] += carall[m][index];
        }

        for(m=0;m<21;m++) allvars[m][1][vv] = zpsi[m];
        vv++;

        /* if current vars are free, compute V_ij and add it to Aop_ij */
        if(i>=0 && tfv[i].iko==R && tfv[i].code/10==1 &&
           j>=0 && tfv[j].iko==S && tfv[j].code==40){
          calcV(i,   j, isi,       vs, tfv, zpsi, Aop);
          calcV(i+1, j, isi+sizex, vs, tfv, zpsi, Aop);
          calcV(i+2, j, isi+sizx2, vs, tfv, zpsi, Aop);
          /* for(m=0;m<21;m++){
               RR[m][isize+j]      = zpsi[m];
               RR[m][isize+size+j] = zpsi[m];
               RR[m][isize+siz2+j] = zpsi[m];
             } */
          j--;
        }


        for(m=0;m<21;m++)                                  /* phi_S (ich1) */
          zphi[m] = 0.0;
        for(ich2=0,index=ch1nv+ivar+2; ich2<nchains; ich2++,index+=ncnv){
          if(ich2==ich1) continue;
          for(m=0;m<21;m++)
            zphi[m] += carall[m][index];
        }

        for(m=0;m<21;m++) allvars[m][1][vv] = zphi[m];
        vv++;

        /* if current vars are free, compute V_ij and add it to Aop_ij */
        if(i>=0 && tfv[i].iko==R && tfv[i].code/10==1 &&
           j>=0 && tfv[j].iko==S && tfv[j].code==30){
          calcV(i,   j, isi,       vs, tfv, zphi, Aop);
          calcV(i+1, j, isi+sizex, vs, tfv, zphi, Aop);
          calcV(i+2, j, isi+sizx2, vs, tfv, zphi, Aop);
          /* for(m=0;m<21;m++){
               RR[m][isize+j]      = zphi[m];
               RR[m][isize+size+j] = zphi[m];
               RR[m][isize+siz2+j] = zphi[m];
             } */
          j--;
        }


        if(props[S].epoc != 1) continue;    /* pass to next S */

        for(m=0;m<21;m++)                                  /* sph_S(ich1) */
          zsph[m] = 0.0;
        for(ich2=0,index=ch1nv+1; ich2<nchains; ich2++,index+=ncnv){
          if(ich2==ich1) continue;
          for(m=0;m<21;m++)
            zsph[m] += carall[m][index];
        }

        for(m=0;m<21;m++) allvars[m][1][vv] = zsph[m];
        vv++;

        /* if current vars are free, compute V_ij and add it to Aop_ij */
        if(i>=0 && tfv[i].iko==R && tfv[i].code/10==1 &&
           j>=0 && tfv[j].iko==S && tfv[j].code/10==2){
          calcV(i,   j,   isi,       vs, tfv, zsph, Aop);
          calcV(i,   j-1, isi,       vs, tfv, zsph, Aop);
          calcV(i+1, j,   isi+sizex, vs, tfv, zsph, Aop);
          calcV(i+1, j-1, isi+sizex, vs, tfv, zsph, Aop);
          calcV(i+2, j,   isi+sizx2, vs, tfv, zsph, Aop);
          calcV(i+2, j-1, isi+sizx2, vs, tfv, zsph, Aop);
          /* for(m=0;m<21;m++){
               RR[m][isize+j]        = zsph[m];
               RR[m][isize+j-1]      = zsph[m];
               RR[m][isize+size+j]   = zsph[m];
               RR[m][isize+size+j-1] = zsph[m];
               RR[m][isize+siz2+j]   = zsph[m];
               RR[m][isize+siz2+j-1] = zsph[m];
             } */
          j -= 2;
        }


        for(m=0;m<21;m++)                                /* car_S(ich1) */
          zcar[m] = 0.0;
        for(ich2=0,index=ch1nv; ich2<nchains; ich2++,index+=ncnv){
          if(ich2==ich1) continue;
          for(m=0;m<21;m++)
            zcar[m] += carall[m][index];
        }

        for(m=0;m<21;m++) allvars[m][1][vv] = zcar[m];
        vv++;

        /* if current vars are free, compute V_ij and add it to Aop_ij (only for "i<=j") */
        if(i>=0 && tfv[i].iko==R && tfv[i].code/10==1 &&
           j>=0 && tfv[j].iko==S && tfv[j].code/10==1){
          calcV(i,   j,   isi,       vs, tfv, zcar, Aop);
          calcV(i,   j-1, isi,       vs, tfv, zcar, Aop);
          calcV(i,   j-2, isi,       vs, tfv, zcar, Aop);
          calcV(i+1, j,   isi+sizex, vs, tfv, zcar, Aop);
          calcV(i+1, j-1, isi+sizex, vs, tfv, zcar, Aop);
          calcV(i+2, j,   isi+sizx2, vs, tfv, zcar, Aop);
          /* for(m=0;m<21;m++){
               RR[m][isize+j]        = zcar[m];
               RR[m][isize+j-1]      = zcar[m];
               RR[m][isize+j-2]      = zcar[m];
               RR[m][isize+size+j]   = zcar[m];
               RR[m][isize+size+j-1] = zcar[m];
               RR[m][isize+siz2+j]   = zcar[m];
             } */
          j -= 3;
        }

      }      /* end of loop over S    */

      if(i>=0 && tfv[i].iko==R && tfv[i].code/10==1){
        i += 3;
        isize += siz3;
        isi += sizx3;
      }


/* sph_R (ich1) ******************************************************/

        vv=0;
        j=chinfo[ich1].ivl;

        for(S=chinfo[ich1].ikl; S>=R; S--){

          k1S = props[S].k1;          /* index of 1st pseudoatom of residue S */
          if(props[S].epoc != 2)
            k1Sp = props[S+1].k1;     /* index of 1st pseudoatom of residue S+1 */


          if(props[S].rt == 1 || props[S].rt == 6)           /* chi_S (ich1) */
            for(m=0;m<21;m++)
              zchi[m] = 0.0;
          else{
            for(m=0;m<21;m++) zchi[m] = allvars[m][1][vv];
            calcL(k1R, k1S+4, 1, dc2, pdb_model, damp_scale, zchi);
            /* for(m=0;m<21;m++)
                zchi[m] = allvars[m][1][vv] + L[m][k1Rna+k1S+4]; */
          }

          for(m=0;m<21;m++) allvars[m][2][vv] = zchi[m];
          vv++;

          /* if current vars are free, compute V_ij and add it to Aop_ij */
          if(i>=0 && tfv[i].iko==R && tfv[i].code/10==2 &&
             j>=0 && tfv[j].iko==S && tfv[j].code==50){
            calcV(i,   j, isi,       vs, tfv, zchi, Aop);
            calcV(i+1, j, isi+sizex, vs, tfv, zchi, Aop);
            /* for(m=0;m<21;m++){
                 RR[m][isize+j]      = zchi[m];
                 RR[m][isize+size+j] = zchi[m];
               } */
            j--;
          }


          if(props[S].epoc == 2)                           /* psi_S (ich1) */
            for(m=0;m<21;m++)
              zpsi[m] = 0.0;
          else{
            for(m=0;m<21;m++)
              zpsi[m] = zphi[m] + allvars[m][1][vv] - allvars[m][1][vv-2];
            calcL(k1R, k1Sp,   1, dc2, pdb_model, damp_scale, zpsi);
            calcL(k1R, k1Sp+1, 1, dc2, pdb_model, damp_scale, zpsi);
            /* for(m=0;m<21;m++)
                zpsi[m] = zphi[m] + allvars[m][1][vv] - allvars[m][1][vv-2]
                         + L[m][k1Rna+k1Sp] + L[m][k1Rna+k1Sp+1]; */
          }

          for(m=0;m<21;m++) allvars[m][2][vv] = zpsi[m];
          vv++;

          /* if current vars are free, compute V_ij and add it to Aop_ij */
          if(i>=0 && tfv[i].iko==R && tfv[i].code/10==2 &&
             j>=0 && tfv[j].iko==S && tfv[j].code==40){
            calcV(i,   j, isi,       vs, tfv, zpsi, Aop);
            calcV(i+1, j, isi+sizex, vs, tfv, zpsi, Aop);
            /* for(m=0;m<21;m++){
                 RR[m][isize+j]      = zpsi[m];
                 RR[m][isize+size+j] = zpsi[m];
               } */
            j--;
          }


          for(m=0;m<21;m++)                                  /* phi_S (ich1) */
            zphi[m] =  zpsi[m] + zchi[m]    + allvars[m][1][vv]
                      - allvars[m][1][vv-1] - allvars[m][1][vv-2];
          calcL(k1R, k1S+2, 1, dc2, pdb_model, damp_scale, zphi);
          if(props[S].rt != 6)
            calcL(k1R, k1S+3, 1, dc2, pdb_model, damp_scale, zphi);
          /* for(m=0;m<21;m++)
              zphi[m] =  zpsi[m] + zchi[m] + allvars[m][1][vv] - allvars[m][1][vv-1]
                        - allvars[m][1][vv-2] + L[m][k1Rna+k1S+2];
          if(props[S].rt != 6)
            for(m=0;m<21;m++)
              zphi[m] += L[m][k1Rna+k1S+3]; */

          for(m=0;m<21;m++) allvars[m][2][vv] = zphi[m];
          vv++;

          /* if current vars are free, compute V_ij and add it to Aop_ij */
          if(i>=0 && tfv[i].iko==R && tfv[i].code/10==2 &&
             j>=0 && tfv[j].iko==S && tfv[j].code==30){
            calcV(i,   j, isi,       vs, tfv, zphi, Aop);
            calcV(i+1, j, isi+sizex, vs, tfv, zphi, Aop);
            /* for(m=0;m<21;m++){
                 RR[m][isize+j]      = zphi[m];
                 RR[m][isize+size+j] = zphi[m];
               } */
            j--;
          }


          if(props[S].epoc != 1) continue;    /* pass to next S */

          for(m=0;m<21;m++)                                  /* sph_S(ich1) */
            zsph[m] = zphi[m] + allvars[m][1][vv] - allvars[m][1][vv-1];
          calcL(k1R, k1S+1, 1, dc2, pdb_model, damp_scale, zsph);
          /*  for(m=0;m<21;m++)
              zsph[m] = zphi[m] + allvars[m][1][vv] - allvars[m][1][vv-1]
                       + L[m][k1Rna+k1S+1]; */

          for(m=0;m<21;m++) allvars[m][2][vv] = zsph[m];
          vv++;

          /* if current vars are free, compute V_ij and add it to Aop_ij (only for "i<=j") */
          if(i>=0 && tfv[i].iko==R && tfv[i].code/10==2 &&
             j>=0 && tfv[j].iko==S && tfv[j].code/10==2){
            calcV(i,   j,   isi,       vs, tfv, zsph, Aop);
            calcV(i,   j-1, isi,       vs, tfv, zsph, Aop);
            calcV(i+1, j,   isi+sizex, vs, tfv, zsph, Aop);
            /* for(m=0;m<21;m++){
                 RR[m][isize+j]        = zsph[m];
                 RR[m][isize+j-1]      = zsph[m];
                 RR[m][isize+size+j]   = zsph[m];
               } */
            j -= 2;
          }

        }     /* end of loop over S */

        if(i>=0 && tfv[i].iko==R && tfv[i].code/10==2){
          i += 2;
          isize += siz2;
          isi += sizx2;
        }


/* phi_R (ich1) ******************************************************/

phiR:   vv=0;
        j=chinfo[ich1].ivl;

        for(S=chinfo[ich1].ikl; S>=R; S--){

          k1S = props[S].k1;          /* index of 1st pseudoatom of residue S */
          if(props[S].epoc != 2)
            k1Sp = props[S+1].k1;     /* index of 1st pseudoatom of residue S+1 */

          if(props[S].rt == 1 || props[S].rt == 6)           /* chi_S (ich1) */
            for(m=0;m<21;m++)
              zchi[m] = 0.0;
          else if(props[R].epoc == 1){
            for(m=0;m<21;m++) zchi[m] = allvars[m][2][vv];
            calcL(k1R+1, k1S+4, 1, dc2, pdb_model, damp_scale, zchi);
            /* for(m=0;m<21;m++)
                zchi[m] = allvars[m][2][vv] + L[m][k1Rna+num_atoms+k1S+4]; */
          }
          else{
            for(m=0;m<21;m++) zchi[m] = allvars[m][4][vv];
            calcL(k1R,   k1S+4, 1, dc2, pdb_model, damp_scale, zchi);
            calcL(k1R+1, k1S+4, 1, dc2, pdb_model, damp_scale, zchi);
            /* for(m=0;m<21;m++)
                 zchi[m] = allvars[m][4][vv] + L[m][k1Rna+k1S+4]
                                             + L[m][k1Rna+num_atoms+k1S+4]; */
          }

          for(m=0;m<21;m++) allvars[m][3][vv] = zchi[m];
          vv++;

          /* if current vars are free, compute V_ij and add it to Aop_ij */
          if(i>=0 && tfv[i].iko==R && tfv[i].code==30 &&
             j>=0 && tfv[j].iko==S && tfv[j].code==50){
            calcV(i, j, isi, vs, tfv, zchi, Aop);
            /* for(m=0;m<21;m++) RR[m][isize+j] = zchi[m]; */
            j--;
          }


          if(props[S].epoc == 2)                           /* psi_S (ich1) */
            for(m=0;m<21;m++)
              zpsi[m] = 0.0;
          else if(props[R].epoc == 1){
            for(m=0;m<21;m++)
              zpsi[m] = zphi[m] + allvars[m][2][vv] - allvars[m][2][vv-2];
            calcL(k1R+1, k1Sp,   1, dc2, pdb_model, damp_scale, zpsi);
            calcL(k1R+1, k1Sp+1, 1, dc2, pdb_model, damp_scale, zpsi);
            /* for(m=0;m<21;m++)
                zpsi[m] = zphi[m] + allvars[m][2][vv] - allvars[m][2][vv-2]
                         + L[m][k1Rna+num_atoms+k1Sp] + L[m][k1Rna+num_atoms+k1Sp+1]; */
          }
          else{
            for(m=0;m<21;m++)
              zpsi[m] = zphi[m] + allvars[m][4][vv] - allvars[m][4][vv-2];
            calcL(k1R,   k1Sp,   1, dc2, pdb_model, damp_scale, zpsi);
            calcL(k1R,   k1Sp+1, 1, dc2, pdb_model, damp_scale, zpsi);
            calcL(k1R+1, k1Sp,   1, dc2, pdb_model, damp_scale, zpsi);
            calcL(k1R+1, k1Sp+1, 1, dc2, pdb_model, damp_scale, zpsi);
            /* for(m=0;m<21;m++)
                zpsi[m] = zphi[m] + allvars[m][4][vv] - allvars[m][4][vv-2]
                         + L[m][k1Rna+k1Sp] + L[m][k1Rna+k1Sp+1]
                         + L[m][k1Rna+num_atoms+k1Sp] + L[m][k1Rna+num_atoms+k1Sp+1]; */
          }

          for(m=0;m<21;m++) allvars[m][3][vv] = zpsi[m];
          vv++;

          /* if current vars are free, compute V_ij and add it to Aop_ij */
          if(i>=0 && tfv[i].iko==R && tfv[i].code==30 &&
             j>=0 && tfv[j].iko==S && tfv[j].code==40){
            calcV(i, j, isi, vs, tfv, zpsi, Aop);
            /* for(m=0;m<21;m++) RR[m][isize+j] = zpsi[m]; */
            j--;
          }


          if(props[R].epoc == 1){                              /* phi_S (ich1) */
            for(m=0;m<21;m++)
              zphi[m] =  zpsi[m] + zchi[m]    + allvars[m][2][vv]
                        - allvars[m][2][vv-1] - allvars[m][2][vv-2];
            calcL(k1R+1, k1S+2, 1, dc2, pdb_model, damp_scale, zphi);
            if(props[S].rt != 6)
              calcL(k1R+1, k1S+3, 1, dc2, pdb_model, damp_scale, zphi);
            
            /* for(m=0;m<21;m++)
                zphi[m] =  zpsi[m] + zchi[m]    + allvars[m][2][vv]
                         - allvars[m][2][vv-1] - allvars[m][2][vv-2]
                         + L[m][k1Rna+num_atoms+k1S+2];
            if(props[S].rt != 6)
              for(m=0;m<21;m++)
                zphi[m] += L[m][k1Rna+num_atoms+k1S+3]; */
          }
          else{
            for(m=0;m<21;m++)
              zphi[m] =  zpsi[m] + zchi[m]    + allvars[m][4][vv]
                        - allvars[m][4][vv-1] - allvars[m][4][vv-2];
            calcL(k1R,   k1S+2, 1, dc2, pdb_model, damp_scale, zphi);
            calcL(k1R+1, k1S+2, 1, dc2, pdb_model, damp_scale, zphi);
            if(props[S].rt != 6){
              calcL(k1R,   k1S+3, 1, dc2, pdb_model, damp_scale, zphi);
              calcL(k1R+1, k1S+3, 1, dc2, pdb_model, damp_scale, zphi);
            }
            /* for(m=0;m<21;m++)
                zphi[m]  =   zpsi[m] + zchi[m] + allvars[m][4][vv]
                          - allvars[m][4][vv-1] - allvars[m][4][vv-2]
                          + L[m][k1Rna+k1S+2] + L[m][k1Rna+num_atoms+k1S+2];
            if(props[S].rt != 6)
              for(m=0;m<21;m++)
                zphi[m] += L[m][k1Rna+k1S+3] + L[m][k1Rna+num_atoms+k1S+3]; */
          }

          for(m=0;m<21;m++) allvars[m][3][vv] = zphi[m];
          vv++;

          /* if current vars are free, compute V_ij and add it to Aop_ij */
          if(i>=0 && tfv[i].iko==R && tfv[i].code==30 &&
             j>=0 && tfv[j].iko==S && tfv[j].code==30){
            calcV(i, j, isi, vs, tfv, zphi, Aop);
            /* for(m=0;m<21;m++) RR[m][isize+j] = zphi[m]; */
            j--;
          }

        }     /* end of loop over S */

        if(i>=0 && tfv[i].iko==R && tfv[i].code==30){
          i++;
          isize += size;
          isi += sizex;
        }


/* psi_R (ich1) ******************************************************/

        vv=0;
        j=chinfo[ich1].ivl;

        for(S=chinfo[ich1].ikl; S>=R; S--){

          k1S = props[S].k1;          /* index of 1st pseudoatom of residue S */
          if(props[S].epoc != 2)
            k1Sp = props[S+1].k1;     /* index of 1st pseudoatom of residue S+1 */

          if(S == R){       /* this case (for chi_S) is done in the chi_R block, below */
            vv++;
            if(i>=0 && tfv[i].iko==R && tfv[i].code==40 &&
               j>=0 && tfv[j].iko==S && tfv[j].code==50)
              j--;
            goto psiS;
          }

          if(props[R].epoc == 2 ||
             props[S].rt == 1 || props[S].rt == 6)       /* chi_S (ich1) */
            for(m=0;m<21;m++)
              zchi[m] = 0.0;
          else{
            for(m=0;m<21;m++) zchi[m] = allvars[m][3][vv];
            calcL(k1R+2, k1S+4, 1, dc2, pdb_model, damp_scale, zchi);
            if(props[R].rt != 6){
              calcL(k1R+3, k1S+4, 1, dc2, pdb_model, damp_scale, zchi);
              if(props[R].rt != 1)
                calcL(k1R+4, k1S+4, 1, dc2, pdb_model, damp_scale, zchi);
            }
            /* for(m=0;m<21;m++)
                zchi[m] = allvars[m][3][vv] + L[m][k1Rna+nat2+k1S+4];
            if(props[R].rt != 6){
              for(m=0;m<21;m++)
                zchi[m] += L[m][k1Rna+nat3+k1S+4];
              if(props[R].rt != 1)
                for(m=0;m<21;m++)
                  zchi[m] += L[m][k1Rna+nat4+k1S+4];
            } */
          }

          for(m=0;m<21;m++) allvars[m][4][vv] = zchi[m];
          vv++;

          /* if current vars are free, compute V_ij and add it to Aop_ij */
          if(i>=0 && tfv[i].iko==R && tfv[i].code==40 &&
             j>=0 && tfv[j].iko==S && tfv[j].code==50){
            calcV(i, j, isi, vs, tfv, zchi, Aop);
            /* for(m=0;m<21;m++) RR[m][isize+j] = zchi[m]; */
            j--;
          }


psiS:     if(props[R].epoc == 2 || props[S].epoc == 2)     /* psi_S (ich1) */
            for(m=0;m<21;m++)
              zpsi[m] = 0.0;
          else{
            for(m=0;m<21;m++)
              zpsi[m] = zphi[m] + allvars[m][3][vv] - allvars[m][3][vv-2];
            calcL(k1R+2, k1Sp,   1, dc2, pdb_model, damp_scale, zpsi);
            calcL(k1R+2, k1Sp+1, 1, dc2, pdb_model, damp_scale, zpsi);
            if(props[R].rt != 6){
              calcL(k1R+3, k1Sp,   1, dc2, pdb_model, damp_scale, zpsi);
              calcL(k1R+3, k1Sp+1, 1, dc2, pdb_model, damp_scale, zpsi);
              if(props[R].rt != 1)
                calcL(k1R+4, k1Sp,   1, dc2, pdb_model, damp_scale, zpsi);
                calcL(k1R+4, k1Sp+1, 1, dc2, pdb_model, damp_scale, zpsi);
            }
            /* for(m=0;m<21;m++)
                zpsi[m] = zphi[m] + allvars[m][3][vv] - allvars[m][3][vv-2]
                                 + L[m][k1Rna+nat2+k1Sp] + L[m][k1Rna+nat2+k1Sp+1];
            if(props[R].rt != 6){
              for(m=0;m<21;m++)
                zpsi[m] += L[m][k1Rna+nat3+k1Sp] + L[m][k1Rna+nat3+k1Sp+1];
              if(props[R].rt != 1)
                for(m=0;m<21;m++)
                  zpsi[m] += L[m][k1Rna+nat4+k1Sp] + L[m][k1Rna+nat4+k1Sp+1];
            } */
          }

          for(m=0;m<21;m++) allvars[m][4][vv] = zpsi[m];
          vv++;

          /* if current vars are free, compute V_ij and add it to Aop_ij */
          if(i>=0 && tfv[i].iko==R && tfv[i].code==40 &&
             j>=0 && tfv[j].iko==S && tfv[j].code==40){
            calcV(i, j, isi, vs, tfv, zpsi, Aop);
            /* for(m=0;m<21;m++) RR[m][isize+j] = zpsi[m]; */
            j--;
          }

          if(S == R) continue;       /* don't do phi_S when S=R */

          if(props[R].epoc == 2)                            /* phi_S (ich1) */
            for(m=0;m<21;m++)
              zphi[m] = 0.0;
          else{
            for(m=0;m<21;m++)
              zphi[m] =  zpsi[m] + zchi[m]    + allvars[m][3][vv]
                        - allvars[m][3][vv-1] - allvars[m][3][vv-2];
            calcL(k1R+2, k1S+2, 1, dc2, pdb_model, damp_scale, zphi);
            if(props[S].rt != 6)
              calcL(k1R+2, k1S+3, 1, dc2, pdb_model, damp_scale, zphi);
            if(props[R].rt != 6)
              calcL(k1R+3, k1S+2, 1, dc2, pdb_model, damp_scale, zphi);
            if(props[R].rt != 6 && props[S].rt != 6)
              calcL(k1R+3, k1S+3, 1, dc2, pdb_model, damp_scale, zphi);
            if(props[R].rt != 1 && props[R].rt != 6){
              calcL(k1R+4, k1S+2, 1, dc2, pdb_model, damp_scale, zphi);
              if(props[S].rt != 6)
                calcL(k1R+4, k1S+3, 1, dc2, pdb_model, damp_scale, zphi);
            }
            /* for(m=0;m<21;m++)
                zphi[m]  =   zpsi[m]  +  zchi[m] + allvars[m][3][vv]
                          - allvars[m][3][vv-1] - allvars[m][3][vv-2]
                          + L[m][k1Rna+nat2+k1S+2];
            if(props[S].rt != 6)
              for(m=0;m<21;m++)
                zphi[m] += L[m][k1Rna+nat2+k1S+3];
            if(props[R].rt != 6)
              for(m=0;m<21;m++)
                zphi[m] += L[m][k1Rna+nat3+k1S+2];
            if(props[R].rt != 6 && props[S].rt != 6)
              for(m=0;m<21;m++)
                zphi[m] += L[m][k1Rna+nat3+k1S+3];
            if(props[R].rt != 1 && props[R].rt != 6){
              for(m=0;m<21;m++)
                zphi[m] += L[m][k1Rna+nat4+k1S+2];
              if(props[S].rt != 6)
                for(m=0;m<21;m++)
                  zphi[m] += L[m][k1Rna+nat4+k1S+3];
            } */
          }

          for(m=0;m<21;m++) allvars[m][4][vv] = zphi[m];
          vv++;

          /* if current vars are free, compute V_ij and add it to Aop_ij */
          if(i>=0 && tfv[i].iko==R && tfv[i].code==40 &&
             j>=0 && tfv[j].iko==S && tfv[j].code==30){
            calcV(i, j, isi, vs, tfv, zphi, Aop);
            /* for(m=0;m<21;m++) RR[m][isize+j] = zphi[m]; */
            j--;
          }

        }    /* end of loop over S */

        if(i>=0 && tfv[i].iko==R && tfv[i].code==40){
          i++;
          isize += size;
          isi += sizex;
        }


/* chi_R (ich1) ******************************************************/

        vv=0;
        j=chinfo[ich1].ivl;

        for(S=chinfo[ich1].ikl; S>=R; S--){

          k1S = props[S].k1;          /* index of 1st pseudoatom of residue S */
          if(props[S].epoc != 2)
            k1Sp = props[S+1].k1;     /* index of 1st pseudoatom of residue S+1 */

          if(props[R].rt == 1 || props[R].rt == 6 ||
             props[S].rt == 1 || props[S].rt == 6)         /* chi_S (ich1) */
            for(m=0;m<21;m++)
              zchi[m] = 0.0;
          else if(S == R){
            for(m=0;m<21;m++) zchi[m] = allvars[m][3][vv];
            calcL(k1R+2, k1S+4, 1, dc2, pdb_model, damp_scale, zchi);
            calcL(k1R+3, k1S+4, 1, dc2, pdb_model, damp_scale, zchi);
            if(props[R].epoc != 2){
              for(m=0;m<21;m++) zchi[m] += zphi[m];
              calcL(k1Rp,   k1S+4, 1, dc2, pdb_model, damp_scale, zchi);
              calcL(k1Rp+1, k1S+4, 1, dc2, pdb_model, damp_scale, zchi);
            }                                                       /* (from Case 14) */
            /* for(m=0;m<21;m++)
                zchi[m] = allvars[m][3][vv] + L[m][k1Rna+nat2+k1S+4]
                                            + L[m][k1Rna+nat3+k1S+4];
            if(props[R].epoc != 2)
              for(m=0;m<21;m++)
                zchi[m] += zphi[m] + L[m][k1Rpna+k1S+4]
                                   + L[m][k1Rpna+num_atoms+k1S+4];   (from Case 14) */
          }
          else{
            calcL(k1R+4, k1S+4, 0, dc2, pdb_model, damp_scale, zchi);
            /* for(m=0;m<21;m++)
                zchi[m] = L[m][k1Rna+nat4+k1S+4]; */
          }

          for(m=0;m<21;m++) allvars[m][5][vv] = zchi[m];
          vv++;

          /* if current vars are free, compute V_ij and add it to Aop_ij */
          if(i>=0 && tfv[i].iko==R && tfv[i].code==50 &&
             j>=0 && tfv[j].iko==S && tfv[j].code==50){
            calcV(i, j, isi, vs, tfv, zchi, Aop);
            /* for(m=0;m<21;m++) RR[m][isize+j] = zchi[m]; */
            j--;
          }


          if(props[R].rt == 1 || props[R].rt == 6 ||
             props[S].epoc == 2)                         /* psi_S (ich1) */
            for(m=0;m<21;m++)
              zpsi[m] = 0.0;
          else{
            for(m=0;m<21;m++) zpsi[m] = zphi[m];
            calcL(k1R+4, k1Sp,   1, dc2, pdb_model, damp_scale, zpsi);
            calcL(k1R+4, k1Sp+1, 1, dc2, pdb_model, damp_scale, zpsi);
            /* for(m=0;m<21;m++)
                zpsi[m] = zphi[m] + L[m][k1Rna+nat4+k1Sp]
                                  + L[m][k1Rna+nat4+k1Sp+1]; */
          }

          for(m=0;m<21;m++) allvars[m][5][vv] = zpsi[m];

          /* if S=R, store the missing case from the psi_R block above
             (this is not used now, but is stored anyway for completeness,
             and in case you find a use for it in the future :-) */
          if(S == R)
            for(m=0;m<21;m++) allvars[m][4][vv-1] = zpsi[m];

          vv++;

          /* if current vars are free, copy zpsi to R_ij */
          if(i>=0 && tfv[i].iko==R && tfv[i].code==50 &&
             j>=0 && tfv[j].iko==S && tfv[j].code==40){
            if(S == R){             /* j<i in this case (missing case from above) */
              /* index = j*size+i; */
              calcV(j, i, j*sizex, vs, tfv, zpsi, Aop);
            }
            else{
              /* index = isize+j; */
              calcV(i, j, isi, vs, tfv, zpsi, Aop);
            }
            /* for(m=0;m<21;m++) RR[m][index] = zpsi[m]; */
            j--;
          }

          if(S == R) continue;       /* don't do phi_S when S=R */

          if(props[R].rt == 1 || props[R].rt == 6)        /* phi_S (ich1) */
            for(m=0;m<21;m++)
              zphi[m] = 0.0;
          else{
            for(m=0;m<21;m++) zphi[m] = zpsi[m] + zchi[m];
            calcL(k1R+4, k1S+2, 1, dc2, pdb_model, damp_scale, zphi);
            if(props[S].rt != 6)
              calcL(k1R+4, k1S+3, 1, dc2, pdb_model, damp_scale, zphi);
            /* for(m=0;m<21;m++)
                zphi[m] = zpsi[m] + zchi[m] + L[m][k1Rna+nat4+k1S+2];
            if(props[S].rt != 6)
              for(m=0;m<21;m++)
                zphi[m] += L[m][k1Rna+nat4+k1S+3]; */
          }

          for(m=0;m<21;m++) allvars[m][5][vv] = zphi[m];
          vv++;

          /* if current vars are free, compute V_ij and add it to Aop_ij */
          if(i>=0 && tfv[i].iko==R && tfv[i].code==50 &&
             j>=0 && tfv[j].iko==S && tfv[j].code==30){
            calcV(i, j, isi, vs, tfv, zphi, Aop);
            /* for(m=0;m<21;m++) RR[m][isize+j] = zphi[m]; */
            j--;
          }

        }      /* end of loop over S */

        if(i>=0 && tfv[i].iko==R && tfv[i].code==50){
          i++;
          isize += size;
          isi += sizex;
        }

    }        /* end of loop over R    */
  }          /* end of loop over ich1 */


/* FINISHED WITH Rij and Vij */



/* add entries corresponding to constraints *********************/

  /* first, constraints from static residues */

  for(p=0;p<srlen;p++){         /* p = index in the 'statres' array */
    i = statres[p].iko;
    k = props[i].k1+1;          /* index of r_b */
    temp  = pow(pdb_model[k-1].x-pdb_model[k].x,2);
    temp += pow(pdb_model[k-1].y-pdb_model[k].y,2);
    temp += pow(pdb_model[k-1].z-pdb_model[k].z,2);
    temp=sqrt(temp);                     /* norm of r_a - r_b  */
    v[0] = (pdb_model[k-1].x-pdb_model[k].x)/temp;
    v[1] = (pdb_model[k-1].y-pdb_model[k].y)/temp;
    v[2] = (pdb_model[k-1].z-pdb_model[k].z)/temp;

    temp  = pow(pdb_model[k+1].x-pdb_model[k].x,2);
    temp += pow(pdb_model[k+1].y-pdb_model[k].y,2);
    temp += pow(pdb_model[k+1].z-pdb_model[k].z,2);
    temp=sqrt(temp);                     /* norm of r_c - r_b  */
    w[0] = (pdb_model[k+1].x-pdb_model[k].x)/temp;
    w[1] = (pdb_model[k+1].y-pdb_model[k].y)/temp;
    w[2] = (pdb_model[k+1].z-pdb_model[k].z)/temp;

    cosg = v[0]*w[0]+v[1]*w[1]+v[2]*w[2];     /* cos(ang(v,w)) */

    for(j=0;j<3;j++)                   /* w - cosg*v */
      wcv[j] = w[j]-cosg*v[j];

    vw[0] = v[1]*w[2]-v[2]*w[1];
    vw[1] = v[2]*w[0]-v[0]*w[2];
    vw[2] = v[0]*w[1]-v[1]*w[0];       /*  v x w  */

    for(i=0,isi=0; i<size; i++,isi+=sizex){
      Wilson(k, i, num_atoms, pdb_model, props, erb, tfv, der);
      Aop[isi+size+6*p+0] = der[k].x;
      Aop[isi+size+6*p+1] = der[k].y;
      Aop[isi+size+6*p+2] = der[k].z;
      Wilson(k-1, i, num_atoms, pdb_model, props, erb, tfv, der);
      sum  = der[k-1].x * wcv[0];
      sum += der[k-1].y * wcv[1];
      sum += der[k-1].z * wcv[2];
      Aop[isi+size+6*p+3] = sum;
      sum  = der[k-1].x * vw[0];
      sum += der[k-1].y * vw[1];
      sum += der[k-1].z * vw[2];
      Aop[isi+size+6*p+4] = sum;
      Wilson(k+1, i, num_atoms, pdb_model, props, erb, tfv, der);
      sum  = der[k+1].x * vw[0];
      sum += der[k+1].y * vw[1];
      sum += der[k+1].z * vw[2];
      Aop[isi+size+6*p+5] = sum;
    }
  }

  /* second, distant constraints     (added on 10/27/2016)      */
  for (c=0; c<dclen; c++){           /* c = index in the 'dconst' array */
    i1 = dconst[c].iko;
    i2 = dconst[c].ikt;
    k  = props[i1].k1+1;    /* index of C-alpha atom of  first residue of the pair */
    l  = props[i2].k1+1;    /* index of C-alpha atom of second residue of the pair */

    for(i=0,isi=0; i<size; i++,isi+=sizex){
      v[0] = pdb_model[k].x - pdb_model[l].x;
      v[1] = pdb_model[k].y - pdb_model[l].y;
      v[2] = pdb_model[k].z - pdb_model[l].z;
      Wilson(k, i, num_atoms, pdb_model, props, erb, tfv, der);
      sum = der[k].x * v[0] + der[k].y * v[1] + der[k].z * v[2];
      Wilson(l, i, num_atoms, pdb_model, props, erb, tfv, der);
      sum -= der[l].x * v[0] + der[l].y * v[1] + der[l].z * v[2];
      Aop[isi+size+nco1+c] = sum;
    }
  }
/***************************************************************/

  /* fill in lower triangular part of Aop */
  for(i=0,isi=0; i<sizex; i++,isi+=sizex)
    for(j=0,jsi=0; j<i; j++,jsi+=sizex)
      Aop[isi+j] = Aop[jsi+i];


/*============================
for(i=0,isi=0; i<sizex; i++,isi+=sizex){
  for(j=0; j<=i; j++)
    printf("%.4e ", Aop[isi+j]);
  printf("\n");
}
============================*/

/* test of constraints
  fca = fopen("test_distance_constraints_matrix.txt", "w");
  for(i=0,isi=0; i<size; i++,isi+=sizex){
    for(j=size+nco1; j<sizex; j++)
      fprintf(fca, "%.4e ", Aop[isi+j]);
    fprintf(fca, "\n");
  }
  fclose(fca);
*/

  /* compute Q^(m)_j, j=1,...,size, and extend it with zeros: */
  for(j=0;j<size;j++){
    qm[j] = 0.0;
    Wilson(-1, j, num_atoms, pdb_model, props, erb, tfv, der);
    for(k=0; k<num_atoms; k++){
      qm[j] += wop[k].x * der[k].x;
      qm[j] += wop[k].y * der[k].y;
      qm[j] += wop[k].z * der[k].z;
    }
  }
  for(j=size;j<sizex;j++)
    qm[j]=0.0;


/* test:
Fnorm = 0.0; Gnorm = 0.0;
for(k=0; k<num_atoms; k++)
  Fnorm += pow(wop[k].x,2)+pow(wop[k].y,2)+pow(wop[k].z,2);
for(j=0;j<size;j++)
  Gnorm += pow(qm[j],2);
Fnorm = sqrt(Fnorm); Gnorm = sqrt(Gnorm);
fprintf(stderr,"Norm of cart. forces = %e   Norm of gen. forces = %e\n",Fnorm,Gnorm);
*/


  /**********************************************/
  /*                                            */
  /*   now solve the linear system:             */
  /*                                            */
  /*             Aop * qdots  =  qm             */
  /*                                            */
  /**********************************************/

  /* In general, to call a Fortran routine
     from C we have to transform the matrix
     from row-major order to column-major order,
     but since in our case the matrix is symmetric,
     this is not necessary.
  */

  /*  OLD: calling dsysv_
  uplo = 'U';
  nrhs=1;
  lda=sizex;
  ldb=sizex;
  lwork = 5*sizex;

  work = (double *) malloc(lwork * sizeof(double));
  ptr_check(work,25);

  dsysv_(&uplo, &sizex, &nrhs, Aop, &lda, ipiv,
         qm, &ldb, work, &lwork, &info);
  */

  int *perm;        // permutation done in the LU decomposition
  double parity;    // parity of the permutation
  perm = (int *) malloc(sizex*sizeof(int));

  /* calling separate functions for LU and back sub (earlier approach):
  1- perform the LU decomp of the matrix Aop:
  info = ludcmp(Aop, sizex, perm, &parity);
  if (info != 0)     // can't do the backsubstitution if matrix is singular
    return info;
  2- if matrix is invertible, do the backsubstitution:
  lubksb(Aop, sizex, perm, qm);  */

  /* calling the single function linsolver: */
  info = linsolver(Aop, sizex, perm, &parity, qm);

  if (info != 0)     // can't do the backsubstitution if matrix is singular
    return info;

  /* copy solution to qdots: */
  for(j=0;j<size;j++)
    qdots[j] = qm[j];

  free(Aop); free(qm);
  /* free(work);
  free(ipiv); */

  for(m=0;m<6;m++)
    free(vs[m]);

  for(m=0;m<10;m++){
    free(Y[m]); free(Z[m]);
  }

  for(m=0;m<21;m++){
    /* free(L[m]); free(RR[m]); */
    free(carall[m]);
  }

  for(ivar=1;ivar<=5;ivar++)
    for(m=0;m<21;m++)
      free(allvars[m][ivar]);

  return info;

}


/*=====================================================================================*/
int calcL(int k, int l, int ep, double dc2, PDB *pdb_model, double damp_scale, double Lp[21])
{
  /* for a pair of pseudo-atoms (k,l) define the damping constant and compute L_{kl} */
  /* dc2 = dist_cut^2 */
  
  /* if ep = 0,  store L[m] in Lp[m] */
  /* if ep != 0,  add  L[m] to Lp[m] */

  double v[3], w[3], r[3], vw[3];
  double temp, temp1, prod1;
  int m, mu, nu;
  int skip, special;

  special = 0;   /* signals a special feature to be considered in the "if" below */

  if (special){
    /* for the particular case of thermosome, the special feature is
       to skip dampers across chains when x^2+y^2<25^2 or |z|>50: */
    skip=0;
    if(pdb_model[k].chain[0] != pdb_model[l].chain[0]){
      if(pow(pdb_model[k].x,2)+pow(pdb_model[k].y,2)<625.0 || fabs(pdb_model[k].z)>50.0)
        skip=1;
      else{
        if(pow(pdb_model[l].x,2)+pow(pdb_model[l].y,2)<625.0 || fabs(pdb_model[l].z)>50.0)
          skip=1;
      }
    }
    if(skip==1){
      if(ep==0)
        for(m=0;m<21;m++)
          Lp[m]=0.0;
      return 0;
    }
    /* end of thermosome case */
  }

  temp = pow(pdb_model[k].x-pdb_model[l].x, 2) +
         pow(pdb_model[k].y-pdb_model[l].y, 2) +
         pow(pdb_model[k].z-pdb_model[l].z, 2);

  if(temp>dc2){
    if(ep==0)
      for(m=0;m<21;m++)
        Lp[m]=0.0;
    return 0;
  }

  temp1 = temp;               /* distance^2 between atoms k and l */
  if(temp<0.25) temp=0.25;
  
  /* D_{kl}: */
  prod1 = damp_scale * 1.22/pow(temp,0.25) / temp1;

  v[0] = pdb_model[k].x; w[0] = pdb_model[l].x;
  v[1] = pdb_model[k].y; w[1] = pdb_model[l].y;
  v[2] = pdb_model[k].z; w[2] = pdb_model[l].z;
  r[0] = v[0]-w[0]; r[1] = v[1]-w[1]; r[2] = v[2]-w[2];
  vw[0] = v[1]*w[2]-v[2]*w[1];
  vw[1] = v[2]*w[0]-v[0]*w[2];
  vw[2] = v[0]*w[1]-v[1]*w[0];

  m=0;
  if(ep==0){
    for(mu=0;mu<3;mu++)
      for(nu=mu;nu<3;nu++){
        Lp[m] = prod1 * vw[mu] * vw[nu];
        m++;
      }
    for(mu=0;mu<3;mu++)
      for(nu=0;nu<3;nu++){
        Lp[m] = prod1 * vw[mu] * r[nu];
        m++;
      }
    for(mu=0;mu<3;mu++)
      for(nu=mu;nu<3;nu++){
        Lp[m] = prod1 * r[mu] * r[nu];
        m++;
      }
  }
  else{
    for(mu=0;mu<3;mu++)
      for(nu=mu;nu<3;nu++){
        Lp[m] += prod1 * vw[mu] * vw[nu];
        m++;
      }
    for(mu=0;mu<3;mu++)
      for(nu=0;nu<3;nu++){
        Lp[m] += prod1 * vw[mu] * r[nu];
        m++;
      }
    for(mu=0;mu<3;mu++)
      for(nu=mu;nu<3;nu++){
        Lp[m] += prod1 * r[mu] * r[nu];
        m++;
      }
  }

  return 1;

}


/*=====================================================================================*/
void calcV(int i, int j, long isi, double *vs[6], fvars *tfv, double RR[21], double *Aop)
{
  /* for a pair of 'variable' indices (i,j), compute Vij and add it to Aop */
  /* isi = i*sizex */

  double sum, vsi[6], vsj[6];
  int m, sigma, xi, eta;

  /* copy the 6-dim factors (e_i, e_i x r_{beta(i)}): */
  for(m=0;m<6;m++){
    vsi[m] = vs[m][i];
    vsj[m] = vs[m][j];
  }


  if( tfv[j].cho  != tfv[i].cho ||
     (tfv[j].code == 30 && tfv[i].code == 50 && tfv[j].iko >  tfv[i].iko) ||
     (tfv[i].code == 40 && tfv[j].code == 50 && tfv[i].iko == tfv[j].iko) ||
     (tfv[j].code == 40 && tfv[i].code == 50 && tfv[j].iko >  tfv[i].iko) ||
     (tfv[j].code == 50 && tfv[i].code == 50 && tfv[j].iko != tfv[i].iko) )
    sigma = -1;   /* M_i inter M_j = empty, hence sigma_{ij} = -1 */
  else
    sigma = 1;    /* M_j is a subset of M_i, hence sigma_{ij} = 1 */


  if(tfv[i].code >= 20 && tfv[j].code >= 20){    /* q_i and q_j angular */

    /* terms with mu=nu: */
    sum  = RR[0]  * vsi[0] * vsj[0] + RR[3]  * vsi[1] * vsj[1];
    sum += RR[5]  * vsi[2] * vsj[2] + RR[15] * vsi[3] * vsj[3];
    sum += RR[18] * vsi[4] * vsj[4] + RR[20] * vsi[5] * vsj[5];

    /* terms with mu<nu: */
    sum += RR[1]  * (vsi[0]*vsj[1] + vsi[1]*vsj[0]);
    sum += RR[2]  * (vsi[0]*vsj[2] + vsi[2]*vsj[0]);
    sum += RR[4]  * (vsi[1]*vsj[2] + vsi[2]*vsj[1]);
    sum += RR[6]  * (vsi[0]*vsj[3] + vsi[3]*vsj[0]);
    sum += RR[7]  * (vsi[0]*vsj[4] + vsi[4]*vsj[0]);
    sum += RR[8]  * (vsi[0]*vsj[5] + vsi[5]*vsj[0]);
    sum += RR[9]  * (vsi[1]*vsj[3] + vsi[3]*vsj[1]);
    sum += RR[10] * (vsi[1]*vsj[4] + vsi[4]*vsj[1]);
    sum += RR[11] * (vsi[1]*vsj[5] + vsi[5]*vsj[1]);
    sum += RR[12] * (vsi[2]*vsj[3] + vsi[3]*vsj[2]);
    sum += RR[13] * (vsi[2]*vsj[4] + vsi[4]*vsj[2]);
    sum += RR[14] * (vsi[2]*vsj[5] + vsi[5]*vsj[2]);
    sum += RR[16] * (vsi[3]*vsj[4] + vsi[4]*vsj[3]);
    sum += RR[17] * (vsi[3]*vsj[5] + vsi[5]*vsj[3]);
    sum += RR[19] * (vsi[4]*vsj[5] + vsi[5]*vsj[4]);

    /* Vij: */
    Aop[isi+j] += sigma * sum;
  }

  if(tfv[i].code >= 20 && tfv[j].code < 20){     /* q_i angular, q_j linear */
    xi = tfv[j].code - 10;   /* 0, 1, 2, for x, y, z, respectively */
    if(xi==0){
      sum  = RR[6]  * vsi[0] + RR[9]  * vsi[1];
      sum += RR[12] * vsi[2] + RR[15] * vsi[3];
      sum += RR[16] * vsi[4] + RR[17] * vsi[5];
    }
    if(xi==1){
      sum  = RR[7]  * vsi[0] + RR[10] * vsi[1];
      sum += RR[13] * vsi[2] + RR[16] * vsi[3];
      sum += RR[18] * vsi[4] + RR[19] * vsi[5];
    }
    if(xi==2){
      sum  = RR[8]  * vsi[0] + RR[11] * vsi[1];
      sum += RR[14] * vsi[2] + RR[17] * vsi[3];
      sum += RR[19] * vsi[4] + RR[20] * vsi[5];
    }
    Aop[isi+j] += -sigma * sum;
  }

  if(tfv[i].code < 20 && tfv[j].code >= 20){     /* q_i linear, q_j angular */
    xi = tfv[i].code - 10;   /* 0, 1, 2, for x, y, z, respectively */
    if(xi==0){
      sum  = RR[6]  * vsj[0] + RR[9]  * vsj[1];
      sum += RR[12] * vsj[2] + RR[15] * vsj[3];
      sum += RR[16] * vsj[4] + RR[17] * vsj[5];
    }
    if(xi==1){
      sum  = RR[7]  * vsj[0] + RR[10] * vsj[1];
      sum += RR[13] * vsj[2] + RR[16] * vsj[3];
      sum += RR[18] * vsj[4] + RR[19] * vsj[5];
    }
    if(xi==2){
      sum  = RR[8]  * vsj[0] + RR[11] * vsj[1];
      sum += RR[14] * vsj[2] + RR[17] * vsj[3];
      sum += RR[19] * vsj[4] + RR[20] * vsj[5];
    }
    Aop[isi+j] += -sigma * sum;
  }

  if(tfv[i].code < 20 && tfv[j].code < 20){      /* q_i and q_j linear */
    xi  = tfv[i].code - 10;   /* 0, 1, 2, for x, y, z, respectively */
    eta = tfv[j].code - 10;   /* 0, 1, 2, for x, y, z, respectively */
    if(xi==0){
      if(eta==0) sum = RR[15];
      if(eta==1) sum = RR[16];
      if(eta==2) sum = RR[17];
    }
    if(xi==1){
      if(eta==0) sum = RR[16];
      if(eta==1) sum = RR[18];
      if(eta==2) sum = RR[19];
    }
    if(xi==2){
      if(eta==0) sum = RR[17];
      if(eta==1) sum = RR[19];
      if(eta==2) sum = RR[20];
    }
    Aop[isi+j] += sigma * sum;
  }

}



/*=====================================================================================*/
int Wilson(int k0, int j, int num_atoms, PDB *pdb_model, tri *props,
           twovec *erb, fvars *tfv, trd *der)
{
  /* compute derivatives of cartesian coordinates of:
     - all pseudo-atoms  if  k0<0,  OR
     - pseudo-atom k0    if  k0>=0,
     with respect to the variable q_j */

  /* returns:  0 if successful;
               1 if tfv[j].code was not valid */


  int i, k, k1, k2, lepoc, nepoc, ilrc, klac;
  double r[3];

  if(k0<0)
    for(k=0; k<num_atoms; k++){
      der[k].x = 0.0;
      der[k].y = 0.0;
      der[k].z = 0.0;
    }
  else
    if(pdb_model[k0].chain[0] != tfv[j].cho){
      der[k0].x = 0.0; der[k0].y = 0.0; der[k0].z = 0.0;
      return 0;
    }


  i  = tfv[j].iko;          /* index of residue to which q_j belongs */
  k1 = props[i].k1;         /* index of 1st atom of residue i */
  lepoc = props[i].epoc;


  if(tfv[j].code == 10){               /* q_j = x of N_i=pdb_model[k1] */
    if(k0<0)
      for(k=k1;k<num_atoms;k++){
        nepoc = props[pdb_model[k].footnote].epoc;
        if(nepoc == 1 && lepoc == 2)     /* finished current chain */
          break;
        lepoc=nepoc;
        der[k].x = 1.0;
      }
    else{
      der[k0].x = 1.0; der[k0].y = 0.0; der[k0].z = 0.0;
    }
    return 0;
  }

  if(tfv[j].code == 11){               /* q_j = y of N_i=pdb_model[k1] */
    if(k0<0)
      for(k=k1;k<num_atoms;k++){
        nepoc = props[pdb_model[k].footnote].epoc;
        if(nepoc == 1 && lepoc == 2)     /* finished current chain */
          break;
        lepoc=nepoc;
        der[k].y = 1.0;
      }
    else{
      der[k0].x = 0.0; der[k0].y = 1.0; der[k0].z = 0.0;
    }
    return 0;
  }

  if(tfv[j].code == 12){               /* q_j = z of N_i=pdb_model[k1] */
    if(k0<0)
      for(k=k1;k<num_atoms;k++){
        nepoc = props[pdb_model[k].footnote].epoc;
        if(nepoc == 1 && lepoc == 2)     /* finished current chain */
          break;
        lepoc=nepoc;
        der[k].z = 1.0;
      }
    else{
      der[k0].x = 0.0; der[k0].y = 0.0; der[k0].z = 1.0;
    }
    return 0;
  }

  if(tfv[j].code == 20){               /* q_j = lambda */
    if(k0<0)
      for(k=k1+1;k<num_atoms;k++){
        nepoc=props[pdb_model[k].footnote].epoc;
        if(nepoc == 1 && lepoc == 2)       /* finished current chain */
          break;
        lepoc=nepoc;
        r[0]=pdb_model[k].x; r[1]=pdb_model[k].y; r[2]=pdb_model[k].z;
        der[k].x = erb[j].e.y * r[2] - erb[j].e.z * r[1];
        der[k].y = erb[j].e.z * r[0] - erb[j].e.x * r[2];
        der[k].z = erb[j].e.x * r[1] - erb[j].e.y * r[0];
      }
    else{
      if(k0>k1){
        r[0]=pdb_model[k0].x; r[1]=pdb_model[k0].y; r[2]=pdb_model[k0].z;
        der[k0].x = erb[j].e.y * r[2] - erb[j].e.z * r[1];
        der[k0].y = erb[j].e.z * r[0] - erb[j].e.x * r[2];
        der[k0].z = erb[j].e.x * r[1] - erb[j].e.y * r[0];
      }
      else{
        der[k0].x = 0.0; der[k0].y = 0.0; der[k0].z = 0.0;
      }
    }
    return 0;
  }

  if(tfv[j].code == 21){               /* q_j = theta */
    if(k0<0)
      for(k=k1+1;k<num_atoms;k++){
        nepoc=props[pdb_model[k].footnote].epoc;
        if(nepoc == 1 && lepoc == 2)     /* finished current chain */
          break;
        lepoc=nepoc;
        r[0]=pdb_model[k].x; r[1]=pdb_model[k].y; r[2]=pdb_model[k].z;
        der[k].x = erb[j].e.y * (r[2]-erb[j].b.z) - erb[j].e.z * (r[1]-erb[j].b.y);
        der[k].y = erb[j].e.z * (r[0]-erb[j].b.x) - erb[j].e.x * (r[2]-erb[j].b.z);
        der[k].z = erb[j].e.x * (r[1]-erb[j].b.y) - erb[j].e.y * (r[0]-erb[j].b.x);
      }
    else{
      if(k0>k1){
        r[0]=pdb_model[k0].x; r[1]=pdb_model[k0].y; r[2]=pdb_model[k0].z;
        der[k0].x = erb[j].e.y * (r[2]-erb[j].b.z) - erb[j].e.z * (r[1]-erb[j].b.y);
        der[k0].y = erb[j].e.z * (r[0]-erb[j].b.x) - erb[j].e.x * (r[2]-erb[j].b.z);
        der[k0].z = erb[j].e.x * (r[1]-erb[j].b.y) - erb[j].e.y * (r[0]-erb[j].b.x);
      }
      else{
        der[k0].x = 0.0; der[k0].y = 0.0; der[k0].z = 0.0;
      }
    }
    return 0;
  }

  if(tfv[j].code == 30){               /* q_j = phi_i */
    if(k0<0)
      for(k=k1+2;k<num_atoms;k++){
        nepoc=props[pdb_model[k].footnote].epoc;
        if(nepoc == 1 && lepoc == 2)     /* finished current chain */
          break;
        lepoc=nepoc;
        r[0]=pdb_model[k].x; r[1]=pdb_model[k].y; r[2]=pdb_model[k].z;
        der[k].x = erb[j].e.y * (r[2]-erb[j].b.z) - erb[j].e.z * (r[1]-erb[j].b.y);
        der[k].y = erb[j].e.z * (r[0]-erb[j].b.x) - erb[j].e.x * (r[2]-erb[j].b.z);
        der[k].z = erb[j].e.x * (r[1]-erb[j].b.y) - erb[j].e.y * (r[0]-erb[j].b.x);
      }
    else{
      if(k0>=k1+2){
        r[0]=pdb_model[k0].x; r[1]=pdb_model[k0].y; r[2]=pdb_model[k0].z;
        der[k0].x = erb[j].e.y * (r[2]-erb[j].b.z) - erb[j].e.z * (r[1]-erb[j].b.y);
        der[k0].y = erb[j].e.z * (r[0]-erb[j].b.x) - erb[j].e.x * (r[2]-erb[j].b.z);
        der[k0].z = erb[j].e.x * (r[1]-erb[j].b.y) - erb[j].e.y * (r[0]-erb[j].b.x);
      }
      else{
        der[k0].x = 0.0; der[k0].y = 0.0; der[k0].z = 0.0;
      }
    }
    return 0;
  }

  if(tfv[j].code == 40){               /* q_j  = psi_i */
    k2 = k1+props[i].nat;              /* index of 1st atom of residue i+1 */
    if(k0<0)
      for(k=k2;k<num_atoms;k++){
        nepoc=props[pdb_model[k].footnote].epoc;
        if(nepoc == 1 && lepoc == 2)     /* finished current chain */
          break;
        lepoc=nepoc;
        r[0]=pdb_model[k].x; r[1]=pdb_model[k].y; r[2]=pdb_model[k].z;
        der[k].x = erb[j].e.y * (r[2]-erb[j].b.z) - erb[j].e.z * (r[1]-erb[j].b.y);
        der[k].y = erb[j].e.z * (r[0]-erb[j].b.x) - erb[j].e.x * (r[2]-erb[j].b.z);
        der[k].z = erb[j].e.x * (r[1]-erb[j].b.y) - erb[j].e.y * (r[0]-erb[j].b.x);
      }
    else{
      if(k0>=k2){
        r[0]=pdb_model[k0].x; r[1]=pdb_model[k0].y; r[2]=pdb_model[k0].z;
        der[k0].x = erb[j].e.y * (r[2]-erb[j].b.z) - erb[j].e.z * (r[1]-erb[j].b.y);
        der[k0].y = erb[j].e.z * (r[0]-erb[j].b.x) - erb[j].e.x * (r[2]-erb[j].b.z);
        der[k0].z = erb[j].e.x * (r[1]-erb[j].b.y) - erb[j].e.y * (r[0]-erb[j].b.x);
      }
      else{
        der[k0].x = 0.0; der[k0].y = 0.0; der[k0].z = 0.0;
      }
    }
    return 0;
  }

  if(tfv[j].code == 50){               /* q_j = chi_i */
    if(k0<0){
      k=k1+4;  /* only pseudo-atom to be moved by chi */
      r[0]=pdb_model[k].x; r[1]=pdb_model[k].y; r[2]=pdb_model[k].z;
      der[k].x = erb[j].e.y * (r[2]-erb[j].b.z) - erb[j].e.z * (r[1]-erb[j].b.y);
      der[k].y = erb[j].e.z * (r[0]-erb[j].b.x) - erb[j].e.x * (r[2]-erb[j].b.z);
      der[k].z = erb[j].e.x * (r[1]-erb[j].b.y) - erb[j].e.y * (r[0]-erb[j].b.x);
    }
    else{
      if(k0==k1+4){
        r[0]=pdb_model[k0].x; r[1]=pdb_model[k0].y; r[2]=pdb_model[k0].z;
        der[k0].x = erb[j].e.y * (r[2]-erb[j].b.z) - erb[j].e.z * (r[1]-erb[j].b.y);
        der[k0].y = erb[j].e.z * (r[0]-erb[j].b.x) - erb[j].e.x * (r[2]-erb[j].b.z);
        der[k0].z = erb[j].e.x * (r[1]-erb[j].b.y) - erb[j].e.y * (r[0]-erb[j].b.x);
      }
      else{
        der[k0].x = 0.0; der[k0].y = 0.0; der[k0].z = 0.0;
      }
    }
    return 0;
  }

  return 1;

}

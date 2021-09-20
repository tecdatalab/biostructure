/*********************************************************************
*                           L I B _ S B A                            *
**********************************************************************
* Library is part of the Situs package URL: situs.biomachina.org     *
* (c) Mirabela Rusu, Stefan Birmanns, and Willy Wriggers, 2011       *
**********************************************************************
*                                                                    *
* Basic support routines for C++ programs derived from Sculptor SVT  *
*                                                                    *
**********************************************************************
* See legal statement for terms of distribution                      *
*********************************************************************/

#include "lib_sba.h"
#include <stdlib.h>


///////////////////////////////////////////////////////////////////////////////
// SVT_SWAP
///////////////////////////////////////////////////////////////////////////////


short int swaptest = 0x0102;

/**
 * are the CPU registers in Big-Endian-Format?
 * \return true if big endian cpu
 */
bool svt_bigEndianMachine(void)
{
  char *p;

  p = (char *)&swaptest;                          /* Zeiger auf Wert 0x0102 */
  return (*p == 0x01);              /* erstes Byte 0x01 => Big-Endian-Format */
}

/**
 * swap a double value
 * \param pValue pointer to the double variable
 */
void svt_swapDouble(double *pValue)
{
  char *cp, temp;

  cp = (char *)pValue;
  svt_swap(cp, cp + 7);
  cp++;
  svt_swap(cp, cp + 5);
  cp++;
  svt_swap(cp, cp + 3);
  cp++;
  svt_swap(cp, cp + 1);
}

/**
 * swap a Real32 value
 * \param pValue pointer to the Real32 variable
 */
void svt_swapReal32(Real32 *pValue)
{
  char *cp, temp;

  cp = (char *)pValue;
  svt_swap(cp, cp + 3);
  cp++;
  svt_swap(cp, cp + 1);
}

/**
 * swap a Int16 value
 * \param pValue pointer to the Int16 variable
 */
void svt_swapInt16(Int16 *pValue)
{
  char *cp, temp;

  cp = (char *)pValue;
  svt_swap(cp, cp + 1);
}

/**
 * swap a Int32 value
 * \param pValue pointer to the Int32 variable
 */
void svt_swapInt32(Int32 *pValue)
{
  char *cp, temp;

  cp = (char *)pValue;
  svt_swap(cp, cp + 3);
  cp++;
  svt_swap(cp, cp + 1);
}


///////////////////////////////////////////////////////////////////////////////
// SVT_TIME
///////////////////////////////////////////////////////////////////////////////

static long beginning;

#ifdef TRUE_WIN32
#  include <sys/timeb.h>
#else
#  include <sys/time.h>
#  include <unistd.h>
#  include <math.h>
#endif

#if defined(SVR4) && !defined(sun)
#define GETTIMEOFDAY(x) gettimeofday(x)
#else
#define GETTIMEOFDAY(x) gettimeofday(x,NULL)
#endif

#ifdef TRUE_WIN32
unsigned long svt_getToD(void)
{
  SYSTEMTIME st;
  GetSystemTime(&st);
  long msec = st.wMilliseconds;

  return GetTickCount();
};
#else
unsigned long svt_getToD(void)
{
  struct timeval now;
  GETTIMEOFDAY(&now);

  return (now.tv_sec * 1000) + (now.tv_usec / 1000);
};
#endif

int svt_getElapsedTime(void)
{
  int elap_time = svt_getToD() - beginning;

  return elap_time;
};

void svt_sleep(unsigned uiMilliSeconds)
{
#ifdef TRUE_WIN32
  Sleep(uiMilliSeconds);
#else
  usleep(uiMilliSeconds * 1000);
#endif

  return;
}

#ifdef TRUE_WIN32
#define strcasecmp stricmp
#endif

///////////////////////////////////////////////////////////////////////////////
// SVT_RND see artistic license below
///////////////////////////////////////////////////////////////////////////////

/* Period parameters */
#define N 624
#define M 397
#define MATRIX_A 0x9908b0df   /* constant vector a */
#define UPPER_MASK 0x80000000 /* most significant w-r bits */
#define LOWER_MASK 0x7fffffff /* least significant r bits */

/* Tempering parameters */
#define TEMPERING_MASK_B 0x9d2c5680
#define TEMPERING_MASK_C 0xefc60000
#define TEMPERING_SHIFT_U(y)  (y >> 11)
#define TEMPERING_SHIFT_S(y)  (y << 7)
#define TEMPERING_SHIFT_T(y)  (y << 15)
#define TEMPERING_SHIFT_L(y)  (y >> 18)

/* global variables */

static unsigned long g_mt[N]; /* the array for the state vector  */
static int g_mti = N + 1; /* g_mti==N+1 means g_mt[N] is not initialized */



/*===================================================================*/
void svt_sgenrand(unsigned long seed)
{
  /* Initialization of array with a seed. Theoretically,               */
  /* there are 2^19937-1 possible states as an intial state.           */
  /* This function allows to choose any of 2^19937-1 ones.             */
  /* Essential bits in "seed_array[]" is following 19937 bits:         */
  /* (seed_array[0]&UPPER_MASK), seed_array[1], ..., seed_array[N-1].  */
  /* (seed_array[0]&LOWER_MASK) is discarded.                          */
  /* Theoretically,                                                    */
  /* (seed_array[0]&UPPER_MASK), seed_array[1], ..., seed_array[N-1]   */
  /* can take any values except all zeros.                             */

  int i;

  for (i = 0; i < N; i++) {
    g_mt[i] = seed & 0xffff0000;
    seed = 69069 * seed + 1;
    g_mt[i] |= (seed & 0xffff0000) >> 16;
    seed = 69069 * seed + 1;
  }
  g_mti = N;
}



double svt_genrand()
{
  /* real random number generator */

  unsigned long y;
  static unsigned long mag01[2] = {0x0, MATRIX_A};

  if (g_mti >= N) { /* generate N words at one time */
    int kk;
    if (g_mti == N + 1) /* if sgenrand() has not been called, */
      svt_sgenrand(4357); /* a default initial seed is used   */

    for (kk = 0; kk < N - M; kk++) {
      y = (g_mt[kk] & UPPER_MASK) | (g_mt[kk + 1] & LOWER_MASK);
      g_mt[kk] = g_mt[kk + M] ^ (y >> 1) ^ mag01[y & 0x1];
    }
    for (; kk < N - 1; kk++) {
      y = (g_mt[kk] & UPPER_MASK) | (g_mt[kk + 1] & LOWER_MASK);
      g_mt[kk] = g_mt[kk + (M - N)] ^ (y >> 1) ^ mag01[y & 0x1];
    }
    y = (g_mt[N - 1] & UPPER_MASK) | (g_mt[0] & LOWER_MASK);
    g_mt[N - 1] = g_mt[M - 1] ^ (y >> 1) ^ mag01[y & 0x1];
    g_mti = 0;
  }

  y = g_mt[g_mti++];
  y ^= TEMPERING_SHIFT_U(y);
  y ^= TEMPERING_SHIFT_S(y) & TEMPERING_MASK_B;
  y ^= TEMPERING_SHIFT_T(y) & TEMPERING_MASK_C;
  y ^= TEMPERING_SHIFT_L(y);

  return ((double)y * 2.3283064370807974e-10);
}
/*  The "Artistic License" relates only to SVT_RND

                                Preamble

The intent of this document is to state the conditions under which a
Package may be copied, such that the Copyright Holder maintains some
semblance of artistic control over the development of the package,
while giving the users of the package the right to use and distribute
the Package in a more-or-less customary fashion, plus the right to make
reasonable modifications.

Definitions:

    - "Package" refers to the collection of files distributed by the
      Copyright Holder, and derivatives of that collection of files
      created through textual modification.

    - "Standard Version" refers to such a Package if it has not been
      modified, or has been modified in accordance with the wishes of
      the Copyright Holder as specified below.

    - "Copyright Holder" is whoever is named in the copyright or
      copyrights for the package.

    - "You" is you, if you're thinking about copying or distributing
      this Package.

    - "Reasonable copying fee" is whatever you can justify on the basis
      of media cost, duplication charges, time of people involved, and
      so on. (You will not be required to justify it to the Copyright
      Holder, but only to the computing community at large as a market
      that must bear the fee.)

    - "Freely Available" means that no fee is charged for the item
      itself, though there may be fees involved in handling the item.
      It also means that recipients of the item may redistribute it
      under the same conditions they received it.

1.  You may make and give away verbatim copies of the source form of the
    Standard Version of this Package without restriction, provided that
    you duplicate all of the original copyright notices and associated
    disclaimers.

2.  You may apply bug fixes, portability fixes and other modifications
    derived from the Public Domain or from the Copyright Holder. A
    Package modified in such a way shall still be considered the
    Standard Version.

3.  You may otherwise modify your copy of this Package in any way,
    provided that you insert a prominent notice in each changed file
    stating how and when you changed that file, and provided that you do
    at least ONE of the following:

    a) place your modifications in the Public Domain or otherwise make
       them Freely Available, such as by posting said modifications to
       Usenet or an equivalent medium, or placing the modifications on a
       major archive site such as uunet.uu.net, or by allowing the
       Copyright Holder to include your modifications in the Standard
       Version of the Package.

    b) use the modified Package only within your corporation or
       organization.

    c) rename any non-standard executables so the names do not conflict
       with standard executables, which must also be provided, and
       provide a separate manual page for each non-standard executable
       that clearly documents how it differs from the Standard Version.

    d) make other distribution arrangements with the Copyright Holder.

4.  You may distribute the programs of this Package in object code or
    executable form, provided that you do at least ONE of the following:

    a) distribute a Standard Version of the executables and library
       files, together with instructions (in the manual page or
       equivalent) on where to get the Standard Version.

    b) accompany the distribution with the machine-readable source of
       the Package with your modifications.

    c) give non-standard executables non-standard names, and clearly
       document the differences in manual pages (or equivalent),
       together with instructions on where to get the Standard Version.

    d) make other distribution arrangements with the Copyright Holder.

5.  You may charge a reasonable copying fee for any distribution of this
    Package. You may charge any fee you choose for support of this
    Package. You may not charge a fee for this Package itself. However,
    you may distribute this Package in aggregate with other (possibly
    commercial) programs as part of a larger (possibly commercial)
    s.ftware distribution provided that you do not advertise this
    Package as a product of your own.

6.  The scripts and library files supplied as input to or produced as
    output from the programs of this Package do not automatically fall
    under the copyright of this Package, but belong to whomever
    generated them, and may be sold commercially, and may be aggregated
    with this Package.

7.  C subroutines (or comparably compiled subroutines in other
    languages) supplied by you and linked into this Package shall not be
    considered part of this Package.

8.  Aggregation of this Package with a commercial distribution is always
    permitted provided that the use of this Package is embedded; that
    is, when no overt attempt is made to make this Package's interfaces
    visible to the end user of the commercial distribution. Such use
    shall not be construed as a distribution of this Package.

9.  The name of the Copyright Holder may not be used to endorse or
    promote products derived from this s.ftware without specific prior
    written permission.

10. THIS PACKAGE IS PROVIDED "AS IS" AND WITHOUT ANY EXPRESS OR IMPLIED
    WARRANTIES, INCLUDING, WITHOUT LIMITATION, THE IMPLIED WARRANTIES OF
    MERCHANTIBILITY AND FITNESS FOR A PARTICULAR PURPOSE.

                                The End

*/

///////////////////////////////////////////////////////////////////////////////
// SVT_RANDOM
///////////////////////////////////////////////////////////////////////////////

//
// some static stuff
//
static svt_randomGenerator s_eRandomGenerator = SPLUS;

static int s_bNativeInitialized = 0;
static int s_bRanduInialized = 0;
static int s_bSPlusInitialized = 0;

static unsigned s_uLastCon_SPLUS = 1803752341u;
static unsigned s_uLastTaus_SPLUS = 3697165728u;
static unsigned s_uLast_RANDU = 1803752341u;


//
// constants
//
const unsigned c_uLargest_SPLUS = 0xffffffff;
const unsigned c_uLargest_RANDU = 0x8fffffff;
const double   c_dPi = 3.141592654;

//
// internal prototypes
//
static unsigned svt_rand_SPLUS();
static unsigned svt_rand_RANDU();
static unsigned svt_rand_NATIVE();

static void svt_ranSeed_SPLUS(unsigned = 0, unsigned = 0);
static void svt_ranSeed_RANDU(unsigned = 0);
static void svt_ranSeed_NATIVE(unsigned = 0);


//
// Implementation
//
void svt_setRandomGenerator(svt_randomGenerator e)
{
  s_eRandomGenerator = e;
}


svt_randomGenerator svt_getRandomGenerator()
{
  return s_eRandomGenerator;

}

void svt_ranSeed(unsigned con, unsigned taus)
{
  switch (svt_getRandomGenerator()) {
    case (SPLUS):
      svt_ranSeed_SPLUS(con, taus);
      break;
    case (RANDU):
      svt_ranSeed_RANDU(con);
      break;
    default:
      svt_ranSeed_NATIVE(con);
      break;
  }
}

void svt_ranSeedAll(unsigned con, unsigned taus)
{
  svt_ranSeed_SPLUS(con, taus);
  svt_ranSeed_RANDU(con);
  svt_ranSeed_NATIVE(con);
}


unsigned svt_ranLargest()
{
  switch (svt_getRandomGenerator()) {
    case (SPLUS):
      return c_uLargest_SPLUS;
    case (RANDU):
      return c_uLargest_RANDU;
    default:
      return RAND_MAX;
  }
}



unsigned svt_rand_SPLUS()
{

  const unsigned lambda = 69069;

  if (!s_bSPlusInitialized)
    svt_ranSeed_SPLUS();

  s_uLastCon_SPLUS = s_uLastCon_SPLUS * lambda;
  s_uLastTaus_SPLUS ^= s_uLastTaus_SPLUS >> 15;
  s_uLastTaus_SPLUS ^= s_uLastTaus_SPLUS << 17;

  return s_uLastTaus_SPLUS ^ s_uLastCon_SPLUS;

}


unsigned svt_rand_RANDU()
{
  if (!s_bRanduInialized)
    svt_ranSeed_RANDU();

  unsigned v = ((1 << 16) + 3) * s_uLast_RANDU;
  v = v & 0x8fffffffu;
  s_uLast_RANDU = v;

  return s_uLast_RANDU;

}

unsigned svt_rand_NATIVE()
{

  if (!s_bNativeInitialized) {
    svt_ranSeed_NATIVE();
  }
  return rand();

}

unsigned svt_rand()
{
  switch (svt_getRandomGenerator()) {
    case (SPLUS):
      return svt_rand_SPLUS();
    case (RANDU):
      return svt_rand_RANDU();
    default:
      return svt_rand_NATIVE();
  }

}


double svt_ranUni()
{
  return svt_rand() / double(svt_ranLargest());
}



double svt_ranNormal()
{

  double r1, r2;
  static int s_iReserve = 0;
  static double s_dReserveValue = 0;

  if (s_iReserve) {
    s_iReserve = 0;
    return s_dReserveValue;
  }

  r1 = sqrt(-2 * log(svt_ranUni()));
  r2 = svt_ranUni();

  s_iReserve = 1;
  s_dReserveValue = r1 * sin(2 * c_dPi * r2);

  return r1 * cos(2 * c_dPi * r2);

}

double svt_ranNormal(double mu, double sigma)
{
  return mu + svt_ranNormal() * sigma;

}

/**
 * generate a random number following a Cauchy Distribution ~ resembles with a normal distribution just wider at the tail
 * the ratio of a two normal distributed variables follow a Cauchy distribution of average 0 and sd 1
 * \return a random number following the Cauchy distribution of average 0 and standard deviation 1
 */
double svt_ranCauchy()
{
  double fN1 = svt_ranNormal();
  double fN2 = svt_ranNormal();

  double fRatio = fN1 / fN2;

  return  fRatio;
}

/**
 * \return a random number following the Cauchy distribution of average mu and standard deviation sigma
 */
double svt_ranCauchy(double mu, double sigma)
{
  return mu + svt_ranCauchy() * sigma;
}

double svt_ranWeibull()
{
  printf("svt_ranWeibull: Not implemented yet!\n");
  return 0;
}


double svt_ranHjorth()
{
  printf("svt_ranHjorth: Not implemented yet!\n");
  return 0;
}



double svt_ranExp()
{
  printf("svt_ranExp: Not implemented yet!\n");
  return 0;
}



unsigned svt_ranBernoulli()
{
  printf("svt_ranBernoulli: Not implemented yet!\n");
  return 0;
}



unsigned svt_ranBinomial()
{
  printf("svt_ranBinomial: Not implemented yet!\n");
  return 0;
}



unsigned svt_ranPoisson()
{
  printf("svt_ranPoisson: Not implemented yet!\n");
  return 0;
}



static void svt_ranSeed_SPLUS(unsigned con, unsigned taus)
{

  if (con == 0)
    con = time(NULL);

  taus = con * 2;

  if (!(con % 2))
    con++;

  s_uLastCon_SPLUS = con;
  s_uLastTaus_SPLUS = taus;
  //printf("Initialized SPLUS with %d, %d\n", con, taus);

  s_bSPlusInitialized = 1;
}


static void svt_ranSeed_RANDU(unsigned start)
{
  if (start == 0)
    start = time(NULL);

  s_uLast_RANDU = start;

  //printf("Initialized RANDU with %d\n", start);
  s_bRanduInialized = 1;
}


static void svt_ranSeed_NATIVE(unsigned start)
{
  if (start == 0)
    start = time(NULL);

  srand(start);

  //printf("Initialized NATIVE with %d\n", start);

  s_bNativeInitialized = 1;
}

///////////////////////////////////////////////////////////////////////////////
// SVT_BOND
///////////////////////////////////////////////////////////////////////////////


/**
 * Constructor
 * create an bond between atom a and b
 * \param pA pointer to first svt_atom object
 * \param pB pointer to second svt_atom object
 */
svt_point_cloud_bond::svt_point_cloud_bond(svt_point_cloud_atom *pA, svt_point_cloud_atom *pB, int iIndexA, int iIndexB) :
  m_pAtomA(pA),
  m_pAtomB(pB),
  m_iIndexA(iIndexA),
  m_iIndexB(iIndexB)
{
}
svt_point_cloud_bond::~svt_point_cloud_bond()
{
}

/**
 * set the first atom
 * \param pA pointer to svt_atom object
 */
void svt_point_cloud_bond::setAtomA(svt_point_cloud_atom *pA)
{
  m_pAtomA = pA;
}
/**
 * get the first atom
 * \return pointer to svt_atom object
 */
svt_point_cloud_atom *svt_point_cloud_bond::getAtomA()
{
  return m_pAtomA;
}
/**
 * set the second atom
 * \param pB pointer to svt_atom object
 */
void svt_point_cloud_bond::setAtomB(svt_point_cloud_atom *pB)
{
  m_pAtomB = pB;
}
/**
 * get the second atom
 * \return pointer to svt_atom object
 */
svt_point_cloud_atom *svt_point_cloud_bond::getAtomB()
{
  return m_pAtomB;
}
/**
 * get the index of atom a
 */
int svt_point_cloud_bond::getIndexA()
{
  return m_iIndexA;
}
/**
 * get the index of atom b
 */
int svt_point_cloud_bond::getIndexB()
{
  return m_iIndexB;
}
/**
 * set the index of atom a
 */
void svt_point_cloud_bond::setIndexA(int iIndexA)
{
  m_iIndexA = iIndexA;
}
/**
 * set the index of atom b
 */
void svt_point_cloud_bond::setIndexB(int iIndexB)
{
  m_iIndexB = iIndexB;
}
///////////////////////////////////////////////////////////////////////////////
// SVT_ATOM
///////////////////////////////////////////////////////////////////////////////

/**
 * Constructor
 */
svt_point_cloud_atom::svt_point_cloud_atom() :
  m_iModel(0),
  m_iPDBIndex(0),
  m_bWater(false),
  m_bHetAtm(false),
  m_aBondList(0),
  m_aBondListI(0)
{
  sprintf(m_aName, "QP");
  // remoteness indicator
  m_cRemoteness = 'D';
  // branch
  m_cBranch = 'B';
  // radius
  m_fRadius = 1.0;
  // alternate location identifier
  m_cAltLoc = ' ';
  // atom resname (i.e. "ALA")
  sprintf(m_aResname, "QPD");
  // short resname name
  m_cShortResname = '-';
  // chain id
  m_cChainID = 'A';
  // residue sequence number
  m_iResSeq = 0;
  // ordinal residue sequence number
  m_iOrdResSeq = 0;
  // iCode
  m_cICode = ' ';
  // the occupancy
  m_fOccupancy = 1.0;
  // the temperature factor
  m_fTempFact = 1.0;
  // note
  m_pNote[0] = 0;
  // segment id
  m_pSegID[0] = 0;
  // element symbol
  sprintf(m_pElement, "QP");
  // charge
  m_pCharge[0] = 0;

  // secondary structure information
  // H   Alpha helix
  // G   3-10 helix
  // I   PI-helix
  // E   Extended conformation
  // B   Isolated bridge
  // T   Turn
  // C   Coil (none of the above)
  // N   Information not available
  m_cSecStruct = 'C';

  // number of residues in the secondary structure this atom belongs to
  m_iSecStructNumResidues = 0;

  // relative atomic mass
  m_fMass = 1.0;

  m_bIsSelected = false;

  if (m_oAAList.size() == 0) { // the list was not created
    svt_resname oAA;
    sprintf(oAA.m_a3LetCode, "ALA");
    oAA.m_c1LetCode = 'A';
    m_oAAList.push_back(oAA);
    sprintf(oAA.m_a3LetCode, "ARG");
    oAA.m_c1LetCode = 'R';
    m_oAAList.push_back(oAA);
    sprintf(oAA.m_a3LetCode, "ASN");
    oAA.m_c1LetCode = 'N';
    m_oAAList.push_back(oAA);
    sprintf(oAA.m_a3LetCode, "ASP");
    oAA.m_c1LetCode = 'D';
    m_oAAList.push_back(oAA);
    sprintf(oAA.m_a3LetCode, "CYS");
    oAA.m_c1LetCode = 'C';
    m_oAAList.push_back(oAA);
    sprintf(oAA.m_a3LetCode, "GLU");
    oAA.m_c1LetCode = 'E';
    m_oAAList.push_back(oAA);
    sprintf(oAA.m_a3LetCode, "GLN");
    oAA.m_c1LetCode = 'Q';
    m_oAAList.push_back(oAA);
    sprintf(oAA.m_a3LetCode, "GLY");
    oAA.m_c1LetCode = 'G';
    m_oAAList.push_back(oAA);
    sprintf(oAA.m_a3LetCode, "HIS");
    oAA.m_c1LetCode = 'H';
    m_oAAList.push_back(oAA);
    sprintf(oAA.m_a3LetCode, "ILE");
    oAA.m_c1LetCode = 'I';
    m_oAAList.push_back(oAA);
    sprintf(oAA.m_a3LetCode, "LEU");
    oAA.m_c1LetCode = 'L';
    m_oAAList.push_back(oAA);
    sprintf(oAA.m_a3LetCode, "LYS");
    oAA.m_c1LetCode = 'K';
    m_oAAList.push_back(oAA);
    sprintf(oAA.m_a3LetCode, "MET");
    oAA.m_c1LetCode = 'M';
    m_oAAList.push_back(oAA);
    sprintf(oAA.m_a3LetCode, "PHE");
    oAA.m_c1LetCode = 'F';
    m_oAAList.push_back(oAA);
    sprintf(oAA.m_a3LetCode, "PRO");
    oAA.m_c1LetCode = 'P';
    m_oAAList.push_back(oAA);
    sprintf(oAA.m_a3LetCode, "SER");
    oAA.m_c1LetCode = 'S';
    m_oAAList.push_back(oAA);
    sprintf(oAA.m_a3LetCode, "THR");
    oAA.m_c1LetCode = 'T';
    m_oAAList.push_back(oAA);
    sprintf(oAA.m_a3LetCode, "TRP");
    oAA.m_c1LetCode = 'W';
    m_oAAList.push_back(oAA);
    sprintf(oAA.m_a3LetCode, "TYR");
    oAA.m_c1LetCode = 'Y';
    m_oAAList.push_back(oAA);
    sprintf(oAA.m_a3LetCode, "VAL");
    oAA.m_c1LetCode = 'V';
    m_oAAList.push_back(oAA);
  }

};

/**
 * Destructor
 */
svt_point_cloud_atom::~svt_point_cloud_atom()
{
};

/**
 * Set the model this atom belongs to
 * \param iModel number of model
 */
void svt_point_cloud_atom::setModel(unsigned int iModel)
{
  m_iModel = iModel;
};

/**
 * Get the model this atom belongs to
 * \return number of model
 */
unsigned int svt_point_cloud_atom::getModel() const
{
  return m_iModel;
};

/**
 * Set the atom type. Automatically strips spaces in front and after the atom name.
 * \param pType pointer to the char array with the name information
 */
void svt_point_cloud_atom::setName(const char *pName)
{
  strcpy(m_aName, pName);

  char cElement = m_aName[1];

  // strip spaces in front
  while (m_aName[0] == ' ') {
    m_aName[0] = m_aName[1];
    m_aName[1] = m_aName[2];
    m_aName[2] = m_aName[3];
    m_aName[3] = m_aName[4];
  }
  // strip spaces at the end
  for (int i = 4; i >= 0; i--)
    if (m_aName[i] == 32)
      m_aName[i] = 0;

  // simple vdw radii. If this is not exact enough uncomment the "correct" vdw radii below. But you will get much more bonds due to
  // the bigger radii.
  m_fRadius = 1.5f;
  switch (cElement) {
    case 'H' :
      m_fRadius = 1.00f;
      break;
    case 'C' :
      m_fRadius = 1.50f;
      break;
    case 'N' :
      m_fRadius = 1.40f;
      break;
    case 'O' :
      m_fRadius = 1.30f;
      break;
    case 'F' :
      m_fRadius = 1.20f;
      break;
    case 'S' :
      m_fRadius = 1.90f;
      break;
  }

  /*
   // standard van der Waals radii are assigned to the common elements. They are taken from Bondi, J.Phys.Chem., 68, 441, 1964. Other elements are assigned van der Waals radii of 2.0A.
   float rad = 2.0f;
   RAD("Ag", 1.72f)
   RAD("Ar", 1.88f)
   RAD("As", 1.85f)
   RAD("Au", 1.66f)
   RAD("Br", 1.85f)
   RAD("C",  1.70f)
   RAD("Cd", 1.58f)
   RAD("Cl", 1.75f)
   RAD("Cu", 1.40f)
   RAD("F",  1.47f)
   RAD("Ga", 1.87f)
   RAD("H",  1.20f)
   RAD("He", 1.40f)
   RAD("Hg", 1.55f)
   RAD("I",  1.98f)
   RAD("In", 1.93f)
   RAD("K",  2.75f)
   RAD("Kr", 2.02f)
   RAD("Li", 1.82f)
   RAD("Mg", 1.73f)
   RAD("N",  1.55f)
   RAD("Na", 2.27f)
   RAD("Ne", 1.54f)
   RAD("Ni", 1.63f)
   RAD("O",  1.52f)
   RAD("P",  1.80f)
   RAD("Pb", 2.02f)
   RAD("Pd", 1.63f)
   RAD("Pt", 1.72f)
   RAD("S",  1.80f)
   RAD("Se", 1.90f)
   RAD("Si", 2.10f)
   RAD("Sn", 2.17f)
   RAD("Te", 2.06f)
   RAD("Tl", 1.96f)
   RAD("U",  1.86f)
   RAD("Xe", 2.16f)
   RAD("Zn", 1.39f)
   */

  m_aName[4] = 0;
};

/**
 * get the atom name
 * \return pointer to the char array with the name information
 */
const char *svt_point_cloud_atom::getName() const
{
  return m_aName;
};

/**
 * Set the PDB-file index of the atom
 * \param iPDBIndex index of the atom
 */
void svt_point_cloud_atom::setPDBIndex(unsigned int iPDBIndex)
{
  m_iPDBIndex = iPDBIndex;
};
/**
 * Get the PDB-file index of the atom
 * \return index of the atom
 */
unsigned int svt_point_cloud_atom::getPDBIndex() const
{
  return m_iPDBIndex;
};

/**
 * get the atom remoteness indicator
 * \return char with the remoteness information (greek letters, A=alpha, B=beta,...)
 */
char svt_point_cloud_atom::getRemoteness() const
{
  return m_cRemoteness;
};

/**
 * set the atom remoteness indicator
 * \param cRemoteness char with the remoteness information (greek letters, A=alpha, B=beta,...)
 */
void svt_point_cloud_atom::setRemoteness(char cRemoteness)
{
  m_cRemoteness = cRemoteness;
};

/**
 * get the atom branch information
 * \return char with the branch information
 */
char svt_point_cloud_atom::getBranch() const
{
  return m_cBranch;
};
/**
 * set the atom branch information
 * \param cBranch char with the branch information
 */
void svt_point_cloud_atom::setBranch(char cBranch)
{
  m_cBranch = cBranch;
};

/**
 * set the atom alternate location indicator
 * \param cAltLoc character with the alternate location indicator
 */
void svt_point_cloud_atom::setAltLoc(char cAltLoc)
{
  m_cAltLoc = cAltLoc;
};

/**
 * get the atom alternate location indicator
 * \return character with the alternate location indicator
 */
char svt_point_cloud_atom::getAltLoc() const
{
  return m_cAltLoc;
};

/**
 * set the atom resname
 * \param pResname pointer to the char array with the residue name information
 */
void svt_point_cloud_atom::setResname(const char *pResname)
{
  strcpy(m_aResname, pResname);
  m_aResname[3] = 0;

  if (strcasecmp(m_aResname, "tip3") == 0 || strcasecmp(m_aResname, "hoh") == 0 || strcasecmp(m_aResname, "h2o") == 0)
    m_bWater = true;
  else
    m_bWater = false;

  unsigned int iSize = m_oAAList.size();
  for (unsigned int iIndex = 0; iIndex < iSize; iIndex++) {
    if (strcasecmp(m_aResname, m_oAAList[iIndex].m_a3LetCode) == 0) {
      m_cShortResname = m_oAAList[iIndex].m_c1LetCode;
      return;
    }
  }
  m_cShortResname = '-'; // is empty if nothing found

};

/**
 * get the atom resname
 * \return pointer to the char array with the residue name information
 */
const char *svt_point_cloud_atom::getResname() const
{
  return m_aResname;
};

/**
 * set the short atom resname
 * \param pResname a char with the short "1letter" residue name information
 */
void svt_point_cloud_atom::setShortResname(const char cResname)
{
  m_cShortResname = cResname;
};

/**
 * get the short atom resname
 * \return a char with the short "1letter" residue name information
 */
char svt_point_cloud_atom::getShortResname() const
{
  return m_cShortResname;
};

/**
 * set the chain id
 * \param cChainID character with the chain id
 */
void svt_point_cloud_atom::setChainID(char cChainID)
{
  m_cChainID = cChainID;
};

/**
 * get the chain id
 * \return character with the chain id
 */
char svt_point_cloud_atom::getChainID() const
{
  return m_cChainID;
};

/**
 * set the ordinal chain id
 * \param iOrdChainID character with the chain id
 */
void svt_point_cloud_atom::setOrdChainID(int iOrdChainID)
{
  m_iOrdChainID = iOrdChainID;
};

/**
 * get the ordinal chain id
 * \return int with the chain id
 */
int svt_point_cloud_atom::getOrdChainID() const
{
  return m_iOrdChainID;
};

/**
 * set the residue sequence number
 * \param iResSeq residue sequence number
 */
void svt_point_cloud_atom::setResidueSeq(int iResSeq)
{
  m_iResSeq = iResSeq;
};

/**
 * get the residue sequence number
 * \return residue sequence number
 */
int svt_point_cloud_atom::getResidueSeq() const
{
  return m_iResSeq;
};

/**
 * set the ordinal residue sequence number
 * \param iOrdResSeq ordinal residue sequence number
 */
void svt_point_cloud_atom::setOrdResidueSeq(int iOrdResSeq)
{
  m_iOrdResSeq = iOrdResSeq;
};
/**
 * get the ordinal residue sequence number
 * \return ordinal residue sequence number
 */
int svt_point_cloud_atom::getOrdResidueSeq() const
{
  return m_iOrdResSeq;
};

/**
 * set the iCode (code for insertion of residues)
 * \param cICode char with the iCode
 */
void svt_point_cloud_atom::setICode(char cICode)
{
  m_cICode = cICode;
};

/**
 * get the iCode (code for insertion of residues)
 * \return char with the iCode
 */
char svt_point_cloud_atom::getICode() const
{
  return m_cICode;
};

/**
 * set the occupancy
 * \param fOccupancy the occupancy
 */
void svt_point_cloud_atom::setOccupancy(float fOccupancy)
{
  m_fOccupancy = fOccupancy;
};

/**
 * get the occupancy
 * \return the occupancy
 */
float svt_point_cloud_atom::getOccupancy() const
{
  return m_fOccupancy;
};

/**
 * set the temperature factor
 * \param fTempFact the temperature factor
 */
void svt_point_cloud_atom::setTempFact(float fTempFact)
{
  m_fTempFact = fTempFact;
};

/**
 * get the temperature factor
 * \return the temperature factor
 */
float svt_point_cloud_atom::getTempFact() const
{
  return m_fTempFact;
};

/**
 * set the note
 * \param pNote pointer to char array with at least 3 chars for the note
 */
void svt_point_cloud_atom::setNote(const char *pNote)
{
  m_pNote[0] = pNote[0];
  m_pNote[1] = pNote[1];
  m_pNote[2] = pNote[2];
  m_pNote[3] = 0;

  if (m_pNote[2] == 0)
    m_pNote[2] = ' ';
  if (m_pNote[1] == 0)
    m_pNote[1] = ' ';
  if (m_pNote[0] == 0)
    m_pNote[0] = ' ';
};

/**
 * get the note
 * \return pointer to char array with the note
 */
const char *svt_point_cloud_atom::getNote() const
{
  return m_pNote;
};

/**
 * set the segment id
 * \param pSegID pointer to char array with at least 4 chars for the segment id
 */
void svt_point_cloud_atom::setSegmentID(const char *pSegID)
{
  m_pSegID[0] = pSegID[0];
  m_pSegID[1] = pSegID[1];
  m_pSegID[2] = pSegID[2];
  m_pSegID[3] = pSegID[3];
  m_pSegID[4] = 0;

  if (m_pSegID[3] == 0)
    m_pSegID[3] = ' ';
  if (m_pSegID[2] == 0)
    m_pSegID[2] = ' ';
  if (m_pSegID[1] == 0)
    m_pSegID[1] = ' ';
  if (m_pSegID[0] == 0)
    m_pSegID[0] = ' ';
};

/**
 * get the segment id
 * \return pointer to char array with the segment id
 */
const char *svt_point_cloud_atom::getSegmentID() const
{
  return m_pSegID;
};

/**
 * set the element symbol
 * \param pElement pointer to char array with at least 2 chars for the element symbol
 */
void svt_point_cloud_atom::setElement(const char *pElement)
{
  m_pElement[0] = pElement[0];
  m_pElement[1] = pElement[1];
  m_pElement[2] = 0;

  //Mirabela added = 13 Carrige return
  if (m_pElement[1] == 0 || m_pElement[1] == 10 || m_pElement[1] == 13)
    m_pElement[1] = ' ';
  if (m_pElement[0] == 0 || m_pElement[0] == 10 || m_pElement[0] == 13)
    m_pElement[0] = ' ';

};

/**
 * get the element symbol
 * \return pointer to char array with the element symbol
 */
const char *svt_point_cloud_atom::getElement() const
{
  return m_pElement;
};

/**
 * set the charge
 * \param pCharge pointer to char array with at least 2 chars for the charge
 */
void svt_point_cloud_atom::setCharge(const char *pCharge)
{
  m_pCharge[0] = pCharge[0];
  m_pCharge[1] = pCharge[1];
  m_pCharge[2] = 0;

  if (m_pCharge[1] == 0 || m_pCharge[1] == 10 || m_pCharge[1] == 13)
    m_pCharge[1] = ' ';
  if (m_pCharge[0] == 0 || m_pCharge[0] == 10 || m_pCharge[0] == 13)
    m_pCharge[0] = ' ';
};

/**
 * get the charge
 * \return pointer to char array with the charge
 */
const char *svt_point_cloud_atom::getCharge() const
{
  return m_pCharge;
};

#define ISNUMBER(X) (X > 47 && X < 58)

/**
 * is the atom a hydrogen atom?
 * \return true if the the atom is an hydrogen atom
 */
bool svt_point_cloud_atom::isHydrogen() const
{
  if (m_aName[0] == 'H' || (ISNUMBER(m_aName[0]) && m_cRemoteness == 'H') || (ISNUMBER(m_aName[0]) && ISNUMBER(m_cRemoteness) && m_cBranch == 'H'))
    return true;
  else
    return false;
}
/**
 * is the atom a QPDB codebook vector?
 * \return true if the the atom is a CV
 */
bool svt_point_cloud_atom::isQPDB() const
{
  if (m_aName[0] == 'Q' && m_cRemoteness == 'P' && m_cBranch == 'D')
    return true;
  else
    return false;
}

/**
 * is the atom part of a water molecule?
 * \return true if the atom os part of a water molecule
 */
bool svt_point_cloud_atom::isWater() const
{
  return m_bWater;
}

/**
 * is the atom a carbon alpha?
 * \return true if the atom a CA
 */
bool svt_point_cloud_atom::isCA() const
{
  if (!m_bHetAtm && (m_aName[0] == 'C' && m_cRemoteness == 'A'))
    return true;
  else
    return false;
}

/**
 * is the atom on the Backbone? (def Backbone: CAlpha (CA), Carbon(C) witout Remoteness with his Oxygen, Nitrogen(N) witout Remoteness and the Oxygen (O)  of the previous
 * \return true if the atom on the Backbone
 */
bool svt_point_cloud_atom::isBackbone() const
{

  if (!m_bHetAtm && ((m_aName[0] == 'C' && m_cRemoteness == 'A')
                     || (m_aName[0] == 'C' && m_cRemoteness == ' ')
                     || (m_aName[0] == 'N' && m_cRemoteness == ' ')
                     || (m_aName[0] == 'O' && m_cRemoteness == ' ')
                     || (m_aName[0] == 'O' && m_cRemoteness == 'X') // COOH of the C terminus
                     || (m_aName[0] == 'P' && m_cRemoteness == ' ') //nucleotide
                     || (m_aName[0] == 'O' && m_cRemoteness == '5' && m_cBranch == '\'')
                     || (m_aName[0] == 'C' && m_cRemoteness == '5' && m_cBranch == '\'')
                     || (m_aName[0] == 'C' && m_cRemoteness == '4' && m_cBranch == '\'')
                     || (m_aName[0] == 'O' && m_cRemoteness == '4' && m_cBranch == '\'')
                     || (m_aName[0] == 'C' && m_cRemoteness == '3' && m_cBranch == '\'')
                     || (m_aName[0] == 'O' && m_cRemoteness == '3' && m_cBranch == '\'')
                     || (m_aName[0] == 'C' && m_cRemoteness == '2' && m_cBranch == '\'')
                     || (m_aName[0] == 'O' && m_cRemoteness == '2' && m_cBranch == '\'')
                     || (m_aName[0] == 'C' && m_cRemoteness == '1' && m_cBranch == '\'')
                    ))
    return true;
  else
    return false;
}

/**
 * is the atom part of a nucleotide?
 * \return true if the atom is part of a nucleotide
 */
bool svt_point_cloud_atom::isNucleotide() const
{
  return (m_aResname[0] == ' ' && m_aResname[1] == ' ');
}

/**
 * set the relative atomic mass
 * \param fMass relative atomic mass
 */
void svt_point_cloud_atom::setMass(Real64 fMass)
{
  m_fMass = fMass;
}
/**
 * get the relative atomic mass
 * \return relative atomic mass
 */
Real64 svt_point_cloud_atom::getMass() const
{
  return m_fMass;
}

/**
 * adjust the atomic mass based on a (simple) periodic table.
 * ToDo: Full periodic table
 */
void svt_point_cloud_atom::adjustMass()
{
  static const char *atom_name[] = {
    "H", "HE", "LI", "BE",  "B",  "C",  "N",  "O",  "F", "NE",
    "NA", "MG", "AL", "SI",  "P",  "S", "CL", "AR",  "K", "CA",
    "SC", "TI",  "V", "CR", "MN", "FE", "CO", "NI", "CU", "ZN"
  };
  static float atom_mass[30] = {
    1.008 ,  4.003,  6.940,  9.012, 10.810, 12.011, 14.007, 16.000,  18.998, 20.170,
    22.989, 24.305, 26.982, 28.086, 30.974, 32.060, 35.453, 39.948,  39.098, 40.078,
    44.956, 47.867, 50.942, 51.996, 54.938, 55.847, 58.933, 58.710,  63.546, 65.380
  };

  for (int i = 0; i < 30; i++) {
    if (strcmp(getName(), atom_name[i]) == 0) {
      setMass(atom_mass[i]);

      //svtout << getName() << " - " << atom_mass[i] << endl;

      return;
    }
  }

  //svtout << getName() << " - NOT FOUND, NO MASS ASSIGNED!" << endl;

};

/**
 * set the secondary structure information for this atom
 * \param cSecStruct secondary structure information
 *
 */
void svt_point_cloud_atom::setSecStruct(char cSecStruct)
{
  m_cSecStruct = cSecStruct;
}

/**
 * get the secondary structure information for this atom
 * \return secondary structure information
 *
 */
char svt_point_cloud_atom::getSecStruct()
{
  return m_cSecStruct;
}

/**
 * set the number of residues in the secondary structure this atom belongs to
 * \param iSecStructNumResidues number of residues
 *
 */
void svt_point_cloud_atom::setSecStructNumResidues(int iSecStructNumResidues)
{
  m_iSecStructNumResidues = iSecStructNumResidues;
}

/**
 * get the number of residues in the secondary structure this atom belongs to
 * \return number of residues
 *
 */
int svt_point_cloud_atom::getSecStructNumResidues()
{
  return m_iSecStructNumResidues;
}

/**
 * Get the van der waals radius
 * \return van der waals radius value (in Angstroem)
 */
Real64 svt_point_cloud_atom::getVDWRadius() const
{
  return m_fRadius;
}

/**
 * Add a bond to the bond list of this atom
 * \param rBond pointer to the svt_point_cloud_bond object which should be added to the bond list
 */
void svt_point_cloud_atom::addToBondList(svt_point_cloud_bond *pBond, unsigned int iIndex)
{
  m_aBondList.push_back(pBond);
  m_aBondListI.push_back(iIndex);
};
/**
 * Remove a bond from the bond list of this atom
 * \param rBond pointer to the svt_point_cloud_bond object which should be added to the bond list
 */
void svt_point_cloud_atom::delFromBondList(unsigned int iIndex)
{
  vector< unsigned int >::iterator pI = find(m_aBondListI.begin(), m_aBondListI.end(), iIndex);
  if (pI == m_aBondListI.end()) {
    return;
  } else {
    m_aBondListI.erase(pI);
    int iIndex = m_aBondListI.end() - pI;
    m_aBondList.erase(m_aBondList.begin() + iIndex);
  }
  return;
};

/**
 * Get bond list
 */
vector< unsigned int > svt_point_cloud_atom::getBondList()
{
  return m_aBondListI;
}

/**
 * Adjust bond list - a bond was erased from the global list and now all bond indexes higher than a certain number have to be decresed by one
 */
void svt_point_cloud_atom::adjustBondList(unsigned int iIndex)
{
  for (unsigned int i = 0; i < m_aBondListI.size(); i++)
    if (m_aBondListI[i] > iIndex)
      m_aBondListI[i] = m_aBondListI[i] - 1;
}

/**
 * deletes all bonds in the bond list
 */
void svt_point_cloud_atom::delBondList()
{
  m_aBondList.clear();
  m_aBondListI.clear();
};

/**
 * set Selected
 */
void svt_point_cloud_atom::setSelected(bool bSelected)
{
  m_bIsSelected = bSelected;
};

/**
* get Selected
*/
bool svt_point_cloud_atom::getSelected()
{
  return m_bIsSelected;
};

/**
 * set the record name : is atom of class hetatm
 */
void svt_point_cloud_atom::setHetAtm(bool bHetAtm)
{
  m_bHetAtm = bHetAtm;
};

/**
 * get the record name : is atom of class hetatm
 */
bool svt_point_cloud_atom::getHetAtm()
{
  return m_bHetAtm;
};


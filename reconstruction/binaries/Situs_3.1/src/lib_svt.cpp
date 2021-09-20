/*********************************************************************
*                           L I B _ S V T                            *
**********************************************************************
* Library is part of the Situs package URL: situs.biomachina.org     *
* (c) Mirabela Rusu, Stefan Birmanns, and Willy Wriggers, 2011       *
**********************************************************************
*                                                                    *
* Genetic algorithm related routines derived from Sculptor SVT       *
*                                                                    *
**********************************************************************
* See legal statement for terms of distribution                      *
*********************************************************************/

#include "lib_svt.h"
#include <stdlib.h>
#include <stdio.h>
#include <pthread.h>


///////////////////////////////////////////////////////////////////////////////
// SVT_CONFIG
///////////////////////////////////////////////////////////////////////////////
#ifdef TRUE_WIN32
#include <direct.h>
#define strcasecmp stricmp
#endif

#ifdef WIN32
#ifndef STATICQT
extern "C" BOOL WINAPI DllMain(HANDLE hModule, DWORD dwFunction, LPVOID lpNot)
{
  svt_registerCout(&std::cout);
  return true;
};
#endif
#endif

/**
 * Constructor
 */
svt_config::svt_config(const char *pFname) :
  m_iItems(0)
{
  for (unsigned int i = 0; i < MAXPARAM; i++) {
    m_aParameter[i].name[0] = 0;
    m_aParameter[i].wert[0] = 0;
  }


  if (pFname == NULL)
    parse(".svtrc");
  else
    parse(pFname);

}
svt_config::~svt_config()
{
}

/**
 * parse the configuration file and store all key value pairs
 * \param pFilename pointer to filename
 */
void svt_config::parse(const char *pFilename)
{
  FILE *pFile = NULL;
  char pTmp[256];
  pFile = fopen(pFilename, "r");

  if (pFile) {
    while (!feof(pFile) && !ferror(pFile) && m_iItems < MAXPARAM - 1) {
      m_iItems++;
      sscanf(pTmp, "%60s = %[^\n]", m_aParameter[m_iItems].name, m_aParameter[m_iItems].wert);
      // strip trailing spaces
      for (int i = strlen(m_aParameter[m_iItems].wert) - 1; i > 0; i--)
        if (m_aParameter[m_iItems].wert[i] == ' ')
          m_aParameter[m_iItems].wert[i] = 0;
        else
          i = 0;

      // comment?
      if (strchr(pTmp, '#') != 0)
        m_iItems--;
      //else
      //    cout << "name: " << m_aParameter[m_iItems].name << " wert: " << m_aParameter[m_iItems].wert << endl;
    }
    fclose(pFile);
  }
  //else
  //    cout << ".svtrc not found! Parameters set to default values!" << endl;
}

// teste ob NAME als Parameter gespeichert ist, wenn ja gib Feldnummer zurck, sonst 0
int svt_config::findItem(const char *name)
{
  if (name == NULL)
    return 0;

  for (int i = 1; i <= m_iItems; i++)
    if (strcasecmp(name, m_aParameter[i].name) == 0) //case
      return i;

  // kein entsprechender Eintrag gefunden
  return 0;
}

// suche Parameter und gib seinen Wert zurck, falls nicht vorhanden DEFAULT_VALUE
bool svt_config::getValue(const char *pParam, bool bDefault)
{
  int iItem = findItem(pParam);

  if (iItem != 0) {
    if (strcasecmp(m_aParameter[iItem].wert, "true") == 0)
      return true;
    else
      return false;
  } else
    return bDefault;
}

int svt_config::getValue(const char *pParam, int default_value)
{
  int iItem = findItem(pParam);

  if (iItem != 0)
    return atoi(m_aParameter[iItem].wert);
  else
    return default_value;
}

float svt_config::getValue(const char *pParam, float default_value)
{
  int iItem = findItem(pParam);

  if (iItem != 0)
    return atof(m_aParameter[iItem].wert);
  else
    return default_value;
}

double svt_config::getValue(const char *pParam, double default_value)
{
  int iItem = findItem(pParam);

  if (iItem != 0)
    return atof(m_aParameter[iItem].wert);
  else
    return default_value;
}

const char *svt_config::getValue(const char *pParam, const char *default_value)
{
  int iItem = findItem(pParam);

  if (iItem != 0)
    return m_aParameter[iItem].wert;
  else
    return default_value;
}


void svt_config::setValue(const char *pParam, bool bNew)
{
  int iItem = findItem(pParam);

  if (iItem == 0) {
    m_iItems++;
    strcpy(m_aParameter[ m_iItems ].name, pParam);
    iItem = m_iItems;
  }

  if (bNew)
    sprintf(m_aParameter[iItem].wert, "true");
  else
    sprintf(m_aParameter[iItem].wert, "false");
}

void svt_config::setValue(const char *pParam, int iValue)
{
  int iItem = findItem(pParam);

  if (iItem == 0) {
    m_iItems++;
    strcpy(m_aParameter[ m_iItems ].name, pParam);
    iItem = m_iItems;
  }

  sprintf(m_aParameter[iItem].wert, "%i", iValue);
}

void svt_config::setValue(const char *pParam, float fValue)
{
  int iItem = findItem(pParam);

  if (iItem == 0) {
    m_iItems++;
    strcpy(m_aParameter[ m_iItems ].name, pParam);
    iItem = m_iItems;
  }

  sprintf(m_aParameter[iItem].wert, "%f", fValue);
}

void svt_config::setValue(const char *pParam, double fValue)
{
  int iItem = findItem(pParam);

  if (iItem == 0) {
    m_iItems++;
    strcpy(m_aParameter[ m_iItems ].name, pParam);
    iItem = m_iItems;
  }

  sprintf(m_aParameter[iItem].wert, "%f", fValue);
}

void svt_config::setValue(const char *pParam, const char *pValue)
{
  int iItem = findItem(pParam);

  if (iItem == 0) {
    m_iItems++;
    strcpy(m_aParameter[ m_iItems ].name, pParam);
    iItem = m_iItems;
  }

  sprintf(m_aParameter[iItem].wert, "%s", pValue);
}


///////////////////////////////////////////////////////////////////////////////
// SVT_SEMAPHORE
///////////////////////////////////////////////////////////////////////////////

void __svt_fatalError(const char *s)
{
  printf("SVT - Fatal Error: %s\n", s);
  exit(1);
}

#ifndef WIN32
svt_semaphore::svt_semaphore(int)
{
  if (pthread_mutex_init(&mux, NULL) < 0)
    __svt_fatalError("svt_system: cannot initialize svt_semaphore.mux!");
}
svt_semaphore::~svt_semaphore()
{
}

void svt_semaphore::P()
{
  if (pthread_mutex_lock(&mux) < 0)
    __svt_fatalError("svt_system: cannot lock svt_semaphore.mux!");
}

void svt_semaphore::V()
{
  if (pthread_mutex_unlock(&mux) < 0)
    __svt_fatalError("svt_system: cannot unlock svt_semaphore.mux!");
}

bool svt_semaphore::tryLock()
{
  if (pthread_mutex_trylock(&mux) == 0)
    return true;
  else
    return false;

}

#else

svt_semaphore::svt_semaphore(int i)
{
  mux = CreateSemaphore(NULL, i, i, NULL);
  if (mux == NULL)
    __svt_fatalError("svt_system: cannot initialize svt_semaphore.mux!");
}
svt_semaphore::~svt_semaphore()
{
  CloseHandle(mux);
}

void svt_semaphore::P()
{
  if (WaitForSingleObject(mux, INFINITE) != WAIT_OBJECT_0)
    __svt_fatalError("svt_system: cannot lock svt_semaphore.mux!");
}

void svt_semaphore::V()
{
  if (!ReleaseSemaphore(mux, 1, NULL))
    __svt_fatalError("svt_system: cannot unlock svt_semaphore.mux!");
}

bool svt_semaphore::tryLock()
{
  return (WaitForSingleObject(mux, 0) == WAIT_OBJECT_0);
}
#endif

///////////////////////////////////////////////////////////////////////////////
// SVT_GA_IND
///////////////////////////////////////////////////////////////////////////////
/**
 * print genes to cout
 */
void svt_ga_ind::printGenes()
{
  char pOut[1024], pOutTmp[256];

  sprintf(pOut, "%8.6f %2d - ", getFitness(), getOrigin());

  for (unsigned int i = 0; i < m_oGenes.size(); i++) {
    sprintf(pOutTmp, "%5.3f ", m_oGenes[i]);
    strcat(pOut, pOutTmp);
  }
  cout << pOut << endl;
}

/**
 * print genes to cout
 */
void svt_ga_ind::printGenes(FILE *file)
{
  if (file != NULL) {
    fprintf(file, "%8.6f - ", getFitness());

    for (unsigned int i = 0; i < m_oGenes.size(); i++)
      fprintf(file, "%8.6f ", m_oGenes[i]);

    fprintf(file, "\n");
  }
}

/**
 * print genes to cout
 */
void svt_ga_ind::printGenesPf()
{
  printf("%8.6f - ", getFitness());

  for (unsigned int i = 0; i < m_oGenes.size(); i++)
    printf("%8.6f ", m_oGenes[i]);

  printf("\n");
}

/**
 * calculate the distance between two individuals
 * \param rOther reference to the other individual
 * \return vector distance between the two gene-vectors
 */
double svt_ga_ind::distance(svt_ga_ind &rOther)
{
  double fDist = 0.0;

  for (unsigned int i = 0; i < m_oGenes.size(); i++)
    fDist += (m_oGenes[i] - rOther.m_oGenes[i]) * (m_oGenes[i] - rOther.m_oGenes[i]);

  return sqrt(fDist) / (double)m_oGenes.size() ;
}



///////////////////////////////////////////////////////////////////////////////
// SVT_GACYLINDER_IND
///////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////
// The Euler angle class
///////////////////////////////////////////////////////////////////////////////

/**
 * Constructor
 */
svt_eulerAngles::svt_eulerAngles()
{
  m_pAngles = NULL;
  m_iAngles = 0;
};

/**
 * Destructor
 */
svt_eulerAngles::~svt_eulerAngles()
{
};

/**
 * Fill the angle table.
 * \param fPsiFrom   lower boundary of the psi angles (in degrees)
 * \param fPsiTo     upper boundary of the psi angles (in degrees)
 * \param fThetaFrom lower boundary of the theta angles (in degrees)
 * \param fThetaTo   upper boundary of the theta angles (in degrees)
 * \param fPhiFrom   lower boundary of the phi angles (in degrees)
 * \param fPhiTo     upper boundary of the phi angles (in degrees)
 * \param fDelta     angular step size in degrees
 */
void svt_eulerAngles::initTable(double fPsiFrom, double fPsiTo, double fThetaFrom, double fThetaTo, double fPhiFrom, double fPhiTo, double fDelta)
{
  if (m_iAngles != 0) {
    free(m_pAngles);
    m_pAngles = NULL;
  }
  //m_iAngles = proportionalEulerAngles( fPsiFrom, fPsiTo, fThetaFrom, fThetaTo, fPhiFrom, fPhiTo, fDelta );
  m_iAngles = proportionalEulerAnglesThetaOrdering(fPsiFrom, fPsiTo, fThetaFrom, fThetaTo, fPhiFrom, fPhiTo, fDelta);
};

/**
 * remove angles from table
 * WARNING: not fully tested
 * \param fPsiFrom   lower boundary of the psi angles (in degrees)
 * \param fPsiTo     upper boundary of the psi angles (in degrees)
 * \param fThetaFrom lower boundary of the theta angles (in degrees)
 * \param fThetaTo   upper boundary of the theta angles (in degrees)
 * \param fPhiFrom   lower boundary of the phi angles (in degrees)
 * \param fPhiTo     upper boundary of the phi angles (in degrees)
 */
void svt_eulerAngles::removeAngles(double fPsiFrom, double fPsiTo, double fThetaFrom, double fThetaTo, double fPhiFrom, double fPhiTo)
{
  if (m_pAngles == NULL) {
    SVTLBO << "ERROR: the table was not yet initialized!" << endl;
  }

  int iNewCount = 0;
  Real64 fPhi, fTheta, fPsi;

  //count how many angles will remain in the table
  for (unsigned iIndex = 0; iIndex < m_iAngles; iIndex++) {
    fPhi = rad2deg(getPhi(iIndex));
    fTheta = rad2deg(getTheta(iIndex));
    fPsi = rad2deg(getPsi(iIndex));
    if ((fPhi <= fPhiFrom || fPhi > fPhiTo) && (fTheta <= fThetaFrom || fTheta > fThetaTo) && (fPsi <= fPsiFrom || fPsi > fPsiTo))
      iNewCount++;
  }

  // allocate memory for the new angles

  float *pAngles = (float *) malloc(iNewCount * 3 * sizeof(float));
  if (pAngles == NULL)
    SVTLBO << "Error: svt_removeAngles, error in memory allocation..." << endl;

  // add am
  unsigned int j = 0;
  for (unsigned iIndex = 0; iIndex < m_iAngles; iIndex++) {
    fPhi = rad2deg(getPhi(iIndex));
    fTheta = rad2deg(getTheta(iIndex));
    fPsi = rad2deg(getPsi(iIndex));

    if ((fPhi <= fPhiFrom || fPhi > fPhiTo) && (fTheta <= fThetaFrom || fTheta > fThetaTo) && (fPsi <= fPsiFrom || fPsi > fPsiTo)) {
      *(pAngles + j * 3 + 0) = (float)deg2rad(fPsi);
      *(pAngles + j * 3 + 1) = (float)deg2rad(fTheta);
      *(pAngles + j * 3 + 2) = (float)deg2rad(fPhi);

      j++;
    }
  }

  //delete all angles
  if (m_iAngles != 0) {
    free(m_pAngles);
    m_pAngles = NULL;
  }

  //set the new table
  m_iAngles = iNewCount;
  m_pAngles = pAngles;



};

/**
 * add the oppsite angles: if angle = (psi, theta, phi)
 * add (-psi, theta, phi) (psi, theta, -phi) (-psi, theta, -phi)
 * eq with (2pi-psi, theta, phi) (psi, theta, 2pi-phi) (2pi-psi, theta, 2pi-phi)
 */
void svt_eulerAngles::addOppositeAngles()
{
  if (m_pAngles == NULL) {
    SVTLBO << "ERROR: the table was not yet initialized!" << endl;
  }


  // allocate memory for the new angles
  int iNewCount = 4 * m_iAngles;
  float *pAngles = (float *) malloc(iNewCount * 3 * sizeof(float));

  if (pAngles == NULL)
    SVTLBO << "Error: svt_removeAngles, error in memory allocation..." << endl;


  Real64 fPhi, fTheta, fPsi;
  unsigned int j = 0;

  for (unsigned iIndex = 0; iIndex < m_iAngles; iIndex++) {
    fPhi = rad2deg(getPhi(iIndex));
    fTheta = rad2deg(getTheta(iIndex));
    fPsi = rad2deg(getPsi(iIndex));

    *(pAngles + j * 3 + 0) = (float)deg2rad(fPsi);
    *(pAngles + j * 3 + 1) = (float)deg2rad(fTheta);
    *(pAngles + j * 3 + 2) = (float)deg2rad(fPhi);
    j++;

    *(pAngles + j * 3 + 0) = (float)deg2rad(2 * PI - fPsi);
    *(pAngles + j * 3 + 1) = (float)deg2rad(fTheta);
    *(pAngles + j * 3 + 2) = (float)deg2rad(fPhi);
    j++;

    *(pAngles + j * 3 + 0) = (float)deg2rad(fPsi);
    *(pAngles + j * 3 + 1) = (float)deg2rad(fTheta);
    *(pAngles + j * 3 + 2) = (float)deg2rad(2 * PI - fPhi);
    j++;

    *(pAngles + j * 3 + 0) = (float)deg2rad(2 * PI - fPsi);
    *(pAngles + j * 3 + 1) = (float)deg2rad(fTheta);
    *(pAngles + j * 3 + 2) = (float)deg2rad(2 * PI - fPhi);
    j++;

  }

  m_iAngles = iNewCount;
  m_pAngles = pAngles;

};


/**
 * How many angles do we have stored in total?
 * \return unsigned long with the number of angles
 */
unsigned long svt_eulerAngles::getAngleCount()
{
  return m_iAngles;
};

/**
 * The angles follow the common PTP (Goldstein) convention. This function returns the Psi angle.
 * \param iIndex index into the table of angles
 * \return psi angle
 */
float svt_eulerAngles::getPsi(unsigned long iIndex)
{
  if (iIndex >= m_iAngles) {
    SVTLBO << "ERROR: iIndex > m_iAngles: " << iIndex << " > " << m_iAngles << endl;
    return *(m_pAngles + (m_iAngles - 1) * 3 + 0);
  } else
    return *(m_pAngles + iIndex * 3 + 0);
};

/**
 * The angles follow the common PTP (Goldstein) convention. This function returns the Theta angle.
 * \param iIndex index into the table of angles
 * \return theta angle
 */
float svt_eulerAngles::getTheta(unsigned long iIndex)
{
  if (iIndex >= m_iAngles) {
    SVTLBO << "ERROR: iIndex > m_iAngles: " << iIndex << " > " << m_iAngles << endl;
    return *(m_pAngles + (m_iAngles - 1) * 3 + 1);
  } else
    return *(m_pAngles + iIndex * 3 + 1);
};

/**
 * The angles follow the common PTP (Goldstein) convention. This function returns the Phi angle.
 * \param iIndex index into the table of angles
 * \return phi angle
 */
float svt_eulerAngles::getPhi(unsigned long iIndex)
{
  if (iIndex >= m_iAngles) {
    SVTLBO << "ERROR: iIndex > m_iAngles: " << iIndex << " > " << m_iAngles << endl;
    return *(m_pAngles + (m_iAngles - 1) * 3 + 2);
  } else
    return *(m_pAngles + iIndex * 3 + 2);
};

/**
 * This function precomputes the angle table. It is called automatically in the constructor.
 */
unsigned long svt_eulerAngles::proportionalEulerAngles(double fPsiFrom, double fPsiTo, double fThetaFrom, double fThetaTo, double fPhiFrom, double fPhiTo, double fDelta)
{
  double psi, theta, phi;
  double psi_ang_dist, psi_real_dist;
  double theta_real_dist, phi_real_dist;
  double psi_steps, theta_steps, phi_steps;
  double psi_range, theta_range, phi_range;
  unsigned long u, j;

  unsigned long iCount;

  if ((fPsiTo   - fPsiFrom) / fDelta < -1  ||
      (fThetaTo - fThetaFrom) / fDelta < -1  ||
      (fPhiTo   - fPhiFrom) / fDelta < -1) {
    SVTLBO << "Error: svt_proportionalEulerAngles(), ranges wrong..." << endl;
    return 0;
  }

  psi_range   = fPsiTo   - fPsiFrom;
  theta_range = fThetaTo - fThetaFrom;
  phi_range   = fPhiTo   - fPhiFrom;

  // Use rounding instead of CEIL to avoid rounding up at x.001
  phi_steps       = rint((phi_range / fDelta) + 0.499);
  phi_real_dist   = phi_range / phi_steps;

  theta_steps     = rint((theta_range / fDelta) + 0.499);
  theta_real_dist = theta_range / theta_steps;

  // Computes the number of angles that will be generated
  u = 0;
  for (phi = fPhiFrom; phi < 360.0 && phi <= fPhiTo;  phi += phi_real_dist) {
    for (theta = fThetaFrom; theta <= 180.0 && theta <= fThetaTo;  theta += theta_real_dist) {
      if (theta == 0.0 || theta == 180.0)
        psi_steps = 1;
      else
        psi_steps = rint(360.0 * cos(deg2rad(90.0 - theta)) / fDelta);

      psi_ang_dist  = 360.0 / psi_steps;
      psi_real_dist = psi_range / (ceil(psi_range / psi_ang_dist));

      for (psi = fPsiFrom; psi < 360.0 && psi <= fPsiTo;  psi += psi_real_dist)
        u++;
    }
  }

  iCount = u;

  // allocate memory
  m_pAngles = (float *) malloc(iCount * 3 * sizeof(float));
  if (m_pAngles == NULL)
    SVTLBO << "Error: svt_proportionalEulerAngles(), error in memory allocation..." << endl;

  j = 0;
  for (phi = fPhiFrom; phi < 360.0 && phi <= fPhiTo;  phi += phi_real_dist) {
    for (theta = fThetaFrom; theta <= 180.0 && theta <= fThetaTo;  theta += theta_real_dist) {
      if (theta == 0.0 || theta == 180.0)
        psi_steps = 1;
      else
        psi_steps = rint(360.0 * cos(deg2rad(90.0 - theta)) / fDelta);

      psi_ang_dist  = 360.0 / psi_steps;
      psi_real_dist = psi_range / (ceil(psi_range / psi_ang_dist));

      for (psi = fPsiFrom; psi < 360.0 && psi <= fPsiTo;  psi += psi_real_dist) {

        *(m_pAngles + j * 3 + 0) = (float)deg2rad(psi);
        *(m_pAngles + j * 3 + 1) = (float)deg2rad(theta);
        *(m_pAngles + j * 3 + 2) = (float)deg2rad(phi);

        //SVTLBO << " [" << phi << " , " << theta << " , " << psi << "]" <<  endl;
        j++;
      }
    }
  }

  return iCount;
};


/**
 * This function precomputes the angle table. It is called automatically in the constructor.
 */
unsigned long svt_eulerAngles::proportionalEulerAnglesThetaOrdering(double fPsiFrom, double fPsiTo, double fThetaFrom, double fThetaTo, double fPhiFrom, double fPhiTo, double fDelta)
{
  double psi, theta, phi;
  double psi_ang_dist, psi_real_dist;
  double theta_real_dist, phi_real_dist;
  double psi_steps, theta_steps, phi_steps;
  double psi_range, theta_range, phi_range;
  unsigned long u, j;

  unsigned long iCount;

  if ((fPsiTo   - fPsiFrom) / fDelta < -1  ||
      (fThetaTo - fThetaFrom) / fDelta < -1  ||
      (fPhiTo   - fPhiFrom) / fDelta < -1) {
    SVTLBO << "Error: svt_proportionalEulerAngles(), ranges wrong..." << endl;
    return 0;
  }

  psi_range   = fPsiTo   - fPsiFrom;
  theta_range = fThetaTo - fThetaFrom;
  phi_range   = fPhiTo   - fPhiFrom;

  // Use rounding instead of CEIL to avoid rounding up at x.001
  phi_steps       = rint((phi_range / fDelta) + 0.499);
  phi_real_dist   = phi_range / phi_steps;

  theta_steps     = rint((theta_range / fDelta) + 0.499);
  theta_real_dist = theta_range / theta_steps;

  // Computes the number of angles that will be generated
  u = 0;
  for (theta = fThetaFrom; theta <= 180.0 && theta <= fThetaTo;  theta += theta_real_dist) {
    if (theta == 0.0 || theta == 180.0)
      psi_steps = 1;
    else
      psi_steps = rint(360.0 * cos(deg2rad(90.0 - theta)) / fDelta);

    for (phi = fPhiFrom; phi < 360.0 && phi <= fPhiTo;  phi += phi_real_dist) {
      psi_ang_dist  = 360.0 / psi_steps;
      psi_real_dist = psi_range / (ceil(psi_range / psi_ang_dist));

      for (psi = fPsiFrom; psi < 360.0 && psi <= fPsiTo;  psi += psi_real_dist)
        u++;
    }
  }

  iCount = u;

  // allocate memory
  m_pAngles = (float *) malloc(iCount * 3 * sizeof(float));
  if (m_pAngles == NULL)
    SVTLBO << "Error: svt_proportionalEulerAngles(), error in memory allocation..." << endl;

  j = 0;
  for (theta = fThetaFrom; theta <= 180.0 && theta <= fThetaTo;  theta += theta_real_dist) {
    if (theta == 0.0 || theta == 180.0)
      psi_steps = 1;
    else
      psi_steps = rint(360.0 * cos(deg2rad(90.0 - theta)) / fDelta);

    for (phi = fPhiFrom; phi < 360.0 && phi <= fPhiTo;  phi += phi_real_dist) {

      psi_ang_dist  = 360.0 / psi_steps;
      psi_real_dist = psi_range / (ceil(psi_range / psi_ang_dist));

      for (psi = fPsiFrom; psi < 360.0 && psi <= fPsiTo;  psi += psi_real_dist) {

        *(m_pAngles + j * 3 + 0) = (float)deg2rad(psi);
        *(m_pAngles + j * 3 + 1) = (float)deg2rad(theta);
        *(m_pAngles + j * 3 + 2) = (float)deg2rad(phi);

        //SVTLBO << " [" << phi << " , " << theta << " , " << psi << "]" <<  endl;
        j++;
      }
    }
  }

  return iCount;
};

/**
 * searches the angles that are within fAngleRange from the angle indicated by iIndex
 * \param iIndex the reference angle around which to search
 * \param fRange how far away from the reference
 * \return a list of indexes that indicates the angles in the angle list that are close to the angle idicated by iIndex
 */
vector<long unsigned int> svt_eulerAngles::getNeighborAngles(unsigned long int iIndex, Real64 fAngleRange)
{

  vector<long unsigned int> oVec;
  if (iIndex >= m_iAngles) {
    SVTLBO << "ERROR: iIndex > m_iAngles: " << iIndex << " > " << m_iAngles << endl;
    return oVec;
  }

  Real64 fPhiRef  = getPhi(iIndex);
  Real64 fThetaRef  = getTheta(iIndex);
  Real64 fPsiRef  = getPsi(iIndex);
  //SVTLBO << "DEBU : " << fPhiRef << " , " <<  fThetaRef << " , " << fPsiRef <<  endl;

  Real64 fPhi, fTheta, fPsi;
  for (long unsigned int iAngle = 0; iAngle < m_iAngles; iAngle++) {
    fPhi    = getPhi(iAngle);
    fTheta    = getTheta(iAngle);
    fPsi    = getPsi(iAngle);

    if (abs(fPhi - fPhiRef) < fAngleRange && abs(fTheta - fThetaRef) < fAngleRange && abs(fPsi - fPsiRef) < fAngleRange) {
      //SVTLBO << "DEBUG: " << fPhi << " , " <<  fTheta << " , " << fPsi << endl;
      oVec.push_back(iAngle);
    }

  }

  return oVec;
};

///////////////////////////////////////////////////////////////////////////////
// svt_gacylinder_ind class
///////////////////////////////////////////////////////////////////////////////

/**
 * Constructor
 */
svt_gacylinder_ind::svt_gacylinder_ind(): svt_ga_ind(),
  m_bWrote(false),
  m_iTurns(10),
  m_fHeightTurn(5.1 / 4.0)
{
  m_oTrans.loadIdentity();
}
/**
     * destructor
     */
svt_gacylinder_ind::~svt_gacylinder_ind()
{

};

/**
 * create the coarse phenotype (equivalent simple pdb) but don't fill yet atomic coordinates
 * \param number of units
 */
void svt_gacylinder_ind::buildCoarsePhenotype()
{
  m_oCoarsePhenotype.clear();
  svt_vector4<Real64> oVec;

  //oVec.x(0.0f); oVec.y(0.0f); oVec.z( -(Real64)m_iTurns*m_fHeightTurn);
  oVec.x(0.0f);
  oVec.y(0.0f);
  oVec.z(0.0);
  m_oCoarsePhenotype.push_back(oVec);

  //oVec.x(0.0f); oVec.y(0.0f); oVec.z( (Real64)m_iTurns*m_fHeightTurn);
  oVec.x(0.0f);
  oVec.y(0.0f);
  oVec.z(1.0);
  m_oCoarsePhenotype.push_back(oVec);
};

/**
 * update the coarse phenotype for the given unit
 * \param transformation matrix for that unit
 */
void svt_gacylinder_ind::updateCoarsePhenotype(svt_ga_mat oMat)
{
  buildCoarsePhenotype();

  for (unsigned int iPoint = 0; iPoint < m_oCoarsePhenotype.size(); iPoint++)
    m_oCoarsePhenotype[iPoint] = oMat * m_oCoarsePhenotype[iPoint]  ;
  setTrans(oMat);
}


/**
 * get the coarse phenotype
 */
svt_point_cloud_pdb< svt_vector4<Real64 > > svt_gacylinder_ind::getCoarsePhenotype()
{
  svt_point_cloud_pdb<svt_vector4<Real64> > oCoarsePhenotype;
  svt_point_cloud_atom oAtom;

  for (unsigned int iIndex = 0; iIndex < (Real64)m_oCoarsePhenotype.size(); iIndex++)
    oCoarsePhenotype.addAtom(oAtom, m_oCoarsePhenotype[iIndex]);

  return oCoarsePhenotype;
};

/**
 * calculate the distance between two individuals
 * \param rOther reference to the other individual
 * \return vector distance between the two gene-vectors
 */
Real64 svt_gacylinder_ind::distance(svt_gacylinder_ind &rOther)
{
  Real64 fDist = 0.0;
  if (m_oCoarsePhenotype.size() != rOther.m_oCoarsePhenotype.size() || m_oCoarsePhenotype.size() == 0) {
    SVTLBO << "Can not compute distance" << m_oCoarsePhenotype.size() << " " << rOther.m_oCoarsePhenotype.size() << " " <<  m_oCoarsePhenotype.size() <<  endl;
    return 0.0;
  }
  /*
      svt_vector4< Real64 > oDiffThis  = m_oCoarsePhenotype[1]-m_oCoarsePhenotype[0];
      svt_vector4< Real64 > oDiffOther = rOther.m_oCoarsePhenotype[1]-rOther.m_oCoarsePhenotype[0];
      svt_vector4< Real64 > oDiff = m_oCoarsePhenotype[0] - rOther.m_oCoarsePhenotype[0];

      svt_vector4< Real64 > oCrossProd;
      oCrossProd.x( oDiffThis.y()*oDiffOther.z() - oDiffThis.z()*oDiffOther.y() );
      oCrossProd.y( oDiffThis.z()*oDiffOther.x() - oDiffThis.x()*oDiffOther.z() );
      oCrossProd.z( oDiffThis.x()*oDiffOther.y() - oDiffThis.y()*oDiffOther.x() );

      fDist = (oDiff.x()*oCrossProd.x()+oDiff.y()*oCrossProd.y()+oDiff.z()*oCrossProd.z() )/oCrossProd.length();

      return fDist;
  */


  for (unsigned int i = 0; i < 1/*m_oCoarsePhenotype.size()*/; i++)
    fDist += m_oCoarsePhenotype[i].distanceSq(rOther.m_oCoarsePhenotype[i]);

  //return sqrt( fDist/(Real64)((Real64)m_oCoarsePhenotype.size()/4.0f) );
  return sqrt(fDist);

};

/**
 * set Wrote on disk
 * \param bWrote whether it was already wrote on disk
 */
void svt_gacylinder_ind::setWrote(bool bWrote)
{
  m_bWrote = bWrote;
};

/**
 * get Wrote on disk
 * \return bWrote whether it was already wrote on disk
 */
bool svt_gacylinder_ind::getWrote()
{
  return m_bWrote;
};

/**
 * get the number of turns
 * \return the number of turns
 */
unsigned int svt_gacylinder_ind::getTurns()
{
  return m_iTurns;
};

/**
* get the number of turns
* \param the number of turns
*/
void svt_gacylinder_ind::setTurns(unsigned int iTurns)
{
  m_iTurns = iTurns;
};


/**
 * get the height of a turn
 * \return the height
 */
Real64 svt_gacylinder_ind::getHeightTurn()
{
  return m_fHeightTurn;
};

/**
* set the height of a turn
* \param the height
*/
void svt_gacylinder_ind::setHeightTurn(Real64 fHeightTurn)
{
  m_fHeightTurn = fHeightTurn;
};

/**
 * get the Transformation
 * \return the transformation
 */
svt_ga_mat svt_gacylinder_ind::getTrans()
{
  return m_oTrans;
};

/**
 * set the Transformation
 * \param the transformation
 */
void svt_gacylinder_ind::setTrans(svt_ga_mat &rTrans)
{
  m_oTrans = rTrans;
};


/**
 * set Fitness Top
 * \param the new fitness
 */
void svt_gacylinder_ind::setFitnessTop(Real64 fFitness)
{
  m_fFitnessTop = fFitness;
};

/**
 * get Fitness Top
 * \return the new fitness
 */
Real64 svt_gacylinder_ind::getFitnessTop()
{
  return m_fFitnessTop;
};

/**
 * set Fitness Bot
 * \param the new fitness
 */
void svt_gacylinder_ind::setFitnessBot(Real64 fFitness)
{
  m_fFitnessBot = fFitness;
};

/**
 * get Fitness Top
 * \return the new fitness
 */
Real64 svt_gacylinder_ind::getFitnessBot()
{
  return m_fFitnessBot;
};

///////////////////////////////////////////////////////////////////////////////
// SVT_TUBE
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
// svt_tube class
///////////////////////////////////////////////////////////////////////////////

/**
 * Constructor
 */
svt_tube::svt_tube():
  m_fAvgFitness(0.0),
  m_fSumFitness(0.0),
  m_iSize(0),
  m_fLength(0),
  m_fMapScore(0.0),
  m_bWasTubeComputed(false),
  m_bWasFineTubeComputed(false),
  m_fPenalty(0.0),
  m_bRecomputeScore(true)
{
};

/**
 * add a new element to tube
 * \param oElem a svt_gacylinder_ind that in considered in this tube
 */
void svt_tube::add2Tube(svt_gacylinder_ind oElem, bool bFlip)
{
  if (m_iSize == 0)
    m_oFirstAddedElem =  oElem.getTrans().translation();

  if (bFlip) {
    vector<svt_gacylinder_ind> oElements;
    while (m_oElements.size() > 0) {
      oElements.push_back(m_oElements[m_oElements.size() - 1]);
      m_oElements.pop_back();
    }
    m_oElements = oElements;
  }

  m_oElements.push_back(oElem);
  m_iSize++;
  m_bWasTubeComputed = false;
  m_bWasFineTubeComputed = false;

  m_fSumFitness += oElem.getFitness();
  if (m_iSize >= 2)
    m_fAvgFitness = (m_fSumFitness / (Real64)m_iSize);// * sqrt(m_iSize);
  else
    m_fAvgFitness = 0.0f; // do not consider yet

};

/**
 * remove the last element of the tube
 */
void svt_tube::pop()
{
  if (m_oElements.size() > 0) {
    m_fSumFitness -= m_oElements[m_oElements.size() - 1].getFitness();
    m_oElements.pop_back();
    m_iSize--;

    if (m_iSize >= 2)
      m_fAvgFitness = m_fSumFitness / (Real64)m_iSize;
    else
      m_fAvgFitness = 0.0f; // do not consider yet
  }

};

/**
 * get the list of individuals added to the tube
 * \return the individuals
 */
vector<svt_gacylinder_ind > svt_tube::getElements()
{
  return m_oElements;
};

/**
 * set the score
 * \param fScore the score assigned
 */
void svt_tube::setMapScore(Real64 fScore)
{
  m_fMapScore = fScore;
};

/**
 * get the score
 * \return the score of the tube
 */
Real64 svt_tube::getMapScore()
{
  return m_fMapScore;
};

/**
 * set the score
 * \param fScore the score assigned
 */
void svt_tube::setAvgScore(Real64 fScore)
{
  m_fAvgFitness = fScore;
};

/**
 * get the score
 * \return the score of the tube
 */
Real64 svt_tube::getAvgScore()
{
  return m_fAvgFitness;
};

/**
 * set the score
 * \param iIndex - which element in vector score
 * \param fScore - the value
 */
void svt_tube::setScore(int iIndex, Real64 fScore)
{
  while (iIndex >= (int)m_oScores.size())
    m_oScores.push_back(0.0f);

  m_oScores[iIndex] = fScore;
};

/**
 * returns the list of scores
 * \param a vector with the different scores
 */
vector<Real64> svt_tube::getScores()
{
  return m_oScores;
};

/**
 * Compute avg score
 */
void svt_tube::computeAvgScore()
{
  m_fSumFitness = 0.0;
  //discard the first two elements and the last two elements as they are the steps that failed
  for (int iIndex = 0; iIndex < (int)m_oElements.size(); iIndex++) {
    m_fSumFitness += m_oElements[iIndex].getFitness();
  }

  if ((int)m_oElements.size() > 0)
    m_fAvgFitness = m_fSumFitness / (Real64)m_oElements.size();
  else
    m_fAvgFitness = 0.0;
};

/**
* get the number of points
* \return the number of points
*/
unsigned int svt_tube::size()
{
  return m_iSize;
};


/**
* get the number of turns
* \return the number of turns
*/
Real64 svt_tube::getTurns(Real64 fRatio)
{
  if (fRatio != 0)
    return m_iSize / (Real64)(fRatio);
  else {
    SVTLBO << "Could not compute the number of turns: fRatio = 0 " << endl;
    return (Real64)m_iSize;
  }
};

/**
 * Compute the length in A0
 */
void svt_tube::computeLength()
{
  //compute points on axis if not yet available
  if (!m_bWasTubeComputed)
    getTube();

  m_fLength = 0.0f;
  for (unsigned int i = 1; i < m_oTube.size(); i++)
    m_fLength += m_oTube[i - 1].distance(m_oTube[i]);
};

/**
* get the length = sum distances between points
* assumes that the points are in order = point 0 is closest 1 and point 2 follows point 1... etc
* \return the length of the tube
*/
Real64 svt_tube::getLength()
{
  if (m_fLength == 0.0 || !m_bWasTubeComputed)
    computeLength();

  return m_fLength;
};

/**
 * get the first element
 */
svt_vector4<Real64> svt_tube::getFirstAddedElem()
{
  return m_oFirstAddedElem;
};

/**
 * set penalty
 * \param the penalty
 */
void svt_tube::setPenalty(Real64 fPenalty)
{
  m_fPenalty = fPenalty;
};

/**
 * increase penalty
 * \param the penalty
 */
void svt_tube::addPenalty(Real64 fPenalty)
{
  m_fPenalty += fPenalty;
};

/**
 * get penalty
 * \return the penalty
 */
Real64 svt_tube::getPenalty()
{
  return m_fPenalty;
};

/**
 * get curvature
 */
Real64 svt_tube::getCurvature()
{
  return m_fCurvature;
}

/**
 * overload < operator
 * \param that another svt_tube element
 */
bool svt_tube::operator<(const svt_tube &that) const
{
  return m_fAvgFitness < that.m_fAvgFitness;
};

/**
 * overload < operator using the max density of the map
 * \param that another svt_tube element
 */
bool svt_tube::lt_mapScore(svt_tube first, svt_tube second)
{
  return first.m_fMapScore < second.m_fMapScore;
};

/**
 * overload < operator using the max density of the map
 * \param that another svt_tube element
 */
bool svt_tube::lt_length(svt_tube first, svt_tube second)
{
  return first.m_iSize < second.m_iSize;
};


/**
 * overload < operator using the max density of the map
 * \param that another svt_tube element
 */
bool svt_tube::lt_wPenatly(svt_tube first, svt_tube second)
{
  return first.m_fAvgFitness - first.getPenalty() < second.m_fAvgFitness - second.getPenalty();
};

/**
 * overload < operator using the max density of the map
 * \param that another svt_tube element
 */
bool  svt_tube::lt_score(svt_tube first, svt_tube second)
{
  if ((int)(first.m_oScores.size()) > 2 && (int)(second.m_oScores.size()) > 2)
    return first.m_oScores[2] < second.m_oScores[2];
  else
    return true;
};



/**
 * get tube as
 * \return the pdb tube
 */
svt_point_cloud_pdb<svt_ga_vec> svt_tube::getTube(svt_point_cloud_pdb<svt_ga_vec> *pTemplate, bool bFine)
{
  //was already computed once just return dont recompute
  if (m_bWasTubeComputed && pTemplate == NULL && !bFine)
    return m_oTube;

  // fine pdb was already computed once just return dont recompute
  if (m_bWasFineTubeComputed && pTemplate == NULL && bFine)
    return m_oFineTube;

  svt_point_cloud_pdb<svt_ga_vec> oPdb;
  svt_point_cloud_atom oAtom;
  oAtom.setName(" C");
  oAtom.setRemoteness('A');
  oAtom.setResname("GLY");
  svt_ga_vec oVec, oNull;
  svt_point_cloud_pdb< svt_ga_vec> oTemplatePdb;
  oNull.x(0.0);
  oNull.y(0.0);
  oNull.z(0.0);
  svt_ga_mat oMat;

  for (unsigned int iIndex = 0; iIndex < m_oElements.size(); iIndex++) {
    oMat = m_oElements[iIndex].getTrans();

    if (pTemplate == NULL) {
      oAtom.setResidueSeq(oPdb.size());
      oVec = oMat * oNull;
      oPdb.addAtom(oAtom, oVec);
    } else {
      oTemplatePdb = oMat * (*pTemplate);
      oPdb.append(oTemplatePdb);
    }
  }
  m_oTube = oPdb;
  m_bWasTubeComputed = true;

  if (pTemplate == NULL) {
    if (bFine  && !m_bWasFineTubeComputed) { // not yet computed
      Real64 fNewX, fNewY, fNewZ, fM, fDist;
      int iCountIntermediates, iCount = 0;
      svt_vector4<Real64> oVec1, oVec2, oVec;

      for (unsigned int iIndex = 0; iIndex < oPdb.size() - 1; iIndex++) {
        iCount++;
        oAtom.setResidueSeq(iCount);
        oAtom.setOrdResidueSeq(iCount);
        oAtom.setPDBIndex(iCount);
        m_oFineTube.addAtom(oAtom, oPdb[iIndex]);

        fDist = oPdb[iIndex].distance(oPdb[iIndex + 1]);
        if (fDist > 6.0) { // should create intermediates
          iCountIntermediates = floor(fDist / 4.0); //how many points should have beween
          for (int iPoint = 1; iPoint < iCountIntermediates; iPoint++) {
            // cout << iIndex << ":" << iPoint << endl;
            oVec1 = oPdb[iIndex];
            oVec2 = oPdb[iIndex + 1];

            fM = 1.0 / iCountIntermediates * iPoint;
            fNewX =  oVec1.x() + fM * (oVec2.x() - oVec1.x());
            fNewY =  oVec1.y() + fM * (oVec2.y() - oVec1.y());
            fNewZ =  oVec1.z() + fM * (oVec2.z() - oVec1.z());

            oVec.x(fNewX);
            oVec.y(fNewY);
            oVec.z(fNewZ);
            iCount++;
            oAtom.setResidueSeq(iCount);
            oAtom.setResidueSeq(iCount);
            oAtom.setOrdResidueSeq(iCount);
            oAtom.setPDBIndex(iCount);
            m_oFineTube.addAtom(oAtom, oVec);
          }
        }
      }
      iCount++;
      oAtom.setResidueSeq(iCount);
      oAtom.setOrdResidueSeq(iCount);
      oAtom.setPDBIndex(iCount);
      m_oFineTube.addAtom(oAtom, oPdb[oPdb.size() - 1]);
      m_bWasFineTubeComputed = true;
    }
  }

  //compute the lenght in A
  computeLength();

  if (!bFine)
    return m_oTube;
  else
    return m_oFineTube;
};

/**
 * get the direction of the tube
 */
svt_ga_vec svt_tube::getDirection()
{
  m_oDirection.x(0.0f);
  m_oDirection.y(0.0f);
  m_oDirection.z(0.0f);

  svt_ga_vec oNull, oPrev, oVec;
  oNull.x(0.0);
  oNull.y(0.0);
  oNull.z(0.0);
  svt_ga_mat oMat;

  if (m_oElements.size() > 0) {
    oMat = m_oElements[0].getTrans();
    oPrev = oMat * oNull;
  }

  unsigned int iCount = 0;
  for (unsigned int iIndex = 1; iIndex < m_oElements.size(); iIndex++) {
    iCount++;
    oMat = m_oElements[iIndex].getTrans();
    oVec = oMat * oNull;

    m_oDirection += (oVec - oPrev);

    oPrev = oVec;
  }

  if (iCount != 0) {
    m_oDirection *= (1.0 / (Real64(iCount)));
    m_oDirection.normalize();
  }

  return m_oDirection;
};

/**
 * print the tube
 */
void svt_tube::print()
{
  SVTLBO << "Tube of " << m_iSize << "  with avg fitness of " << m_fAvgFitness << endl;

  for (unsigned int iIndex = 0; iIndex < m_oElements.size(); iIndex++) {
    SVTLBO;
    m_oElements[iIndex].printGenes();
  }
};

/**
 * compute the volume underneath the tube
 */
void svt_tube::fillExplored(svt_ga_vol &rVol, svt_ga_vol *pVol)
{
  svt_point_cloud_pdb< svt_vector4<Real64> > oPdb;
  oPdb = getTube();

  /*
      svt_ga_vol oExtractedVol;
      oExtractedVol = rVol;
      oExtractedVol.setValue(0.0);

      svt_ga_mat oMat;

      oPdb.projectMass( &oExtractedVol,oMat, false );
  */
  svt_ga_vec oVecVol;
  Real64 fWidth = rVol.getWidth();
  Real64 fGridX = rVol.getGridX();
  Real64 fGridY = rVol.getGridY();
  Real64 fGridZ = rVol.getGridZ();
  Real64 fDist, fDensity;

  for (unsigned int iX = 0; iX < rVol.getSizeX() - 1; iX++) {
    oVecVol.x(fGridX + iX * fWidth + fWidth / 2.0);
    for (unsigned int iY = 0; iY < rVol.getSizeY() - 1; iY++) {
      oVecVol.y(fGridY + iY * fWidth  + fWidth / 2.0);
      for (unsigned int iZ = 0; iZ < rVol.getSizeZ() - 1; iZ++) {
        oVecVol.z(fGridZ + iZ * fWidth  + fWidth / 2.0);
        for (unsigned int iIndex = 0; iIndex < oPdb.size(); iIndex++) {
          fDist = oPdb[iIndex].distance(oVecVol);
          if (fDist <= 6.0 / 2.0) {
            fDensity = rVol.at(iX, iY, iZ) ;
            pVol->setAt(iX, iY, iZ, fDensity);
            //oExtractedVol.setAt(iX, iY, iZ, fDensity );
          }
        }
      }
    }
  }

  //create blurring kernel
  //    svt_ga_vol oKernel;
  //    oKernel.create1DBlurringKernel(1.0, 5.0);
  //   oExtractedVol.convolve1D3D(oKernel, false);

  //get the matrix
  //  m_fMapScore =  oExtractedVol.correlation(rVol);

};

/**
 *get the high resolution version of the tube
 */
svt_point_cloud_pdb<svt_ga_vec> svt_tube::getHRTube()
{
  return m_oHRTube;
};

/**
 * creates a highresolution version of the helix
 */
void svt_tube::createHighResTube(svt_ga_vol &rVol, svt_point_cloud_pdb<svt_ga_vec> *pTemplate)
{
  svt_point_cloud_pdb< svt_ga_vec > oHelix, oBestHelix;
  Real64 fBestCorr, fCorr;

  svt_ga_mat oMat, oConstRot;
  oConstRot.loadIdentity();

  int iStep = 5;
  Real64 fAngleStep = 1.5, fAngleMin, fAngleMax, fAngleBest;

  vector<svt_point_cloud_pdb <svt_ga_vec> > oPdbs;
  for (int iIndex = m_oElements.size() - 1; iIndex >= 0; iIndex -= iStep) {
    // refine rotation
    if (iIndex == (int)m_oElements.size() - 1) {
      fAngleMin = -180;
      fAngleMax = 180;
    } else {
      fAngleMin = -45;
      fAngleMax = 45;
    }

    oMat.loadIdentity();
    oMat.rotate(2, deg2rad(fAngleMin));

    fBestCorr = pTemplate->projectMassCorr(&rVol, m_oElements[iIndex].getTrans() * oConstRot * oMat, false);
    fAngleBest = fAngleMin;
    for (Real64 fAngle = fAngleMin + fAngleStep; fAngle < fAngleMax; fAngle += fAngleStep) {
      oMat.loadIdentity();
      oMat.rotate(2, deg2rad(fAngle));

      fCorr = pTemplate->projectMassCorr(&rVol, m_oElements[iIndex].getTrans() * oConstRot * oMat, false);
      if (fCorr > fBestCorr) {
        fBestCorr = fCorr;
        fAngleBest = fAngle;
      }
    }

    oMat.loadIdentity();
    oMat.rotate(2, deg2rad(fAngleBest));

    oPdbs.push_back(m_oElements[iIndex].getTrans()*oConstRot * oMat * (*pTemplate));
    oConstRot = oMat * oConstRot;

    oConstRot.rotate(2, deg2rad(130));
  }

  for (unsigned int iIndex = 0; iIndex < oPdbs.size(); iIndex ++) {
    oHelix = oPdbs[iIndex];
    if (iIndex == 0)
      m_oHRTube =  oHelix;
    else
      m_oHRTube.append(oHelix);
  }
};

/**
 * Estimate curvature of the axis
 */
void svt_tube::estimate_curvature()
{
  svt_point_cloud_pdb<svt_ga_vec> oPdb = getTube();
  svt_ga_vec oRef, oDiff;

  m_fCurvature = 0.0;
  if (oPdb.size() >= 3) {
    oRef = oPdb[1] - oPdb[0];

    for (unsigned int iIndex = 2; iIndex < oPdb.size(); iIndex++) {
      oDiff = oPdb[iIndex] - oPdb[iIndex - 1];
      m_fCurvature += oRef.distance(oDiff);
      oRef = oDiff;
    }

    m_fCurvature /= (Real64)(oPdb.size() - 2);
  }
};

/**
 * discard the points that are at the ends of the tube is their score is < fScore
 */
void svt_tube::discardPointsAtEnds(Real64 fScore)
{
  print();

  //at the end
  while (m_oElements.size() > 0 && m_oElements[m_oElements.size() - 1].getFitness() < fScore)
    m_oElements.pop_back();

  //at the begining
  while (m_oElements.size() > 0 && m_oElements[0].getFitness() < fScore)
    m_oElements.erase(m_oElements.begin());

  //recompute Scores
  m_iSize = m_oElements.size();
  computeAvgScore();

  print();
};

///////////////////////////////////////////////////////////////////////////////
// SVT_CREATETHREAD
///////////////////////////////////////////////////////////////////////////////
#ifdef WIN32
/**
 * Create a seperate thread with the function func
 * \param pFunc pointer to the function which should be executed in the seperate thread
 * \param pArg  pointer to the arguments for the function
 * \param iPriority priority of the thread (e.g. SVT_THREAD_PRIORITY_HIGH)
 */
void svt_createThread(void *(*func)(void *), void *arg, int iPriority)
{
  unsigned long nThreadID;

  HANDLE handle = CreateThread(0, 0, (LPTHREAD_START_ROUTINE)func, (LPVOID)arg, 0, &nThreadID);

  if (handle)
    switch (iPriority) {
      case SVT_THREAD_PRIORITY_NORMAL:
        break;
      case SVT_THREAD_PRIORITY_HIGH:
        //SetPriorityClass(handle, HIGH_PRIORITY_CLASS);
        SetThreadPriority(handle, THREAD_PRIORITY_ABOVE_NORMAL);
        break;
      case SVT_THREAD_PRIORITY_LOW:
        SetThreadPriority(handle, THREAD_PRIORITY_BELOW_NORMAL);
        break;
      default:
        break;
    }
};

/**
 * Terminate thread. Has to be called by the thread itself!
 */
void svt_terminateThread()
{

};

#else

#include <pthread.h>
#include <unistd.h>

pthread_t recalc_it;

/**
 * Create a seperate thread with the function func
 * \param pFunc pointer to the function which should be executed in the seperate thread
 * \param pArg  pointer to the arguments for the function
 * \param iPriority priority of the thread (e.g. SVT_THREAD_PRIORITY_HIGH)
 */
void svt_createThread(void *(*func)(void *), void *arg, int)
{
  pthread_create(&recalc_it, NULL, func, arg);
  pthread_detach(recalc_it);
};

/**
 * Terminate thread. Has to be called by the thread itself!
 */
void svt_terminateThread()
{
  pthread_exit(NULL);
};

#endif


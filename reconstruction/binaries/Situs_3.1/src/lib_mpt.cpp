/*********************************************************************
*                           L I B _ M P T                            *
**********************************************************************
* Library is part of the Situs package URL: situs.biomachina.org     *
* (c) Stefan Birmanns and Willy Wriggers, 2009-2012                  *
**********************************************************************
*                                                                    *
* Support routines for matchpoint (matchpt) tool.                    *
*                                                                    *
**********************************************************************
* See legal statement for terms of distribution                      *
*********************************************************************/

///////////////////////////////////////////////////////////////////////////////
// INCLUDES AND CONSTANTS
///////////////////////////////////////////////////////////////////////////////

#include "situs.h"
#include "lib_mpt.h"
#include "lib_err.h"
#include "lib_vio.h"
#include "lib_pio.h"
#include <math.h>
#include <iostream>
#include <string.h>
#include <limits>


bool isPositive(const double &value)
{
  return value > std::numeric_limits<double>::epsilon();
};

///////////////////////////////////////////////////////////////////////////////
// VEC4
///////////////////////////////////////////////////////////////////////////////

/**
 * Constructor
 * \param fX initial x coordinate
 * \param fY initial y coordinate
 * \param fZ initial z coordinate
 * \param fW initial w coordinate
 */
vec4::vec4(double fX, double fY, double fZ, double fW)
{
  x(fX);
  y(fY);
  z(fZ);
  w(fW);
}

vec4::vec4(double fValue, double fW)
{
  x(fValue);
  y(fValue);
  z(fValue);
  w(fW);
}

vec4::~vec4()
{
};

vec4 &vec4::operator=(const vec4 &that)
{
  x(that.x());
  y(that.y());
  z(that.z());
  w(that.w());
  return *this;
};

double &vec4::operator[](unsigned i)
{
  return m_aData[i];
};

const double &vec4::operator[](unsigned i) const
{
  return m_aData[i];
};

//
// arithmetic operators
//
vec4 &vec4::operator+=(const vec4 &p)
{
  (*this)[0] += p[0];
  (*this)[1] += p[1];
  (*this)[2] += p[2];
  return *this;
}

vec4 &vec4::operator+=(const double &f)
{
  (*this)[0] += f;
  (*this)[1] += f;
  (*this)[2] += f;
  return *this;
}

vec4 &vec4::operator-=(const vec4 &p)
{
  (*this)[0] -= p[0];
  (*this)[1] -= p[1];
  (*this)[2] -= p[2];
  return *this;
}

vec4 &vec4::operator-=(const double &f)
{
  (*this)[0] -= f;
  (*this)[1] -= f;
  (*this)[2] -= f;
  return *this;
}

vec4 &vec4::operator*=(const double &f)
{
  (*this)[0] *= f;
  (*this)[1] *= f;
  (*this)[2] *= f;
  return *this;
}

vec4 &vec4::operator/=(const double &f)
{
  (*this)[0] /= f;
  (*this)[1] /= f;
  (*this)[2] /= f;
  return *this;
}

vec4 vec4::operator-(const vec4 &p)
{
  vec4 v = (*this);
  v[0] -= p[0];
  v[1] -= p[1];
  v[2] -= p[2];
  v[3] =  p[3];
  return v;
};

vec4 vec4::operator+(const vec4 &p)
{
  vec4 v = (*this);
  v[0] += p[0];
  v[1] += p[1];
  v[2] += p[2];
  v[3] =  p[3];
  return v;
};

/**
 * performs oMat * this, and stores result in this
 */
vec4 &vec4::operator*=(const mat4 &oMat)
{
  double fSum;
  vec4 oResult;

  //
  // loop over each element of the new mat4
  //
  for (unsigned iRow = 0; iRow < 4; iRow++) {
    fSum = 0;
    for (unsigned i = 0; i < 4; i++)
      fSum += oMat[iRow][i] * (*this)[i];
    oResult[iRow] = fSum;
  }

  (*this) = oResult;

  return (*this);
};

/**
 * get/set methods to manipulate the coordinate
 */
double vec4::x() const
{
  return (*this)[0];
}
void vec4::x(double value)
{
  (*this)[0] = value;
}

double vec4::y() const
{
  return (*this)[1];
}
void vec4::y(double value)
{
  (*this)[1] = value;
}

double vec4::z() const
{
  return (*this)[2];
}
void vec4::z(double value)
{
  (*this)[2] = value;
}

double vec4::w() const
{
  return (*this)[3];
}
void vec4::w(double value)
{
  (*this)[3] = value;
}

/**
 * set all three coords of the vec4 at once
 * \param fX x coord
 * \param fY y coord
 * \param fZ z coord
 */
void vec4::set(double fX, double fY, double fZ, double fW)
{
  x(fX);
  y(fY);
  z(fZ);
  w(fW);
}

void vec4::set(double value, double fW)
{
  *this = value;
  w(fW);
}

void vec4::set(const double *p)
{
  memcpy(m_aData, p, 4 * sizeof(double));
}

/**
 * get the squares length of the vec4
 * \return length^2
 */
double vec4::lengthSq() const
{
  return x() * x() + y() * y() + z() * z();
}

/**
 * get the length of the vec4
 * \return length
 */
double vec4::length() const
{
  return sqrt(lengthSq());
}

/**
 * calculate the distance between two vec4s
 * \param oVec the other vec4
 * \return distance
 */
double vec4::distance(const vec4 &oVec) const
{
  double fDiffX = x() - oVec.x();
  double fDiffY = y() - oVec.y();
  double fDiffZ = z() - oVec.z();

  double fLengthSq = fDiffX * fDiffX + fDiffY * fDiffY + fDiffZ * fDiffZ;
  double fLength =  sqrt(fLengthSq);

  return fLength;
}

/**
 * calculate the squared distance between two vec4s
 * \param oVec the other vec4
 * \return distance
 */
double vec4::distanceSq(const vec4 &oVec) const
{
  return (*this - oVec).lengthSq();
}

/**
 * normalize the vec4, return *this to allow daisy chaining
 */
vec4 &vec4::normalize()
{
  double fLength = length();

  if (isPositive(fLength))
    (*this) /= fLength;

  return *this;
}

/**
 * Direct access to the memory
 */
double *vec4::c_data()
{
  return m_aData;
};

/**
 * Direct access to the memory
 */
const double *vec4::c_data() const
{
  return m_aData;
};

/**
 * Print content of vector to stdout
 */
void vec4::print()
{
  MLBO << "(" << x() << ", " << y() << ", " << z() << ", " << w() << ")" << endl;
};

vec4 operator-(const vec4 &p1, const vec4 &p2)
{
  return vec4(p1.x() - p2.x(), p1.y() - p2.y(), p1.z() - p2.z(), p1.w() - p2.w());
};

vec4 operator*(const mat4 &oM,  const vec4 &oV)
{
  double fSum;
  vec4 oResult;

  //
  // loop over each element of the new mat4
  //
  for (unsigned iRow = 0; iRow < 4; iRow++) {
    fSum = 0;
    for (unsigned i = 0; i < 4; i++)
      fSum += oM[iRow][i] * oV[i];
    oResult[iRow] = fSum;
  }

  return oResult;
};

vec4 operator/(const vec4 &V,  const double &F)
{
  return vec4(V.x() / F, V.y() / F, V.z() / F, V.w());
};

///////////////////////////////////////////////////////////////////////////////
// MAT4
///////////////////////////////////////////////////////////////////////////////

/**
 * Constructor
 */
mat4::mat4()
{
  loadIdentity();
};

mat4::mat4(const mat4 &rThat)
{
  memcpy(m_aData, rThat.m_aData, 16 * sizeof(double));
};

mat4::~mat4()
{
};

/**
 * Operators
 */
mat4 &mat4::operator=(const mat4 &rThat)
{
  memcpy(m_aData, rThat.m_aData, 16 * sizeof(double));
  return (*this);
};

mat4 &mat4::operator=(const double &fValue)
{
  for (unsigned int i = 0; i < 16; i++)
    m_aData[i] = fValue;

  return (*this);
};

/**
 * mat4 multiplication.
 */
//mat4& mat4::operator*=(const mat4& B);

/**
 * sets the mat4 to the identity mat4
 */
void mat4::loadIdentity()
{
  (*this) = double(0);
  (*this)[0][0] = (*this)[1][1] = (*this)[2][2] = (*this)[3][3] = double(1);
};

/**
 * adds a translation (from right)
 * \param fX x translation
 * \param fY y translation
 * \param fZ z translation
 */
mat4 &mat4::translate(double fX, double fY, double fZ)
{
  (*this)[0][3] += fX;
  (*this)[1][3] += fY;
  (*this)[2][3] += fZ;
  return *this;
};

/**
 * adds a translation
 * \param rVec reference to svt_vec44
 */
mat4 &mat4::translate(const vec4 &rVec)
{
  (*this)[0][3] += rVec.x();
  (*this)[1][3] += rVec.y();
  (*this)[2][3] += rVec.z();
  return *this;
};

/**
 * get the translation component
 * \return svt_vec44 object
 */
vec4 mat4::translation() const
{
  vec4 oVec((*this)[0][3], (*this)[1][3], (*this)[2][3]);
  return oVec;
};

/**
 * get the x translation
 * \return x translation
 */
double mat4::translationX() const
{
  return (*this)[0][3];
};

/**
 * set the x translation
 * \param fX the new x translation
 */
mat4 &mat4::setTranslationX(double fX)
{
  (*this)[0][3] = fX;
  return *this;
};

/**
 * get the y translation
 * \return y translation
 */
double mat4::translationY() const
{
  return (*this)[1][3];
};

/**
 * set the y translation
 * \param fY the new y translation
 */
mat4 &mat4::setTranslationY(double fY)
{
  (*this)[1][3] = fY;
  return *this;
};

/**
 * get the z translation
 * \return z translation
 */
double mat4::translationZ() const
{
  return (*this)[2][3];
};

/**
 * set the z translation
 * \param fZ the new z translation
 */
mat4 &mat4::setTranslationZ(double fZ)
{
  (*this)[2][3] = fZ;
  return *this;
};

/**
 * set translation
 * \param fX x component
 * \param fY y component
 * \param fZ z component
 */
mat4 &mat4::setTranslation(double fX, double fY, double fZ)
{
  setTranslationX(fX);
  setTranslationY(fY);
  return setTranslationZ(fZ);
};

/**
 * set the translation component
 * \param rVec reference to svt_vec44 object
 */
mat4 &mat4::setTranslation(const vec4 &rVec)
{
  return setTranslation(rVec.x(), rVec.y(), rVec.z());
};

/**
 * Range-unchecked Dereference Operator
 * intented to be used as matrix[iRow][iColumn]
 * this method returns a pointer to the first element of the iRow´th Row
 * the second [] is done by c
 * range-unchecked
 */
double *mat4::operator[](unsigned iRow)
{
  return &(m_aData[iRow * 4]);
};

/**
 * Range-unchecked Dereference Operator
 * intented to be used as matrix[iRow][iColumn]
 * this method returns a pointer to the first element of the iRow´th Row
 * the second [] is done by c
 * range-unchecked
 */
const double *mat4::operator[](unsigned iRow) const
{
  return &(m_aData[iRow * 4]);
};

/**
 * Direct access to the memory
 */
double *mat4::c_data()
{
  return m_aData;
};

/**
 * Jacobi transformation
 * \param rEigenvectors svt_matrix object to store the eigenvectors (columns)
 * \param rEigenvalues svt_vector object to store the eigenvalues
 */
bool mat4::jacobi(mat4 &rEigenvectors, vec4 &rEigenvalues)
{
  double  fSM;
  double  fTheta;
  double  fC, fS, fT;
  double  fTau;
  double  fH, fG;
  double  fThresh;
  double  aB[4];
  double  aZ[4];
  unsigned int  iP, iQ, i, j;
  double  oA[4][4];
  int iRots;

  // initializations
  for (i = 0; i < 4; i++) {
    aB[i] = rEigenvalues[i] = (*this)[i][i];
    aZ[i] = 0.0;

    for (j = 0; j < 4; j++) {
      rEigenvectors[i][j] = (i == j) ? 1.0f : 0.0f;
      oA[i][j] = (*this)[i][j];
    }
  }

  iRots = 0;

  for (i = 0; i < 50; i++) {
    fSM = 0.0;

    for (iP = 0; iP < 4 - 1; iP++) {
      for (iQ = iP + 1; iQ < 4; iQ++) {
        fSM += fabs(oA[iP][iQ]);
      }
    }

    if (fSM == 0.0)
      return false;

    fThresh = (i < 3 ? (.2 * fSM / (4 * 4)) : 0.0);

    for (iP = 0; iP < 4 - 1; iP++) {
      for (iQ = iP + 1; iQ < 4; iQ++) {
        fG = 100.0 * fabs(oA[iP][iQ]);

        if (i > 3 &&
            (fabs(rEigenvalues[iP]) + fG == fabs(rEigenvalues[iP])) &&
            (fabs(rEigenvalues[iQ]) + fG == fabs(rEigenvalues[iQ]))) {
          oA[iP][iQ] = 0.0;
        } else if (fabs(oA[iP][iQ]) > fThresh) {
          fH = rEigenvalues[iQ] - rEigenvalues[iP];

          if (fabs(fH) + fG == fabs(fH)) {
            fT = oA[iP][iQ] / fH;
          } else {
            fTheta = .5 * fH / oA[iP][iQ];
            fT = 1.0 / (fabs(fTheta) + sqrt(1 + fTheta * fTheta));
            if (fTheta < 0.0)  fT = -fT;
          }
          // End of computing tangent of rotation angle

          fC = 1.0 / sqrt(1.0 + fT * fT);
          fS = fT * fC;

          fTau = fS / (1.0 + fC);
          fH   = fT * oA[iP][iQ];

          aZ[iP]    -= fH;
          aZ[iQ]    += fH;

          rEigenvalues[iP] -= double(fH);
          rEigenvalues[iQ] += double(fH);

          oA[iP][iQ] = 0.0;

          for (j = 0; j < iP; j++) {
            fG = oA[j][iP];
            fH = oA[j][iQ];

            oA[j][iP] = fG - fS * (fH + fG * fTau);
            oA[j][iQ] = fH + fS * (fG - fH * fTau);
          }

          for (j = iP + 1; j < iQ; j++) {
            fG = oA[iP][j];
            fH = oA[j][iQ];

            oA[iP][j] = fG - fS * (fH + fG * fTau);
            oA[j][iQ] = fH + fS * (fG - fH * fTau);
          }

          for (j = iQ + 1; j < 4; j++) {
            fG = oA[iP][j];
            fH = oA[iQ][j];

            oA[iP][j] = fG - fS * (fH + fG * fTau);
            oA[iQ][j] = fH + fS * (fG - fH * fTau);
          }

          for (j = 0; j < 4; j++) {
            fG = rEigenvectors[j][iP];
            fH = rEigenvectors[j][iQ];

            rEigenvectors[j][iP] = double(fG - fS * (fH + fG * fTau));
            rEigenvectors[j][iQ] = double(fH + fS * (fG - fH * fTau));
          }
        }
        iRots++;
      }
    }
    for (iP = 0; iP < 4; iP++) {
      rEigenvalues[iP] = double(aB[iP] += aZ[iP]);
      aZ[iP] = 0;
    }
  }

  return true;

};

mat4 &mat_mult(const mat4 &A, const mat4 &B, mat4 &C)
{
  double sum;

  //
  // loop over each element of the new matrix
  //
  for (unsigned iRow = 0; iRow < 4; iRow++)
    for (unsigned iCol = 0; iCol < 4; iCol++)
      // compute each element as product of A(iRow,*) * B(*,iCol)
    {
      sum = 0;
      for (unsigned i = 0; i < 4; i++)
        sum += A[iRow][i] * B[i][iCol];
      C[iRow][iCol] = sum;
    }

  return C;
}

/**
 * Print content of matrix to stdout
 */
void mat4::print()
{
  for (unsigned int y = 0; y < 4; y++)
    printf("  %05.2f %05.2f %05.2f %05.2f\n", (*this)[y][0], (*this)[y][1], (*this)[y][2], (*this)[y][3]);
};

mat4 operator*(const mat4 &A, const mat4 &B)
{
  mat4 C;
  return mat_mult(A, B, C);
}

///////////////////////////////////////////////////////////////////////////////
// POINTCLOUD
///////////////////////////////////////////////////////////////////////////////

pointcloud::pointcloud() : sampled()
{
};

/**
 * Get all points in point cloud in an stl vector.
 * \return reference to vector of svt_vector4 objects
 */
vector< vec4 > &pointcloud::getPoints()
{
  return m_oPoints;
};

/**
 * Add a point to point cloud.
 * \param rVec vec4 object
 */
void pointcloud::addPoint(vec4 &rVec)
{
  m_oPoints.push_back(rVec);
};

/**
 * Add a point to point cloud.
 * \param fX x coord
 * \param fY y coord
 * \param fZ z coord
 */
void pointcloud::addPoint(double fX, double fY, double fZ)
{
  vec4 oVec(fX, fY, fZ);
  m_oPoints.push_back(oVec);
};

/**
 * Get a point out of point cloud.
 * \param iIndex index of point
 * \return reference to vec4 object
 */
vec4 &pointcloud::getPoint(unsigned int iIndex)
{
  return m_oPoints[iIndex];
};

/**
 * Size of point cloud.
 * \return size of pc
 */
unsigned int pointcloud::size() const
{
  return m_oPoints.size();
};

/**
 * Dereference operator (not range checked!).
 * \param iIndex index of point in point cloud
 */
vec4 &pointcloud::operator[](unsigned int iIndex)
{
  return m_oPoints[iIndex];
};

/**
 * product of matrix and point cloud
 */
pointcloud operator*(const mat4 &rM, pointcloud &rPC)
{
  pointcloud oPC;
  vec4 oP;

  for (unsigned int i = 0; i < rPC.size(); i++) {
    oP = rM * rPC[i];
    oPC.addPoint(oP);
    oPC.addMass(rPC.getMass(i));
  }

  return oPC;
};

/**
 * Delete all points in point cloud
 */
void pointcloud::delAllPoints()
{
  m_oPoints.clear();
};

/**
 * Calculate center of atoms (COA with all masses = 1).
 * \return svt_vector4 with the COA
 */
vec4 pointcloud::coa()
{
  unsigned int iNum = size();
  unsigned int i = 0;
  vec4 oCOAVec(0.0,  0.0,  0.0);

  for (i = 0; i < iNum; i++)
    oCOAVec =  oCOAVec + (*this)[i];

  oCOAVec /= (double) iNum;

  return oCOAVec;
};

/**
 * calculate RMSD between this and another PC. The PCs must be matched already! To get the minimal RMSD please use align() first!
 * \param rPC second point cloud
 */
double pointcloud::rmsd(pointcloud &rPC)
{
  if (size() != rPC.size()) {
    MLBO << "Cannot calculate rmsd - the point-clouds have different size!" << endl;
    exit(1);
  }

  unsigned int iNumA = size();
  unsigned int i;

  double fRMSD = 0.0;
  double fDist;

  for (i = 0; i < iNumA; i++) {
    fDist = (*this)[i].distanceSq(rPC[i]);

    fRMSD += fDist;
  }

  fRMSD *= 1.0 / (double)(iNumA);
  fRMSD = sqrt(fRMSD);

  return fRMSD;
};

/**
 * Calculate RMSD between this and another PC. The points are matched using the nearest neighbor relationship.
 * All the nearest neighbor distances are calculated and then sorted (slow!) and only the first N-percent are used for the matching and the following RMSD calculation.
 * The idea is that approximately N-percent of the points are outliers which would increase the RMSD significantly, although the overall deviation of the two point clouds is
 * actually small.
 * \param rPC second point cloud
 * \param fPercent percentage of neighbors that should be used for the rmsd calculation (value between 0 and 1!)
 */
double pointcloud::rmsd_NN_Outliers(pointcloud &rPC, double fPercent, mat4 *pMatrix)
{
  unsigned int iNumA = size();
  unsigned int i;

  double fRMSD = 0.0;
  double fDist;
  unsigned int iNeighbor;
  vector< neighbor > *pNeighbors = new vector< neighbor >;
  pointcloud oTmp;

  if (pMatrix != NULL)
    oTmp = (*pMatrix) * (*this);
  else
    oTmp = (*this);

  for (i = 0; i < iNumA; i++) {
    iNeighbor = rPC.nearestNeighbor(oTmp[i]);
    fDist = oTmp[i].distanceSq(rPC[iNeighbor]);

    pNeighbors->push_back(neighbor(fDist, i, iNeighbor));
  }

  sort(pNeighbors->begin(), pNeighbors->end());

  unsigned int iPercent = (unsigned int)((double)(iNumA) * fPercent);

  vector<int> aModelMatch;
  vector<int> aSceneMatch;

  for (i = 0; i < iPercent; i++) {
    fRMSD += (*pNeighbors)[i].getScore();

    aModelMatch.push_back((*pNeighbors)[i].getIndexA());
    aSceneMatch.push_back((*pNeighbors)[i].getIndexB());
  }

  fRMSD /= (double)(iPercent);
  fRMSD = sqrt(fRMSD);

  if (pMatrix != NULL) {
    mat4 oOrig = (*pMatrix);
    mat4 oMat;
    oMat = oTmp.kearsley(aModelMatch, aSceneMatch, oTmp, rPC);
    (*pMatrix) = oOrig * oMat;
  }

  delete pNeighbors;

  return fRMSD;
}

/**
 * Write pdb file.
 * \param pFilename pointer to array of char with the filename
 * \param bAppend if true, the pdb structure as append at the end of an existing structure file
 */
void pointcloud::writePDB(const char *pFilename, bool bAppend)
{
  unsigned int i;
  unsigned int iNum = size();
  FILE *pFile = NULL;

  if (!bAppend)
    pFile = fopen(pFilename, "w");
  else {
    pFile = fopen(pFilename, "a");
    fprintf(pFile, "END\n");
  }

  for (i = 0; i < iNum; i++)
    fprintf(
      pFile, "ATOM  %5i %2s%c%c%c%3s %c%4i%c  %8.*f%8.*f%8.*f%6.2f%6.2f %3s  %4s%2s%2s\n",
      i + 1,
      "QP",
      'D',
      'B',
      ' ',
      "QPDB",
      ' ',
      0,
      'Q',
      coord_precision(this->getPoint(i).x()), this->getPoint(i).x(),
      coord_precision(this->getPoint(i).y()), this->getPoint(i).y(),
      coord_precision(this->getPoint(i).z()), this->getPoint(i).z(),
      0.0,
      0.0,
      "QPDB",
      "QPDB",
      "QPDB",
      "QPDB"
    );

  fclose(pFile);
};

/**
 * Load a pdb file.
 * \param pointer to array of char with the filename
 */
void pointcloud::loadPDB(const char *pFilename)
{
  // open file
  columnreader oReader(pFilename);
  char *pBuffer;

  // clear all old points
  delAllPoints();

  if (!oReader.eof())
    while (1) {
      // read until next atom record reached
      while (1) {
        if (oReader.readLine() == false)
          return;

        pBuffer = oReader.extractString(0, 5);

        //is it an ATOM record?
        if (pBuffer[0] == 'A' && pBuffer[1] == 'T' &&  pBuffer[2] == 'O' &&  pBuffer[3] == 'M')
          break;
        // is it an HETATM record?
        if (pBuffer[0] == 'H' && pBuffer[1] == 'E' &&  pBuffer[2] == 'T' &&  pBuffer[3] == 'A' &&  pBuffer[4] == 'T' &&  pBuffer[5] == 'M')
          break;
      }

      //now parse the coordinates
      double fX = oReader.extractdouble(30, 37);
      double fY = oReader.extractdouble(38, 45);
      double fZ = oReader.extractdouble(46, 53);

      addPoint(fX, fY, fZ);

      adjustMass(oReader.extractString(12, 13));
    }
};

/**
 * Replace coordinates in a pdb file.
 * \param pointer to array of char with the input filename
 * \param pointer to array of char with the output filename
 */
void pointcloud::replacePDB(const char *pFilename_IN, const char *pFilename_OUT)
{
  // open input file
  columnreader oReader(pFilename_IN);
  char *pBuffer;
  char pCoords[256];
  unsigned int i = 0;

  // open output file
  FILE *pFile = fopen(pFilename_OUT, "w");

  if (!oReader.eof())
    while (1) {
      // read until next atom record reached
      while (1) {
        if (oReader.readLine() == false)
          return;

        pBuffer = oReader.getLine();

        //is it an ATOM record?
        if (pBuffer[0] == 'A' && pBuffer[1] == 'T' &&  pBuffer[2] == 'O' &&  pBuffer[3] == 'M')
          break;
        // is it an HETATM record?
        if (pBuffer[0] == 'H' && pBuffer[1] == 'E' &&  pBuffer[2] == 'T' &&  pBuffer[3] == 'A' &&  pBuffer[4] == 'T' &&  pBuffer[5] == 'M')
          break;
      }

      // now replace only the coordinate columns
      sprintf(pCoords, "%8.3f%8.3f%8.3f", this->getPoint(i).x(), this->getPoint(i).y(), this->getPoint(i).z());
      memcpy(pBuffer + 30, pCoords, 23);

      fprintf(pFile, "%s", pBuffer);

      i++;
      if (i > size())
        break;
    }

  fclose(pFile);
};

/**
 * calculate the average nearest neighbor distance
 * in order to reduce the complexity a random test is done and once the average stabilizes, the search is stopped.
 * \param fPrecision if average does not change by more than fPrecision between two iterations the calculation is stopped
 * \return average nearest neighbor distance
 */
double pointcloud::averageNNDistance(double fPrecision)
{
  unsigned int iNum = size();
  double fNum = (double)(iNum);
  unsigned int iTest = 0;
  double fCount = 0.0;
  double fAvg = 0.0;
  double fAvgDiv = 0.0;
  double fAvgDivOld = 1.0E10;
  double fMinDist;
  double fDist;
  unsigned int iOldTest = 0;

  if (fNum > 0) {
    // loop until desired precision is reached...
    while (fabs(fAvgDiv - fAvgDivOld) > fPrecision && fCount < 1000000.0) {
      while (iTest == iOldTest) iTest = (unsigned int)((rand() / (RAND_MAX + 1.0)) * fNum);
      iOldTest = iTest;

      // determine nearest neighbor
      fMinDist = 1.0E10;
      for (unsigned int i = 0; i < iNum; i++) {
        fDist = (*this)[i].distance((*this)[iTest]);
        if (fDist < fMinDist && fDist != 0.0) fMinDist = fDist;
      }

      // calculate average
      fAvg += fMinDist;
      fAvgDivOld = fAvgDiv;
      fAvgDiv = fAvg / fCount;
      fCount += 1.0;
    }
  }

  return fAvgDiv;
};

/**
 * Set tree pruning parameter
 * \param fEpsilon epsilon
 */
void pointcloud::setEpsilon(double fEpsilon)
{
  m_fEpsilon = fEpsilon;
};
/**
 * Get tree pruning parameter.
 * \param fEpsilon epsilon
 */
double pointcloud::getEpsilon() const
{
  return m_fEpsilon;
};
/**
 * set tolorance distance for the anchor determination
 * \param fLambda lambda
 */
void pointcloud::setLambda(double fLambda)
{
  m_fLambda = fLambda;
};
/**
 * get tolorance distance for the anchor determination
 * \return lambda
 */
double pointcloud::getLambda() const
{
  return m_fLambda;
};
/**
 * set nearest neighbor matching zone size
 * \param fGamma gamma
 */
void pointcloud::setGamma(double fGamma)
{
  m_fGamma = fGamma;
};
/**
 * get nearest neighbor matching zone size
 * \return gamma
 */
double pointcloud::getGamma() const
{
  return m_fGamma;
};

/**
 * set the maximal size of the matching zone
 * \param iZoneSize
 */
void pointcloud::setZoneSize(unsigned int iZoneSize)
{
  m_iZoneSize = iZoneSize;
};
/**
 * get the maximal size of the matching zone
 * \return maximal size of matching zone
 */
unsigned int pointcloud::getZoneSize() const
{
  return m_iZoneSize;
};

/**
 * set the maximal number of wildcard matches
 * \param iMaxNoMatch maximal number of wildcard matches
 */
void pointcloud::setWildcards(unsigned int iMaxNoMatch)
{
  m_iMaxNoMatch = iMaxNoMatch;
};
/**
 * get the maximal number of wildcard matches
 * \return maximal number of wildcard matches
 */
unsigned int pointcloud::getWildcards() const
{
  return m_iMaxNoMatch;
};

/**
 * set the penalty for wildcard matches
 * \param fSkipPenalty penalty for a single wildcard
 */
void pointcloud::setSkipPenalty(double fSkipPenalty)
{
  m_fSkipPenalty = fSkipPenalty;
}
/**
 * get the penalty for wildcard matches
 * \return penalty for a single wildcard
 */
double pointcloud::getSkipPenalty() const
{
  return m_fSkipPenalty;
}

/**
 * if two solutions are very close only the one with the higher score is considered. The other solutions are removed.
 * \param minimal distance between solutions
 */
void pointcloud::setUnique(double fUnique)
{
  m_fUnique = fUnique;
}
/**
 * if two solutions are very close only the one with the higher score is considered. The other solutions are removed.
 * \return minimal distance between solutions
 */
double pointcloud::getUnique() const
{
  return m_fUnique;
};

/**
 * set the number of runs. Each time a different set of anchor points get selected from the probe structure, beginning with the three points furthest away from each other and the COA.
 * \param iRuns number of runs
 */
void pointcloud::setRuns(unsigned int iRuns)
{
  m_iRuns = iRuns;
};
/**
 * get the number of runs. Each time a different set of anchor points get selected from the probe structure, beginning with the three points furthest away from each other and the COA.
 * \return number of runs
 */
unsigned int pointcloud::getRuns() const
{
  return m_iRuns;
};

/**
 * Set the next point selection scheme to COA
 */
void pointcloud::setNextPointCOA(bool bNextPointCOA)
{
  m_bNextPointCOA = bNextPointCOA;
};

/**
 * match two point clouds.
 * This function will perform a full N->M matching based on an achor-based search.
 * \param rPC second point cloud
 * \param rMatch vector of unsigned ints with the indices of the points of second cloud to which the points of this cloud are matched. This vector will get erased and then filled during the matching.
 * \param rMatrices vector of svt_matrix4 objects with the transformations according to the rMatch
 */
void pointcloud::match(pointcloud &rPC, vector< matchresult > &rMatch)
{
  // calculate COA
  vec4 oCOA = coa();
  unsigned int i, j;

  unsigned int iFstIndex = 0;
  unsigned int iSndIndex = 0;
  unsigned int iTrdIndex = 0;

  vector< int  > aOldAnchors;
  vector< mat4 > aBestSolution;
  vector< int  > aBestMatch;

  // some status output
  //MLBO << endl;
  //MLBO << "Point-cloud matching..." << endl;
  //MLBO << "Subcomponent has " << size() << " points " << endl;
  //MLBO << "Gets docked into a structure with " << rPC.size() << " points " << endl;
  //MLBO << endl;

  // size of the two pointclouds
  unsigned int iNumA = size();

  // sum of all seed atoms evaluated
  unsigned int iSumSeeds = 0;

  // array with all the results from all runs
  vector< matchresult > aResults;

  for (unsigned int iRun = 0; iRun < m_iRuns; iRun++) {
    // calc best anchor of model point cloud
    // 1st: which one is furthest away from the center
    double fLength = 0.0;
    double fDist = 0.0;
    for (i = 0; i < iNumA; i++) {
      fDist = (*this)[i].distance(oCOA);

      if (fDist > fLength) {
        if (find(aOldAnchors.begin(), aOldAnchors.end(), (int)(i)) == aOldAnchors.end()) {
          fLength = fDist;
          iFstIndex = i;
        }
      }
    }

    // 2nd: which one is furthest away from the first point?
    double fTmpDist = 0.0;
    fDist = 0.0;
    for (i = 0; i < iNumA; i++) {
      fTmpDist = (*this)[i].distance((*this)[iFstIndex]);

      if (i != iFstIndex && fTmpDist > fDist) {
        if (find(aOldAnchors.begin(), aOldAnchors.end(), (int)(i)) == aOldAnchors.end()) {
          fDist = fTmpDist;
          iSndIndex = i;
        }
      }
    }

    // 3rd: which one is furthest away from the first two?
    double fSumDist = 0.0;
    for (i = 0; i < iNumA; i++) {
      if (i != iFstIndex && i != iSndIndex) {
        fDist = (*this)[i].distance((*this)[iFstIndex]) + (*this)[i].distance((*this)[iSndIndex]);
        if (fDist > fSumDist) {
          if (find(aOldAnchors.begin(), aOldAnchors.end(), (int)(i)) == aOldAnchors.end()) {
            fSumDist = fDist;
            iTrdIndex = i;
          }
        }
      }
    }

    if (find(aOldAnchors.begin(), aOldAnchors.end(), (int)(iFstIndex)) != aOldAnchors.end() &&
        find(aOldAnchors.begin(), aOldAnchors.end(), (int)(iSndIndex)) != aOldAnchors.end() &&
        find(aOldAnchors.begin(), aOldAnchors.end(), (int)(iTrdIndex)) != aOldAnchors.end()) {
      //MLBO << "Max. runs is higher, but no more anchor points available to try." << endl;
      break;
    }

    //MLBO << "Anchors: (" << iFstIndex << ", " << iSndIndex << ", " << iTrdIndex << ")" << endl;

    // run the matching
    vector< mat4 > aMatrix;
    vector< matchresult > aMatch;
    match(rPC, aMatch, (unsigned int)(iFstIndex), (unsigned int)(iSndIndex), (unsigned int)(iTrdIndex));

    // append all the results to the global result array
    for (unsigned int k = 0; k < aMatch.size(); k++)
      aResults.push_back(aMatch[k]);

    iSumSeeds += m_iNumSeeds;

    // add the anchors to the old anchors list
    aOldAnchors.push_back(iFstIndex);
    aOldAnchors.push_back(iSndIndex);
    aOldAnchors.push_back(iTrdIndex);
  }

  // analyze the found solutions
  sort(aResults.begin(), aResults.end());

  // filter
  rMatch.clear();
  vector< matchresult > aFinal;
  vector< vec4        > aComs;
  vec4 oTransVec;
  unsigned int iNum = aResults.size();
  for (i = 0; i < iNum; i++) {
    oTransVec = aResults[i].getMatrix() * oCOA;

    for (j = 0; j < aComs.size(); j++) {
      if (aComs[j].distance(oTransVec) < m_fUnique)
        break;
    }

    if (j >= aComs.size()) {
      aComs.push_back(oTransVec);
      rMatch.push_back(aResults[i]);
    }
  }

  m_iNumSeeds = iSumSeeds;
}

/**
 * match two point clouds.
 * This function will perform a full N->M matching based on an achor-based search.
 * \param rPC second point cloud
 * \param pMatch vector of unsigned ints with the indices of the points of second cloud to which the points of this cloud are matched. This vector will get erased and then filled during the matching.
 * \param iAnchorA index of first anchor point (if not provided the routine will find a set of appropriate anchor points)
 * \param iAnchorB index of second anchor point (if not provided the routine will find a set of appropriate anchor points)
 * \param iAnchorC index of third anchor point (if not provided the routine will find a set of appropriate anchor points)
 * \return copy of this cloud matched to rPC
 */
void pointcloud::match(pointcloud &rPC, vector< matchresult > &rMatch, unsigned int iAnchorA, unsigned int iAnchorB, unsigned int iAnchorC)
{
  unsigned int i, j;

  // start time
  //long iStartTime = svt_getToD();

  // empty the match
  rMatch.clear();

  // average nn distance
  double fEpsilon = averageNNDistance(0.1);
  setEpsilon(fEpsilon + (1.0 * fEpsilon));

  // size of the two pointclouds
  unsigned int iNumB = rPC.size();

  // anchors
  unsigned int iFstIndex = iAnchorA;
  unsigned int iSndIndex = iAnchorB;
  unsigned int iTrdIndex = iAnchorC;

  // distances between anchors
  double fDistAB = (*this)[iFstIndex].distance((*this)[iSndIndex]);
  double fDistBC = (*this)[iSndIndex].distance((*this)[iTrdIndex]);
  double fDistAC = (*this)[iFstIndex].distance((*this)[iTrdIndex]);
  vector< vec4 > aSeeds;

  // now find appropriate seed anchors in the scene, based on the distances...
  for (i = 0; i < iNumB; i++)
    for (j = 0; j < iNumB; j++)
      if (i != j) {
        if (fabs(rPC[i].distance(rPC[j]) - fDistAB) < m_fLambda) {
          vec4 oVec(i, j, -1);
          aSeeds.push_back(oVec);

        }
      }

  m_iNumSeeds = aSeeds.size();
  vector< vec4 > aFinalSeeds;

  for (i = 0; i < m_iNumSeeds; i++)
    for (j = 0; j < iNumB; j++) {
      if (aSeeds[i].x() != (int)(j) && aSeeds[i].y() != (int)(j)) {
        if (
          (fabs(rPC[j].distance(rPC[(unsigned int)(aSeeds[i].x()) ]) - fDistAC) < m_fLambda &&
           fabs(rPC[j].distance(rPC[(unsigned int)(aSeeds[i].y()) ]) - fDistBC) < m_fLambda) ||
          (fabs(rPC[j].distance(rPC[(unsigned int)(aSeeds[i].y()) ]) - fDistAC) < m_fLambda &&
           fabs(rPC[j].distance(rPC[(unsigned int)(aSeeds[i].x()) ]) - fDistBC) < m_fLambda)
        ) {
          vec4 oVec = aSeeds[i];
          oVec.z(j);
          aFinalSeeds.push_back(oVec);
        }
      }
    }

  m_iNumSeeds = aFinalSeeds.size();

  //MLBO << "Number of potential anchor matches: " << aFinalSeeds.size() << endl;

  // investigate the quality of the anchor matches
  mat4 oMatrix;
  vector<int> aMatches;
  vector< matchresult > aResults;
  double fRMSD = 0.0;

  if (m_iNumSeeds > 0) {
    //MLBO << "Searching" << endl;
    //MLBO << ".";
  }

  for (i = 0; i < m_iNumSeeds; i++) {
    vector<int> aModelAnchor;
    aModelAnchor.push_back(iFstIndex);
    aModelAnchor.push_back(iSndIndex);
    aModelAnchor.push_back(iTrdIndex);
    vector<int> aSceneAnchor;
    aSceneAnchor.push_back((unsigned int)(aFinalSeeds[i].x()));
    aSceneAnchor.push_back((unsigned int)(aFinalSeeds[i].y()));
    aSceneAnchor.push_back((unsigned int)(aFinalSeeds[i].z()));

    oMatrix = optimize(NULL, &aModelAnchor, &aSceneAnchor, *this, rPC, true);

    if (aModelAnchor.size() > 3) {
      fRMSD = calcRMSD(*this, rPC, oMatrix, &aModelAnchor, &aSceneAnchor);
      aResults.push_back(matchresult(fRMSD, oMatrix, aModelAnchor, aSceneAnchor));
    }
  }
  //if (m_iNumSeeds > 0)
  //    cout << endl;

  rMatch.clear();

  if (aResults.size() == 0) {
    //MLBO << "  Error - there was no valid result in this run!" << endl;
    return;
  }

  // sort the resulting matrices
  sort(aResults.begin(), aResults.end());
  // copy the 10 best results into a matrix array
  for (i = 0; i < aResults.size() && i < 1000; i++)
    rMatch.push_back(aResults[i]);
};

/**
 * internal convenience function - no real difference to rmsd - only here we don't need to create/copy a pointcloud
 * plus here we take the matching into account! If NOMATCH this point will not be used.
 */
double pointcloud::calcRMSD(pointcloud &rModel, pointcloud &rScene, mat4 &rMatrix, vector<int> *pModelMatch, vector<int> *pSceneMatch)
{
  unsigned int iNum = pModelMatch->size();
  unsigned int i;

  if (iNum != (*pSceneMatch).size()) {
    MLBO << "ERROR: calcRMSD: pSceneMatch does not contain enough point matches!!!" << endl;
    exit(1);
  }

  double fRMSD = 0.0;
  int iCount = 0;
  vec4 oVec;

  for (i = 0; i < iNum; i++) {
    if ((*pSceneMatch)[i] != NOMATCH) {
      oVec = rMatrix * rModel[(*pModelMatch)[i]];
      fRMSD += oVec.distanceSq(rScene[(*pSceneMatch)[i]]);
      iCount++;
    }
  }

  fRMSD /= iCount;
  fRMSD = sqrt(fRMSD);

  return fRMSD;
};

/**
 * optimize - subroutine for the match() procedure
 */
mat4 pointcloud::optimize(optState *pState, vector<int> *pModelMatch, vector<int> *pSceneMatch, pointcloud &rModel, pointcloud &rScene, bool bInit)
{
  if (pState == NULL)
    pState = new optState;

  if (bInit == true) {
    pState->fBestRMSD = 1.0E10;
    pState->iNumM = rModel.size();
    pState->iNumS = rScene.size();
    pState->iLevel = 0;

    if (m_bNextPointCOA == true) {
      // calculate distances of all model vectors to COA
      vec4 oCOA = rModel.coa();
      pState->aDistToCOA.clear();
      for (unsigned int i = 0; i < pState->iNumM; i++)
        pState->aDistToCOA.push_back(pc_dist(rModel[i].distance(oCOA), i));

      sort(pState->aDistToCOA.begin(), pState->aDistToCOA.end());
    }

  } else
    pState->iLevel++;

  // find an unmatched model vector...
  vec4 oVec;
  unsigned int iIndex = 0, j = 0, k = 0;

  if (pModelMatch->size() < rModel.size()) {
    // get optimal transformation for the old match
    mat4 oMatrix = kearsley(*pModelMatch, *pSceneMatch, rModel, rScene);

    // look for an unmatched point by moving towards the COA...
    if (m_bNextPointCOA) {
      unsigned int i;
      for (i = pState->aDistToCOA.size(); i > 0; i--)
        if (find(pModelMatch->begin(), pModelMatch->end(), (int)(pState->aDistToCOA[i - 1].getIndex())) == pModelMatch->end())
          break;

      iIndex = pState->aDistToCOA[i - 1].getIndex();

    } else {

      // look for an unmatched point that is closest to a scenepoint
      double fMinDist = 1.0E10;
      double fDist = 0.0;
      for (unsigned int i = 0; i < rModel.size(); i++) {
        // which model point is unmatched?
        if (find(pModelMatch->begin(), pModelMatch->end(), (int)(i)) == pModelMatch->end()) {
          // now loop over the scene...
          for (k = 0; k < rScene.size(); k++) {
            // which scene point is unmatched?
            if (find(pSceneMatch->begin(), pSceneMatch->end(), (int)(k)) == pSceneMatch->end()) {
              // calculate distance to unmatched model vector i
              fDist = rModel[i].distance(rScene[k]);

              // closest so far?
              if (fDist < fMinDist) {
                fMinDist = fDist;
                iIndex = i;
              }
            }
          }
        }
      }
    }

    // transform
    oVec = oMatrix * rModel[iIndex];

    // now lets find the nearest neighbors in the scene...
    vector< pc_dist > aNeighbors;
    double fDist = 0.0;
    for (j = 0; j < pState->iNumS; j++) {
      // is this scene point unmatched?
      for (k = 0; k < pSceneMatch->size(); k++) {
        //MLBO << "pSceneMatch[" << k << "]:" << (*pSceneMatch)[k] << endl;

        if ((*pSceneMatch)[k] == (int)j)
          break;
      }
      // yes...
      if (k >= pSceneMatch->size()) {
        fDist = rScene[j].distance(oVec);

        //MLBO << j << " - fDist: " << fDist << " m_fGamma: " << m_fGamma << endl;

        if (fDist < m_fGamma)
          aNeighbors.push_back(pc_dist(fDist, j));
      }
    }

    // test all the points
    unsigned int iNumN = aNeighbors.size();
    sort(aNeighbors.begin(), aNeighbors.end());
    unsigned int iCount = 0;

    // do we have potential matching points?
    if (iNumN > 0) {
      for (j = 0; j < iNumN; j++) {
        if (iCount < m_iZoneSize) {
          vector<int> *pTmpMatchModel = new vector<int>;
          *pTmpMatchModel = *pModelMatch;
          vector<int> *pTmpMatchScene = new vector<int>;
          *pTmpMatchScene = *pSceneMatch;
          pTmpMatchModel->push_back(iIndex);
          pTmpMatchScene->push_back(aNeighbors[j].getIndex());
          optimize(pState, pTmpMatchModel, pTmpMatchScene, rModel, rScene, false);
          delete pTmpMatchModel;
          delete pTmpMatchScene;
          iCount++;
        }
      }
      // no, no potential matches, so we can try a NOMATCH, but not too many...
    }
    if (iNumN == 0 || aNeighbors[0].getScore() > 3.0) {
      unsigned int iNoMatch = 0;
      for (j = 0; j < pSceneMatch->size(); j++)
        if ((*pSceneMatch)[j] == NOMATCH)
          iNoMatch++;

      // if there are no more than iMaxNoMatch blind matches, try them
      if (iNoMatch < m_iMaxNoMatch) {
        vector<int> *pTmpMatchModel = new vector<int>;
        *pTmpMatchModel = *pModelMatch;
        vector<int> *pTmpMatchScene = new vector<int>;
        *pTmpMatchScene = *pSceneMatch;
        pTmpMatchModel->push_back(iIndex);
        pTmpMatchScene->push_back(NOMATCH);
        optimize(pState, pTmpMatchModel, pTmpMatchScene, rModel, rScene, false);
        delete pTmpMatchModel;
        delete pTmpMatchScene;
      }
    }

    // no unmatched vector left, so lets evaluate this solution...
  } else {

    if (pModelMatch->size() == rModel.size()) {
      mat4 oMatrix = kearsley(*pModelMatch, *pSceneMatch, rModel, rScene);
      double fRMSD = 0.0;
      if (m_iMaxNoMatch == 0) {
        fRMSD = calcRMSD(rModel, rScene, oMatrix, pModelMatch, pSceneMatch);

      } else {

        double fPercent = (double)(rModel.size()) / (double)(m_iMaxNoMatch);
        if (fPercent > 1.0)
          fPercent = 1.0;
        fRMSD = rModel.rmsd_NN_Outliers(rScene, fPercent);
      }

      // do we have a winner?
      if (fRMSD < pState->fBestRMSD) {
        pState->fBestRMSD = fRMSD;
        pState->oBestMatrix = oMatrix;
        pState->aBestModelMatch = *pModelMatch;
        pState->aBestSceneMatch = *pSceneMatch;

        unsigned int iNoMatch = 0;
        for (j = 0; j < pSceneMatch->size(); j++)
          if ((*pSceneMatch)[j] == NOMATCH)
            iNoMatch++;
      }
    }

  }

  // OK, do we exit this now? Then copy the best match...
  if (bInit == true) {
    if (pState->aBestModelMatch.size() == rModel.size()) {
      (*pModelMatch) = pState->aBestModelMatch;
      (*pSceneMatch) = pState->aBestSceneMatch;
    }
  }

  return pState->oBestMatrix;
};

/**
 * find the nearest neighbor to a query point in the point cloud
 * \param rVec reference to svt_vector4 object - the query point
 * \return index to nearest point in point cloud
 */
unsigned int pointcloud::nearestNeighbor(vec4 &rVec)
{
  unsigned int iNumA = size();

  double fMinDist = 1.0E10;
  unsigned int iMinDist = 0;
  double fDist;

  for (unsigned int i = 0; i < iNumA; i++) {
    fDist = (*this)[i].distance(rVec);

    if (fDist < fMinDist) {
      iMinDist = i;
      fMinDist = fDist;
    }
  }

  return iMinDist;
};

/**
 * Least-squares fit of a model point set to a scene point set.
 *
 * The method is based on Kearsley's quaternion method. The method by Kabsch produces
 * in some cases wrong results, e.g. try to align the following pdb onto itself:
 *
 * ATOM   1421  CA  ASP   195      23.907   2.453  21.715  1.00  4.10
 * ATOM   1627  CA  GLU   221      20.828   3.631  18.679  1.00  3.82
 * ATOM    295  CA  GLY    44      48.823  29.146  -1.483  1.00 22.85
 * ATOM    313  CA  GLY    47      42.862  22.490   3.028  1.00  8.42
 *
 * This should of course yield the identity operator and an RMSD of
 * 0.0 A, but the Kabsch routine actually yields ... 1.368 A, plus an
 * operator that isn't anywhere near identity !
 * Found on the O mailing list.
 *
 * \param rM vector of ints with the indices of the first (model) vectors that are used
 * \param rS vector of ints with the indices of the second (scene) vectors that are used
 * \param rModel vector of svt_vector4 with the model vectors
 * \param rScene vector of svt_vector4 with the scene vectors
 * \return matrix with optimal transformation of the model into the scene
 */

mat4 pointcloud::kearsley(vector< int > &rMM, vector< int > &rSM, pointcloud &rModel, pointcloud &rScene)
{
  unsigned int i;
  unsigned int iCount = rMM.size();

  vector< int > rM;
  vector< int > rS;

  for (i = 0; i < iCount; i++) {
    if (rSM[i] != NOMATCH) {
      rM.push_back(rMM[i]);
      rS.push_back(rSM[i]);
    }
  }
  iCount = rM.size();

  // precomp.: calc COA for both model and scene and transformation
  vec4 oModelCOA;
  for (i = 0; i < iCount; i++)
    oModelCOA = oModelCOA + rModel[rM[i]];
  oModelCOA = oModelCOA / (double)(iCount);

  vec4 oSceneCOA;
  for (i = 0; i < iCount; i++)
    oSceneCOA = oSceneCOA + rScene[rS[i]];
  oSceneCOA = oSceneCOA / (double)(iCount);

  // first step: setup Q matrix (see kearsley acta cryst a45)

  // distance plus and minus
  vec4 oDP;
  vec4 oDM;

  // setup matrix
  mat4 oQ;
  // as svt_matrix creates an identity matrix, we have to set the diagonal elements to 0
  oQ[0][0] = 0.0;
  oQ[1][1] = 0.0;
  oQ[2][2] = 0.0;
  oQ[3][3] = 0.0;

  // now construct matrix
  for (i = 0; i < iCount; i++) {
    oDM = (rModel[rM[i]] - oModelCOA) - (rScene[rS[i]] - oSceneCOA);
    oDP = (rModel[rM[i]] - oModelCOA) + (rScene[rS[i]] - oSceneCOA);

    oQ[0][0] += (oDM.x() * oDM.x()) + (oDM.y() * oDM.y()) + (oDM.z() * oDM.z());
    oQ[1][1] += (oDM.x() * oDM.x()) + (oDP.y() * oDP.y()) + (oDP.z() * oDP.z());
    oQ[2][2] += (oDP.x() * oDP.x()) + (oDM.y() * oDM.y()) + (oDP.z() * oDP.z());
    oQ[3][3] += (oDP.x() * oDP.x()) + (oDP.y() * oDP.y()) + (oDM.z() * oDM.z());

    oQ[0][1] += (oDP.y() * oDM.z()) - (oDM.y() * oDP.z());
    oQ[0][2] += (oDM.x() * oDP.z()) - (oDP.x() * oDM.z());
    oQ[0][3] += (oDP.x() * oDM.y()) - (oDM.x() * oDP.y());

    oQ[1][0] += (oDP.y() * oDM.z()) - (oDM.y() * oDP.z());
    oQ[1][2] += (oDM.x() * oDM.y()) - (oDP.x() * oDP.y());
    oQ[1][3] += (oDM.x() * oDM.z()) - (oDP.x() * oDP.z());

    oQ[2][0] += (oDM.x() * oDP.z()) - (oDP.x() * oDM.z());
    oQ[2][1] += (oDM.x() * oDM.y()) - (oDP.x() * oDP.y());
    oQ[2][3] += (oDM.y() * oDM.z()) - (oDP.y() * oDP.z());

    oQ[3][0] += (oDP.x() * oDM.y()) - (oDM.x() * oDP.y());
    oQ[3][1] += (oDM.x() * oDM.z()) - (oDP.x() * oDP.z());
    oQ[3][2] += (oDM.y() * oDM.z()) - (oDP.y() * oDP.z());
  }

  // second step: solve matrix using Jacobi transformation
  mat4 oEigenvectors;
  vec4 oEigenvalues;
  oQ.jacobi(oEigenvectors, oEigenvalues);

  // third step: sort the eigenvectors according to their eigenvalues
  eigensort(oEigenvectors, oEigenvalues);

  // forth step: generate transformation matrix
  mat4 oFinal;
  mat4 oRot;

  oFinal[0][3] -= oModelCOA[0];
  oFinal[1][3] -= oModelCOA[1];
  oFinal[2][3] -= oModelCOA[2];

  oRot[0][0] = (oEigenvectors[0][3] * oEigenvectors[0][3]) + (oEigenvectors[1][3] * oEigenvectors[1][3]) - (oEigenvectors[2][3] * oEigenvectors[2][3]) - (oEigenvectors[3][3] * oEigenvectors[3][3]);
  oRot[1][0] = 2.0 * (oEigenvectors[1][3] * oEigenvectors[2][3] + oEigenvectors[0][3] * oEigenvectors[3][3]);
  oRot[2][0] = 2.0 * (oEigenvectors[1][3] * oEigenvectors[3][3] - oEigenvectors[0][3] * oEigenvectors[2][3]);

  oRot[0][1] = 2.0 * (oEigenvectors[1][3] * oEigenvectors[2][3] - oEigenvectors[0][3] * oEigenvectors[3][3]);
  oRot[1][1] = oEigenvectors[0][3] * oEigenvectors[0][3] + oEigenvectors[2][3] * oEigenvectors[2][3] - oEigenvectors[1][3] * oEigenvectors[1][3] - oEigenvectors[3][3] * oEigenvectors[3][3];
  oRot[2][1] = 2.0 * (oEigenvectors[2][3] * oEigenvectors[3][3] + oEigenvectors[0][3] * oEigenvectors[1][3]);

  oRot[0][2] = 2.0 * (oEigenvectors[1][3] * oEigenvectors[3][3] + oEigenvectors[0][3] * oEigenvectors[2][3]);
  oRot[1][2] = 2.0 * (oEigenvectors[2][3] * oEigenvectors[3][3] - oEigenvectors[0][3] * oEigenvectors[1][3]);
  oRot[2][2] = oEigenvectors[0][3] * oEigenvectors[0][3] + oEigenvectors[3][3] * oEigenvectors[3][3] - oEigenvectors[1][3] * oEigenvectors[1][3] - oEigenvectors[2][3] * oEigenvectors[2][3];

  oFinal = oRot * oFinal;

  oFinal[0][3] += oSceneCOA[0];
  oFinal[1][3] += oSceneCOA[1];
  oFinal[2][3] += oSceneCOA[2];

  return oFinal;
};

/**
 * sort eigenvectors according to their eigenvalues
 * \param pEigenvectors pointer to svt_matrix with the eigenvectors as column vectors
 * \param pEigenvalues pointer to svt_vector with the eigenvalues
 */
void pointcloud::eigensort(mat4 &rEigenvectors, vec4 &rEigenvalues)
{
  double fP;
  int i, j, k;

  for (i = 0; i < 4; i++) {
    fP = rEigenvalues[k = i];

    for (j = i + 1; j < 4; j++)
      if (rEigenvalues[j] >= fP)
        fP = rEigenvalues[k = j];

    if (k != i) {
      rEigenvalues[k] = rEigenvalues[i];
      rEigenvalues[i] = fP;
      for (j = 0; j < 4; j++) {
        fP = rEigenvectors[j][i];
        rEigenvectors[j][i] = rEigenvectors[j][k];
        rEigenvectors[j][k] = fP;
      }
    }
  }
};

/**
 * blur the pdb structure and thereby create an artificial low-resolution map
 * \param fWidth voxel width of the target map
 * \param fResolution resolution of the target map
 */
volume *pointcloud::blur(double fWidth, double fResolution)
{
  // create gaussian kernel
  volume oKernel;
  oKernel.createKernel(fWidth, fResolution);

  // bring lattice into register with origin
  vec4 oMin = getMinCoord();
  vec4 oMax = getMaxCoord();

  double fMinx = (fWidth * floor(oMin.x() / fWidth));
  double fMaxx = (fWidth * ceil(oMax.x() / fWidth));
  double fMiny = (fWidth * floor(oMin.y() / fWidth));
  double fMaxy = (fWidth * ceil(oMax.y() / fWidth));
  double fMinz = (fWidth * floor(oMin.z() / fWidth));
  double fMaxz = (fWidth * ceil(oMax.z() / fWidth));

  // allocate protein density map
  int iMargin = (int) ceil((double)(oKernel.getSizeX()) / 2.0);
  int iExtx = (int)(ceil((fMaxx - fMinx) / fWidth)) + (2 * iMargin) + 1;
  int iExty = (int)(ceil((fMaxy - fMiny) / fWidth)) + (2 * iMargin) + 1;
  int iExtz = (int)(ceil((fMaxz - fMinz) / fWidth)) + (2 * iMargin) + 1;

  volume *pMap = new volume(iExtx, iExty, iExtz);
  pMap->setValue(0.0);
  pMap->setWidth(fWidth);

  // interpolate structure to protein map and keep track of variability - i.e. create a volumetric map with peaks at the atomic positions...
  double fVarp = 0.0;
  double fGx, fGy, fGz, fA, fB, fC;
  int iX0, iY0, iZ0, iX1, iY1, iZ1;
  unsigned int i;
  unsigned int iAtomsNum = size();
  double fSumMass = 0.0;

  for (i = 0; i < iAtomsNum; i++) {
    // compute position within grid
    fGx = iMargin + (((*this)[i].x() - fMinx) / fWidth);
    fGy = iMargin + (((*this)[i].y() - fMiny) / fWidth);
    fGz = iMargin + (((*this)[i].z() - fMinz) / fWidth);

    iX0 = (int)(floor(fGx));
    iY0 = (int)(floor(fGy));
    iZ0 = (int)(floor(fGz));
    iX1 = iX0 + 1;
    iY1 = iY0 + 1;
    iZ1 = iZ0 + 1;

    // interpolate
    fA = (double)(iX1) - fGx;
    fB = (double)(iY1) - fGy;
    fC = (double)(iZ1) - fGz;

    pMap->setValue(iX0, iY0, iZ0, 
                   pMap->getValue(iX0, iY0, iZ0) + fA * fB * fC * getMass(i));
    fVarp += getMass(i) * fA * fB * fC * ((1.0 - fA) * (1.0 - fA) + 
                                          (1.0 - fB) * (1.0 - fB) + 
                                          (1.0 - fC) * (1.0 - fC));
    pMap->setValue(iX0, iY0, iZ1, 
                   pMap->getValue(iX0, iY0, iZ1) + fA * fB * (1.0 - fC) * 
                   getMass(i));
    fVarp += getMass(i) * fA * fB * (1.0 - fC) * 
             ((1.0 - fA) * (1.0 - fA) + (1.0 - fB) * (1.0 - fB) + fC * fC);
    pMap->setValue(iX0, iY1, iZ0, 
                   pMap->getValue(iX0, iY1, iZ0) + fA * (1 - fB) * fC *
                   getMass(i));
    fVarp += getMass(i) * fA * (1.0 - fB) * fC * 
             ((1.0 - fA) * (1.0 - fA) + fB * fB + (1.0 - fC) * (1.0 - fC));
    pMap->setValue(iX1, iY0, iZ0, 
                   pMap->getValue(iX1, iY0, iZ0) + (1.0 - fA) * fB * fC * 
                   getMass(i));
    fVarp += getMass(i) * (1.0 - fA) * fB * fC * 
             (fA * fA + (1.0 - fB) * (1.0 - fB) + (1.0 - fC) * (1.0 - fC));
    pMap->setValue(iX0, iY1, iZ1, 
                   pMap->getValue(iX0, iY1, iZ1) + fA * (1 - fB) * (1.0 - fC) *
                   getMass(i));
    fVarp += getMass(i) * fA * (1.0 - fB) * (1.0 - fC) * 
             ((1.0 - fA) * (1.0 - fA) + fB * fB + fC * fC);
    pMap->setValue(iX1, iY1, iZ0, 
                   pMap->getValue(iX1, iY1, iZ0) + (1.0 - fA) * (1 - fB) * fC *
                   getMass(i));
    fVarp += getMass(i) * (1.0 - fA) * (1.0 - fB) * fC * 
             (fA * fA + fB * fB + (1.0 - fC) * (1.0 - fC));
    pMap->setValue(iX1, iY0, iZ1, 
                   pMap->getValue(iX1, iY0, iZ1) + (1.0 - fA) * fB * 
                   (1.0 - fC) * getMass(i));
    fVarp += getMass(i) * (1.0 - fA) * fB * (1.0 - fC) * 
             (fA * fA + (1.0 - fB) * (1.0 - fB) + fC * fC);
    pMap->setValue(iX1, iY1, iZ1, 
                   pMap->getValue(iX1, iY1, iZ1) + (1.0 - fA) * (1 - fB) * 
                   (1.0 - fC) * getMass(i));
    fVarp += getMass(i) * (1.0 - fA) * (1.0 - fB) * (1.0 - fC) * 
             (fA * fA + fB * fB + fC * fC);
    fSumMass += getMass(i);
  }
  fVarp /= fSumMass;

  // convolve
  volume oKernelVarp;
  oKernelVarp.createKernel(fWidth, fResolution, fVarp);
  pMap->convolve(oKernelVarp);

  // set correct position of map relative to pdb
  pMap->setOrig(fMinx - ((double)(iMargin)*fWidth), fMiny - ((double)(iMargin)*fWidth), fMinz - ((double)(iMargin)*fWidth));

  // return
  return pMap;
};

/**
 * Get the minimal coordinates of the point cloud - it will return a vector that has in each dimension the information about the minimal
 * coordinate it has found in the cloud.
 */
vec4 pointcloud::getMinCoord()
{
  double fMinX = 1.0E10;
  double fMinY = 1.0E10;
  double fMinZ = 1.0E10;

  unsigned int i;

  for (i = 0; i < size(); i++) {
    if ((*this)[i].x() < fMinX)
      fMinX = (*this)[i].x();
    if ((*this)[i].y() < fMinY)
      fMinY = (*this)[i].y();
    if ((*this)[i].z() < fMinZ)
      fMinZ = (*this)[i].z();
  }

  vec4 oVec;
  oVec.x(fMinX);
  oVec.y(fMinY);
  oVec.z(fMinZ);

  return oVec;
};

/**
 * Get the maximal coordinates of the point cloud - it will return a vector that has in each dimension the information about the maximal
 * coordinate it has found in the cloud.
 */
vec4 pointcloud::getMaxCoord()
{
  double fMaxX = -1.0E10;
  double fMaxY = -1.0E10;
  double fMaxZ = -1.0E10;

  unsigned int i;

  for (i = 0; i < size(); i++) {
    if ((*this)[i].x() > fMaxX)
      fMaxX = (*this)[i].x();
    if ((*this)[i].y() > fMaxY)
      fMaxY = (*this)[i].y();
    if ((*this)[i].z() > fMaxZ)
      fMaxZ = (*this)[i].z();
  }

  vec4 oVec;
  oVec.x(fMaxX);
  oVec.y(fMaxY);
  oVec.z(fMaxZ);

  return oVec;
};

/**
 * Get the mass of a atom number i
 */
double pointcloud::getMass(unsigned int i)
{
  return m_aMass[ i ];
};

/**
 * adjust the atomic mass based on a (simple) periodic table.
 * ToDo: Full periodic table
 */
void pointcloud::adjustMass(char *pName)
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
    if (strcmp(pName, atom_name[i]) == 0) {
      m_aMass.push_back(atom_mass[i]);
      return;
    }
  }

  m_aMass.push_back(40.0);
};

/**
 * Add a mass
 */
void pointcloud::addMass(double fMass)
{
  m_aMass.push_back(fMass);
};

/**
 * sample the object randomly and return a vector that refrects the probability distribution of the object
 */
vec4 pointcloud::sample()
{
  unsigned int iIndex = (unsigned int)((rand() / (RAND_MAX + 1.0)) * (double)(size()));
  return (*this)[iIndex];
};

///////////////////////////////////////////////////////////////////////////////
// COLUMNREADER
///////////////////////////////////////////////////////////////////////////////

/**
 * Constructor
 */
columnreader::columnreader(const char *pFilename) : m_pFile(NULL)
{
  m_pFile = fopen(pFilename, "r");
};
columnreader::~columnreader()
{
  if (m_pFile != NULL)
    fclose(m_pFile);
};

/**
 * Read next line.
 */
bool columnreader::readLine()
{
  if (m_pFile == NULL)
    return false;

  // read the complete line. return if EOF.
  if (fgets(m_pLine, MAXLINE, m_pFile) == NULL)
    return false;

  return true;
};

/**
 * Extract a string.
 * \param iStart first column (starts with 0!)
 * \param iEnd last column (this column is still read!)
 * \return string - please delete[] after use!
 */
char *columnreader::extractString(unsigned int iStart, unsigned int iEnd)
{
  unsigned int i = 0;
  unsigned int iLength = iEnd - iStart;

  if (iEnd > strlen(m_pLine))
    iLength = strlen(m_pLine) - iStart - 1;

  if (iStart >= strlen(m_pLine) || iLength <= 0) {
    m_pString[0] = 0;
    return m_pString;
  }

  for (i = 0; i <= iLength; i++)
    m_pString[i] = m_pLine[i + iStart];

  m_pString[i] = 0;
  return m_pString;
};

/**
 * Extract a char.
 * \param iCol column where the char resides.
 * \return the char.
 */
char columnreader::extractChar(unsigned int iCol)
{
  if (iCol < strlen(m_pLine))
    return m_pLine[iCol];
  else
    return ' ';
};

/**
 * Extract a real32.
 * \param iStart first column (starts with 0!)
 * \param iEnd last column (this column is still read!)
 * \return Real32 value
 */
double columnreader::extractdouble(unsigned int iStart, unsigned int iEnd)
{
  double fVal = 0.0f;

  char *pString = extractString(iStart, iEnd);
  fVal = atof(pString);

  return fVal;
};

/**
 * Extract an int.
 * \param iStart first column (starts with 0!)
 * \param iEnd last column (this column is still read!)
 * \return Int value
 */
int columnreader::extractInt(unsigned int iStart, unsigned int iEnd)
{
  int iVal = 0;

  char *pString = extractString(iStart, iEnd);
  iVal = atoi(pString);

  return iVal;
};

/**
 * Get length of line
 */
unsigned int columnreader::getLength() const
{
  return strlen(m_pLine);
};

/**
 * Get the entire line
 */
char *columnreader::getLine()
{
  return m_pLine;
};

/**
 * EOF test.
 */
bool columnreader::eof()
{
  if (m_pFile != NULL)
    return (bool)(feof(m_pFile));
  else
    return true;
};

///////////////////////////////////////////////////////////////////////////////
// MATCHRESULT
///////////////////////////////////////////////////////////////////////////////


matchresult::matchresult(double fScore, mat4 oMatrix, vector<int> &aModelMatch, vector<int> &aSceneMatch)
{
  m_fScore = fScore;
  m_oMatrix = oMatrix;
  m_aModelMatch = aModelMatch;
  m_aSceneMatch = aSceneMatch;
};

double matchresult::getScore() const
{
  return m_fScore;
};
void matchresult::setScore(double fScore)
{
  m_fScore = fScore;
};
mat4 matchresult::getMatrix() const
{
  return m_oMatrix;
};
void matchresult::setMatrix(mat4 oMatrix)
{
  m_oMatrix = oMatrix;
};
vector<int> &matchresult::getModelMatch()
{
  return m_aModelMatch;
};
vector<int> &matchresult::getSceneMatch()
{
  return m_aSceneMatch;
};
vector<int> &matchresult::getMatch()
{
  if (m_aMatch.size() == 0)
    makeSorted();

  return m_aMatch;
};
unsigned int matchresult::compareMatch(matchresult &rOther)
{
  if (m_aMatch.size() == 0)
    makeSorted();

  vector<int> &rOtherM = rOther.getMatch();
  unsigned int iNum = m_aMatch.size();
  unsigned int iDiff = 0;

  if (rOtherM.size() != m_aMatch.size()) {
    MLBO << "ERROR: rOtherM.size(): " << rOtherM.size() << " m_aMatch.size(): " << m_aMatch.size() << endl;
    exit(1);
  }

  for (unsigned int i = 0; i < iNum; i++)
    if (rOtherM[i] != m_aMatch[i] && rOtherM[i] != NOMATCH && m_aMatch[i] != NOMATCH)
      iDiff++;

  return iDiff;
};

void matchresult::printMatch()
{
  if (m_aMatch.size() == 0)
    makeSorted();

  printf("(");
  for (unsigned int i = 0; i < m_aMatch.size(); i++)
    if (m_aMatch[i] != NOMATCH)
      if (i != m_aMatch.size() - 1)
        printf("%2i,", m_aMatch[i] + 1);
      else
        printf("%2i)", m_aMatch[i] + 1);
    else
      printf(" - ");
};

bool matchresult::operator<(const matchresult &rR) const
{
  return m_fScore < rR.m_fScore;
};

/**
 * Generate the sorted match
 */
void matchresult::makeSorted()
{
  if (m_aMatch.size() == 0) {
    m_aMatch = vector<int>(m_aModelMatch.size());
    for (unsigned int i = 0; i < m_aModelMatch.size(); i++)
      m_aMatch[m_aModelMatch[i]] = m_aSceneMatch[i];
  }
};

///////////////////////////////////////////////////////////////////////////////
// VOLUME
///////////////////////////////////////////////////////////////////////////////

/**
 * Constructor
 */
volume::volume(const char *pFilename) :
  m_fWidth(0.0),
  m_iSizeX(0),
  m_iSizeY(0),
  m_iSizeZ(0),
  m_fOrigX(0.0),
  m_fOrigY(0.0),
  m_fOrigZ(0.0),
  m_pPhi(0),
  m_bChanged(true),
  m_fCutoff(-1.0E30)
{
  if (pFilename != 0) {
    loadVolume(pFilename);
  }
};
volume::volume(unsigned int iSizeX, unsigned int iSizeY, unsigned int iSizeZ) :
  m_fWidth(0.0),
  m_iSizeX(iSizeX),
  m_iSizeY(iSizeY),
  m_iSizeZ(iSizeZ),
  m_fOrigX(0.0),
  m_fOrigY(0.0),
  m_fOrigZ(0.0),
  m_pPhi(0),
  m_bChanged(true),
  m_fCutoff(-1.0E30)
{
  m_pPhi = new double[ m_iSizeX * m_iSizeY * m_iSizeZ ];
  for (unsigned int i = 0; i < m_iSizeX * m_iSizeY * m_iSizeZ; i++)
    m_pPhi[i] = 0.0;
};
volume::~volume()
{
  if (m_pPhi != 0)
    free(m_pPhi);
}

/**
 * read density file
 * \param pFilename pointer to array of char with the filename
 */
void volume::loadVolume(const char *pFilename)
{
  read_vol((char *)pFilename, &m_fWidth, &m_fOrigX, &m_fOrigY, &m_fOrigZ, 
           &m_iSizeX, &m_iSizeY, &m_iSizeZ, &m_pPhi);
  m_bChanged = true;
};


/**
 * Get the value at a position inside the volume
 * Does boundary checking!
 * \param iX x coordinate
 * \param iY y coordinate
 * \param iZ z coordinate
 * \return value
 */
double volume::getValue(unsigned int iX, unsigned int iY, unsigned int iZ) const
{
  if (iX < m_iSizeX && iY < m_iSizeY && iZ < m_iSizeZ)
    return m_pPhi[iX + (iY * m_iSizeX) + (iZ * m_iSizeX * m_iSizeY)];
  else
    return 0.0;
};

/**
 * Changes one voxel value.
 * Does boundary checking!
 * \param iX x coordinate
 * \param iY y coordinate
 * \param iZ z coordinate
 * \param fValue new value
 */
void volume::setValue(unsigned int iX, unsigned int iY, unsigned int iZ, double fValue)
{
  if (iX < m_iSizeX && iY < m_iSizeY && iZ < m_iSizeZ)
    m_pPhi[iX + iY * m_iSizeX + iZ * m_iSizeX * m_iSizeY] = fValue;

  m_bChanged = true;
};
/**
 * Changes all voxel values.
 * \param fValue new value
 */
void volume::setValue(double fValue)
{
  unsigned int iNum = m_iSizeX * m_iSizeY * m_iSizeZ;
  unsigned int i;
  for (i = 0; i < iNum; i++)
    m_pPhi[i] = fValue;

  m_bChanged = true;
};

/**
 * get the value at a position inside the volume
 * \param iCount counter
 * \return value
 */
double volume::getValue(unsigned int iCount) const
{
  if (iCount < m_iSizeX * m_iSizeY * m_iSizeZ)
    return m_pPhi[iCount];
  else
    return 0.0;
};

/**
 * Calculate correlation with other volume object
 * \param rVol reference to other volume object
 */
double volume::correlation(volume &rVol)
{
  // is there actually data in the two objects?
  if (size() == 0 || rVol.size() == 0)
    return -1.0e9;

  // calculate correlation
  int iX, iY, iZ;
  double fCorrelation = 0.0;
  double fOrgX = this->getOrigX() - rVol.getOrigX();
  double fOrgY = this->getOrigY() - rVol.getOrigY();
  double fOrgZ = this->getOrigZ() - rVol.getOrigZ();
  // use of (int)floor to make the rounding correct for negative numbers
  int iOrgX = (int)floor(fOrgX / m_fWidth + 0.5);
  int iOrgY = (int)floor(fOrgY / m_fWidth + 0.5);
  int iOrgZ = (int)floor(fOrgZ / m_fWidth + 0.5);
  double fValueA, fValueB;
  double fValueA_Sum = 0.0,  fValueB_Sum = 0.0, fValueA_Sq_Sum = 0.0, fValueB_Sq_Sum = 0.0, fCoSum = 0.0;
  unsigned int iAddedVoxels = 0;

  for (iZ = 0; iZ < (int)(m_iSizeZ); iZ++)
    for (iY = 0; iY < (int)(m_iSizeY); iY++)
      for (iX = 0; iX < (int)(m_iSizeX); iX++)
        if (iX + iOrgX >= 0 && iY + iOrgY >= 0 && iZ + iOrgZ >= 0 && 
            iX + iOrgX < (int)(rVol.getSizeX()) && 
            iY + iOrgY < (int)(rVol.getSizeY()) && iZ + iOrgZ < 
            (int)(rVol.getSizeZ())) {
          fValueA = (double)(this->getValue(iX, iY, iZ));
          fValueB = (double)(rVol.getValue(iX + iOrgX, iY + iOrgY, iZ + iOrgZ));

          if (fValueA != 0.0) {
            fCoSum    += fValueA * fValueB;
            fValueA_Sq_Sum  += fValueA * fValueA;
            fValueB_Sq_Sum  += fValueB * fValueB;
            fValueA_Sum   += fValueA;
            fValueB_Sum   += fValueB;
            iAddedVoxels++;
          }
        }

  double fAddedVoxels = (double)(iAddedVoxels);
  fCorrelation = (fAddedVoxels * fCoSum - fValueA_Sum * fValueB_Sum) / 
                 sqrt((fAddedVoxels * fValueA_Sq_Sum - fValueA_Sum * 
                       fValueA_Sum) * (fAddedVoxels * fValueB_Sq_Sum - 
                                       fValueB_Sum * fValueB_Sum));

  return fCorrelation;

};

/**
 * Calculate correlation with pointcloud object
 * \param rPDB reference to pointcloud object
 */
double volume::correlation(pointcloud &rPDB, double fResolution)
{
  volume *pVol = rPDB.blur(m_fWidth, fResolution);

  double fCorr = pVol->correlation(*this);

  delete pVol;

  return fCorr;
};

/**
 * Get the x size of the volume.
 * \return x size of the volume
 */
unsigned int volume::getSizeX() const
{
  return m_iSizeX;
};
/**
 * Get the y size of the volume.
 * \return y size of the volume
 */
unsigned int volume::getSizeY() const
{
  return m_iSizeY;
};
/**
 * Get the z size of the volume.
 * \return z size of the volume
 */
unsigned int volume::getSizeZ() const
{
  return m_iSizeZ;
};
/**
 * Get the number of voxels of the volume.
 * \return number of voxels of the volume
 */
unsigned int volume::size() const
{
  return m_iSizeX * m_iSizeY * m_iSizeZ;
};

/**
 * get position of first voxel
 */
double volume::getOrigX() const
{
  return m_fOrigX;
};
/**
 * get position of first voxel
 */
double volume::getOrigY() const
{
  return m_fOrigY;
};
/**
 * get position of first voxel
 */
double volume::getOrigZ() const
{
  return m_fOrigZ;
};

/**
 * Set voxel width of the map
 */
void volume::setWidth(double fWidth)
{
  m_fWidth = fWidth;
};

/**
 * set position of first voxel
 */
void volume::setOrig(double fOrigX, double fOrigY, double fOrigZ)
{
  m_fOrigX = fOrigX;
  m_fOrigY = fOrigY;
  m_fOrigZ = fOrigZ;
};

/**
 * Create a Gaussian blurring kernel volume (Situs scheme)
 * \param fWidth the voxel width of the target map one wants to convolve with the kernel
 * \param fResolution the target resolution
 * \param fVarp variance of map (if 0 no correction for lattice interpolation smoothing effects = default)
 */
void volume::createKernel(double fWidth, double fResolution, double fVarp)
{
  double fSigma1 = fResolution / 2.0;
  double fKmsd = fSigma1 * fSigma1 / (fWidth * fWidth);

  fResolution = 2.0 * sqrt((fSigma1 * fSigma1) + (fVarp * fWidth * fWidth));

  double fVarmap = fKmsd;
  fVarmap -= fVarp;

  if (fVarmap < 0.0) {
    MLBO << "Error: lattice smoothing exceeds kernel size" << endl;
    exit(1);
  }

  // sigma-1D
  double fSigmaMap = sqrt(fVarmap / 3.0);

  // truncate at 3 sigma-1D == sqrt(3) sigma-3D
  unsigned int iSizeH = (int) ceil(3.0 * fSigmaMap);
  unsigned int iSize = 2 * iSizeH + 1;

  m_pPhi = new double[ iSize * iSize * iSize ];
  m_iSizeX = iSize;
  m_iSizeY = iSize;
  m_iSizeZ = iSize;
  setValue(0.0);

  // kernel width
  setWidth(fWidth);

  // write Gaussian within 3 sigma-1D to map
  double fBValue = -1.0 / (2.0 * fSigmaMap * fSigmaMap);
  double fCValue = 9.0 * fSigmaMap * fSigmaMap;
  double fScale = 0;
  double fDSqu;

  unsigned int iX, iY, iZ;
  for (iZ = 0; iZ < iSize; iZ++)
    for (iY = 0; iY < iSize; iY++)
      for (iX = 0; iX < iSize; iX++) {

        fDSqu = (iX - iSizeH) * (iX - iSizeH) +
                (iY - iSizeH) * (iY - iSizeH) +
                (iZ - iSizeH) * (iZ - iSizeH);

        if (fDSqu <= fCValue)
          setValue(iX, iY, iZ, exp(fDSqu * fBValue));

        fScale += getValue(iX, iY, iZ);
      }

  for (iZ = 0; iZ < iSize; iZ++)
    for (iY = 0; iY < iSize; iY++)
      for (iX = 0; iX < iSize; iX++)
        setValue(iX, iY, iZ, getValue(iX, iY, iZ) / fScale);
};

/**
 * Convolve this volume with another one (typically a 3D kernel filter)
 * \param rKernel reference to kernel volume
 */
void volume::convolve(volume &rKernel)
{
  unsigned int iSize = rKernel.getSizeX();
  unsigned int iSizeS = iSize * iSize;
  unsigned int iX, iY, iZ;
  int iKX, iKY, iKZ;
  double fVal;
  int iDim = (int)((double)(rKernel.getSizeX()) * 0.5f);
  int iStart = -iDim;
  int iEnd = iDim + 1;
  unsigned int iSizeXY = m_iSizeX * m_iSizeY;

  double *pTmp = new double[ m_iSizeX * m_iSizeY * m_iSizeZ ];
  for (unsigned int i = 0; i < m_iSizeX * m_iSizeY * m_iSizeZ; i++)
    pTmp[i] = 0.0;

  for (iZ = 0; iZ < m_iSizeZ; iZ++) {
    for (iY = 0; iY < m_iSizeY; iY++) {
      for (iX = 0; iX < m_iSizeX; iX++) {
        fVal = getValue((iX) + (iY * m_iSizeX) + (iZ * m_iSizeX * m_iSizeY));

        if (fVal != 0.0) {
          for (iKZ = iStart; iKZ < iEnd; iKZ++)
            for (iKY = iStart; iKY < iEnd; iKY++)
              for (iKX = iStart; iKX < iEnd; iKX++) {
                if ((((int)(iX) + iKX) < (int)(m_iSizeX) && 
                     ((int)(iX) + iKX) >= 0) && 
                    (((int)(iY) + iKY) < (int)(m_iSizeY) && 
                     ((int)(iY) + iKY) >= 0) && 
                    (((int)(iZ) + iKZ) < (int)(m_iSizeZ) && 
                     ((int)(iZ) + iKZ) >= 0))
                  pTmp[(iX + iKX) + ((iY + iKY)*m_iSizeX) + ((iZ + iKZ)*iSizeXY) ] += 
                  (rKernel.getValue(((iKX - iStart) * iSizeS) + ((iKY - iStart) * iSize) + iKZ - iStart) * fVal);
              }
        }

      }
    }
  }

  delete m_pPhi;
  m_pPhi = pTmp;
  m_bChanged = true;
};

/**
 * Set cutoff for sampling
 */
void volume::setCutoff(double fCutoff)
{
  m_fCutoff = fCutoff;
};
/**
 * Get cutoff for sampling
 */
double volume::getCutoff() const
{
  return m_fCutoff;
};

/**
 * Calculate/update the minimum and maximum density values.
 */
void volume::calcMinMaxDensity()
{
  unsigned int iCount;
  unsigned int iNVox = m_iSizeX * m_iSizeY * m_iSizeZ;

  m_fMaxDensity = getValue(0, 0, 0);
  m_fMinDensity = getValue(0, 0, 0);

  double fValue;

  for (iCount = 0; iCount < iNVox; iCount++) {
    fValue = getValue(iCount);

    if (m_fMaxDensity < fValue)
      m_fMaxDensity = fValue;

    if (m_fMinDensity > fValue)
      m_fMinDensity = fValue;
  }

  m_bChanged = false;
};

/**
 * get the minimum density.
 * \return the minimum density
 */
double volume::getMinDensity()
{
  if (m_bChanged)
    calcMinMaxDensity();

  return m_fMinDensity;
};
/**
 * get the maximum density.
 *  \return the maximum density
 */
double volume::getMaxDensity()
{
  if (m_bChanged)
    calcMinMaxDensity();

  return m_fMaxDensity;
};

/**
 * sample the object randomly and return a vector that refrects the probability distribution of the object
 */
vec4 volume::sample()
{
  double fMax = getMaxDensity();
  double fNumV = (double)(m_iSizeX * m_iSizeY * m_iSizeZ);
  bool bFound = false;
  unsigned int iIndex;
  double fVoxel, fDensity;

  while (bFound == false) {
    iIndex = (unsigned int)((rand() / (RAND_MAX + 1.0)) * fNumV);
    fVoxel = getValue(iIndex);
    fDensity = (double)((rand() / (RAND_MAX + 1.0)) * fMax);

    if (fVoxel > fDensity && fVoxel > m_fCutoff)
      bFound = true;
  }

  vec4 oVec;

  unsigned int iSizeXY = m_iSizeX * m_iSizeY;
  unsigned int iIndz = iIndex / iSizeXY;
  iIndex -= iIndz * iSizeXY;
  unsigned int iIndy = iIndex / m_iSizeX;
  unsigned int iIndx = iIndex - iIndy * m_iSizeX;

  oVec.x(getOrigX() + (iIndx * m_fWidth));
  oVec.y(getOrigY() + (iIndy * m_fWidth));
  oVec.z(getOrigZ() + (iIndz * m_fWidth));

  return oVec;
};

///////////////////////////////////////////////////////////////////////////////
// CLUSTERING
///////////////////////////////////////////////////////////////////////////////

/**
 * train the network
 */
void clustering::train(sampled &rObject)
{
  m_iCount = 0;

  m_fLi *= this->m_aW.size();
  m_fTi *= this->m_iMaxstep;
  m_fTf *= this->m_iMaxstep;

  // train
  vec4 oData;
  unsigned int i;

  m_aDistance.clear();
  for (i = 0; i < this->m_aW.size(); i++)
    m_aDistance.push_back(_cluster_dist(i, 0.0));

  // calculate the average training time and update the busy indicator every second
  for (i = 0; i < this->m_iMaxstep; i++) {
    oData = rObject.sample();
    train(oData);
  }

  m_fTi /= this->m_iMaxstep;
  m_fTf /= this->m_iMaxstep;
  m_fLi /= this->m_aW.size();
};

/**
 * training, input a new vector from the data and the codebook vectors will get adjusted according to the TRN algo
 * \param fData the data vector
 */
void clustering::train(vec4 &fData)
{
  unsigned int i, j;

  // calculate and sort the distances
  for (i = 0; i < this->m_aW.size(); i++)
    m_aDistance[i] = _cluster_dist(i, this->m_aW[i].distance(fData));
  sort(m_aDistance.begin(), m_aDistance.end());

  double fCountMaxstep = (double)(this->m_iCount) / (double)(this->m_iMaxstep);

  // update lambda and epsilon
  double fEpsilon = m_fEi * pow(m_fEf / m_fEi, fCountMaxstep);
  double fLambda =  m_fLi * pow(m_fLf / m_fLi, fCountMaxstep);

  // update the codebook vectors
  unsigned int iIndex;
  double fFact, fTmp;

  // ok, now update..
  unsigned int iNum = this->m_aW.size();
  for (i = 0; i < iNum; i++) {
    iIndex = m_aDistance[i].getIndex();
    fFact = fEpsilon * exp(-(double)(i) / fLambda);

    for (j = 0; j < 4; j++) {
      fTmp = fData[j] - this->m_aW[iIndex][j];
      this->m_aW[iIndex][j] = this->m_aW[iIndex][j] + (fFact * fTmp);
    }
  }

  // increase counter
  this->m_iCount++;
};

/**
 * Set the number of steps the algorithm should run
 * \param iMaxstep maximum number of iterations
 */
void clustering::setMaxstep(unsigned int iMaxstep)
{
  m_iMaxstep = iMaxstep;
};
/**
 * Get the number of steps the algorithm should run
 * \return maximum number of iterations
 */
unsigned int clustering::getMaxstep() const
{
  return m_iMaxstep;
};

/**
 * Get the current iteration number
 * \return current iteration
 */
unsigned int clustering::getStep() const
{
  return m_iCount;
};

/**
 * Add a codebook vector. Each call will add a single codebook vector and initialize it with the provided svt_vector object
 * \param fData
 */
void clustering::addVector(vec4 &fData)
{
  m_aW.addPoint(fData);
  m_aW.addMass(1.0);
};

/**
 * Get all codebook vectors
 * \return stl vector with all the svt_vector objects
 */
pointcloud clustering::getCodebook()
{
  return m_aW;
};

/**
 * Delete all codebook vectors
 */
void clustering::delCodebook()
{
  m_iCount = 0;
  m_aW.delAllPoints();
};

/**
 * Get number of codebook vectors
 * \return number of codebook vectors
 */
unsigned int clustering::getCodebookSize() const
{
  return m_aW.size();
};

/**
 * Get a single codebook vector
 * \param i index of codebook vector
 * \return reference to svt_multivec object
 */
vec4 &clustering::getVector(int i)
{
  return m_aW[i];
};

/**
 * Get average variability
 */
double clustering::getVariability()
{
  return variability;
};

/**
 * Output the codebook vectors as a PDB file
 * \param pFname pointer to char with the filename
 */
void clustering::writePDB(char *pFname)
{
  FILE *pFile = fopen(pFname, "w");

  if (pFile == NULL)
    exit(1);

  // now transform the model and output the vectors as pdb file
  vec4 oVec(2);
  for (unsigned int i = 0; i < m_aW.size(); i++) {
    fprintf(pFile, "ATOM  %5i %-4s%c%3s %c%4i%c   %8.*f%8.*f%8.*f%6.2f%6.2f %3s  %4s%2s%2s\n",
            0,
            "QPDB",
            'Q',
            "001",
            'A',
            0,
            'Q',
            coord_precision(m_aW[i].x()), m_aW[i].x(),
            coord_precision(m_aW[i].y()), m_aW[i].y(),
            coord_precision(m_aW[i].z()), m_aW[i].z(),
            0.0f,
            0.0f,
            "QPDB",
            "QPDB",
            "QPDB",
            "QPDB"
           );
  }

  fclose(pFile);
};

/**
 * Cluster a data set - the data set must be derived of sampled (e.g. volume).
 * \param iCV number of codebook vectors the routine should use for the clustering
 * \param iRuns how many runs of the clustering should be clustered together?
 * \param rPDB the object derived from svt_sampled that should analyzed.
 */
pointcloud clustering::cluster(unsigned int iCV, unsigned int iRuns, sampled &rPDB)
{
  vector< pointcloud > oVectors;
  pointcloud oCVMAP;
  unsigned int i, j, k;
  double fMaxClusterVar = 0.0;
  double *aClusterVar = new double[iCV];

  //MLBO << "  Runs: ";

  // run the clustering algo
  for (k = 0; k < iRuns; k++) {
    delCodebook();

    for (i = 0; i < iCV; i++) {
      vec4 oVec = rPDB.sample();
      addVector(oVec);
    }

    train(rPDB);

    oCVMAP = getCodebook();
    oVectors.push_back(oCVMAP);

    //cout << ".";
    //fflush( 0 );
  }
  //cout << endl;
  double fAvgClusterVar = 0.0;

  oCVMAP = oVectors[(unsigned int)(((rand() / (RAND_MAX + 1.0)) * (double)(iRuns))) ];
  // how are the clusters populated?
  vector< vector< unsigned int > > iClusterPop;
  // where did every cv ended up being?
  vector< vector< unsigned int > > aClusterCV;

  double fDiff = 10.0;
  double fOldDiff = 0.0;
  unsigned int iR;
  unsigned int iClusterSize;
  vec4 oVec;
  unsigned int iMinInd = 0;

  // compute cluster centers
  while (fabs(fDiff - fOldDiff) > 1.0E-20) {
    // initialize the cluster population vector
    iClusterPop.clear();
    aClusterCV.clear();
    for (j = 0; j < iCV; j++) {
      vector< unsigned int > iCluster;
      iClusterPop.push_back(iCluster);
    }

    // add every codebook vector in every run to a cluster
    for (i = 0; i < iRuns; i++) {
      vector< unsigned int > aCluster;
      aClusterCV.push_back(aCluster);

      for (j = 0; j < iCV; j++) {
        iMinInd = oCVMAP.nearestNeighbor(oVectors[i][j]);

        iClusterPop[iMinInd].push_back(i * iCV + j);
        aClusterCV[i].push_back(iMinInd);
      }
    }

    // compute the average position = cluster center for every cluster
    fOldDiff = fDiff;
    fDiff = 0.0;
    for (j = 0; j < iCV; j++) {
      iClusterSize = iClusterPop[j].size();

      oVec.set(0.0);
      for (i = 0; i < iClusterSize; i++) {
        iR = (unsigned int)(((double)(iClusterPop[j][i]) / (double)(iCV)));

        oVec = oVec + oVectors[iR][iClusterPop[j][i] - (iR * iCV)];
      }
      oVec = oVec / (double)(iClusterSize);

      fDiff += oVec.distance(oCVMAP[j]);

      oCVMAP[j] = oVec;
    }
  }

  // compute variability within clusters: average rmsd between cluster centers and cv's within cluster
  double fClusterVar;
  fAvgClusterVar = 0.0;
  fMaxClusterVar = 0.0;
  double fDist = 0.0;

  for (j = 0; j < iCV; j++) {
    fClusterVar = 0.0;

    iClusterSize = iClusterPop[j].size();

    if (iClusterSize == 0) {
      cout << endl;
      MLBO << "Warning! Cluster does not contain any data points!" << endl;
    }

    // calculate average
    for (i = 0; i < iClusterSize; i++) {
      iR = (unsigned int)(((double)(iClusterPop[j][i]) / (double)(iCV)));
      fDist = oCVMAP[j].distanceSq(oVectors[iR][iClusterPop[j][i] - (iR * iCV)]);
      fClusterVar += fDist;
    }
    fClusterVar = fClusterVar / (double)(iClusterSize);
    fClusterVar = sqrt(fClusterVar);
    aClusterVar[j] = fClusterVar;

    // some stats
    if (fClusterVar > fMaxClusterVar)
      fMaxClusterVar = fClusterVar;

    fAvgClusterVar += fClusterVar;

  }
  fAvgClusterVar /= iCV;

  char pInfo[256];
  char pInfo2[256];
  sprintf(pInfo, ": %5.3f A", fAvgClusterVar);
  sprintf(pInfo2, "%2i", iCV);
  MLBO << "Variability of " << pInfo2 << " codebook vectors in " << iRuns << " statistically independent runs" << pInfo << endl;

  delete [] aClusterVar;

  for (j = 0; j < iCV; j++)
    this->m_aW[j] = oCVMAP[j];

  this->variability = fAvgClusterVar;

  return this->m_aW;
};

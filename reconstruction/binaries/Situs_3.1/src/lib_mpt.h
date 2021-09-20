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

#ifndef __SITUS_LIB_MPT
#define __SITUS_LIB_MPT

///////////////////////////////////////////////////////////////////////////////
// INCLUDES AND CONSTANTS
///////////////////////////////////////////////////////////////////////////////

// constants
#define MAXLINE 2048
#define NOMATCH 100000000

// text output
#define MLBO cout << "lib_mpt> "
#define MPTO cout << "matchpt> "

// includes
#include <vector>
#include <algorithm>
using namespace std;
#include <stdio.h>

// forward declarations
class mat4;
class matchresult;
class volume;

// functions
bool isPositive(const double &value);

///////////////////////////////////////////////////////////////////////////////
// VEC4
///////////////////////////////////////////////////////////////////////////////

/**
 * Basic vec4 class
 */
class vec4
{
  private:

    double m_aData[4];

  public:

    /**
     * Constructor
     * \param fX initial x coordinate
     * \param fY initial y coordinate
     * \param fZ initial z coordinate
     * \param fW initial w coordinate
     */
    vec4(double fX, double fY, double fZ, double fW = double(1));
    vec4(double fValue = double(0), double fW = double(1));
    virtual ~vec4();

    //
    // arithmetic operators
    //
    vec4        &operator= (const vec4 &that);
    double       &operator[](unsigned i);
    const double &operator[](unsigned i) const;
    vec4        &operator+=(const vec4 &p);
    vec4        &operator+=(const double &f);
    vec4        &operator-=(const vec4 &p);
    vec4        &operator-=(const double &f);
    vec4        &operator*=(const double &f);
    vec4        &operator/=(const double &f);
    vec4         operator- (const vec4 &p);
    vec4         operator+ (const vec4 &p);
    vec4        &operator*=(const mat4 &oMat);

    /**
     * get/set methods to manipulate the coordinate
     */
    double x() const;
    void x(double value);

    double y() const;
    void y(double value);

    double z() const;
    void z(double value);

    double w() const;
    void w(double value);

    /**
     * set all three coords of the vec4 at once
     * \param fX x coord
     * \param fY y coord
     * \param fZ z coord
     */
    void set(double fX, double fY, double fZ, double fW = double(1));

    void set(double value, double fW = double(1));

    void set(const double *p);

    /**
     * get the squares length of the vec4
     * \return length^2
     */
    double lengthSq() const;

    /**
     * get the length of the vec4
     * \return length
     */
    double length() const;

    /**
     * calculate the distance between two vec4s
     * \param oVec the other vec4
     * \return distance
     */
    double distance(const vec4 &oVec) const;

    /**
     * calculate the squared distance between two vec4s
     * \param oVec the other vec4
     * \return distance
     */
    double distanceSq(const vec4 &oVec) const;

    /**
     * normalize the vec4, return *this to allow daisy chaining
     */
    vec4 &normalize();

    /**
     * Direct access to the memory
     */
    double *c_data();
    /**
     * Direct access to the memory
     */
    const double *c_data() const;

    /**
     * Print content of vector to stdout
     */
    void print();
};
vec4 operator-(const vec4 &p1, const vec4 &p2);
vec4 operator*(const mat4 &M,  const vec4 &V);
vec4 operator/(const vec4 &V,  const double &F);

///////////////////////////////////////////////////////////////////////////////
// MAT4
///////////////////////////////////////////////////////////////////////////////

/**
 * Basic 4x4 mat4
 */
class mat4
{
  private:

    double m_aData[16];

  public:

    /**
     * Constructor
     */
    mat4();
    mat4(const mat4 &rThat);
    virtual ~mat4();

    /**
     * Operators
     */
    mat4 &operator=(const mat4 &rThat);
    mat4 &operator=(const double &fValue);

    /**
     * mat4 multiplication.
     */
    //mat4& operator*=(const mat4& B);

    /**
     * sets the mat4 to the identity mat4
     */
    void loadIdentity();

    /**
     * adds a translation (from right)
     * \param fX x translation
     * \param fY y translation
     * \param fZ z translation
     */
    mat4 &translate(double fX, double fY, double fZ);

    /**
     * adds a translation
     * \param rVec reference to svt_vec44
     */
    mat4 &translate(const vec4 &rVec);

    /**
     * get the translation component
     * \return svt_vec44 object
     */
    vec4 translation() const;

    /**
     * get the x translation
     * \return x translation
     */
    double translationX() const;

    /**
     * set the x translation
     * \param fX the new x translation
     */
    mat4 &setTranslationX(double fX);

    /**
     * get the y translation
     * \return y translation
     */
    double translationY() const;

    /**
     * set the y translation
     * \param fY the new y translation
     */
    mat4 &setTranslationY(double fY);

    /**
     * get the z translation
     * \return z translation
     */
    double translationZ() const;

    /**
     * set the z translation
     * \param fZ the new z translation
     */
    mat4 &setTranslationZ(double fZ);

    /**
     * set translation
     * \param fX x component
     * \param fY y component
     * \param fZ z component
     */
    mat4 &setTranslation(double fX, double fY, double fZ);

    /**
     * set the translation component
     * \param rVec reference to svt_vec44 object
     */
    mat4 &setTranslation(const vec4 &rVec);

    /**
     * Range-unchecked Dereference Operator
     * intented to be used as matrix[iRow][iColumn]
     * this method returns a pointer to the first element of the iRow´th Row
     * the second [] is done by c
     * range-unchecked
     */
    double *operator[](unsigned iRow);

    /**
     * Range-unchecked Dereference Operator
     * intented to be used as matrix[iRow][iColumn]
     * this method returns a pointer to the first element of the iRow´th Row
     * the second [] is done by c
     * range-unchecked
     */
    const double *operator[](unsigned iRow) const;

    /**
     * Direct access to the memory
     */
    double *c_data();

    /**
     * Jacobi transformation
     * \param rEigenvectors svt_matrix object to store the eigenvectors (columns)
     * \param rEigenvalues svt_vector object to store the eigenvalues
     */
    bool jacobi(mat4 &rEigenvectors, vec4 &rEigenvalues);

    /**
     * Print content of matrix to stdout
     */
    void print();
};

mat4 operator*(const mat4 &A, const mat4 &B);

///////////////////////////////////////////////////////////////////////////////
// NEIGHBOR HELPER CLASS
///////////////////////////////////////////////////////////////////////////////

/**
 * Helper class
 */
class neighbor
{
  protected:

    double m_fScore;
    unsigned int m_iIndexA;
    unsigned int m_iIndexB;

  public:

    neighbor(double fScore, unsigned int iIndexA, unsigned int iIndexB)
    {
      m_fScore = fScore;
      m_iIndexA = iIndexA;
      m_iIndexB = iIndexB;
    };

    double getScore() const
    {
      return m_fScore;
    };
    unsigned int getIndexA() const
    {
      return m_iIndexA;
    };
    unsigned int getIndexB() const
    {
      return m_iIndexB;
    };

    bool operator<(const neighbor &rR) const
    {
      return m_fScore < rR.m_fScore;
    };
};

///////////////////////////////////////////////////////////////////////////////
// DISTANCE MATCHING HELPER CLASS
///////////////////////////////////////////////////////////////////////////////

/**
 * matching helper class
 */
class pc_dist
{
  protected:

    double m_fScore;
    unsigned int m_iIndex;

  public:

    pc_dist(double fScore, unsigned int iIndex)
    {
      m_fScore = fScore;
      m_iIndex = iIndex;
    };

    double getScore() const
    {
      return m_fScore;
    };
    unsigned int getIndex() const
    {
      return m_iIndex;
    };

    bool operator<(const pc_dist &rR) const
    {
      return m_fScore < rR.m_fScore;
    };
};

/**
 * State of the optimize function
 */
typedef struct {
  double        fBestRMSD;
  vector<int>  aBestModelMatch;
  vector<int>  aBestSceneMatch;
  mat4         oBestMatrix;
  unsigned int iNumM;
  unsigned int iNumS;
  unsigned int iLevel;
  vector< pc_dist > aDistToCOA;

} optState;

///////////////////////////////////////////////////////////////////////////////
// SAMPLED
///////////////////////////////////////////////////////////////////////////////

/**
 * Pure abstract base class of an object that can be sampled
 */
class sampled
{
  public:

    /**
     */
    virtual ~sampled()
    {
    };

    /**
     * sample the object randomly and return a vector that reflects the probability distribution of the object
     */
    virtual vec4 sample() = 0;
};

///////////////////////////////////////////////////////////////////////////////
// POINTCLOUD
///////////////////////////////////////////////////////////////////////////////

/**
 * Basic point cloud class
 */
class pointcloud : public sampled
{
  private:

    vector< vec4 >  m_oPoints;
    vector< double > m_aMass;

    // tree pruning parameter
    double m_fEpsilon;

    // tolorance distance
    double m_fLambda;

    // ranking distance
    double m_fGamma;

    // size of the matching zone
    unsigned int m_iZoneSize;

    // number of runs with different anchor points
    unsigned int m_iRuns;

    // maximum number of wildcard matches
    unsigned int m_iMaxNoMatch;

    // number of seed/anchor matches evaluated during last matching run
    unsigned int m_iNumSeeds;

    // uniqueness rmsd threshold for the solutions
    double m_fUnique;

    // next point selection scheme
    bool m_bNextPointCOA;

    // timestep information for time-varying pointclouds
    unsigned int m_iTimestep;

    // penalty for the wildcards in the matching algo
    double m_fSkipPenalty;

  public:

    pointcloud();

    /**
     * Get all points in point cloud in an stl vector.
     * \return reference to vector of svt_vector4 objects
     */
    vector< vec4 > &getPoints();

    /**
     * Add a point to point cloud.
     * \param rVec vec4 object
     */
    void addPoint(vec4 &rVec);

    /**
     * Add a point to point cloud.
     * \param fX x coord
     * \param fY y coord
     * \param fZ z coord
     */
    void addPoint(double fX, double fY, double fZ);

    /**
     * Get a point out of point cloud.
     * \param iIndex index of point
     * \return reference to vec4 object
     */
    vec4 &getPoint(unsigned int iIndex);

    /**
     * Size of point cloud.
     * \return size of pc
     */
    unsigned int size() const;

    /**
     * Dereference operator (not range checked!).
     * \param iIndex index of point in point cloud
     */
    vec4 &operator[](unsigned int iIndex);

    /**
     * Delete all points in point cloud
     */
    void delAllPoints();

    /**
     * Calculate center of atoms (COA with all masses = 1).
     * \return svt_vector4 with the COA
     */
    vec4 coa();

    /**
     * calculate RMSD between this and another PC. The PCs must be matched already! To get the minimal RMSD please use align() first!
     * \param rPC second point cloud
     */
    double rmsd(pointcloud &rPC);
    /**
     * Calculate RMSD between this and another PC. The points are matched using the nearest neighbor relationship.
     * All the nearest neighbor distances are calculated and then sorted (slow!) and only the first N-percent are used for the matching and the following RMSD calculation.
     * The idea is that approximately N-percent of the points are outliers which would increase the RMSD significantly, although the overall deviation of the two point clouds is
     * actually small.
     * \param rPC second point cloud
     * \param fPercent percentage of neighbors that should be used for the rmsd calculation
     * \param rMatrix reference to svt_matrix4 object
     */
    double rmsd_NN_Outliers(pointcloud &rPC, double fPercent, mat4 *pMatrix = NULL);

    /**
     * Write pdb file.
     * \param pFilename pointer to array of char with the filename
     * \param bAppend if true, the pdb structure as append at the end of an existing structure file
     */
    void writePDB(const char *pFilename, bool bAppend = false);
    /**
     * Load a pdb file.
     * \param pointer to array of char with the filename
     */
    void loadPDB(const char *pFilename);
    /**
     * Replace coordinates in a pdb file.
     * \param pointer to array of char with the input filename
     * \param pointer to array of char with the output filename
     */
    void replacePDB(const char *pFilename_IN, const char *pFilename_OUT);

    /**
     * calculate the average nearest neighbor distance
     * in order to reduce the complexity a random test is done and once the average stabilizes, the search is stopped.
     * \param fPrecision if average does not change by more than fPrecision between two iterations the calculation is stopped
     * \return average nearest neighbor distance
     */
    double averageNNDistance(double fPrecision);

    /**
     * Set tree pruning parameter
     * \param fEpsilon epsilon
     */
    void setEpsilon(double fEpsilon);
    /**
     * Get tree pruning parameter.
     * \param fEpsilon epsilon
     */
    double getEpsilon() const;
    /**
     * set tolorance distance for the anchor determination
     * \param fLambda lambda
     */
    void setLambda(double fLambda);
    /**
     * get tolorance distance for the anchor determination
     * \return lambda
     */
    double getLambda() const;
    /**
     * set nearest neighbor matching zone size
     * \param fGamma gamma
     */
    void setGamma(double fGamma);
    /**
     * get nearest neighbor matching zone size
     * \return gamma
     */
    double getGamma() const;

    /**
     * set the maximal size of the matching zone
     * \param iZoneSize
     */
    void setZoneSize(unsigned int iZoneSize);
    /**
     * get the maximal size of the matching zone
     * \return maximal size of matching zone
     */
    unsigned int getZoneSize() const;

    /**
     * set the maximal number of wildcard matches
     * \param iMaxNoMatch maximal number of wildcard matches
     */
    void setWildcards(unsigned int iMaxNoMatch);
    /**
     * get the maximal number of wildcard matches
     * \return maximal number of wildcard matches
     */
    unsigned int getWildcards() const;

    /**
     * set the penalty for wildcard matches
     * \param fSkipPenalty penalty for a single wildcard
     */
    void setSkipPenalty(double fSkipPenalty);
    /**
     * get the penalty for wildcard matches
     * \return penalty for a single wildcard
     */
    double getSkipPenalty() const;

    /**
     * if two solutions are very close only the one with the higher score is considered. The other solutions are removed.
     * \param minimal distance between solutions
     */
    void setUnique(double fUnique);
    /**
     * if two solutions are very close only the one with the higher score is considered. The other solutions are removed.
     * \return minimal distance between solutions
     */
    double getUnique() const;

    /**
     * set the number of runs. Each time a different set of anchor points get selected from the probe structure, beginning with the three points furthest away from each other and the COA.
     * \param iRuns number of runs
     */
    void setRuns(unsigned int iRuns);
    /**
     * get the number of runs. Each time a different set of anchor points get selected from the probe structure, beginning with the three points furthest away from each other and the COA.
     * \return number of runs
     */
    unsigned int getRuns() const;

    /**
     * Full or simple version of the matching algorithm to be performed.
     * \param bSimple if true only a reduced version is performed (faster but less accurate)
     */
    void setSimple(bool bSimple);
    /**
     * Full or simple version of the matching algorithm to be performed.
     * \return if true only a reduced version is performed (faster but less accurate)
     */
    bool getSimple() const;

    /**
     * Set the next point selection scheme to COA
     */
    void setNextPointCOA(bool bNextPointCOA);

    /**
     * match two point clouds.
     * This function will perform a full N->M matching based on an achor-based search.
     * \param rPC second point cloud
     * \param rMatch vector of unsigned ints with the indices of the points of second cloud to which the points of this cloud are matched. This vector will get erased and then filled during the matching.
     * \param rMatrices vector of svt_matrix4 objects with the transformations according to the rMatch
     */
    void match(pointcloud &rPC, vector< matchresult > &rMatch);

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
    void match(pointcloud &rPC, vector< matchresult > &rMatch, unsigned int iAnchorA, unsigned int iAnchorB, unsigned int iAnchorC);

    /**
     * internal convenience function - no real difference to rmsd - only here we don't need to create/copy a pointcloud
     * plus here we take the matching into account! If NOMATCH this point will not be used.
     */
    double calcRMSD(pointcloud &rModel, pointcloud &rScene, mat4 &rMatrix, vector<int> *pModelMatch, vector<int> *pSceneMatch);

    /**
     * optimize - subroutine for the match() procedure
     */
    mat4 optimize(optState *pState, vector<int> *pModelMatch, vector<int> *pSceneMatch, pointcloud &rModel, pointcloud &rScene, bool bInit);

    /**
     * find the nearest neighbor to a query point in the point cloud
     * \param rVec reference to svt_vector4 object - the query point
     * \return index to nearest point in point cloud
     */
    unsigned int nearestNeighbor(vec4 &rVec);

    /**
     * least-squares fit of a model point set to a scene point set
     * \param rM vector of ints with the indices of the first (model) vectors that are used
     * \param rS vector of ints with the indices of the second (scene) vectors that are used
     * \param rModel vector of svt_vector4 with the model vectors
     * \param rScene vector of svt_vector4 with the scene vectors
     * \return matrix with optimal transformation of the model into the scene
     */
    mat4 kearsley(vector< int > &rM, vector< int > &rS, pointcloud &rModel, pointcloud &rScene);

    /**
     * sort eigenvectors according to their eigenvalues
     * \param pEigenvectors pointer to svt_matrix with the eigenvectors as column vectors
     * \param pEigenvalues pointer to svt_vector with the eigenvalues
     */
    void eigensort(mat4 &rEigenvectors, vec4 &rEigenvalues);

    /**
     * blur the pdb structure and thereby create an artificial low-resolution map
     * \param fWidth voxel width of the target map
     * \param fResolution resolution of the target map
     */
    volume *blur(double fWidth, double fResolution);

    /**
     * Get the minimal coordinates of the point cloud - it will return a vector that has in each dimension the information about the minimal
     * coordinate it has found in the cloud.
     */
    vec4 getMinCoord();
    /**
     * Get the maximal coordinates of the point cloud - it will return a vector that has in each dimension the information about the maximal
     * coordinate it has found in the cloud.
     */
    vec4 getMaxCoord();

    /**
     * Get the mass of a atom number i
     */
    double getMass(unsigned int i);

    /**
     * adjust the atomic mass based on a (simple) periodic table.
     * ToDo: Full periodic table
     */
    void adjustMass(char *pName);

    /**
     * Add a mass
     */
    void addMass(double fMass);

    /**
     * sample the object randomly and return a vector that refrects the probability distribution of the object
     */
    vec4 sample();
};

/**
 * product of matrix and point cloud
 */
pointcloud operator*(const mat4 &rM, pointcloud &rPC);

///////////////////////////////////////////////////////////////////////////////
// COLUMNREADER
///////////////////////////////////////////////////////////////////////////////

/**
 * Column-based file format reader.
 */
class columnreader
{
  protected:

    char m_pLine[MAXLINE];
    FILE *m_pFile;

    char m_pString[MAXLINE];

  public:

    /**
     * Constructor
     */
    columnreader(const char *pFilename);
    ~columnreader();

    /**
     * Read next line.
     */
    bool readLine();

    /**
     * Extract a string.
     * \param iStart first column (starts with 0!)
     * \param iEnd last column (this column is still read!)
     * \return string - please delete[] after use!
     */
    char *extractString(unsigned int iStart, unsigned int iEnd);

    /**
     * Extract a char.
     * \param iCol column where the char resides.
     * \return the char.
     */
    char extractChar(unsigned int iCol);

    /**
     * Extract a real32.
     * \param iStart first column (starts with 0!)
     * \param iEnd last column (this column is still read!)
     * \return Real32 value
     */
    double extractdouble(unsigned int iStart, unsigned int iEnd);

    /**
     * Extract an int.
     * \param iStart first column (starts with 0!)
     * \param iEnd last column (this column is still read!)
     * \return Int value
     */
    int extractInt(unsigned int iStart, unsigned int iEnd);

    /**
     * Get length of line
     */
    unsigned int getLength() const;

    /**
     * Get the entire line
     */
    char *getLine();

    /**
     * EOF test.
     */
    bool eof();
};

///////////////////////////////////////////////////////////////////////////////
// MATCHRESULT
///////////////////////////////////////////////////////////////////////////////

/**
 * Helper class that encapsulates a single result of the matching process
 */
class matchresult
{
  protected:

    double m_fScore;
    mat4 m_oMatrix;
    vector<int> m_aModelMatch;
    vector<int> m_aSceneMatch;
    vector<int> m_aMatch;

  public:

    matchresult(double fScore, mat4 oMatrix, vector<int> &aModelMatch, vector<int> &aSceneMatch);

    double getScore() const;
    void setScore(double fScore);
    mat4 getMatrix() const;
    void setMatrix(mat4 oMatrix);
    vector<int> &getModelMatch();
    vector<int> &getSceneMatch();
    vector<int> &getMatch();
    unsigned int compareMatch(matchresult &rOther);

    void printMatch();
    bool operator<(const matchresult &rR) const;

  private:

    /**
     * Generate the sorted match
     */
    void makeSorted();
};

///////////////////////////////////////////////////////////////////////////////
// VOLUME
///////////////////////////////////////////////////////////////////////////////

class volume : public sampled
{
  protected:

    double m_fWidth;

    unsigned int m_iSizeX;
    unsigned int m_iSizeY;
    unsigned int m_iSizeZ;

    double m_fOrigX;
    double m_fOrigY;
    double m_fOrigZ;

    double *m_pPhi;

    bool m_bChanged;

    double m_fCutoff;

    double m_fMinDensity;
    double m_fMaxDensity;

  public:

    /**
     * Constructor
     */
    volume(const char *pFilename = 0);
    volume(unsigned int iSizeX, unsigned int iSizeY, unsigned int iSizeZ);
    ~volume();

    /**
     * read density file
     * \param pFilename pointer to array of char with the filename
     */
    void loadVolume(const char *pFilename);


    /**
     * Get the value at a position inside the volume.
     * Does boundary checking!
     * \param iX x coordinate
     * \param iY y coordinate
     * \param iZ z coordinate
     * \return value
     */
    double getValue(unsigned int iX, unsigned int iY, unsigned int iZ) const;

    /**
     * Get the value at a position inside the volume.
     * Does boundary checking!
     * \param iCount counter
     * \return value
     */
    double getValue(unsigned int iCount) const;

    /**
     * Changes one voxel value.
     * Does boundary checking!
     * \param iX x coordinate
     * \param iY y coordinate
     * \param iZ z coordinate
     * \param fValue new value
     */
    void setValue(unsigned int iX, unsigned int iY, unsigned int iZ, double fValue);
    /**
     * Changes all voxel values.
     * \param fValue new value
     */
    void setValue(double fValue);

    /**
     * Calculate correlation with other volume object
     * \param rVol reference to other volume object
     */
    double correlation(volume &rVol);
    /**
     * Calculate correlation with pointcloud object
     * \param rPDB reference to pointcloud object
     * \param fResolution resolution of target map
     */
    double correlation(pointcloud &rPDB, double fResolution);

    /**
     * Get the x size of the volume.
     * \return x size of the volume
     */
    unsigned int getSizeX() const;
    /**
     * Get the y size of the volume.
     * \return y size of the volume
     */
    unsigned int getSizeY() const;
    /**
     * Get the z size of the volume.
     * \return z size of the volume
     */
    unsigned int getSizeZ() const;
    /**
     * Get the number of voxels of the volume.
     * \return number of voxels of the volume
     */
    unsigned int size() const;

    /**
     * get position of first voxel
     */
    double getOrigX() const;
    /**
     * get position of first voxel
     */
    double getOrigY() const;
    /**
     * get position of first voxel
     */
    double getOrigZ() const;

    /**
     * Set voxel width of the map
     */
    void setWidth(double fWidth);
    /**
     * set position of first voxel
     */
    void setOrig(double fOrigX, double fOrigY, double fOrigZ);

    /**
     * Create a Gaussian blurring kernel volume (Situs scheme)
     * \param fWidth the voxel width of the target map one wants to convolve with the kernel
     * \param fResolution the target resolution
     * \param fVarp variance of map (if 0 no correction for lattice interpolation smoothing effects = default)
     */
    void createKernel(double fWidth, double fResolution, double fVarp = 0);
    /**
     * Convolve this volume with another one (typically a 3D kernel filter)
     * \param rKernel reference to kernel volume
     */
    void convolve(volume &rKernel);

    /**
     * Set cutoff for sampling
     */
    void setCutoff(double fCutoff);
    /**
     * Get cutoff for sampling
     */
    double getCutoff() const;

    /**
     * sample the object randomly and return a vector that refrects the probability distribution of the object
     */
    vec4 sample();

    /**
     * Calculate/update the minimum and maximum density values.
     */
    void calcMinMaxDensity();

    /**
     * get the minimum density.
     * \return the minimum density
     */
    double getMinDensity();
    /**
     * get the maximum density.
     *  \return the maximum density
     */
    double getMaxDensity();
};

///////////////////////////////////////////////////////////////////////////////
// CLUSTERING
///////////////////////////////////////////////////////////////////////////////

/**
 * distance sorting helper class
 */
class _cluster_dist
{
  protected:

    double m_fDist;
    int m_iIndex;

  public:

    _cluster_dist(int iIndex, double fDist)
    {
      m_iIndex = iIndex;
      m_fDist = fDist;
    };

    double getDist() const
    {
      return m_fDist;
    };
    int getIndex() const
    {
      return m_iIndex;
    };

    bool operator<(const _cluster_dist &rD) const
    {
      return m_fDist < rD.m_fDist;
    };

};

/**
 * NG algorithm
 */
class clustering
{
  protected:

    // codebook vectors
    pointcloud m_aW;

    // average variability
    double variability;

    // number of updates so far
    unsigned int m_iCount;

    // number of maximal updates
    unsigned int m_iMaxstep;

    // lambda
    double m_fLi;
    double m_fLf;
    // epsilon
    double m_fEi;
    double m_fEf;
    // time
    double m_fTf;
    double m_fTi;

    // array of distancies
    vector<_cluster_dist> m_aDistance;

  public:

    clustering() :
      m_iCount(1),
      m_iMaxstep(500000),
      m_fLi(0.2),
      m_fLf(0.02),
      m_fEi(0.1),
      m_fEf(0.001)
    {
    };
    virtual ~clustering()
    {
    };

    /**
     * Train the network
     */
    void train(sampled &rObject);
    /**
     * Training, input a new vector from the data and the codebook vectors will get adjusted according to the NG algo
     * \param fData the data vector
     */
    void train(vec4 &fData);

    /**
     * Set the number of steps the algorithm should run
     * \param iMaxstep maximum number of iterations
     */
    void setMaxstep(unsigned int iMaxstep);
    /**
     * Get the number of steps the algorithm should run
     * \return maximum number of iterations
     */
    unsigned int getMaxstep() const;

    /**
     * get the current iteration number
     * \return current iteration
     */
    unsigned int getStep() const;

    /**
     * add a codebook vector. Each call will add a single codebook vector and initialize it with the provided svt_vector object
     * \param fData
     */
    void addVector(vec4 &fData);

    /**
     * get all codebook vectors
     */
    pointcloud getCodebook();

    /**
     * delete all codebook vectors
     */
    void delCodebook();

    /**
     * get number of codebook vectors
     * \return number of codebook vectors
     */
    unsigned int getCodebookSize() const;

    /**
     * get a single codebook vector
     * \param i index of codebook vector
     * \return reference to svt_multivec object
     */
    vec4 &getVector(int i);

    /**
     * get average Variability
     */
    double getVariability();

    /**
     * Output the codebook vectors as a PDB file
     * \param pFname pointer to char with the filename
     */
    void writePDB(char *pFname);

    /**
     * Cluster a data set - the data set must be derived of svt_sampled (e.g. svt_volume).
     * \param iCV number of codebook vectors the routine should use for the clustering
     * \param iRuns how many runs of the clustering should be clustered together?
     * \param rPDB the object derived from svt_sampled that should analyzed.
     */
    pointcloud cluster(unsigned int iCV, unsigned int iRuns, sampled &rPDB);
};

#endif

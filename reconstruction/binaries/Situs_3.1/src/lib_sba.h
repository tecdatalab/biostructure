/*********************************************************************
*                           L I B _ S B A                            *
**********************************************************************
* Library is part of the Situs package URL: situs.biomachina.org     *
* (c) Mirabela Rusu, Stefan Birmanns, and Willy Wriggers, 2011- 2012 *
**********************************************************************
*                                                                    *
* Basic support routines for C++ programs derived from Sculptor SVT  *
*                                                                    *
**********************************************************************
* See legal statement for terms of distribution                      *
*********************************************************************/

#ifndef __SITUS_LIB_SVTBASICS
#define __SITUS_LIB_SVTBASICS

#include "situs.h"
#include "lib_err.h"
#include "lib_vio.h"
#include "lib_pio.h"
#include <cmath>
#include <stack>
#include <vector>
#include <math.h>
#include <stdio.h>
#include <iostream>
#include <string.h>
#include <float.h>
#include <pthread.h>
#include <algorithm>

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace std;

typedef char Int8;
typedef short Int16;
typedef int Int32;
typedef float Real32;
typedef double Real64;

#define PI 3.1415926535
#define EQ_EPS 0.00000000001

// text output
#define SVTLBBO cout << "lib_sba> "

//forward declarations

template<class T> class svt_vector3;
template<class T> class svt_vector4;
template<class T> class svt_matrix4;
class svt_point_cloud_atom;
template<class T> class svt_volume;

/**
 * PI as float constant
 */

const Real32 fPI = 3.1415926535f;

/**
 * PI as double constant
 */
const Real64 dPI = 3.141592653589793116;

/**
 * EPSILON
 */
const Real64 EPSILON = 1.0E-6;

///////////////////////////////////////////////////////////////////////////////
// SVT_SWAP
///////////////////////////////////////////////////////////////////////////////
#define svt_swap(a, b) (temp)=*(a);\
  *(a)=*(b);\
  *(b)=(temp);

#define SwapInt(a)  SwapFloat((float *)(a))
#define SwapLong(a) SwapFloat((float *)(a))

/**
 * are the CPU registers in Big-Endian-Format?
 * \return true if big endian cpu
 */
bool svt_bigEndianMachine(void);
/**
 * swap a double value
 * \param pValue pointer to the double variable
 */
void svt_swapDouble(double *pValue);
/**
 * swap a Real32 value
 * \param pValue pointer to the Real32 variable
 */
void svt_swapReal32(Real32 *pValue);
/**
 * swap a Int16 value
 * \param pValue pointer to the Int16 variable
 */
void svt_swapInt16(Int16 *value);
/**
 * swap a Int32 value
 * \param pValue pointer to the Int32 variable
 */
void svt_swapInt32(Int32 *value);


///////////////////////////////////////////////////////////////////////////////
// SVT_TIME
///////////////////////////////////////////////////////////////////////////////
/// return the time of day in millisecond (milliseconds since Jan. 1, 1970)
unsigned long svt_getToD(void);

/// return the elapsed time since starting the program in milliseconds
int svt_getElapsedTime(void);

/// sleep for the given amount of milliseconds
void svt_sleep(unsigned uiMilliSeconds);

///////////////////////////////////////////////////////////////////////////////
// SVT_RANDOM
///////////////////////////////////////////////////////////////////////////////

/**
 * seed for the random number generator
 * \param iSeed seed value
 */
void svt_sgenrand(unsigned long);

/**
 * generate a random number
 * \return the new random number
 */
Real64 svt_genrand();

/*
 ** enums of supported random number generators.
 */
typedef enum { SPLUS, RANDU, NATIVE } svt_randomGenerator;

void svt_setRandomGenerator(svt_randomGenerator);
svt_randomGenerator svt_getRandomGenerator();

void svt_ranSeed(unsigned con = 0, unsigned taus = 0);
void svt_ranSeedAll(unsigned con = 0, unsigned taus = 0);

unsigned svt_ranLargest();

/// Gleichverteilung auf [0, svt_ranLargest()]
unsigned svt_rand();

/// Gleichverteilung auf [0,1]
Real64 svt_ranUni();

/// Normalverteilung zu (0,1) nach Box-Muller-Algorithmus
Real64 svt_ranNormal();

/// Normalverteilung zu (mu, sigma^2)
Real64 svt_ranNormal(Real64 mu, Real64 sigma);

/// Cauchyverteilung
/**
 * generate a random number following a Cauchy Distribution ~ resembles with a normal distribution just wider at the tail
 * the ratio of a two normal distributed variables follow a Cauchy distribution of average 0 and sd 1
 * \return a random number following the Cauchy distribution of average 0 and standard deviation 1
 */
Real64 svt_ranCauchy();

/**
 * \return a random number following the Cauchy distribution of average mu and standard deviation sigma
 */
Real64 svt_ranCauchy(Real64 mu, Real64 sigma);

/// Weibull-Verteilung
Real64 svt_ranWeibull();

/// Hjorth-Verteilung
Real64 svt_ranHjorth();

/// Exponential-Verteilung
Real64 svt_ranExp();

/// Bernoulli-Verteilung
unsigned svt_ranBernoulli();

/// Binomial-Verteilung
unsigned svt_ranBinomial();

/// Poisson-Verteilung
unsigned svt_ranPoisson();

///////////////////////////////////////////////////////////////////////////////
// svt_cmath
///////////////////////////////////////////////////////////////////////////////

/**
 * convert degrees to radian (float version)
 */
inline Real32 deg2rad(Real32 f)
{
  return fPI * (f / 180.0f);
}

/**
 * convert degrees to radian (double version)
 */
inline Real64 deg2rad(Real64 f)
{
  return dPI * (f / 180.0);
}
inline Real64 deg2rad(int f)
{
  return dPI * (f / 180.0);
}


/**
 * convert radian to degrees (float version)
 */
inline Real32 rad2deg(Real32 f)
{
  return 180.0f * (f / fPI);
}

/**
 * convert radian to degrees (double version)
 */
inline Real64 rad2deg(Real64 f)
{
  return 180.0  * (f / dPI);
}
inline Real64 rad2deg(int f)
{
  return 180.0  * (f / dPI);
}

/**
 * swap content of to values as template function
 * (same as STL swap(), which does not seem to work everywhere)
 */
template <class T>
inline void svt_swap_values(T &a, T &b)
{
  T tmp = a;
  a = b;
  b = tmp;
}
/**
 * Check if a float value is significantly greater that 0
 * (seems like all supported compilers can now handle numeric_limits,
 *  but anyhow...)
 */
#if !defined(WIN32_MSVC) && !defined(IRIX) && !defined(GPLUSPLUS2)

#include <limits>



template <class T>
inline bool isPositive(const T &value)
{
  return value > numeric_limits<T>::epsilon();
}

#else
#  include <float.h>

inline bool isPositive(float f)
{
  return f > FLT_EPSILON;
}

inline bool isPositive(double f)
{
  return f > DBL_EPSILON;
}

inline bool isPositive(long double f)
{
  return f > LDBL_EPSILON;
}

#endif

///////////////////////////////////////////////////////////////////////////////
// svt_matrix
///////////////////////////////////////////////////////////////////////////////

/********************************************************************
 *                                                                  *
 *  file: svt_matrix.h                                              *
 *                                                                  *
 *  specification for class: svt_matrix                             *
 *  (pure inline template class)                                    *
 *                                                                  *
 *  f.delonge                                                       *
 *                                                                  *
 ********************************************************************/


///////////////////////////////////////////////////////////////////////////////
// SVT_MATRIX4
///////////////////////////////////////////////////////////////////////////////
/** An abstract matrix base class as template.
  defines a lot of usfil operators, such as +,-,*,...
  \author Frank Delonge
*/


template<class T>
class svt_matrix
{
  public:

    enum ResizeMode {Uninitialized,
                     Clear,
                     Save,
                     SaveClamp,
                     SaveCenter,
                     SaveCenterClamp,
                     Interpolate
                    };

    // All operators are now based on public methods because of
    // trouble with some compilers at friend decalaration.
    // Therefore, the matrix class does not have any friends, but
    // operators shall be documented that way.
    // (Therewith, these operators appear in the svt_matrix class documentation
    //  generated by doxygen).

#ifdef DOXYGEN

    /** Compare 2 matrixes.
     */
    friend bool operator==<T> (const svt_matrix<T> &A, const svt_matrix<T> &B);

    /** Same as !(A==B).
     */
    friend bool operator!=<T> (const svt_matrix<T> &A, const svt_matrix<T> &B);


    /** Add 2 matrixes elementwise.
    */
    friend svt_matrix<T> operator+<T> (const svt_matrix<T> &A, const svt_matrix<T> &B);

    /** Add matrix and scalar -> matrix.
      Value is added to each element of A.
    */
    friend svt_matrix<T> operator+<T> (const svt_matrix<T> &A, const T &value);


    /** Add scalar and matrix -> matrix.
      value is added to each element of A
    */
    friend svt_matrix<T> operator+<T> (const T &value, const svt_matrix<T> &B);

    /** substract 2 matrixes elementwise
    */
    friend svt_matrix<T> operator-<T> (const svt_matrix<T> &A, const svt_matrix<T> &B);



    /** Subtract matrix a scalar from matrix A -> matrix.
      Value is substracted from each element of A.
    */
    friend svt_matrix<T> operator-<T> (const svt_matrix<T> &A, const T &value);

    /** Subtract scalar and matrix -> matrix.
      Same as -(A-value).
    */
    friend svt_matrix<T> operator-<T> (const T &value, const svt_matrix<T> &B);



    /** Negates each element.
      Same as 0-A.
    */
    friend svt_matrix<T> operator-<T> (const svt_matrix<T> &A);

    /** Matrix multiplication.
      A(h1,b1)*B(h2,b2) -> matrix(h1, b2)
    */
    friend svt_matrix<T> operator*<T> (const svt_matrix<T> &A, const svt_matrix<T> &B);

    /// scalar multiplication (elementwise)
    friend svt_matrix<T> operator*<T> (const svt_matrix<T> &A, const T &value);

    /// scalar multiplication (elementwise)
    friend svt_matrix<T> operator*<T> (const T &value, const svt_matrix<T> &A);

    /// same as A*(1/value)
    friend svt_matrix<T> operator/<T> (const svt_matrix<T> &A, const T &value);

    ///compaires
    // friend svt_matrix<T> operator<<T> (const svt_matrix<T>& A, const svt_matrix<T>& B);

    /// get the maximum value of A
    friend T max<T> (const svt_matrix<T> &A);

    /// get the minimum value of A
    friend T min<T> (const svt_matrix<T> &A);

#endif


    //////////////////////////////////////////////////////////////////////
    //                                                                  //
    //                P U B L I C   M E T H O D S                       //
    //                                                                  //
    //////////////////////////////////////////////////////////////////////

  public:


    //
    // constructors
    //

    /** Create a matrix A sized iHeight x iWidth.
      If args are omitted, an empty matrix is created
      values are uninitialized.
    */
    svt_matrix(unsigned iHeight = 0, unsigned iWidth = 0);

    /** create a matrix A sized iHeight x iWidth with default values initialValue
      same as svt_matrix A(h,w); A=initialValue;
    */
    svt_matrix(unsigned iHeight, unsigned iWidth, const T &initialValue);

    /**
     * copy-constructor, e.g. svt_matrix B(A);
     */
    svt_matrix(const svt_matrix &that);


    /// destructor
    virtual ~svt_matrix();


    //
    // assign operators
    //


    /** assign a scalar value, e.g. A=1.
     * The size of A will be unchanged, value is assigned to all elements.
     */
    svt_matrix<T> &operator=(const T &value);

    /** Assign a matrix B to a matrix A, e.g. A=B.
      Size and values will become same as B.
    */
    svt_matrix<T> &operator=(const svt_matrix<T> &that);


    //
    // dereference operators
    //

    /** dereference a matrix A: A[row][column] (unchecked).
        (the second [] is evaluated by c itself.)
    */
    T *operator[](unsigned iRow);

    /** dereference a constant matrix A: A[row][column] (unchecked).
      (the second [] is evaluated by c itself.)
      In this case, the expression is not an l-value!
    */
    const T *operator[](unsigned iRow) const;

    /// dereference a matrix A: A[row][column] (range-checked)
    T &at(unsigned iRow, unsigned iColumn);

    /** Dereference a constant matrix A: A[row][column] (range-checked).
      In this case, the expression is not an l-value!
    */
    const T &at(unsigned iRow, unsigned iColumn) const;

    //
    // modification operators
    //

    /// performs a ++ on each element (prefix)
    svt_matrix<T>  &operator++();

    /// performs a ++ on each element (postfix)
    svt_matrix<T>  &operator++(int);


    /// adds value to each element
    svt_matrix<T> &operator+=(const T &value);

    /** Adds the matrix B.
     if matrixes are not size-conform, error is thrown
    */
    svt_matrix<T> &operator+=(const svt_matrix<T> &B);

    /// performs a -- on each element (prefix)
    svt_matrix<T> &operator--();

    /// performs a -- on each element (postfix)
    svt_matrix<T> &operator--(int);

    /// substracts value from each element
    svt_matrix<T> &operator-=(const T &value);

    /** substracts the matrix B.
      \error, if matrixes are not size-conform.
    */
    svt_matrix<T> &operator-=(const svt_matrix<T> &B);

    /// scalar multiplication with value
    svt_matrix<T> &operator*=(const T &value);

    /** matrix multiplication.
      \error, if matrixes are not size-conform.
    */
    svt_matrix<T> &operator*=(const svt_matrix<T> &B);

    /// same as A*=(1/value)
    svt_matrix<T> &operator/=(const T &value);

    /**
     * nxn matrix inversion by Gauss-Jordan elimination
     */
    void gaussjordan();

    /**
     * Filles the matrix with the values of the argument matrix. It places them as a block
     * \param oSrc provides the values to copy
     * \param iRow the row position where the top right cornet of oSrc will be placed
     * \param iCol the column position where the top right cornet of oSrc will be placed
     */
    void fill(svt_matrix<T> &oSrc, unsigned int iRow, unsigned int iCol);


    /**
     * min/max
     */
    //@{
    T min() const;
    T max() const;
    T minOnRow(unsigned int iIndexRow) const;
    T maxOnRow(unsigned int iIndexRow) const;
    T nthMaxOnRow(unsigned int iIndexRow, unsigned int iIndex) const;
    //@}
    //
    // query size
    //

    /// get amount of columns
    unsigned width() const;

    /// get amount of rows
    unsigned height() const;

    //
    // modify size
    //


    /** Resize to new dimensions.
     * The Parameter eResizeMode steers handling of former content: <br>
     *   Uninitialized   - The new content does not get initialized at all (default).<br>
     *   Clear           - The new content is initialized with the value tClearValue. <br>
     *   Save            - The former content will get copied into the new matrix, alway starting at the top left corner.
     *                     If the new size is greater than the old one, all new entries become tInitialValue. <br>
     *   SaveClamp       - Same as above, but new values are copies of the last old value in their row or column. <br>
     *   SaveCenter      - The old content is centered in the new matrix, new values at the borders become tInitialValue. <br>
     *   SaveCenterClamp - Same as above, but border values are clamped to first inner old value. <br>
     *   Interpolate     - Bilinear Interpolation to new dimensions. <br>
     *
     * If old size == new size, nothing happens no matter which other arguments are specified.
     * If the matrix is of fixed size (and thus uses stack memory), attempting
     * to resize the matrix will yield an error.
     */


    void resize(unsigned iHeight, unsigned iWidth, ResizeMode eResizeMode = Uninitialized, const T &tClearValue = T(0));

    /** get data as c-type vector.
     */
    T *c_data() const;



    //////////////////////////////////////////////////////////////////////
    //                                                                  //
    //             P R O T E C T E D   S T U F F                        //
    //                                                                  //
    //////////////////////////////////////////////////////////////////////

  protected:

    /** Constructor to pass in memory pointer.
     *  This constructor can be used by inherited classes
     *  that use fix stack memory.
    */
    svt_matrix(unsigned iHeight, unsigned iWidth, T *data);

    /** create a matrix A sized iHeight x iWidth with default values initialValue
      same as svt_matrix A(h,w); A=initialValue;
    */
    svt_matrix(unsigned iHeight, unsigned iWidth, const T &initialValue, T *data);

    /**
     * copy-constructor, e.g. svt_matrix B(A);
     */
    svt_matrix(const svt_matrix &that, T *data);

    /// amount of columns
    unsigned m_iWidth;

    /// amount of rows
    unsigned m_iHeight;

    /// data
    T *m_pData;


  private:

    /// init method used by contructors
    void init(unsigned iHeight, unsigned iWidth);
    bool m_bUsesDynMem;

};

/////////////////////////////////////////////////////////////////////////////




//////////////////////////////////////////////////////////////////////
//                                                                  //
//             I N L I N E   D E F I N I T I O N S                  //
//                                                                  //
//////////////////////////////////////////////////////////////////////




///////////////////////////////////////////////////////////////////////////
// init methods used by construcots                           (private)  //
//   void svt_matrx::init()                                              //
///////////////////////////////////////////////////////////////////////////
template <class T>
inline void svt_matrix<T>::init(unsigned iHeight, unsigned iWidth)
{
  if (m_bUsesDynMem) {
    if (iWidth * iHeight != 0) {
      m_iWidth = iWidth;
      m_iHeight = iHeight;
      m_pData = new T[m_iHeight * m_iWidth];
    } else {
      m_iWidth = 0;
      m_iHeight = 0;
      m_pData = NULL;
    }
  }
  return;
}



///////////////////////////////////////////////////////////////////////////
//   Public Constructors                                                 //
///////////////////////////////////////////////////////////////////////////



///////////////////////////////////////////////////////////////////////////
// Class Constructor:                                                    //
//   svt_matrix::svt_matrix(unsigned, unsigned)                          //
///////////////////////////////////////////////////////////////////////////
template <class T>
inline svt_matrix<T>::svt_matrix(unsigned iHeight, unsigned iWidth)
  : m_bUsesDynMem(true)
{
  init(iHeight, iWidth);
  return;
}


///////////////////////////////////////////////////////////////////////////
// Class Constructor:                                                    //
//   svt_matrix::svt_matrix(unsigned, unsigned, T)                       //
///////////////////////////////////////////////////////////////////////////
template <class T>
inline svt_matrix<T>::svt_matrix(unsigned iHeight, unsigned iWidth, const T &initialValue)
  : m_bUsesDynMem(true)
{
  init(iHeight, iWidth);

  for (unsigned i = 0; i < m_iHeight * m_iWidth; i++)
    m_pData[i] = initialValue;

}


///////////////////////////////////////////////////////////////////////////
// Class Copy Constructor:                                               //
//   svt_matrix::svt_matrix(const svt_matrix&)                           //
///////////////////////////////////////////////////////////////////////////
template <class T>
inline svt_matrix<T>::svt_matrix(const svt_matrix &that)
  : m_iWidth(that.m_iWidth), m_iHeight(that.m_iHeight),
    m_bUsesDynMem(true)
{
  if (m_iHeight * m_iWidth != 0) {
    m_pData = new T[m_iHeight * m_iWidth];
    memcpy(m_pData, that.m_pData, m_iHeight * m_iWidth * sizeof(T));
  } else
    m_pData = NULL;

}

///////////////////////////////////////////////////////////////////////////
//   Protected Constructors                                              //
///////////////////////////////////////////////////////////////////////////



///////////////////////////////////////////////////////////////////////////
// Class Constructor:                                                    //
//   svt_matrix::svt_matrix(unsigned, unsigned)                          //
///////////////////////////////////////////////////////////////////////////
template <class T>
inline svt_matrix<T>::svt_matrix(unsigned iHeight, unsigned iWidth, T *data)
  : m_iWidth(iWidth), m_iHeight(iHeight), m_pData(data), m_bUsesDynMem(false)
{}


///////////////////////////////////////////////////////////////////////////
// Class Constructor:                                                    //
//   svt_matrix::svt_matrix(unsigned, unsigned, T)                       //
///////////////////////////////////////////////////////////////////////////
template <class T>
inline svt_matrix<T>::svt_matrix(unsigned iHeight, unsigned iWidth, const T &initialValue, T *data)
  : m_iWidth(iWidth), m_iHeight(iHeight), m_pData(data), m_bUsesDynMem(false)
{
  for (unsigned i = 0; i < m_iHeight * m_iWidth; i++)
    m_pData[i] = initialValue;
}


///////////////////////////////////////////////////////////////////////////
// Class Copy Constructor:                                               //
//   svt_matrix::svt_matrix(const svt_matrix&)                           //
///////////////////////////////////////////////////////////////////////////
template <class T>
inline svt_matrix<T>::svt_matrix(const svt_matrix &that, T *data)
  : m_iWidth(that.m_iWidth), m_iHeight(that.m_iHeight),
    m_pData(data), m_bUsesDynMem(false)
{
  memcpy(m_pData, that.m_pData, m_iHeight * m_iWidth * sizeof(T));
}


///////////////////////////////////////////////////////////////////////////
// Class Destructor:                                      (virtual)      //
//   svt_matrix::~svt_matrix()                                           //
///////////////////////////////////////////////////////////////////////////
template <class T>
inline svt_matrix<T>::~svt_matrix()
{
  if (m_pData && m_bUsesDynMem)
    delete[] m_pData;
}


///////////////////////////////////////////////////////////////////////////
// Assign Operator                                                       //
//    svt_matrix::operator=(const T& value)                              //
// assign value to each element                                          //
///////////////////////////////////////////////////////////////////////////
template <class T>
inline svt_matrix<T> &svt_matrix<T>::operator=(const T &value)
{

  for (unsigned i = 0; i < m_iHeight * m_iWidth; i++)
    m_pData[i] = value;

  return *this;

}


///////////////////////////////////////////////////////////////////////////
// Assign Operator                                                       //
//    svt_matrix::operator=(const svt_matrix<T>& that)                   //
// assign that to this                                                   //
///////////////////////////////////////////////////////////////////////////
template <class T>
inline svt_matrix<T> &svt_matrix<T>::operator=(const svt_matrix<T> &that)
{
  //
  // self-test
  //
  if (this == &that)
    return *this;

  //
  // resize to shape of that
  // (this may throw an error...)
  //
  resize(that.height(), that.width());

  //
  // copy that's data
  //
  memcpy(m_pData, that.m_pData, m_iHeight * m_iWidth * sizeof(T));
  return *this;

}


///////////////////////////////////////////////////////////////////////////
// range-unchecked Dereference Operator                                  //
//    svt_matrix::operator[](unsigned iRow)                              //
// intented to be used as matrix[iRow][iColumn]                          //
// this method returns a pointer to the first element of the iRow´th Row //
// the second [] is done by c                                            //
// range-unchecked                                                       //
///////////////////////////////////////////////////////////////////////////
template <class T>
inline T *svt_matrix<T>::operator[](unsigned iRow)
{
  return &(m_pData[iRow * m_iWidth]);
}


///////////////////////////////////////////////////////////////////////////
// range-unchecked Dereference Operator                                  //
//    svt_matrix::operator[](unsigned iRow) const                        //
// intented to be used as matrix[iRow][iColumn]                          //
// same as above, but for const matrixes                                 //
///////////////////////////////////////////////////////////////////////////
template <class T>
inline const T *svt_matrix<T>::operator[](unsigned iRow) const
{
  return &(m_pData[iRow * m_iWidth]);
}


///////////////////////////////////////////////////////////////////////////
// range-checked element access method                                   //
//    svt_matrix::at(unsigned iRow, unsigned iColumn)                    //
// thows svt_error if iRow or iColumn is out of range                //
///////////////////////////////////////////////////////////////////////////
template <class T>
inline T &svt_matrix<T>::at(unsigned iRow, unsigned iCol)
{
  if (iRow >= m_iHeight || iCol >= m_iWidth)
    error_sba(85010,  "svt_matrix:: at():size ot of range!");
  return (*this)[iRow][iCol];

}


///////////////////////////////////////////////////////////////////////////
// range-checked element access method for constant matrixes             //
//    svt_matrix::at(unsigned iRow, unsigned iColumn)                    //
// thows svt_error if iRow or iColumn is out of range                //
///////////////////////////////////////////////////////////////////////////
template <class T>
inline const T &svt_matrix<T>::at(unsigned iRow, unsigned iCol) const
{
  if (iRow >= m_iHeight || iCol >= m_iWidth)
    error_sba(85010, "svt_matrix:: at():size ot of range!");
  return (*this)[iRow][iCol];
}


///////////////////////////////////////////////////////////////////////////
// prefix ++ operator: performs ++ on each element                       //
//    svt_matrix::operator++                                             //
///////////////////////////////////////////////////////////////////////////
template <class T>
inline svt_matrix<T> &svt_matrix<T>::operator++()
{
  for (unsigned i = 0; i < m_iHeight * m_iWidth; i++)
    ++m_pData[i];

  return *this;
}


///////////////////////////////////////////////////////////////////////////
// postfix ++ operator: performs ++ on each element                      //
//    svt_matrix::operator++                                             //
///////////////////////////////////////////////////////////////////////////
template <class T>
inline svt_matrix<T> &svt_matrix<T>::operator++(int)
{
  for (unsigned i = 0; i < m_iHeight * m_iWidth; i++)
    m_pData[i]++;

  return *this;
}


///////////////////////////////////////////////////////////////////////////
// += operator for scalar values                                         //
//    svt_matrix::operator+= (const T& value)                            //
// adds value to each matrix element                                     //
///////////////////////////////////////////////////////////////////////////
template <class T>
inline svt_matrix<T> &svt_matrix<T>::operator+=(const T &value)
{
  for (unsigned i = 0; i < m_iHeight * m_iWidth; i++)
    m_pData[i] += value;

  return *this;
}


///////////////////////////////////////////////////////////////////////////
// += operator for second matrix B                                       //
//    svt_matrix::operator+= (const svt_matrix<T>& B)                    //
// adds matix B to this matrix, throws svt_error if dimensions       //
// are not conform                                                       //
///////////////////////////////////////////////////////////////////////////
template <class T>
inline svt_matrix<T> &svt_matrix<T>::operator+=(const svt_matrix<T> &B)
{

  if (m_iWidth != B.m_iWidth || m_iHeight != B.m_iHeight)
    error_sba(85010, "svt_matrix:: matrix addition of 2 matrixes of non-conform size!");

  for (unsigned i = 0; i < m_iHeight * m_iWidth; i++)
    m_pData[i] += B.m_pData[i];

  return *this;
}


///////////////////////////////////////////////////////////////////////////
// prefix -- operator: performs -- on each element                       //
//    svt_matrix::operator--                                             //
///////////////////////////////////////////////////////////////////////////
template <class T>
inline svt_matrix<T> &svt_matrix<T>::operator--()
{
  for (unsigned i = 0; i < m_iHeight * m_iWidth; i++)
    --m_pData[i];

  return *this;
}


///////////////////////////////////////////////////////////////////////////
// postix -- operator: performs -- on each element                       //
//    svt_matrix::operator--                                             //
///////////////////////////////////////////////////////////////////////////
template <class T>
inline svt_matrix<T> &svt_matrix<T>::operator--(int)
{
  for (unsigned i = 0; i < m_iHeight * m_iWidth; i++)
    m_pData[i]--;

  return *this;
}


///////////////////////////////////////////////////////////////////////////
// -= operator for scalar values                                         //
//    svt_matrix::operator-= (const T& value)                            //
// substracts value from each matrix element                             //
///////////////////////////////////////////////////////////////////////////
template <class T>
inline svt_matrix<T> &svt_matrix<T>::operator-=(const T &value)
{
  for (unsigned i = 0; i < m_iHeight * m_iWidth; i++)
    m_pData[i] -= value;

  return *this;
}


///////////////////////////////////////////////////////////////////////////
// -= operator for second matrix B                                       //
//    svt_matrix::operator-= (const svt_matrix<T>& B)                    //
// substracts matix B from this matrix, throws svt_error if          //
// dimensions are not conform                                            //
///////////////////////////////////////////////////////////////////////////
template <class T>
inline svt_matrix<T> &svt_matrix<T>::operator-=(const svt_matrix<T> &B)
{

  if (m_iWidth != B.m_iWidth || m_iHeight != B.m_iHeight)
    error_sba(85010, "svt_matrix:: matrix subtraction of 2 matrixes of non-conform size!");

  for (unsigned i = 0; i < m_iHeight * m_iWidth; i++)
    m_pData[i] -= B.m_pData[i];

  return *this;
}


///////////////////////////////////////////////////////////////////////////
// *= operator for scalar values (scalar multiplication)                 //
//    svt_matrix::operator*= (const T& value)                            //
// multiplicates each matrix element with value                          //
///////////////////////////////////////////////////////////////////////////
template <class T>
inline svt_matrix<T> &svt_matrix<T>::operator*=(const T &value)
{
  for (unsigned i = 0; i < m_iHeight * m_iWidth; i++)
    m_pData[i] *= value;

  return *this;
}


///////////////////////////////////////////////////////////////////////////
// *= operator for a second matrix B (matrix multiplication)             //
//    svt_matrix::operator*= (const svt_matrix<T>& B)                    //
// performs matrix multiplacation this * B                               //
// note the the dimension of this usually changes:                       //
//    this(h1,b1) * B(h2,b2) -> this(h1,b2)    if b1==h2                 //
//    throws svt_error                     otherwise                 //
///////////////////////////////////////////////////////////////////////////
template <class T>
inline svt_matrix<T> &svt_matrix<T>::operator*=(const svt_matrix<T> &B)
{
  //
  // use *(matrix, matrix)-operator
  //
  *this = (*this) * B;
  return *this;
}


///////////////////////////////////////////////////////////////////////////
// /= operator for scalar value                                          //
//    svt_matrix::operator/= (const T& B)                                //
// same as this *= (1/value)                                             //
///////////////////////////////////////////////////////////////////////////
template <class T>
inline svt_matrix<T> &svt_matrix<T>::operator/=(const T &value)
{
  for (unsigned i = 0; i < m_iHeight * m_iWidth; i++)
    m_pData[i] /= value;

  return *this;
}


///////////////////////////////////////////////////////////////////////////
// query amount of columns                                               //
//    svt_matrix::width() const                                          //
///////////////////////////////////////////////////////////////////////////
template <class T>
inline unsigned svt_matrix<T>::width() const
{
  return m_iWidth;
}


///////////////////////////////////////////////////////////////////////////
// query amount of rows                                                  //
//    svt_matrix::height() const                                         //
///////////////////////////////////////////////////////////////////////////
template <class T>
inline unsigned svt_matrix<T>::height() const
{
  return m_iHeight;
}


///////////////////////////////////////////////////////////////////////////
// change dimensions                                                     //
//    svt_matrix::resize( unsigned iHeight, unsigned iWidth,             //
//                        bool bSaveData, const T& tInitialValue )       //
// change dimensions to new height and width                             //
// if sizes are same as old sizes, nothing happens                       //
// if bSaveData is true (default), old content will be resotred          //
// new elements (or all elements, if bSaveData is set to false)          //
// will become iInitialValue                                             //
///////////////////////////////////////////////////////////////////////////
template <class T>
inline void svt_matrix<T>::resize(unsigned iHeight, unsigned iWidth,
                                  ResizeMode eResizeMode, const T &tFillValue)
{

  //
  // nothing to to if we are already in desired shape
  //
  if ((iWidth == m_iWidth) && (iHeight == m_iHeight))
    return;


  //
  // is static memory is used, we cannot resize it
  // -> throw error
  //
  if (!m_bUsesDynMem)
    error_sba(85010, "svt_matrix:: attempting to resize a fixed-size matrix!");

  if (iWidth == 0 || iHeight == 0) {
    m_iWidth  = 0;
    m_iHeight = 0;

    if (m_pData)
      delete [] m_pData;
    m_pData = NULL;
    return;
  }


  if ((m_iWidth * m_iHeight == iWidth * iHeight) && (eResizeMode == Uninitialized))
    return;

  T *data_tmp = new T [iWidth * iHeight];

  if (!m_pData) eResizeMode = Uninitialized;

  switch (eResizeMode) {
      int i, j, iRef, jRef, iRef2, jRef2;
      int off_i, off_j;
      Real64 rel_x, rel_y;

    case Clear:
      for (i = 0; i < int(iWidth * iHeight); i++) data_tmp[i] = tFillValue;
      break;

    case Save:
      for (i = 0; i < int(iHeight); i++)
        for (j = 0; j < int(iWidth); j++)
          data_tmp[i * iWidth + j] =
            (i < int(m_iHeight) && j < int(m_iWidth)) ?
            m_pData[i * m_iWidth + j] : tFillValue;
      break;

    case SaveClamp:
      for (i = 0; i < int(iHeight); i++) {
        iRef = (i >= int(m_iHeight)) ? m_iHeight - 1 : i;
        for (j = 0; j < int(iWidth); j++) {
          jRef = (j >= int(m_iWidth)) ? m_iWidth - 1 : j;
          data_tmp[i * iWidth + j] = m_pData[iRef * m_iWidth + jRef];
        }
      }
      break;

    case SaveCenter:
      off_i = (int(iHeight) - int(m_iHeight)) / 2;
      off_j = (int(iWidth) - int(m_iWidth)) / 2;

      for (iRef = -off_i, i = 0; i < int(iHeight); iRef++, i++)
        for (jRef = -off_j, j = 0; j < int(iWidth); jRef++, j++)
          data_tmp[i * iWidth + j] =
            (iRef >= 0 && iRef < int(m_iHeight) &&
             jRef >= 0 && jRef < int(m_iWidth)) ?
            m_pData[iRef * m_iWidth + jRef] : tFillValue;
      break;

    case SaveCenterClamp:
      off_i = (int(iHeight) - int(m_iHeight)) / 2;
      off_j = (int(iWidth) - int(m_iWidth)) / 2;

      for (iRef = -off_i, i = 0; i < int(iHeight); iRef++, i++) {
        if (iRef < 0) iRef2 = 0;
        else if (iRef >= int(m_iHeight)) iRef2 = m_iHeight - 1;
        else iRef2 = iRef;

        for (jRef = -off_j, j = 0; j < int(iWidth); jRef++, j++) {
          if (jRef < 0) jRef2 = 0;
          else if (jRef >= int(m_iWidth)) jRef2 = m_iWidth - 1;
          else jRef2 = jRef;
          data_tmp[i * iWidth + j] = m_pData[iRef2 * m_iWidth + jRef2];
        }
      }
      break;

    case Interpolate:

      rel_x = (iWidth > 1)  ? Real64(m_iWidth - 1) / Real64(iWidth - 1) : Real64(m_iWidth - 1) / 2;
      rel_y = (iHeight > 1) ? Real64(m_iHeight - 1) / Real64(iHeight - 1) : Real64(m_iHeight - 1) / 2;

      for (i = 0; i < int(iHeight); i++) {
        Real64 doy = (iHeight > 1) ? i * rel_y : rel_y;
        iRef = (int)(doy + DBL_EPSILON);
        doy = doy - iRef;
        iRef2 = (iRef == int(m_iHeight - 1)) ? iRef : iRef + 1;
        for (j = 0; j < int(iWidth); j++) {
          Real64 dox = (iWidth > 1) ?  j * rel_x : rel_x;
          jRef = (int)(dox + DBL_EPSILON);
          dox = dox - jRef;
          jRef2 = (jRef == int(m_iWidth - 1)) ? jRef : jRef + 1;
          Real64 val = (1.0 - doy) * (1.0 - dox) * m_pData[iRef * m_iWidth + jRef ];
          val += (1.0 - doy) * (dox)    * m_pData[iRef * m_iWidth + jRef2];
          val += (doy)     * (1.0 - dox) * m_pData[iRef2 * m_iWidth + jRef ];
          val += (doy)     * (dox)    * m_pData[iRef2 * m_iWidth + jRef2];
          data_tmp[i * iWidth + j] = T(val);
        }
      }

      break;

    default:
      break;
  }

  if (m_pData)
    delete [] m_pData;

  m_pData = data_tmp;
  m_iWidth = iWidth;
  m_iHeight = iHeight;
}


template <class T>
inline T  *svt_matrix<T>::c_data() const
{
  return m_pData;
}


template <class T>
inline T svt_matrix<T>::min() const
{
  if (!m_pData)
    error_sba(85010, "min:: matrix has no elements");
  T tMin = *m_pData;
  for (unsigned i = 1; i < height()*width(); i++)
    if (m_pData[i] < tMin)
      tMin = m_pData[i];
  return tMin;

}

/**
 * min on the row
 */
template <class T>
inline T svt_matrix<T>::minOnRow(unsigned int iIndexRow) const
{
  if (!m_pData)
    error_sba(85010, "min:: matrix has no elements");
  T tMin = m_pData[iIndexRow * width()];
  for (unsigned i = iIndexRow * width() + 1; i < (1 + iIndexRow)*width(); i++) {
    if (m_pData[i] < tMin)
      tMin = m_pData[i];
  }
  return tMin;

}

/**
 * max on the row
 */
template <class T>
inline T svt_matrix<T>::maxOnRow(unsigned int iIndexRow) const
{
  if (!m_pData)
    error_sba(85010, "maxOnRow:: matrix has no elements");
  T tMax = m_pData[iIndexRow * width()];
  for (unsigned i = iIndexRow * width() + 1; i < (1 + iIndexRow)*width(); i++) {
    if (m_pData[i] > tMax)
      tMax = m_pData[i];
  }
  return tMax;

}
/**
 * The n-th greater element on the rows (index n-1 in the ascending sorted list)
 *\param iIndexRow the index of the rows
 *\param iIndex = n
 */
template <class T>
inline T svt_matrix<T>::nthMaxOnRow(unsigned int iIndexRow, unsigned int iIndex) const
{
  if (!m_pData)
    error_sba(85010, "max:: matrix has no elements");

  vector<T> oVec;

  for (unsigned int iIndexCol = iIndexRow * width(); iIndexCol < (iIndexRow + 1)*width(); iIndexCol++)
    oVec.push_back(m_pData[iIndexCol]);

  sort(oVec.begin(), oVec.end()); // ascending

  if (oVec.size() >= iIndex) {
    return oVec[iIndex - 1];
  } else {
    SVTLBBO << "Warning: the index is out of bounds. The max value is returned!" << endl;
    return oVec[oVec.size() - 1];
  }
}

template <class T>
inline T svt_matrix<T>::max() const
{
  if (!m_pData)
    error_sba(85010, "max:: matrix has no elements");
  T tMax = *m_pData;
  for (unsigned i = 1; i < height()*width(); i++)
    if (m_pData[i] > tMax)
      tMax = m_pData[i];
  return tMax;

}

/**
 * nxn matrix inversion by Gauss-Jordan elimination
 */
template <class T>
void svt_matrix<T>::gaussjordan()
{
  unsigned int i, j, k, l, m, ic = 0, ir = 0;
  int q;
  int *coli, *rowi, *pivi;
  T pivinv, currmax, dummy;
  T tmp;

  // is this a nxn matrix?
  if (width() != height())
    return;

  // memory allocation
  coli = new int[m_iWidth];
  if (coli == NULL) {
    cout << "svt_matrix> Error: Could not satisfy memory allocation request" << endl;
    exit(52910);
  }

  rowi = new int[m_iWidth];
  if (rowi == NULL) {
    cout << "svt_matrix> Error: Could not satisfy memory allocation request" << endl;
    exit(52920);
  }

  pivi = new int[m_iWidth];
  if (pivi == NULL) {
    cout << "svt_matrix> Error: Could not satisfy memory allocation request" << endl;
    exit(52930);
  }

  for (j = 0; j < m_iWidth; j++)
    pivi[j] = 0;

  for (i = 0; i < m_iWidth; i++) {
    currmax = 0;
    for (j = 0; j < m_iWidth; j++)
      if (pivi[j] != 1)
        for (k = 0; k < m_iWidth; k++) {
          if (pivi[k] == 0) {
            if (fabs((long Real64)(*this)[j][k]) >= currmax) {
              currmax = (T)fabs((long Real64)(*this)[j][k]);
              ir = j;
              ic = k;
            }
          } else if (pivi[k] > 1) {
            cout << "svt_matrix> Error: Matrix singular - pivot" << endl;
            exit(52940);
          }
        }

    ++(pivi[ic]);

    if (ir != ic)
      for (l = 0; l < m_iWidth; l++) {
        tmp = (*this)[ir][l];
        (*this)[ir][l] = (*this)[ic][l];
        (*this)[ic][l] = tmp;
      }

    rowi[i] = ir;
    coli[i] = ic;
    if ((*this)[ic][ic] == 0) {
      cout << "svt_matrix> Error: Matrix singular" << endl;
      exit(52950);
    }

    pivinv = (T)(1.0) / (*this)[ic][ic];
    (*this)[ic][ic] = 1;
    for (l = 0; l < m_iWidth; l++)
      (*this)[ic][l] *= pivinv;

    for (m = 0; m < m_iWidth; m++)
      if (m != ic) {
        dummy = (*this)[m][ic];
        (*this)[m][ic] = 0;
        for (l = 0; l < m_iWidth; l++)
          (*this)[m][l] -= (*this)[ic][l] * dummy;
      }
  }

  for (q = m_iWidth - 1; q >= 0; q--)
    if (rowi[q] != coli[q])
      for (k = 0; k < m_iWidth; k++) {
        tmp = (*this)[k][rowi[q]];
        (*this)[k][rowi[q]] = (*this)[k][coli[q]];
        (*this)[k][coli[q]] = tmp;
      }

  // deallocate memory
  delete[] coli;
  delete[] rowi;
  delete[] pivi;
}

/**
 * Filles the matrix with the values of the argument matrix. It places them as a block
 * \param oSrc provides the values to copy
 * \param iRow the row position where the top right cornet of oSrc will be placed
 * \param iCol the column position where the top right cornet of oSrc will be placed
 */
template<class T>
void svt_matrix<T>::fill(svt_matrix<T> &oSrc, unsigned int iRow, unsigned int iCol)
{
  // error check - check range
  if (oSrc.width() + iCol > m_iWidth)
    error_sba(85010, "Fill error: the Destination matrix is too small! Insuficient columns.");

  if (oSrc.height() + iRow > m_iHeight)
    error_sba(85010, "Fill error: the Destination matrix is too small! Insuficient rows.");


  //fill in the Dest matric the src as a bloct starting with position (iRow, iCol)
  for (unsigned int iIndexRow = 0; iIndexRow < oSrc.height(); iIndexRow++)
    for (unsigned int iIndexCol = 0; iIndexCol < oSrc.width(); iIndexCol++)
      if (m_pData[(iIndexRow + iRow)*m_iWidth + (iIndexCol + iCol)] != 0) {
        char pErrorMsg[80];
        sprintf(pErrorMsg, "Fill error at %d:%d the src matrix overlaps the values in Destination matrix", iRow, iCol);
        error_sba(85010, pErrorMsg);
      } else
        m_pData[(iIndexRow + iRow)*m_iWidth + (iIndexCol + iCol)] = oSrc[iIndexRow][iIndexCol];
}


/////////////////////////////////////////////////////////////////////////////



//////////////////////////////////////////////////////////////////////
//                                                                  //
//                F R I E N D   O P E R A T O R S                   //
//                                                                  //
//////////////////////////////////////////////////////////////////////



///////////////////////////////////////////////////////////////////////////
// == : compare 2 matrixes                                               //
//    bool operator== (const svt_matrix<T>& A, const svt_matrix<T>& B)   //
///////////////////////////////////////////////////////////////////////////


template <class T>
inline bool operator== (const svt_matrix<T> &A, const svt_matrix<T> &B)
{
  if (&A == &B)
    return true;

  if (A.width() != B.width()  ||  A.height() != B.height())
    return false;

  for (unsigned i = 0; i < A.width()*A.height(); i++)
    if (A.c_data()[i] != B.c_data()[i])
      return false;


  return true;
}


///////////////////////////////////////////////////////////////////////////
// != : compare 2 matrixes                       (same as !(A==B)        //
//    bool operator== (const svt_matrix<T>& A, const svt_matrix<T>& B)   //
///////////////////////////////////////////////////////////////////////////
template <class T>
inline bool operator!= (const svt_matrix<T> &A, const svt_matrix<T> &B)
{
  return !(A == B);
}


///////////////////////////////////////////////////////////////////////////
// + : add 2 matrixes                                                    //
//   svt_matrix operator+ (const svt_matrix<T>& A,                       //
//                         const svt_matrix<T>& B )                      //
// throws svt_error if dimesions are not conform                     //
///////////////////////////////////////////////////////////////////////////
template <class T>
inline svt_matrix<T> operator+ (const svt_matrix<T> &A, const svt_matrix<T> &B)
{
  svt_matrix<T> C(A);
  C += B;
  return C;
}


///////////////////////////////////////////////////////////////////////////
// + : add matrix and scalar -> matrix                                   //
//   svt_matrix operator+ (const svt_matrix<T>& A, const T& value)       //
///////////////////////////////////////////////////////////////////////////
template <class T>
inline svt_matrix<T> operator+ (const svt_matrix<T> &A, const T &value)
{
  svt_matrix<T> C(A);
  C += value;
  return C;
}


///////////////////////////////////////////////////////////////////////////
// + : add scalar and matrix -> matrix                                   //
//   svt_matrix operator+ (const T& value, const svt_matrix<T>& A)       //
///////////////////////////////////////////////////////////////////////////
template <class T>
inline svt_matrix<T> operator+ (const T &value, const svt_matrix<T> &A)
{
  return A + value;
}


///////////////////////////////////////////////////////////////////////////
// - : substract 2 matrixes                                              //
//   svt_matrix operator- (const svt_matrix<T>& A,                       //
//                         const svt_matrix<T>& B )                      //
// throws svt_error if dimesions are not conform                     //
///////////////////////////////////////////////////////////////////////////
template <class T>
inline svt_matrix<T> operator- (const svt_matrix<T> &A, const svt_matrix<T> &B)
{
  svt_matrix<T> C(A);
  C -= B;
  return C;
}


///////////////////////////////////////////////////////////////////////////
// - : substract scalar from matrix -> matrix                            //
//   svt_matrix operator- (const svt_matrix<T>& A, const T& value)       //
///////////////////////////////////////////////////////////////////////////
template <class T>
inline svt_matrix<T> operator- (const svt_matrix<T> &A, const T &value)
{
  svt_matrix<T> C(A);
  C -= value;
  return C;
}


///////////////////////////////////////////////////////////////////////////
// - : substract matrix from scalar                ( same as -(A-value) )//
//   svt_matrix operator+ (const T& value, const svt_matrix<T>& A)       //
///////////////////////////////////////////////////////////////////////////
template <class T>
inline svt_matrix<T> operator- (const T &value, const svt_matrix<T> &A)
{
  return -(A - value);
}


///////////////////////////////////////////////////////////////////////////
// - : negate each matrix element                                        //
//   svt_matrix operator- (const svt_matrix<T>& A)                       //
// throws svt_error if dimesions are not conform                     //
///////////////////////////////////////////////////////////////////////////
template <class T>
inline svt_matrix<T> operator- (const svt_matrix<T> &A)
{

  svt_matrix<T> C(A.height(), A.width());

  for (unsigned i = 0; i < C.height()*C.width(); i++)
    C.c_data()[i] = -A.c_data()[i];

  return C;
}



///////////////////////////////////////////////////////////////////////////
// * : multiplication of 2 matrices                                      //
//   svt_matrix operator* (const svt_matrix<T>& A,                       //
//                         const svt_matrix<T>& B )                      //
// performs matrix multiplacation A(h1,b1) * B(h2,b2) -> matrix(h1,b2)   //
// throws svt_error if b1 != h2                                      //
///////////////////////////////////////////////////////////////////////////

template <class T>
inline svt_matrix<T> &mat_mult
(const svt_matrix<T> &A, const svt_matrix<T> &B, svt_matrix<T> &C)
{

  T sum;

  //
  // loop over each element of the new matrix
  //
  for (unsigned iRow = 0; iRow < C.height(); iRow++)
    for (unsigned iCol = 0; iCol < C.width(); iCol++)
      // compute each element as product of A(iRow,*) * B(*,iCol)
    {
      sum = 0;
      for (unsigned i = 0; i < A.width(); i++)
        sum += A[iRow][i] * B[i][iCol];
      C[iRow][iCol] = sum;
    }

  return C;
}

template <class T>
inline svt_matrix<T> operator*(const svt_matrix<T> &A, const svt_matrix<T> &B)
{
  if (A.width() != B.height())
    error_sba(85010, "matrix multiplication of 2 matrixes of non-conform size!");

  svt_matrix<T> C(A.height(), B.width());
  return mat_mult(A, B, C);
}


///////////////////////////////////////////////////////////////////////////
// * : matrix scalar multiplication                                      //
//   svt_matrix operator* (const svt_matrix<T>& A, const T& value        //
///////////////////////////////////////////////////////////////////////////
template <class T>
inline svt_matrix<T> operator*(const svt_matrix<T> &A, const T &value)
{
  svt_matrix<T> C(A);
  C *= value;
  return C;
}


///////////////////////////////////////////////////////////////////////////
// * : matrix scalar multiplication                                      //
//   svt_matrix operator* (const T& value, const svt_matrix<T>& A        //
///////////////////////////////////////////////////////////////////////////
template <class T>
inline svt_matrix<T> operator*(const T &value, const svt_matrix<T> &A)
{
  svt_matrix<T> C(A);
  C *= value;
  return C;
}


///////////////////////////////////////////////////////////////////////////
// / : same as A * (1/value)                                             //
//   svt_matrix operator/ (const svt_matrix<T>& A, const T& value        //
///////////////////////////////////////////////////////////////////////////
template <class T>
inline svt_matrix<T> operator/ (const svt_matrix<T> &A, const T &value)
{
  svt_matrix<T> C(A);
  C /= value;
  return C;
}


///////////////////////////////////////////////////////////////////////////
// get the maximum value of A                                            //
//   const T& max(const svt_matrix<T>& A)                                //
///////////////////////////////////////////////////////////////////////////
template <class T>
inline T max(const svt_matrix<T> &A)
{
  if (!A.c_data())
    error_sba(85010, "max:: matrix has no elements");

  T tMax = A[0][0];

  for (unsigned i = 1; i < A.height() * A.width(); i++)
    if (A.c_data()[i] > tMax)
      tMax = A.c_data()[i];

  return tMax;

}


///////////////////////////////////////////////////////////////////////////
// get the minimum value of A                                            //
//   const T& min(const svt_matrix<T>& A)                                //
///////////////////////////////////////////////////////////////////////////
template <class T>
inline T min(const svt_matrix<T> &A)
{
  if (!A.c_data())
    error_sba(85010, "min:: matrix has no elements");

  T tMin = A[0][0];

  for (unsigned i = 1; i < A.height() * A.width(); i++)
    if (A.c_data()[i] < tMin)
      tMin = A.c_data()[i];

  return tMin;

}



///////////////////////////////////////////////////////////////////////////////
// SVT_MatrIX4
///////////////////////////////////////////////////////////////////////////////
class svt_scaleInfo
{
  private:
    svt_scaleInfo() {};
    static bool sm_bEverScaled;

  public:
    static bool everScaled();
    static void setEverScaled();
};

typedef svt_matrix4<float> svt_mat4real32;
typedef svt_matrix4<Real64> svt_mat4real64;

// Sum / Difference / Product
template<class T>
inline svt_matrix4<T> operator+(const svt_matrix4<T> &A, const svt_matrix4<T> &B);

template<class T>
inline svt_matrix4<T> operator+(const svt_matrix4<T> &A, const T &value);

template<class T>
inline svt_matrix4<T> operator+(const T &value, const svt_matrix4<T> &A);


template<class T>
inline svt_matrix4<T> operator-(const svt_matrix4<T> &A, const svt_matrix4<T> &B);

template<class T>
inline svt_matrix4<T> operator-(const svt_matrix4<T> &A, const T &value);

template<class T>
inline svt_matrix4<T> operator-(const T &value, const svt_matrix4<T> &A);

template<class T>
inline svt_matrix4<T> operator-(const svt_matrix4<T> &A);

template<class T>
inline svt_matrix4<T> operator*(const svt_matrix4<T> &A, const T &value);

template<class T>
inline svt_matrix4<T> operator*(const T &value, const svt_matrix4<T> &A);


template<class T>
inline svt_matrix4<T> operator/(const svt_matrix4<T> &A, const T &value);

// matrix product
template<class T>
inline svt_matrix4<T> operator*(const svt_matrix4<T> &A, const svt_matrix4<T> &B);

/**
 * A template 4x4 matrix
 *@author Stefan Birmanns, Frank Delonge
 */
template<class T>
class svt_matrix4 : public svt_matrix<T>
{
  private:

    T stack_data[16];

    bool m_bScaled;
    void inner_invert();

    void tmp_dump(T *rs, int *);

  public:

    /**
     * Constructor
     */
    svt_matrix4() : svt_matrix<T>(4, 4, stack_data), m_bScaled(false)
    {
      loadIdentity();
    }

    svt_matrix4(const svt_matrix4<T> &that)
      : svt_matrix<T>(that, stack_data), m_bScaled(that.m_bScaled) {};

    svt_matrix4(const svt_matrix<T> &that)
      : svt_matrix<T>(that, stack_data), m_bScaled(true) {};


    virtual ~svt_matrix4() { };

    /* operators */
    svt_matrix4<T> &operator=(const svt_matrix<T> &that)
    {
      m_bScaled = true;
      return static_cast<svt_matrix4<T>&>(svt_matrix<T>::operator=(that));
    }

    svt_matrix4<T> &operator=(const T &value)
    {
      m_bScaled = true;
      return static_cast<svt_matrix4<T>&>(svt_matrix<T>::operator=(value));
    }

    /* special assign operators
     * note: it is assumed that matrix4 will only be instanciated with
             float or Real64
    */
    svt_matrix4<T> &operator=(const svt_matrix4<Real64> &that);
    svt_matrix4<T> &operator=(const svt_matrix4<float> &that);


    /** matrix multiplication.
     */
    svt_matrix4<T> &operator*=(const svt_matrix4<T> &B);


    /**
     * sets the matrix to the identity matrix
     */
    void loadIdentity()
    {
      (*this) = T(0);
      (*this)[0][0] = (*this)[1][1] = (*this)[2][2] = (*this)[3][3] = T(1);
      m_bScaled = false;
    }


    /**
     * sets the content of the matrix from a string with 16 real
     */
    void setFromString(const char *str)
    {
      unsigned int i = 0;
      char tmpstr[256];
      int tmpptr = 0;
      int valcnt = 0;

      while (i < strlen(str) && valcnt < 16) {
        if (str[i] != ' ' && i != strlen(str) - 1)
          tmpstr[tmpptr++] = str[i];
        else {
          if (i == strlen(str) - 1) {
            tmpstr[tmpptr++] = str[i];
          }

          if (tmpptr != 0) {
            tmpstr[tmpptr] = 0;
            (*this)[valcnt / 4][valcnt % 4] = static_cast<T>(atof(tmpstr));
            valcnt++;
          }
          tmpptr = 0;
        }
        i++;
      }
      m_bScaled = true;
    }

    /**
     * adds a translation (from right)
     * \param fX x translation
     * \param fY y translation
     * \param fZ z translation
     */
    svt_matrix4<T> &translate(T fX, T fY, T fZ)
    {
      (*this)[0][3] += fX;
      (*this)[1][3] += fY;
      (*this)[2][3] += fZ;
      return *this;
    }

    /**
     * adds a translation
     * \param rVec reference to svt_vector4
     */
    svt_matrix4<T> &translate(const svt_vector4<T> &rVec);

    /**
     * get the translation component
     * \return svt_vector4 object
     */
    svt_vector4<T> translation() const;

    /**
     * get the x translation
     * \return x translation
     */
    T translationX() const
    {
      return (*this)[0][3];
    }

    /**
     * set the x translation
     * \param fX the new x translation
     */
    svt_matrix4<T> &setTranslationX(T fX)
    {
      (*this)[0][3] = fX;
      return *this;
    }

    /**
     * get the y translation
     * \return y translation
     */
    T translationY() const
    {
      return (*this)[1][3];
    }

    /**
     * set the y translation
     * \param fY the new y translation
     */
    svt_matrix4<T> &setTranslationY(T fY)
    {
      (*this)[1][3] = fY;
      return *this;
    }

    /**
     * get the z translation
     * \return z translation
     */
    T translationZ() const
    {
      return (*this)[2][3];
    }

    /**
     * set the z translation
     * \param fZ the new z translation
     */
    svt_matrix4<T> &setTranslationZ(T fZ)
    {
      (*this)[2][3] = fZ;
      return *this;
    }

    /**
     * set translation
     * \param fX x component
     * \param fY y component
     * \param fZ z component
     */
    svt_matrix4<T> &setTranslation(T fX, T fY, T fZ)
    {
      setTranslationX(fX);
      setTranslationY(fY);
      return setTranslationZ(fZ);
    }

    /**
     * set the translation component
     * \param rVec reference to svt_vector4 object
     */
    svt_matrix4<T> &setTranslation(const svt_vector4<T> &rVec)
    {
      return setTranslation(rVec.x(), rVec.y(), rVec.z());
    }

    /**
     * adds a rotation around a vector
     * \param rVec reference to Vector3f
     * \param fAlpha rotation angle in rad
     */
    svt_matrix4<T> &rotate(T x, T y, T z, T fAlpha)
    {
      svt_matrix4<T> oTmp; // == identity after creation

      T c  = cos(fAlpha);
      T c_ = 1 - c;
      T s  = sin(fAlpha);
      T tmp;

      oTmp[0][0] = c + x * x * c_;
      oTmp[1][1] = c + y * y * c_;
      oTmp[2][2] = c + z * z * c_;

      tmp = x * y * c_;
      oTmp[1][0] = tmp + z * s;
      oTmp[0][1] = tmp - z * s;
      tmp = x * z * c_;
      oTmp[2][0] = tmp - y * s;
      oTmp[0][2] = tmp + y * s;
      tmp = z * y * c_;
      oTmp[2][1] = tmp + x * s;
      oTmp[1][2] = tmp - x * s;

      *this = oTmp * (*this);

      return *this;
    }

    /**
     * adds a rotation around a vector
     * \param rVec reference to svt_vector4<T>
     * \param fAlpha rotation angle in rad
     */
    svt_matrix4<T> &rotate(const svt_vector4<T> &rVec, T fAlpha);

    /**
     * adds a rotation around a vector
     * \param rVec reference to svt_vector4<T>
     *             wherethe w-coordinate is taken as rotation angle (in rad)
     *
     */
    svt_matrix4<T> &rotate(const svt_vector4<T> &rVec);


    // meaning of axis values: (0,1,2) = (x, y, z)
    svt_matrix4<T>  &rotate(int iAxis, T fAlpha)
    {
      svt_matrix4<T> oTmp;  // oTmp == identity after creation
      T c = cos(fAlpha);
      T s = sin(fAlpha);

      switch (iAxis) {
        case 0:
          oTmp[1][1] = oTmp[2][2] = c;
          oTmp[2][1] = s;
          oTmp[1][2] = -s;
          break;
        case 1:
          oTmp[0][0] = oTmp[2][2] = c;
          oTmp[0][2] = s;
          oTmp[2][0] = -s;
          break;
        case 2:
          oTmp[0][0] = oTmp[1][1] = c;
          oTmp[1][0] = s;
          oTmp[0][1] = -s;
          break;
      }
      *this = oTmp * (*this);
      return *this;
    }

    /**
     * add rotation around phi, theta, psi following the Goldstein convention (rotation around z,x,z)
     * \param fPhi phi angle
     * \param fTheta theta angle
     * \param fPsi psi angle
     */
    svt_matrix4<T> &rotatePTP(T fPhi, T fTheta, T fPsi)
    {
      svt_matrix4<T> oTmp;  // oTmp == identity after creation

      Real64 fSin_psi   = sin(fPsi);
      Real64 fCos_psi   = cos(fPsi);
      Real64 fSin_theta = sin(fTheta);
      Real64 fCos_theta = cos(fTheta);
      Real64 fSin_phi   = sin(fPhi);
      Real64 fCos_phi   = cos(fPhi);

      oTmp[0][0] = fCos_psi * fCos_phi - fCos_theta * fSin_phi * fSin_psi;
      oTmp[0][1] = fCos_psi * fSin_phi + fCos_theta * fCos_phi * fSin_psi;
      oTmp[0][2] = fSin_psi * fSin_theta;
      oTmp[1][0] = -fSin_psi * fCos_phi  - fCos_theta * fSin_phi * fCos_psi;
      oTmp[1][1] = -fSin_psi * fSin_phi  + fCos_theta * fCos_phi * fCos_psi;
      oTmp[1][2] =  fCos_psi * fSin_theta;
      oTmp[2][0] =  fSin_theta * fSin_phi;
      oTmp[2][1] = -fSin_theta * fCos_phi;
      oTmp[2][2] =  fCos_theta;

      *this = oTmp * (*this);
      return *this;
    };

    //
    // scaling
    //

    /**
     * scale one axis by given factor
     */
    svt_matrix4<T> &scale(int iAxis, T factor)
    {
      m_bScaled = true;
      svt_scaleInfo::setEverScaled();
      (*this)[0][iAxis] *= factor;
      (*this)[1][iAxis] *= factor;
      (*this)[2][iAxis] *= factor;
      return *this;
    }


    /**
     * scale all 3 axis uniformly by factor value
     */
    svt_matrix4<T> &scale(T value)
    {
      return scale(0, value).scale(1, value).scale(2, value);
    }

    /**
     * scale all 3 axis differently by given factors
     */
    svt_matrix4<T> &scale(T x, T y, T z)
    {
      return scale(0, x).scale(1, y).scale(2, z);
    }


    /**
     * scale all axis by differenttly by factors of given vector
     */
    svt_matrix4<T> &scale(const svt_vector4<T> &vec);


    /**
     * compute the inverse matrix
     */
    svt_matrix4<T> &invert(void);

    /**
     * Jacobi transformation
     * \param rEigenvectors svt_matrix object to store the eigenvectors (columns)
     * \param rEigenvalues svt_vector object to store the eigenvalues
     */
    bool jacobi(svt_matrix4<T> &rEigenvectors, svt_vector4<T> &rEigenvalues);

    bool scaled() const
    {
      return m_bScaled;
    }

    /**
     * Transpose the matrix
     */
    void transpose();
};



///////////////////////////////////////////////////////////////////////////////
// svt_matrix4
///////////////////////////////////////////////////////////////////////////////



template<class T>
inline svt_matrix4<T> &svt_matrix4<T>::operator=(const svt_matrix4<Real64> &that)
{
  m_bScaled = that.scaled();
  for (int i = 0; i < 4; i++)
    for (int j = 0; j < 4; j++)
      (*this)[i][j] = (T) that[i][j];
  return *this;
}

#ifndef WIN32
template<>
inline svt_matrix4<Real64> &svt_matrix4<Real64>::operator=(const svt_matrix4<Real64> &that)
{
  m_bScaled = that.m_bScaled;
  memcpy(m_pData, that.c_data(), 16 * sizeof(Real64));
  return *this;
}
#endif

template<class T>
inline svt_matrix4<T> &svt_matrix4<T>::operator=(const svt_matrix4<float> &that)
{
  m_bScaled = that.scaled();
  for (int i = 0; i < 4; i++)
    for (int j = 0; j < 4; j++)
      (*this)[i][j] = (T) that[i][j];
  return *this;
}

#ifndef WIN32
template<>
inline svt_matrix4<float> &svt_matrix4<float>::operator=(const svt_matrix4<float> &that)
{
  m_bScaled = that.m_bScaled;
  memcpy(m_pData, that.c_data(), 16 * sizeof(float));
  return *this;
}
#endif

/**
 * adds a rotation around a vector
 * \param rVec reference to svt_vector4<T>
 * \param fAlpha rotation angle
 */
template<class T>
inline svt_matrix4<T> &svt_matrix4<T>::rotate(const svt_vector4<T> &rVec, T fAlpha)
{
  return rotate(rVec[0], rVec[1], rVec[2], fAlpha);
}

/**
 * adds a translation (from right)
 * \param rVec reference to svt_vector4
 */
template<class T>
inline svt_matrix4<T> &svt_matrix4<T>::rotate(const svt_vector4<T> &rVec)
{
  return rotate(rVec[0], rVec[1], rVec[2], rVec[3]);
}

/**
 * get the translation component
 * \return svt_vector4 object
 */
template<class T>
inline svt_vector4<T> svt_matrix4<T>::translation() const
{
  return svt_vector4<T>(translationX(), translationY(), translationZ());
}

/**
 * set the translation component
 * \param rVec reference to svt_vector4 object
 */
template<class T>
inline svt_matrix4<T> &svt_matrix4<T>::translate(const svt_vector4<T> &rVec)
{
  return translate(rVec[0], rVec[1], rVec[2]);
}


template<class T>
inline svt_matrix4<T> &scale(const svt_vector4<T> &rVec)
{
  return scale(rVec[0], rVec[1], rVec[2]);
}

/**
 * Performs the Lu decomposition for Gaussian elimination with partial
 * pivoting.
 * Changed by fd: only examine first 3 rows since the last row is always (0 0 0 1)
 * \param pivot - contains the pivot history...the original position numbers of the pivoted rows.
 */


template <class T>
void svt_matrix4<T>::inner_invert()
{

  if (!m_bScaled) {
    // just transpose if no scaling has been performed
    svt_swap_values((*this)[0][1], (*this)[1][0]);
    svt_swap_values((*this)[0][2], (*this)[2][0]);
    svt_swap_values((*this)[1][2], (*this)[2][1]);
  } else {
    // invert upper submatrix with gauss-algorithm and pivot-search
    T   rMat[9]  = {1, 0, 0, 0, 1, 0, 0, 0, 1};
    int pivot[3] = {0, 1, 2};
    int iRow, iCol, i;

    // forward substitution
    for (iRow = 0; iRow < 3; iRow++) {
      // find pivot element
      int iPR = iRow; // pivot row
      T tPivotValue = fabs((*this)[pivot[iRow]][iRow]);
      for (i = iRow + 1; i < 3; i++)
        if (fabs((*this)[pivot[i]][iRow]) > tPivotValue) {
          iPR = i;
          tPivotValue = fabs((*this)[pivot[i]][iRow]);
        }
      if (iPR != iRow)
        svt_swap_values(pivot[iRow], pivot[iPR]);
      iPR = pivot[iRow];

      // devide current pivot row + right side by pivot elemtent
      tPivotValue = (*this)[iPR][iRow];
      for (iCol = iRow; iCol < 3; iCol++)
        (*this)[iPR][iCol] /= tPivotValue;
      for (iCol = 0; iCol < 3; iCol++)
        rMat[3 * iPR + iCol] /= tPivotValue;

      // rest rows: row operations to eliminate leading entry
      for (i = iRow + 1; i < 3; i++) {
        T factor = (*this)[pivot[i]][iRow];
        for (iCol = iRow; iCol < 3; iCol++)
          (*this)[pivot[i]][iCol] -= factor * (*this)[iPR][iCol];
        for (iCol = 0; iCol < 3; iCol++)
          rMat[3 * pivot[i] + iCol] -= factor * rMat[3 * iPR + iCol];
      }
    }

    // backward substituation
    for (iRow = 2; iRow > 0; iRow--)
      for (i = iRow - 1; i >= 0; i--) {
        T factor = (*this)[pivot[i]][iRow];
        (*this)[pivot[i]][iRow] = T(0);
        for (int iCol = 0; iCol < 3; iCol++)
          rMat[3 * pivot[i] + iCol] -= factor * rMat[3 * pivot[iRow] + iCol];
      }

    // copy solution from right side into *this
    for (iRow = 0; iRow < 3; iRow++)
      for (iCol = 0; iCol < 3; iCol++)
        (*this)[iRow][iCol] = rMat[3 * pivot[iRow] + iCol];

  }

}

/**
 * compute the inverse matrix
 */
template<class T>
inline svt_matrix4<T> &svt_matrix4<T>::invert()
{
  // compute and store inverse of upper 3x3 submatrix R
  inner_invert();

  // the new translation vector can be calculated as -R^(-1)*T
  // where T is the existing translation vector
  T vec[3];
  int i;
  for (i = 0; i < 3; i++) {
    vec[i] = T(0);
    for (int j = 0; j < 3; j++)
      vec[i] -= (*this)[i][j] * (*this)[j][3];
  }
  for (i = 0; i < 3; i++)
    (*this)[i][3] = vec[i];

  return *this;
}

/**
 * Jacobi transformation
 * \param rEigenvectors svt_matrix object to store the eigenvectors (columns)
 * \param rEigenvalues svt_vector object to store the eigenvalues
 */
template<class T>
bool svt_matrix4<T>::jacobi(svt_matrix4<T> &rEigenvectors, svt_vector4<T> &rEigenvalues)
{
  Real64  fSM;
  Real64  fTheta;
  Real64  fC, fS, fT;
  Real64  fTau;
  Real64  fH, fG;
  Real64  fThresh;
  Real64  aB[4];
  Real64  aZ[4];
  int  iP, iQ, i, j;
  Real64  oA[4][4];
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

          rEigenvalues[iP] -= T(fH);
          rEigenvalues[iQ] += T(fH);

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

            rEigenvectors[j][iP] = T(fG - fS * (fH + fG * fTau));
            rEigenvectors[j][iQ] = T(fH + fS * (fG - fH * fTau));
          }
        }
        iRots++;
      }
    }
    for (iP = 0; iP < 4; iP++) {
      rEigenvalues[iP] = T(aB[iP] += aZ[iP]);
      aZ[iP] = 0;
    }
  }

  return true;
}

template<class T>
inline void svt_matrix4<T>::transpose()
{
  int iX, iY;

  T gl_mat[16];

  for (iX = 0; iX < 4; iX++)
    for (iY = 0; iY < 4; iY++)
      gl_mat[4 * iY + iX] = (Real64)((*this)[iX][iY]);

  for (iX = 0; iX < 4; iX++)
    for (iY = 0; iY < 4; iY++)
      (*this)[iY][iX] = gl_mat[4 * iY + iX];
}

template <class T>
inline svt_matrix4<T> &svt_matrix4<T>::operator*=(const svt_matrix4<T> &B)
{
  //
  // use *(matrix4, matrix4)-operator
  //
  *this = (*this) * B;
  return *this;
}





//
// binary operators

// Sum / Difference / Product
template<class T>
inline svt_matrix4<T> operator+(const svt_matrix4<T> &A, const svt_matrix4<T> &B)
{
  svt_matrix4<T> C(A);
  C += B;
  return C;
}



template<class T>
inline svt_matrix4<T> operator+(const svt_matrix4<T> &A, const T &value)
{
  svt_matrix4<T> C(A);
  C += value;
  return C;
}

template<class T>
inline svt_matrix4<T> operator+(const T &value, const svt_matrix4<T> &A)
{
  svt_matrix4<T> C(A);
  C += value;
  return C;
}


template<class T>
inline svt_matrix4<T> operator-(const svt_matrix4<T> &A, const svt_matrix4<T> &B)
{
  svt_matrix4<T> C(A);
  C -= B;
  return C;
}


template<class T>
inline svt_matrix4<T> operator-(const svt_matrix4<T> &A, const T &value)
{
  svt_matrix4<T> C(A);
  C -= value;
  return C;
}


template<class T>
inline svt_matrix4<T> operator-(const T &value, const svt_matrix4<T> &A)
{
  return -A + value;
}


template<class T>
inline svt_matrix4<T> operator-(const svt_matrix4<T> &A)
{
  svt_matrix4<T> C;
  for (unsigned i = 0; i < C.height()*C.width(); i++)
    C.c_data()[i] = -A.c_data()[i];
  return C;
}


template<class T>
inline svt_matrix4<T> operator*(const svt_matrix4<T> &A, const T &value)
{
  svt_matrix4<T> C(A);
  C *= value;
  return C;
}


template<class T>
inline svt_matrix4<T> operator*(const T &value, const svt_matrix4<T> &A)
{
  svt_matrix4<T> C(A);
  C *= value;
  return C;
}


template<class T>
inline svt_matrix4<T> operator/(const svt_matrix4<T> &A, const T &value)
{
  svt_matrix4<T> C(A);
  C /= value;
  return C;
}


// matrix product
template<class T>
inline svt_matrix4<T> operator*(const svt_matrix4<T> &A, const svt_matrix4<T> &B)
{
  svt_matrix4<T> C;
  return mat_mult(A, B, C);
}



///////////////////////////////////////////////////////////////////////////////
// SVT_Vector4
///////////////////////////////////////////////////////////////////////////////

#if !defined(NOLIMITS)
#include <limits>
#else
#include <float.h>
#endif

typedef svt_vector4<float> svt_vec4real32;
typedef svt_vector4<Real64> svt_vec4real64;

//
// binary operators
//
// NOTE: all arithmetic operators only correspond to the first 3 vector entries,
//       the w-coordinate always stays untouched


// Sum / Difference / Product

template<class T>
inline svt_vector4<T> operator+(const svt_vector4<T> &p1, const svt_vector4<T> &p2);

template<class T>
inline svt_vector4<T> operator+(const svt_vector4<T> &p, const T &f);

template<class T>
inline svt_vector4<T> operator+(const T &f, const svt_vector4<T> &p);

template<class T>
inline svt_vector4<T> operator-(const svt_vector4<T> &p1, const svt_vector4<T> &p2);

template<class T>
inline svt_vector4<T> operator-(const svt_vector4<T> &p, const T &f);

template<class T>
inline svt_vector4<T> operator-(const T &f, const svt_vector4<T> &p);

template<class T>
inline svt_vector4<T> operator-(const svt_vector4<T> &p);


template<class T>
inline svt_vector4<T> operator*(const svt_vector4<T> &p, const T &f);

template<class T>
inline svt_vector4<T> operator*(const T &f, const svt_vector4<T> &p);

template<class T>
inline svt_vector4<T> operator*(const svt_matrix4<T> &M, const svt_vector4<T> &V);

template<class T>
inline svt_vector4<T> operator*(const svt_vector4<T> &V, const svt_matrix4<T> &M);

template<class T>
inline svt_vector4<T> operator/(const svt_vector4<T> &p1, const T &f);

// Scalar Product
template<class T>
inline T operator*(const svt_vector4<T> &p1, const svt_vector4<T> &p2);


// Vector Product
template<class T>
inline svt_vector4<T>
vectorProduct(const svt_vector4<T> &p1, const svt_vector4<T> &p2);

// Vector Product of Quaternions
template<class T>
inline svt_vector4<T> quaternionVectorProduct(const svt_vector4<T> &p1, const svt_vector4<T> &p2);



///////////////////////////////////////////////////////////////////////////////
// declaration
///////////////////////////////////////////////////////////////////////////////



/** A 4 value template vector
  *@author Stefan Birmanns, Frank Delonge
  */
template<class T> class svt_vector4 : public svt_matrix<T>
{
  private:

    T stack_data[4];

  public:

    /**
     * Constructor
     * \param fX initial x coordinate
     * \param fY initial y coordinate
     * \param fZ initial z coordinate
     * \param fW initial w coordinate
     */
    svt_vector4(T fX, T fY, T fZ, T fW = T(1)) : svt_matrix<T>(4, 1, stack_data)
    {
      x(fX);
      y(fY);
      z(fZ);
      w(fW);
    }

    svt_vector4(T fValue = T(0), T fW = T(1)) : svt_matrix<T>(4, 1, stack_data)
    {
      x(fValue);
      y(fValue);
      z(fValue);
      w(fW);
    }

    svt_vector4(const svt_matrix<T> &that) : svt_matrix<T>(that, stack_data)
    {}

    virtual ~svt_vector4()
    {
    };

    svt_vector4<T> &operator=(const svt_vector4<T> &that);

    inline T &operator[](unsigned i)
    {
      return svt_matrix<T>::m_pData[i];
    }

    inline const T &operator[](unsigned i) const
    {
      return svt_matrix<T>::m_pData[i];
    }

    //
    // arithmetic operators
    // These operators need to be redefined because
    // only the first 3 coordinates shall be considered
    // all operators return *this to allow daisy chaining
    //
    svt_vector4<T> &operator+=(const svt_vector4<T> &p)
    {
      (*this)[0] += p[0];
      (*this)[1] += p[1];
      (*this)[2] += p[2];
      return *this;
    }

    svt_vector4<T> &operator+=(const T &f)
    {
      (*this)[0] += f;
      (*this)[1] += f;
      (*this)[2] += f;
      return *this;
    }

    svt_vector4<T> &operator-=(const svt_vector4<T> &p)
    {
      (*this)[0] -= p[0];
      (*this)[1] -= p[1];
      (*this)[2] -= p[2];
      return *this;
    }

    svt_vector4<T> &operator-=(const T &f)
    {
      (*this)[0] -= f;
      (*this)[1] -= f;
      (*this)[2] -= f;
      return *this;
    }

    svt_vector4<T> &operator*=(const T &f)
    {
      (*this)[0] *= f;
      (*this)[1] *= f;
      (*this)[2] *= f;
      return *this;
    }

    svt_vector4<T> &operator/=(const T &f)
    {
      (*this)[0] /= f;
      (*this)[1] /= f;
      (*this)[2] /= f;
      return *this;
    }

    /**
     * performs oMat * this, and stores result in this
     */
    svt_vector4<T> &operator*=(const svt_matrix4<T> &oMat);

    /**
     * get/set methods
     */
    inline T x() const
    {
      return (*this)[0];
    }
    inline void x(T value)
    {
      (*this)[0] = value;
    }

    inline T y() const
    {
      return (*this)[1];
    }
    inline void y(T value)
    {
      (*this)[1] = value;
    }

    inline T z() const
    {
      return (*this)[2];
    }
    inline void z(T value)
    {
      (*this)[2] = value;
    }

    inline T w() const
    {
      return (*this)[3];
    }
    inline void w(T value)
    {
      (*this)[3] = value;
    }

    /**
     * set all three coords of the vector at once
     * \param fX x coord
     * \param fY y coord
     * \param fZ z coord
     */
    inline void set(T fX, T fY, T fZ, T fW = T(1))
    {
      x(fX);
      y(fY);
      z(fZ);
      w(fW);
    }

    inline void set(T value, T fW = T(1))
    {
      *this = value;
      w(fW);
    }

    inline void set(const T *p)
    {
      memcpy(svt_matrix<T>::m_pData, p, 4 * sizeof(T));
    }

    /**
     * initialize the quaternion with a random uniform distribution rotation. Note that x^2+y^2+z^2+w^2 = 1
     * see Shoemake, Graphics Gems III.6, pp.124-132, "Uniform Random Rotations",
     * the vectors describe a ball with radius 1
     *
     */
    inline void initUniformQuaternion()
    {
      T fAlpha1, fAlpha2, fR1, fR2, fX0;

      fAlpha1 = svt_genrand() * 2.0f * PI; // no between 0 and 2PI
      fX0  = svt_genrand();
      fR1  = sqrt(1.0f - fX0);

      x(sin(fAlpha1) * fR1);
      y(cos(fAlpha1) * fR1);

      fAlpha2 = svt_genrand() * 2.0f * PI; // no between 0 and 2PI
      fR2  = sqrt(fX0);

      z(sin(fAlpha2) * fR2);
      w(cos(fAlpha2) * fR2);

    }

    /**
     * convert the unit quaternion (x^2+y^2+z^2+w^2 = 1) into one with a normal vector and an angle  between -pi and pi
     */
    inline void setQuat2Rot()
    {
      T fW = w();
      if (fabs(fW) > 1.0f) {
        error_sba(85010, "w should be between -1 and 1");
        return;
      }

      Real64 fAngle = 2.0f * acos(fW);   // value between 0 and 360
      if (fW == 1.0f) {
        x(1.0);
        y(0.0);
        z(0.0);
      } else {
        Real64 fInvSinHalfAngle = 1.0 / sin(fAngle / 2.0);

        x(x()*fInvSinHalfAngle);
        y(y()*fInvSinHalfAngle);
        z(z()*fInvSinHalfAngle);

        normalize();
      }

      //convert angle to -180, 180 interval
      /*      if (fAngle>PI)
                fAngle -= 2.0*PI;
            else if (fAngle<-PI)
                fAngle += 2.0*PI;*/

      w(fAngle);
    }


    /**
     * get the squares length of the vector
     * \return length^2
     */
    inline T lengthSq() const
    {
      return x() * x() + y() * y() + z() * z();
    }

    /**
     * get the length of the vector
     * \return length
     */
    inline T length() const
    {
      return sqrt(lengthSq());
    }

    /**
    * get the length of the quaternion
    * \return length
    */
    inline T quaternionLengthSq() const
    {
      return x() * x() + y() * y() + z() * z() + w() * w() ;
    }

    /**
     * calculate the distance between two vectors
     * \param oVec the other vector
     * \return distance
     */
    inline T distance(const svt_vector4<T> &oVec) const
    {
      T *oTmp = oVec.c_data();
      T *oThis = (*this).c_data();

      T fDiffX = oThis[0] - oTmp[0];
      T fDiffY = oThis[1] - oTmp[1];
      T fDiffZ = oThis[2] - oTmp[2];

      T fLengthSq = fDiffX * fDiffX + fDiffY * fDiffY + fDiffZ * fDiffZ;
      T fLength =  sqrt(fLengthSq);

      return fLength;
    }

    /**
     * calculate the squared distance between two vectors
     * \param oVec the other vector
     * \return distance
     */
    inline T distanceSq(const svt_vector4<T> &oVec) const
    {
      return (*this - oVec).lengthSq();
    }

    /**
     * normalize the vector, return *this to allow daisy chaining
     */
    inline svt_vector4<T> &normalize()
    {
      T fLength = length();

      if (isPositive(fLength))
        (*this) /= fLength;

      return *this;
    }

    /**
     * normalize the vector, include w-coordinate in scaling
     * return *this to allow daisy chaining
     */
    svt_vector4<T> &normalize4()
    {
      T fLength = length();
      if (isPositive(fLength)) {
        *this /= fLength;
        (*this)[3] /= fLength;
      }
      return *this;
    }

};




///////////////////////////////////////////////////////////////////////////////
// definition
///////////////////////////////////////////////////////////////////////////////

template<class T>
svt_vector4<T> &svt_vector4<T>::operator=(const svt_vector4<T> &that)
{
  x(that.x());
  y(that.y());
  z(that.z());
  w(that.w());
  return *this;
}

template<class T>
inline svt_vector4<T> operator-(const svt_vector4<T> &p)
{
  svt_vector4<T> v;
  v[0] = -p[0];
  v[1] = -p[1];
  v[2] = -p[2];
  v[3] =  p[3];
  return v;
}

template<class T>
inline svt_vector4<T> &svt_vector4<T>::operator*=(const svt_matrix4<T> &oMat)
{
  *this = oMat * (*this);
  return *this;
}


//
// arithmetic operators
// (need to be redefined because only first 3 components are taken into consideration)
//

// add 2 points
template<class T>
inline svt_vector4<T> operator+(const svt_vector4<T> &p1, const svt_vector4<T> &p2)
{
  return svt_vector4<T>(p1.x() + p2.x(), p1.y() + p2.y(), p1.z() + p2.z());
}


// add point and scalar
template<class T>
inline svt_vector4<T> operator+(const svt_vector4<T> &p, const T &f)
{
  return svt_vector4<T>(p.x() + f, p.y() + f, p.z() + f);
}


// add scalar and point
template<class T>
inline svt_vector4<T> operator+(const T &f, const svt_vector4<T> &p)
{
  return svt_vector4<T>(p.x() + f, p.y() + f, p.z() + f);
}

// substract 2 points
template<class T>
inline svt_vector4<T> operator-(const svt_vector4<T> &p1, const svt_vector4<T> &p2)
{
  return svt_vector4<T>(p1.x() - p2.x(), p1.y() - p2.y(), p1.z() - p2.z());
}


// substract point and scalar
template<class T>
inline svt_vector4<T> operator-(const svt_vector4<T> &p, const T &f)
{
  return svt_vector4<T>(p.x() - f, p.y() - f, p.z() - f);
}


// substract scalar and point
template<class T>
inline svt_vector4<T> operator-(const T &f, const svt_vector4<T> &p)
{
  return svt_vector4<T>(f - p.x(), f - p.y(), f - p.z());
}


// multiply/devide point with scalar
template<class T>
inline svt_vector4<T> operator*(const svt_vector4<T> &p, const T &f)
{
  return svt_vector4<T>(f * p.x(), f * p.y(), f * p.z());
}


template<class T>
inline svt_vector4<T> operator*(const T &f, const svt_vector4<T> &p)
{
  return svt_vector4<T>(f * p.x(), f * p.y(), f * p.z());
}


template<class T>
inline svt_vector4<T> operator/(const svt_vector4<T> &p, const T &f)
{
  return svt_vector4<T>(p.x() / f, p.y() / f, p.z() / f);
}


// Scalar Product
template<class T>
inline T operator*(const svt_vector4<T> &p1, const svt_vector4<T> &p2)
{
  return p1.x() * p2.x() + p1.y() * p2.y() + p1.z() * p2.z();
}


// Vector Product
template<class T>
inline svt_vector4<T> vectorProduct(const svt_vector4<T> &p1, const svt_vector4<T> &p2)
{
  return svt_vector4<T>(p1.y() * p2.z() - p1.z() * p2.y(),
                        p1.z() * p2.x() - p1.x() * p2.z(),
                        p1.x() * p2.y() - p1.y() * p2.x());
}

// Vector Product of Quaternions
template<class T>
inline svt_vector4<T> quaternionVectorProduct(const svt_vector4<T> &p1, const svt_vector4<T> &p2)
{
  return svt_vector4<T>(p1.w() * p2.x() + p1.x() * p2.w() + p1.y() * p2.z() - p1.z() * p2.y(),
                        p1.w() * p2.y() + p1.y() * p2.w() + p1.z() * p2.x() - p1.x() * p2.z(),
                        p1.w() * p2.z() + p1.z() * p2.w() + p1.x() * p2.y() - p1.y() * p2.x(),
                        p1.w() * p2.w() - p1.x() * p2.x() - p1.y() * p2.y() - p1.z() * p2.z());
}

// Product of matrix and vector
template<class T>
inline svt_vector4<T> operator*(const svt_matrix4<T> &M, const svt_vector4<T> &V)
{
  svt_vector4<T> oVec(V);
  return mat_mult(M, V, oVec);
}


// Product of matrix and vector
// ...does the same as the method above, although not aritmnetically correct...
// ...but usefull :)
template<class T>
inline svt_vector4<T> operator*(const svt_vector4<T> &V, const svt_matrix4<T> &M)
{
  svt_vector4<T> oVec(V);
  return mat_mult(M, V, oVec);
}

///////////////////////////////////////////////////////////////////////////////
// SVT_BOND
///////////////////////////////////////////////////////////////////////////////
/**
 * A bond class (bond between two atoms).
 *@author Stefan Birmanns
 */
class svt_point_cloud_bond
{

  protected:

    // bond between atom a and b
    svt_point_cloud_atom *m_pAtomA;
    svt_point_cloud_atom *m_pAtomB;

    int m_iIndexA;
    int m_iIndexB;

  public:

    /**
     * Constructor
     * create an bond between atom a and b
     * \param pA pointer to first svt_atom object
     * \param pB pointer to second svt_atom object
     */
    svt_point_cloud_bond(svt_point_cloud_atom *pA, svt_point_cloud_atom *pB, int iIndexA = -1, int iIndexB = -1);
    ~svt_point_cloud_bond();

    /**
     * set the first atom
     * \param pA pointer to svt_atom object
     */
    void setAtomA(svt_point_cloud_atom *pA);
    /**
     * get the first atom
     * \return pointer to svt_atom object
     */
    svt_point_cloud_atom *getAtomA();
    /**
     * get the index of atom a
     */
    int getIndexA();
    /**
     * get the index of atom b
     */
    int getIndexB();
    /**
     * set the index of atom a
     */
    void setIndexA(int iIndexA);
    /**
     * set the index of atom b
     */
    void setIndexB(int iIndexB);
    /**
     * set the second atom
     * \param pB pointer to svt_atom object
     */
    void setAtomB(svt_point_cloud_atom *pB);
    /**
     * get the second atom
     * \return pointer to svt_atom object
     */
    svt_point_cloud_atom *getAtomB();
};

///////////////////////////////////////////////////////////////////////////////
// SVT_SAMPLED
///////////////////////////////////////////////////////////////////////////////


/**
 * Pure abstract base class of an object that can be sampled
 * \author Stefan Birmanns
 */
template<class T> class svt_sampled
{
  public:

    /**
     */
    virtual ~svt_sampled()
    {
    };

    /**
     * sample the object randomly and return a vector that reflects the probability distribution of the object
     */
    virtual T sample() = 0;
};



///////////////////////////////////////////////////////////////////////////////
// SVT_VOLUME
///////////////////////////////////////////////////////////////////////////////

/**
 * A container class for volumetric data.
 *@author Stefan Birmanns
 */
template<class T> class svt_volume : public svt_sampled< svt_vector4<Real64> >
{

  public:

    enum GradientMode {
      CentralDistance,
      Sobel
    };

  protected:

    T *m_pData;

    unsigned int m_iSizeX;
    unsigned int m_iSizeY;
    unsigned int m_iSizeZ;

    T m_fMaxDensity;
    T m_fMinDensity;
    T m_fAvgDensity;

    Real64 m_fWidth; // voxel width in angstroem

    Real64 m_fGridX; // x position of first voxel
    Real64 m_fGridY; // y position of first voxel
    Real64 m_fGridZ; // z position of first voxel

    bool m_bInternal; // if true, the data is allocated internally and will also get deleted when the object is destroyed!

    bool m_bChanged; // true signals that the content of the volume has changed and there a calcMinMaxDensity has to be called.

    Real64 m_fCutoff; // cutoff voxel value for the sampling

    Real64 m_fNorm;  // Norm of the volume = sqrt( sum( at(i) * at(i) ) );

  public:

    /**
     * Constructor
     */
    svt_volume();

    /**
     * Constructor.
     * This constructor will create a volume of a certain size.
     * Attention: The memory will get delete when the object gets destroyed!
     * \param iSizeX size in x direction of the new volume
     * \param iSizeY size in y direction of the new volume
     * \param iSizeZ size in z direction of the new volume
     * \param fInit the value the voxels get initialized to (default: 0)
     */
    svt_volume(unsigned int iSizeX, unsigned int iSizeY, unsigned int iSizeZ, const T fInit = 0.0);

    /**
     * Constructor.
     * This constructor will create a volume of a certain size.
     * Attention: The memory will get delete when the object gets destroyed!
     * \param iSizeX size in x direction of the new volume
     * \param iSizeY size in y direction of the new volume
     * \param iSizeZ size in z direction of the new volume
     * \param pData pointer to memory block that gets copied into the new object.
     */
    svt_volume(unsigned int iSizeX, unsigned int iSizeY, unsigned int iSizeZ, const T *pData);

    /**
     * Copy constructor
     */
    svt_volume(svt_volume<T> &rVol);

    /**
     * Destructor
     */
    virtual ~svt_volume();

    ///////////////////////////////////////////////////////////////////////////
    // Operators
    ///////////////////////////////////////////////////////////////////////////

    /**
     * Assign a scalar value, e.g. A=1.
     * The size of A will be unchanged, value is assigned to all elements.
     */
    inline svt_volume<T> &operator=(const T &fValue);

    /**
     * Assign a volume B to a volume A, e.g. A=B.
     * Size and values will become same as B.
     * Attention: This makes a deep copy if the memory was allocated by the object B itself. If
     * the memory was just linked to B via a pointer (see setData()), then only the pointer
     * is copied but not the data itself. Please use deepCopy() if you wish to copy the data in such case.
     */
    svt_volume<T> &operator=(const svt_volume<T> &rThat);

    /**
     * Multiply with scalar.
     * \param fScalar the scalar the voxels of the volume are getting multiplied with.
     */
    svt_volume<T> &operator*(const T fScalar);

    /**
     * Add another volume.
     * \param rVol reference to the other volume
     */
    void operator+=(const svt_volume<T> &rVol);

    /**
     * Multiply with matrix - resulting volume has different dimensions that the original
     * \param rMat a transformation matrix that will translate/rotate the volume -
     */
    svt_volume<T> operator*(svt_matrix4< T > oMat);

    /**
     *  apply the transformation of the current volume (this) and project onto
     *  \param rDest the volume where to project
     *  \param the tranformation to apply
     *  ATTENTION: rDest should be already allocated
     */
    void applyTransformation(svt_volume<Real64> &rDest, svt_matrix4<Real64> oMat);

    /**
     * get dimension of the volume after appling transformation
     * \param the transformation to apply
     * \param rNewGrid the new origin
     * \return the new size
     */
    svt_vector4<unsigned int> getDimensions(svt_matrix4< T > oMat, svt_vector4<Real64> &rNewGrid);

    /**
     * update the dimensions, size and origin
     * \param rNewSize the size up to now, it will be updated if needed
     * \param rNewGrid the origin up to now, it will be updated
     * \param rDir gives the corner that will be consider
     */
    void updateDimensions(svt_matrix4< T > oMat, svt_vector4<Real64> &rNewGrid,  svt_vector4<unsigned int> &rNewSize, svt_vector4<unsigned int> oDir);


    ///////////////////////////////////////////////////////////////////////////
    // Public Methods
    ///////////////////////////////////////////////////////////////////////////

    /**
     * Allocate memory for a certain sized volume.
     * Attention: If there is already a dataset stored internally, then the memory will get deleted!!
     * \param iSizeX x size of the new volume
     * \param iSizeY y size of the new volume
     * \param iSizeZ z size of the new volume
     * \param fInit the value the voxels get initialized to (default: 0.0)
     * \param bInit if true the voxels get initialized to fInit, if false, the data stays uninitialized (default: true)
     */
    void allocate(unsigned int iSizeX, unsigned int iSizeY, unsigned int iSizeZ, const T fInit = 0.0, bool bInit = true);

    /**
     * Do a deep copy from another svt_volume object - regardless if internal memory or link to external memory.
     * \param rVolume reference to other object
     */
    void deepCopy(svt_volume<T> &rVolume);
    /**
     * Do a deep copy from a memory block
     * \param iSizeX x size of the new volume data
     * \param iSizeY y size of the new volume data
     * \param iSizeZ z size of the new volume data
     * \param pData pointer to data
     */
    void deepCopy(unsigned int iSizeX, unsigned int iSizeY, unsigned int iSizeZ, T *pData);

    /**
     * Set volume memory to internally allocated memory or not
     * \param bInternal if true the destructor of the object will delete memory!
     */
    void setInternal(bool bInternal);

    /**
     * Print the content of the volume
     */
    void print();

    /**
     * Links the complete volume data set to an external memory block.
     * \param iSizeX x size of the new volume data
     * \param iSizeY y size of the new volume data
     * \param iSizeZ z size of the new volume data
     * \param pData pointer to the external volume data
     */
    void setData(unsigned int iSizeX, unsigned int iSizeY, unsigned int iSizeZ, T *pData);

    /**
     * Returns the pointer to the memory location of the data.
     * Attention: If this is internally allocated memory, do not delete! The memory will get deleted automatically when the object is destroyed.
     * \return pointer to T
     */
    T *getData() const;

    /**
     * Changes one voxel value.
     * Does boundary checking!
     * \param iX x coordinate
     * \param iY y coordinate
     * \param iZ z coordinate
     * \param fValue new value
     */
    inline void setValue(unsigned int iX, unsigned int iY, unsigned int iZ, T fValue);
    /**
     * Changes all voxel values.
     * \param fValue new value
     */
    inline void setValue(T fValue);

    /**
     * Get the value at a position inside the volume.
     * Does boundary checking!
     * \param iX x coordinate
     * \param iY y coordinate
     * \param iZ z coordinate
     * \return value
     */
    virtual T getValue(unsigned int iX, unsigned int iY, unsigned int iZ) const;

    /**
     * Get the value at a position inside the volume.
     * Does boundary checking!
     * \param iCount counter
     * \return value
     */
    virtual T getValue(unsigned int iCount) const;

    /**
     * Get the value at a position inside the volume - real space positions, considers origin, trilinear interpolation.
     * \param fX x coordinate
     * \param fY y coordinate
     * \param fZ z coordinate
     * \return value
     */
    T getRealSpaceValue(Real64 fX, Real64 fY, Real64 fZ) const;

    /**
     * Get the value at a position inside the volume.
     * No boundary checking!
     * \param iX x coordinate
     * \param iY y coordinate
     * \param iZ z coordinate
     * \return value
     */
    T at(unsigned int iX, unsigned int iY, unsigned int iZ) const;
    /**
     * Get the value at a position inside the volume.
     * No boundary checking!
     * \param iCount counter
     * \return value
     */
    T at(unsigned int iCount) const;
    /**
     * Set the value at a position inside the volume.
     * No boundary checking!
     * \param iX x coordinate
     * \param iY y coordinate
     * \param iZ z coordinate
     * \param fValue new voxel value
     */
    void setAt(unsigned int iX, unsigned int iY, unsigned int iZ, T fValue);
    /**
     * Set the value at a position inside the volume.
     * No boundary checking!
     * \param iCount counter
     * \param fValue new voxel value
     */
    void setAt(unsigned int iCount, T fValue);

    /**
     * set the size
     */
    void setSize(unsigned int iSizeX, unsigned int iSizeY, unsigned int iSizeZ);
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
     * set position of first voxel
     */
    inline void setGrid(Real64 fGridX, Real64 fGridY, Real64 fGridZ);
    /**
     * set X position of first voxel
     */
    inline void setGridX(Real64 fGridX);
    /**
     * set Y position of first voxel
     */
    inline void setGridY(Real64 fGridY);
    /**
     * set Z position of first voxel
     */
    inline void setGridZ(Real64 fGridZ);
    /**
     * get position of first voxel
     */
    inline Real64 getGridX() const;
    /**
     * get position of first voxel
     */
    inline Real64 getGridY() const;
    /**
     * get position of first voxel
     */
    inline Real64 getGridZ() const;

    /**
     * get the interpolated value at a position inside the volume
     * \param fX x coordinate
     * \param fY y coordinate
     * \param fZ z coordinate
     * \return value
     */
    virtual T getIntValue(Real64 fX, Real64 fY, Real64 fZ) const;

    /**
     * get the minimum density.
     * \return the minimum density
     */
    T getMinDensity();
    /**
     * get the maximum density.
     *  \return the maximum density
     */
    T getMaxDensity();
    /**
     * get the average density.
     *  \return the average density
     */
    T getAvgDensity();

    /**
     * get norm of the volume = sqrt( sum( at(i) * at(i) ) );
     * \return norm
     */
    Real64 getNorm();

    /**
     * set voxel width (in Angstroem)
     * \param fWidth voxel width
     */
    inline void setWidth(Real64 fWidth);
    /**
     * get voxel width (in Angstroem)
     * \return voxel width
     */
    inline Real64 getWidth() const;

    /**
     * normalize volume
     */
    void normalize();

    /**
     * normalize volume
     * \param fSigma -  blurring kernel following local normalization (default fSigma=1.0 - don't blur)
     */
    void locallyNormalize(T fSigma, bool bProgress = false);

    /**
     * normalize volume
     * \param fSigma -  blurring kernel following local normalization (default fSigma=1.0 - don't blur)
     */
    void locallyNormalizeCorrectBorders(T fSigma, bool bProgress = false);

    /**
     * Threshold volume. All voxel values above a threshold value are cutoff and set to the threshold.
     * \param fMaxDensity new maximum density
     * \param fMinDensity new minimum density
     */
    void threshold(T fMinDensity, T fMaxDensity);

    /**
     * Crop volume - cuts out and returns a subvolume
     * \param iNewMinX new minimum x voxel index
     * \param iNewMaxX new maximal x voxel index
     * \param iNewMinY new minimum y voxel index
     * \param iNewMaxY new maximal y voxel index
     * \param iNewMinZ new minimum z voxel index
     * \param iNewMaxZ new maximal z voxel index
     */
    svt_volume<T> &crop(
      unsigned int iNewMinX, unsigned int iNewMaxX,
      unsigned int iNewMinY, unsigned int iNewMaxY,
      unsigned int iNewMinZ, unsigned int iNewMaxZ
    );

    /**
     * Delete a spherical subregion
     * \param fCenterX center voxel coordinate for the spherical region
     * \param fCenterY center voxel coordinate for the spherical region
     * \param fCenterZ center voxel coordinate for the spherical region
     * \param fRadius radius of the sphere
     * \return returns a new svt_volume object with the spherical region cut out
     */
    svt_volume<T> &cutSphere(Real64 fCenterX, Real64 fCenterY, Real64 fCenterZ, Real64 fRadius);

    /**
     * cuts out a spherical subregion
     * \param fCenterX center voxel coordinate for the spherical region
     * \param fCenterY center voxel coordinate for the spherical region
     * \param fCenterZ center voxel coordinate for the spherical region
     * \param fRadius radius of the sphere
     * \return returns a new svt_volume object with the cutted out spherical region
     */
    svt_volume<T> &copySphere(Real64 fCenterX, Real64 fCenterY, Real64 fCenterZ, Real64 fRadius);


    /**
     * Mask with another svt_volume object - all the voxels in this vol are multiplied with the mask volume voxels (multiplied by 0 for not overlapping voxels).
     * \param rMask reference to other object
     */
    svt_volume<T> &mask(svt_volume<T> &rMask);
    /**
     * Make a mask volume - if voxel > threshold it is set to 1, otherwise 0
     * \param fThreshold threshold value
     */
    svt_volume<T> &makeMask(T fThreshold);
    /**
     * Invert mask - all voxels have to be between 0 and 1. The voxel values are inverted.
     */
    svt_volume<T> &invertMask();

    /**
    * shrinks the original volume to the occupied Volume - redefine m_fGrid and m_iSize such that volume is not padded by zeros
    */
    void shrinkToOccupiedVolume();

    /**
     * get the new border
     * \param iMinX the min value of X
     * \param iMaxX the max value of X
     * \param iMinY the min value of Y
     * \param iMaxY the max value of Y
     * \param iMinZ the min value of Z
     * \param iMaxZ the max value of Z
     * \return the new index (voxel value)
     */
    unsigned int getNewBorder(unsigned int iMinX, unsigned int iMaxX, unsigned int iMinY, unsigned int iMaxY, unsigned int iMinZ, unsigned int iMaxZ);

    /**
     * interpolate the map to another lattice
     */
    void interpolate_map(Real64 fWidthX, Real64 fWidthY = 0.0, Real64 fWidthZ = 0.0);

    /**
     * Calc gradient map
     * \param eMode gradient method selector: CentralDistance or Sobel
     * \param rGradientX reference to an svt_volume object for the gradient map in x direction
     * \param rGradientY reference to an svt_volume object for the gradient map in y direction
     * \param rGradientZ reference to an svt_volume object for the gradient map in z direction
     */
    void calcGradient(const GradientMode eMode, svt_volume<T> &rGradientX, svt_volume<T> &rGradientY, svt_volume<T> &rGradientZ) const;

    /**
     * Set densities below limit to zero
     * \param fLimit the threshold
     */
    void threshold(T fLimit);

    /**
     * Convolve this volume with another one (typically a 3D kernel filter)
     * \param rKernel reference to kernel volume
     * \param bNormalize normalizes during convolution
     * \param bProgress show a progress bar
     */
    void convolve(svt_volume<T> &rKernel, bool bProgress);
    /**
     * Convolve this volume with another one, which is a 1D volume (only x axis). The 1D kernel will get convolved in all three directions, so this function should be
     * used with linear separable kernels.
     * \param rKernel reference to kernel volume
     * \param bNormalize normalizes during convolution
     * \param bProgress show a progress bar
     */
    void convolve1D3D(svt_volume<T> &rKernel, bool bProgress);

    /**
     * Create a Laplacian kernel volume.
     * Attention: This will overwrite the current content of the volume object with the filter kernel. It will allocate the memory internally.
     */
    void createLaplacian();

    /**
     * Create a Gaussian kernel volume within SigmaFactor*fSigma
     * Attention: This will overwrite the current content of the volume object with the filter kernel. It will allocate the memory internally.
     * \param fSigma1D sigma of map
     * \param fSigmaFactor sigma factor
     */
    void createGaussian(double fSigma1D, double fSigmaFactor);

    /**
     * Create a one-dimensional Gaussian blurring kernel volume (Situs scheme)
     * \param fWidth the voxel width of the target map one wants to convolve with the kernel
     * \param fResolution the target resolution
     * \param fVarp variance of map (if 0 no correction for lattice interpolation smoothing effects = default)
     * Attention: This will overwrite the current content of the volume object with the filter kernel. It will allocate the memory internally.
     */
    void create1DBlurringKernel(Real64 fWidth, Real64 fResolution, Real64 fVarp = 0.0);

    /**
     * Create a Gaussian blurring kernel volume (Situs scheme)
     * \param fWidth the voxel width of the target map one wants to convolve with the kernel
     * \param fResolution the target resolution
     * \param fVarp variance of map (if 0 no correction for lattice interpolation smoothing effects = default)
     */
    void createSitusBlurringKernel(Real64 fWidth, Real64 fResolution, Real64 fVarp = 0.0);

    /**
     * Create a Laplacian of a Gaussian kernel volume.
     * \param fWidth the voxel width of the target map one wants to convolve with the kernel
     * \param fResolution the target resolution
     * \param fVarp variance of map (if 0 no correction for lattice interpolation smoothing effects = default)
     * Attention: This will overwrite the current content of the volume object with the filter kernel. It will allocate the memory internally.
     */
    void createLaplacianOfGaussian(Real64 fWidth, Real64 fResolution, Real64 fVarp = 0.0);

    /**
      * calculate Fischer Determinant
      * \param fVectorVal_VarianceOnePath reference vector with variance
      * \param fVectorVal_VarianceOnePath_Order reference for vector with descending order of fVectorVal_VarianceOnePath
      * \param iK_star reference for iK_star that maximalizes Functional F(K) and divided fVectorVal_VarianceOnePath for two classes.
      */
    void FisherDiscriminant(vector <Real64> &fVectorVal_VarianceOnePath, vector <unsigned int > &iVectorVal_VarianceOnePath_Order, unsigned int &iK_star);


    /**
     * Filtrate Volume with Bilateral Filter
     * Attention: This will overwrite the current content of the volume object with the filter kernel. It will allocate the memory internally.
     * \param fSigma1D1 sigma of Gaussuian part of filter
     * \param fsigma1D2 sigma of intensity part of filter
     * \param iSize size of the kernel
     * \param bProgress show a progress bar
     */
    void filtrBilateral(double fSigma1D1, double fSigma1D2, unsigned int iSize,  bool bProgress);

    /**
     * calculate One Geodesic Step in path based on Volume
     * \param oPathVector store the output Geodesic path
     * \param oPathSoFar  store current Geodesic path
     * \param oToPos next (x,y,z) position of theoretical jumping particle
     * \param iFullLength the required length of path
     * \param iCubeSize the size of Kernel
     */

    void stepGeodesic(vector < vector<svt_vector4< int> > > &oPathVector, vector<svt_vector4< int> > &oPathSoFar, svt_vector4< int> oToPos, int iStepNumber, int iFullLength, int iCubeSize, int iNeigboorModel);

    /**
      * calculate all Geodesic Path for Volume
     * \param oPathVector store the output Geodesic path
     * \param oPathSoFar  store current Geodesic path
     * \param oToPos next (x,y,z) position of theoretical jumping particle
     * \param iFullLength the required length of path
     * \param iCubeSize the size of Kernel
     */

    void findGeodesicPaths(vector < vector<svt_vector4< int> > > &oPathVector, svt_vector4< int> oToPos, int iFullLength, int iCubeSize, int iNeigboorModel);


    /**
      * Filtrate Volume with Geodesic Path Filter
      * Attention: This will overwrite the current content of the volume object with the filter kernel. It will allocate the memory internally.
      * \param fBetha the filter equivalent of Gaussian Sigma
      * \param iMaskDim kernel size
      * \param iPathLength requred path length
      * \param bProgress show a progress bar
     */
    void filtrGeodesicPath(double fBetha, int iNeigboorModel, unsigned int iMaskDim,  unsigned int iPathLength,  bool bProgress);

    /**
     * Create an Identity kernel volume
     * Attention: This will overwrite the current content of the volume object with the filter kernel. It will allocate the memory internally.
     */
    void createIdentity();

    /**
     * Create sphere at the center of the volume
     * \param fRadius radius of the sphere
     */
    void createSphere(Real64 fRadius);

    /**
     * Applies a Laplacian filter.
     * \param bRelax if true the poisson relaxation procedure is applied
     */
    void applyLaplacian(bool bRelax = true)
    {
      svt_volume<T> oLaplace;
      oLaplace.createLaplacian();

      // relaxation procedure
      if (bRelax) {
        //unsigned int iIgnored[3];  /* zero margin ignored in fast kernel convolution */
        //iIgnored[0] = (m_iSizeX+1)/4; iIgnored[1] = (m_iSizeY+1)/4; iIgnored[2] = (m_iSizeZ+1)/4;
        //pphiLaplace = new MYFLT[iNvox];
        //int x;
        //unsigned int iNVox = m_iSizeX*m_iSizeY*m_iSizeZ;
        //for (x=0;x<=iNVox;x++)
        //    pphiLaplace[x] = pphi[x];
        //relax_laplacian(&pphiLaplace, m_iExtX, m_iExtY, m_iExtZ, iIgnored, 5.0);
      }

      convolve(oLaplace, false);
    };

    /**
     * Calculate correlation with other svt_volume object
     * \param rVol reference to other svt_volume object
     */
    Real64 correlation(svt_volume<T> &rVol, bool bMask = true);

    /**
     * Calculate correlation with other svt_volume object
     * \param rVol reference to other svt_volume object
     */
    Real64 correlationColores(svt_volume<T> &rVol);

    /**
     * Calculates the internal correlation of the voxel values.
     */
    svt_volume<Real64> *internalCorr(unsigned int iWidth);

    /**
     * Grow/Shrink - this function blows up/interpolates the volume at a different size (not a padding routine!).
     * \param iNewSizeX new size (x dimension)
     * \param iNewSizeY new size (y dimension)
     * \param iNewSizeZ new size (z dimension)
     */
    void resize(unsigned int iNewSeizeX, unsigned int iNewSizeY, unsigned int iNewSizeZ);

    /**
     * Calculate number of occupied voxels
     * \param fThreshold threshold value - only voxels higher than the threshold are counted.
     * \return number of voxels with a density higher than fThreshold
     */
    unsigned long getOccupied(T fThreshold) const;

    /**
     * Load a file. This function looks at the extension to determine which function has to be used to actually load the file.
     * \param pFname pointer to array of char with the filename
     */
    svt_matrix4<T> load(const char *pFname);
    /**
     * Save a file. This function looks at the extension to determine which function has to be used to actually save the file.
     * \param pFname pointer to array of char with the filename
     */
    void save(const char *pFname);

    /**
     * Set cutoff for sampling
     * \param fCutoff a voxel value lower than this value is not considered for the sampling
     */
    void setCutoff(Real64 fCutoff);
    /**
     * Get cutoff for sampling
     * \return a voxel value lower than this value is not considered for the sampling
     */
    Real64 getCutoff() const;

    /**
     * sample the object randomly and return a vector that refrects the probability distribution of the object
     */
    svt_vector4<Real64> sample();

    /**
     * sample the object in the sphere randomly and return a vector that refrects the probability distribution of the object
     * \param oCenter keeps coordinates of center of sphere
     * \param fRadius is radius of sphere
     */
    svt_vector4<Real64> sampleSphere(svt_vector4<Real64> oCenter, Real64 fRadius);

    /**
     * Non-recursive flood-fill segmentation algorithm. All voxels that are connected to a start voxel and are above the threshold are kept, the others are removed.
     * The algorithm creates a mask that is later blurred by a Gaussian kernel. The sigma of the gaussian can be specified.
     * \param iStartX x index of the start voxel
     * \param iStartY y index of the start voxel
     * \param iStartZ z index of the start voxel
     * \param fTreshold threshold for the floodfill
     * \param fGaussian sigma of the gaussian the mask is convoluted with (if 0, no blurring of the mask is applied)
     */
    void floodfill_segmentation(unsigned int iStartX, unsigned int iStartY, unsigned int iStartZ, T fThreshold, Real64 fGaussian);

    /**
     * Non-recursive flood-fill algorithm. The algorithm fills the voxels that are above the threshold (and that are connected to the start voxel) with the specified value.
     * \param iStartX x index of the start voxel
     * \param iStartY y index of the start voxel
     * \param iStartZ z index of the start voxel
     * \param fTreshold threshold for the floodfill
     * \param fGaussian sigma of the gaussian the mask is convoluted with (if 0, no blurring of the mask is applied)
     */
    void floodfill(unsigned int iStartX, unsigned int iStartY, unsigned int iStartZ, T fThreshold, T fFillValue);

    /**
    * gives the best ISO threshold for fitting the object with a given volume
    */
    T bestISO(svt_volume<T> &rTarget, T fThreshold);

    /**
     * Calculates the corresponding isosurface threshold value for this volume. The algorithm tries to cover exactly the voxel values of the other volume with the ones from here, so that
     * both isosurfaces have a similar size. rVol_A would typically have a higher resolution but show the exact same system.
     *
     * \param rVol_A other volume
     * \param fThresh good known threshold value of rVol_A
     */
    T correspondingISO(svt_volume<T> &rVol_A, T fThreshold);

    /**
    * scales the entire map by a given factor
    */
    void scale(T fScale);

    /**
    * shifts the entire map by a given factor
    */
    void shift(T fShift);

    /**
    * scales by the slope
    */
    void scaleBySlope(T fSlope);

    /**
    * returns the slope given two points
    */
    T getSlope(T fX1, T fY1, T fX2, T fY2);

    /**
    * does an interpolated scaling of the object to match another
    */
    void interpolatedScale(T fThreshold1, T fThreshold2, T fNew1, T fNew2);

    /**
    * calculates the root mean square of the volume
    */
    T getRMS();

    /**
     * Remove every structure that has more than two neighbors
     * \param fThreshold voxel above that value are considered occupied
     */
    void removeNeighbors(T fThreshold);

    /**
     * How many direct neighbors does a voxel have?
     * \return number of occupied neighboring voxels
     */
    unsigned int numNeighbors(unsigned int iX, unsigned int iY, unsigned int iZ, T fThreshold);

  private:

    /**
     * Calculate/update the minimum and maximum density values.
     */
    void calcMinMaxDensity();

    /**
     * Simple internal routine that calculates the index of a voxel, depending on the data order
     * \param iIndex number of the voxel
     * \param iDataOrder data order
     */
    unsigned int calc_xyz_index(unsigned int iIndex, unsigned int iDataOrder);

    /**
     * Simple internal routine that calculates the x,y,z index of a voxel
     * \param iIndex number of the voxel
     * \param iDataOrder data order
     * \param rX reference to unsigned int
     * \param rY reference to unsigned int
     * \param rZ reference to unsigned int
     */
    void calc_xyz(unsigned int iIndex, unsigned int iDataOrder, unsigned int iMaxX, unsigned int iMaxY, unsigned int iMaxZ, unsigned int &rX, unsigned int &rY, unsigned int &rZ);
};

///////////////////////////////////////////////////////////////////////////////
// Implementation
///////////////////////////////////////////////////////////////////////////////

/**
 * Constructor
 */
template<class T>
inline svt_volume<T>::svt_volume() :
  m_pData(NULL),
  m_fWidth(1.0),
  m_fGridX(0.0),
  m_fGridY(0.0),
  m_fGridZ(0.0),
  m_bInternal(false),
  m_bChanged(true),
  m_fCutoff(-1.0E30),
  m_fNorm(0.0)
{
  setData(0, 0, 0, (T *)NULL);
};

/**
 * Constructor.
 * This constructor will create a volume of a certain size. Attention: The memory will get delete when the object gets destroyed!
 * \param iSizeX size in x direction of the new volume
 * \param iSizeY size in y direction of the new volume
 * \param iSizeZ size in z direction of the new volume
 * \param fInit the value the voxels get initialized to
 */
template<class T>
inline svt_volume<T>::svt_volume(unsigned int iSizeX, unsigned int iSizeY, unsigned int iSizeZ, const T fInit) :
  m_pData(NULL),
  m_fWidth(1.0),
  m_fGridX(0.0),
  m_fGridY(0.0),
  m_fGridZ(0.0),
  m_bInternal(true),
  m_bChanged(true),
  m_fCutoff(-1.0E30),
  m_fNorm(0.0)
{
  allocate(iSizeX, iSizeY, iSizeZ, fInit);
};

/**
 * Constructor.
 * This constructor will create a volume of a certain size.
 * Attention: The memory will get delete when the object gets destroyed!
 * \param iSizeX size in x direction of the new volume
 * \param iSizeY size in y direction of the new volume
 * \param iSizeZ size in z direction of the new volume
 * \param pData pointer to memory block that gets copied into the new object.
 */
template<class T>
inline svt_volume<T>::svt_volume(unsigned int iSizeX, unsigned int iSizeY, unsigned int iSizeZ, const T *pData) :
  m_pData(NULL),
  m_fWidth(1.0),
  m_fGridX(0.0),
  m_fGridY(0.0),
  m_fGridZ(0.0),
  m_bInternal(true),
  m_bChanged(true),
  m_fCutoff(-1.0E30),
  m_fNorm(0.0)
{
  allocate(iSizeX, iSizeY, iSizeZ);

  memcpy(m_pData, pData, sizeof(T)*iSizeX * iSizeY * iSizeZ);
};

/**
 * Copy constructor
 */
template<class T>
inline svt_volume<T>::svt_volume(svt_volume<T> &rVol) : svt_sampled< svt_vector4<Real64> >(),
  m_pData(NULL),
  m_fWidth(1.0),
  m_fGridX(0.0),
  m_fGridY(0.0),
  m_fGridZ(0.0),
  m_bInternal(true),
  m_bChanged(true),
  m_fCutoff(-1.0E30),
  m_fNorm(0.0)
{
  *this = rVol;
};

/**
 * Destructor
 */
template<class T>
inline svt_volume<T>::~svt_volume()
{
  if (m_bInternal)
    delete[] m_pData;
};

///////////////////////////////////////////////////////////////////////////////
// Operators
///////////////////////////////////////////////////////////////////////////////

/**
 * Assign a scalar value, e.g. A=1.
 * The size of A will be unchanged, value is assigned to all elements.
 */
template<class T>
inline svt_volume<T> &svt_volume<T>::operator=(const T &fValue)
{
  setValue(fValue);

  return *this;
};

/**
 * Assign a volume B to a volume A, e.g. A=B.
 * Size and values will become same as B.
 * Attention: This makes a deep copy if the memory was allocated by the object B itself. If
 * the memory was just linked to B via a pointer (see setData()), then only the pointer
 * is copied but not the data itself. Please use deepCopy() if you wish to copy the data in such case.
 */
template<class T>
svt_volume<T> &svt_volume<T>::operator=(const svt_volume<T> &rVol)
{
  m_iSizeX = rVol.m_iSizeX;
  m_iSizeY = rVol.m_iSizeY;
  m_iSizeZ = rVol.m_iSizeZ;

  m_fMaxDensity = rVol.m_fMaxDensity;
  m_fMinDensity = rVol.m_fMinDensity;

  m_fWidth = rVol.m_fWidth;

  m_fGridX = rVol.m_fGridX;
  m_fGridY = rVol.m_fGridY;
  m_fGridZ = rVol.m_fGridZ;

  m_bInternal = rVol.m_bInternal;

  if (m_pData != NULL && m_bInternal) {
    //SVTLBBO << "Attention: memory for new volume gets allocated, " << endl;
    //SVTLBBO << "   but object had already content (which gets destroyed)" << endl;
    delete[] m_pData;
  }

  // deep copy or just pointer copy
  if (rVol.m_bInternal == true) {
    unsigned int iNum = m_iSizeX * m_iSizeY * m_iSizeZ;
    m_pData = new T[iNum];
    memcpy(m_pData, rVol.m_pData, iNum * sizeof(T));

  } else
    m_pData = rVol.m_pData;

  m_bChanged = rVol.m_bChanged;
  m_fCutoff = rVol.m_fCutoff;
  m_fNorm = rVol.m_fNorm;

  return *this;
};

/**
 * Multiply with scalar.
 * \param fScalar the scalar the voxels of the volume are getting multiplied with.
 */
template<class T>
svt_volume<T> &svt_volume<T>::operator*(const T fScalar)
{
  unsigned int iNum = m_iSizeX * m_iSizeY * m_iSizeZ;

  for (unsigned int i = 0; i < iNum; i++)
    setAt(i, at(i) * fScalar);

  m_bChanged = true;

  return *this;
};
/**
 * Add another volume.
 * \param rVol reference to the other volume
 */
template<class T>
void svt_volume<T>::operator+=(const svt_volume<T> &rVol)
{
  if (m_fGridX != rVol.getGridX() || m_fGridY != rVol.getGridY() || m_fGridZ != rVol.getGridZ() ||
      m_iSizeX != rVol.getSizeX() || m_iSizeY != rVol.getSizeY() || m_iSizeZ != rVol.getSizeZ() ||
      m_fWidth != rVol.getWidth()) {
    error_sba(85010, "Can not execute += operator as the maps have different size!");
    return;
  }

  unsigned int iNum = m_iSizeX * m_iSizeY * m_iSizeZ;

  for (unsigned int i = 0; i < iNum; i++)
    setAt(i, at(i) + rVol.at(i));

  m_bChanged = true;
};

/**
 * Multiply with matrix - resulting volume has different dimensions that the original
 * \param rMat a transformation matrix that will translate/rotate the volume -
 */
template< class T>
svt_volume<T> svt_volume<T>::operator*(svt_matrix4< T > oMat)
{
  svt_vector4<Real64> oGrid;
  svt_vector4<unsigned int> oSize = getDimensions(oMat, oGrid);

  //first get the size of the new volume
  svt_volume<T> oVol(oSize.x(), oSize.y(), oSize.z());
  oVol.setGridX(oGrid.x());
  oVol.setGridY(oGrid.y());
  oVol.setGridZ(oGrid.z());
  oVol.setWidth(m_fWidth);

  //apply transformation to this and project on oVol
  applyTransformation(oVol, oMat);

  oVol.shrinkToOccupiedVolume();

  return oVol;
};

/**
 *  apply the transformation of the current volume (this) and project onto
 *  \param rDest the volume where to project
 *  \param the tranformation to apply
 *  ATTENTION: rDest should be already allocated
 */
template< class T>
void svt_volume<T>::applyTransformation(svt_volume<Real64> &rDest, svt_matrix4<Real64> oMat)
{
  oMat.invert(); // get the inverse transformation

  svt_vector4<Real64> oVec, oInvVec;
  int iX0, iY0, iZ0, iX1, iY1, iZ1;
  Real64 fGx, fGy, fGz, fA, fB, fC;
  Real64 fNewDensity;
  for (unsigned int iX = 0; iX < rDest.getSizeX(); iX++) {
    oVec.x(rDest.getGridX() + iX * m_fWidth);
    for (unsigned int iY = 0; iY < rDest.getSizeY(); iY++) {
      oVec.y(rDest.getGridY() + iY * m_fWidth);
      for (unsigned int iZ = 0; iZ < rDest.getSizeZ(); iZ++) {
        oVec.z(rDest.getGridZ() + iZ * m_fWidth);
        oInvVec = oMat * oVec; // this vector is in the non rotated spate

        fGx = (oInvVec.x() - m_fGridX) / m_fWidth;
        fGy = (oInvVec.y() - m_fGridY) / m_fWidth;
        fGz = (oInvVec.z() - m_fGridZ) / m_fWidth;

        iX0 = floor(fGx);
        iY0 = floor(fGy);
        iZ0 = floor(fGz);

        iX1 = ceil(fGx);
        iY1 = ceil(fGy);
        iZ1 = ceil(fGz);

        fA = fGx - iX0;
        fB = fGy - iY0;
        fC = fGz - iZ0;

        /* interpolate */

        if (iX0 >= 0 && iX0 <= (int)m_iSizeX && iY0 >= 0 && iY0 <= (int)m_iSizeY && iZ0 >= 0 && iZ0 <= (int)m_iSizeZ &&
            iX1 >= 0 && iX1 <= (int)m_iSizeX && iY1 >= 0 && iY1 <= (int)m_iSizeY && iZ1 >= 0 && iZ1 <= (int)m_iSizeZ) { // inside the original space
          fNewDensity =
            fA  * fB    * fC    * this->getValue(iX1, iY1, iZ1) +
            (1 - fA)  * fB    * fC    * this->getValue(iX0, iY1, iZ1) +
            fA  * (1 - fB)  * fC    * this->getValue(iX1, iY0, iZ1) +
            fA  * fB    * (1 - fC)  * this->getValue(iX1, iY1, iZ0) +
            fA  * (1 - fB)  * (1 - fC)  * this->getValue(iX1, iY0, iZ0) +
            (1 - fA)  * fB    * (1 - fC)  * this->getValue(iX0, iY1, iZ0) +
            (1 - fA)  * (1 - fB)  * fC    * this->getValue(iX0, iY0, iZ1) +
            (1 - fA)  * (1 - fB)  * (1 - fC)  * this->getValue(iX0, iY0, iZ0);
        } else
          fNewDensity = 0.0;

        rDest.setValue(iX, iY, iZ, rDest.at(iX, iY, iZ) + fNewDensity);
      }
    }
  }

}

/**
 * get dimension of the volume after appling transformation
 * \param the transformation to apply
 * \param rNewGrid the new origin
 * \return the new size
 */
template< class T>
svt_vector4<unsigned int> svt_volume<T>::getDimensions(svt_matrix4< T > oMat, svt_vector4<Real64> &rNewGrid)
{
  svt_vector4<unsigned int> oSize;
  oSize.x(0);
  oSize.y(0);
  oSize.z(0);

  //init the first origin
  rNewGrid.x(m_fGridX);
  rNewGrid.y(m_fGridY);
  rNewGrid.z(m_fGridZ);
  rNewGrid = oMat * rNewGrid;

  //update grid and size
  updateDimensions(oMat, rNewGrid, oSize, svt_vector4<unsigned int>(1, 0, 0));
  updateDimensions(oMat, rNewGrid, oSize, svt_vector4<unsigned int>(0, 1, 0));
  updateDimensions(oMat, rNewGrid, oSize, svt_vector4<unsigned int>(0, 0, 1));
  updateDimensions(oMat, rNewGrid, oSize, svt_vector4<unsigned int>(1, 1, 0));
  updateDimensions(oMat, rNewGrid, oSize, svt_vector4<unsigned int>(0, 1, 1));
  updateDimensions(oMat, rNewGrid, oSize, svt_vector4<unsigned int>(1, 0, 1));
  updateDimensions(oMat, rNewGrid, oSize, svt_vector4<unsigned int>(1, 1, 1));

  //rNewGrid is now the origin given by the rotated bounding box; not necessary on the same grid as this

  svt_vector4<Real64> oFarCorner;
  oFarCorner.x(rNewGrid.x() + m_fWidth * oSize.x());
  oFarCorner.y(rNewGrid.y() + m_fWidth * oSize.y());
  oFarCorner.z(rNewGrid.z() + m_fWidth * oSize.z());

  //put volume on the original grid
  rNewGrid.x(m_fGridX + floor((rNewGrid.x() - m_fGridX) / m_fWidth)*m_fWidth);
  rNewGrid.y(m_fGridY + floor((rNewGrid.y() - m_fGridY) / m_fWidth)*m_fWidth);
  rNewGrid.z(m_fGridZ + floor((rNewGrid.z() - m_fGridZ) / m_fWidth)*m_fWidth);

  oSize.x(ceil((oFarCorner.x() - rNewGrid.x()) / m_fWidth));
  oSize.y(ceil((oFarCorner.y() - rNewGrid.y()) / m_fWidth));
  oSize.z(ceil((oFarCorner.z() - rNewGrid.z()) / m_fWidth));

  return oSize;
}

/**
 * update the dimensions, size and origin
 * \param rNewSize the size up to now, it will be updated if needed
 * \param rNewGrid the origin up to now, it will be updated
 * \param rDir gives the corner that will be consider
 */
template< class T>
void svt_volume<T>::updateDimensions(svt_matrix4< T > oMat, svt_vector4<Real64> &rNewGrid,  svt_vector4<unsigned int> &rNewSize, svt_vector4<unsigned int> oDir)
{
  svt_vector4<Real64> oNewCorner, oCorner;
  oCorner.x(m_fGridX + m_fWidth * oDir.x()*m_iSizeX);
  oCorner.y(m_fGridY + m_fWidth * oDir.y()*m_iSizeY);
  oCorner.z(m_fGridZ + m_fWidth * oDir.z()*m_iSizeZ);
  //get the coordinates of the new corner
  oNewCorner = oMat * oCorner;

  //origin update
  if (rNewGrid.x() > oNewCorner.x())
    rNewGrid.x(oNewCorner.x());

  if (rNewGrid.y() > oNewCorner.y())
    rNewGrid.y(oNewCorner.y());

  if (rNewGrid.z() > oNewCorner.z())
    rNewGrid.z(oNewCorner.z());

  //size update
  unsigned int iNewSize;
  iNewSize =  ceil(abs(rNewGrid.x() - oNewCorner.x()) / m_fWidth);
  if (iNewSize > rNewSize.x())
    rNewSize.x(iNewSize);

  iNewSize =  ceil(abs(rNewGrid.y() - oNewCorner.y()) / m_fWidth);
  if (iNewSize > rNewSize.y())
    rNewSize.y(iNewSize);

  iNewSize =  ceil(abs(rNewGrid.z() - oNewCorner.z()) / m_fWidth);
  if (iNewSize > rNewSize.z())
    rNewSize.z(iNewSize);
}




///////////////////////////////////////////////////////////////////////////////
// Public functions
///////////////////////////////////////////////////////////////////////////////

/**
 * Allocate memory for a certain sized volume.
 * Attention: If there is already a dataset stored internally, then the memory will get deleted!!
 * \param iSizeX x size of the new volume
 * \param iSizeY y size of the new volume
 * \param iSizeZ z size of the new volume
 * \param fInit the value the voxels get initialized to
 * \param bInit if true the voxels get initialized to fInit, if false, the data stays uninitialized.
 */
template<class T>
inline void svt_volume<T>::allocate(unsigned int iSizeX, unsigned int iSizeY, unsigned int iSizeZ, const T fInit, bool bInit)
{
  if (m_pData != NULL && m_bInternal) {
    //SVTLBBO << "Attention: memory for new volume gets allocated, " << endl;
    //SVTLBBO << "   but object had already content (which gets destroyed)" << endl;
    delete[] m_pData;
  } else if (m_pData != NULL && !m_bInternal) {
    SVTLBBO << "Attention: memory for new volume gets allocated, " << endl;
    SVTLBBO << "   but object is still linked to something else..." << endl;
  }

  unsigned int iNum = iSizeX * iSizeY * iSizeZ;
  m_pData = new T[iNum];
  m_bInternal = true;

  m_iSizeX = iSizeX;
  m_iSizeY = iSizeY;
  m_iSizeZ = iSizeZ;

  if (bInit) {
    m_fMaxDensity = fInit;
    m_fMinDensity = fInit;

    unsigned int i;
    for (i = 0; i < iNum; i++)
      m_pData[i] = fInit;
  }
};

/**
 * Do a deep copy from another svt_volume object - regardless if internal memory or link to external memory.
 * \param rVolume reference to other object
 */
template<class T>
inline void svt_volume<T>::deepCopy(svt_volume<T> &rVolume)
{
  // use assign operator to copy all the member vars
  *this = rVolume;

  // copy data anyway even if external
  if (m_bInternal == false) {
    unsigned int iNum = m_iSizeX * m_iSizeY * m_iSizeZ;
    m_pData = new T[iNum];
    memcpy(m_pData, rVolume.m_pData, iNum * sizeof(T));
    m_bInternal = true;
  }

  m_bChanged = true;
};

/**
 * Do a deep copy from a memory block
 * \param iSizeX x size of the new volume data
 * \param iSizeY y size of the new volume data
 * \param iSizeZ z size of the new volume data
 * \param pData pointer to data
 */
template<class T>
inline void svt_volume<T>::deepCopy(unsigned int iSizeX, unsigned int iSizeY, unsigned int iSizeZ, T *pData)
{
  allocate(iSizeX, iSizeY, iSizeZ);
  unsigned int iNum = iSizeX * iSizeY * iSizeZ;
  memcpy(m_pData, pData, iNum * sizeof(T));

  m_bChanged = true;
}

/**
 * Set volume memory to internally allocated memory or not
 * \param bInternal if true the destructor of the object will delete memory!
 */
template<class T>
void svt_volume<T>::setInternal(bool bInternal)
{
  m_bInternal = bInternal;
};

/**
 * Print the content of the volume
 */
template<class T>
inline void svt_volume<T>::print()
{

  unsigned int iX, iY, iZ;

  SVTLBBO << endl;
  for (iZ = 0; iZ < m_iSizeZ; iZ++) {
    for (iY = 0; iY < m_iSizeY; iY++) {
      SVTLBBO << " ";
      for (iX = 0; iX < m_iSizeX; iX++) {
        cout << at(iX, iY, iZ);
        if (iX < m_iSizeX - 1)
          cout << ", ";
      }
      cout << endl;
    }
    SVTLBBO << endl;
  }
};

/**
 * Links the complete volume data set to an external memory block.
 * \param iSizeX x size of the new volume data
 * \param iSizeY y size of the new volume data
 * \param iSizeZ z size of the new volume data
 * \param pData pointer to the external volume data
 */
template<class T>
inline void svt_volume<T>::setData(unsigned int iSizeX, unsigned int iSizeY, unsigned int iSizeZ, T *pData)
{
  m_iSizeX = iSizeX;
  m_iSizeY = iSizeY;
  m_iSizeZ = iSizeZ;
  m_pData = pData;

  m_bChanged = true;
  m_bInternal = false;
};

/**
 * Returns the pointer to the memory location of the data
 * \return pointer to T
 */
template<class T>
inline T *svt_volume<T>::getData() const
{
  return m_pData;
};

/**
 * Changes one voxel value.
 * Does boundary checking!
 * \param iX x coordinate
 * \param iY y coordinate
 * \param iZ z coordinate
 * \param fValue new value
 */
template<class T>
inline void svt_volume<T>::setValue(unsigned int iX, unsigned int iY, unsigned int iZ, T fValue)
{
  if (iX < m_iSizeX && iY < m_iSizeY && iZ < m_iSizeZ) {
    m_pData[iX + iY * m_iSizeX + iZ * m_iSizeX * m_iSizeY] = fValue;
    m_bChanged = true;
  }
};
/**
 * Changes all voxel values.
 * \param fValue new value
 */
template<class T>
inline void svt_volume<T>::setValue(T fValue)
{
  unsigned int iNum = m_iSizeX * m_iSizeY * m_iSizeZ;
  unsigned int i;
  for (i = 0; i < iNum; i++)
    m_pData[i] = fValue;

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
template<class T>
inline T svt_volume<T>::getValue(unsigned int iX, unsigned int iY, unsigned int iZ) const
{
  if (iX < m_iSizeX && iY < m_iSizeY && iZ < m_iSizeZ)
    return m_pData[iX + (iY * m_iSizeX) + (iZ * m_iSizeX * m_iSizeY)];
  else
    return (T)(0);
};

/**
 * get the value at a position inside the volume
 * \param iCount counter
 * \return value
 */
template<class T>
inline T svt_volume<T>::getValue(unsigned int iCount) const
{
  if (iCount < m_iSizeX * m_iSizeY * m_iSizeZ)
    return m_pData[iCount];
  else
    return (T)(0.0);
};

/**
 * Get the value at a position inside the volume - real space positions, considers origin, trilinearly interpolated.
 * \param fX x coordinate
 * \param fY y coordinate
 * \param fZ z coordinate
 * \return value
 */
template<class T>
inline T svt_volume<T>::getRealSpaceValue(Real64 fX, Real64 fY, Real64 fZ) const
{
  fX -= m_fGridX;
  fY -= m_fGridY;
  fZ -= m_fGridZ;

  fX /= m_fWidth;
  fY /= m_fWidth;
  fZ /= m_fWidth;

  return getIntValue(fX, fY, fZ);
};


/**
 * Get the value at a position inside the volume.
 * No boundary checking!
 * \param iX x coordinate
 * \param iY y coordinate
 * \param iZ z coordinate
 * \return value
 */
template<class T>
inline T svt_volume<T>::at(unsigned int iX, unsigned int iY, unsigned int iZ) const
{
  return m_pData[iX + (iY * m_iSizeX) + (iZ * m_iSizeX * m_iSizeY)];
};
/**
 * Get the value at a position inside the volume.
 * No boundary checking!
 * \param iCount counter
 * \return value
 */
template<class T>
inline T svt_volume<T>::at(unsigned int iCount) const
{
  return m_pData[iCount];
};

/**
 * Set the value at a position inside the volume.
 * No boundary checking!
 * \param iX x coordinate
 * \param iY y coordinate
 * \param iZ z coordinate
 * \param fValue new voxel value
 */
template<class T>
inline void svt_volume<T>::setAt(unsigned int iX, unsigned int iY, unsigned int iZ, T fValue)
{
  m_pData[iX + (iY * m_iSizeX) + (iZ * m_iSizeX * m_iSizeY)] = fValue;
  m_bChanged = true;
};
/**
 * Set the value at a position inside the volume.
 * No boundary checking!
 * \param iCount counter
 * \param fValue new voxel value
 */
template<class T>
inline void svt_volume<T>::setAt(unsigned int iCount, T fValue)
{
  m_pData[iCount] = fValue;
  m_bChanged = true;
};

/**
 * get the x size of the volume
 * \return x size of the volume
 */
template<class T>
inline unsigned int svt_volume<T>::getSizeX() const
{
  return m_iSizeX;
};
/**
 * get the y size of the volume
 * \return y size of the volume
 */
template<class T>
inline unsigned int svt_volume<T>::getSizeY() const
{
  return m_iSizeY;
};
/**
 * get the z size of the volume
 * \return z size of the volume
 */
template<class T>
inline unsigned int svt_volume<T>::getSizeZ() const
{
  return m_iSizeZ;
};
/**
 * Get the number of voxels of the volume.
 * \return number of voxels of the volume
 */
template<class T>
inline unsigned int svt_volume<T>::size() const
{
  return m_iSizeX * m_iSizeY * m_iSizeZ;
};

/**
 * set position of first voxel
 */
template<class T>
inline void svt_volume<T>::setGrid(Real64 fGridX, Real64 fGridY, Real64 fGridZ)
{
  m_fGridX = fGridX;
  m_fGridY = fGridY;
  m_fGridZ = fGridZ;
};
/**
 * set X position of first voxel
 */
template<class T>
inline void svt_volume<T>::setGridX(Real64 fGridX)
{
  m_fGridX = fGridX;
};
/**
 * set Y position of first voxel
 */
template<class T>
inline void svt_volume<T>::setGridY(Real64 fGridY)
{
  m_fGridY = fGridY;
};
/**
 * set Z position of first voxel
 */
template<class T>
inline void svt_volume<T>::setGridZ(Real64 fGridZ)
{
  m_fGridZ = fGridZ;
};
/**
 * get position of first voxel
 */
template<class T>
inline Real64 svt_volume<T>::getGridX() const
{
  return m_fGridX;
};
/**
 * get position of first voxel
 */
template<class T>
inline Real64 svt_volume<T>::getGridY() const
{
  return m_fGridY;
};
/**
 * get position of first voxel
 */
template<class T>
inline Real64 svt_volume<T>::getGridZ() const
{
  return m_fGridZ;
};

/**
 * get the interpolated value at a position inside the volume
 * \param fX x coordinate
 * \param fY y coordinate
 * \param fZ z coordinate
 * \return value
 */
template<class T>
inline T svt_volume<T>::getIntValue(Real64 fX, Real64 fY, Real64 fZ) const
{
  if (fX > (Real64)(m_iSizeX - 1) || fY > (Real64)(m_iSizeY - 1) || fZ > (Real64)(m_iSizeZ - 1) || fX < 0 || fY < 0 || fZ < 0)
    return (T)(0);

  int iV0X = (int)(floor(fX));
  int iV0Y = (int)(floor(fY));
  int iV0Z = (int)(floor(fZ));
  int iV1X = (int)(ceil(fX));
  int iV1Y = (int)(ceil(fY));
  int iV1Z = (int)(ceil(fZ));
  T fA = (T)(iV1X) - (T)(fX);
  T fB = (T)(iV1Y) - (T)(fY);
  T fC = (T)(iV1Z) - (T)(fZ);
  T fValue = (T)(0);

  // 0 0 0
  fValue += getValue(iV0X, iV0Y, iV0Z) * (fA * fB * fC);
  // 0 0 1
  fValue += getValue(iV1X, iV0Y, iV0Z) * (((T)(1) - fA) * fB * fC);
  // 0 1 0
  fValue += getValue(iV0X, iV1Y, iV0Z) * (fA * ((T)(1) - fB) * fC);
  // 0 1 1
  fValue += getValue(iV1X, iV1Y, iV0Z) * (((T)(1) - fA) * ((T)(1) - fB) * fC);
  // 1 0 0
  fValue += getValue(iV0X, iV0Y, iV1Z) * (fA * fB * ((T)(1) - fC));
  // 1 0 1
  fValue += getValue(iV1X, iV0Y, iV1Z) * (((T)(1) - fA) * fB * ((T)(1) - fC));
  // 1 1 0
  fValue += getValue(iV0X, iV1Y, iV1Z) * (fA * ((T)(1) - fB) * ((T)(1) - fC));
  // 1 1 1
  fValue += getValue(iV1X, iV1Y, iV1Z) * (((T)(1) - fA) * ((T)(1) - fB) * ((T)(1) - fC));

  return fValue;
};

/**
 * get the minimum density.
 * \return the minimum density
 */
template<class T>
inline T svt_volume<T>::getMinDensity()
{
  if (m_bChanged)
    calcMinMaxDensity();

  return m_fMinDensity;
};
/**
 * get the maximum density.
 *  \return the maximum density
 */
template<class T>
inline T svt_volume<T>::getMaxDensity()
{
  if (m_bChanged)
    calcMinMaxDensity();

  return m_fMaxDensity;
};
/**
 * get the average density.
 *  \return the average density
 */
template<class T>
inline T svt_volume<T>::getAvgDensity()
{
  if (m_bChanged)
    calcMinMaxDensity();

  return m_fAvgDensity;
};

/**
 * get norm of the volume = sqrt( sum( at(i) * at(i) ) );
 * \return norm
 */
template<class T>
inline Real64 svt_volume<T>::getNorm()
{
  if (m_bChanged)
    calcMinMaxDensity();

  return m_fNorm;
};

/**
 * set voxel width (in Angstroem)
 * \param fWidth voxel width
 */
template<class T>
inline void svt_volume<T>::setWidth(Real64 fWidth)
{
  m_fWidth = fWidth;
};
/**
 * get voxel width (in Angstroem)
 * \return voxel width
 */
template<class T>
inline Real64 svt_volume<T>::getWidth() const
{
  return m_fWidth;
};

/**
 * normalize volume
 */
template<class T>
inline void svt_volume<T>::normalize()
{
  unsigned long iNVox, iCount;
  iNVox = m_iSizeX * m_iSizeY * m_iSizeZ;

  T fMax = getMaxDensity();
  T fMin = getMinDensity();
  T fLength = fMax - fMin;

  if (fLength == (T)(0.0))
    return;

  for (iCount = 0; iCount < iNVox; iCount++) {
    m_pData[iCount] = (m_pData[iCount] - fMin) / fLength;
  };

  m_fMinDensity = (T)(0.0);
  m_fMaxDensity = (T)(1.0);

  m_bChanged = true;
};

/**
 * normalize volume
 * \param fSigma1 - sigma used in the gaussian calculation
 * \param fSigma2 - second  sigma used in the gaussian calculation
 * \param fSigma -  blurring kernel following local normalization (default fSigma=1.0 - don't blur)
 */
template<class T>
inline void svt_volume<T>::locallyNormalize(T fSigma, bool bProgress)
{
  svt_volume<T> oGaussian;

  if (fSigma <= 2.5)
    oGaussian.createGaussian(fSigma, fSigma);
  else
    oGaussian.createGaussian(fSigma, 2.0);

  SVTLBBO << "Applying local normalization using a gaussian of size (" << oGaussian.getSizeX() << " x " << oGaussian.getSizeY() << " x " << oGaussian.getSizeZ() << ") " <<  endl;

  svt_volume<T> oAverage, oSd;
  oAverage    = *this;
  oSd         = *this;

  //compute the average
  oAverage.convolve(oGaussian, bProgress);
  oAverage.scale(-1.0);



  *this += oAverage;

  //compute the standard deviation
  unsigned int iNum = oSd.getSizeX() * oSd.getSizeY() * oSd.getSizeZ();
  for (unsigned int i = 0; i < iNum; i++)
    oSd.setAt(i, oSd.at(i) * oSd.at(i));

  oSd.convolve(oGaussian, bProgress);

  for (unsigned int i = 0; i < iNum; i++) {
    oSd.setAt(i, sqrt(oSd.at(i)  - oAverage.at(i)*oAverage.at(i)));
    if (oSd.at(i) != 0)
      setAt(i, at(i) / oSd.at(i));
  }

};

/**
 * normalize volume
 * \param fSigma -  blurring kernel following local normalization (default fSigma=1.0 - don't blur)
 */
template<class T>
inline void svt_volume<T>::locallyNormalizeCorrectBorders(T fSigma, bool bProgress)
{

  svt_volume<T> oGaussian;

  if (fSigma <= 2.5)
    oGaussian.createGaussian(fSigma, fSigma);
  else
    oGaussian.createGaussian(fSigma, 2.0);

  //fsigma = 1 - the map remains the same
  if (oGaussian.getSizeX() == 1 || oGaussian.getSizeY() == 1 || oGaussian.getSizeZ() == 1)
    return;

  SVTLBBO << "Applying local normalization using a gaussian of size (" << oGaussian.getSizeX() << " x " << oGaussian.getSizeY() << " x " << oGaussian.getSizeZ() << ") " <<  endl;

  Real64 fHW = int (oGaussian.getSizeX() / 2.0);

#ifdef _OPENMP
  int iThreads = omp_get_max_threads();
  omp_set_num_threads(iThreads);
  SVTLBBO << "Starting parallel local normalization on " << iThreads << " cores" <<  endl;
#endif

  svt_volume<Real64> oLocNormVol;
  oLocNormVol = *this;
  Real64 fCount;
  Real64 fAvg, fSqAvg, fVal, fW;

#ifdef _OPENMP
  #pragma omp parallel for \
  shared(oLocNormVol, oGaussian) \
  private( fCount, fAvg, fSqAvg, fVal, fW ) \
  schedule(dynamic, 1)
#endif

  for (unsigned int iIndexZ = 0; iIndexZ < m_iSizeZ; iIndexZ++) {
    //run block in separate threads
    for (unsigned int iIndexX = 0; iIndexX < m_iSizeX; iIndexX++) {
      for (unsigned int iIndexY = 0; iIndexY < m_iSizeY; iIndexY++) {
#ifdef _OPENMP
        #pragma omp critical
#endif
        //process in window for voxel iIndexX, iIndexY, iIndexZ
        fCount = 0;
        fAvg = 0;
        fSqAvg = 0;
        for (int i = iIndexX - fHW; i <= iIndexX + fHW; i++) {
          for (int j = iIndexY - fHW; j <= iIndexY + fHW; j++) {
            for (int k = iIndexZ - fHW; k <= iIndexZ + fHW; k++)
              if (i >= 0 && i < m_iSizeX && j >= 0 && j < m_iSizeY && k >= 0 && k < m_iSizeZ) {
                fW = oGaussian.at(i - iIndexX + fHW, j - iIndexY + fHW, k - iIndexZ + fHW);
                fVal = at(i, j, k);
                fAvg += fVal * fW;
                fSqAvg += fVal * fVal * fW;
                fCount += fW;
              }
          }
        }
        fAvg /= (Real64)fCount;
        fSqAvg /= (Real64)fCount;
        if (abs(sqrt(fSqAvg - fAvg * fAvg)) > 1e-6)
          oLocNormVol.setAt(iIndexX, iIndexY, iIndexZ, (at(iIndexX, iIndexY, iIndexZ) - fAvg) / sqrt(fSqAvg - fAvg * fAvg)); //the sd
        else
          oLocNormVol.setAt(iIndexX, iIndexY, iIndexZ, (at(iIndexX, iIndexY, iIndexZ) - fAvg)); //the sd
      }
    }
  }

  *this = oLocNormVol;

};


/**
 * Threshold volume. All voxel values above a threshold value are cutoff and set to the threshold.
 * \param fMaxDensity new maximum density
 * \param fMinDensity new minimum density
 */
template<class T>
inline void svt_volume<T>::threshold(T fMinDensity, T fMaxDensity)
{
  unsigned long iNVox, iCount;
  iNVox = m_iSizeX * m_iSizeY * m_iSizeZ;

  for (iCount = 0; iCount < iNVox; iCount++) {
    if (m_pData[iCount] > fMaxDensity)
      m_pData[iCount] = fMaxDensity;

    if (m_pData[iCount] < fMinDensity)
      m_pData[iCount] = fMinDensity;
  };

  m_fMinDensity = (T)(fMinDensity);
  m_fMaxDensity = (T)(fMaxDensity);

  m_bChanged = true;
};

/**
 * Crop volume - cuts out and returns a subvolume
 * \param iNewMinX new minimum x voxel index
 * \param iNewMaxX new maximal x voxel index
 * \param iNewMinY new minimum y voxel index
 * \param iNewMaxY new maximal y voxel index
 * \param iNewMinZ new minimum z voxel index
 * \param iNewMaxZ new maximal z voxel index
 */
template<class T>
inline svt_volume<T> &svt_volume<T>::crop(
  unsigned int iNewMinX, unsigned int iNewMaxX,
  unsigned int iNewMinY, unsigned int iNewMaxY,
  unsigned int iNewMinZ, unsigned int iNewMaxZ
)
{
  svt_volume<T> oVol;
  oVol.allocate(
    iNewMaxX - iNewMinX + 1,
    iNewMaxY - iNewMinY + 1,
    iNewMaxZ - iNewMinZ + 1
  );

  for (unsigned int iZ = iNewMinZ; iZ <= iNewMaxZ; iZ++)
    for (unsigned int iY = iNewMinY; iY <= iNewMaxY; iY++)
      for (unsigned int iX = iNewMinX; iX <= iNewMaxX; iX++)
        oVol.setAt(iX - iNewMinX, iY - iNewMinY, iZ - iNewMinZ, this->at(iX, iY, iZ));

  Real64 fNewGridX = m_fGridX + (iNewMinX * m_fWidth);
  Real64 fNewGridY = m_fGridY + (iNewMinY * m_fWidth);
  Real64 fNewGridZ = m_fGridZ + (iNewMinZ * m_fWidth);

  // delete old data
  if (m_bInternal)
    delete[] m_pData;

  m_pData = oVol.getData();
  m_bInternal = true;

  m_fGridX = fNewGridX;
  m_fGridY = fNewGridY;
  m_fGridZ = fNewGridZ;

  m_iSizeX = oVol.getSizeX();
  m_iSizeY = oVol.getSizeY();
  m_iSizeZ = oVol.getSizeZ();

  m_bChanged = true;

  oVol.m_bInternal = false;

  return *this;
};

/**
 * Delete a spherical subregion
 * \param fCenterX center voxel coordinate for the spherical region
 * \param fCenterY center voxel coordinate for the spherical region
 * \param fCenterZ center voxel coordinate for the spherical region
 * \param fRadius radius of the sphere
 * \return returns a new svt_volume object with the spherical region cut out
 */
template<class T>
inline svt_volume<T> &svt_volume<T>::cutSphere(Real64 fCenterX, Real64 fCenterY, Real64 fCenterZ, Real64 fRadius)
{
  svt_volume<T> *pVol = new svt_volume<T>();
  pVol->deepCopy(*this);
  Real64 fX, fY, fZ;
  Real64 fRadSq = fRadius * fRadius;

  for (unsigned int iZ = 0; iZ < m_iSizeZ; iZ++) {
    fZ = (Real64)(iZ);

    for (unsigned int iY = 0; iY < m_iSizeY; iY++) {
      fY = (Real64)(iY);

      for (unsigned int iX = 0; iX < m_iSizeX; iX++) {
        fX = (Real64)(iX);

        if (

          ((fX - fCenterX) * (fX - fCenterX)) +
          ((fY - fCenterY) * (fY - fCenterY)) +
          ((fZ - fCenterZ) * (fZ - fCenterZ))

          < fRadSq)

          pVol->setAt(iX, iY, iZ, 0.0);
      }
    }
  }

  return *pVol;
};

/**
 * copy a spherical subregion
 * \param fCenterX center voxel coordinate for the spherical region
 * \param fCenterY center voxel coordinate for the spherical region
 * \param fCenterZ center voxel coordinate for the spherical region
 * \param fRadius radius of the sphere
 * \return returns a new svt_volume object with the cutted out spherical region
 */
template<class T>
inline svt_volume<T> &svt_volume<T>::copySphere(Real64 fCenterX, Real64 fCenterY, Real64 fCenterZ, Real64 fRadius)
{
  svt_volume<T> *pVol = new svt_volume<T>();
  pVol->deepCopy(*this);

  Real64 fDist;

  svt_vector4<Real64> oCenter;
  svt_vector4<Real64> oVoxel;

  oCenter.x((fCenterX - m_fGridX) / m_fWidth);
  oCenter.y((fCenterY - m_fGridY) / m_fWidth);
  oCenter.z((fCenterZ - m_fGridZ) / m_fWidth);


  for (unsigned int iZ = 0; iZ < m_iSizeZ; iZ++) {
    for (unsigned int iY = 0; iY < m_iSizeY; iY++) {
      for (unsigned int iX = 0; iX < m_iSizeX; iX++) {
        oVoxel.x((Real64) iX);
        oVoxel.y((Real64) iY);
        oVoxel.z((Real64) iZ);
        fDist = oVoxel.distance(oCenter);

        if (fDist > (fRadius / m_fWidth))
          pVol->setAt(iX, iY, iZ, 0.0);
      }
    }
  }

  unsigned int iMinX = 0;
  unsigned int iMaxX = m_iSizeX - 1;

  unsigned int iMinY = 0;
  unsigned int iMaxY = m_iSizeY - 1;

  unsigned int iMinZ = 0;
  unsigned int iMaxZ = m_iSizeZ - 1;


  if (floor(fCenterX + m_fGridX * m_fWidth - fRadius) >= 0 &&  floor(fCenterX + m_fGridX * m_fWidth - fRadius) < m_iSizeX)
    iMinX =  floor((fCenterX - fRadius - m_fGridX) / m_fWidth) ;

  if (ceil(fCenterX + m_fGridX * m_fWidth + fRadius) >= 0 &&  ceil(fCenterX + m_fGridX * m_fWidth + fRadius) < m_iSizeX)
    iMaxX =  ceil((fCenterX + fRadius - m_fGridX) / m_fWidth)  ;

  if (floor(fCenterY + m_fGridY * m_fWidth - fRadius) >= 0 &&  floor(fCenterY + m_fGridY * m_fWidth - fRadius) < m_iSizeY)
    iMinY =  floor((fCenterY - fRadius - m_fGridY) / m_fWidth) ;

  if (ceil(fCenterY + m_fGridY * m_fWidth + fRadius) >= 0 &&  ceil(fCenterY + m_fGridY * m_fWidth + fRadius) < m_iSizeY)
    iMaxY =  ceil((fCenterY + fRadius - m_fGridY) / m_fWidth)  ;

  if (floor(fCenterZ + m_fGridZ * m_fWidth - fRadius) >= 0 &&  floor(fCenterZ + m_fGridZ * m_fWidth - fRadius) < m_iSizeZ)
    iMinZ =  floor((fCenterZ - fRadius - m_fGridZ) / m_fWidth) ;

  if (ceil(fCenterZ + m_fGridZ * m_fWidth + fRadius) >= 0 &&  ceil(fCenterZ + m_fGridZ * m_fWidth + fRadius) < m_iSizeZ)
    iMaxZ =  ceil((fCenterZ + fRadius - m_fGridZ) / m_fWidth)  ;

  pVol->crop(iMinX, iMaxX, iMinY, iMaxY, iMinZ, iMaxZ);

  return *pVol;
};


/**
 * Mask with another svt_volume object - all the voxels in this vol are multiplied with the mask volume voxels (multiplied by 0 for not overlapping voxels).
 * \param rMask reference to other object
 */
template<class T>
svt_volume<T> &svt_volume<T>::mask(svt_volume<T> &rMask)
{
  Real64 fX, fY, fZ;
  T fVal;

  if (abs(m_fWidth - rMask.getWidth()) > 1e-6) { // || fmod((m_fGridX - rMask.getGridX()) , m_fWidth) != 0.0 || fmod((m_fGridY - rMask.getGridY()), m_fWidth) != 0.0 || fmod((m_fGridZ - rMask.getGridZ()), m_fWidth) != 0.0)
    //if (m_fWidth != rMask.getWidth())
    //{
    SVTLBBO << "mask(): Interpolation necessary as voxel width is different!" << endl;
    SVTLBBO << "        Voxel width: " << m_fWidth << " Mask: " << rMask.getWidth() << endl;
    //}
    //else
    //{
    //    SVTLBBO << "mask(): Interpolation necessary as origin does not lie on the same lattice!" << endl;
    //    SVTLBBO << "        Origin: " << m_fGridX << ", " << m_fGridY << ", " << m_fGridZ << " Mask: " << rMask.getGridX() << ", " << rMask.getGridY() << ", " << rMask.getGridZ() << endl;
    //    SVTLBBO << "DEBUG: " << fmod((m_fGridX - rMask.getGridX()) , m_fWidth) << ", " << fmod((m_fGridY - rMask.getGridY()) , m_fWidth) << ", " << fmod((m_fGridZ - rMask.getGridZ()) , m_fWidth) << endl;
    //}

    for (unsigned int iZ = 0; iZ < m_iSizeZ; iZ++) {
      fZ = (Real64)(iZ);
      fZ *= m_fWidth;
      fZ += m_fGridZ;

      for (unsigned int iY = 0; iY < m_iSizeY; iY++) {
        fY = (Real64)(iY);
        fY *= m_fWidth;
        fY += m_fGridY;

        for (unsigned int iX = 0; iX < m_iSizeX; iX++) {
          fX = (Real64)(iX);
          fX *= m_fWidth;
          fX += m_fGridX;

          fVal = this->at(iX, iY, iZ) * rMask.getRealSpaceValue(fX, fY, fZ);
          this->setAt(iX, iY, iZ, fVal);
        }
      }
    }

    shrinkToOccupiedVolume();

  } else {

    int iStartX = (int)((rMask.getGridX() - m_fGridX) / m_fWidth);
    int iStartY = (int)((rMask.getGridY() - m_fGridY) / m_fWidth);
    int iStartZ = (int)((rMask.getGridZ() - m_fGridZ) / m_fWidth);

    int iMinX = (int)m_iSizeX;
    int iMaxX = 0;
    int iMinY = (int)m_iSizeY;
    int iMaxY = 0;
    int iMinZ = (int)m_iSizeZ;
    int iMaxZ = 0;

    for (int iZ = 0; iZ < (int)m_iSizeZ; iZ++) {
      for (int iY = 0; iY < (int)m_iSizeY; iY++) {
        for (int iX = 0; iX < (int)m_iSizeX; iX++) {
          if (iZ < iStartZ || iY < iStartY || iX < iStartX || iX - iStartX >= (int)rMask.getSizeX() || iY - iStartY >= (int)rMask.getSizeY() || iZ - iStartZ >= (int)rMask.getSizeZ()) {
            fVal = 0.0;
            this->setAt(iX, iY, iZ, 0.0);
          } else {
            fVal = this->at(iX, iY, iZ) * rMask.getValue(iX - iStartX, iY - iStartY, iZ - iStartZ);
            this->setAt(iX, iY, iZ, fVal);

            if (fVal != 0.0) {
              if (iX > iMaxX)
                iMaxX = iX;
              if (iX < iMinX)
                iMinX = iX;

              if (iY > iMaxY)
                iMaxY = iY;
              if (iY < iMinY)
                iMinY = iY;

              if (iZ > iMaxZ)
                iMaxZ = iZ;
              if (iZ < iMinZ)
                iMinZ = iZ;
            }
          }
        }
      }
    }

    iMinX = iMinX - 1 >= 0 ? iMinX - 1 : 0;
    iMinY = iMinY - 1 >= 0 ? iMinY - 1 : 0;
    iMinZ = iMinZ - 1 >= 0 ? iMinZ - 1 : 0;

    iMaxX = iMaxX + 1 < (int)getSizeX() ? iMaxX + 1 : iMaxX;
    iMaxY = iMaxY + 1 < (int)getSizeY() ? iMaxY + 1 : iMaxY;
    iMaxZ = iMaxZ + 1 < (int)getSizeZ() ? iMaxZ + 1 : iMaxZ;

    SVTLBBO << "New volume size: " << iMaxX - iMinX << " x " << iMaxY - iMinY << " x " << iMaxZ - iMinZ << endl;
    SVTLBBO << "New volume size: " << iMaxX  << " " <<   iMinX << " x " << iMaxY  << " " <<  iMinY << " x " << iMaxZ  << " " <<  iMinZ << endl;
    crop(iMinX, iMaxX, iMinY, iMaxY, iMinZ, iMaxZ);
  }

  m_bChanged = true;

  return *this;
};

/**
 * Make a mask volume - if voxel > threshold it is set to 1, otherwise 0
 * \param fThreshold threshold value
 */
template<class T>
svt_volume<T> &svt_volume<T>::makeMask(T fThreshold)
{
  unsigned int iNum = m_iSizeX * m_iSizeY * m_iSizeZ;

  for (unsigned int i = 0; i < iNum; i++) {
    if (this->at(i) > fThreshold)
      this->setAt(i, 1.0);
    else
      this->setAt(i, 0.0);
  }

  return *this;
};

/**
 * Invert mask - all voxels have to be between 0 and 1. The voxel values are inverted.
 */
template<class T>
svt_volume<T> &svt_volume<T>::invertMask()
{
  unsigned int iNum = m_iSizeX * m_iSizeY * m_iSizeZ;
  T fValue;

  for (unsigned int i = 0; i < iNum; i++) {
    fValue = this->at(i);

    this->setAt(i, (T)(1.0) - fValue);
  }

  return *this;
};

/**
 * shrinks the original volume to the occupied Volume - redefine m_fGrid and m_iSize such that volume is not padded by zeros
 */
template<class T>
void svt_volume<T>::shrinkToOccupiedVolume()
{
  unsigned int iNewMinX = 0;
  unsigned int iNewMaxX = m_iSizeX;
  unsigned int iNewMinY = 0;
  unsigned int iNewMaxY = m_iSizeY;
  unsigned int iNewMinZ = 0;
  unsigned int iNewMaxZ = m_iSizeZ;

  iNewMinX = getNewBorder(0           , 1         , 0         , m_iSizeY  , 0         , m_iSizeZ);
  iNewMaxX = getNewBorder(m_iSizeX - 1  , m_iSizeX  , 0         , m_iSizeY  , 0         , m_iSizeZ);

  iNewMinY = getNewBorder(0           , m_iSizeX  , 0         , 1         , 0         , m_iSizeZ);
  iNewMaxY = getNewBorder(0           , m_iSizeX  , m_iSizeY - 1, m_iSizeY  , 0         , m_iSizeZ);

  iNewMinZ = getNewBorder(0           , m_iSizeX  , 0         , m_iSizeY  , 0         , 1);
  iNewMaxZ = getNewBorder(0           , m_iSizeX  , 0         , m_iSizeY  , m_iSizeZ - 1, m_iSizeZ);

  crop(iNewMinX, iNewMaxX, iNewMinY, iNewMaxY, iNewMinZ, iNewMaxZ);
}

/**
 * get the new border
 * \param iMinX the min value of X
 * \param iMaxX the max value of X
 * \param iMinY the min value of Y
 * \param iMaxY the max value of Y
 * \param iMinZ the min value of Z
 * \param iMaxZ the max value of Z
 * \return the new index (voxel value)
 */
template<class T>
unsigned int svt_volume<T>::getNewBorder(
  unsigned int iMinX, unsigned int iMaxX,
  unsigned int iMinY, unsigned int iMaxY,
  unsigned int iMinZ, unsigned int iMaxZ
)
{
  bool bIsNull;
  unsigned int iNewBorder = 0;
  int iIncrement = +1;

  //should go downwards
  if ((iMinX + 1 == iMaxX && iMinX + 1 == m_iSizeX) || (iMinY + 1 == iMaxY && iMinY + 1 == m_iSizeY) || (iMinZ + 1 == iMaxZ && iMinZ + 1 == m_iSizeZ))
    iIncrement = -1;

  do {
    bIsNull = true;
    for (unsigned int iX = iMinX; iX < iMaxX; iX++)
      for (unsigned int iY = iMinY; iY < iMaxY; iY++)
        for (unsigned int iZ = iMinZ; iZ < iMaxZ; iZ++)
          if (this->getValue(iX, iY, iZ) != 0.0)
            bIsNull = false;

    if (iMinX + 1 == iMaxX) {
      iNewBorder = iMinX;
      iMinX += iIncrement;
      iMaxX += iIncrement;
    }

    if (iMinY + 1 == iMaxY) {
      iNewBorder = iMinY;
      iMinY += iIncrement;
      iMaxY += iIncrement;
    }

    if (iMinZ + 1 == iMaxZ) {
      iNewBorder = iMinZ;
      iMinZ += iIncrement;
      iMaxZ += iIncrement;
    }
    //SVTLBBO << "DEBUG: " <<  iMinX  << "-" << iMaxX  << " " <<  iMinY  << "-" << iMaxY  <<  " " << iMinZ  << "-" << iMaxZ << " " << iNewBorder << endl;
  } while (bIsNull && (int)iMinX >= 0 && iMaxX <= m_iSizeX && (int)iMinY >= 0 && iMaxY <= m_iSizeY && (int)iMinZ >= 0 && iMaxZ <= m_iSizeZ);

  return iNewBorder;
}

/**
 * change voxel spacings and bring map origin into register with coordinate system origin
 */
template<class T>
void svt_volume<T>::interpolate_map(Real64 fWidthX, Real64 fWidthY, Real64 fWidthZ)
{
  if (fWidthY == 0)
    fWidthY = fWidthX;
  if (fWidthZ == 0)
    fWidthZ = fWidthX;

  Real64 fDeltaX, fDeltaY, fDeltaZ;
  Real64 fXpos, fYpos, fZpos, fGx, fGy, fGz, fA, fB, fC;
  int iX0, iY0, iZ0, iX1, iY1, iZ1;
  int iSx, iSy, iSz, iEx, iEy, iEz;

  /* output start index rel. to coordinate system origin, asserting that outmap is fully embedded in inmap */
  iSx = ceil(m_fGridX / fWidthX);
  iSy = ceil(m_fGridY / fWidthY);
  iSz = ceil(m_fGridZ / fWidthZ);

  /* output end index rel. to coordinate system origin, asserting that outmap is fully embedded in inmap */
  iEx = floor((m_fGridX + m_fWidth * (m_iSizeX - 1)) / fWidthX);
  iEy = floor((m_fGridY + m_fWidth * (m_iSizeY - 1)) / fWidthY);
  iEz = floor((m_fGridZ + m_fWidth * (m_iSizeZ - 1)) / fWidthZ);

  /* assign output grid size */
  int iSizeX = iEx - iSx + 1;
  int iSizeY = iEy - iSy + 1;
  int iSizeZ = iEz - iSz + 1;
  if (iSizeX < 2 || iSizeY < 2 || iSizeZ < 2) {
    error_sba(85010, "Interpolation output map size underflow!");
    exit(1);
  }

  svt_volume<T> oVol(iSizeX, iSizeY, iSizeZ);
  oVol.setWidth(fWidthX);

  /* save origin shift */
  fDeltaX = iSx * fWidthX - m_fGridX;
  fDeltaY = iSy * fWidthY - m_fGridY;
  fDeltaZ = iSz * fWidthZ - m_fGridZ;

  Real64 fNewDensity;
  for (int iZ = 0; iZ < iSizeZ; iZ++) {
    for (int iY = 0; iY < iSizeY; iY++) {
      for (int iX = 0; iX < iSizeX; iX++) {

        /* determine position of outmap voxel relative to start of inmap */
        fXpos = fDeltaX + iX * fWidthX;
        fYpos = fDeltaY + iY * fWidthY;
        fZpos = fDeltaZ + iZ * fWidthZ;

        /* compute position in inmap voxel units */
        fGx = (fXpos / m_fWidth);
        fGy = (fYpos / m_fWidth);
        fGz = (fZpos / m_fWidth);

        /* compute bounding box voxel indices and linear distances */
        iX0 = floor(fGx);
        iY0 = floor(fGy);
        iZ0 = floor(fGz);
        iX1 = ceil(fGx);
        iY1 = ceil(fGy);
        iZ1 = ceil(fGz);
        fA = fGx - iX0;
        fB = fGy - iY0;
        fC = fGz - iZ0;

        /* interpolate */
        fNewDensity =
          fA    * fB    * fC    * this->getValue(iX1, iY1, iZ1) +
          (1 - fA)  * fB    * fC    * this->getValue(iX0, iY1, iZ1) +
          fA    * (1 - fB)  * fC    * this->getValue(iX1, iY0, iZ1) +
          fA    * fB    * (1 - fC)  * this->getValue(iX1, iY1, iZ0) +
          fA    * (1 - fB)  * (1 - fC)  * this->getValue(iX1, iY0, iZ0) +
          (1 - fA)  * fB    * (1 - fC)  * this->getValue(iX0, iY1, iZ0) +
          (1 - fA)  * (1 - fB)  * fC    * this->getValue(iX0, iY0, iZ1) +
          (1 - fA)  * (1 - fB)  * (1 - fC)  * this->getValue(iX0, iY0, iZ0);

        oVol.setAt(iX, iY, iZ, fNewDensity);
      }
    }
  }

  oVol.setGrid(m_fGridX + fDeltaX, m_fGridY + fDeltaY , m_fGridZ + fDeltaZ);
  *this = oVol;

  m_bChanged = true;
}


/**
 * Calc gradient map.
 * \param eMode gradient method selector: CentralDistance or Sobel
 * \param rGradientX reference to an svt_volume object for the gradient map in x direction
 * \param rGradientY reference to an svt_volume object for the gradient map in y direction
 * \param rGradientZ reference to an svt_volume object for the gradient map in z direction
 */
template<class T>
inline void svt_volume<T>::calcGradient(const GradientMode eMode, svt_volume<T> &rGradientX, svt_volume<T> &rGradientY, svt_volume<T> &rGradientZ) const
{
  rGradientX.allocate(m_iSizeX, m_iSizeY, m_iSizeZ);
  rGradientY.allocate(m_iSizeX, m_iSizeY, m_iSizeZ);
  rGradientZ.allocate(m_iSizeX, m_iSizeY, m_iSizeZ);

  unsigned int iX, iY, iZ;
  T fValX, fValY, fValZ;

  int iTime = svt_getElapsedTime();
  SVTLBBO << "Gradient calculation..." << endl;

  // central difference gradient
  if (eMode == CentralDistance) {
    for (iZ = 0; iZ <= m_iSizeZ; iZ++)
      for (iY = 0; iY <= m_iSizeY; iY++)
        for (iX = 0; iX <= m_iSizeX; iX++) {
          if (iX != 0 && iX != m_iSizeX)
            fValX = (getValue(iX + 1, iY, iZ) - getValue(iX - 1, iY, iZ)) / 2;
          else
            fValX = 0.0f;
          rGradientX.setValue(iX, iY, iZ, fValX);

          if (iY != 0 && iY != m_iSizeY)
            fValY = (getValue(iX, iY + 1, iZ) - getValue(iX, iY - 1, iZ)) / 2;
          else
            fValY = 0.0f;
          rGradientY.setValue(iX, iY, iZ, fValY);

          if (iZ != 0 && iZ != m_iSizeZ)
            fValZ = (getValue(iX, iY, iZ + 1) - getValue(iX, iY, iZ - 1)) / 2;
          else
            fValZ = 0.0f;
          rGradientZ.setValue(iX, iY, iZ, fValZ);
        }

  }

  // sobel
  if (eMode == Sobel) {
    svt_vector4<Real64> oVec;
    int iKX, iKY, iKZ;
    // sobel kernels (fastest running index = z, then y and then x!!)
    Real64 xKernel[27] = {
      1, 3, 1,
      3, 6, 3,
      1, 3, 1,
      0, 0, 0,
      0, 0, 0,
      0, 0, 0,
      -1, -3, -1,
      -3, -6, -3,
      -1, -3, -1
    };
    Real64 yKernel[27] = {
      1, 3, 1,
      0, 0, 0,
      -1, -3, -1,
      3, 6, 3,
      0, 0, 0,
      -3, -6, -3,
      1, 3, 1,
      0, 0, 0,
      -1, -3, -1
    };
    Real64 zKernel[27] = {
      1, 0, -1,
      3, 0, -3,
      1, 0, -1,
      3, 0, -3,
      6, 0, -6,
      3, 0, -3,
      1, 0, -1,
      3, 0, -3,
      1, 0, -1
    };

    for (iZ = 0; iZ < m_iSizeZ; iZ++) {
      for (iX = 0; iX < m_iSizeX; iX++) {
        for (iY = 0; iY < m_iSizeY; iY++) {
          // apply x kernel
          Real64 xVal = 0.0;
          for (iKZ = -1; iKZ < 2; iKZ++)
            for (iKY = -1; iKY < 2; iKY++)
              for (iKX = -1; iKX < 2; iKX++) {
                if ((((int)(iX) + iKX) < (int)(m_iSizeX) && ((int)(iX) + iKX) >= 0) && (((int)(iY) + iKY) < (int)(m_iSizeY) && ((int)(iY) + iKY) >= 0) && (((int)(iZ) + iKZ) < (int)(m_iSizeZ) && ((int)(iZ) + iKZ) >= 0))
                  xVal += xKernel[((iKX + 1) * 9) + ((iKY + 1) * 3) + iKZ + 1] * Real64(m_pData[(iX + iKX) + ((iY + iKY) * m_iSizeX) + ((iZ + iKZ) * m_iSizeX * m_iSizeY)]);
              }

          // apply y kernel
          Real64 yVal = 0.0;
          for (iKZ = -1; iKZ < 2; iKZ++)
            for (iKY = -1; iKY < 2; iKY++)
              for (iKX = -1; iKX < 2; iKX++) {
                if ((((int)(iX) + iKX) < (int)(m_iSizeX) && ((int)(iX) + iKX) >= 0) && (((int)(iY) + iKY) < (int)(m_iSizeY) && ((int)(iY) + iKY) >= 0) && (((int)(iZ) + iKZ) < (int)(m_iSizeZ) && ((int)(iZ) + iKZ) >= 0))
                  yVal += yKernel[((iKX + 1) * 9) + ((iKY + 1) * 3) + iKZ + 1] * Real64(m_pData[(iX + iKX) + ((iY + iKY) * m_iSizeX) + ((iZ + iKZ) * m_iSizeX * m_iSizeY)]);
              }

          // apply z kernel
          Real64 zVal = 0.0;
          for (iKX = -1; iKX < 2; iKX++)
            for (iKY = -1; iKY < 2; iKY++)
              for (iKZ = -1; iKZ < 2; iKZ++) {
                if ((((int)(iX) + iKX) < (int)(m_iSizeX) && ((int)(iX) + iKX) >= 0) && (((int)(iY) + iKY) < (int)(m_iSizeY) && ((int)(iY) + iKY) >= 0) && (((int)(iZ) + iKZ) < (int)(m_iSizeZ) && ((int)(iZ) + iKZ) >= 0))
                  zVal += zKernel[((iKX + 1) * 9) + ((iKY + 1) * 3) + iKZ + 1] * Real64(m_pData[(iX + iKX) + ((iY + iKY) * m_iSizeX) + ((iZ + iKZ) * m_iSizeX * m_iSizeY)]);
              }

          oVec.x(xVal);
          oVec.y(yVal);
          oVec.z(zVal);
          oVec.w(0.0);

          if (iX == 0 || iY == 0 || iZ == 0 || iX == m_iSizeX - 1 || iY == m_iSizeY - 1 || iZ == m_iSizeZ - 1) {
            oVec.x(1.0);
            oVec.y(1.0);
            oVec.z(1.0);
          }

          rGradientX.setValue(iX, iY, iZ, oVec.x());
          rGradientY.setValue(iX, iY, iZ, oVec.y());
          rGradientZ.setValue(iX, iY, iZ, oVec.z());
        }
      }
    }
  }

  iTime = svt_getElapsedTime() - iTime;
  SVTLBBO << "Gradient runtime: " << iTime << endl;
};

/**
 * Set densities below limit to zero
 * \param fLimit the threshold
 */
template<class T>
inline void svt_volume<T>::threshold(T fLimit)
{
  unsigned int i;
  unsigned int iNVox = m_iSizeX * m_iSizeY * m_iSizeZ;

  for (i = 0; i < iNVox; i++)
    if (m_pData[i] < fLimit)
      m_pData[i] = (T)(0.0);

  m_bChanged = true;
};

/**
 * Convolve this volume with another one, which is a 1D volume (only x axis). The 1D kernel will get convolved in all three directions, so this function should be
 * used with linear separable kernels.
 * \param rKernel reference to kernel volume
 * \param bNormalize normalizes during convolution
 * \param bProgress show a progress bar
 */
template<class T>
inline void svt_volume<T>::convolve1D3D(svt_volume<T> &rKernel, bool bProgress)
{
  unsigned int i, iX, iY, iZ;
  int iKX;
  T fVal;

  // allocate memory for temporary volume
  svt_volume<T> oTmp;
  oTmp.allocate(m_iSizeX, m_iSizeY, m_iSizeZ, 0.0);
  // set the current volume to changed
  m_bChanged = true;

  // calculate dimension of kernel volume
  int iDim = (int)((Real32)(rKernel.getSizeX()) * 0.5f);
  int iStart = -iDim;
  int iEnd = iDim + 1;

  try {
    // in x direction
    for (iZ = 0; iZ < m_iSizeZ; iZ++) {
      for (iX = 0; iX < m_iSizeX; iX++) {
        for (iY = 0; iY < m_iSizeY; iY++) {
          fVal = m_pData[(iX) + (iY * m_iSizeX) + (iZ * m_iSizeX * m_iSizeY)];

          if (fVal != (T)(0)) {
            for (iKX = iStart; iKX < iEnd; iKX++)
              if (((int)(iX) + iKX) < (int)(m_iSizeX) && ((int)(iX) + iKX) >= 0)
                oTmp.setAt(iX + iKX, iY, iZ, oTmp.at(iX + iKX, iY, iZ) + (rKernel.at(iKX - iStart) * fVal));
          }
        }
      }


    }

    // set target volume to 0
    for (i = 0; i < m_iSizeX * m_iSizeY * m_iSizeZ; i++)
      m_pData[i] = 0;

    // in y direction
    for (iZ = 0; iZ < m_iSizeZ; iZ++) {
      for (iX = 0; iX < m_iSizeX; iX++) {
        for (iY = 0; iY < m_iSizeY; iY++) {
          fVal = oTmp.at((iX) + (iY * m_iSizeX) + (iZ * m_iSizeX * m_iSizeY));

          if (fVal != (T)(0)) {
            for (iKX = iStart; iKX < iEnd; iKX++)
              if (((int)(iY) + iKX) < (int)(m_iSizeY) && ((int)(iY) + iKX) >= 0)
                this->setAt(iX, iY + iKX, iZ, this->at(iX, iY + iKX, iZ) + (rKernel.at(iKX - iStart) * fVal));
          }
        }
      }

    }

    oTmp.setValue(0.0);

    // in z direction
    for (iZ = 0; iZ < m_iSizeZ; iZ++) {
      for (iX = 0; iX < m_iSizeX; iX++) {
        for (iY = 0; iY < m_iSizeY; iY++) {
          fVal = this->at((iX) + (iY * m_iSizeX) + (iZ * m_iSizeX * m_iSizeY));

          if (fVal != (T)(0)) {
            for (iKX = iStart; iKX < iEnd; iKX++)
              if (((int)(iZ) + iKX) < (int)(m_iSizeZ) && ((int)(iZ) + iKX) >= 0)
                oTmp.setAt(iX, iY, iZ + iKX, oTmp.at(iX, iY, iZ + iKX) + (rKernel.at(iKX - iStart) * fVal));
          }
        }
      }
    }

    memcpy(this->m_pData, oTmp.getData(), sizeof(T)*m_iSizeX * m_iSizeY * m_iSizeZ);


  } catch (int e) {
  }
};

/**
 * Convolve this volume with another one (typically a 3D kernel filter)
 * \param rKernel reference to kernel volume
 * \param bNormalize normalizes during convolution
 * \param bProgress show a progress bar
 */
template<class T>
inline void svt_volume<T>::convolve(svt_volume<T> &rKernel, bool bProgress)
{
  unsigned int iSize = rKernel.getSizeX();
  unsigned int iSizeS = iSize * iSize;
  unsigned int iX, iY, iZ;
  int iKX, iKY, iKZ;
  Real64 fVal;

  Real64 fOldWidth = m_fWidth;
  Real64 fOldGridX = m_fGridX;
  Real64 fOldGridY = m_fGridY;
  Real64 fOldGridZ = m_fGridZ;

  T *pData = m_pData;
  bool bInternal = m_bInternal;
  m_pData = NULL;
  m_bInternal = false;

  svt_volume<T> oTmp;
  oTmp.allocate(m_iSizeX, m_iSizeY, m_iSizeZ, 0.0);

  m_pData = pData;
  bInternal = bInternal;
  m_bChanged = true;

  // calculate dimension of kernel volume
  int iDim = (int)((Real32)(rKernel.getSizeX()) * 0.5f);
  int iStart = -iDim;
  int iEnd = iDim + 1;

  try {

    for (iZ = 0; iZ < m_iSizeZ; iZ++) {
      for (iX = 0; iX < m_iSizeX; iX++) {
        for (iY = 0; iY < m_iSizeY; iY++) {
          fVal = m_pData[(iX) + (iY * m_iSizeX) + (iZ * m_iSizeX * m_iSizeY)];

          if (fVal != (T)(0)) {
            for (iKZ = iStart; iKZ < iEnd; iKZ++)
              for (iKY = iStart; iKY < iEnd; iKY++)
                for (iKX = iStart; iKX < iEnd; iKX++) {
                  if ((((int)(iX) + iKX) < (int)(m_iSizeX) && ((int)(iX) + iKX) >= 0) && (((int)(iY) + iKY) < (int)(m_iSizeY) && ((int)(iY) + iKY) >= 0) && (((int)(iZ) + iKZ) < (int)(m_iSizeZ) && ((int)(iZ) + iKZ) >= 0))
                    oTmp.setAt(iX + iKX, iY + iKY, iZ + iKZ, oTmp.at(iX + iKX, iY + iKY, iZ + iKZ) + (rKernel.at(((iKX - iStart)*iSizeS) + ((iKY - iStart)*iSize) + iKZ - iStart) * fVal));
                }
          }

        }
      }
    }

    *this = oTmp;

    m_bChanged = true;

    // store map properties of old map in oTmp
    m_fWidth = fOldWidth;
    m_fGridX = fOldGridX;
    m_fGridY = fOldGridY;
    m_fGridZ = fOldGridZ;

  } catch (int e) {
  }

};

/**
 * Create a Laplacian kernel volume
 */
template<class T>
void svt_volume<T>::createLaplacian()
{
  T aLap[3][3][3] = {
    { { 0,  0,      0},  {      0,  1 / 12.,      0},  { 0,      0, 0} },
    { { 0,  1 / 12.,  0},  {  1 / 12., -6 / 12., 1 / 12.0},  { 0,  1 / 12., 0} },
    { { 0,  0,      0},  {      0,  1 / 12.,      0},  { 0,      0, 0} }
  };

  deepCopy(3, 3, 3, (T *)aLap);
};

/**
 * Create a Gaussian kernel volume within SigmaFactor*fSigma
 * \param fSigma1D sigma of map
 * \param fSigmaFactor sigma factor
 */
template<class T>
void svt_volume<T>::createGaussian(double fSigma1D, double fSigmaFactor)
{
  // truncate at fSigmaFactor * fSigma1D
  unsigned int iSizeH = (int) ceil(fSigmaFactor * fSigma1D);
  unsigned int iSize = 2 * iSizeH - 1;

  allocate(iSize, iSize, iSize);

  // write Gaussian within fSigmaFactor * fSigma1D
  Real64 fBValue = -1 / (2.0 * fSigma1D * fSigma1D);
  Real64 fCValue = fSigmaFactor * fSigmaFactor * fSigma1D * fSigma1D;
  Real64 fScale = 0;
  Real64 fDSqu;

  unsigned int iX, iY, iZ;
  for (iZ = 0; iZ < iSize; iZ++)
    for (iY = 0; iY < iSize; iY++)
      for (iX = 0; iX < iSize; iX++) {
        fDSqu = (iX - iSizeH + 1) * (iX - iSizeH + 1) +
                (iY - iSizeH + 1) * (iY - iSizeH + 1) +
                (iZ - iSizeH + 1) * (iZ - iSizeH + 1);

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
 * Create a Gaussian blurring kernel volume (Situs scheme)
 * \param fWidth the voxel width of the target map one wants to convolve with the kernel
 * \param fResolution the target resolution
 * \param fVarp variance of map (if 0 no correction for lattice interpolation smoothing effects = default)
 * Attention: This will overwrite the current content of the volume object with the filter kernel. It will allocate the memory internally.
 */
template<class T>
void svt_volume<T>::createSitusBlurringKernel(Real64 fWidth, Real64 fResolution, Real64 fVarp)
{
  Real64 fSigma1 = fResolution / 2.0;
  Real64 fKmsd = fSigma1 * fSigma1 / (fWidth * fWidth);

  fResolution = 2.0 * sqrt((fSigma1 * fSigma1) + (fVarp * fWidth * fWidth));

  Real64 fVarmap = fKmsd;
  fVarmap -= fVarp;

  if (fVarmap < 0.0) {
    SVTLBBO << "Error: lattice smoothing exceeds kernel size" << endl;
    exit(1);
  }

  // sigma-1D
  Real64 fSigmaMap = sqrt(fVarmap / 3.0);

  // truncate at 3 sigma-1D == sqrt(3) sigma-3D
  unsigned int iSizeH = (int) ceil(3.0 * fSigmaMap);
  unsigned int iSize = 2 * iSizeH + 1;

  allocate(iSize, iSize, iSize);
  setValue(0.0);

  // kernel width
  setWidth(fWidth);

  // write Gaussian within 3 sigma-1D to map
  Real64 fBValue = -1.0 / (2.0 * fSigmaMap * fSigmaMap);
  Real64 fCValue = 9.0 * fSigmaMap * fSigmaMap;
  Real64 fScale = 0;
  Real64 fDSqu;

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
 * Create a one-dimensional Gaussian blurring kernel volume (Situs scheme)
 * \param fWidth the voxel width of the target map one wants to convolve with the kernel
 * \param fResolution the target resolution
 * \param fVarp variance of map (if 0 no correction for lattice interpolation smoothing effects = default)
 * Attention: This will overwrite the current content of the volume object with the filter kernel. It will allocate the memory internally.
 */
template<class T>
void svt_volume<T>::create1DBlurringKernel(Real64 fWidth, Real64 fResolution, Real64 fVarp)
{
  Real64 fSigma1 = fResolution / 2.0;
  Real64 fKmsd = fSigma1 * fSigma1 / (fWidth * fWidth);

  fResolution = 2.0 * sqrt((fSigma1 * fSigma1) + (fVarp * fWidth * fWidth));

  Real64 fVarmap = fKmsd;
  fVarmap -= fVarp;

  if (fVarmap < 0.0) {
    SVTLBBO << "Error: lattice smoothing exceeds kernel size" << endl;
    exit(1);
  }

  // sigma-1D
  Real64 fSigmaMap = sqrt(fVarmap / 3.0);

  // truncate at 3 sigma-1D == sqrt(3) sigma-3D
  unsigned int iSizeH = (int) ceil(3.0 * fSigmaMap);
  unsigned int iSize = 2 * iSizeH + 1;

  allocate(iSize, 1, 1);
  setValue(0.0);

  // kernel width
  setWidth(fWidth);

  // write Gaussian within 3 sigma-1D to map
  Real64 fBValue = -1.0 / (2.0 * fSigmaMap * fSigmaMap);
  Real64 fCValue = 9.0 * fSigmaMap * fSigmaMap;
  Real64 fScale = 0;
  Real64 fDSqu;

  unsigned int iX;
  for (iX = 0; iX < iSize; iX++) {
    fDSqu = (iX - iSizeH) * (iX - iSizeH);

    if (fDSqu <= fCValue)
      setValue(iX, 0, 0, exp(fDSqu * fBValue));

    fScale += getValue(iX, 0, 0);
  }

  for (iX = 0; iX < iSize; iX++)
    setValue(iX, 0, 0, getValue(iX, 0, 0) / fScale);
};

/**
 * Create a Laplacian of a Gaussian kernel volume.
 * \param fWidth the voxel width of the target map one wants to convolve with the kernel
 * \param fResolution the target resolution
 * \param fVarp variance of map (if 0 no correction for lattice interpolation smoothing effects = default)
 * Attention: This will overwrite the current content of the volume object with the filter kernel. It will allocate the memory internally.
 */
template<class T>
void svt_volume<T>::createLaplacianOfGaussian(Real64 fWidth, Real64 fResolution, Real64 fVarp)
{
  Real64 fSigma1 = fResolution / 2.0;
  Real64 fKmsd = fSigma1 * fSigma1 / (fWidth * fWidth);

  fResolution = 2.0 * sqrt((fSigma1 * fSigma1) + (fVarp * fWidth * fWidth));

  Real64 fVarmap = fKmsd;
  fVarmap -= fVarp;

  if (fVarmap < 0.0) {
    SVTLBBO << "Error: lattice smoothing exceeds kernel size" << endl;
    exit(1);
  }

  // sigma-1D
  Real64 fSigmaMap = sqrt(fVarmap / 3.0);

  // truncate at 3 sigma-1D == sqrt(3) sigma-3D
  unsigned int iSizeH = (int) ceil(3.0 * fSigmaMap);
  unsigned int iSize = 2 * iSizeH + 1;

  allocate(iSize, iSize, iSize);
  setValue(0.0);

  // kernel width
  setWidth(fWidth);

  // write Gaussian within 3 sigma-1D to map
  Real64 fBValue = -1.0 / (2.0 * fSigmaMap * fSigmaMap);
  Real64 fCValue = 9.0 * fSigmaMap * fSigmaMap;
  Real64 fScale = 0;
  Real64 fDSqu;
  Real64 fLapl = -1.0 / (PI * fSigmaMap * fSigmaMap * fSigmaMap * fSigmaMap);
  Real64 fTerm;

  unsigned int iX, iY, iZ;
  for (iZ = 0; iZ < iSize; iZ++)
    for (iY = 0; iY < iSize; iY++)
      for (iX = 0; iX < iSize; iX++) {

        fDSqu = (iX - iSizeH) * (iX - iSizeH) +
                (iY - iSizeH) * (iY - iSizeH) +
                (iZ - iSizeH) * (iZ - iSizeH);

        fTerm = exp(fDSqu * fBValue);

        if (fDSqu <= fCValue)
          setValue(iX, iY, iZ, fLapl * ((1.0 - fTerm)*fTerm));

        fScale += getValue(iX, iY, iZ);
      }

  for (iZ = 0; iZ < iSize; iZ++)
    for (iY = 0; iY < iSize; iY++)
      for (iX = 0; iX < iSize; iX++)
        setValue(iX, iY, iZ, getValue(iX, iY, iZ) / fScale);
};

/* /\***/
/*  * Filtrate Volume with Bilateral Filter */
/*  * Attention: This will overwrite the current content of the volume object with the filter kernel. It will allocate the memory internally. */
/*  * \param fSigma1D1 sigma of Gaussuian part of filter */
/*  * \param fsigma1D2 sigma of Intensity part of filter */
/*  * \param iSize size of the kernel */
/*  * \param bProgress show a progress bar */
/*  *\/ */
/* template<class T> */
/* inline void svt_volume<T>::filtrBilateral(double fSigma1D1, double fSigma1D2, unsigned int iSize,  bool bProgress) */
/* { */
/*     //unsigned int iSize = rKernel.getSizeX(); */
/*     //unsigned int iSizeS = iSize * iSize; */
/*     unsigned int iX, iY, iZ; */
/*     int iKX, iKY, iKZ; */
/*     Real64 fVal, fValInc, fValVol,fDeltaInc, fDeltaVol,  fDeltaIncSq, fDeltaVolSq, fKernel ; */

/*     Real64 fOldWidth = m_fWidth; */
/*     Real64 fOldGridX = m_fGridX; */
/*     Real64 fOldGridY = m_fGridY; */
/*     Real64 fOldGridZ = m_fGridZ; */

/*     T* pData = m_pData; */
/*     bool bInternal = m_bInternal; */
/*     m_pData = NULL; */
/*     m_bInternal = false; */

/*     svt_volume<T> oTmp; */
/*     oTmp.allocate( m_iSizeX, m_iSizeY, m_iSizeZ, 0.0 ); */

/*     m_pData = pData; */
/*     bInternal = bInternal; */
/*     m_bChanged = true; */

/*     Real64 fScale = 0; */

/*     // calculate dimension of kernel volume */

/*     // truncate at fSigmaFactor * fSigma1D */
/*     //unsigned int iSizeH1 = (int) ceil( fSigmaFactor * fSigma1D1 ); */
/*     //unsigned int iSizeH2 = (int) ceil( fSigmaFactor * fSigma1D2 ); */
/*     //max of iSizeH1 iSizeH2 */
/*     //unsigned int iSizeH = max(iSizeH1,iSizeH2); */


/*     iSize = 2 * iSize - 1; */

/*     int iDim = (int)((Real32)(iSize) * 0.5f); */
/*     int iStart = -iDim; */
/*     int iEnd = iDim+1; */

/*     try */
/*     { */

/*  for(iZ=0; iZ<m_iSizeZ; iZ++) */
/*  { */
/*      for(iX=0; iX<m_iSizeX; iX++) */
/*      { */
/*    for(iY=0; iY<m_iSizeY; iY++) */
/*    { */
/*        fVal = m_pData[(iX)+(iY*m_iSizeX)+(iZ*m_iSizeX*m_iSizeY)]; */

/*                     if (fVal != (T)(0)) */
/*                     { */
/*                         for(iKZ = iStart; iKZ < iEnd; iKZ++) */
/*                         for(iKY = iStart; iKY < iEnd; iKY++) */
/*                         for(iKX = iStart; iKX < iEnd; iKX++) */
/*                         { */
/*                            if ((((int)(iX)+iKX) < (int)(m_iSizeX) && ((int)(iX)+iKX) >= 0) && (((int)(iY)+iKY) < (int)(m_iSizeY) && ((int)(iY)+iKY) >= 0) && (((int)(iZ)+iKZ) < (int)(m_iSizeZ) && ((int)(iZ)+iKZ) >= 0)) */
/*                               { */
/*                                      fValInc   =   m_pData[(iX+iKX)+((iY+iKY)*m_iSizeX)+((iZ+iKZ)*m_iSizeX*m_iSizeY)]; */
/*                                      fDeltaInc = sqrt(iKX*iKX+iKY*iKY+iKZ*iKZ); */
/*                                      fDeltaIncSq = fDeltaInc * fDeltaInc; */

/*                                      fValVol   =   m_pData[(iX)+((iY)*m_iSizeX)+((iZ)*m_iSizeX*m_iSizeY)]; */
/*                                      fDeltaVol = sqrt((fValVol-fValInc)*(fValVol-fValInc)); */
/*                                      fDeltaVolSq = fDeltaVol * fDeltaVol; */

/*                                      fKernel= exp(-1/(2.0*fSigma1D1*fSigma1D1) * fDeltaIncSq) *exp(-1/(2.0*fSigma1D2*fSigma1D2)   * fDeltaVolSq ); */
/*                                      oTmp.setAt( iX, iY, iZ, oTmp.at( iX, iY, iZ ) + fKernel*fValInc); */
/*                                         fScale += fKernel*fValInc; */

/*                                } */

/*                        } */
/*                     } */

/*                 } */
/*             } */
/*  } */

/*     for (iZ=0; iZ<iSize; iZ++) */
/*  for (iY=0; iY<iSize; iY++) */
/*      for (iX=0; iX<iSize; iX++) */
/*        oTmp.setAt( iX, iY, iZ,  oTmp.at( iX, iY, iZ ) /fScale); */

/*  *this = oTmp; */

/*  m_bChanged = true; */

/*  // store map properties of old map in oTmp */
/*  m_fWidth = fOldWidth; */
/*  m_fGridX = fOldGridX; */
/*  m_fGridY = fOldGridY; */
/*  m_fGridZ = fOldGridZ; */

/*     } catch (int e) */
/*     { */
/*     } */

/* }; */




/**
 * Create an Identity kernel volume
 */
template<class T>
void svt_volume<T>::createIdentity()
{
  allocate(3, 3, 3, 0.0);
  setAt(1, 1, 1, 1.0);
};

/**
 * Create sphere at the center of the volume
 * \param fRadius radius of the sphere
 */
template<class T>
void svt_volume<T>::createSphere(Real64 fRadius)
{
  svt_vector4<Real64> oCenter;
  svt_vector4<Real64> oVoxel;
  Real64 fDist;
  oCenter.x((Real64)(m_iSizeX) * 0.5 * m_fWidth);
  oCenter.y((Real64)(m_iSizeY) * 0.5 * m_fWidth);
  oCenter.z((Real64)(m_iSizeZ) * 0.5 * m_fWidth);

  unsigned int iX, iY, iZ;

  for (iZ = 0; iZ < m_iSizeZ; iZ++)
    for (iY = 0; iY < m_iSizeY; iY++)
      for (iX = 0; iX < m_iSizeX; iX++) {
        oVoxel.x((Real64)(iX) * m_fWidth);
        oVoxel.y((Real64)(iY) * m_fWidth);
        oVoxel.z((Real64)(iZ) * m_fWidth);

        fDist = oVoxel.distance(oCenter);

        if (fDist < fRadius)
          setValue(iX, iY, iZ, 1.0 - (fDist / fRadius));
        else
          setValue(iX, iY, iZ, 0.0);
      }

};

///////////////////////////////////////////////////////////////////////////////
// Sampling functions
///////////////////////////////////////////////////////////////////////////////

/**
 * Set cutoff for sampling
 * \param fCutoff a voxel value lower than this value is not considered for the sampling
 */
template<class T>
void svt_volume<T>::setCutoff(Real64 fCutoff)
{
  m_fCutoff = fCutoff;
};
/**
 * Get cutoff for sampling
 * \return a voxel value lower than this value is not considered for the sampling
 */
template<class T>
Real64 svt_volume<T>::getCutoff() const
{
  return m_fCutoff;
};

/**
 * sample the object randomly and return a vector that refrects the probability distribution of the object
 */
template<class T>
svt_vector4<Real64> svt_volume<T>::sample()
{
  T fMax = getMaxDensity();
  Real64 fNumV = (Real64)(m_iSizeX * m_iSizeY * m_iSizeZ);
  bool bFound = false;
  unsigned int iIndex;
  T fVoxel, fDensity;

  while (bFound == false) {
    iIndex = (unsigned int)(svt_genrand() * fNumV);
    fVoxel = at(iIndex);
    fDensity = (T)(svt_genrand() * fMax);

    if (fVoxel > fDensity && fVoxel > m_fCutoff)
      bFound = true;
  }

  svt_vector4<Real64> oVec;

  unsigned int iSizeXY = m_iSizeX * m_iSizeY;
  unsigned int iIndz = iIndex / iSizeXY;
  iIndex -= iIndz * iSizeXY;
  unsigned int iIndy = iIndex / m_iSizeX;
  unsigned int iIndx = iIndex - iIndy * m_iSizeX;

  oVec.x(getGridX() + (iIndx * m_fWidth));
  oVec.y(getGridY() + (iIndy * m_fWidth));
  oVec.z(getGridZ() + (iIndz * m_fWidth));

  //oVec.x( (iIndx * m_fWidth) - ( m_iSizeX * m_fWidth * 0.5 ) );
  //oVec.y( (iIndy * m_fWidth) - ( m_iSizeY * m_fWidth * 0.5 ) );
  //oVec.z( (iIndz * m_fWidth) - ( m_iSizeZ * m_fWidth * 0.5 ) );

  return oVec;
};

/**
 * Sample inside sphere the object randomly and return a vector that refrects the probability distribution of the object
 * \param oCenter object class svt_vector4 with xyz coordinates of center of sphere
 * \param fRadius radius of the sphere
 */
template<class T>
svt_vector4<Real64> svt_volume<T>::sampleSphere(svt_vector4<Real64> oCenter, Real64 fRadius)
{
  T fMax = getMaxDensity();
  Real64 fNumV = (Real64)(m_iSizeX * m_iSizeY * m_iSizeZ);
  Real64 fDist, fGridX, fGridY, fGridZ;

  bool bFound = false;
  unsigned int iIndex;
  T fVoxel, fDensity;

  unsigned int iSizeXY = m_iSizeX * m_iSizeY;
  unsigned int iIndx, iIndy, iIndz;

  fGridX = getGridX();
  fGridY = getGridY();
  fGridZ = getGridZ();


  while (bFound == false) {
    iIndex = (unsigned int)(svt_genrand() * fNumV);
    fVoxel = at(iIndex);
    fDensity = (T)(svt_genrand() * fMax);


    iIndz = (iIndex / iSizeXY);
    iIndex -= iIndz * iSizeXY;
    iIndy = (iIndex / m_iSizeX);
    iIndx = (iIndex - iIndy * m_iSizeX);


    fDist = ((fGridX + (iIndx * m_fWidth)) - oCenter.x()) * ((fGridX + (iIndx * m_fWidth)) - oCenter.x()) +
            ((fGridY + (iIndy * m_fWidth)) - oCenter.y()) * ((fGridY + (iIndy * m_fWidth)) - oCenter.y()) +
            ((fGridZ + (iIndz * m_fWidth)) - oCenter.z()) * ((fGridZ + (iIndz * m_fWidth)) - oCenter.z());
    fDist = sqrt(fDist);

    if (fVoxel > fDensity && fVoxel > m_fCutoff && fDist < fRadius)
      bFound = true;
  }

  svt_vector4<Real64> oVec;

  oVec.x(fGridX + (iIndx * m_fWidth));
  oVec.y(fGridY + (iIndy * m_fWidth));
  oVec.z(fGridZ + (iIndz * m_fWidth));

  return oVec;
};

/**
* gives the best ISO threshold for fitting the object with a given volume
* this outputs the corresponding threshold for the data set
* this algorithm adjusts for the voxel spacing of each object
*/
template<class T>
T svt_volume<T>::bestISO(svt_volume<T> &rTarget, T fThreshold)
{
  double fTargetVol = 0; //Holds target volume
  vector<T> oVoxels; //Holds the voxel data

  for (unsigned int i = 0; i < rTarget.getSizeX(); i++) //Get target voxel amount
    for (unsigned int j = 0; j < rTarget.getSizeY(); j++)
      for (unsigned int k = 0; k < rTarget.getSizeZ(); k++)
        if (rTarget.at(i, j, k) >= fThreshold)
          fTargetVol++;

  T fTargetCubed = rTarget.getWidth() * rTarget.getWidth() * rTarget.getWidth();
  T fCubed = m_fWidth * m_fWidth * m_fWidth;
  fTargetVol *= fTargetCubed / fCubed; //Adjust for voxel spacing to get volume

  for (unsigned int i = 0; i < getSizeX(); i++) //Organize voxel data
    for (unsigned int j = 0; j < getSizeY(); j++)
      for (unsigned int k = 0; k < getSizeZ(); k++)
        if (getValue(i, j, k) > 0) oVoxels.push_back(getValue(i, j, k));
  sort(oVoxels.begin(), oVoxels.end()); //Sort the voxel information

  int fIndex = oVoxels.size() - fTargetVol;
  if (fIndex < 0) fIndex = 0;
  if (fIndex > oVoxels.size() - 1) fIndex = oVoxels.size() - 1;
  return oVoxels.at(fIndex); //Return the ISO threshold
};

/**
 * Calculates the corresponding isosurface threshold value for this volume. The algorithm tries to cover exactly the voxel values of the other volume with the ones from here, so that
 * both isosurfaces have a similar size. rVol_A would typically have a higher resolution but show the exact same system.
 *
 * \param rVol_A other volume
 * \param fThresh good known threshold value of rVol_A
 */
template<class T>
T svt_volume<T>::correspondingISO(svt_volume<T> &rVol_A, T fThresh_A)
{
  Real64 fX, fY, fZ;
  T fThresh_B = (T)(1.0E10);

  for (unsigned int iX = 0; iX < this->getSizeX(); iX++)
    for (unsigned int iY = 0; iY < this->getSizeY(); iY++)
      for (unsigned int iZ = 0; iZ < this->getSizeZ(); iZ++) {
        fX = ((Real64)(iX) * m_fWidth) + m_fGridX;
        fY = ((Real64)(iY) * m_fWidth) + m_fGridY;
        fZ = ((Real64)(iZ) * m_fWidth) + m_fGridZ;

        if (rVol_A.getRealSpaceValue(fX, fY, fZ) >= fThresh_A && this->getValue(iX, iY, iZ) < fThresh_B)
          fThresh_B = this->getValue(iX, iY, iZ);
      }

  return fThresh_B;
};

/**
* scales the map by a given factor
*/
template<class T>
void svt_volume<T>::scale(T fScale)
{
  for (unsigned int i = 0; i < getSizeX(); i++)
    for (unsigned int j = 0; j < getSizeY(); j++)
      for (unsigned int k = 0; k < getSizeZ(); k++)
        setAt(i, j, k, at(i, j, k) * fScale);
  return;
};

/**
* shifts the map by a given offset
*/
template<class T>
void svt_volume<T>::shift(T fShift)
{
  for (unsigned int i = 0; i < getSizeX(); i++)
    for (unsigned int j = 0; j < getSizeY(); j++)
      for (unsigned int k = 0; k < getSizeZ(); k++)
        setAt(i, j, k, at(i, j, k) + fShift);
  return;
};

/**
* scales using a given slope (for testing)
*/
template<class T>
void svt_volume<T>::scaleBySlope(T fSlope)
{
  for (unsigned int i = 0; i < getSizeX(); i++)
    for (unsigned int j = 0; j < getSizeY(); j++)
      for (unsigned int k = 0; k < getSizeZ(); k++)
        setAt(i, j, k, at(i, j, k) * at(i, j, k) * fSlope);
  return;
};

/**
* returns the slope given two points
*/
template<class T>
T svt_volume<T>::getSlope(T fX1, T fY1, T fX2, T fY2)
{
  return ((fY2 - fY1) / (fX2 - fX1));
};

/**
* interpolated scaling of the volume based on given information
* takes in two old thresholds and two new thresholds (from bestISO method)
*/
template<class T>
void svt_volume<T>::interpolatedScale(T fThreshold1, T fThreshold2, T fNew1, T fNew2)
{
  T fScale1 = fNew1 / fThreshold1;
  T fScale2 = fNew2 / fThreshold2;
  T fSlope = getSlope(fThreshold1, fScale1, fThreshold2, fScale2);
  T fShift = fNew1 - fSlope * fThreshold1;
  for (unsigned int i = 0; i < getSizeX(); i++)
    for (unsigned int j = 0; j < getSizeY(); j++)
      for (unsigned int k = 0; k < getSizeZ(); k++)
        if (at(i, j, k) > 0) {
          setAt(i, j, k, (sqrt(at(i, j, k)) / sqrt(fSlope)));
        }
  shift(fShift); //Adjust for shift value
  return;
};

/**
* calculates the root mean square of the volume
*/
template<class T>
T svt_volume<T>::getRMS()
{
  T fSum = 0;
  for (unsigned int i = 0; i < getSizeX(); i++)
    for (unsigned int j = 0; j < getSizeY(); j++)
      for (unsigned int k = 0; k < getSizeZ(); k++)
        fSum += pow(at(i, j, k), 2);
  fSum /= (getSizeX() * getSizeY() * getSizeZ());
  return sqrt(fSum);
};

/**
 * Remove every structure that has more than two neighbors
 * \param fThreshold voxel above that value are considered occupied
 */
template<class T>
void svt_volume<T>::removeNeighbors(T fThreshold)
{
  for (unsigned int iX = 0; iX < (unsigned int)getSizeX(); iX++)
    for (unsigned int iY = 0; iY < (unsigned int)getSizeY(); iY++)
      for (unsigned int iZ = 0; iZ < (unsigned int)getSizeZ(); iZ++)
        if (getValue(iX, iY, iZ) > fThreshold  && numNeighbors(iX, iY, iZ, fThreshold) == 1)
          setValue(iX, iY, iZ, 1);
        else
          setValue(iX, iY, iZ, 0);
};

/**
 * How many direct neighbors does a voxel have?
 * \return number of occupied neighboring voxels
 */
template<class T>
unsigned int svt_volume<T>::numNeighbors(unsigned int iX, unsigned int iY, unsigned int iZ, T fThreshold)
{
  int iVX = (int)iX;
  int iVY = (int)iY;
  int iVZ = (int)iZ;

  unsigned int iNeighb = 0;

  for (int iKX = -1; iKX < 2; iKX++)
    for (int iKY = -1; iKY < 2; iKY++)
      for (int iKZ = -1; iKZ < 2; iKZ++) {
        if (
          (iVX + iKX >= 0 && iVX + iKX < getSizeX()) &&
          (iVY + iKY >= 0 && iVY + iKY < getSizeY()) &&
          (iVZ + iKZ >= 0 && iVZ + iKZ < getSizeZ()) &&
          (getValue(iVX + iKX, iVY + iKY, iVZ + iKZ) > fThreshold)
        )
          iNeighb++;
      }

  return iNeighb;
};

///////////////////////////////////////////////////////////////////////////////
// Internal functions
///////////////////////////////////////////////////////////////////////////////

/**
 * Calculate/update the minimum and maximum density values.
 */
template<class T>
inline void svt_volume<T>::calcMinMaxDensity()
{
  unsigned int iCount;
  unsigned int iNVox = m_iSizeX * m_iSizeY * m_iSizeZ;

  m_fMaxDensity = getValue(0, 0, 0);
  m_fMinDensity = getValue(0, 0, 0);

  double fAvgDensity = (double)(0.0);
  m_fNorm = 0.0;
  T fValue;

  for (iCount = 0; iCount < iNVox; iCount++) {
    fValue = getValue(iCount);

    if (m_fMaxDensity < fValue)
      m_fMaxDensity = fValue;

    if (m_fMinDensity > fValue)
      m_fMinDensity = fValue;

    m_fNorm += fValue * fValue;
    fAvgDensity += (double)(getValue(iCount));
  }

  m_fNorm = sqrt(m_fNorm);
  fAvgDensity /= (double)(iNVox);

  m_fAvgDensity = (T)(fAvgDensity);

  m_bChanged = false;
};

/**
 * Calculate correlation with other svt_volume object
 * \param rVol reference to other svt_volume object
 */
template<class T>
Real64 svt_volume<T>::correlation(svt_volume<T> &rVol, bool bMask)
{
  // is there actually data in the two objects?
  if (size() == 0 || rVol.size() == 0)
    return -1.0e9;

  // is the voxelwidth indeed identical?
  if (fabs(m_fWidth - rVol.getWidth()) > EPSILON) {
    SVTLBBO << "Voxel width different: " << m_fWidth << " != " << rVol.getWidth() << endl;
    error_sba(85010, "The cross-correlation coefficient could not be computed due to different voxel width!");
    return -1.0e9;
  }

  // calculate correlation
  int iX, iY, iZ;
  Real64 fCorrelation = 0.0;
  Real64 fOrgX = this->getGridX() - rVol.getGridX();
  Real64 fOrgY = this->getGridY() - rVol.getGridY();
  Real64 fOrgZ = this->getGridZ() - rVol.getGridZ();
  // use of (int)floor to make the rounding correct for negative numbers
  int iOrgX = (int)floor(fOrgX / m_fWidth + 0.5f);
  int iOrgY = (int)floor(fOrgY / m_fWidth + 0.5f);
  int iOrgZ = (int)floor(fOrgZ / m_fWidth + 0.5f);
  Real64 fValueA, fValueB;
  Real64 fValueA_Sum = 0,  fValueB_Sum = 0, fValueA_Sq_Sum = 0, fValueB_Sq_Sum = 0, fCoSum = 0;
  unsigned int iAddedVoxels = 0;

  for (iZ = 0; iZ < (int)(m_iSizeZ); iZ++)
    for (iY = 0; iY < (int)(m_iSizeY); iY++)
      for (iX = 0; iX < (int)(m_iSizeX); iX++)
        if (iX + iOrgX >= 0 && iY + iOrgY >= 0 && iZ + iOrgZ >= 0 && iX + iOrgX < (int)(rVol.getSizeX()) && iY + iOrgY < (int)(rVol.getSizeY()) && iZ + iOrgZ < (int)(rVol.getSizeZ())) {
          fValueA = (Real64)(this->at(iX, iY, iZ));
          fValueB = (Real64)(rVol.at(iX + iOrgX, iY + iOrgY, iZ + iOrgZ));

          if (!bMask || (bMask && fValueA != 0)) {
            fCoSum    += fValueA * fValueB;
            fValueA_Sq_Sum  += fValueA * fValueA;
            fValueB_Sq_Sum  += fValueB * fValueB;
            fValueA_Sum   += fValueA;
            fValueB_Sum   += fValueB;
            iAddedVoxels ++;
          }
        }
  if (sqrt((iAddedVoxels * fValueA_Sq_Sum - fValueA_Sum * fValueA_Sum) * (iAddedVoxels * fValueB_Sq_Sum - fValueB_Sum * fValueB_Sum)) != 0)
    fCorrelation = (iAddedVoxels * fCoSum - fValueA_Sum * fValueB_Sum) / sqrt((iAddedVoxels * fValueA_Sq_Sum - fValueA_Sum * fValueA_Sum) * (iAddedVoxels * fValueB_Sq_Sum - fValueB_Sum * fValueB_Sum));
  else
    fCorrelation = 0.0;

  return fCorrelation;
};

/**
 * Calculate correlation with other svt_volume object
 * \param rVol reference to other svt_volume object
 */
template<class T>
Real64 svt_volume<T>::correlationColores(svt_volume<T> &rVol)
{
  // calculate norm
  Real64 fNormA = this->getNorm();
  Real64 fNormB = rVol.getNorm();

  if (fNormA == 0.0 || fNormB == 0.0)
    return 0.0;

  // calculate correlation
  int iX, iY, iZ;
  Real64 fCorrelation = 0.0;
  Real64 fOrgX = this->getGridX() - rVol.getGridX();
  Real64 fOrgY = this->getGridY() - rVol.getGridY();
  Real64 fOrgZ = this->getGridZ() - rVol.getGridZ();
  // use of (int)floor to make the rounding correct for negative numbers
  int iOrgX = (int)floor(fOrgX / m_fWidth + 0.5f);
  int iOrgY = (int)floor(fOrgY / m_fWidth + 0.5f);
  int iOrgZ = (int)floor(fOrgZ / m_fWidth + 0.5f);

  for (iZ = 0; iZ < (int)(m_iSizeZ); iZ++)
    for (iY = 0; iY < (int)(m_iSizeY); iY++)
      for (iX = 0; iX < (int)(m_iSizeX); iX++)
        if (iX + iOrgX >= 0 && iY + iOrgY >= 0 && iZ + iOrgZ >= 0 && iX + iOrgX < (int)(rVol.getSizeX()) && iY + iOrgY < (int)(rVol.getSizeY()) && iZ + iOrgZ < (int)(rVol.getSizeZ()))
          fCorrelation += this->at(iX, iY, iZ) * rVol.at(iX + iOrgX, iY + iOrgY, iZ + iOrgZ);

  fCorrelation /= fNormA * fNormB;

  return fCorrelation;
};

/**
 * Calculates the internal correlation of the voxel values.
 */
template<class T>
svt_volume<Real64> *svt_volume<T>::internalCorr(unsigned int iWidth)
{
  Real64 fCorr;
  svt_volume<Real64> *pCorrVol = new svt_volume<Real64>(this->getSizeX(), this->getSizeY(), this->getSizeZ());

  try {

    for (int iX = 0; iX < (int)this->getSizeX(); iX++) {

      for (int iY = 0; iY < (int)this->getSizeY(); iY++)
        for (int iZ = 0; iZ < (int)this->getSizeZ(); iZ++) {
          fCorr = 0.0;

          for (int iKX = -iWidth; iKX < (int)iWidth; iKX++)
            for (int iKY = -iWidth; iKY < (int)iWidth; iKY++)
              for (int iKZ = -iWidth; iKZ < (int)iWidth; iKZ++)

                if (
                  (iKX != 0 || iKY != 0 || iKZ != 0) &&
                  (iX + iKX >= 0 && iX + iKX < (int)this->getSizeX()) &&
                  (iY + iKY >= 0 && iY + iKY < (int)this->getSizeY()) &&
                  (iZ + iKZ >= 0 && iZ + iKZ < (int)this->getSizeZ())
                )

                  fCorr += getValue(iX, iY, iZ) * getValue(iX + iKX, iY + iKY, iZ + iKZ);

          pCorrVol->setValue(iX, iY, iZ, fCorr);
        }
    }

  } catch (int e) {
    SVTLBBO << "Warning: Internal correlation calculation was aborted, volume data might be incomplete..." << endl;
  }

  return pCorrVol;
};


/**
 * Grow - this function blows up the volume at a different size ( only padding routine!).
 * \param iNewSizeX new size (x dimension)
 * \param iNewSizeY new size (y dimension)
 * \param iNewSizeZ new size (z dimension)
 */
template<class T>
void svt_volume<T>::resize(unsigned int iNewSizeX, unsigned int iNewSizeY, unsigned int iNewSizeZ)
{


  unsigned int iOldSizeX = getSizeX();
  unsigned int iOldSizeY = getSizeY();
  unsigned int iOldSizeZ = getSizeZ();

  if ((iOldSizeX <= iNewSizeX) && (iOldSizeY <= iNewSizeY) && (iOldSizeZ <= iNewSizeZ)) {

    svt_volume<T> *pVol = new svt_volume<T>(iNewSizeX, iNewSizeY, iNewSizeZ, 0.0);

    Real64 fVoxelValue;
    Real64 fWidth = getWidth();

    Real64 fOldGridX = getGridX();
    Real64 fOldGridY = getGridY();
    Real64 fOldGridZ = getGridZ();

    Real64 fNewGridX = (fOldGridX / fWidth - abs(double(iOldSizeX) - double(iNewSizeX)) / (2.0)) * fWidth;
    Real64 fNewGridY = (fOldGridY / fWidth - abs(double(iOldSizeY) - double(iNewSizeY)) / (2.0)) * fWidth;
    Real64 fNewGridZ = (fOldGridZ / fWidth - abs(double(iOldSizeZ) - double(iNewSizeZ)) / (2.0)) * fWidth;

    unsigned int iCurrX, iCurrY, iCurrZ;


    for (unsigned int iX = 0; iX < iNewSizeX; iX++)
      for (unsigned int iY = 0; iY < iNewSizeY; iY++)
        for (unsigned int iZ = 0; iZ < iNewSizeZ; iZ++)
          if ((iX >= floor(iNewSizeX / 2 - iOldSizeX / 2)) && (iX < floor(iNewSizeX / 2 + iOldSizeX / 2)) &&
              (iY >= floor(iNewSizeY / 2 - iOldSizeY / 2)) && (iY < floor(iNewSizeY / 2 + iOldSizeY / 2)) &&
              (iZ >= floor(iNewSizeZ / 2 - iOldSizeZ / 2)) && (iZ < floor(iNewSizeZ / 2 + iOldSizeZ / 2))) {
            iCurrX = iX - floor(iNewSizeX / 2 - iOldSizeX / 2);
            iCurrY = iY - floor(iNewSizeY / 2 - iOldSizeY / 2);
            iCurrZ = iZ - floor(iNewSizeZ / 2 - iOldSizeZ / 2);

            fVoxelValue = getValue(iCurrX, iCurrY, iCurrZ);
            pVol->setValue(iX, iY, iZ, fVoxelValue);
          }



    pVol -> setWidth(fWidth);
    pVol -> setGrid(fNewGridX, fNewGridY, fNewGridZ);

    *this = (*pVol);
    m_bChanged = true;

  } else {
    SVTLBBO << "Cannot resize to smaller volume size" << endl;
  }

};

/**
 * Calculate number of occupied voxels
 * \param fThreshold threshold value - only voxels higher than the threshold are counted.
 * \return number of voxels with a density higher than fThreshold
 */
template<class T>
unsigned long svt_volume<T>::getOccupied(T fThreshold) const
{
  unsigned long iCount = 0;
  unsigned long iNVox = m_iSizeX * m_iSizeY * m_iSizeZ;
  unsigned long i;

  for (i = 0; i < iNVox; i++)
    if (at(i) >= fThreshold)
      iCount++;

  return iCount;
}

///////////////////////////////////////////////////////////////////////////////
// File formats
///////////////////////////////////////////////////////////////////////////////

/**
 * Load a file. This function looks at the extension to determine which function has to be used to actually load the file.
 * \param pFname pointer to array of char with the filename
 */
template<class T>
inline svt_matrix4<T> svt_volume<T>::load(const char *pFname)
{
  T *pPhi = NULL;
  read_vol((char *)pFname, &m_fWidth, &m_fGridX, &m_fGridY, &m_fGridZ, &m_iSizeX, &m_iSizeY, &m_iSizeZ, &pPhi);
  setData(m_iSizeX, m_iSizeY, m_iSizeZ, pPhi);

  m_bInternal = true;
  svt_matrix4<T> oMatrix;
  return oMatrix;
};

/**
 * Save a file. This function looks at the extension to determine which function has to be used to actually save the file.
 * \param pFname pointer to array of char with the filename
 */
template<class T>
inline void svt_volume<T>::save(const char *pFname)
{
  write_vol((char *)pFname, m_fWidth, m_fGridX, m_fGridY, m_fGridZ, m_iSizeX, m_iSizeY, m_iSizeZ, m_pData);
};



/**
 * Non-recursive flood-fill segmentation algorithm. All voxels that are connected to a start voxel and are above the threshold are kept, the others are removed.
 * The algorithm creates a mask that is later blurred by a Gaussian kernel. The sigma of the gaussian can be specified.
 * \param iStartX x index of the start voxel
 * \param iStartY y index of the start voxel
 * \param iStartZ z index of the start voxel
 * \param fTreshold threshold for the floodfill
 * \param fGaussian sigma of the gaussian the mask is convoluted with (if 0, no blurring of the mask is applied)
 */
template<class T>
void svt_volume<T>::floodfill_segmentation(unsigned int iStartX, unsigned int iStartY, unsigned int iStartZ, T fThreshold, Real64 fGaussian)
{
  svt_volume<T> oVisited(this->m_iSizeX, this->m_iSizeY, this->m_iSizeZ);
  oVisited = (T)(0);

  stack< svt_vector4<int> > oQueue;
  stack< svt_vector4<int> > oQueueNotIn;
  svt_vector4<int> oStart(iStartX, iStartY, iStartZ);
  oQueue.push(oStart);

  if (this->getValue(iStartX, iStartY, iStartZ) < fThreshold)
    SVTLBBO << "Warning: Start voxel is lower than threshold in floodfill!" << endl;

  svt_vector4<int> oCurrent;
  svt_vector4<int> oTemp;

  SVTLBBO << "Floodfill: Generating mask..." << endl;

  unsigned int iRegion = 0;
  unsigned int iNVox = m_iSizeX * m_iSizeY * m_iSizeZ;
  unsigned int iTotalRegion = iNVox;
  unsigned int iRadius = 10;

  try {
    while ((oQueue.size() > 0 || oQueueNotIn.size() > 0)) {
      // select new region to be manipulated, if no voxels are left in the current region
      if (iRegion == 0 || oQueue.size() == 0) {
        for (unsigned int i = 0; i < oQueueNotIn.size(); i++) {
          oQueue.push(oQueueNotIn.top());
          oQueueNotIn.pop();
        }

        unsigned int iSX = 0;
        unsigned int iEX = m_iSizeX;
        unsigned int iSY = 0;
        unsigned int iEY = m_iSizeY;
        unsigned int iSZ = 0;
        unsigned int iEZ = m_iSizeZ;

        if (iStartX > iRadius)
          iSX = iStartX - iRadius;
        if (iStartY > iRadius)
          iSY = iStartY - iRadius;
        if (iStartZ > iRadius)
          iSZ = iStartZ - iRadius;
        if (iStartX + iRadius < m_iSizeX)
          iEX = iStartX + iRadius;
        if (iStartY + iRadius < m_iSizeY)
          iEY = iStartY + iRadius;
        if (iStartZ + iRadius < m_iSizeZ)
          iEZ = iStartZ + iRadius;

        for (unsigned int iX = iSX; iX < iEX; iX++)
          for (unsigned int iY = iSY; iY < iEY; iY++)
            for (unsigned int iZ = iSZ; iZ < iEZ; iZ++)
              if (sqrt((iX - iStartX) * (iX - iStartX) + (iY - iStartY) * (iY - iStartY) + (iZ - iStartZ) * (iZ - iStartZ)) < iRadius)
                if (oVisited.getValue(iX, iY, iZ) == 0) {
                  oVisited.setValue(iX, iY, iZ, 3);
                  iRegion++;
                  iTotalRegion--;
                }
        iRadius += 10;

      }

      oCurrent = oQueue.top();
      oQueue.pop();

      if (oVisited.getValue(oCurrent.x(), oCurrent.y(), oCurrent.z()) == 0) {
        oQueueNotIn.push(oCurrent);

      } else if (oVisited.getValue(oCurrent.x(), oCurrent.y(), oCurrent.z()) == 3) {
        iRegion--;

        if (this->getValue(oCurrent.x(), oCurrent.y(), oCurrent.z()) > fThreshold) {
          oVisited.setValue(oCurrent.x(), oCurrent.y(), oCurrent.z(), 1);

          // z = 0

          // x-1 y+1 z+0
          if (oCurrent.x() - 1 > 0 && oCurrent.y() + 1 < (int)this->m_iSizeY && (oVisited.getValue(oCurrent.x() - 1, oCurrent.y() + 1, oCurrent.z()) == 0 || oVisited.getValue(oCurrent.x() - 1, oCurrent.y() + 1, oCurrent.z()) == 3)) {
            oTemp = oCurrent;
            oTemp.x(oTemp.x() - 1);
            oTemp.y(oTemp.y() + 1);
            oQueue.push(oTemp);
          }
          // x+0 y+1 z+0
          if (oCurrent.y() + 1 < (int)this->m_iSizeY && (oVisited.getValue(oCurrent.x(), oCurrent.y() + 1, oCurrent.z()) == 0 || oVisited.getValue(oCurrent.x(), oCurrent.y() + 1, oCurrent.z()) == 3)) {
            oTemp = oCurrent;
            oTemp.y(oTemp.y() + 1);
            oQueue.push(oTemp);
          }
          // x+1 y+1 z+0
          if (oCurrent.x() + 1 < (int)this->m_iSizeX && oCurrent.y() + 1 < (int)this->m_iSizeY && (oVisited.getValue(oCurrent.x() + 1, oCurrent.y() + 1, oCurrent.z()) == 0 || oVisited.getValue(oCurrent.x() + 1, oCurrent.y() + 1, oCurrent.z()) == 3)) {
            oTemp = oCurrent;
            oTemp.x(oTemp.x() + 1);
            oTemp.y(oTemp.y() + 1);
            oQueue.push(oTemp);
          }

          // x-1 y+0 z+0
          if (oCurrent.x() - 1 > 0 && (oVisited.getValue(oCurrent.x() - 1, oCurrent.y(), oCurrent.z()) == 0 || oVisited.getValue(oCurrent.x() - 1, oCurrent.y(), oCurrent.z()) == 3)) {
            oTemp = oCurrent;
            oTemp.x(oTemp.x() - 1);
            oQueue.push(oTemp);
          }
          // x+1 y+0 z+0
          if (oCurrent.x() + 1 < (int)this->m_iSizeX && (oVisited.getValue(oCurrent.x() + 1, oCurrent.y(), oCurrent.z()) == 0 || oVisited.getValue(oCurrent.x() + 1, oCurrent.y(), oCurrent.z()) == 3)) {
            oTemp = oCurrent;
            oTemp.x(oTemp.x() + 1);
            oQueue.push(oTemp);
          }

          // x-1 y-1 z+0
          if (oCurrent.x() - 1 > 0 && oCurrent.y() - 1 > 0 && (oVisited.getValue(oCurrent.x() - 1, oCurrent.y(), oCurrent.z()) == 0 || oVisited.getValue(oCurrent.x() - 1, oCurrent.y(), oCurrent.z()) == 3)) {
            oTemp = oCurrent;
            oTemp.x(oTemp.x() - 1);
            oTemp.y(oTemp.y() - 1);
            oQueue.push(oTemp);
          }
          // x+0 y-1 z+0
          if (oCurrent.y() - 1 > 0 && (oVisited.getValue(oCurrent.x(), oCurrent.y() - 1, oCurrent.z()) == 0 || oVisited.getValue(oCurrent.x(), oCurrent.y() - 1, oCurrent.z()) == 3)) {
            oTemp = oCurrent;
            oTemp.y(oTemp.y() - 1);
            oQueue.push(oTemp);
          }
          // x+1 y-1 z+0
          if (oCurrent.x() + 1 < (int)this->m_iSizeX && oCurrent.y() - 1 > 0 && (oVisited.getValue(oCurrent.x() + 1, oCurrent.y(), oCurrent.z()) == 0 || oVisited.getValue(oCurrent.x() + 1, oCurrent.y(), oCurrent.z()) == 3)) {
            oTemp = oCurrent;
            oTemp.x(oTemp.x() + 1);
            oTemp.y(oTemp.y() - 1);
            oQueue.push(oTemp);
          }

          // z = 1
          if (oCurrent.z() + 1 < (int)this->m_iSizeZ) {

            // x-1 y+1 z+1
            if (oCurrent.x() - 1 > 0 && oCurrent.y() + 1 < (int)this->m_iSizeY && (oVisited.getValue(oCurrent.x() - 1, oCurrent.y() + 1, oCurrent.z() + 1) == 0 || oVisited.getValue(oCurrent.x() - 1, oCurrent.y() + 1, oCurrent.z() + 1) == 3)) {
              oTemp = oCurrent;
              oTemp.x(oTemp.x() - 1);
              oTemp.y(oTemp.y() + 1);
              oTemp.z(oTemp.z() + 1);
              oQueue.push(oTemp);
            }
            // x+0 y+1 z+1
            if (oCurrent.y() + 1 < (int)this->m_iSizeY && (oVisited.getValue(oCurrent.x(), oCurrent.y() + 1, oCurrent.z() + 1) == 0 || oVisited.getValue(oCurrent.x(), oCurrent.y() + 1, oCurrent.z() + 1) == 3)) {
              oTemp = oCurrent;
              oTemp.y(oTemp.y() + 1);
              oTemp.z(oTemp.z() + 1);
              oQueue.push(oTemp);
            }
            // x+1 y+1 z+1
            if (oCurrent.x() + 1 < (int)this->m_iSizeX && oCurrent.y() + 1 < (int)this->m_iSizeY && (oVisited.getValue(oCurrent.x() + 1, oCurrent.y() + 1, oCurrent.z() + 1) == 0 || oVisited.getValue(oCurrent.x() + 1, oCurrent.y() + 1, oCurrent.z() + 1) == 3)) {
              oTemp = oCurrent;
              oTemp.x(oTemp.x() + 1);
              oTemp.y(oTemp.y() + 1);
              oTemp.z(oTemp.z() + 1);
              oQueue.push(oTemp);
            }

            // x-1 y+0 z+1
            if (oCurrent.x() - 1 > 0 && (oVisited.getValue(oCurrent.x() - 1, oCurrent.y(), oCurrent.z() + 1) == 0 || oVisited.getValue(oCurrent.x() - 1, oCurrent.y(), oCurrent.z() + 1) == 3)) {
              oTemp = oCurrent;
              oTemp.x(oTemp.x() - 1);
              oTemp.z(oTemp.z() + 1);
              oQueue.push(oTemp);
            }
            // x+0 y+0 z+1
            if (oVisited.getValue(oCurrent.x(), oCurrent.y(), oCurrent.z() + 1) == 0 || oVisited.getValue(oCurrent.x(), oCurrent.y(), oCurrent.z() + 1) == 3) {
              oTemp = oCurrent;
              oTemp.z(oTemp.z() + 1);
              oQueue.push(oTemp);
            }
            // x+1 y+0 z+1
            if (oCurrent.x() + 1 < (int)this->m_iSizeX && (oVisited.getValue(oCurrent.x() + 1, oCurrent.y(), oCurrent.z() + 1) == 0 || oVisited.getValue(oCurrent.x() + 1, oCurrent.y(), oCurrent.z() + 1) == 3)) {
              oTemp = oCurrent;
              oTemp.x(oTemp.x() + 1);
              oTemp.z(oTemp.z() + 1);
              oQueue.push(oTemp);
            }

            // x-1 y-1 z+1
            if (oCurrent.x() - 1 > 0 && oCurrent.y() - 1 > 0 && (oVisited.getValue(oCurrent.x() - 1, oCurrent.y() - 1, oCurrent.z() + 1) == 0 || oVisited.getValue(oCurrent.x() - 1, oCurrent.y() - 1, oCurrent.z() + 1) == 3)) {
              oTemp = oCurrent;
              oTemp.x(oTemp.x() - 1);
              oTemp.y(oTemp.y() - 1);
              oTemp.z(oTemp.z() + 1);
              oQueue.push(oTemp);
            }
            // x+0 y-1 z+1
            if (oCurrent.y() - 1 > 0 && (oVisited.getValue(oCurrent.x(), oCurrent.y() - 1, oCurrent.z() + 1) == 0 || oVisited.getValue(oCurrent.x(), oCurrent.y() - 1, oCurrent.z() + 1) == 3)) {
              oTemp = oCurrent;
              oTemp.y(oTemp.y() - 1);
              oTemp.z(oTemp.z() + 1);
              oQueue.push(oTemp);
            }
            // x+1 y-1 z+1
            if (oCurrent.x() + 1 < (int)this->m_iSizeX && oCurrent.y() - 1 > 0 && (oVisited.getValue(oCurrent.x() + 1, oCurrent.y() - 1, oCurrent.z() + 1) == 0 || oVisited.getValue(oCurrent.x() + 1, oCurrent.y() - 1, oCurrent.z() + 1) == 3)) {
              oTemp = oCurrent;
              oTemp.x(oTemp.x() + 1);
              oTemp.y(oTemp.y() - 1);
              oTemp.z(oTemp.z() + 1);
              oQueue.push(oTemp);
            }

          }

          // z = -1
          if (oCurrent.z() - 1 > 0) {

            // x-1 y+1 z-1
            if (oCurrent.x() - 1 > 0 && oCurrent.y() + 1 < (int)this->m_iSizeY && (oVisited.getValue(oCurrent.x() - 1, oCurrent.y() + 1, oCurrent.z() - 1) == 0 || oVisited.getValue(oCurrent.x() - 1, oCurrent.y() + 1, oCurrent.z() - 1) == 3)) {
              oTemp = oCurrent;
              oTemp.x(oTemp.x() - 1);
              oTemp.y(oTemp.y() + 1);
              oTemp.z(oTemp.z() - 1);
              oQueue.push(oTemp);
            }
            // x+0 y+1 z-1
            if (oCurrent.y() + 1 < (int)this->m_iSizeY && (oVisited.getValue(oCurrent.x(), oCurrent.y() + 1, oCurrent.z() - 1) == 0 || oVisited.getValue(oCurrent.x(), oCurrent.y() + 1, oCurrent.z() - 1) == 3)) {
              oTemp = oCurrent;
              oTemp.y(oTemp.y() + 1);
              oTemp.z(oTemp.z() - 1);
              oQueue.push(oTemp);
            }
            // x+1 y+1 z-1
            if (oCurrent.x() + 1 < (int)this->m_iSizeX && oCurrent.y() + 1 < (int)this->m_iSizeY && (oVisited.getValue(oCurrent.x() + 1, oCurrent.y() + 1, oCurrent.z() - 1) == 0 || oVisited.getValue(oCurrent.x() + 1, oCurrent.y() + 1, oCurrent.z() - 1) == 3)) {
              oTemp = oCurrent;
              oTemp.x(oTemp.x() + 1);
              oTemp.y(oTemp.y() + 1);
              oTemp.z(oTemp.z() - 1);
              oQueue.push(oTemp);
            }

            // x-1 y+0 z-1
            if (oCurrent.x() - 1 > 0 && (oVisited.getValue(oCurrent.x() - 1, oCurrent.y(), oCurrent.z() - 1) == 0 || oVisited.getValue(oCurrent.x() - 1, oCurrent.y(), oCurrent.z() - 1) == 3)) {
              oTemp = oCurrent;
              oTemp.x(oTemp.x() - 1);
              oTemp.z(oTemp.z() - 1);
              oQueue.push(oTemp);
            }
            // x+0 y+0 z-1
            if (oVisited.getValue(oCurrent.x(), oCurrent.y(), oCurrent.z() - 1) == 0 || oVisited.getValue(oCurrent.x(), oCurrent.y(), oCurrent.z() - 1) == 3) {
              oTemp = oCurrent;
              oTemp.z(oTemp.z() - 1);
              oQueue.push(oTemp);
            }
            // x+1 y+0 z-1
            if (oCurrent.x() + 1 < (int)this->m_iSizeX && (oVisited.getValue(oCurrent.x() + 1, oCurrent.y(), oCurrent.z() - 1) == 0 || oVisited.getValue(oCurrent.x() + 1, oCurrent.y(), oCurrent.z() - 1) == 3)) {
              oTemp = oCurrent;
              oTemp.x(oTemp.x() + 1);
              oTemp.z(oTemp.z() - 1);
              oQueue.push(oTemp);
            }

            // x-1 y-1 z-1
            if (oCurrent.x() - 1 > 0 && oCurrent.y() - 1 > 0  && (oVisited.getValue(oCurrent.x() - 1, oCurrent.y(), oCurrent.z() - 1) == 0 || oVisited.getValue(oCurrent.x() - 1, oCurrent.y(), oCurrent.z() - 1) == 3)) {
              oTemp = oCurrent;
              oTemp.x(oTemp.x() - 1);
              oTemp.y(oTemp.y() - 1);
              oTemp.z(oTemp.z() - 1);
              oQueue.push(oTemp);
            }
            // x+0 y-1 z-1
            if ((oCurrent.y() - 1 > 0 && oVisited.getValue(oCurrent.x(), oCurrent.y() - 1, oCurrent.z() - 1) == 0)
                || (oCurrent.y() - 1 > 0 && oVisited.getValue(oCurrent.x(), oCurrent.y() - 1, oCurrent.z() - 1) == 3)) {
              oTemp = oCurrent;
              oTemp.y(oTemp.y() - 1);
              oTemp.z(oTemp.z() - 1);
              oQueue.push(oTemp);
            }
            // x+1 y-1 z-1
            if (oCurrent.x() + 1 < (int)this->m_iSizeX && oCurrent.y() - 1 > 0 && (oVisited.getValue(oCurrent.x() + 1, oCurrent.y() - 1, oCurrent.z() - 1) == 0 || oVisited.getValue(oCurrent.x() + 1, oCurrent.y() - 1, oCurrent.z() - 1) == 3)) {
              oTemp = oCurrent;
              oTemp.x(oTemp.x() + 1);
              oTemp.y(oTemp.y() - 1);
              oTemp.z(oTemp.z() - 1);
              oQueue.push(oTemp);
            }

          }

        } else
          oVisited.setValue(oCurrent.x(), oCurrent.y(), oCurrent.z(), 2);
      }
    }
  } catch (int e) {
    SVTLBBO << "Floodfill canceled..." << endl;
    return;
  }

  unsigned int iNum = this->m_iSizeX * this->m_iSizeY * this->m_iSizeZ;
  bool bOnlyNull = true;
  for (unsigned i = 0; i < iNum; i++) {
    if (oVisited.at(i) != 0) {
      bOnlyNull = false;
      if (oVisited.at(i) == 2 || oVisited.at(i) == 3)
        oVisited.setAt(i, T(0));
    }
  }

  if (bOnlyNull == false && fGaussian > 0.0) {
    SVTLBBO << "   Convolving mask with Gaussian kernel..." << endl;

    svt_volume<T> oGaussian;
    oGaussian.createGaussian(fGaussian, 3.0);
    oVisited.convolve(oGaussian, true);
    oVisited.normalize();
  }

  SVTLBBO << "   Applying mask..." << endl;

  for (unsigned i = 0; i < iNum; i++)
    this->setAt(i, this->at(i) * oVisited.at(i));

  this->m_bChanged = true;

  SVTLBBO << "   Done." << endl;
};

/**
 * Non-recursive flood-fill algorithm. The algorithm only fills the voxels that are above the threshold (and that are connected to the start voxel) with the specified value.
 * \param iStartX x index of the start voxel
 * \param iStartY y index of the start voxel
 * \param iStartZ z index of the start voxel
 * \param fTreshold threshold for the floodfill
 * \param fGaussian sigma of the gaussian the mask is convoluted with (if 0, no blurring of the mask is applied)
 */
template<class T>
void svt_volume<T>::floodfill(unsigned int iStartX, unsigned int iStartY, unsigned int iStartZ, T fThreshold, T fFillValue)
{
  svt_volume<T> oVisited(this->m_iSizeX, this->m_iSizeY, this->m_iSizeZ);
  oVisited = (T)(0);

  stack< svt_vector4<int> > oQueue;
  svt_vector4<int> oStart(iStartX, iStartY, iStartZ);
  oQueue.push(oStart);

  if (this->getValue(iStartX, iStartY, iStartZ) < fThreshold)
    SVTLBBO << "Warning: Start voxel is lower than threshold in floodfill!" << endl;

  svt_vector4<int> oCurrent;
  svt_vector4<int> oTemp;

  SVTLBBO << "Floodfill: Generating mask..." << endl;

  while (oQueue.size() > 0) {
    oCurrent = oQueue.top();
    oQueue.pop();

    if (oVisited.getValue(oCurrent.x(), oCurrent.y(), oCurrent.z()) == 0) {
      if (this->getValue(oCurrent.x(), oCurrent.y(), oCurrent.z()) > fThreshold) {
        oVisited.setValue(oCurrent.x(), oCurrent.y(), oCurrent.z(), 1);

        // z = 0

        // x-1 y+1 z+0
        if (oCurrent.x() - 1 >= 0 && oCurrent.y() + 1 < (int)this->m_iSizeY) {
          oTemp = oCurrent;
          oTemp.x(oTemp.x() - 1);
          oTemp.y(oTemp.y() + 1);
          oQueue.push(oTemp);
        }
        // x+0 y+1 z+0
        if (oCurrent.y() + 1 < (int)this->m_iSizeY) {
          oTemp = oCurrent;
          oTemp.y(oTemp.y() + 1);
          oQueue.push(oTemp);
        }
        // x+1 y+1 z+0
        if (oCurrent.x() + 1 < (int)this->m_iSizeX && oCurrent.y() + 1 < (int)this->m_iSizeY) {
          oTemp = oCurrent;
          oTemp.x(oTemp.x() + 1);
          oTemp.y(oTemp.y() + 1);
          oQueue.push(oTemp);
        }

        // x-1 y+0 z+0
        if (oCurrent.x() - 1 >= 0) {
          oTemp = oCurrent;
          oTemp.x(oTemp.x() - 1);
          oQueue.push(oTemp);
        }
        // x+1 y+0 z+0
        if (oCurrent.x() + 1 < (int)this->m_iSizeX) {
          oTemp = oCurrent;
          oTemp.x(oTemp.x() + 1);
          oQueue.push(oTemp);
        }

        // x-1 y-1 z+0
        if (oCurrent.x() - 1 >= 0 && oCurrent.y() - 1 >= 0) {
          oTemp = oCurrent;
          oTemp.x(oTemp.x() - 1);
          oTemp.y(oTemp.y() - 1);
          oQueue.push(oTemp);
        }
        // x+0 y-1 z+0
        if (oCurrent.y() - 1 >= 0) {
          oTemp = oCurrent;
          oTemp.y(oTemp.y() - 1);
          oQueue.push(oTemp);
        }
        // x+1 y-1 z+0
        if (oCurrent.x() + 1 < (int)this->m_iSizeX && oCurrent.y() - 1 >= 0) {
          oTemp = oCurrent;
          oTemp.x(oTemp.x() + 1);
          oTemp.y(oTemp.y() - 1);
          oQueue.push(oTemp);
        }

        // z = 1
        if (oCurrent.z() + 1 < (int)this->m_iSizeZ) {

          // x-1 y+1 z+1
          if (oCurrent.x() - 1 >= 0 && oCurrent.y() + 1 < (int)this->m_iSizeY) {
            oTemp = oCurrent;
            oTemp.x(oTemp.x() - 1);
            oTemp.y(oTemp.y() + 1);
            oTemp.z(oTemp.z() + 1);
            oQueue.push(oTemp);
          }
          // x+0 y+1 z+1
          if (oCurrent.y() + 1 < (int)this->m_iSizeY) {
            oTemp = oCurrent;
            oTemp.y(oTemp.y() + 1);
            oTemp.z(oTemp.z() + 1);
            oQueue.push(oTemp);
          }
          // x+1 y+1 z+1
          if (oCurrent.x() + 1 < (int)this->m_iSizeX && oCurrent.y() + 1 < (int)this->m_iSizeY) {
            oTemp = oCurrent;
            oTemp.x(oTemp.x() + 1);
            oTemp.y(oTemp.y() + 1);
            oTemp.z(oTemp.z() + 1);
            oQueue.push(oTemp);
          }

          // x-1 y+0 z+1
          if (oCurrent.x() - 1 >= 0) {
            oTemp = oCurrent;
            oTemp.x(oTemp.x() - 1);
            oTemp.z(oTemp.z() + 1);
            oQueue.push(oTemp);
          }
          // x+0 y+0 z+1
          {
            oTemp = oCurrent;
            oTemp.z(oTemp.z() + 1);
            oQueue.push(oTemp);
          }
          // x+1 y+0 z+1
          if (oCurrent.x() + 1 < (int)this->m_iSizeX) {
            oTemp = oCurrent;
            oTemp.x(oTemp.x() + 1);
            oTemp.z(oTemp.z() + 1);
            oQueue.push(oTemp);
          }

          // x-1 y-1 z+1
          if (oCurrent.x() - 1 >= 0 && oCurrent.y() - 1 >= 0) {
            oTemp = oCurrent;
            oTemp.x(oTemp.x() - 1);
            oTemp.y(oTemp.y() - 1);
            oTemp.z(oTemp.z() + 1);
            oQueue.push(oTemp);
          }
          // x+0 y-1 z+1
          if (oCurrent.y() - 1 >= 0) {
            oTemp = oCurrent;
            oTemp.y(oTemp.y() - 1);
            oTemp.z(oTemp.z() + 1);
            oQueue.push(oTemp);
          }
          // x+1 y-1 z+1
          if (oCurrent.x() + 1 < (int)this->m_iSizeX && oCurrent.y() - 1 >= 0) {
            oTemp = oCurrent;
            oTemp.x(oTemp.x() + 1);
            oTemp.y(oTemp.y() - 1);
            oTemp.z(oTemp.z() + 1);
            oQueue.push(oTemp);
          }

        }

        // z = -1
        if (oCurrent.z() - 1 >= 0) {

          // x-1 y+1 z-1
          if (oCurrent.x() - 1 >= 0 && oCurrent.y() + 1 < (int)this->m_iSizeY) {
            oTemp = oCurrent;
            oTemp.x(oTemp.x() - 1);
            oTemp.y(oTemp.y() + 1);
            oTemp.z(oTemp.z() - 1);
            oQueue.push(oTemp);
          }
          // x+0 y+1 z-1
          if (oCurrent.y() + 1 < (int)this->m_iSizeY) {
            oTemp = oCurrent;
            oTemp.y(oTemp.y() + 1);
            oTemp.z(oTemp.z() - 1);
            oQueue.push(oTemp);
          }
          // x+1 y+1 z-1
          if (oCurrent.x() + 1 < (int)this->m_iSizeX && oCurrent.y() + 1 < (int)this->m_iSizeY) {
            oTemp = oCurrent;
            oTemp.x(oTemp.x() + 1);
            oTemp.y(oTemp.y() + 1);
            oTemp.z(oTemp.z() - 1);
            oQueue.push(oTemp);
          }

          // x-1 y+0 z-1
          if (oCurrent.x() - 1 >= 0) {
            oTemp = oCurrent;
            oTemp.x(oTemp.x() - 1);
            oTemp.z(oTemp.z() - 1);
            oQueue.push(oTemp);
          }
          // x+0 y+0 z-1
          {
            oTemp = oCurrent;
            oTemp.z(oTemp.z() - 1);
            oQueue.push(oTemp);
          }
          // x+1 y+0 z-1
          if (oCurrent.x() + 1 < (int)this->m_iSizeX) {
            oTemp = oCurrent;
            oTemp.x(oTemp.x() + 1);
            oTemp.z(oTemp.z() - 1);
            oQueue.push(oTemp);
          }

          // x-1 y-1 z-1
          if (oCurrent.x() - 1 >= 0 && oCurrent.y() - 1 >= 0) {
            oTemp = oCurrent;
            oTemp.x(oTemp.x() - 1);
            oTemp.y(oTemp.y() - 1);
            oTemp.z(oTemp.z() - 1);
            oQueue.push(oTemp);
          }
          // x+0 y-1 z-1
          if (oCurrent.y() - 1 >= 0) {
            oTemp = oCurrent;
            oTemp.y(oTemp.y() - 1);
            oTemp.z(oTemp.z() - 1);
            oQueue.push(oTemp);
          }
          // x+1 y-1 z-1
          if (oCurrent.x() + 1 < (int)this->m_iSizeX && oCurrent.y() - 1 >= 0) {
            oTemp = oCurrent;
            oTemp.x(oTemp.x() + 1);
            oTemp.y(oTemp.y() - 1);
            oTemp.z(oTemp.z() - 1);
            oQueue.push(oTemp);
          }

        }

      } else
        oVisited.setValue(oCurrent.x(), oCurrent.y(), oCurrent.z(), 2);
    }
  }

  SVTLBBO << "   Applying mask..." << endl;

  unsigned int iNum = this->m_iSizeX * this->m_iSizeY * this->m_iSizeZ;
  for (unsigned i = 0; i < iNum; i++)
    if (oVisited.at(i) == 1)
      this->setAt(i, fFillValue);

  this->m_bChanged = true;

  SVTLBBO << "   Done." << endl;
};

///////////////////////////////////////////////////////////////////////////////
// MRC File Format
///////////////////////////////////////////////////////////////////////////////

/* Header for MRC file format */
/* Adapted from EMAN2 nightly build 12/05/05 by Steve Ludtke */

/* Label sizes */
enum { MRC_LABEL_SIZE = 80 };
enum { MRC_USER       = 25 };
enum { MRC_NUM_LABELS = 10 };

/*  The different modes supported by the MRC format. */
enum { MRC_MODE_char          = 0 };
enum { MRC_MODE_short         = 1 };
enum { MRC_MODE_float         = 2 };
enum { MRC_MODE_short_COMPLEX = 3 };
enum { MRC_MODE_float_COMPLEX = 4 };
enum { MRC_MODE_long          = 6 };

/* Data order */
enum { MRC_DATA_ORDER_XYZ,
       MRC_DATA_ORDER_XZY,
       MRC_DATA_ORDER_YXZ,
       MRC_DATA_ORDER_YZX,
       MRC_DATA_ORDER_ZXY,
       MRC_DATA_ORDER_ZYX
     };

#define EPSILON 0.02

/**
 * Simple internal routine that calculates the x,y,z index of a voxel
 * \param iIndex number of the voxel
 * \param iDataOrder data order
 * \param rX reference to unsigned int
 * \param rY reference to unsigned int
 * \param rZ reference to unsigned int
 */
template<class T>
void svt_volume<T>::calc_xyz(unsigned int iIndex, unsigned int iDataOrder, unsigned int iMaxX, unsigned int iMaxY, unsigned int, unsigned int &rX, unsigned int &rY, unsigned int &rZ)
{
  int xo = 0;
  int yo = 0;
  int zo = 0;

  zo = iIndex / (iMaxX * iMaxY);
  yo = (iIndex - zo * iMaxX * iMaxY) / iMaxX;
  xo = iIndex - zo * iMaxX * iMaxY - yo * iMaxX;

  switch (iDataOrder) {

    case MRC_DATA_ORDER_XYZ:
      rX  = xo;
      rY  = yo;
      rZ  = zo;
      break;

    case MRC_DATA_ORDER_XZY:
      rX  = xo;
      rY  = zo;
      rZ  = yo;
      break;

    case MRC_DATA_ORDER_YXZ:
      rX  = yo;
      rY  = xo;
      rZ  = zo;
      break;

    case MRC_DATA_ORDER_YZX:
      rX  = yo;
      rY  = zo;
      rZ  = xo;
      break;

    case MRC_DATA_ORDER_ZXY:
      rX  = zo;
      rY  = xo;
      rZ  = yo;
      break;

    case MRC_DATA_ORDER_ZYX:
      rX  = zo;
      rY  = yo;
      rZ  = xo;
      break;

    default:
      SVTLBBO << "calc_xyz_index_map: Unknown data order" << endl;
  }
};


/**
 * Simple internal routine that calculates the index of a voxel, depending on the data order
 * \param iIndex number of the voxel
 * \param iDataOrder data order
 */
template<class T>
unsigned int svt_volume<T>::calc_xyz_index(unsigned int iIndex, unsigned int iDataOrder)
{
  int xo = 0;
  int yo = 0;
  int zo = 0;
  int x = 0;
  int y = 0;
  int z = 0;
  int idx;
  int nx = 0;
  int ny = 0;
  int nz = 0;

  zo = iIndex / (m_iSizeX * m_iSizeY);
  yo = (iIndex - zo * m_iSizeX * m_iSizeY) / m_iSizeX;
  xo = iIndex - zo * m_iSizeX * m_iSizeY - yo * m_iSizeX;

  switch (iDataOrder) {

    case MRC_DATA_ORDER_XYZ:
      x  = xo;
      y  = yo;
      z  = zo;
      nx = m_iSizeX;
      ny = m_iSizeY;
      nz = m_iSizeZ;
      break;

    case MRC_DATA_ORDER_XZY:
      x  = xo;
      y  = zo;
      z  = yo;
      nx = m_iSizeX;
      ny = m_iSizeZ;
      nz = m_iSizeY;
      break;

    case MRC_DATA_ORDER_YXZ:
      x  = yo;
      y  = xo;
      z  = zo;
      nx = m_iSizeY;
      ny = m_iSizeX;
      nz = m_iSizeZ;
      break;

    case MRC_DATA_ORDER_YZX:
      x  = yo;
      y  = zo;
      z  = xo;
      nx = m_iSizeY;
      ny = m_iSizeZ;
      nz = m_iSizeX;
      break;

    case MRC_DATA_ORDER_ZXY:
      x  = zo;
      y  = xo;
      z  = yo;
      nx = m_iSizeZ;
      ny = m_iSizeX;
      nz = m_iSizeY;
      break;

    case MRC_DATA_ORDER_ZYX:
      x  = zo;
      y  = yo;
      z  = xo;
      nx = m_iSizeZ;
      ny = m_iSizeY;
      nz = m_iSizeX;
      break;

    default:
      SVTLBBO << "calc_xyz_index_map: Unknown data order" << endl;
  }

  idx = x + nx * y + nx * ny * z;

  return idx;
};

/*******************************************************************************
Protected Data Structures
*******************************************************************************/



///////////////////////////////////////////////////////////////////////////////
// SVT_COLUMN_READER
///////////////////////////////////////////////////////////////////////////////

#define MAXLINE 2048

/**
 * Column-based file format reader.
 */
class svt_column_reader
{
  protected:

    char m_pLine[MAXLINE];
    FILE *m_pFile;

    char m_pString[MAXLINE];

  public:

    /**
     * Constructor
     */
    svt_column_reader(const char *pFilename) :
      m_pFile(NULL)
    {
      m_pFile = fopen(pFilename, "r");
    };
    ~svt_column_reader()
    {
      if (m_pFile != NULL)
        fclose(m_pFile);
    };

    /**
     * Read next line.
     */
    inline bool readLine()
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
    inline char *extractString(unsigned int iStart, unsigned int iEnd)
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
    inline char extractChar(unsigned int iCol)
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
    inline Real32 extractReal32(unsigned int iStart, unsigned int iEnd)
    {
      Real32 fVal = 0.0f;

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
    inline int extractInt(unsigned int iStart, unsigned int iEnd)
    {
      int iVal = 0;

      char *pString = extractString(iStart, iEnd);
      iVal = atoi(pString);

      return iVal;
    };

    /**
     * Get length of line
     */
    inline unsigned int getLength() const
    {
      return strlen(m_pLine);
    };

    /**
     * EOF test.
     */
    inline bool eof()
    {
      if (m_pFile != NULL)
        return (bool)(feof(m_pFile));
      else
        return true;
    };

};


///////////////////////////////////////////////////////////////////////////////
// SVT_ATOM
///////////////////////////////////////////////////////////////////////////////

typedef struct {
  //the three letter codes for the amino acids
  char m_a3LetCode[5];
  // the one letter code
  char m_c1LetCode;

} svt_resname;

// the list with the aa long and short names
static vector<svt_resname> m_oAAList;

/**
* basic class containing the atom information
*/
class svt_point_cloud_atom
{
  protected:

    // to which model does this atom belong
    unsigned int m_iModel;
    // PDB-file index of the atom
    unsigned int m_iPDBIndex;

    // atom name (i.e. "C")
    char   m_aName[5];
    // remoteness indicator
    char   m_cRemoteness;
    // branch
    char   m_cBranch;
    // radius
    Real32 m_fRadius;
    // alternate location identifier
    char   m_cAltLoc;
    // atom resname (i.e. "ALA")
    char   m_aResname[5];
    // atom short resname (i.e. "A")
    char   m_cShortResname;
    // chain id
    char   m_cChainID;
    // ordinal chain id (note it's an integer!)
    int    m_iOrdChainID;
    // residue sequence number
    int    m_iResSeq;
    // ordinal res sequence number
    int    m_iOrdResSeq;
    // iCode
    char   m_cICode;
    // the occupancy
    Real32 m_fOccupancy;
    // the temperature factor
    Real32 m_fTempFact;
    // note
    char   m_pNote[4];
    // segment id
    char   m_pSegID[5];
    // element symbol
    char   m_pElement[3];
    // charge
    char   m_pCharge[3];

    // secondary structure information
    // H   Alpha helix
    // G   3-10 helix
    // I   PI-helix
    // E   Extended conformation
    // B   Isolated bridge
    // T   Turn
    // C   Coil (none of the above)
    // N   Information not available
    char   m_cSecStruct;

    // number of residues in the secondary structure this atom belongs to
    int    m_iSecStructNumResidues;

    // relative atomic mass
    Real64 m_fMass;

    // is this atom part of a water molecule?
    bool   m_bWater;

    //is this atom an "hetatm" record in the pdb
    bool m_bHetAtm;

    // bond list
    vector< svt_point_cloud_bond * > m_aBondList;
    // index-based bond list
    vector< unsigned int > m_aBondListI;

    // Is Atom selected
    bool m_bIsSelected;


  public:

    /**
     * Constructor
     */
    svt_point_cloud_atom();
    /**
     * Destructor
     */
    virtual ~svt_point_cloud_atom();

    /**
     * Set the model this atom belongs to
     * \param iModel number of model
     */
    void setModel(unsigned int iModel);
    /**
     * Get the model this atom belongs to
     * \return number of model
     */
    unsigned int getModel() const;

    /**
     * Set the PDB-file index of the atom
     * \param iPDBIndex index of the atom
     */
    void setPDBIndex(unsigned int iPDBIndex);
    /**
     * Get the PDB-file index of the atom
     * \return index of the atom
     */
    unsigned int getPDBIndex() const;

    /**
     * set the atom name
     * \param pName pointer to the char array with the name information
     */
    void setName(const char *pType);
    /**
     * get the atom name
     * \return pointer to the char array with the name information
     */
    const char *getName() const;
    /**
     * get the atom remoteness indicator
     * \return char with the remoteness information (greek letters, A=alpha, B=beta,...)
     */
    char getRemoteness() const;
    /**
     * set the atom remoteness indicator
     * \param cRemoteness char with the remoteness information (greek letters, A=alpha, B=beta,...)
     */
    void setRemoteness(char cRemoteness);
    /**
     * get the atom branch information
     * \return char with the branch information
     */
    char getBranch() const;
    /**
     * set the atom branch information
     * \param cBranch char with the branch information
     */
    void setBranch(char cBranch);
    /**
     * set the atom alternate location indicator
     * \param cAltLoc character with the alternate location indicator
     */
    void setAltLoc(char cAltLoc);
    /**
     * get the atom alternate location indicator
     * \return character with the alternate location indicator
     */
    char getAltLoc() const;
    /**
     * set the atom resname
     * \param pResname pointer to the char array with the residue name information
     */
    void setResname(const char *pResname);
    /**
     * get the atom resname
     * \return pointer to the char array with the residue name information
     */
    const char *getResname() const;
    /**
    * set the short atom resname
    * \param pResname a char with the short "1letter" residue name information
    */
    void setShortResname(const char cResname);
    /**
     * get the short atom resname
     * \return a char with the short "1letter" residue name information
     */
    char getShortResname() const;
    /**
      * set the chain id
      * \param cChainID character with the chain id
      */
    void setChainID(char cChainID);
    /**
     * get the chain id
     * \return character with the chain id
     */
    char getChainID() const;
    /**
     * set the chain id
     * \param cOrdChainID character with the chain id
     */
    void setOrdChainID(int iOrdChainID);
    /**
     * get the chain id
     * \return int with the chain id
     */
    int getOrdChainID() const;
    /**
     * set the residue sequence number
     * \param iResSeq residue sequence number
     */
    void setResidueSeq(int iResSeq);
    /**
     * get the residue sequence number
     * \return residue sequence number
     */
    int getResidueSeq() const;
    /**
     * set the ordinal residue sequence number
     * \param iOrdResSeq ordinal residue sequence number
     */
    void setOrdResidueSeq(int iOrdResSeq);
    /**
     * get the ordinal residue sequence number
     * \return ordinal residue sequence number
     */
    int getOrdResidueSeq() const;
    /**
     * set the iCode (code for insertion of residues)
     * \param cICode char with the iCode
     */
    void setICode(char cICode);
    /**
     * get the iCode (code for insertion of residues)
     * \return char with the iCode
     */
    char getICode() const;
    /**
     * set the occupancy
     * \param fOccupancy the occupancy
     */
    void setOccupancy(float fOccupancy);
    /**
     * get the occupancy
     * \return the occupancy
     */
    float getOccupancy() const;
    /**
     * set the temperature factor
     * \param fTempFact the temperature factor
     */
    void setTempFact(float fTempFact);
    /**
     * get the temperature factor
     * \return the temperature factor
     */
    float getTempFact() const;
    /**
     * set the note
     * \param pNote pointer to char array with at least 3 chars for the note
     */
    void setNote(const char *pNote);
    /**
     * get the note
     * \return pointer to char array with the note
     */
    const char *getNote() const;
    /**
     * set the segment id
     * \param pSegID pointer to char array with at least 4 chars for the segment id
     */
    void setSegmentID(const char *pSegID);
    /**
     * get the segment id
     * \return pointer to char array with the segment id
     */
    const char *getSegmentID() const;
    /**
     * set the element symbol
     * \param pElement pointer to char array with at least 2 chars for the element symbol
     */
    void setElement(const char *pElement);
    /**
     * get the element symbol
     * \return pointer to char array with the element symbol
     */
    const char *getElement() const;
    /**
     * set the charge
     * \param pCharge pointer to char array with at least 2 chars for the charge
     */
    void setCharge(const char *pCharge);
    /**
     * get the charge
     * \return pointer to char array with the charge
     */
    const char *getCharge() const;

    /**
     * is the atom a hydrogen atom?
     * \return true if the the atom is an hydrogen atom
     */
    bool isHydrogen() const;
    /**
     * is the atom a QPDB codebook vector?
     * \return true if the the atom is a CV
     */
    bool isQPDB() const;
    /**
     * is the atom part of a water molecule?
     * \return true if the atom os part of a water molecule
     */
    bool isWater() const;
    /**
     * is the atom a carbon alpha?
     * \return true if the atom a CA
     */
    bool isCA() const;
    /**
     * is the atom on the Backbone?
     * \return true if the atom on the Backbone
     */
    bool isBackbone() const;
    /**
     * is the atom part of a nucleotide?
     * \return true if the atom is part of a nucleotide
     */
    bool isNucleotide() const;

    /**
     * set the relative atomic mass
     * \param fMass relative atomic mass
     */
    void setMass(Real64 fMass);
    /**
     * get the relative atomic mass
     * \return relative atomic mass
     */
    Real64 getMass() const;

    /**
     * adjust the atomic mass based on a (simple) periodic table.
     * ToDo: Full periodic table
     */
    void adjustMass();

    /**
     * set the secondary structure information for this atom
     * \param cSecStruct secondary structure information
     */
    void setSecStruct(char cSecStruct);

    /**
     * get the secondary structure information for this atom
     * \return secondary structure information
     */
    char getSecStruct();

    /**
     * set the length in number of residues of the secondary structure element to which this atom belongs
     * \param iSecStructNumResidues number of residues
     */
    void setSecStructNumResidues(int iSecStructNumResidues);

    /**
     * get the length in number of residues of the secondary structure element to which this atom belongs
     * \return number of residues
     */
    int getSecStructNumResidues();

    /**
     * Get the van der waals radius
     * \return van der waals radius value (in Angstroem)
     */
    Real64 getVDWRadius() const;

    /**
     * Add a bond to the bond list of this atom
     * \param rBond pointer to the svt_point_cloud_bond object which should be added to the bond list
     */
    void addToBondList(svt_point_cloud_bond *pBond, unsigned int iIndex);
    /**
     * Remove a bond from the bond list of this atom
     * \param iIndex index of the svt_point_cloud_bond object which should be deleted from the bond list
     */
    void delFromBondList(unsigned int iIndex);

    /**
     * Get bond list
     */
    vector< unsigned int > getBondList();
    /**
     * Adjust bond list - a bond was erased from the global list and now all bond indexes higher than a certain number have to be decresed by one
     */
    void adjustBondList(unsigned int iIndex);

    /**
     * deletes all bonds in the bond list
     */
    void delBondList();


    /**
     * set Selected
     */
    void setSelected(bool bSelected);

    /**
     * get Selected
     */
    bool getSelected();

    /**
     * set the record name : is atom of class hetatm
     */
    void setHetAtm(bool bHetAtm);

    /**
     * get the record name : is atom of class hetatm
     */
    bool getHetAtm();


};



///////////////////////////////////////////////////////////////////////////////
// SVT_POINT_CLOUD
///////////////////////////////////////////////////////////////////////////////
/***************************************************************************
                          svt_point_cloud
                          ---------------
    begin                : 04/12/2005
    author               : Stefan Birmanns
    email                : Stefan.Birmanns@uth.tmc.edu
 ***************************************************************************/

#define NOMATCH 100000000

/**
 * Helper class
 */
template<class T> class svt_neighbor
{
  protected:

    T m_fScore;
    unsigned int m_iIndexA;
    unsigned int m_iIndexB;

  public:

    svt_neighbor(T fScore, unsigned int iIndexA, unsigned int iIndexB)
    {
      m_fScore = fScore;
      m_iIndexA = iIndexA;
      m_iIndexB = iIndexB;
    };

    inline T getScore() const
    {
      return m_fScore;
    };
    inline unsigned int getIndexA() const
    {
      return m_iIndexA;
    };
    inline unsigned int getIndexB() const
    {
      return m_iIndexB;
    };

    inline bool operator<(const svt_neighbor<T> &rR) const
    {
      return m_fScore < rR.m_fScore;
    };
};


/**
 * A matching helper class that encapsulates a single result of the matching process
 */
template<class T> class svt_matchResult
{
  protected:

    T m_fScore;
    svt_matrix4<T> m_oMatrix;
    vector<int> m_aModelMatch;
    vector<int> m_aSceneMatch;
    vector<int> m_aMatch;

  public:

    svt_matchResult(T fScore, svt_matrix4<T> oMatrix, vector<int> &aModelMatch, vector<int> &aSceneMatch)
    {
      m_fScore = fScore;
      m_oMatrix = oMatrix;
      m_aModelMatch = aModelMatch;
      m_aSceneMatch = aSceneMatch;
    };

    inline T getScore() const
    {
      return m_fScore;
    };
    inline void setScore(T fScore)
    {
      m_fScore = fScore;
    };
    inline svt_matrix4<T> getMatrix() const
    {
      return m_oMatrix;
    };
    inline void setMatrix(svt_matrix4<Real64> oMatrix)
    {
      m_oMatrix = oMatrix;
    };
    inline vector<int> &getModelMatch()
    {
      return m_aModelMatch;
    };
    inline vector<int> &getSceneMatch()
    {
      return m_aSceneMatch;
    };
    vector<int> &getMatch()
    {
      if (m_aMatch.size() == 0)
        makeSorted();

      return m_aMatch;
    };
    inline unsigned int compareMatch(svt_matchResult<T> &rOther)
    {
      if (m_aMatch.size() == 0)
        makeSorted();

      vector<int> &rOtherM = rOther.getMatch();
      unsigned int iNum = m_aMatch.size();
      unsigned int iDiff = 0;

      if (rOtherM.size() != m_aMatch.size()) {
        SVTLBBO << "ERROR: rOtherM.size(): " << rOtherM.size() << " m_aMatch.size(): " << m_aMatch.size() << endl;
        exit(1);
      }

      for (unsigned int i = 0; i < iNum; i++)
        if (rOtherM[i] != m_aMatch[i] && rOtherM[i] != NOMATCH && m_aMatch[i] != NOMATCH)
          iDiff++;

      return iDiff;
    };

    inline void printMatch()
    {
      if (m_aMatch.size() == 0)
        makeSorted();

      SVTLBBO;
      for (unsigned int i = 0; i < m_aMatch.size(); i++)
        if (m_aMatch[i] != NOMATCH)
          if (i != m_aMatch.size() - 1)
            printf("%02i,", m_aMatch[i]);
          else
            printf("%02i", m_aMatch[i]);
        else
          printf(" - ");
      cout << endl;
    };

    inline bool operator<(const svt_matchResult<T> &rR) const
    {
      return m_fScore < rR.m_fScore;
    };

  private:

    /**
     * Generate the sorted match
     */
    void makeSorted()
    {
      if (m_aMatch.size() == 0) {
        m_aMatch = vector<int>(m_aModelMatch.size());
        for (unsigned int i = 0; i < m_aModelMatch.size(); i++)
          m_aMatch[m_aModelMatch[i]] = m_aSceneMatch[i];
      }
    };
};


/**
 * A point cloud class.
 * \author Stefan Birmanns
 */
template<class T> class svt_point_cloud : public svt_sampled< T >
{
  protected:

    vector< vector< T > > m_aPoints;

    // tree pruning parameter
    Real64 m_fEpsilon;
    // tolorance distance
    Real64 m_fLambda;
    // ranking distance
    Real64 m_fGamma;
    // size of the matching zone
    unsigned int m_iZoneSize;
    // number of runs with different anchor points
    unsigned int m_iRuns;

    // maximum number of wildcard matches
    unsigned int m_iMaxNoMatch;

    // number of seed/anchor matches evaluated during last matching run
    unsigned int m_iNumSeeds;

    // in case of simple==true only a very simple version of the matching algorithm is performed...
    bool m_bSimple;

    // uniqueness rmsd threshold for the solutions
    Real64 m_fUnique;

    // next point selection scheme
    bool m_bNextPointCOA;

    // timestep information for time-varying pointclouds
    unsigned int m_iTimestep;

    // penalty for the wildcards in the matching algo
    Real64 m_fSkipPenalty;

  public:

    /**
     * Constructor
     */
    svt_point_cloud() :
      m_aPoints(1),
      m_fEpsilon(30.0),
      m_fLambda(5.0),
      m_fGamma(1.0),
      m_iZoneSize(3),
      m_iRuns(1),
      m_iMaxNoMatch(3),
      m_iNumSeeds(0),
      m_bSimple(false),
      m_fUnique(1.0),
      m_bNextPointCOA(false),
      m_iTimestep(0),
      m_fSkipPenalty(1.0)
    {
    };

    ~svt_point_cloud();

    /**
     * Get all points in point cloud in an stl vector.
     * \return reference to vector of svt_vector4 objects
     */
    inline vector< T > &getPoints();

    /**
     * Add a point to point cloud.
     * \param rVec svt_vector4 object
     */
    inline void addPoint(T &rVec);
    /**
     * Add a point to point cloud.
     * \param fX x coord
     * \param fY y coord
     * \param fZ z coord
     */
    inline void addPoint(Real64 fX, Real64 fY, Real64 fZ);

    /**
     * Get a point out of point cloud.
     * \param iIndex index of point
     * \return reference to svt_vector4 object
     */
    inline T &getPoint(unsigned int iIndex);

    /**
     * Replace point in point cloud.
     * \param iIndex index of point
     * \param rVec svt_vector4 object
     */
    inline void setPoint(unsigned int iIndex, T &rVec);

    /**
     * Delete point in point cloud.
     * \param iIndex index of point
     */
    inline void delPoint(unsigned int iIndex);

    /**
     * Delete all points in point cloud
     */
    inline void delAllPoints();

    /**
     * Size of point cloud.
     * \return size of pc
     */
    inline unsigned int size() const;

    /**
     * Dereference operator (not range checked!).
     * \param iIndex index of point in point cloud
     */
    inline T &operator[](unsigned int iIndex);

    /**
     * Calculate center of atoms (COA with all masses = 1).
     * \return svt_vector4 with the COA
     */
    T coa();

    /**
     * Calculate the geometric center of the point cloud.
     * \return svt_vector4 with the geometric center
     */
    T geometricCenter();

    /**
     * calculate RMSD between this and another PC. The PCs must be matched already! To get the minimal RMSD please use align() first!
     * \param rPC second point cloud
     */
    Real64 rmsd(svt_point_cloud<T> &rPC);
    /**
     * calculate RMSD between this and another PC. The points are matched using the nearest neighbor relationship.
     * \param rPC second point cloud
     */
    Real64 rmsd_NN(svt_point_cloud<T> &rPC);
    /**
     * Calculate RMSD between this and another PC. The points are matched using the nearest neighbor relationship.
     * All the nearest neighbor distances are calculated and then sorted (slow!) and only the first N-percent are used for the matching and the following RMSD calculation.
     * The idea is that approximately N-percent of the points are outliers which would increase the RMSD significantly, although the overall deviation of the two point clouds is
     * actually small.
     * \param rPC second point cloud
     * \param fPercent percentage of neighbors that should be used for the rmsd calculation
     * \param rMatrix reference to svt_matrix4 object
     */
    Real64 rmsd_NN_Outliers(svt_point_cloud<T> &rPC, Real64 fPercent, svt_matrix4<Real64> *pMatrix = NULL);

    /**
     * calculate hausdorff distance between this and another PC. The PCs must be matched and aligned already!
     * \param rPC second point cloud
     */
    Real64 hausdorff(svt_point_cloud<T> &rPC);

    /**
     * find the nearest neighbor to a query point in the point cloud
     * \param rVec reference to svt_vector4 object - the query point
     * \return index to nearest point in point cloud
     */
    unsigned int nearestNeighbor(T &rVec);

    /**
     * calculate the average nearest neighbor distance
     * in order to reduce the complexity a random test is done and once the average stabilizes, the search is stopped.
     * \param fPrecision if average does not change by more than fPrecision between two iterations the calculation is stopped
     * \return average nearest neighbor distance
     */
    Real64 averageNNDistance(Real64 fPrecision);

    /**
     * calculate the maximum nearest neighbor distance
     * \return maximum nearest neighbor distance
     */
    Real64 maxNNDistance();

    /**
     * Align point cloud with another point cloud.
     * This function does not perform a matching procedure, it assumes that the ith point in the first
     * cloud is associated to the ith point in the second point cloud. It aligns both clouds so that they have a minimal RMSD.
     * \param rPC second point cloud
     */
    //void align( svt_point_cloud<T>& rPC );

    /**
     * Get the minimal coordinates of the point cloud - it will return a vector that has in each dimension the information about the minimal
     * coordinate it has found in the cloud.
     */
    T getMinCoord();
    /**
     * Get the maximal coordinates of the point cloud - it will return a vector that has in each dimension the information about the maximal
     * coordinate it has found in the cloud.
     */
    T getMaxCoord();

    //
    // clustering
    //

    /**
     * sample the object randomly and return a vector that refrects the probability distribution of the object
     */
    T sample();

    //
    // file i/o
    //

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
     * Write CSV file - comma separated values x,y,z. It will write only the coordinates of the current frame.
     * \param pFilename pointer to array of char with the filename
     */
    void writeCSV(const char *pFilename);

    //
    // matching functions
    //

    /**
     * Set tree pruning parameter
     * \param fEpsilon epsilon
     */
    void setEpsilon(Real64 fEpsilon);
    /**
     * Get tree pruning parameter.
     * \param fEpsilon epsilon
     */
    Real64 getEpsilon() const;
    /**
     * set tolorance distance for the anchor determination
     * \param fLambda lambda
     */
    void setLambda(Real64 fLambda);
    /**
     * get tolorance distance for the anchor determination
     * \return lambda
     */
    Real64 getLambda() const;
    /**
     * set nearest neighbor matching zone size
     * \param fGamma gamma
     */
    void setGamma(Real64 fGamma);
    /**
     * get nearest neighbor matching zone size
     * \return gamma
     */
    Real64 getGamma() const;

    /**
     * set the maximal size of the matching zone
     * \param iZoneSize
     */
    inline void setZoneSize(unsigned int iZoneSize);
    /**
     * get the maximal size of the matching zone
     * \return maximal size of matching zone
     */
    inline unsigned int getZoneSize() const;

    /**
     * set the maximal number of wildcard matches
     * \param iMaxNoMatch maximal number of wildcard matches
     */
    inline void setWildcards(unsigned int iMaxNoMatch);
    /**
     * get the maximal number of wildcard matches
     * \return maximal number of wildcard matches
     */
    inline unsigned int getWildcards() const;

    /**
     * set the penalty for wildcard matches
     * \param fSkipPenalty penalty for a single wildcard
     */
    inline void setSkipPenalty(Real64 fSkipPenalty);
    /**
     * get the penalty for wildcard matches
     * \return penalty for a single wildcard
     */
    inline Real64 getSkipPenalty() const;

    /**
     * if two solutions are very close only the one with the higher score is considered. The other solutions are removed.
     * \param minimal distance between solutions
     */
    inline void setUnique(Real64 fUnique);
    /**
     * if two solutions are very close only the one with the higher score is considered. The other solutions are removed.
     * \return minimal distance between solutions
     */
    inline Real64 getUnique() const;

    /**
     * set the number of runs. Each time a different set of anchor points get selected from the probe structure, beginning with the three points furthest away from each other and the COA.
     * \param iRuns number of runs
     */
    inline void setRuns(unsigned int iRuns);
    /**
     * get the number of runs. Each time a different set of anchor points get selected from the probe structure, beginning with the three points furthest away from each other and the COA.
     * \return number of runs
     */
    inline unsigned int getRuns() const;

    /**
     * Full or simple version of the matching algorithm to be performed.
     * \param bSimple if true only a reduced version is performed (faster but less accurate)
     */
    inline void setSimple(bool bSimple);
    /**
     * Full or simple version of the matching algorithm to be performed.
     * \return if true only a reduced version is performed (faster but less accurate)
     */
    inline bool getSimple() const;

    /**
     * Set the next point selection scheme to COA
     */
    inline void setNextPointCOA(bool bNextPointCOA);

    /**
     * match two point clouds.
     * This function will perform a full N->M matching based on an achor-based search.
     * \param rPC second point cloud
     * \param pMatch vector of unsigned ints with the indices of the points of second cloud to which the points of this cloud are matched. This vector will get erased and then filled during the matching.
     * \return copy of this cloud matched to rPC
     */
    void match(svt_point_cloud<T> &rPC, vector< svt_matchResult<Real64> > &rMatch);

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
    void match(svt_point_cloud<T> &rPC, vector< svt_matchResult<Real64> > &rMatch, unsigned int iAnchorA, unsigned int iAnchorB, unsigned int iAnchorC);

    /**
     * Get number of seed/anchor matches evaluated during last matching run.
     */
    unsigned int getNumSeeds() const;

    /**
     * Calculate the cross-correlation between this pointcloud and a volumetric object. The idea is to use the points in the pointcloud to sample the volume and to
     * form thereby ad reduce cross-correlation coefficient equation.
     * See also:
     * Stefan Birmanns and Willy Wriggers. Interactive Fitting Augmented by Force-Feedback and Virtual Reality. J. Struct. Biol., 2003, Vol. 144, pp. 123-131.
     * \param rPointMat reference to svt_matrix4 object with the transformation matrix of the pointcloud object
     * \param rVol reference to svt_volume object
     * \param rVolMat reference to svt_matrix4 object with the transformation matrix of the volumetric object
     * \return the cross-correlation
     */
    Real64 calcCC(svt_matrix4<Real64> &rPointMat, svt_volume<Real64> &rVol, svt_matrix4<Real64> &rVolMat);

  private:

    /**
     * sort eigenvectors according to their eigenvalues
     * \param pEigenvectors pointer to svt_matrix with the eigenvectors as column vectors
     * \param pEigenvalues pointer to svt_vector with the eigenvalues
     */
    void eigensort(svt_matrix4<Real64> &rEigenvectors, T &rEigenvalues);

  public:

    /**
     * least-squares fit of a model point set to a scene point set
     * \param rM vector of ints with the indices of the first (model) vectors that are used
     * \param rS vector of ints with the indices of the second (scene) vectors that are used
     * \param rModel vector of svt_vector4 with the model vectors
     * \param rScene vector of svt_vector4 with the scene vectors
     * \return matrix with optimal transformation of the model into the scene
     */
    svt_matrix4<Real64> kearsley(vector< int > &rM, vector< int > &rS, svt_point_cloud<T> &rModel, svt_point_cloud<T> &rScene);

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
     * \param rModel vector of svt_vector4 with the model vectors
     * \param rScene vector of svt_vector4 with the scene vectors
     * \return matrix with optimal transformation of the model into the scene
     */
    svt_matrix4<Real64> kearsley(svt_point_cloud<T> &rModel, svt_point_cloud<T> &rScene);

    //private:

    /**
     * optimize_simple - subroutine for the match() procedure
     */
    svt_matrix4<Real64> optimize_simple(vector<int> *pModelMatch, vector<int> *pSceneMatch, svt_point_cloud<T> &rModel, svt_point_cloud<T> &rScene, bool bInit);
    /**
     * optimize - subroutine for the match() procedure
     */
    svt_matrix4<Real64> optimize(vector<int> *pModelMatch, vector<int> *pSceneMatch, svt_point_cloud<T> &rModel, svt_point_cloud<T> &rScene, bool bInit);

    /**
     * internal convenience function - no real difference with hausdorff - only here we don't need to create/copy a svt_point_cloud
     */
    Real64 calcHausdorff(svt_point_cloud<T> &rModel, svt_point_cloud<T> &rScene, svt_matrix4<Real64> &rMatrix);
    /**
     * internal convenience function - no real difference with rmsd - only here we don't need to create/copy a svt_point_cloud
     */
    Real64 calcRMSD(svt_point_cloud<T> &rModel, svt_matrix4<Real64> oMatA, svt_matrix4<Real64> oMatB);
    /**
     * internal convenience function - no real difference with rmsd - only here we don't need to create/copy a svt_point_cloud
     */
    Real64 calcRMSD(svt_point_cloud<T> &rModel, svt_point_cloud<T> &rScene, svt_matrix4<Real64> &rMatrix);
    /**
     * internal convenience function - no real difference with rmsd - only here we don't need to create/copy a svt_point_cloud
     * plus here we take the matching into account! If NOMATCH this point will not used
     */
    Real64 calcRMSD(svt_point_cloud<T> &rModel, svt_point_cloud<T> &rScene, svt_matrix4<Real64> &rMatrix, vector<int> *pModelMatch);
    /**
     * internal convenience function - no real difference with rmsd - only here we don't need to create/copy a svt_point_cloud
     * plus here we take the matching into account! If NOMATCH this point will not be used.
     */
    Real64 calcRMSD(svt_point_cloud<T> &rModel, svt_point_cloud<T> &rScene, svt_matrix4<Real64> &rMatrix, vector<int> *pModelMatch, vector<int> *pSceneMatch);

    /**
     * Calculate the Pair Distribution Function - improved version without multiple function calls
     * \param fWidth is width of a bin (held constant between the range Dmin and Dmax)
     */
    void calcPairDistribution(Real64 fWidth = 1.0, bool bShowProgress = true);

    /**
     * Calculate the Pair Distribution Function
     * \param fWidth is width of a bin (held constant between the range Dmin and Dmax)
     */
    void calcPairDistribution2(Real64 fWidth = 1.0, bool bShowProgress = true);



    /**
     * Compute the distance matrix between all points in the pointcloud
     * \return svt_matrix object with the distances
     */
    svt_matrix<Real64> getDistanceMat();

    /**
     * \name Timesteps for time-varying point clouds
     */
    //@{

    /**
     * set the current timestep
     * \param iTimestep the new timestep
     */
    void setTimestep(unsigned int iTimestep);
    /**
     * get the current timestep
     * \return the timestep
     */
    int getTimestep();
    /**
     * get the maximum timestep
     * \return the maximum timestep
     */
    int getMaxTimestep();
    /**
     * Add a new timestep.
     */
    void addTimestep();

    //@}
};
/***************************************************************************
                          svt_point_cloud_match
                          ---------------------
    begin                : 02/10/2006
    author               : Stefan Birmanns
    email                : Stefan.Birmanns@uth.tmc.edu
 ***************************************************************************/

///////////////////////////////////////////////////////////////////////////////
// Implementation of point cloud matching functionality of svt_point_cloud
///////////////////////////////////////////////////////////////////////////////
// DO NOT INCLUDE DIRECTLY! File gets included by other header file!
///////////////////////////////////////////////////////////////////////////////

/**
 * matching helper class
 */
template<class T> class _pc_dist
{
  protected:

    T m_fScore;
    unsigned int m_iIndex;

  public:

    _pc_dist(T fScore, unsigned int iIndex)
    {
      m_fScore = fScore;
      m_iIndex = iIndex;
    };

    inline T getScore() const
    {
      return m_fScore;
    };
    inline unsigned int getIndex() const
    {
      return m_iIndex;
    };

    inline bool operator<(const _pc_dist &rR) const
    {
      return m_fScore < rR.m_fScore;
    };
};

/**
 * Destructor
 */
template<class T>
svt_point_cloud<T>::~svt_point_cloud()
{
};


/**
 * Set tree pruning parameter
 * \param fEpsilon epsilon
 */
template<class T>
inline void svt_point_cloud<T>::setEpsilon(Real64 fEpsilon)
{
  m_fEpsilon = fEpsilon;
};
/**
 * Get tree pruning parameter.
 * \param fEpsilon epsilon
 */
template<class T>
inline Real64 svt_point_cloud<T>::getEpsilon() const
{
  return m_fEpsilon;
};
/**
 * set tolorance distance for the anchor determination
 * \param fLambda lambda
 */
template<class T>
inline void svt_point_cloud<T>::setLambda(Real64 fLambda)
{
  m_fLambda = fLambda;
};
/**
 * get tolorance distance for the anchor determination
 * \return lambda
 */
template<class T>
inline Real64 svt_point_cloud<T>::getLambda() const
{
  return m_fLambda;
};
/**
 * set nearest neighbor matching zone size
 * \param fGamma gamma
 */
template<class T>
inline void svt_point_cloud<T>::setGamma(Real64 fGamma)
{
  m_fGamma = fGamma;
};
/**
 * get nearest neighbor matching zone size
 * \return gamma
 */
template<class T>
inline Real64 svt_point_cloud<T>::getGamma() const
{
  return m_fGamma;
};

/**
 * set the maximal size of the matching zone
 * \param iZoneSize
 */
template<class T>
inline void svt_point_cloud<T>::setZoneSize(unsigned int iZoneSize)
{
  m_iZoneSize = iZoneSize;
};
/**
 * get the maximal size of the matching zone
 * \return maximal size of matching zone
 */
template<class T>
inline unsigned int svt_point_cloud<T>::getZoneSize() const
{
  return m_iZoneSize;
};

/**
 * set the maximal number of wildcard matches
 * \param iMaxNoMatch maximal number of wildcard matches
 */
template<class T>
inline void svt_point_cloud<T>::setWildcards(unsigned int iMaxNoMatch)
{
  m_iMaxNoMatch = iMaxNoMatch;
};
/**
 * get the maximal number of wildcard matches
 * \return maximal number of wildcard matches
 */
template<class T>
inline unsigned int svt_point_cloud<T>::getWildcards() const
{
  return m_iMaxNoMatch;
};

/**
 * set the penalty for wildcard matches
 * \param fSkipPenalty penalty for a single wildcard
 */
template<class T>
inline void svt_point_cloud<T>::setSkipPenalty(Real64 fSkipPenalty)
{
  m_fSkipPenalty = fSkipPenalty;
}
/**
 * get the penalty for wildcard matches
 * \return penalty for a single wildcard
 */
template<class T>
inline Real64 svt_point_cloud<T>::getSkipPenalty() const
{
  return m_fSkipPenalty;
}

/**
 * if two solutions are very close only the one with the higher score is considered. The other solutions are removed.
 * \param minimal distance between solutions
 */
template<class T>
inline void svt_point_cloud<T>::setUnique(Real64 fUnique)
{
  m_fUnique = fUnique;
}
/**
 * if two solutions are very close only the one with the higher score is considered. The other solutions are removed.
 * \return minimal distance between solutions
 */
template<class T>
inline Real64 svt_point_cloud<T>::getUnique() const
{
  return m_fUnique;
};

/**
 * set the number of runs. Each time a different set of anchor points get selected from the probe structure, beginning with the three points furthest away from each other and the COA.
 * \param iRuns number of runs
 */
template<class T>
inline void svt_point_cloud<T>::setRuns(unsigned int iRuns)
{
  m_iRuns = iRuns;
};
/**
 * get the number of runs. Each time a different set of anchor points get selected from the probe structure, beginning with the three points furthest away from each other and the COA.
 * \return number of runs
 */
template<class T>
inline unsigned int svt_point_cloud<T>::getRuns() const
{
  return m_iRuns;
};

/**
 * Full or simple version of the matching algorithm to be performed.
 * \param bSimple if true only a reduced version is performed (faster but less accurate)
 */
template<class T>
inline void svt_point_cloud<T>::setSimple(bool bSimple)
{
  m_bSimple = bSimple;
};
/**
 * Full or simple version of the matching algorithm to be performed.
 * \return if true only a reduced version is performed (faster but less accurate)
 */
template<class T>
inline bool svt_point_cloud<T>::getSimple() const
{
  return m_bSimple;
};

/**
 * Set the next point selection scheme to COA
 */
template<class T>
inline void svt_point_cloud<T>::setNextPointCOA(bool bNextPointCOA)
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
template<class T>
void svt_point_cloud<T>::match(svt_point_cloud<T> &rPC, vector< svt_matchResult<Real64> > &rMatch)
{
  // calculate COA
  T oCOA = coa();
  unsigned int i, j;

  unsigned int iFstIndex = 0;
  unsigned int iSndIndex = 0;
  unsigned int iTrdIndex = 0;

  vector< int > aOldAnchors;
  vector< svt_matrix4<Real64> > aBestSolution;
  vector< int > aBestMatch;

  // some status output
  SVTLBBO << endl;
  SVTLBBO << "Point-cloud matching..." << endl;
  SVTLBBO << "  Subcomponent has " << size() << " points " << endl;
  SVTLBBO << "  Gets docked into a structure with " << rPC.size() << " points " << endl;

  // size of the two pointclouds
  unsigned int iNumA = size();

  // sum of all seed atoms evaluated
  unsigned int iSumSeeds = 0;

  // array with all the results from all runs
  vector< svt_matchResult<Real64> > aResults;

  for (unsigned int iRun = 0; iRun < m_iRuns; iRun++) {
    if (m_iRuns != 1) {
      SVTLBBO << endl;
      SVTLBBO << "  Run: " << iRun << endl;
    }

    // calc best anchor of model point cloud
    // 1st: which one is furthest away from the center
    Real64 fLength = 0.0;
    Real64 fDist = 0.0;
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
    Real64 fTmpDist = 0.0;
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
    Real64 fSumDist = 0.0;
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

    SVTLBBO << "  First Anchor: " << iFstIndex << endl;
    SVTLBBO << "  Second Anchor: " << iSndIndex << endl;
    SVTLBBO << "  Third Anchor: " << iTrdIndex << endl;

    // run the matching
    vector< svt_matrix4<Real64> > aMatrix;
    vector< svt_matchResult<Real64> > aMatch;
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
  vector< svt_matchResult< Real64 > > aFinal;
  vector< svt_vector4<Real64> > aComs;
  svt_vector4<Real64> oTransVec;
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

  SVTLBBO << endl;

  if (rMatch.size() > 0)
    SVTLBBO << "  Best Result: " << rMatch[0].getScore() << endl;
  else
    SVTLBBO << "  Best Result: No result found!" << endl;

  m_iNumSeeds = iSumSeeds;

  // end
  SVTLBBO << "  Done." << endl;
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
template<class T>
void svt_point_cloud<T>::match(svt_point_cloud<T> &rPC, vector< svt_matchResult<Real64> > &rMatch, unsigned int iAnchorA, unsigned int iAnchorB, unsigned int iAnchorC)
{
  unsigned int i, j;

  // start time
  //long iStartTime = svt_getToD();

  // empty the match
  rMatch.clear();

  // average nn distance
  Real64 fEpsilon = averageNNDistance(0.1);
  setEpsilon(fEpsilon + (1.0 * fEpsilon));

  // size of the two pointclouds
  unsigned int iNumB = rPC.size();

  // anchors
  unsigned int iFstIndex = iAnchorA;
  unsigned int iSndIndex = iAnchorB;
  unsigned int iTrdIndex = iAnchorC;

  // distances between anchors
  Real64 fDistAB = (*this)[iFstIndex].distance((*this)[iSndIndex]);
  Real64 fDistBC = (*this)[iSndIndex].distance((*this)[iTrdIndex]);
  Real64 fDistAC = (*this)[iFstIndex].distance((*this)[iTrdIndex]);
  vector< svt_vector4<int> > aSeeds;

  // now find appropriate seed anchors in the scene, based on the distances...
  for (i = 0; i < iNumB; i++)
    for (j = 0; j < iNumB; j++)
      if (i != j) {
        if (fabs(rPC[i].distance(rPC[j]) - fDistAB) < m_fLambda) {
          svt_vector4<int> oVec(i, j, -1);
          aSeeds.push_back(oVec);

        }
      }

  m_iNumSeeds = aSeeds.size();
  vector< svt_vector4<int> > aFinalSeeds;

  for (i = 0; i < m_iNumSeeds; i++)
    for (j = 0; j < iNumB; j++) {
      if (aSeeds[i].x() != (int)(j) && aSeeds[i].y() != (int)(j)) {
        if (
          (fabs(rPC[j].distance(rPC[ aSeeds[i].x() ]) - fDistAC) < m_fLambda &&
           fabs(rPC[j].distance(rPC[ aSeeds[i].y() ]) - fDistBC) < m_fLambda) ||
          (fabs(rPC[j].distance(rPC[ aSeeds[i].y() ]) - fDistAC) < m_fLambda &&
           fabs(rPC[j].distance(rPC[ aSeeds[i].x() ]) - fDistBC) < m_fLambda)
        ) {
          svt_vector4<int> oVec = aSeeds[i];
          oVec.z(j);
          aFinalSeeds.push_back(oVec);
        }
      }
    }

  m_iNumSeeds = aFinalSeeds.size();

  SVTLBBO << "  Number of found anchor matches: " << aFinalSeeds.size() << endl;

  // investigate the quality of the anchor matches
  svt_matrix4<Real64> oMatrix;
  vector<int> aMatches;
  vector< svt_matchResult<Real64> > aResults;
  Real64 fRMSD = 0.0;

  try {
    for (i = 0; i < m_iNumSeeds; i++) {

      // run A.1
      {
        vector<int> aModelAnchor;
        aModelAnchor.push_back(iFstIndex);
        aModelAnchor.push_back(iSndIndex);
        aModelAnchor.push_back(iTrdIndex);
        vector<int> aSceneAnchor;
        aSceneAnchor.push_back(aFinalSeeds[i].x());
        aSceneAnchor.push_back(aFinalSeeds[i].y());
        aSceneAnchor.push_back(aFinalSeeds[i].z());

        oMatrix = optimize(&aModelAnchor, &aSceneAnchor, *this, rPC, true);

        if (aModelAnchor.size() > 3) {
          fRMSD = calcRMSD(*this, rPC, oMatrix, &aModelAnchor, &aSceneAnchor);
          aResults.push_back(svt_matchResult<Real64>(fRMSD, oMatrix, aModelAnchor, aSceneAnchor));
        }
      }
      // run A.2
      /*      {
          vector<int> aModelAnchor;
          aModelAnchor.push_back( iFstIndex );
          aModelAnchor.push_back( iSndIndex );
          aModelAnchor.push_back( iTrdIndex );
          vector<int> aSceneAnchor;
          aSceneAnchor.push_back( aFinalSeeds[i].x() );
          aSceneAnchor.push_back( aFinalSeeds[i].z() );
          aSceneAnchor.push_back( aFinalSeeds[i].y() );

          oMatrix = optimize( &aModelAnchor, &aSceneAnchor, *this, rPC, true );

          if (aModelAnchor.size() > 3)
          {
              fRMSD = calcRMSD( *this, rPC, oMatrix, &aModelAnchor, &aSceneAnchor);
              aResults.push_back( svt_matchResult<Real64>( fRMSD, oMatrix, aModelAnchor, aSceneAnchor ) );
          }
            }
            // run A.3
            {
          vector<int> aModelAnchor;
          aModelAnchor.push_back( iFstIndex );
          aModelAnchor.push_back( iSndIndex );
          aModelAnchor.push_back( iTrdIndex );
          vector<int> aSceneAnchor;
          aSceneAnchor.push_back( aFinalSeeds[i].y() );
          aSceneAnchor.push_back( aFinalSeeds[i].x() );
          aSceneAnchor.push_back( aFinalSeeds[i].z() );

          oMatrix = optimize( &aModelAnchor, &aSceneAnchor, *this, rPC, true );

          if (aModelAnchor.size() > 3)
          {
              fRMSD = calcRMSD( *this, rPC, oMatrix, &aModelAnchor, &aSceneAnchor);
              aResults.push_back( svt_matchResult<Real64>( fRMSD, oMatrix, aModelAnchor, aSceneAnchor ) );
          }
            }
            // run A.4
            {
          vector<int> aModelAnchor;
          aModelAnchor.push_back( iFstIndex );
          aModelAnchor.push_back( iSndIndex );
          aModelAnchor.push_back( iTrdIndex );
          vector<int> aSceneAnchor;
          aSceneAnchor.push_back( aFinalSeeds[i].y() );
          aSceneAnchor.push_back( aFinalSeeds[i].z() );
          aSceneAnchor.push_back( aFinalSeeds[i].x() );

          oMatrix = optimize( &aModelAnchor, &aSceneAnchor, *this, rPC, true );

          if (aModelAnchor.size() > 3)
          {
              fRMSD = calcRMSD( *this, rPC, oMatrix, &aModelAnchor, &aSceneAnchor);
              aResults.push_back( svt_matchResult<Real64>( fRMSD, oMatrix, aModelAnchor, aSceneAnchor ) );
          }
            }
            // run A.5
            {
          vector<int> aModelAnchor;
          aModelAnchor.push_back( iFstIndex );
          aModelAnchor.push_back( iSndIndex );
          aModelAnchor.push_back( iTrdIndex );
          vector<int> aSceneAnchor;
          aSceneAnchor.push_back( aFinalSeeds[i].z() );
          aSceneAnchor.push_back( aFinalSeeds[i].x() );
          aSceneAnchor.push_back( aFinalSeeds[i].y() );

          oMatrix = optimize( &aModelAnchor, &aSceneAnchor, *this, rPC, true );

          if (aModelAnchor.size() > 3)
          {
              fRMSD = calcRMSD( *this, rPC, oMatrix, &aModelAnchor, &aSceneAnchor);
              aResults.push_back( svt_matchResult<Real64>( fRMSD, oMatrix, aModelAnchor, aSceneAnchor ) );
          }
            }
            // run A.6
            {
          vector<int> aModelAnchor;
          aModelAnchor.push_back( iFstIndex );
          aModelAnchor.push_back( iSndIndex );
          aModelAnchor.push_back( iTrdIndex );
          vector<int> aSceneAnchor;
          aSceneAnchor.push_back( aFinalSeeds[i].z() );
          aSceneAnchor.push_back( aFinalSeeds[i].y() );
          aSceneAnchor.push_back( aFinalSeeds[i].x() );

          oMatrix = optimize( &aModelAnchor, &aSceneAnchor, *this, rPC, true );

          if (aModelAnchor.size() > 3)
          {
              fRMSD = calcRMSD( *this, rPC, oMatrix, &aModelAnchor, &aSceneAnchor);
              aResults.push_back( svt_matchResult<Real64>( fRMSD, oMatrix, aModelAnchor, aSceneAnchor ) );
          }
                  }
                  */

    }

  } catch (int e) {
    return;
  }


  // FIXME: Debug code...
  /*cout << endl;
  for (unsigned int i=0; i<aBestModelM.size(); i++)
  printf( "%02i,", aBestModelM[i] );
  cout << endl;

  for (unsigned int i=0; i<aBestSceneM.size(); i++)
  if (aBestSceneM[i] != NOMATCH)
    printf( "%02i,", aBestSceneM[i] );
  else
    printf( " - " );
  cout << endl;
  getchar();
  */

  rMatch.clear();

  // elapsed time
  //long iETime = svt_getToD() - iStartTime;
  //unsigned int iHours = (unsigned int)((double)(iETime) / 3.6E9);
  //iETime -= (long)(iHours * 3.6E9);
  //unsigned int iMins = (unsigned int)((double)(iETime) / 6.0E7);
  //iETime -= (long)(iMins * 6.0E7);
  //unsigned int iSecs = (unsigned int)((double)(iETime) / 1.0E6);
  //SVTLBBO << "Elapsed time: " << iHours << ":" << iMins << ":" << iSecs << endl;

  if (aResults.size() == 0) {
    SVTLBBO << "  Error - there was no valid result in this run!" << endl;
    return;
  }

  // sort the resulting matrices
  sort(aResults.begin(), aResults.end());
  // copy the 10 best results into a matrix array
  for (i = 0; i < aResults.size() && i < 1000; i++)
    rMatch.push_back(aResults[i]);

  //SVTLBBO << "  Best Result: " << aResults[0].getScore() << endl;
};

/**
 * Get number of seed/anchor matches evaluated during last matching run.
 */
template<class T>
unsigned int svt_point_cloud<T>::getNumSeeds() const
{
  return m_iNumSeeds;
};

/**
 * Calculate the cross-correlation between this pointcloud and a volumetric object. The idea is to use the points in the pointcloud to sample the volume and to
 * form thereby a reduce cross-correlation coefficient equation.
 * See also:
 * Stefan Birmanns and Willy Wriggers. Interactive Fitting Augmented by Force-Feedback and Virtual Reality. J. Struct. Biol., 2003, Vol. 144, pp. 123-131.
 * \param rPointMat reference to svt_matrix4 object with the initial transformation matrix of the volumetric object
 * \param rVol reference to svt_volume object
 * \param rVolMat reference to svt_matrix4 object with the transformation matrix of the volumetric object
 * \return the cross-correlation
 */
template<class T>
Real64 svt_point_cloud<T>::calcCC(svt_matrix4<Real64> &rPointMat, svt_volume<Real64> &rVol, svt_matrix4<Real64> &rVolMat)
{
  Real64 fCC = 0.0;
  int iVoxelX, iVoxelY, iVoxelZ;
  Real64 fWidth = rVol.getWidth();
  Real64 fVoxelValue = 0.0;
  int iCount = 0;

  // inner = -1, if codebook vector is inner vector from laplace data
  //Real64 fInnerFac = 1; //(inner ? -1. : 1.);

  // calculate the center of the volume data
  T oOrigin;
  oOrigin.x((fWidth * rVol.getSizeX()) * 0.5f);
  oOrigin.y((fWidth * rVol.getSizeY()) * 0.5f);
  oOrigin.z((fWidth * rVol.getSizeZ()) * 0.5f);

  // transformation of point cloud into volumetric data
  svt_matrix4<Real64> oTransMatrix = rVolMat;
  oTransMatrix.invert();
  oTransMatrix *= rPointMat;

  // now calculate the force for each codebook vector
  for (unsigned int i = 0; i < size(); i++) {
    // get original position
    T oVec = (*this)[i];
    // move it into the situs coordinate system
    oVec = oTransMatrix * oVec;
    //oVec += oOrigin;

    // get the index of the voxel lying under the cv
    iVoxelX = (int)(floor(oVec.x() / fWidth));
    iVoxelY = (int)(floor(oVec.y() / fWidth));
    iVoxelZ = (int)(floor(oVec.z() / fWidth));

    // calculate the gradient at the voxel position - negativ weight for inner codebookvectors
    if (iVoxelX < 0 || iVoxelY < 0 || iVoxelZ < 0 || iVoxelX >= (int)(rVol.getSizeX()) || iVoxelY >= (int)(rVol.getSizeY()) || iVoxelZ >= (int)(rVol.getSizeZ()))
      continue;

    // position of voxel inside the volume
    Real64 fVoxelPosX = iVoxelX * fWidth;
    Real64 fVoxelPosY = iVoxelY * fWidth;
    Real64 fVoxelPosZ = iVoxelZ * fWidth;

    // position of cv inside the voxel
    Real64 fCVXV = (oVec.x() - fVoxelPosX) / fWidth;
    fCVXV += iVoxelX;
    Real64 fCVYV = (oVec.y() - fVoxelPosY) / fWidth;
    fCVYV += iVoxelY;
    Real64 fCVZV = (oVec.z() - fVoxelPosZ) / fWidth;
    fCVZV += iVoxelZ;

    // get interpolated value from volume
    fVoxelValue = rVol.getIntValue(fCVXV, fCVYV, fCVZV);
    if (fVoxelValue != 0.0)
      iCount++;
    fCC += fVoxelValue;
  }

  if (iCount != 0)
    fCC /= (Real64)(iCount);
  else
    fCC = 0;

  return fCC;
}

///////////////////////////////////////////////////////////////////////////////
// Private functions
///////////////////////////////////////////////////////////////////////////////

/**
 * sort eigenvectors according to their eigenvalues
 * \param pEigenvectors pointer to svt_matrix with the eigenvectors as column vectors
 * \param pEigenvalues pointer to svt_vector with the eigenvalues
 */
template<class T>
void svt_point_cloud<T>::eigensort(svt_matrix4<Real64> &rEigenvectors, T &rEigenvalues)
{
  Real64 fP;
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
template<class T>
svt_matrix4<Real64> svt_point_cloud<T>::kearsley(vector< int > &rMM, vector< int > &rSM, svt_point_cloud<T> &rModel, svt_point_cloud<T> &rScene)
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
  T oModelCOA;
  for (i = 0; i < iCount; i++)
    oModelCOA = oModelCOA + rModel[rM[i]];
  oModelCOA = oModelCOA / (double)(iCount);

  T oSceneCOA;
  for (i = 0; i < iCount; i++)
    oSceneCOA = oSceneCOA + rScene[rS[i]];
  oSceneCOA = oSceneCOA / (double)(iCount);

  // first step: setup Q matrix (see kearsley acta cryst a45)

  // distance plus and minus
  T oDP;
  T oDM;

  // setup matrix
  svt_matrix4<Real64> oQ;
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
  svt_matrix4<Real64> oEigenvectors;
  T oEigenvalues;
  oQ.jacobi(oEigenvectors, oEigenvalues);

  // third step: sort the eigenvectors according to their eigenvalues
  eigensort(oEigenvectors, oEigenvalues);

  // forth step: generate transformation matrix
  svt_matrix4<Real64> oFinal;
  svt_matrix4<Real64> oRot;

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
 * \param rModel vector of svt_vector4 with the model vectors
 * \param rScene vector of svt_vector4 with the scene vectors
 * \return matrix with optimal transformation of the model into the scene
 */
template<class T>
svt_matrix4<Real64> svt_point_cloud<T>::kearsley(svt_point_cloud<T> &rModel, svt_point_cloud<T> &rScene)
{
  unsigned int i;
  unsigned int iCount = rModel.size();

  // precomp.: calc COA for both model and scene and transformation
  T oModelCOA;
  for (i = 0; i < iCount; i++)
    oModelCOA = oModelCOA + rModel[i];
  oModelCOA = oModelCOA / (double)(iCount);

  T oSceneCOA;
  for (i = 0; i < iCount; i++)
    oSceneCOA = oSceneCOA + rScene[i];
  oSceneCOA = oSceneCOA / (double)(iCount);

  T oTransCOA;
  oTransCOA = oSceneCOA - oModelCOA;

  // first step: setup Q matrix (see kearsley acta cryst a45)

  // distance plus and minus
  T oDP;
  T oDM;

  // setup matrix
  svt_matrix4<Real64> oQ;
  // as svt_matrix creates an identity matrix, we have to set the diagonal elements to 0
  oQ[0][0] = 0.0;
  oQ[1][1] = 0.0;
  oQ[2][2] = 0.0;
  oQ[3][3] = 0.0;

  // now construct matrix
  for (i = 0; i < iCount; i++) {
    oDM = (rModel[i] - oModelCOA) - (rScene[i] - oSceneCOA);
    oDP = (rModel[i] - oModelCOA) + (rScene[i] - oSceneCOA);

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
  svt_matrix4<Real64> oEigenvectors;
  T oEigenvalues;
  oQ.jacobi(oEigenvectors, oEigenvalues);

  // third step: sort the eigenvectors according to their eigenvalues
  eigensort(oEigenvectors, oEigenvalues);

  // forth step: generate transformation matrix
  svt_matrix4<Real64> oFinal;
  svt_matrix4<Real64> oRot;

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
 * optimize_simple - subroutine for the match() procedure
 */
template<class T>
svt_matrix4<Real64> svt_point_cloud<T>::optimize_simple(vector<int> *pModelMatch, vector<int> *pSceneMatch, svt_point_cloud<T> &rModel, svt_point_cloud<T> &rScene, bool)
{
  SVTLBBO << "SIMPLE!" << endl;

  svt_matrix4<Real64> oMatrix = kearsley(*pModelMatch, *pSceneMatch, rModel, rScene);

  // transform...
  unsigned int i, j;
  unsigned int iIndex;
  svt_vector4<Real64> oVec;

  for (i = 0; i < rModel.size(); i++) {
    // is vector i unmatched?
    for (j = 0; j < pModelMatch->size(); j++)
      if ((*pModelMatch)[j] == (int)(i))
        break;

    // yes...
    if (j == pModelMatch->size()) {
      oVec = oMatrix * rModel[i];

      iIndex = rScene.nearestNeighbor(oVec);

      // is vector iIndex unmatched?
      for (j = 0; j < pSceneMatch->size(); j++)
        if ((*pSceneMatch)[j] == (int)(iIndex))
          break;

      // yes...
      if (j == pSceneMatch->size()) {
        pModelMatch->push_back(i);
        pSceneMatch->push_back(iIndex);
      }
    }
  }

  oMatrix = kearsley(*pModelMatch, *pSceneMatch, rModel, rScene);

  return oMatrix;
};

/**
 * optimize - subroutine for the match() procedure
 */
template<class T>
svt_matrix4<Real64> svt_point_cloud<T>::optimize(vector<int> *pModelMatch, vector<int> *pSceneMatch, svt_point_cloud<T> &rModel, svt_point_cloud<T> &rScene, bool bInit)
{
  static Real64 fBestRMSD = 1.0E10;
  static vector<int> aBestModelMatch;
  static vector<int> aBestSceneMatch;
  static svt_matrix4<Real64> oBestMatrix;
  static unsigned int iNumM;
  static unsigned int iNumS;
  static unsigned int iLevel;
  static vector< _pc_dist<Real64> > aDistToCOA;

  if (bInit == true) {
    fBestRMSD = 1.0E10;
    iNumM = rModel.size();
    iNumS = rScene.size();
    iLevel = 0;

    if (m_bNextPointCOA == true) {
      // calculate distances of all model vectors to COA
      T oCOA = rModel.coa();
      aDistToCOA.clear();
      for (unsigned int i = 0; i < iNumM; i++)
        aDistToCOA.push_back(_pc_dist<Real64>(rModel[i].distance(oCOA), i));

      sort(aDistToCOA.begin(), aDistToCOA.end());
    }

  } else
    iLevel++;

  // FIXME: Debug output
  //SVTLBBO << "DEBUG: iLevel:" << iLevel << endl;
  //getchar();

  // find an unmatched model vector...
  T oVec;
  unsigned int iIndex = 0, j = 0, k = 0;

  if (pModelMatch->size() < rModel.size()) {
    // get optimal transformation for the old match
    svt_matrix4<Real64> oMatrix = kearsley(*pModelMatch, *pSceneMatch, rModel, rScene);

    // look for an unmatched point by moving towards the COA...
    if (m_bNextPointCOA) {
      unsigned int i;
      for (i = aDistToCOA.size(); i > 0; i--)
        if (find(pModelMatch->begin(), pModelMatch->end(), (int)(aDistToCOA[i - 1].getIndex())) == pModelMatch->end())
          break;

      iIndex = aDistToCOA[i - 1].getIndex();

    } else {

      // look for an unmatched point that is closest to a scenepoint
      Real64 fMinDist = 1.0E10;
      Real64 fDist = 0.0;
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

    // FIXME: Debug output
    /*cout << endl;
    for (unsigned int i=0; i<pModelMatch->size(); i++)
        printf( "%03i,", (*pModelMatch)[i] );
          printf("%03i", iIndex );
    cout << endl;

    for (unsigned int i=0; i<pSceneMatch->size(); i++)
        if ((*pSceneMatch)[i] != NOMATCH)
      printf( "%03i,", (*pSceneMatch)[i] );
        else
      printf( "--- " );
    cout << endl;

          SVTLBBO << oMatrix << endl;

    svt_point_cloud< svt_vector4<Real64> > oWrite;
    for(unsigned int iI=0; iI<rModel.size(); iI++)
    {
              oVec = oMatrix * rModel[iI];
        oWrite.addPoint( oVec );
    }
          oWrite.writePDB( "optimize.pdb" );
    getchar();*/

    // transform
    oVec = oMatrix * rModel[iIndex];

    // now lets find the nearest neighbors in the scene...
    vector< _pc_dist<Real64> > aNeighbors;
    Real64 fDist = 0.0;
    for (j = 0; j < iNumS; j++) {
      // is this scene point unmatched?
      for (k = 0; k < pSceneMatch->size(); k++) {
        //cout << "pSceneMatch[" << k << "]:" << (*pSceneMatch)[k] << endl;

        if ((*pSceneMatch)[k] == (int)j)
          break;
      }
      // yes...
      if (k >= pSceneMatch->size()) {
        fDist = rScene[j].distance(oVec);

        //cout << j << " - fDist: " << fDist << " m_fGamma: " << m_fGamma << endl;

        if (fDist < m_fGamma)
          aNeighbors.push_back(_pc_dist<Real64>(fDist, j));
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
          // FIXME: Debug...
          //SVTLBBO << "DEBUG: Add Match..." << endl;

          vector<int> *pTmpMatchModel = new vector<int>;
          *pTmpMatchModel = *pModelMatch;
          vector<int> *pTmpMatchScene = new vector<int>;
          *pTmpMatchScene = *pSceneMatch;
          pTmpMatchModel->push_back(iIndex);
          pTmpMatchScene->push_back(aNeighbors[j].getIndex());
          optimize(pTmpMatchModel, pTmpMatchScene, rModel, rScene, false);
          delete pTmpMatchModel;
          delete pTmpMatchScene;
          iCount++;
        }
      }
      // no, no potential matches, so we can try a NOMATCH, but not too many...
    }
    if (iNumN == 0 || aNeighbors[0].getScore() > 3.0) {
      // as we sort by nearest neighbor distance, it is strange that we did not find any more neighbors.
      /*if (rModel.size() - pModelMatch->size() < m_iMaxNoMatch)
      {
      svt_matrix4<Real64> oMatrix = kearsley( *pModelMatch, *pSceneMatch, rModel, rScene );
      Real64 fRMSD = calcRMSD( rModel, rScene, oMatrix, pModelMatch, pSceneMatch );

                // fill the rest with NOMATCH
      for(i=0; i<rModel.size(); i++)
      {
        // which model point is unmatched?
        if (find(pModelMatch->begin(), pModelMatch->end(), (int)(i) ) == pModelMatch->end())
        {
      pModelMatch->push_back( iIndex );
      pSceneMatch->push_back( NOMATCH );
        }
      }

      // do we have a winner?
      if (fRMSD < fBestRMSD)
      {
        fBestRMSD = fRMSD;
        oBestMatrix = oMatrix;
        aBestModelMatch = *pModelMatch;
        aBestSceneMatch = *pSceneMatch;
      }
      };*/

      unsigned int iNoMatch = 0;
      for (j = 0; j < pSceneMatch->size(); j++)
        if ((*pSceneMatch)[j] == NOMATCH)
          iNoMatch++;

      // if there are no more than iMaxNoMatch blind matches, try them
      if (iNoMatch < m_iMaxNoMatch) {
        // FIXME: Debug...
        //SVTLBBO << "DEBUG: Add NoMatch..." << endl;

        vector<int> *pTmpMatchModel = new vector<int>;
        *pTmpMatchModel = *pModelMatch;
        vector<int> *pTmpMatchScene = new vector<int>;
        *pTmpMatchScene = *pSceneMatch;
        pTmpMatchModel->push_back(iIndex);
        pTmpMatchScene->push_back(NOMATCH);
        optimize(pTmpMatchModel, pTmpMatchScene, rModel, rScene, false);
        delete pTmpMatchModel;
        delete pTmpMatchScene;

        /*} else {
        // FIXME: Debug...
        SVTLBBO << "DEBUG: No more wildcards available" << endl;
        }*/
      }
    }

    // no unmatched vector left, so lets evaluate this solution...
  } else {

    if (pModelMatch->size() == rModel.size()) {
      //SVTLBBO << "DEBUG: OK, complete!" << endl;

      svt_matrix4<Real64> oMatrix = kearsley(*pModelMatch, *pSceneMatch, rModel, rScene);
      Real64 fRMSD = 0.0;
      if (m_iMaxNoMatch == 0) {
        fRMSD = calcRMSD(rModel, rScene, oMatrix, pModelMatch, pSceneMatch);

      } else {

        Real64 fPercent = (Real64)(rModel.size()) / (Real64)(m_iMaxNoMatch);
        if (fPercent > 1.0)
          fPercent = 1.0;
        fRMSD = rModel.rmsd_NN_Outliers(rScene, fPercent);
      }

      //SVTLBBO << "DEBUG: fRMSD: " << fRMSD << endl;

      // do we have a winner?
      if (fRMSD < fBestRMSD) {
        fBestRMSD = fRMSD;
        oBestMatrix = oMatrix;
        aBestModelMatch = *pModelMatch;
        aBestSceneMatch = *pSceneMatch;

        unsigned int iNoMatch = 0;
        for (j = 0; j < pSceneMatch->size(); j++)
          if ((*pSceneMatch)[j] == NOMATCH)
            iNoMatch++;
      }
    }

  }

  // OK, do we exit this now? Then copy the best match...
  if (bInit == true) {
    if (aBestModelMatch.size() == rModel.size()) {
      //SVTLBBO << "DEBUG: End of recursion, plus matching complete!" << endl;

      /*svt_matrix4<Real64> oMatrix = kearsley( *pModelMatch, *pSceneMatch, rModel, rScene );
      Real64 fRMSD = calcRMSD( rModel, rScene, oMatrix, pModelMatch, pSceneMatch );

      SVTLBBO << "DEBUG: fRMSD: " << fRMSD << endl;
      SVTLBBO << "DEBUG: fBestRMSD: " << fBestRMSD << endl;

      // do we have a winner?
      if (fRMSD < fBestRMSD)
      {
      SVTLBBO << "DEBUG: Winnerrrrr!" << endl;

      fBestRMSD = fRMSD;
      oBestMatrix = oMatrix;
      aBestModelMatch = *pModelMatch;
      aBestSceneMatch = *pSceneMatch;

      unsigned int iNoMatch = 0;
      for(j=0; j<pSceneMatch->size(); j++)
        if ((*pSceneMatch)[j] == NOMATCH)
      iNoMatch++;

      if (iNoMatch != 0)
        SVTLBBO << endl << "Wildcards used: " << iNoMatch << endl;

      }*/
      (*pModelMatch) = aBestModelMatch;
      (*pSceneMatch) = aBestSceneMatch;
    }
  }

  return oBestMatrix;
};

/**
 * internal convenience function - no real difference with hausdorff - only here we don't need to create/copy a svt_point_cloud
 */
template<class T>
Real64 svt_point_cloud<T>::calcHausdorff(svt_point_cloud<T> &rModel, svt_point_cloud<T> &rScene, svt_matrix4<Real64> &rMatrix)
{
  unsigned int iNumA = rModel.size();
  unsigned int iNumB = rScene.size();
  unsigned int i, j;

  Real64 fMinDist;
  Real64 fHD = 0.0;
  T oVec;
  int iIndexA = 0, iIndexB = 0;
  int iMinIndex = 0;

  for (i = 0; i < iNumA; i++) {
    oVec = rMatrix * rModel[i];
    fMinDist = 1.0E10;
    iMinIndex = -1;

    for (j = 0; j < iNumB; j++) {
      if (oVec.distance(rScene[j]) < fMinDist) {
        fMinDist = oVec.distance(rScene[j]);
        iMinIndex = j;
      }
    }

    if (fMinDist > fHD) {
      fHD = fMinDist;
      iIndexA = i;
      iIndexB = iMinIndex;
    }
  }

  return fHD;
};

/**
 * internal convenience function - no real difference with rmsd - only here we don't need to create/copy a svt_point_cloud
 */
template<class T>
Real64 svt_point_cloud<T>::calcRMSD(svt_point_cloud<T> &rModel, svt_matrix4<Real64> oMatA, svt_matrix4<Real64> oMatB)
{
  unsigned int iNum = rModel.size();
  unsigned int i, j;

  Real64 fMinDist;
  Real64 fDist;
  Real64 fRMSD = 0.0;
  T oVecA;
  T oVecB;

  for (i = 0; i < iNum; i++) {
    oVecA = oMatA * rModel[i];
    fMinDist = 1.0E10;

    for (j = 0; j < iNum; j++) {
      oVecB = oMatB * rModel[j];
      fDist = oVecA.distanceSq(oVecB);

      if (fDist < fMinDist)
        fMinDist = fDist;
    }

    fRMSD += fMinDist;
  }

  fRMSD /= iNum;
  fRMSD = sqrt(fRMSD);

  return fRMSD;
};

/**
 * internal convenience function - no real difference with rmsd - only here we don't need to create/copy a svt_point_cloud
 */
template<class T>
Real64 svt_point_cloud<T>::calcRMSD(svt_point_cloud<T> &rModel, svt_point_cloud<T> &rScene, svt_matrix4<Real64> &rMatrix)
{
  unsigned int iNumA = rModel.size();
  unsigned int iNumB = rScene.size();
  unsigned int i, j;

  Real64 fMinDist;
  Real64 fRMSD = 0.0;
  T oVec;

  for (i = 0; i < iNumA; i++) {
    oVec = rMatrix * rModel[i];
    fMinDist = 1.0E10;

    for (j = 0; j < iNumB; j++) {
      if (oVec.distanceSq(rScene[j]) < fMinDist)
        fMinDist = oVec.distanceSq(rScene[j]);
    }

    fRMSD += fMinDist;
  }

  fRMSD /= iNumA;
  fRMSD = sqrt(fRMSD);

  return fRMSD;
};

/**
 * internal convenience function - no real difference with rmsd - only here we don't need to create/copy a svt_point_cloud
 * plus here we take the matching into account! If NOMATCH this point will not be used.
 */
template<class T>
Real64 svt_point_cloud<T>::calcRMSD(svt_point_cloud<T> &rModel, svt_point_cloud<T> &rScene, svt_matrix4<Real64> &rMatrix, vector<int> *pSceneMatch)
{
  unsigned int iNumA = rModel.size();
  unsigned int iNumB = rScene.size();
  unsigned int i, j;

  if (iNumA != (*pSceneMatch).size()) {
    SVTLBBO << "DEBUG: calcRMSD: pSceneMatch does not contain enough point matches!!!" << endl;
    exit(1);
  }

  Real64 fMinDist;
  Real64 fRMSD = 0.0;
  int iCount = 0;
  T oVec;

  for (i = 0; i < iNumA; i++) {
    if ((*pSceneMatch)[i] != NOMATCH) {

      oVec = rMatrix * rModel[i];
      fMinDist = 1.0E10;

      for (j = 0; j < iNumB; j++) {
        if (oVec.distanceSq(rScene[j]) < fMinDist)
          fMinDist = oVec.distanceSq(rScene[j]);
      }

      fRMSD += fMinDist;
      iCount++;

    } else

      fRMSD += m_fSkipPenalty;
  }

  fRMSD /= iCount;
  fRMSD = sqrt(fRMSD);

  return fRMSD;
};

/**
 * internal convenience function - no real difference with rmsd - only here we don't need to create/copy a svt_point_cloud
 * plus here we take the matching into account! If NOMATCH this point will not be used.
 */
template<class T>
Real64 svt_point_cloud<T>::calcRMSD(svt_point_cloud<T> &rModel, svt_point_cloud<T> &rScene, svt_matrix4<Real64> &rMatrix, vector<int> *pModelMatch, vector<int> *pSceneMatch)
{

  /**svt_point_cloud<T> oTmp = rMatrix * rModel;

  return rScene.rmsd_NN( oTmp );*/

  unsigned int iNum = pModelMatch->size();
  unsigned int i;

  if (iNum != (*pSceneMatch).size()) {
    SVTLBBO << "DEBUG: calcRMSD: pSceneMatch does not contain enough point matches!!!" << endl;
    exit(1);
  }

  Real64 fRMSD = 0.0;
  int iCount = 0;
  T oVec;

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



/***************************************************************************
                          svt_point_cloud_basic
                          ---------------------
    begin                : 02/10/2006
    author               : Stefan Birmanns
    email                : Stefan.Birmanns@uth.tmc.edu
 ***************************************************************************/

///////////////////////////////////////////////////////////////////////////////
// Implementation of basic functionality of svt_point_cloud
///////////////////////////////////////////////////////////////////////////////
// DO NOT INCLUDE DIRECTLY! File gets included by other header file!
///////////////////////////////////////////////////////////////////////////////

/**
 * Get all points in point cloud in an stl vector
 * \return reference to vector of svt_vector4 objects
 */
template<class T>
inline vector< T > &svt_point_cloud<T>::getPoints()
{
  return m_aPoints[m_iTimestep];
};

/**
 * Add a point to point cloud.
 * \param rVec svt_vector4 object
 */
template<class T>
inline void svt_point_cloud<T>::addPoint(T &rVec)
{
  m_aPoints[m_iTimestep].push_back(rVec);
};

/**
 * Add a point to point cloud.
 * \param fX x coord
 * \param fY y coord
 * \param fZ z coord
 */
template<class T>
inline void svt_point_cloud<T>::addPoint(Real64 fX, Real64 fY, Real64 fZ)
{
  T oVec;
  oVec.x(fX);
  oVec.y(fY);
  oVec.z(fZ);
  m_aPoints[m_iTimestep].push_back(oVec);
};

/**
 * Get a point out of point cloud.
 * \param iIndex index of point
 * \return reference to svt_vector4 object
 */
template<class T>
inline T &svt_point_cloud<T>::getPoint(unsigned int iIndex)
{
  return m_aPoints[m_iTimestep][iIndex];
};

/**
 * Replace point in point cloud.
 * \param iIndex index of point
 * \param rVec svt_vector4 object
 */
template<class T>
inline void svt_point_cloud<T>::setPoint(unsigned int iIndex, T &rVec)
{
  m_aPoints[m_iTimestep][iIndex] = rVec;
};

/**
* Delete point in point cloud.
* \param iIndex index of point
*/
template<class T>
inline void svt_point_cloud<T>::delPoint(unsigned int iIndex)
{
  m_aPoints[m_iTimestep].erase(m_aPoints[m_iTimestep].begin() + iIndex);
};


/**
 * Delete all points in point cloud
 */
template<class T>
inline void svt_point_cloud<T>::delAllPoints()
{
  m_aPoints[m_iTimestep].clear();
};

/**
 * Size of point cloud
 * \return size of pc
 */
template<class T>
inline unsigned int svt_point_cloud<T>::size() const
{
  return m_aPoints[m_iTimestep].size();
};

/**
 * Dereference operator (not range checked!).
 * \param iIndex index of point in point cloud
 */
template<class T>
inline T &svt_point_cloud<T>::operator[](unsigned int iIndex)
{
  return m_aPoints[m_iTimestep][iIndex];
};

/**
 * Calculate the geometric center of the point cloud.
 * \return svt_vector4 with the geometric center
 */
template<class T>
T svt_point_cloud<T>::geometricCenter()
{
  unsigned int iNum = size();
  unsigned int i = 0;
  T oMinVec(1.0E10,  1.0E10,  1.0E10);
  T oMaxVec(-1.0E10, -1.0E10, -1.0E10);

  for (i = 0; i < iNum; i++) {
    if ((*this)[i].x() < oMinVec.x())
      oMinVec.x((*this)[i].x());
    if ((*this)[i].y() < oMinVec.y())
      oMinVec.y((*this)[i].y());
    if ((*this)[i].z() < oMinVec.z())
      oMinVec.z((*this)[i].z());

    if ((*this)[i].x() > oMaxVec.x())
      oMaxVec.x((*this)[i].x());
    if ((*this)[i].y() > oMaxVec.y())
      oMaxVec.y((*this)[i].y());
    if ((*this)[i].z() > oMaxVec.z())
      oMaxVec.z((*this)[i].z());
  }

  T oTrans = (oMaxVec - oMinVec) * 0.5;
  oTrans = oTrans + oMinVec;

  return oTrans;
};

/**
 * Calculate center of atoms (COM with all masses = 1).
 * \return svt_vector4 with the COA
 */
template<class T>
T svt_point_cloud<T>::coa()
{
  unsigned int iNum = size();
  unsigned int i = 0;
  T oCOAVec(0.0,  0.0,  0.0);

  for (i = 0; i < iNum; i++)
    oCOAVec =  oCOAVec + (*this)[i];

  oCOAVec /= (double) iNum;

  return oCOAVec;
};

/**
 * calculate RMSD between this and another PC. The PCs must be matched already! To get the minimal RMSD please use align() first!
 * \param rPC second point cloud
 */
template<class T>
Real64 svt_point_cloud<T>::rmsd(svt_point_cloud<T> &rPC)
{
  if (size() != rPC.size()) {
    SVTLBBO << "Cannot calculate rmsd - the point_clouds have different size!" << endl;
    exit(1);
  }

  unsigned int iNumA = size();
  unsigned int i;

  Real64 fRMSD = 0.0;
  Real64 fDist;

  for (i = 0; i < iNumA; i++) {
    fDist = (*this)[i].distanceSq(rPC[i]);

    fRMSD += fDist;
  }

  fRMSD *= 1.0 / (Real64)(iNumA);
  fRMSD = sqrt(fRMSD);

  return fRMSD;
};

/**
 * calculate RMSD between this and another PC. The points are matched using the nearest neighbor relationship.
 * \param rPC second point cloud
 */
template<class T>
Real64 svt_point_cloud<T>::rmsd_NN(svt_point_cloud<T> &rPC)
{
  unsigned int iNumA = size();
  unsigned int i;

  Real64 fRMSD = 0.0;
  Real64 fDist;
  unsigned int iNeighbor;

  for (i = 0; i < iNumA; i++) {
    iNeighbor = rPC.nearestNeighbor((*this)[i]);
    fDist = (*this)[i].distanceSq(rPC[iNeighbor]);

    fRMSD += fDist;
  }

  fRMSD *= 1.0 / (Real64)(iNumA);
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
template<class T>
Real64 svt_point_cloud<T>::rmsd_NN_Outliers(svt_point_cloud<T> &rPC, Real64 fPercent, svt_matrix4<Real64> *pMatrix)
{
  unsigned int iNumA = size();
  unsigned int i;

  Real64 fRMSD = 0.0;
  Real64 fDist;
  unsigned int iNeighbor;
  vector< svt_neighbor<Real64> > *pNeighbors = new vector< svt_neighbor<Real64> >;
  svt_point_cloud<T> oTmp;

  if (pMatrix != NULL)
    oTmp = (*pMatrix) * (*this);
  else
    oTmp = (*this);

  for (i = 0; i < iNumA; i++) {
    iNeighbor = rPC.nearestNeighbor(oTmp[i]);
    fDist = oTmp[i].distanceSq(rPC[iNeighbor]);

    pNeighbors->push_back(svt_neighbor<Real64>(fDist, i, iNeighbor));
  }

  sort(pNeighbors->begin(), pNeighbors->end());

  unsigned int iPercent = (unsigned int)((Real64)(iNumA) * fPercent);

  vector<int> aModelMatch;
  vector<int> aSceneMatch;

  for (i = 0; i < iPercent; i++) {
    fRMSD += (*pNeighbors)[i].getScore();

    aModelMatch.push_back((*pNeighbors)[i].getIndexA());
    aSceneMatch.push_back((*pNeighbors)[i].getIndexB());
  }

  fRMSD /= (Real64)(iPercent);
  fRMSD = sqrt(fRMSD);

  if (pMatrix != NULL) {
    svt_matrix4<Real64> oOrig = (*pMatrix);
    svt_matrix4<Real64> oMat;
    oMat = oTmp.kearsley(aModelMatch, aSceneMatch, oTmp, rPC);
    (*pMatrix) = oOrig * oMat;
  }

  delete pNeighbors;

  return fRMSD;
}


/**
 * calculate hausdorff distance between this and another PC. The PCs must be matched and aligned already!
 * \param rPC second point cloud
 */
template<class T>
Real64 svt_point_cloud<T>::hausdorff(svt_point_cloud<T> &rPC)
{
  unsigned int iNumA = size();
  unsigned int i;

  Real64 fDist;
  Real64 fHD = 0.0;
  unsigned int iNeighbor;

  for (i = 0; i < iNumA; i++) {
    iNeighbor = rPC.nearestNeighbor((*this)[i]);
    fDist = (*this)[i].distance(rPC[iNeighbor]);

    if (fDist > fHD)
      fHD = fDist;
  }

  return fHD;
};

/**
 * find the nearest neighbor to a query point in the point cloud
 * \param rVec reference to svt_vector4 object - the query point
 * \return index to nearest point in point cloud
 */
template<class T>
unsigned int svt_point_cloud<T>::nearestNeighbor(T &rVec)
{
  unsigned int iNumA = size();

  Real64 fMinDist = 1.0E10;
  unsigned int iMinDist = 0;
  Real64 fDist;

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
 * calculate the average nearest neighbor distance
 * in order to reduce the complexity a random test is done and once the average stabilizes, the search is stopped.
 * \param fPrecision if average does not change by more than fPrecision between two iterations the calculation is stopped
 * \return average nearest neighbor distance
 */
template<class T>
Real64 svt_point_cloud<T>::averageNNDistance(Real64 fPrecision)
{
  unsigned int iNum = size();
  Real64 fNum = (Real64)(iNum);
  unsigned int iTest = 0;
  Real64 fCount = 0.0;
  Real64 fAvg = 0.0;
  Real64 fAvgDiv = 0.0;
  Real64 fAvgDivOld = 1.0E10;
  Real64 fMinDist;
  unsigned int iMinDist;
  Real64 fDist;
  unsigned int iOldTest = 0;

  if (fNum > 0) {
    // loop until desired precision is reached...
    while (fabs(fAvgDiv - fAvgDivOld) > fPrecision && fCount < 1000000.0) {
      while (iTest == iOldTest)
        iTest = (unsigned int)(svt_genrand() * fNum);
      iOldTest = iTest;

      // determine nearest neighbor
      fMinDist = 1.0E10;
      iMinDist = 0;

      for (unsigned int i = 0; i < iNum; i++) {
        fDist = (*this)[i].distance((*this)[iTest]);

        if (fDist < fMinDist && fDist != 0.0) {
          iMinDist = i;
          fMinDist = fDist;
        }
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
 * calculate the maximum nearest neighbor distance
 * \return maximum nearest neighbor distance
 */
template<class T>
Real64 svt_point_cloud<T>::maxNNDistance()
{
  unsigned int iNum = size();
  Real64 fMinDist;
  Real64 fMaxDist = 0.0;
  Real64 fDist;

  // loop until desired precision is reached...
  for (unsigned int i = 0; i < iNum; i++) {
    // determine nearest neighbor
    fMinDist = 1.0E10;

    for (unsigned int j = 0; j < iNum; j++) {
      fDist = (*this)[i].distance((*this)[j]);

      if (fDist < fMinDist && fDist != 0.0)
        fMinDist = fDist;
    }

    // is this the maximum?
    if (fMinDist > fMaxDist)
      fMaxDist = fMinDist;
  }

  return fMaxDist;
};

/**
 * Get the minimal coordinates of the point cloud - it will return a vector that has in each dimension the information about the minimal
 * coordinate it has found in the cloud.
 */
template<class T>
T svt_point_cloud<T>::getMinCoord()
{
  Real64 fMinX = 1.0E10;
  Real64 fMinY = 1.0E10;
  Real64 fMinZ = 1.0E10;

  unsigned int i;

  for (i = 0; i < size(); i++) {
    if (getPoint(i).x() < fMinX)
      fMinX = getPoint(i).x();
    if (getPoint(i).y() < fMinY)
      fMinY = getPoint(i).y();
    if (getPoint(i).z() < fMinZ)
      fMinZ = getPoint(i).z();
  }

  T oVec;
  oVec.x(fMinX);
  oVec.y(fMinY);
  oVec.z(fMinZ);

  return oVec;
};

/**
 * Get the maximal coordinates of the point cloud - it will return a vector that has in each dimension the information about the maximal
 * coordinate it has found in the cloud.
 */
template<class T>
T svt_point_cloud<T>::getMaxCoord()
{
  Real64 fMaxX = -1.0E10;
  Real64 fMaxY = -1.0E10;
  Real64 fMaxZ = -1.0E10;

  unsigned int i;

  for (i = 0; i < size(); i++) {
    if (getPoint(i).x() > fMaxX)
      fMaxX = getPoint(i).x();
    if (getPoint(i).y() > fMaxY)
      fMaxY = getPoint(i).y();
    if (getPoint(i).z() > fMaxZ)
      fMaxZ = getPoint(i).z();
  }

  T oVec;
  oVec.x(fMaxX);
  oVec.y(fMaxY);
  oVec.z(fMaxZ);

  return oVec;
};

/**
 * product of matrix and point cloud
 */
template<class T>
inline svt_point_cloud<T> operator*(const svt_matrix4<Real64> &rM, svt_point_cloud<T> &rPC)
{
  svt_point_cloud<T> oPC;
  T oP;

  for (unsigned int i = 0; i < rPC.size(); i++) {
    oP = rM * rPC[i];
    oPC.addPoint(oP);
  }

  return oPC;
};

/**
 * sample the object randomly and return a vector that refrects the probability distribution of the object
 */
template<class T>
T svt_point_cloud<T>::sample()
{
  unsigned int iIndex = (unsigned int)(svt_genrand() * (double)(size()));

  return this->getPoint(iIndex);
};

/**
 * Write pdb file.
 * \param pFilename pointer to array of char with the filename
 * \param bAppend if true, the pdb structure as append at the end of an existing structure file
 */
template<class T>
void svt_point_cloud<T>::writePDB(const char *pFilename, bool bAppend)
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

  int iTS;
  for (iTS = 0; iTS < this->getMaxTimestep(); iTS++) {
    if (this->getMaxTimestep() > 0)
      SVTLBBO << "  Writing frame " << iTS << endl;

    this->setTimestep(iTS);

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

    if (this->getMaxTimestep() > 0)
      fprintf(pFile, "END\n");
  }

  if (this->getMaxTimestep() > 0)
    this->setTimestep(0);

  fclose(pFile);

};

/**
 * Write CSV file - comma separated values x,y,z. It will write only the coordinates of the current frame.
 * \param pFilename pointer to array of char with the filename
 */
template<class T>
void svt_point_cloud<T>::writeCSV(const char *pFilename)
{
  unsigned int i;
  unsigned int iNum = size();
  FILE *pFile = fopen(pFilename, "w");

  for (i = 0; i < iNum; i++)
    fprintf(
      pFile, "%8.3f,%8.3f,%8.3f\n",
      this->getPoint(i).x(),
      this->getPoint(i).y(),
      this->getPoint(i).z()
    );

  fclose(pFile);
};

/**
 * Load a pdb file.
 * \param pointer to array of char with the filename
 */
template<class T>
void svt_point_cloud<T>::loadPDB(const char *pFilename)
{
  // open file
  svt_column_reader oReader(pFilename);
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
      double fX = oReader.extractReal32(30, 37);
      double fY = oReader.extractReal32(38, 45);
      double fZ = oReader.extractReal32(46, 53);

      this->addPoint(fX, fY, fZ);
    }
};

/**
 * set the current timestep
 * \param iTimestep the new timestep
 */
template<class T>
void svt_point_cloud<T>::setTimestep(unsigned int iTimestep)
{
  if (iTimestep < m_aPoints.size())
    m_iTimestep = iTimestep;
};

/**
 * get the current timestep
 * \return the timestep
 */
template<class T>
int svt_point_cloud<T>::getTimestep()
{
  return m_iTimestep;
};

/**
 * get the maximum timestep
 * \return the maximum timestep
 */
template<class T>
int svt_point_cloud<T>::getMaxTimestep()
{
  return m_aPoints.size();
};

/**
 * Add a new timestep
 */
template<class T>
void svt_point_cloud<T>::addTimestep()
{
  m_aPoints.push_back(m_aPoints[ m_aPoints.size() - 1 ]);
  m_iTimestep = m_aPoints.size() - 1;
};

/**
 * Compute the distance matrix between all points in the pointcloud
 * \return svt_matrix object with the distances
 */
template<class T>
svt_matrix<Real64> svt_point_cloud<T>::getDistanceMat()
{
  unsigned int iSize = size();
  svt_matrix<Real64> oDist(iSize, iSize);
  Real64 fDist;
  for (unsigned int iIndex1 = 0; iIndex1 < iSize; iIndex1++) {
    for (unsigned int iIndex2 = 0; iIndex2 < iSize; iIndex2++) {
      fDist = (*this)[iIndex1].distance((*this)[iIndex2]);
      oDist[iIndex1][iIndex2] = fDist;
    }
  }

  return oDist;
}


///////////////////////////////////////////////////////////////////////////////
// SVT_POINT_CLOUD_PDB
///////////////////////////////////////////////////////////////////////////////
/***************************************************************************
                          svt_point_cloud_pdb
                          -------------------
    begin                : 02/10/2006
    author               : Stefan Birmanns
    email                : Stefan.Birmanns@uth.tmc.edu
 ***************************************************************************/
#include <stdlib.h>
#include <limits.h>
//added libraries for fetchPDB
#ifndef WIN32
#include <sys/types.h>
#include <sys/socket.h>
#include <netdb.h>
#endif

#ifdef WIN32
#include <string.h>
#include <windows.h>
#include <winsock.h>
#include <winsock2.h>
#include <ws2tcpip.h>
#define strcasecmp _stricmp
#endif

#ifdef MACOSX
#include <CFURL.h>
#endif

#define MAXLINE 2048

///////////////////////////////////////////////////////////////////////////////
// svt_point_cloud_pdb
///////////////////////////////////////////////////////////////////////////////

enum SymmetryType {
  SYMMETRY_C,
  SYMMETRY_D,
  SYMMETRY_H
};

typedef enum {
  SSE_NOTAVAILABLE,
  SSE_PDB, // from the pdb
  SSE_STRIDE, //stride
  SSE_OTHER
} svt_secStructSource;

enum Select {
  ALL,
  BACKBONE,
  TRACE
};

/**
 * class containing secondary structure information as read from Pdb
 */
class svt_sse_pdb
{
  public:
    char m_aType[7];
    int m_iNum;

    char m_aID[4];

    char m_aInitialResname[5];
    char m_cInitialResChainID;
    int m_iInitialResSeq;
    char m_cInitialICode;

    char m_aTerminalResname[5];
    char m_cTerminalResChainID;
    int m_iTerminalResSeq;
    char m_cTerminalICode;

    //helix information
    int m_iClass;
    char m_aComment[31];
    int m_iLen;

    // sheet information
    int m_iNumStrands;
    int m_iSense;

    char m_aCurAtom[5];
    char m_aCurResname[5];
    char m_cCurResChainID;
    char m_aCurResSeq[5];
    char m_cCurICode;

    char m_aPrevAtom[5];
    char m_aPrevResname[5];
    char m_cPrevResChainID;
    char m_aPrevResSeq[5];
    char m_cPrevICode;

};

/**
 * This is a special form of a point cloud that encapsulates a molecular pdb file
 * \author Stefan Birmanns
 */
template<class T> class svt_point_cloud_pdb : public svt_point_cloud<T>
{
  protected:

    vector< svt_point_cloud_atom > m_aAtoms;
    vector< svt_point_cloud_bond > m_aBonds;


    vector< const char * >  m_aAtomNames;
    vector< const char * >  m_aAtomResnames;
    vector< char >         m_aAtomChainIDs;
    vector< char >         m_aAtomSecStructureIDs;
    vector< const char * >  m_aAtomSegmentIDs;
    vector< int >          m_aAtomResSeqs;
    vector< unsigned int > m_aAtomModels;

    // should the atoms belonging to a water molecule be ignored in clustering this structure?
    bool m_bIgnoreWater;

    // shall we load HETATM records?
    bool m_bLoadHetAtm;

    // enumeration functions;
    unsigned int m_iAtomEnum;
    unsigned int m_iBondEnum;

    //the number of carbon alpha in the molecule, the number of backbone atoms in the molecule
    unsigned int m_iNumCAAtoms;
    unsigned int m_iNumBBAtoms;

    // was the secondary structure information for this molecule already specified?
    bool m_bSecStructAvailable;

    // were does the secondary structure come from: 0 not available (not yet computed/read); 1, from pdb, 2, from stride
    svt_secStructSource m_eSecStructSource;

    // distance cutoff value that determines the maximum distance for which consecutive residues are connected
    float m_fTraceCutoff;

    // filename of pdb structure
    string m_pFilename;

    //graph distance in the psf file
    //at position i,j is placed the graph distance between atom i and atom j
    svt_matrix < Real64 > m_oGraphDists;

    //was the psf file been read? without the bonds information the m_oGraphDist cannot be computed
    bool m_bPSFRead;

    // Remarks relative to the pdb file (read form the original file, also added durring execution)
    vector< string > m_aRemarks;

    //list of modified residues
    vector< vector<string> > m_aModResidues;

    // list of all the secondary Structure elements
    vector <svt_sse_pdb> m_aSsePdb;

  public:

    /**
     * Constructor
     */
    svt_point_cloud_pdb();
    /**
     * Copy Constructor for svt_point_cloud (not pdb) object
     */
    svt_point_cloud_pdb(svt_point_cloud<T> &rOther);
    /**
     * Destructor
     */
    ~svt_point_cloud_pdb();



    /**
     * \name File format functions
     */
    //@{


    /**
     * Fetch PDB data from www.rcsb.org
     * \param pPdbID pointer to array of 4 char PDBID
     */
    bool fetchPDB(const char *pPdbID);

    /**
     * Load a pdb file.
     * \param pointer to array of char with the filename
     */
    void loadPDB(const char *pFilename);

    /**
     * Write pdb file.
     * \param pFilename pointer to array of char with the filename
     * \param bAppend if true, the pdb structure as append at the end of an existing structure file
     * \param bFirst if true only the first frame is written to the pdb file
     */
    void writePDB(const char *pFilename, bool bAppend = false, bool bAddRemarks = true, bool bFirst = false);

    /**
     * Opens and reads a psf file.
     * \param pFilename pointer to array of char with the filename
     */
    bool loadPSF(const char *pFilename);

    /**
     * Write psf file
     * \param pFilename pointer to array of char with the filename
     */
    void writePSF(const char *pFilename);

    /**
     * loads a dcd file
     * \param oFilename filename
     */
    bool loadDCD(const char *pFilename);

    //@}

    /**
     * \name File Format Helper Functions
     */
    //@{

    /**
     * setupHTTPGet is a method that constructs the HTTP/1.1 GET command from it's parameters.
     * \param cDirectory the exact directory from the server
     * \param cFilename the file to retrieve
     * \param cHost the server to connect to
     * \param cGetstr the string that will contain the entire GET command
     */
    void setupHTTPGet(const char *cDirectory, const char *cFilename, const char *cHost, char *cGetstr);

    /**
     * hexValue is a method that converts a single hex digit into decimal.
     * \param cHexchar the character to convert
     * \return decimal base 10 of hex character
     */
    int hexValue(char cHexchar);

    /**
     * hextoint is a method that converts an entire hex string into decimal.
     * \param pHexstr pointer to the string to convert
     * \return decimal base 10 of entire hex string
     */
    int hextoint(char *pHexstr);

    /**
     * peekPDB is a method that retrieves the HTTP header of the chosen PDBID. It checks the first
     * line of the HTTP header to determine whether the PDBID is valid.  If it is not, it closes the
     * connection, and returns false.  Otherwise, the method retrieves the entire header and returns true.
     * \param iSocket the socket used to communicate/recv with the server
     */
    bool peekPDB(int iSocket);

    /**
     * parseLine is a method that reads one line from pdb data, parses the information, and
     * integrates it into the pdb object.
     * \param pLine pointer to one line of data
     * \param iModel integer of which molecular model
     * \param iCount integer of number of atoms
     * \param iPrevSeq integer of previous residue sequence
     * \param iOrdSeq integer of ordinal sequence
     * \param cPrevChain character of previous chain ID
     */
    void parseLine(const char *pLine, unsigned int &iModel, unsigned int &iCount, int &iPrevSeq, int &iOrdSeq, char &cPrevChain);

    //@}

    /**
     * \name Data management
     * With this functions you can add or delete atoms, etc.
     */
    //@{

    /**
     * add atom
     */
    void addAtom(svt_point_cloud_atom &rAtom, svt_vector4<Real64> &rVec);

    /**
     * get atom with a special index (attention this is not the index in the pdb file, because this index is not unique in a lot of pdb files. svt_pdb simply renumbers all the atoms to get a unique index for each atom).
     * \param iIndex index of the atom you want to get
     * \return reference to the svt_atom object you want to get.
     */
    svt_point_cloud_atom &atom(unsigned int i);
    /**
     * get atom with a special index (attention this is not the index in the pdb file, because this index is not unique in a lot of pdb files. svt_pdb simply renumbers all the atoms to get a unique index for each atom).
     * \param iIndex index of the atom you want to get
     * \return pointer to the svt_atom object you want to get.
     */
    svt_point_cloud_atom *getAtom(unsigned int i);

    /**
     * Delete all atoms.
     * Attention: Please keep in mind that the internal data structures are using stl vectors. In most stl implementations the vectors cannot shrink - they will increase
     * the memory allocation every time one stores an element, but will not shrink if one deletes elements. In addition, to keep the operating system calls for the
     * memory management at a minimum, they will allocate more memory than what they would actually need to store the current elements. Most implementations will double the amount
     * of memory that they currently hold, whenever they exhaust their current capacity.
     * For the function here, that means that if you store 2mio atoms in the point_cloud_pdb and then delete all atoms, the memory consumption will still be for 2mio atoms. If you
     * from now on store only 30 atoms, it does not matter, the memory consumption will not go down. If you would like to store fluctuating numbers of atoms, please dynamically
     * create new objects on the heap with new and delete them again.
     */
    void deleteAllAtoms();

    /**
     * add a bond
     * \param iAtomA index of atom A
     * \param iAtomB index of atom B
     */
    void addBond(unsigned int iAtomA, unsigned int iAtomB);
    /**
     * is there a bond between two atoms?
     * \param iAtomA index of atom A
     * \param iAtomB index of atom B
     */
    bool isBond(unsigned int iAtomA, unsigned int iAtomB);
    /**
     * remove bond between two atoms?
     * \param iAtomA index of atom A
     * \param iAtomB index of atom B
     */
    void delBond(unsigned int iAtomA, unsigned int iAtomB);


    /**
     * find the bond object connecting two atoms
     * \param pA pointer to atom a
     * \param pB pointer to atom b
     * \return pointer to svt_point_cloud_bond object. If there is no such bond between these atoms, the return value will be NULL
     */
    svt_point_cloud_bond *findBond(svt_point_cloud_atom *pA, svt_point_cloud_atom *pB);

    /**
     * get number of bonds
     * \return number of bonds
     */
    unsigned int getBondsNum();
    /**
     * get bonds
     * \return pointer to array of bonds
     */
    vector< svt_point_cloud_bond > &getBonds();


    /**
     * delete all bonds
     */
    void deleteAllBonds();

    /**
     * Append another pdb structure
     * \param rOther reference to other pdb point cloud
     * \param iModel the model number of the other, appended pdb structure (default: 0)
     * \return reference to point_cloud_pdb object
     */
    svt_point_cloud_pdb<T> &append(svt_point_cloud_pdb<T> &rOther, int iModel = -1);

    /**
     * Add a remark
     * \param oRemark the remark
     */
    void addRemark(string oRemark);

    //@}

    /**
     * \name Enumeration functions
     * With this functions you can search the svt_pdb object for atoms or bonds with certain properties.
     */
    //@{

    /**
     * reset atom enumeration functions. This functions set the atom enumeration index back to the first atom.
     */
    void resetAtomEnum();
    /**
     * look for all atoms of a certain type. This function will return the next atom of this type, starting by a enumeration index. This index will move on every time you use an enumeration function. You can reset the index with resetAtomEnum().
     * \param pType pointer to an char array with the name of the atom
     */
    int enumAtomType(const char *pType);
    /**
     * look for all atoms of a certain residue name. This function will return the next atom with this residue name, starting by a enumeration index. This index will move on every time you use an enumeration function. You can reset the index with resetAtomEnum().
     * \param pResname pointer to an char array with the residue name
     */
    int enumAtomResname(const char *pResname);
    /**
     * look for all atoms with a certain chain id. This function will return the next atom of this type, starting by a enumeration index. This index will move on every time you use an enumeration function. You can reset the index with resetAtomEnum().
     * \param pType pointer to an char array with the name of the atom
     */
    int enumAtomChainID(char cChainID);
    /**
     * look for all atoms with a certain segment id. This function will return the next atom of this type, starting by a enumeration index. This index will move on every time you use an enumeration function. You can reset the index with resetAtomEnum().
     * \param pType pointer to an char array with the name of the atom
     */
    int enumAtomSegmentID(const char *pSegmentID);
    /**
     * look for all atoms with a certain residue sequence number. This function will return the next atom of this type, starting by a enumeration index. This index will move on every time you use an enumeration function. You can reset the index with resetAtomEnum().
     * \param pType pointer to an char array with the name of the atom
     */
    int enumAtomResidueSeq(int iResidueSeq);
    /**
     * look for all atoms with a certain model number. This function will return the next atom of this type, starting by a enumeration index. This index will move on every time you use an enumeration function. You can reset the index with resetAtomEnum().
     * \param int with the model number
     */
    int enumAtomModel(unsigned int iModel);
    /**
     * look for all atoms. This function will return the next atom, starting by a enumeration index. This index will move on every time you use an enumeration function. You can reset the index with resetAtomEnum().
     */
    int enumAtomAll();
    /**
     * look for all atoms part of a nucleotide. This function will return the next atom, starting by a enumeration index. This index will move on every time you use an enumeration function. You can reset the index with resetAtomEnum().
     */
    int enumAtomNucleotide();
    /**
     * look for all HET atoms. This function will return the next atom, starting by a enumeration index. This index will move on every time you use an enumeration function. You can reset the index with resetAtomEnum().
     */
    int enumAtomHet();
    /**
     * look for all water atoms. This function will return the next atom, starting by a enumeration index. This index will move on every time you use an enumeration function. You can reset the index with resetAtomEnum().
     */
    int enumAtomWater();
    /**
     * reset bond enumeration functions. This functions set the bond enumeration index back to the first bond.
     */
    void resetBondEnum();
    /**
     * look for all bonds. This function will return the next bond, starting by a enumeration index. This index will move on every time you use an enumeration function. You can reset the index with resetBondEnum().
     */
    int enumBondAll();
    /**
     * look for all bonds which connect a certain atom with other atoms. This function will return the next bond connecting this atom, starting by a enumeration index. This index will move on every time you use an enumeration function. You can reset the index with resetBondEnum().
     * \param pAtom pointer to an svt_point_cloud_atom object
     */
    int enumBondAtom(svt_point_cloud_atom *pAtom);

    //@}

    /**
     * \name Information
     * With these functions one can get basic statistics about the atoms in the pdb file.
     */
    //@{

    /**
     * get an array with all the different atom types (names) in the pdb file
     * \return pointer to an vector with strings of all different atom types (names) used in the pdb. Dont delete this pointer! svt_pdb will take care of it.
     */
    vector<const char *> &getAtomNames();
    /**
     * get an array with all the different residue names in the pdb file
     * \return pointer to an vector with strings of all different residue names used in the pdb. Dont delete this pointer! svt_pdb will take care of it.
     */
    vector<const char *> &getAtomResnames();
    /**
     * get an array with all the different chain id's in the pdb file
     * \return pointer to an vector with strings of all different residue names used in the pdb. Dont delete this pointer! svt_pdb will take care of it.
     */
    vector<char> &getAtomChainIDs();
    /**
     * get an array with all the different segment id's in the pdb file
     * \return pointer to an vector with strings of all different residue names used in the pdb. Dont delete this pointer! svt_pdb will take care of it.
     */
    vector<const char *> &getAtomSegmentIDs();
    /**
     * get an array with all the different residue sequence numbers in the pdb file
     * \return pointer to an vector with strings of all different residue names used in the pdb. Dont delete this pointer! svt_pdb will take care of it.
     */
    vector<int> &getAtomResidueSeqs();
    /**
     * get an array with all the different secondary structure ids in the pdb file
     * \return pointer to an vector with strings of all different secondary structure ids used in the pdb. Dont delete this pointer! svt_pdb will take care of it.
     */
    vector<char> &getAtomSecStructIDs();
    /**
     * get an array with all the different models
     * \return reference to vector with ints
     */
    vector<unsigned int> &getAtomModels();
    /**
     * get a copy of a certain model
     * \param iModel number of model
     */
    svt_point_cloud_pdb<T> *getModel(unsigned int iModel);

    /**
     * get the number of Carbon Alpha in the pdb
     * \return the number of CA in the pdb
     */
    unsigned int getNumCAAtoms();

    /**
     * get the number of Backbone atoms in the pdb
     * \return the number of backbone atoms in the pdb
     */
    unsigned int getNumBBAtoms();

    /**
     * get the minimum x coord of all atoms
     * \return minimum x coord
     */
    Real32 getMinXCoord();
    /**
     * get the minimum y coord of all atoms
     * \return minimum y coord
     */
    Real32 getMinYCoord();
    /**
     * get the minimum z coord of all atoms
     * \return minimum z coord
     */
    Real32 getMinZCoord();
    /**
     * get the maximum x coord of all atoms
     * \return maximum x coord
     */
    Real32 getMaxXCoord();
    /**
     * get the maximum y coord of all atoms
     * \return maximum y coord
     */
    Real32 getMaxYCoord();
    /**
     * get the maximum z coord of all atoms
     * \return maximum z coord
     */
    Real32 getMaxZCoord();

    /**
     * get the maximum bfactor of all atoms
     * \return maximum bfactor
     */
    Real32 getMaxTempFact();
    /**
     * get the minimum bfactor of all atoms
     * \return minimum bfactor
     */
    Real32 getMinTempFact();

    /**
     * indicates if the psf_file has been reads
     * \return a bool true(1) if the psf has been read; false(0)  in the other case
     */
    bool getPSFRead();

    /**
     * Get the distance in the graph between 2 atoms(codebook vectors)
     * \param iIndex1 an int representing the row
     * \param iIndex2 an int representing the column
     * \return a float - the distance (inside the graph) between the atom iIndex1, iIndex2
     */
    Real64 getGraphDist(unsigned int iIndex1, unsigned int iIndex2);

    /**
     * Get the pointer to the graph distance matrix
     * \return svt_matrix& object m_oGraphDists = distances (inside the graph) between any atoms
     */
    svt_matrix<Real64> &getGraphDists();

    /**
     * calculate the sphericity of the pdb
     */
    Real32 getSphericity();

    /**
     * calculate the bonds between the atoms (according to their distance the bonds are guessed)
     * \param bShowProgress indicates rather the progress bar is shown
     */
    void calcBonds(bool bShowProgress = true);
    /**
     * calculate all lists (atom types, residues, ...)
     */
    void calcAllLists();
    /**
     * calculate array with the different atom type
     */
    void calcAtomNames();
    /**
     * calculate array with the different atom resnames
     */
    void calcAtomResnames();
    /**
     * calculate array with the different atom chain ids
     */
    void calcAtomChainIDs();
    /**
     * calculate array with the different atom secondary structure ids
     */
    void calcAtomSecStructureIDs();
    /**
     * calculate array with the different atom segment ids
     */
    void calcAtomSegmentIDs();
    /**
     * calculate array with the different atom residue sequence numbers
     */
    void calcAtomResidueSeqs();
    /**
     * calculate array with the different atom model numbers (Attention: This function is NOT CALLED in calcAllLists as typically this array is automatically build during loadPDB. If a pdb is build by hand, call this function after assembly of the structure!).
     */
    void calcAtomModels();

    //@}

    /**
     * \name Secondary structure information
     */
    //@{
    /**
     * compute the ordinal chainIDs. Those identify chains in the pdb based only on a C-alpha to C-alpha
     * distance criterion and are supposed to be calculated independently from the pdb chainIDs
     */
    void calcOrdinalChainIDs();
    /**
     * calculate the secondary structure information
     */
    void calcSecStruct();
    /**
     * runs stride on one pdb containing one model or chain
     */
    void runStride();
    /**
    * Compress the secondary structure information into the ssePdb object
    * \ATTENTION: erases the existant compressed information and recompresses it from the pdb atom entry
    */
    void compressSecStruct();
    /**
     * compute and set the length of the secondary structure elements
     */
    void calcSecStructLengths();
    /**
     * is the secondary structure information already specified?
     * \param bSecStruct true if the atoms know their secondary structure membership
     */
    void setSecStructAvailable(bool bSecStruct);
    /**
     * is the secondary structure information already specified?
     * return true if the atoms know their secondary structure membership
     */
    bool getSecStructAvailable();
    /**
    * set the secondary structure source
    * \param 0 - not available, 1 - pdb, 2 - stride
    */
    void setSecStructSource(int eSecStructSource);
    /**
     * set the secondary structure source
     * \return 0 - not available, 1 - pdb, 2 - stride
     */
    int getSecStructSource();
    /**
     * get the compressed list of secondary structure- the one used to hold the information wrote in pdb in HELIX and SHEET entry
     */
    vector <svt_sse_pdb> &getSecStructCompressedList();
    /**
    * set the compressed list of secondary structure- the one used to hold the information wrote in pdb in HELIX and SHEET entry
    */
    void setSecStructCompressedList(vector <svt_sse_pdb> aSse);
    /**
     * get the sse from the list for this atom
     * \param oAtom the atom
     */
    char getSecStructFromCompressedList(svt_point_cloud_atom oAtom);
    /**
     * set the distance cutoff value that determines the maximum distance for which consecutive residues are connected
     * \param fTubeCutoff the distance cutoff value
     */
    void setTraceCutoff(float fTraceCutoff);
    /**
     * set the distance cutoff value that determines the maximum distance for which consecutive residues are connected
     * \return the distance cutoff value
     */
    float getTraceCutoff();
    /**
     * get the secondary structure
     */

    //@}

    /**
     * \name Misc functions
     */
    //@{

    /**
     * Get the maximum b-factor
     */
    Real32 getMaxTempFact() const;

    /**
     * Get the average b-factor
     */
    Real32 getAvgTempFact() const;

    /**
     * should the atoms belonging to a water molecule be ignored in clustering this structure? Default: false.
     */
    void setIgnoreWater(bool bIgnoreWater);
    /**
     * shall we read HETATM records?
     * \param bHetAtm if true hetatm atoms are read (default: true)
     */
    void setHetAtm(bool bHetAtm);


    /**
     * sample the object randomly and return a vector that refrects the probability distribution of the object
     */
    T sample();

    /**
     * Project mass
     * \param fWidth voxel width of the target map
     * \param bCA if true, only the CA atoms are evaluated
     */
    svt_volume<Real64> *projectMass(Real64 fWidth, bool bCA = false);

    /**
     * Project mass
     * \param pMap pointer to map that will hold projected structure
     * \param pTransformation pointer to transformation to apply (in Sculptor coordinates)
     * \param bCA if true, only the CA atoms are evaluated
     */
    void projectMass(svt_volume<Real64> *pMap, svt_matrix4<Real64> oTransformation = svt_matrix4<Real64>(), bool bCA = false, int iAtoms = -1);

    /**
     * Project distance
     * \param pMap pointer to map that will hold projected distances
     * \param pTransformation pointer to transformation to apply
     */
    void projectDistance(svt_volume<Real64> *pMap, svt_matrix4<Real64> oTransformation, bool bBackbone);
    /**
     * blur the pdb structure and thereby create an artificial low-resolution map
     * \param fWidth voxel width of the target map
     * \param fResolution resolution of the target map
     * \param fAdjX adjust x coordinates by this value (e.g. to uncenter) default: 0
     * \param fAdjY adjust y coordinates by this value (e.g. to uncenter) default: 0
     * \param fAdjZ adjust z coordinates by this value (e.g. to uncenter) default: 0
     * \param bProgress if true a progress bar is shown
     */
    svt_volume<Real64> *blur(Real64 fWidth, Real64 fResolution, Real64 fAdjX = 0.0, Real64 fAdjY = 0.0, Real64 fAdjZ = 0.0, bool bProgress = false);
    /**
     * blur the pdb structure and thereby create an artificial low-resolution map
     * \param fWidth voxel width of the target map
     * \param fResolution resolution of the target map
     * \param fAdjX adjust x coordinates by this value (e.g. to uncenter) default: 0
     * \param fAdjY adjust y coordinates by this value (e.g. to uncenter) default: 0
     * \param fAdjZ adjust z coordinates by this value (e.g. to uncenter) default: 0
     * \param bProgress if true a progress bar is shown
     */
    svt_volume<Real64> *blur1D(Real64 fWidth, Real64 fResolution, Real64 fAdjX = 0.0, Real64 fAdjY = 0.0, Real64 fAdjZ = 0.0, bool bProgress = false);

    /**
     * Project mass and correlation - no new volume will get created, the correlation is calculated on the fly
     * \param pMap pointer to map that holds the target map
     * \param pTransformation pointer to transformation to apply
     * \param bCA if true, only the CA atoms are evaluated
     */
    Real64 projectMassCorr(svt_volume<Real64> *pMap, svt_matrix4<Real64> oTransformation, bool bCA);

    /**
     * Blur and correlation - no new volume will get created, the correlation is calculated on the fly
     * \param pMap pointer to map that holds the target map
     * \param pTransformation pointer to transformation to apply
     * \param bCA if true, only the CA atoms are evaluated
     */
    Real64 blurCorr(svt_volume<Real64> *pKernel, svt_volume<Real64> *pMap, svt_matrix4<Real64> oTransformation, bool bCA);

    /**
     * Calculate rmsd between this and another structure. The other structure should have equal or larger size. If larger, an oligomer is assumed,
     * and the overlap with a chain is determined. The most overlapping chain is then used as reference and the rmsd is computed.
     * \param rOlig reference to other structure
     * \param bAlign align the two structures first (default: false) - valid only for the case if the two structures have the same size, ignored otherwise!
     */
    Real64 rmsd(svt_point_cloud_pdb<T> &rOlig, bool bAlign = false, Select iSelection = ALL, bool bShowProgress = true);

    /**
    * Calculate dRMSD (distance RMSD - intramolecular distances) between this and another structure.
    * \param rOlig reference to other structure
    */
    Real64 drmsd(svt_point_cloud_pdb<T> &rOlig);


    /**
     * Get a chain
     * \param cChainID the chain ID
     * \return returns the atoms of the chain as pointer to another svt_point_cloud_pdb object
     */
    svt_point_cloud_pdb<T> getChain(char cChainID);

    /**
    * Get the trace
    * \return returns the CA as svt_point_cloud_pdb object
    */
    svt_point_cloud_pdb<T> getTrace(unsigned int iSkip = 0);

    /**
     * gets the backbone of the structure
     * \return svt_point_cloud_pdb with the backbone
     */
    svt_point_cloud_pdb<T> getBackbone();

    /**
     * Get the number of atoms in
     * \param iSelection can be ALL, BACKBONE, or CA; with ALL- function equivalent with this->size()
     */
    unsigned int getAtomsNumber(Select iSelection);

    /**
     * Calculates and returns the Carbon Alpha (CA) of the atom's residue
     * \param iIndexAtom an integer representing the atom's index in the pdb structure
     * \return unsigned int - the index of the Carbon Alpha
     */
    unsigned int getCA(unsigned int iIndexAtom);

    /**
     * Determine it iIndexAtom is on the backbone
     * \param iIndexAtom an integer representing the atom's index in the pdb structure
     * \return unsigned int - iIndexAtom if the atom is on the backbone or the CA if not
     */
    unsigned int getBackbone(unsigned int iIndexAtom);

    /**
     * Compute the Distances between any atom of the structure using the topology - graph connections given in the PSF file
     */
    bool computePSFGraphDistMat();

    /**
     * \param bPSFRead a bool indicating if the psf files has been read;
     */
    void setPSFRead(bool bPSFRead);
    //@}


    /**
    * \name Select functions
    */
    //@{

    /**
     * select atoms that are in the residue
     * \param iResidueSeq the number of the residue
     */
    void selectAtomResidueSeq(int iResidueSeq);

    /**
     * get the index of the selected/deselected atoms
     * \param bIsSelected the selection status = true - is selected; false  - is not selected
     * \return the index of the atoms with selected status indicated by bIsSelected
     */
    vector<unsigned int> getSelection(bool bIsSelected = true);

    //@}

    /**
     *\name Symmetry
     */
    //@{
    /**
     * create symmetric oligomers
     * \param order is equivalent with the number of monomers in the symmetric unit
     * \param symmetry type
     * \param symmetry axis
     * \param fOffsetAxis1 is the offset form the first axis
     * \param fOffsetAxis2 is the offset form the second axis
     * \return pdb containing the symmetric oligomere
     */
    svt_point_cloud_pdb<T> applySymmetry(unsigned int iOrder = 6, const SymmetryType eType = SYMMETRY_C, char cAxis = 'z', Real64 fOffsetAxis1 = 0, Real64 fOffsetAxis2 = 0);
    //@}

};

///////////////////////////////////////////////////////////////////////////////
// Implementation
///////////////////////////////////////////////////////////////////////////////

/**
 * Constructor
 */
template<class T>
svt_point_cloud_pdb<T>::svt_point_cloud_pdb() : svt_point_cloud<T>(),
  m_bIgnoreWater(false),
  m_bLoadHetAtm(true),
  m_iAtomEnum(0),
  m_iBondEnum(0),
  m_iNumCAAtoms(0),
  m_iNumBBAtoms(0),
  m_bSecStructAvailable(false),
  m_eSecStructSource(SSE_NOTAVAILABLE),
  m_fTraceCutoff(8.0f),
  m_pFilename(""),
  m_bPSFRead(false)
{
};
/**
 * Copy Constructor for svt_point_cloud (not pdb) object
 */
template<class T>
svt_point_cloud_pdb<T>::svt_point_cloud_pdb(svt_point_cloud<T> &rOther) : svt_point_cloud<T>(),
  m_bIgnoreWater(false),
  m_bLoadHetAtm(true),
  m_iAtomEnum(0),
  m_iBondEnum(0),
  m_iNumCAAtoms(0),
  m_iNumBBAtoms(0),
  m_bSecStructAvailable(false),
  m_eSecStructSource(SSE_NOTAVAILABLE),
  m_fTraceCutoff(8.0f),
  m_pFilename(""),
  m_bPSFRead(false)
{
  svt_point_cloud_atom oAtom;

  for (unsigned int i = 0; i < rOther.size(); i++)
    addAtom(oAtom, rOther[i]);
};


/**
 * Destructor
 */
template<class T>
svt_point_cloud_pdb<T>::~svt_point_cloud_pdb()
{
};


/**
 * product of matrix and point cloud
 */
template<class T>
inline svt_point_cloud_pdb<T> operator*(const svt_matrix4<Real64> &rM, svt_point_cloud_pdb<T> &rPC)
{
  svt_point_cloud_pdb<T> oPC(rPC);

  for (unsigned int i = 0; i < rPC.size(); i++)
    oPC[i] = rM * rPC[i];

  return oPC;
};

///////////////////////////////////////////////////////////////////////////////
// File format functions
///////////////////////////////////////////////////////////////////////////////

/**
 * Fetch a pdb file from www.rcsb.org
 * \param pPdbID pointer to array of complete file name to fetch
 */
template<class T>
bool svt_point_cloud_pdb<T>::fetchPDB(const char *pPdbID)
{
  SVTLBBO << "Fetching " << pPdbID << endl;

  //used to hold the HTTP GET command
  char cHTTPGetstr[60] = {};
  //local variable that holds whether PDB exists
  bool bPdbExist;
  //variables for the socket and data transmission
#ifdef WIN32
  SOCKET iSockfd;
#else
  int iSockfd;
#endif
  int iHexval, iStatus, iNumbytes, iVal, iNum;
  //addrinfo sets variables given by the server
  struct addrinfo oHints, *pRes;
  //hexarray stores the number of bytes (in hex format) of each chunk to be received
  char cHexarray[9] = {};
  //one character buffer used in recv commands
  char cEachchar[1];
  //buffer gap to eliminate new lines and display correct hex values and data
  char cGap[2];
  //server to connect to
  const char cHost[] = "www.rcsb.org";

  // clear all old atoms
  deleteAllAtoms();
  //variables used for parseLine and passing parsestr into parseLine
  unsigned int iModel = 0;
  unsigned int iCount = 0;
  int iPrevResSeq = -1;
  int iOrdResSeq = 0;
  char cPrevChainID = '-';
  char cParsestr[MAXLINE] = {};

  //setup socket info
  memset(&oHints, 0, sizeof oHints);
  oHints.ai_family = AF_UNSPEC;
  oHints.ai_socktype = SOCK_STREAM;
  oHints.ai_protocol = 0;
  try {

#ifdef WIN32
    WSADATA wsadata;
    if (WSAStartup(MAKEWORD(1, 1), &wsadata) != 0) {
      SVTLBBO << "fetchPDB: couldn't initialize winsockets" << endl;
      return false;
    }
#endif

    //get addr info and setup/connect to socket with checks
    if ((iStatus = getaddrinfo(cHost, "80", &oHints, &pRes)) != 0) {
      SVTLBBO << "fetchPDB: getaddrinfo failure: " << gai_strerror(iStatus) << endl;
      return false;
    }

    iSockfd = socket(pRes->ai_family, pRes->ai_socktype, pRes->ai_protocol);
#ifdef WIN32
    if (iSockfd == INVALID_SOCKET) {
      SVTLBBO << "WSAGetLastError(): " << WSAGetLastError() << endl;
      SVTLBBO << "fetchPDB: error during socket() call." << endl;
      return false;
    }
#else
    if (iSockfd < 0) {
      SVTLBBO << "fetchPDB: error during socket() call." << endl;
      return false;
    }
#endif

    if ((iStatus = connect(iSockfd, pRes->ai_addr, pRes->ai_addrlen)) != 0) {
      SVTLBBO << "fetchPDB: connect error." << endl;
      return false;
    }

    freeaddrinfo(pRes);

    //setup string for the get command using HTTP/1.1 format, ie: GET <directory><filename> HTTP/1.1\nHost: <host>\n\n
    setupHTTPGet("/pdb/files/", pPdbID, cHost, cHTTPGetstr);

    //send the GET command to the server
    iStatus = send(iSockfd, cHTTPGetstr, sizeof cHTTPGetstr, 0);

    //method to determine if PDB exists, sets global variable to false if it does not, otherwise eliminates header
    bPdbExist = peekPDB(iSockfd);

    //if PDB does not exist, this is the end of the fetchPDB method.  Report failure.
    if (!bPdbExist) {
      SVTLBBO << pPdbID << " is not in the database.  Connection closing." << endl;
      return false;
    }

    //loop to iteratively retrieve data if the PDBID exists
    while (true) {

      //loop counter for hexarray location and termination of loop
      iVal = 0;
      //get the hex value to serve as counter in data transmission - eliminate hex number in data
      while (true) {

        iNumbytes = recv(iSockfd, cEachchar, sizeof cEachchar, 0);
        //check to see what is received.  If it is a letter or number, store it.  Otherwise, its the end of of the hex chunk
        if (isalnum(*cEachchar)) {
          //if the first number retrieved is 0, we're finished, otherwise store the variable and keep going.
          if (iVal == 0 && cEachchar[0] == '0') {
#ifdef WIN32
            WSACleanup();
#else
            close(iSockfd);
#endif
            SVTLBBO << pPdbID << " fetched." << endl;
            if (iCount == 0 && this->getMaxTimestep() > 1)
              this->m_aPoints.pop_back();

            if (m_aAtomModels.size() == 0)
              m_aAtomModels.push_back(0);
            else
              sort(m_aAtomModels.begin(), m_aAtomModels.end());

            calcSecStructLengths();
            SVTLBBO << "Loaded " << m_aSsePdb.size() << " secondary structure elements." << endl;
            SVTLBBO << "Loaded " << this->size() << " atoms from " << pPdbID << "." << endl;
            this->setTimestep(0);
            return true;
          } else
            cHexarray[iVal] = cEachchar[0];

        } else {
          cHexarray[iVal] = '\0';
          break;
        }
        iVal++;
      }

      //gap of one after hex number for proper formatting of pdb data
      recv(iSockfd, cEachchar, sizeof cEachchar, 0);

      //convert hexarray (which stores the 4 or less digit hex number) into decimal and use as counter.
      iHexval = hextoint(cHexarray);

      //Determine where in the line we are.  If this is a new line, then setup new line with iNum = 0. Otherwise,
      //the previous chunk ended midline, so we must resume finishing the line from the next chunk.  This correctly
      //sets the placeholder (iNum) for the next data value in the line (cParsestr).
      if (strlen(cParsestr) == 0)
        iNum = 0;
      else
        iNum = strlen(cParsestr);

      //Using the counter (iHexval) - ie, the size of the chunk, continue receiving data until the chunk is finished.
      //Store data into a line (cParsestr).  When we finish the line, send to parseLine and reset placeholder (iNum) and string.  If
      //the chunk is ended, add a null character to the end of the string.  This will be used as indicated above to
      //determine the length and therefore the placeholder of the line.
      while (iHexval > 0) {
        recv(iSockfd, cEachchar, sizeof cEachchar, 0);
        iHexval--;
        if (*cEachchar == '\n') {
          cParsestr[iNum] = '\n';
          parseLine(cParsestr, iModel, iCount, iPrevResSeq, iOrdResSeq, cPrevChainID);
          iNum = 0;
          cParsestr[0] = '\0';
        } else {
          cParsestr[iNum] = *cEachchar;
          iNum++;
          if (iHexval == 0)
            cParsestr[iNum] = '\0';
        }
      }
      //space gap after chunk has been sent for proper formatting of the pdb data
      recv(iSockfd, cGap, sizeof cGap, 0);
    }

  } catch (int e) {
    deleteAllAtoms();
#ifdef WIN32
    WSACleanup();
#else
    close(iSockfd);
#endif
    return true;
  }
  return true;
};

/**
 * Load a pdb file.
 * \param pFilename pointer to array of char with the filename
 */
template<class T>
void svt_point_cloud_pdb<T>::loadPDB(const char *pFilename)
{
  char m_pLine[MAXLINE];
  FILE *m_pFile;
  bool bEof;

  SVTLBBO << "Opening " << pFilename << endl;
  //open file
  m_pFile = fopen(pFilename, "r");

  if (m_pFile != NULL)
    bEof = (bool)(feof(m_pFile));
  else
    bEof = true;

  if (!bEof)
    m_pFilename.assign(pFilename);
  else {
    SVTLBBO << "Error opening file: " << pFilename << endl;
    return;
  }

  // clear all old atoms
  deleteAllAtoms();
  unsigned int iModel = 0;
  unsigned int iCount = 0;
  int iPrevResSeq = -1;
  int iOrdResSeq = 0;
  char cPrevChainID = '-';

  // read until next atom record reached
  while (true) {

    if (fgets(m_pLine, MAXLINE, m_pFile) == NULL) {
      if (iCount == 0 && this->getMaxTimestep() > 1)
        this->m_aPoints.pop_back();

      if (m_aAtomModels.size() == 0)
        m_aAtomModels.push_back(0);
      else
        sort(m_aAtomModels.begin(), m_aAtomModels.end());

      this->setTimestep(0);

      calcOrdinalChainIDs();
      calcSecStructLengths();

      SVTLBBO << "Loaded " << m_aSsePdb.size() << " secondary structure elements." << endl;
      SVTLBBO << "Loaded " << this->size() << " atoms from file " << pFilename << "." << endl;
      return;
    } else
      parseLine(m_pLine, iModel, iCount, iPrevResSeq, iOrdResSeq, cPrevChainID);

  }

};


/**
 * Write pdb file.
 * \param pFilename pointer to array of char with the filename
 * \param bAppend if true, the pdb structure as append at the end of an existing structure file
 * \param bFirst if true only the first frame is written to the pdb file
 */
template<class T>
void svt_point_cloud_pdb<T>::writePDB(const char *pFilename, bool bAppend, bool bAddRemarks, bool bFirst)
{
  unsigned int iNum = m_aAtoms.size();

  FILE *pFile = NULL;
  if (!bAppend)
    pFile = fopen(pFilename, "wt");
  else {
    pFile = fopen(pFilename, "a");
  }

  if (pFile == NULL) {
    SVTLBBO << "Error: Unable to write to file " << pFilename << "! Please check permissions, disk space, etc." << endl;
    return;
  }

  if (bAddRemarks) {
    for (unsigned int iIndexRemark = 0; iIndexRemark < m_aRemarks.size(); iIndexRemark++) {
      fprintf(pFile, "REMARK %s", m_aRemarks[iIndexRemark].c_str());
    }
  }
  unsigned int iModel = m_aAtoms[0].getModel();

  //code could be used if the modified res are to be outputed into the PDB
  //     for (unsigned int iModRes=0; iModRes<m_aModResidues.size(); iModRes++ )
  //     {
  //  fprintf(pFile, "MODRES %s %s %s%s%s%s  %s\n",m_aModResidues[iModRes][0].c_str(),m_aModResidues[iModRes][1].c_str(), m_aModResidues[iModRes][2].c_str(),
  //    m_aModResidues[iModRes][3].c_str(),m_aModResidues[iModRes][4].c_str(),m_aModResidues[iModRes][5].c_str(),m_aModResidues[iModRes][6].c_str() );
  //     }
  //

  if (m_eSecStructSource == SSE_OTHER)
    compressSecStruct();

  int iNoHelices = 0;
  int iNoSheets = 0;
  for (unsigned int iSse = 0; iSse < m_aSsePdb.size(); iSse++) {
    if (strcmp(m_aSsePdb[iSse].m_aType, "HELIX ") == 0) { // they are the same
      iNoHelices++;
      fprintf(pFile, "%s %3i %3s %3s %c %4i%c %3s %c %4i%c%2i%30s %5i%4s\n",
              m_aSsePdb[iSse].m_aType,
              m_aSsePdb[iSse].m_iNum % 1000,
              m_aSsePdb[iSse].m_aID,
              m_aSsePdb[iSse].m_aInitialResname,
              m_aSsePdb[iSse].m_cInitialResChainID,
              m_aSsePdb[iSse].m_iInitialResSeq,
              m_aSsePdb[iSse].m_cInitialICode,
              m_aSsePdb[iSse].m_aTerminalResname,
              m_aSsePdb[iSse].m_cTerminalResChainID,
              m_aSsePdb[iSse].m_iTerminalResSeq,
              m_aSsePdb[iSse].m_cTerminalICode,
              m_aSsePdb[iSse].m_iClass,
              m_aSsePdb[iSse].m_aComment,
              m_aSsePdb[iSse].m_iLen,
              " "
             );
    }

  }

  for (unsigned int iSse = 0; iSse < m_aSsePdb.size(); iSse++) {
    if (strcmp(m_aSsePdb[iSse].m_aType, "SHEET ") == 0) { // they are the same
      iNoSheets++;
      fprintf(pFile, "%s %3i %3s%2i %3s %c%4i%c %3s %c%4i%c%2d %4s%3s %c%4s%c %4s%3s %c%4s%c%10s\n",
              m_aSsePdb[iSse].m_aType,
              m_aSsePdb[iSse].m_iNum,
              m_aSsePdb[iSse].m_aID,
              m_aSsePdb[iSse].m_iNumStrands,
              m_aSsePdb[iSse].m_aInitialResname,
              m_aSsePdb[iSse].m_cInitialResChainID,
              m_aSsePdb[iSse].m_iInitialResSeq,
              m_aSsePdb[iSse].m_cInitialICode,
              m_aSsePdb[iSse].m_aTerminalResname,
              m_aSsePdb[iSse].m_cTerminalResChainID,
              m_aSsePdb[iSse].m_iTerminalResSeq,
              m_aSsePdb[iSse].m_cTerminalICode,
              m_aSsePdb[iSse].m_iSense,
              m_aSsePdb[iSse].m_aCurAtom,
              m_aSsePdb[iSse].m_aCurResname,
              m_aSsePdb[iSse].m_cCurResChainID,
              m_aSsePdb[iSse].m_aCurResSeq,
              m_aSsePdb[iSse].m_cCurICode,
              m_aSsePdb[iSse].m_aPrevAtom,
              m_aSsePdb[iSse].m_aPrevResname,
              m_aSsePdb[iSse].m_cPrevResChainID,
              m_aSsePdb[iSse].m_aPrevResSeq,
              m_aSsePdb[iSse].m_cPrevICode,
              " "
             );
    }
  }


  int iTS;
  char cChainID;
  char pRecordName[7];

  for (iTS = 0; iTS < this->getMaxTimestep(); iTS++) {
    //  SVTLBBO << "iTS:  " << iTS  << endl;

    if (this->getMaxTimestep() > 1)
      SVTLBBO << "  Writing frame " << iTS << " / " << this->getMaxTimestep() << endl;

    this->setTimestep(iTS);

    if (m_aAtomModels.size() > 1)
      fprintf(pFile, "MODEL        %i\n", iModel);


    for (unsigned int i = 0; i < iNum; i++) {

      if (m_aAtoms[i].getModel() != iModel) {
        fprintf(pFile, "ENDMDL\n");
        fprintf(pFile, "MODEL        %i\n", m_aAtoms[i].getModel());
        iModel = m_aAtoms[i].getModel();
      }

      cChainID = m_aAtoms[i].getChainID();
      if (cChainID == '-')
        cChainID = ' ';

      if (m_aAtoms[i].getHetAtm())
        sprintf(pRecordName, "%s", "HETATM\0");
      else
        sprintf(pRecordName, "%s", "ATOM  \0");

      fprintf(pFile, "%s%5i %2s%c%c%c%-2s %c%4i%c   %8.*f%8.*f%8.*f%6.2f%6.2f %3s  %4s%2s%2s\n",
              pRecordName,
              m_aAtoms[i].getPDBIndex(),
              m_aAtoms[i].getName(),
              m_aAtoms[i].getRemoteness(),
              m_aAtoms[i].getBranch(),
              m_aAtoms[i].getAltLoc(),
              m_aAtoms[i].getResname(),
              cChainID,
              m_aAtoms[i].getResidueSeq(),
              m_aAtoms[i].getICode(),
              coord_precision(this->getPoint(i).x()), this->getPoint(i).x(),
              coord_precision(this->getPoint(i).y()), this->getPoint(i).y(),
              coord_precision(this->getPoint(i).z()), this->getPoint(i).z(),
              m_aAtoms[i].getOccupancy(),
              m_aAtoms[i].getTempFact(),
              m_aAtoms[i].getNote(),
              m_aAtoms[i].getSegmentID(),
              m_aAtoms[i].getElement(),
              m_aAtoms[i].getCharge()
             );

    }

    if (m_aAtomModels.size() > 1) {
      fprintf(pFile, "ENDMDL\n");
    }

    if (bFirst)
      break;

    if (this->getMaxTimestep() > 0)
      fprintf(pFile, "END\n");
  }


  this->setTimestep(0);

  fclose(pFile);
};

/**
 * opens and reads a psf file
 */
template<class T>
bool svt_point_cloud_pdb<T>::loadPSF(const char *pFilename)
{
  unsigned int iNrAtoms, iNrBonds;
  unsigned int i, j;

  //yes the psf file has been read
  m_bPSFRead = true;


  // do we have to create also new coordinate entries?
  bool bNewAtoms;
  if (m_aAtoms.size() > 0)
    bNewAtoms = false;
  else
    bNewAtoms = true;

  // open file
  svt_column_reader oReader(pFilename);
  if (oReader.eof()) {
    cout << "loadPSF: cannot load psf file!" << endl;
    return false;
  }

  // load the header
  oReader.readLine();
  char *pBuffer = oReader.extractString(0, 3);
  if (pBuffer[0] != 'P' && pBuffer[1] != 'S' && pBuffer[2] != 'F') {
    SVTLBBO << "loadPSF: cannot load psf file - first string != PSF " << endl;
    return false;
  }
  // read empty line after PSF
  oReader.readLine();

  // number of lines in the header
  unsigned int iNrLines;
  oReader.readLine();
  iNrLines = oReader.extractInt(0, 7);
  for (i = 0; i < iNrLines; i++)
    oReader.readLine();
  // read empty line after header
  oReader.readLine();

  // read the number of atoms in the file
  oReader.readLine();
  iNrAtoms = oReader.extractInt(0, 7);
  unsigned int iAtomNr = 0;
  char pName[3];
  pName[0] = 0;
  pName[1] = 0;
  pName[2] = 0;
  // read the atoms
  for (i = 0; i < iNrAtoms; ++i) {
    oReader.readLine();

    svt_point_cloud_atom oAtom;

    // atom number
    iAtomNr = oReader.extractInt(0, 7);
    oAtom.setPDBIndex(iAtomNr);
    // segment id
    oAtom.setSegmentID(oReader.extractString(9, 12));
    // residue number
    oAtom.setResidueSeq(oReader.extractInt(14, 17));
    // residue name
    oAtom.setResname(oReader.extractString(19, 21));
    // atom name
    pName[0] = oReader.extractChar(23);
    pName[1] = oReader.extractChar(24);
    oAtom.setName(pName);
    // atom remoteness
    oAtom.setRemoteness(oReader.extractChar(25));
    oAtom.setBranch(oReader.extractChar(26));
    // atom type
    //oReader.extractString( 30, 32 );
    // atom charge
    oAtom.setCharge(oReader.extractString(35, 48));
    // atom mass
    oAtom.setMass(oReader.extractReal32(49, 59));

    // some defaults
    oAtom.setAltLoc(' ');
    oAtom.setChainID('A');
    oAtom.setICode(' ');
    oAtom.setOccupancy(1.0);
    oAtom.setTempFact(1.0);

    //create atom
    if (bNewAtoms) {
      svt_vector4<Real64> oNull(0.0, 0.0, 0.0);

      // no coordinates in a psf file
      addAtom(oAtom, oNull);

    } else {
      //do not update information cause the pdb already has an entry relative to this atom
    }
  }
  // cut the blank line after the atom record
  oReader.readLine();

  // get the number of bonds
  oReader.readLine();
  iNrBonds = oReader.extractInt(0, 7);

  unsigned int iA, iB, iBonds;
  iBonds = 0;
  while (iBonds < iNrBonds) {

    oReader.readLine();

    for (j = 0; j < oReader.getLength() - 16; j += 16) {
      iA = oReader.extractInt(j + 0, j + 7);
      iB = oReader.extractInt(j + 8, j + 15);

      addBond(iA - 1, iB - 1);

      iBonds++;
    }
  }




  return true;
};

/**
 * Write psf file
 * \param pFilename pointer to array of char with the filename
 */
template<class T>
void svt_point_cloud_pdb<T>::writePSF(const char *pFilename)
{
  unsigned int i, k = 0;

  FILE *pFout = fopen(pFilename, "w");

  fprintf(pFout, "PSF \n");
  fprintf(pFout, " \n");
  fprintf(pFout, "       1 !NTITLE\n");
  fprintf(pFout, " REMARK Connectivity file created by SVT/Sculptor - Attention: Atm, only the connectivity information is valid in this file!\n");
  fprintf(pFout, " \n");
  fprintf(pFout, "%8d !NATOM\n", this->size());

  for (i = 0; i < this->size(); i++)
    fprintf(pFout, "%8d %4s %-4d %4s QVOL QVOL   0.000000       0.00000           0\n", this->atom(i).getPDBIndex(), this->atom(i).getSegmentID(), this->atom(i).getResidueSeq(), this->atom(i).getResname());

  //fprintf(pFout, "%8d QVOL %-4d QVOL QVOL QVOL   0.000000       0.00000           0\n", i, i+1);

  fprintf(pFout, " \n");

  fprintf(pFout, "%8d !NBOND: bonds\n", (int)m_aBonds.size());

  for (i = 0; i < m_aBonds.size(); i++) {
    if (k == 4) {
      k = 0;
      fprintf(pFout, " \n");
    }
    k++;
    fprintf(pFout, "%8d%8d", m_aBonds[i].getIndexA() + 1, m_aBonds[i].getIndexB() + 1);
  }
  if (k > 0) fprintf(pFout, "\n");

  fprintf(pFout, "\n");
  fprintf(pFout, "       0 !NTHETA: angles\n");
  fprintf(pFout, "\n");
  fprintf(pFout, "       0 !NPHI: dihedrals\n");
  fprintf(pFout, "\n");
  fprintf(pFout, "       0 !NIMPHI: impropers\n");
  fprintf(pFout, "\n");
  fprintf(pFout, "       0 !NDON: donors\n");
  fprintf(pFout, "\n");
  fprintf(pFout, "       0 !NACC: acceptors\n");
  fprintf(pFout, "\n");
  fprintf(pFout, "       0 !NNB\n");
  fprintf(pFout, "\n");
  fprintf(pFout, "       0       0 !NGRP\n");
  fprintf(pFout, "\n");

  fclose(pFout);

  SVTLBBO << "Connectivity data written to PSF file " << pFilename << endl;

};


/**
 * loads a dcd file
 * \param oFilename filename
 */
template<class T>
bool svt_point_cloud_pdb<T>::loadDCD(const char *pFilename)
{
  /*  int j;
    int i=0;

    //
    // Open file
    //
    FILE* pFile = fopen( pFilename, "rb" );
    if (pFile == NULL)
    {
        SVTLBBO << "Cannot open dcd file: " << pFilename << endl;
        return false;
    }
    m_pFilename = pFilename;

    SVTLBBO << "Load DCD: opening " << pFilename << " for binary reading" << endl;

    //
    // First Block
    //
    Int32 iBlockSize = svt_readInt32( pFile );
    if (iBlockSize != 84)
    {
        SVTLBBO << "Load DCD:" << "\n\tDCD file may be corrupt" << "\n\tThe size of the first block is not 84" << endl;
        return false;
    }

    // CORD
    if (svt_readInt8( pFile ) != 'C' || svt_readInt8( pFile ) != 'O' || svt_readInt8( pFile ) != 'R' || svt_readInt8( pFile ) != 'D')
    {
        SVTLBBO << "Load DCD: This is not a coordinate DCD file! Skipping..." << endl;
        return false;
    }

    // HEADER
    Int32 iNrSets = svt_readInt32( pFile );
    SVTLBBO << "  NSET    = " << iNrSets << endl;
    SVTLBBO << "  ISTRT   = " << svt_readInt32( pFile ) << endl;
    SVTLBBO << "  NSAVC   = " << svt_readInt32( pFile ) << endl;
    SVTLBBO << "  5-ZEROS = " << svt_readInt32( pFile ) << ", " << svt_readInt32(pFile) << ", " << svt_readInt32(pFile) << ", " << svt_readInt32(pFile) << ", " << svt_readInt32(pFile) << endl;
    SVTLBBO << "  NAMNF   = " << svt_readInt32( pFile ) << endl;
    SVTLBBO << "  DELTA   = " << svt_readReal64( pFile ) << endl;
    SVTLBBO << "  9-ZEROS = " << svt_readInt32( pFile ) << ", " << svt_readInt32(pFile) << ", " << svt_readInt32(pFile) << ", " << svt_readInt32(pFile) << ", " << svt_readInt32(pFile) << ", " << svt_readInt32(pFile) << ", " << svt_readInt32(pFile) << ", " << svt_readInt32(pFile) << ", " << svt_readInt32(pFile) << endl;

    // End blocksize of first block
    //fseek(pFile, (11*sizeof(Int32) + sizeof(Real64) + 9*sizeof(Int32)), SEEK_SET);
    iBlockSize = svt_readInt32( pFile );
    if (iBlockSize != 84)
    {
        SVTLBBO << "Load DCD: DCD file may be corrupt, the end-size of the first block is not 84" << endl;
        return false;
    }


    //
    // Second block, title block
    //
    iBlockSize = svt_readInt32( pFile );

    // number of strings in title block
    Int32 iNrStrings = svt_readInt32( pFile );
    SVTLBBO << "  NTITLE  = " << iNrStrings << endl;
    if (iBlockSize != (4 + (80 * iNrStrings) ) )
    {
  SVTLBBO << "Load DCD: DCD file may be corrupt, the size of the title record is not (4 + 80*n)" << endl;
  return false;
    }

    char pTitle[100];
    for(i=0; i<iNrStrings; i++)
    {
  svt_readBuffer( pFile, 80, pTitle );
        pTitle[80] = 0;
  SVTLBBO << "            " << i << " : " << pTitle << endl;
    }

    //fseek( pFile, (91 + 8 + (iNrStrings * 80) + 1), SEEK_SET );
    //
    // End blocksize of second block
    iBlockSize = svt_readInt32( pFile );
    if (iBlockSize != (4 + (80 * iNrStrings) ) )
    {
        SVTLBBO << "Load DCD: DCD file may be corrupt, the end-size of the title record is not (4 + 80*n)" << endl;
        return false;
    }


    //
    // Third block, number of atoms
    //
    iBlockSize = svt_readInt32( pFile );
    if (iBlockSize != 4)
    {
        SVTLBBO << "Load DCD: DCD file may be corrupt" << endl;
        return false;
    }

    // read the number of atoms
    int iAtomNr = svt_readInt32( pFile );
    SVTLBBO << "  NATM    = " << iAtomNr << endl;
    if (iAtomNr != (int)(m_aAtoms.size()))
    {
        SVTLBBO << "Load DCD: DCD file has different number of atoms compared to current molecule?" << endl;
        return false;
    }

    // End blocksize of third block
    iBlockSize = svt_readInt32( pFile );
    if (iBlockSize != 4)
    {
        SVTLBBO << "Load DCD: DCD file may be corrupt" << endl;
        return false;
    }


    //
    // Next blocks, coordinates
    //
    for(j = 0; j < iAtomNr; j++)
    {
  this->m_aPoints[0][j].x( 0 );
  this->m_aPoints[0][j].y( 0 );
  this->m_aPoints[0][j].z( 0 );
    }

    i=0;

    // load in the coords
    while(!feof(pFile))
    {
  SVTLBBO << "  Reading timestep: " << i << "\r";

        //
  // X column
        //
  iBlockSize = svt_readInt32( pFile );
  if (iBlockSize == 0)
            break;
  if (iBlockSize != 4*iAtomNr)
  {
      SVTLBBO << "Load DCD: DCD file may be corrupt" << endl;
      SVTLBBO << "  Blocksize start of X column: " << iBlockSize << endl;
      return false;
  }
  for(j = 0; j < iAtomNr; j++)
      this->m_aPoints[i][j].x( svt_readReal32( pFile ) );
  iBlockSize = svt_readInt32( pFile );
  if (iBlockSize != 4*iAtomNr)
  {
      SVTLBBO << "Load DCD: DCD file may be corrupt" << endl;
      SVTLBBO << "  Blocksize end of X column: " << iBlockSize << endl;
      return false;
  }

  //
  // Y column
        //
        iBlockSize = svt_readInt32( pFile );
  if (iBlockSize != 4*iAtomNr)
  {
      SVTLBBO << "Load DCD: DCD file may be corrupt" << endl;
      SVTLBBO << "  Blocksize start of Y column: " << iBlockSize << endl;
      return false;
  }
        for(j = 0; j < iAtomNr; j++)
      this->m_aPoints[i][j].y( svt_readReal32( pFile ) );
  iBlockSize = svt_readInt32( pFile );
  if (iBlockSize != 4*iAtomNr)
  {
      SVTLBBO << "Load DCD: DCD file may be corrupt" << endl;
      SVTLBBO << "  Blocksize end of Y column: " << iBlockSize << endl;
      return false;
  }

  //
  // Z column
        //
  iBlockSize = svt_readInt32( pFile );
  if (iBlockSize != 4*iAtomNr)
  {
      SVTLBBO << "Load DCD: DCD file may be corrupt" << endl;
      SVTLBBO << "  Blocksize start of Z column: " << iBlockSize << endl;
      return false;
  }
  for(j = 0; j < iAtomNr; j++)
      this->m_aPoints[i][j].z( svt_readReal32( pFile ) );
  iBlockSize = svt_readInt32( pFile );
  if (iBlockSize != 4*iAtomNr)
  {
      SVTLBBO << "Load DCD: DCD file may be corrupt" << endl;
      SVTLBBO << "  Blocksize end of Z column: " << iBlockSize << endl;
      return false;
  }

        // next timestep
  this->addTimestep();
        i++;
    }

    this->setTimestep( 0 );

    fclose( pFile );

    cout << endl;
    SVTLBBO << "Load DCD: Finished loading " << i << " timesteps" << endl;
  */
  return true;
};

///////////////////////////////////////////////////////////////////////////////
// File Format Helper Functions
///////////////////////////////////////////////////////////////////////////////

/**
 * setupHTTPGet is a method that constructs the HTTP/1.1 GET command from it's parameters.
 * \param cDirectory the exact directory from the server
 * \param cFilename the file to retrieve
 * \param cHost the server to connect to
 * \param cGetstr the string that will contain the entire GET command
 */
template<class T>
void svt_point_cloud_pdb<T>::setupHTTPGet(const char *cDirectory, const char *cFilename, const char *cHost, char *cGetstr)
{
  strcat(cGetstr, "GET ");
  strcat(cGetstr, cDirectory);
  strcat(cGetstr, cFilename);
  strcat(cGetstr, " HTTP/1.1\nHost: ");
  strcat(cGetstr, cHost);
  strcat(cGetstr, "\n\n");
};

/**
 * hexValue is a method that converts a single hex digit into decimal.
 * \param cHexchar the character to convert
 * \return decimal base 10 of hex character
 */
template<class T>
int svt_point_cloud_pdb<T>::hexValue(char cHexchar)
{
  const char *pDigits = "0123456789ABCDEF";
  int iNum;
  int iLen = strlen(pDigits);
  for (iNum = 0; iNum < iLen; iNum++) {
    if (toupper((unsigned char) cHexchar) == pDigits[iNum])
      break;
  }

  return iNum;
};

/**
 * hextoint is a method that converts an entire hex string into decimal.
 * \param pHexstr pointer to the string to convert
 * \return decimal base 10 of entire hex string
 */
template<class T>
int svt_point_cloud_pdb<T>::hextoint(char *pHexstr)
{
  int iLen = strlen(pHexstr);
  int iPower = 0;
  int iResult = 0;
  int iNum;
  for (iNum = iLen - 1; iNum >= 0; iNum--) {
    iResult = iResult + (int)pow(16.0, iPower) * hexValue(pHexstr[iNum]);
    iPower++;
  }
  return iResult;
};

/**
 * peekPDB is a method that retrieves the HTTP header of the chosen PDBID. It checks the first
 * line of the HTTP header to determine whether the PDBID is valid.  If it is not, it closes the
 * connection, and returns false.  Otherwise, the method retrieves the entire header and returns true.
 * \param iSocket the socket used to communicate/recv with the server
 */
template<class T>
bool svt_point_cloud_pdb<T>::peekPDB(int iSocket)
{
  char cSingle[1];
  int iReturnChar = 0;
  //do while loop to check the first line of http header - specifically checks for any 4's in the first line
  do {
    recv(iSocket, cSingle, sizeof cSingle, 0);
    //first line contains a 4, indicative of any number of 400 server commands (all failures) so return false
    if (*cSingle == '4') {
#ifdef WIN32
      WSACleanup();
#else
      close(iSocket);
#endif
      return false;
    }

  } while (*cSingle != '\r' && *cSingle != '\n');
  //if we made it this far, then the PDB exists, so let's remove the rest of the http header before passing it on by checking for 2
  //consecutive '\r''\n' chars  - indicating the end of the header and return true.
  while (iReturnChar < 4) {
    recv(iSocket, cSingle, sizeof cSingle, 0);
    if (*cSingle == '\r' || *cSingle == '\n')
      iReturnChar++;
    else
      iReturnChar = 0;
  }
  return true;
};

/**
 * parseLine is a method that reads one line from pdb data, parses the information, and
 * integrates it into the pdb object.
 * \param pLine pointer to one line of data
 * \param iModel integer of which molecular model
 * \param iCount integer of number of atoms
 * \param iPrevSeq integer of previous residue sequence
 * \param iOrdSeq integer of ordinal sequence
 * \param cPrevChain character of previous chain ID
 */
template<class T>
void svt_point_cloud_pdb<T>::parseLine(const char *pLine, unsigned int &iModel, unsigned int &iCount, int &iPrevSeq, int &iOrdSeq, char &cPrevChain)
{
  string oRemark;
  Real64 fX, fY, fZ;
  string oBuffer = pLine;

  //if it does not have 80 characters just add empty spaces
  while (oBuffer.size() < 80)
    oBuffer += " ";

  //make sure input line is long enough for substring commands...if not, tell user.
  if (oBuffer.length() < 3) {
    SVTLBBO << "PDB file corrupt or missing data." << endl;
    return;
  }

  string pBuffer = oBuffer.substr(0, 6);
  string oEndframeA = oBuffer.substr(0, 3);
  string oEndframeB = oBuffer.substr(3, 1);

  //is it an ATOM record?
  if (pBuffer == "ATOM  " || (pBuffer == "HETATM" && m_bLoadHetAtm)) {
    svt_point_cloud_atom oAtom;

    // index
    int index = atoi((oBuffer.substr(6, 5)).c_str());
    oAtom.setPDBIndex(index);
    // model number
    oAtom.setModel(iModel);

    //now parse the coordinates
    fX = atof((oBuffer.substr(30, 8)).c_str());
    fY = atof((oBuffer.substr(38, 8)).c_str());
    fZ = atof((oBuffer.substr(46, 8)).c_str());

    if (this->getMaxTimestep() == 1) {
      this->addPoint(fX, fY, fZ);

      //type information - str
      oAtom.setName((oBuffer.substr(12, 2)).c_str());
      //type char
      oAtom.setRemoteness((oBuffer.substr(14, 1)).c_str()[0]);
      //type char
      oAtom.setBranch((oBuffer.substr(15, 1)).c_str()[0]);
      //alternate location identifier information type char
      oAtom.setAltLoc((oBuffer.substr(16, 1)).c_str()[0]);
      //resname information type str
      oAtom.setResname((oBuffer.substr(17, 4)).c_str());
      //chain id type char
      char cChainID = (oBuffer.substr(21, 1)).c_str()[0];
      if (cChainID == ' ')
        cChainID = '-';
      oAtom.setChainID(cChainID);
      if (cChainID != cPrevChain)
        iOrdSeq = 0;
      cPrevChain = cChainID;
      // residue sequence number type int
      oAtom.setResidueSeq(atoi((oBuffer.substr(22, 4)).c_str()));
      // ordinal sequence number
      if (iPrevSeq != oAtom.getResidueSeq())
        iOrdSeq++;
      iPrevSeq = oAtom.getResidueSeq();
      oAtom.setOrdResidueSeq(iOrdSeq);
      // iCode type char
      oAtom.setICode((oBuffer.substr(26, 1)).c_str()[0]);
      // fOccupancy type Real32
      oAtom.setOccupancy(atof((oBuffer.substr(54, 6)).c_str()));
      // temperature b factor type Real32
      oAtom.setTempFact(atof((oBuffer.substr(60, 6)).c_str()));
      // note type string
      oAtom.setNote((oBuffer.substr(67, 3)).c_str());
      // segment id type string
      oAtom.setSegmentID((oBuffer.substr(72, 4)).c_str());
      // element symbol type string
      oAtom.setElement((oBuffer.substr(76, 2)).c_str());
      // charge type string
      oAtom.setCharge((oBuffer.substr(78, 2)).c_str());

      // adjust mass
      oAtom.adjustMass();

      m_aAtoms.push_back(oAtom);

      if (oAtom.isCA())
        m_iNumCAAtoms++;

      if (oAtom.isBackbone())
        m_iNumBBAtoms++;
    } else {

      if (iCount >= m_aAtoms.size()) {
        SVTLBBO << "Error, file contains a trajectory with frames consisting of a varying number of atoms. Loading aborted!" << endl;
        return;
      }

      T oVec;
      oVec.x(fX);
      oVec.y(fY);
      oVec.z(fZ);

      this->setPoint(iCount, oVec);

      // increase local counter
      iCount++;

    }

    // is it an HETATM record?
    if (m_bLoadHetAtm && pBuffer == "HETATM") {
      int iModResNo;
      bool bIsModRes = false;
      int iIndexModRes = 0;
      for (unsigned int iModRes = 0; iModRes < m_aModResidues.size(); iModRes++) {
        sscanf(m_aModResidues[iModRes][3].c_str(), "%i", &iModResNo);
        if (m_aAtoms[ m_aAtoms.size() - 1].getResname() == m_aModResidues[iModRes][1] && m_aAtoms[ m_aAtoms.size() - 1].getChainID() == m_aModResidues[iModRes][2][0] &&
            m_aAtoms[ m_aAtoms.size() - 1].getResidueSeq() == iModResNo) {
          bIsModRes = true;
          iIndexModRes = iModRes;
        }
      }
      m_aAtoms[ m_aAtoms.size() - 1].setHetAtm(!bIsModRes);

      if (bIsModRes)
        m_aAtoms[ m_aAtoms.size() - 1].setResname(m_aModResidues[iIndexModRes][5].c_str());
    }

    // see if atom is in secondary structure element
    if (m_aSsePdb.size() > 0) {
      //default secondary structure element is coil
      char cSse = getSecStructFromCompressedList(m_aAtoms[ m_aAtoms.size() - 1]);
      m_aAtoms[ m_aAtoms.size() - 1].setSecStruct(cSse);
    }
  }


  // new frame?
  if (oEndframeA == "END" && oEndframeB != "M") {
    iCount = 0;
    this->addTimestep();

    SVTLBBO << "  Found frame " << this->getMaxTimestep() << endl;
  }

  // each model is treated as chain!
  if (pBuffer == "MODEL ") {
    string pBufferInt = oBuffer.substr(12, 4);
    iModel = (unsigned int)(atoi(pBufferInt.c_str()));

    SVTLBBO << "  Found Model: " << iModel << endl;

    if (find(m_aAtomModels.begin(), m_aAtomModels.end(), iModel) == m_aAtomModels.end())
      m_aAtomModels.push_back(iModel);
  }

  if (pBuffer == "REMARK") {
    string oRemark = oBuffer.substr(7, MAXLINE - 7);
    addRemark(oRemark);
  }

  if (pBuffer == "MODRES") {
    vector<string> oModResidue;
    oModResidue.push_back((oBuffer.substr(7, 4)).c_str());   //0 - idcode
    oModResidue.push_back((oBuffer.substr(12, 3)).c_str());  //1 - residue name
    oModResidue.push_back((oBuffer.substr(16, 2)).c_str());  //2 - chainid
    oModResidue.push_back((oBuffer.substr(18, 4)).c_str());  //3 - sequence number
    oModResidue.push_back((oBuffer.substr(22, 2)).c_str());  //4 -
    oModResidue.push_back((oBuffer.substr(24, 3)).c_str());  //5 - standard residues name
    oModResidue.push_back((oBuffer.substr(29, 41)).c_str());

    m_aModResidues.push_back(oModResidue);
  }

  if (pBuffer == "HELIX ") {
    svt_sse_pdb oHelix;

    strcpy(oHelix.m_aType, (oBuffer.substr(0, 6)).c_str());

    oHelix.m_iNum = atoi((oBuffer.substr(7, 3)).c_str());
    strcpy(oHelix.m_aID, (oBuffer.substr(11 , 3)).c_str());

    strcpy(oHelix.m_aInitialResname, (oBuffer.substr(15, 3)).c_str());
    oHelix.m_cInitialResChainID = (oBuffer.substr(19, 1)).c_str()[0];
    oHelix.m_iInitialResSeq = atoi((oBuffer.substr(21, 4)).c_str());
    oHelix.m_cInitialICode = (oBuffer.substr(25, 1)).c_str()[0];

    strcpy(oHelix.m_aTerminalResname, (oBuffer.substr(27, 3)).c_str());
    oHelix.m_cTerminalResChainID = (oBuffer.substr(31, 1)).c_str()[0];
    oHelix.m_iTerminalResSeq = atoi((oBuffer.substr(33, 4)).c_str());
    oHelix.m_cTerminalICode = (oBuffer.substr(37, 1)).c_str()[0];

    oHelix.m_iClass = atoi((oBuffer.substr(38, 2)).c_str());
    strcpy(oHelix.m_aComment, (oBuffer.substr(40, 30)).c_str());
    oHelix.m_iLen = atoi((oBuffer.substr(71, 5)).c_str());

    m_aSsePdb.push_back(oHelix);

    //also set the secStructureAvailable(true);
    setSecStructSource(SSE_PDB);
  }

  if (pBuffer == "SHEET ") {
    svt_sse_pdb oSheet;

    sprintf(oSheet.m_aType, "%s", (oBuffer.substr(0, 6)).c_str());

    oSheet.m_iNum = atoi((oBuffer.substr(7, 3)).c_str());
    strcpy(oSheet.m_aID, (oBuffer.substr(11 , 3)).c_str());

    oSheet.m_iNumStrands = atoi((oBuffer.substr(14, 2)).c_str());

    strcpy(oSheet.m_aInitialResname, (oBuffer.substr(17, 3)).c_str());
    oSheet.m_cInitialResChainID = (oBuffer.substr(21 , 1)).c_str()[0];
    oSheet.m_iInitialResSeq = atoi((oBuffer.substr(22, 4)).c_str());
    oSheet.m_cInitialICode = (oBuffer.substr(26, 1)).c_str()[0];

    strcpy(oSheet.m_aTerminalResname, (oBuffer.substr(28, 3)).c_str());
    oSheet.m_cTerminalResChainID = (oBuffer.substr(32, 1)).c_str()[0];
    oSheet.m_iTerminalResSeq = atoi((oBuffer.substr(33, 4)).c_str());
    oSheet.m_cTerminalICode = (oBuffer.substr(37, 1)).c_str()[0];

    oSheet.m_iSense = atoi((oBuffer.substr(38, 2)).c_str());

    strcpy(oSheet.m_aCurAtom, (oBuffer.substr(41, 4)).c_str());
    strcpy(oSheet.m_aCurResname, (oBuffer.substr(45, 3)).c_str());
    oSheet.m_cCurResChainID = (oBuffer.substr(49 , 1)).c_str()[0];
    strcpy(oSheet.m_aCurResSeq, (oBuffer.substr(50, 4)).c_str());
    oSheet.m_cCurICode = (oBuffer.substr(54, 1)).c_str()[0];

    strcpy(oSheet.m_aPrevAtom, (oBuffer.substr(56, 4)).c_str());
    strcpy(oSheet.m_aPrevResname, (oBuffer.substr(60, 3)).c_str());
    oSheet.m_cPrevResChainID = (oBuffer.substr(64 , 1)).c_str()[0];
    strcpy(oSheet.m_aPrevResSeq, (oBuffer.substr(65, 4)).c_str());
    oSheet.m_cPrevICode = (oBuffer.substr(69, 1)).c_str()[0];

    m_aSsePdb.push_back(oSheet);

    //also set the secStructureAvailable(true);
    setSecStructSource(SSE_PDB);
  }


};


///////////////////////////////////////////////////////////////////////////////
// Data management functions
///////////////////////////////////////////////////////////////////////////////

/**
 * add atom
 */
template<class T>
void svt_point_cloud_pdb<T>::addAtom(svt_point_cloud_atom &rAtom, svt_vector4<Real64> &rVec)
{
  this->addPoint(rVec.x(), rVec.y(), rVec.z());
  m_aAtoms.push_back(rAtom);
  if (m_aAtoms[m_aAtoms.size() - 1].getPDBIndex() == 0) {
    m_aAtoms[m_aAtoms.size() - 1].setPDBIndex(m_aAtoms.size());
  }

  if (rAtom.isCA())
    m_iNumCAAtoms++;

  if (rAtom.isBackbone())
    m_iNumBBAtoms++;

};

/**
 * get atom with a special index (attention this is not the index in the pdb file, because this index is not unique in a lot of pdb files. svt_pdb simply renumbers all the atoms to get a unique index for each atom).
 * \param iIndex index of the atom you want to get
 * \return reference to the svt_atom object you want to get.
 */
template<class T>
svt_point_cloud_atom &svt_point_cloud_pdb<T>::atom(unsigned int i)
{
  return m_aAtoms[i];
};

/**
 * get atom with a special index (attention this is not the index in the pdb file, because this index is not unique in a lot of pdb files. svt_pdb simply renumbers all the atoms to get a unique index for each atom).
 * \param iIndex index of the atom you want to get
 * \return pointer to the svt_atom object you want to get.
 */
template<class T>
svt_point_cloud_atom *svt_point_cloud_pdb<T>::getAtom(unsigned int i)
{
  return &m_aAtoms[i];
};

/**
 * delete all atoms
 */
template<class T>
void svt_point_cloud_pdb<T>::deleteAllAtoms()
{
  this->delAllPoints();
  m_aAtoms.clear();
  m_aBonds.clear();

  m_aAtomNames.clear();
  m_aAtomResnames.clear();
  m_aAtomChainIDs.clear();
  m_aAtomSecStructureIDs.clear();
  m_aAtomSegmentIDs.clear();
  m_aAtomResSeqs.clear();
  m_aAtomModels.clear();

  m_aRemarks.clear();
  m_aModResidues.clear();

  m_iNumCAAtoms = 0;
  m_iNumBBAtoms = 0;

  m_oGraphDists.resize(0, 0);
};

/**
 * add a bond
 * \param iAtomA index of atom A
 * \param iAtomB index of atom B
 */
template<class T>
void svt_point_cloud_pdb<T>::addBond(unsigned int iAtomA, unsigned int iAtomB)
{
  // is there already a bond?
  vector< unsigned int > aList = m_aAtoms[iAtomA].getBondList();
  for (unsigned int i = 0; i < aList.size(); i++)
    if ((m_aBonds.size() > aList[i] && (m_aBonds[aList[i]].getIndexB() == (int)(iAtomB) || m_aBonds[aList[i]].getIndexA() == (int)(iAtomB)))) {
      SVTLBBO << "DEBUG: addBond called, but there was already a bond in place?!" << endl;
      SVTLBBO << iAtomA << "\t" << iAtomB << "\t" << aList.size() << "\t" << aList[i] << "\t" << m_aBonds[aList[i]].getIndexB() << "\t" << m_aBonds[aList[i]].getIndexA()  << endl;
      //return;
    }

  m_aBonds.push_back(svt_point_cloud_bond(&m_aAtoms[iAtomA], &m_aAtoms[iAtomB], iAtomA, iAtomB));
  m_aAtoms[iAtomA].addToBondList(&m_aBonds.back(), m_aBonds.size() - 1);
  m_bPSFRead = true;
};
/**
 * is there a bond between two atoms?
 * \param iAtomA index of atom A
 * \param iAtomB index of atom B
 */
template<class T>
bool svt_point_cloud_pdb<T>::isBond(unsigned int iAtomA, unsigned int iAtomB)
{
  vector< unsigned int > aList = m_aAtoms[iAtomA].getBondList();

  if (aList.size() == 0)
    return false;

  for (unsigned int i = 0; i < aList.size(); i++)
    if (m_aBonds[aList[i]].getIndexB() == (int)(iAtomB) || m_aBonds[aList[i]].getIndexA() == (int)(iAtomB))
      return true;

  return false;
};
/**
 * remove bond between two atoms?
 * \param iAtomA index of atom A
 * \param iAtomB index of atom B
 */
template<class T>
void svt_point_cloud_pdb<T>::delBond(unsigned int iAtomA, unsigned int iAtomB)
{
  vector< unsigned int > aList = m_aAtoms[iAtomA].getBondList();

  if (aList.size() == 0) {
    SVTLBBO << "DEBUG: delBond was called, but there are no bonds there?!" << endl;
    return;
  }
  //for(unsigned int i=0; i<aList.size(); i++)
  //    SVTLBBO << "DEBUG: CV: " << iAtomA << " has a bond: " << m_aBonds[aList[i]].getIndexA() << " - " << m_aBonds[aList[i]].getIndexB() << endl;

  for (unsigned int i = 0; i < aList.size(); i++) {
    if (m_aBonds[aList[i]].getIndexB() == (int)(iAtomB) || m_aBonds[aList[i]].getIndexA() == (int)(iAtomB)) {
      m_aAtoms[iAtomA].delFromBondList(aList[i]);

      //vector< unsigned int > aNewList = m_aAtoms[iAtomA].getBondList();
      //for(unsigned int j=0; j<aNewList.size(); j++)
      //    SVTLBBO << "DEBUG: After CV: " << iAtomA << " has a bond: " << m_aBonds[aNewList[j]].getIndexA() << " - " << m_aBonds[aNewList[j]].getIndexB() << endl;

      // now correct the other bonds
      for (unsigned int j = 0; j < m_aAtoms.size(); j++)
        m_aAtoms[j].adjustBondList(aList[i]);

      m_aBonds.erase(m_aBonds.begin() + aList[i]);

      return;
    }
  }

  SVTLBBO << "DEBUG: delBond was called, but I cannot find that bond?!" << endl;
  return;
};

/**
 * find the bond object connecting two atoms
 * \param pA pointer to atom a
 * \param pB pointer to atom b
 * \return pointer to svt_point_cloud_bond object. If there is no such bond between these atoms, the return value will be NULL
 */
template<class T>
svt_point_cloud_bond *svt_point_cloud_pdb<T>::findBond(svt_point_cloud_atom *pA, svt_point_cloud_atom *pB)
{
  for (unsigned int i = 0; i < m_aBonds.size(); i++)
    if ((m_aBonds[i].getAtomA() == pA && m_aBonds[i].getAtomB() == pB) || (m_aBonds[i].getAtomA() == pB && m_aBonds[i].getAtomB() == pA))
      return &(m_aBonds[i]);

  return NULL;
};

/**
 * get number of bonds
 * \return number of bonds
 */
template<class T>
unsigned int svt_point_cloud_pdb<T>::getBondsNum()
{
  return m_aBonds.size();
}

/**
 * get bonds
 * \return pointer to array of bonds
 */
template<class T>
vector< svt_point_cloud_bond > &svt_point_cloud_pdb<T>::getBonds()
{
  return m_aBonds;
}

/**
 * delete all bonds
 */
template<class T>
void svt_point_cloud_pdb<T>::deleteAllBonds()
{
  m_aBonds.clear();

  unsigned int iAtomsNum = m_aAtoms.size();
  for (unsigned int iIndex = 0; iIndex < iAtomsNum; iIndex++) {
    //printf("Atom: %d delete %d bonds\n", iIndex, m_aAtoms[iIndex].getBondList().size() );
    m_aAtoms[iIndex].delBondList();
  }
};

/**
 * Append another pdb structure
 * \param rOther reference to other pdb point cloud
 * \param iModel the model number of the other, appended pdb structure (default: 0)
 * \return reference to point_cloud_pdb object
 */
template<class T>
svt_point_cloud_pdb<T> &svt_point_cloud_pdb<T>::append(svt_point_cloud_pdb<T> &rOther, int iModel)
{
  unsigned int iSize = this->size();
  unsigned int iAtomStart = m_aAtoms.size();

  for (unsigned int i = 0; i < rOther.size(); i++) {
    this->addAtom(rOther.atom(i), rOther[i]);
    if (iModel >= 0)
      this->atom(iSize + i).setModel(iModel);
  }

  for (unsigned int i = 0; i < rOther.m_aBonds.size(); i++)
    addBond(rOther.m_aBonds[i].getIndexA() + iAtomStart, rOther.m_aBonds[i].getIndexB() + iAtomStart);

  return *this;
};

/**
 * Add a remark
 * \param oRemark the remark as a string
 */
template<class T>
void svt_point_cloud_pdb<T>::addRemark(string oRemark)
{
  m_aRemarks.push_back(oRemark);
}


///////////////////////////////////////////////////////////////////////////////
// Enumeration functions
///////////////////////////////////////////////////////////////////////////////

/**
 * reset atom enumeration functions. This functions set the atom enumeration index back to the first atom.
 */
template<class T>
void svt_point_cloud_pdb<T>::resetAtomEnum()
{
  m_iAtomEnum = 0;
}

/**
 * look for all atoms of a certain type. This function will return the next atom of this type, starting by a enumeration index. This index will move on every time you use an enumeration function. You can reset the index with resetAtomEnum().
 * \param pType pointer to an char array with the name of the atom
 */
template<class T>
int svt_point_cloud_pdb<T>::enumAtomType(const char *pType)
{
  for (unsigned int i = m_iAtomEnum; i < m_aAtoms.size(); i++)
    if (strcmp(pType, m_aAtoms[i].getName()) == 0) {
      m_iAtomEnum = i + 1;
      return i;
    }

  return -1;
}

/**
 * look for all atoms of a certain residue name. This function will return the next atom with this residue name, starting by a enumeration index. This index will move on every time you use an enumeration function. You can reset the index with resetAtomEnum().
 * \param pResname pointer to an char array with the residue name
 */
template<class T>
int svt_point_cloud_pdb<T>::enumAtomResname(const char *pResname)
{
  for (unsigned int i = m_iAtomEnum; i < m_aAtoms.size(); i++)
    if (strcmp(pResname, m_aAtoms[i].getResname()) == 0) {
      m_iAtomEnum = i + 1;
      return i;
    }

  return -1;
}

/**
 * look for all atoms with a certain chain id. This function will return the next atom of this type, starting by a enumeration index. This index will move on every time you use an enumeration function. You can reset the index with resetAtomEnum().
 * \param pType pointer to an char array with the name of the atom
 */
template<class T>
int svt_point_cloud_pdb<T>::enumAtomChainID(char cChainID)
{
  for (unsigned int i = m_iAtomEnum; i < m_aAtoms.size(); i++)
    if (m_aAtoms[i].getChainID() == cChainID) {
      m_iAtomEnum = i + 1;
      return i;
    }

  return -1;
}

/**
 * look for all atoms with a certain segment id. This function will return the next atom of this type, starting by a enumeration index. This index will move on every time you use an enumeration function. You can reset the index with resetAtomEnum().
 * \param pType pointer to an char array with the name of the atom
 */
template<class T>
int svt_point_cloud_pdb<T>::enumAtomSegmentID(const char *pSegmentID)
{
  for (unsigned int i = m_iAtomEnum; i < m_aAtoms.size(); i++)
    if (strcmp(pSegmentID, m_aAtoms[i].getSegmentID()) == 0) {
      m_iAtomEnum = i + 1;
      return i;
    }

  return -1;
}
/**
 * look for all atoms with a certain model number. This function will return the next atom of this type, starting by a enumeration index. This index will move on every time you use an enumeration function. You can reset the index with resetAtomEnum().
 * \param int with the model number
 */
template<class T>
int svt_point_cloud_pdb<T>::enumAtomModel(unsigned int iModel)
{
  for (unsigned int i = m_iAtomEnum; i < m_aAtoms.size(); i++)
    if (m_aAtoms[i].getModel() == iModel) {
      m_iAtomEnum = i + 1;
      return i;
    }

  return -1;
};

/**
 * look for all atoms with a certain residue sequence number. This function will return the next atom of this type, starting by a enumeration index. This index will move on every time you use an enumeration function. You can reset the index with resetAtomEnum().
 * \param pType pointer to an char array with the name of the atom
 */
template<class T>
int svt_point_cloud_pdb<T>::enumAtomResidueSeq(int iResidueSeq)
{
  for (unsigned int i = m_iAtomEnum; i < m_aAtoms.size(); i++)
    if (m_aAtoms[i].getResidueSeq() == iResidueSeq) {
      m_iAtomEnum = i + 1;
      return i;
    }

  return -1;
}

/**
 * look for all atoms. This function will return the next atom, starting by a enumeration index. This index will move on every time you use an enumeration function. You can reset the index with resetAtomEnum().
 */
template<class T>
int svt_point_cloud_pdb<T>::enumAtomAll()
{
  if (m_iAtomEnum < m_aAtoms.size()) {
    m_iAtomEnum++;
    return m_iAtomEnum - 1;
  } else
    return -1;
}

/**
 * look for all atoms part of a nucleotide. This function will return the next atom, starting by a enumeration index. This index will move on every time you use an enumeration function. You can reset the index with resetAtomEnum().
 */
template<class T>
int svt_point_cloud_pdb<T>::enumAtomNucleotide()
{
  for (unsigned int i = m_iAtomEnum; i < m_aAtoms.size(); i++) {
    if (m_aAtoms[i].isNucleotide()) {
      m_iAtomEnum = i + 1;
      return (int) i;
    }
  }

  return -1;
}

/**
 * look for all HET atoms. This function will return the next atom, starting by a enumeration index. This index will move on every time you use an enumeration function. You can reset the index with resetAtomEnum().
 */
template<class T>
int svt_point_cloud_pdb<T>::enumAtomHet()
{
  for (unsigned int i = m_iAtomEnum; i < m_aAtoms.size(); i++) {
    if (m_aAtoms[i].getHetAtm()) {
      m_iAtomEnum = i + 1;
      return (int) i;
    }
  }

  return -1;
}

/**
 * look for all water atoms. This function will return the next atom, starting by a enumeration index. This index will move on every time you use an enumeration function. You can reset the index with resetAtomEnum().
 */
template<class T>
int svt_point_cloud_pdb<T>::enumAtomWater()
{
  for (unsigned int i = m_iAtomEnum; i < m_aAtoms.size(); i++) {
    if (m_aAtoms[i].isWater()) {
      m_iAtomEnum = i + 1;
      return (int) i;
    }
  }

  return -1;
}

/**
 * reset bond enumeration functions. This functions set the bond enumeration index back to the first bond.
 */
template<class T>
void svt_point_cloud_pdb<T>::resetBondEnum()
{
  m_iBondEnum = 0;
}

/**
 * look for all bonds. This function will return the next bond, starting by a enumeration index. This index will move on every time you use an enumeration function. You can reset the index with resetBondEnum().
 */
template<class T>
int svt_point_cloud_pdb<T>::enumBondAll()
{
  if (m_iBondEnum < m_aBonds.size()) {
    m_iBondEnum++;
    return m_iBondEnum - 1;
  } else
    return -1;
}

/**
 * look for all bonds which connect a certain atom with other atoms. This function will return the next bond connecting this atom, starting by a enumeration index. This index will move on every time you use an enumeration function. You can reset the index with resetBondEnum().
 * \param pAtom pointer to an svt_point_cloud_atom object
 */
template<class T>
int svt_point_cloud_pdb<T>::enumBondAtom(svt_point_cloud_atom *pAtom)
{
  for (unsigned int i = m_iBondEnum; i < m_aBonds.size(); i++)
    if (m_aBonds[i].getAtomA() == pAtom || m_aBonds[i].getAtomB() == pAtom) {
      m_iBondEnum = i + 1;
      return i;
    }

  return -1;
}

///////////////////////////////////////////////////////////////////////////////
// Info functions
///////////////////////////////////////////////////////////////////////////////

/**
 * Get an array with all the different atom types (names) in the pdb file
 * \return pointer to an vector with strings of all different atom types (names) used in the pdb. Dont delete this pointer! svt_pdb will take care of it.
 */
template<class T>
vector<const char *> &svt_point_cloud_pdb<T>::getAtomNames()
{
  return m_aAtomNames;
}
/**
 * Get an array with all the different residue names in the pdb file
 * \return pointer to an vector with strings of all different residue names used in the pdb. Dont delete this pointer! svt_pdb will take care of it.
 */
template<class T>
vector<const char *> &svt_point_cloud_pdb<T>::getAtomResnames()
{
  return m_aAtomResnames;
}
/**
 * Get an array with all the different chain id's in the pdb file
 * \return pointer to an vector with strings of all different residue names used in the pdb. Dont delete this pointer! svt_pdb will take care of it.
 */
template<class T>
vector<char> &svt_point_cloud_pdb<T>::getAtomChainIDs()
{
  return m_aAtomChainIDs;
}
/**
 * Get an array with all the different segment id's in the pdb file
 * \return pointer to an vector with strings of all different residue names used in the pdb. Dont delete this pointer! svt_pdb will take care of it.
 */
template<class T>
vector<const char *> &svt_point_cloud_pdb<T>::getAtomSegmentIDs()
{
  return m_aAtomSegmentIDs;
}
/**
 * Get an array with all the different residue sequence numbers in the pdb file
 * \return pointer to an vector with strings of all different residue names used in the pdb. Dont delete this pointer! svt_pdb will take care of it.
 */
template<class T>
vector<int> &svt_point_cloud_pdb<T>::getAtomResidueSeqs()
{
  return m_aAtomResSeqs;
}
/**
 * Get an array with all the different secondary structure ids in the pdb file
 * \return pointer to an vector with strings of all different secondary structure ids used in the pdb. Dont delete this pointer! svt_pdb will take care of it.
 */
template<class T>
vector<char> &svt_point_cloud_pdb<T>::getAtomSecStructIDs()
{
  return m_aAtomSecStructureIDs;
}
/**
 * get an array with all the different models
 * \return reference to vector with ints
 */
template<class T>
vector<unsigned int> &svt_point_cloud_pdb<T>::getAtomModels()
{
  return m_aAtomModels;
};

/**
 * get a copy of a certain model
 * \param iModel number of model
 */
template<class T>
svt_point_cloud_pdb<T> *svt_point_cloud_pdb<T>::getModel(unsigned int iModel)
{
  svt_point_cloud_pdb<T> *pNew = new svt_point_cloud_pdb<T>();

  for (unsigned int i = 0; i < this->size(); i++) {
    if (this->atom(i).getModel() == iModel)
      pNew->addAtom(this->atom(i), (*this)[i]);
  }

  return pNew;
};

/**
 * get the number of Carbon Alpha in the pdb
 * \return the number of CA in the pdb
 */
template<class T>
unsigned int svt_point_cloud_pdb<T>::getNumCAAtoms()
{
  return m_iNumCAAtoms;
};

/**
 * get the number of Backbone atoms in the pdb
 * \return the number of backbone atoms in the pdb
 */
template<class T>
unsigned int svt_point_cloud_pdb<T>::getNumBBAtoms()
{
  return m_iNumBBAtoms;
};

/**
 * Get the minimum x coord of all atoms
 * \return minimum x coord
 */
template<class T>
Real32 svt_point_cloud_pdb<T>::getMinXCoord()
{
  Real32 min_x = 10000.0f;
  unsigned int i;
  for (i = 0; i < m_aAtoms.size(); i++)
    if ((*this)[i].x() < min_x)
      min_x = (*this)[i].x();

  return min_x;
}

/**
 * Get the minimum y coord of all atoms
 * \return minimum y coord
 */
template<class T>
Real32 svt_point_cloud_pdb<T>::getMinYCoord()
{
  Real32 min_y = 10000.0f;
  unsigned int i;
  for (i = 0; i < m_aAtoms.size(); i++)
    if ((*this)[i].y() < min_y)
      min_y = (*this)[i].y();

  return min_y;
}

/**
 * Get the minimum z coord of all atoms
 * \return minimum z coord
 */
template<class T>
Real32 svt_point_cloud_pdb<T>::getMinZCoord()
{
  Real32 min_z = 10000.0f;
  unsigned int i;
  for (i = 0; i < m_aAtoms.size(); i++)
    if ((*this)[i].z() < min_z)
      min_z = (*this)[i].z();

  return min_z;
}

/**
 * Get the maximum x coord of all atoms
 * \return maximum x coord
 */
template<class T>
Real32 svt_point_cloud_pdb<T>::getMaxXCoord()
{
  Real32 max_x = -10000.0f;
  unsigned int i;
  for (i = 0; i < m_aAtoms.size(); i++)
    if ((*this)[i].x() > max_x)
      max_x = (*this)[i].x();

  return max_x;
}

/**
 * Get the maximum y coord of all atoms
 * \return maximum y coord
 */
template<class T>
Real32 svt_point_cloud_pdb<T>::getMaxYCoord()
{
  Real32 max_y = -10000.0f;
  unsigned int i;
  for (i = 0; i < m_aAtoms.size(); i++)
    if ((*this)[i].y() > max_y)
      max_y = (*this)[i].y();

  return max_y;
}

/**
 * Get the maximum z coord of all atoms
 * \return maximum z coord
 */
template<class T>
Real32 svt_point_cloud_pdb<T>::getMaxZCoord()
{
  Real32 max_z = -10000.0f;
  unsigned int i;
  for (i = 0; i < m_aAtoms.size(); i++)
    if ((*this)[i].z() > max_z)
      max_z = (*this)[i].z();

  return max_z;
}

/**
 * get the maximum bfactor of all atoms
 * \return maximum bfactor
 */
template<class T>
Real32 svt_point_cloud_pdb<T>::getMaxTempFact()
{
  Real32 fMax = -10000.0f;
  unsigned int i;
  for (i = 0; i < m_aAtoms.size(); i++)
    if (m_aAtoms[i].getTempFact() > fMax)
      fMax = m_aAtoms[i].getTempFact();

  return fMax;
}

/**
 * get the minimum bfactor of all atoms
 * \return minimum bfactor
 */
template<class T>
Real32 svt_point_cloud_pdb<T>::getMinTempFact()
{
  Real32 fMin = 10000.0f;
  unsigned int i;
  for (i = 0; i < m_aAtoms.size(); i++)
    if (m_aAtoms[i].getTempFact() > fMin)
      fMin = m_aAtoms[i].getTempFact();

  return fMin;
}

/**
 * indicates if the psf_file has been reads
 * \return a bool true(1) if the psf has been read; false(0)  in the other case
 */
template<class T>
bool svt_point_cloud_pdb<T>::getPSFRead()
{
  return m_bPSFRead;
}

/**
 * Get the distance in the graph between 2 atoms(codebook vectors)
 * \param iIndex1 an int representing the row
 * \param iIndex2 an int representing the column
 * \return a float - the distance (inside the graph) between the atom iIndex1, iIndex2
 */
template<class T>
Real64 svt_point_cloud_pdb<T>::getGraphDist(unsigned int iIndex1, unsigned int iIndex2)
{
  return m_oGraphDists[iIndex1][iIndex2];
}

/**
 * Get the pointer to the graph distance matrix
 * \return svt_matrix& object m_oGraphDists = distances (inside the graph) between any atoms
 */
template<class T>
svt_matrix<Real64> &svt_point_cloud_pdb<T>::getGraphDists()
{
  return m_oGraphDists;
}


/**
 * Calculate the sphericity of the pdb
 */
template<class T>
Real32 svt_point_cloud_pdb<T>::getSphericity()
{
  float cumx = 0;
  float cumy = 0;
  float cumz = 0;
  float cum;
  int i;

  int iAtomNum = m_aAtoms.size();

  for (i = 0; i < iAtomNum; ++i) {
    cumx += (*this)[i].x() * 1.0E10;
    cumy += (*this)[i].y() * 1.0E10;
    cumz += (*this)[i].z() * 1.0E10;
  }
  cumx /= (float)iAtomNum;
  cumy /= (float)iAtomNum;
  cumz /= (float)iAtomNum;
  cum = 0;
  for (i = 0; i < iAtomNum; ++i) {
    cum += ((*this)[i].x() * 1.0E10 - cumx) * ((*this)[i].x() * 1.0E10 - cumx) + ((*this)[i].y() * 1.0E10 - cumy) * ((*this)[i].y() * 1.0E10 - cumy) + ((*this)[i].z() * 1.0E10 - cumz) * ((*this)[i].z() * 1.0E10 - cumz);
  }
  cum /= (float)iAtomNum;
  cum = sqrt(cum);

  return iAtomNum / (0.60 * cum * cum * cum);
};

/**
 * calculate all lists (atom types, residues, ...)
 */
template<class T>
void svt_point_cloud_pdb<T>::calcAllLists()
{
  calcAtomNames();
  calcAtomResnames();
  calcAtomChainIDs();
  calcAtomSegmentIDs();
  calcAtomResidueSeqs();
  calcAtomSecStructureIDs();
};

/**
 * calculate array with the different atom type
 */
template<class T>
void svt_point_cloud_pdb<T>::calcAtomNames()
{
  bool bAdd;
  unsigned int iL1, iL2, i, j, k;
  m_aAtomNames.clear();

  for (i = 0; i < m_aAtoms.size(); i++) {
    bAdd = true;

    // search list
    for (j = 0; j < m_aAtomNames.size(); j++) {
      iL1 = strlen(m_aAtoms[i].getName());
      iL2 = strlen(m_aAtomNames[j]);

      // already in list?
      if (iL1 == iL2) {
        for (k = 0; k < iL1; k++) {
          if (m_aAtoms[i].getName()[k] != m_aAtomNames[j][k])
            break;
          if (k == (iL1 - 1))
            bAdd = false;
        }
      }
    }

    // the atom type is not in the list -> add it.
    if (bAdd)
      m_aAtomNames.push_back(m_aAtoms[i].getName());
  }
}

/**
 * calculate array with the different atom resnames
 */
template<class T>
void svt_point_cloud_pdb<T>::calcAtomResnames()
{
  bool bAdd;
  unsigned int iL1, iL2, i, j, k;
  m_aAtomResnames.clear();

  for (i = 0; i < m_aAtoms.size(); i++) {
    bAdd = true;

    // search list
    for (j = 0; j < m_aAtomResnames.size(); j++) {
      iL1 = strlen(m_aAtoms[i].getResname());
      iL2 = strlen(m_aAtomResnames[j]);

      // already in list?
      if (iL1 == iL2) {
        for (k = 0; k < iL1; k++) {
          if (m_aAtoms[i].getResname()[k] != m_aAtomResnames[j] [k])
            break;
          if (k == (iL1 - 1))
            bAdd = false;
        }
      }
    }

    // the atom type is not in the list -> add it.
    if (bAdd)
      m_aAtomResnames.push_back(m_aAtoms[i].getResname());
  }
}

/**
 * calculate array with the different atom chain ids
 */
template<class T>
void svt_point_cloud_pdb<T>::calcAtomChainIDs()
{
  bool bAdd;
  unsigned int i, j;
  m_aAtomChainIDs.clear();

  for (i = 0; i < m_aAtoms.size(); i++) {
    bAdd = true;

    // search list
    for (j = 0; j < m_aAtomChainIDs.size(); j++)
      // already in list?
      if (m_aAtoms[i].getChainID() == m_aAtomChainIDs[j])
        bAdd = false;

    // the atom type is not in the list -> add it.
    if (bAdd)
      m_aAtomChainIDs.push_back(m_aAtoms[i].getChainID());
  }
}

/**
 * calculate array with the different atom secondary structure ids
 */
template<class T>
void svt_point_cloud_pdb<T>::calcAtomSecStructureIDs()
{
  bool bAdd;
  unsigned int i, j;
  m_aAtomSecStructureIDs.clear();

  for (i = 0; i < m_aAtoms.size(); i++) {
    bAdd = true;

    // search list
    for (j = 0; j < m_aAtomSecStructureIDs.size(); j++)
      // already in list?
      if (m_aAtoms[i].getSecStruct() == m_aAtomSecStructureIDs[j])
        bAdd = false;

    // the atom type is not in the list -> add it.
    if (bAdd)
      m_aAtomSecStructureIDs.push_back(m_aAtoms[i].getSecStruct());
  }
}

/**
 * calculate array with the different atom segment ids
 */
template<class T>
void svt_point_cloud_pdb<T>::calcAtomSegmentIDs()
{
  bool bAdd;
  unsigned int iL1, iL2, i, j, k;
  m_aAtomSegmentIDs.clear();

  for (i = 0; i < m_aAtoms.size(); i++) {
    bAdd = true;

    // search list
    for (j = 0; j < m_aAtomSegmentIDs.size(); j++) {
      iL1 = strlen(m_aAtoms[i].getSegmentID());
      iL2 = strlen(m_aAtomSegmentIDs[j]);

      // already in list?
      if (iL1 == iL2) {
        for (k = 0; k < iL1; k++) {
          if (m_aAtoms[i].getSegmentID()[k] != m_aAtomSegmentIDs[j] [k])
            break;
          if (k == (iL1 - 1))
            bAdd = false;
        }
      }
    }

    // the atom type is not in the list -> add it.
    if (bAdd)
      m_aAtomSegmentIDs.push_back(m_aAtoms[i].getSegmentID());
  }
}

/**
 * calculate array with the different atom residue sequence numbers
 */
template<class T>
void svt_point_cloud_pdb<T>::calcAtomResidueSeqs()
{
  bool bAdd;
  unsigned int i, j;
  m_aAtomResSeqs.clear();

  for (i = 0; i < m_aAtoms.size(); i++) {
    bAdd = true;

    // search list
    for (j = 0; j < m_aAtomResSeqs.size(); j++)
      // already in list?
      if (m_aAtoms[i].getResidueSeq() == m_aAtomResSeqs[j])
        bAdd = false;

    // the atom type is not in the list -> add it.
    if (bAdd)
      m_aAtomResSeqs.push_back(m_aAtoms[i].getResidueSeq());
  }
}

/**
 * calculate array with the different atom model numbers (Attention: This function is NOT CALLED in calcAllLists as typically this array is automatically build during loadPDB. If a pdb is build by hand, call this function after assembly of the structure!).
 */
template<class T>
void svt_point_cloud_pdb<T>::calcAtomModels()
{
  m_aAtomModels.clear();

  for (unsigned int i = 0; i < m_aAtoms.size(); i++) {
    if (find(m_aAtomModels.begin(), m_aAtomModels.end(), m_aAtoms[i].getModel()) == m_aAtomModels.end())
      m_aAtomModels.push_back(m_aAtoms[i].getModel());
  }

  sort(m_aAtomModels.begin(), m_aAtomModels.end());
};

///////////////////////////////////////////////////////////////////////////////
// Misc functions
///////////////////////////////////////////////////////////////////////////////

/**
 * should the atoms belonging to a water molecule be ignored in clustering this structure? Default: false.
 */
template<class T>
void svt_point_cloud_pdb<T>::setIgnoreWater(bool bIgnoreWater)
{
  m_bIgnoreWater = bIgnoreWater;
};

/**
 * shall we read HETATM records?
 * \param bHetAtm if true hetatm atoms are read (default: true)
 */
template<class T>
void svt_point_cloud_pdb<T>::setHetAtm(bool bHetAtm)
{
  m_bLoadHetAtm = bHetAtm;
};

/**
 * Get the maximum b-factor
 */
template<class T>
Real32 svt_point_cloud_pdb<T>::getMaxTempFact() const
{
  unsigned int i;
  unsigned int iNum = m_aAtoms.size();
  Real32 fMax = 0.0;

  for (i = 0; i < iNum; i++)
    if (m_aAtoms[i].getTempFact() > fMax)
      fMax = m_aAtoms[i].getTempFact();

  return fMax;
}

/**
 * Get the average b-factor
 */
template<class T>
Real32 svt_point_cloud_pdb<T>::getAvgTempFact() const
{
  unsigned int i;
  unsigned int iNum = m_aAtoms.size();
  Real32 fAvg = 0.0;

  for (i = 0; i < iNum; i++)
    fAvg += m_aAtoms[i].getTempFact();

  return fAvg / (Real64)(iNum);
}

/**
 * sample the object randomly and return a vector that refrects the probability distribution of the object
 */
template<class T>
T svt_point_cloud_pdb<T>::sample()
{
  unsigned int iCount = 0;
  unsigned int iIndex = 0;

  while (iCount++ < 1000000) {
    iIndex = (unsigned int)(svt_genrand() * (double)(m_aAtoms.size()));

    if (m_bIgnoreWater) {
      if (!m_aAtoms[iIndex].isWater() && m_aAtoms[iIndex].getMass() > (Real64)(svt_genrand() * 55.0)) {
        break;
      }
    } else {
      if (m_aAtoms[iIndex].getMass() > (Real64)(svt_genrand() * 55.0))
        break;
    }
  }

  //iIndex = (unsigned int)(svt_genrand() * (double)(m_aAtoms.size()));
  return this->getPoint(iIndex);
};

/**
 * Project mass
 * \param fWidth voxel width of the target map
 */
template<class T>
svt_volume<Real64> *svt_point_cloud_pdb<T>::projectMass(Real64 fWidth, bool bCA)
{
  // bring lattice into register with origin
  T oMin = this->getMinCoord();
  T oMax = this->getMaxCoord();

  Real64 fMinx = fWidth * floor(oMin.x() / fWidth);
  Real64 fMaxx = fWidth * ceil(oMax.x() / fWidth);
  Real64 fMiny = fWidth * floor(oMin.y() / fWidth);
  Real64 fMaxy = fWidth * ceil(oMax.y() / fWidth);
  Real64 fMinz = fWidth * floor(oMin.z() / fWidth);
  Real64 fMaxz = fWidth * ceil(oMax.z() / fWidth);

  //SVTLBBO << "blur: minimal coords: " << oMin << endl;
  //SVTLBBO << "blur: maximal coords: " << oMax << endl;

  // allocate protein density map
  int iMargin = 0;
  int iExtx = (int)(ceil((fMaxx - fMinx) / fWidth)) + (2 * iMargin) + 1;
  int iExty = (int)(ceil((fMaxy - fMiny) / fWidth)) + (2 * iMargin) + 1;
  int iExtz = (int)(ceil((fMaxz - fMinz) / fWidth)) + (2 * iMargin) + 1;

  //SVTLBBO << "blur: target map size: " << iExtx << " x " << iExty << " x " << iExtz << endl;

  svt_volume<Real64> *pMap = new svt_volume<Real64>(iExtx, iExty, iExtz);
  pMap->setValue(0.0);
  pMap->setWidth(fWidth);

  //SVTLBBO << "blur: project mass" << endl;

  // interpolate structure to protein map and keep track of variability - i.e. create a volumetric map with peaks at the atomic positions...
  Real64 fVarp = 0.0;
  Real64 fGx, fGy, fGz, fA, fB, fC;
  int iX0, iY0, iZ0, iX1, iY1, iZ1;
  unsigned int i;
  unsigned int iAtomsNum = this->size();

  for (i = 0; i < iAtomsNum; i++) {
    if (!bCA || this->getAtom(i)->isCA()) {
      // compute position within grid
      fGx = iMargin + ((this->getPoint(i).x() - fMinx) / fWidth);
      fGy = iMargin + ((this->getPoint(i).y() - fMiny) / fWidth);
      fGz = iMargin + ((this->getPoint(i).z() - fMinz) / fWidth);

      iX0 = (int)(floor(fGx));
      iY0 = (int)(floor(fGy));
      iZ0 = (int)(floor(fGz));
      iX1 = iX0 + 1;
      iY1 = iY0 + 1;
      iZ1 = iZ0 + 1;

      // interpolate
      fA = (Real64)(iX1) - fGx;
      fB = (Real64)(iY1) - fGy;
      fC = (Real64)(iZ1) - fGz;

      pMap->setValue(iX0, iY0, iZ0, pMap->getValue(iX0, iY0, iZ0) + fA * fB * fC                 * m_aAtoms[i].getMass());
      fVarp += fA * fB * fC * ((1.0 - fA) * (1.0 - fA) + (1.0 - fB) * (1.0 - fB) + (1.0 - fC) * (1.0 - fC));
      pMap->setValue(iX0, iY0, iZ1, pMap->getValue(iX0, iY0, iZ1) + fA * fB * (1.0 - fC)           * m_aAtoms[i].getMass());
      fVarp += fA * fB * (1.0 - fC) * ((1.0 - fA) * (1.0 - fA) + (1.0 - fB) * (1.0 - fB) + fC * fC);
      pMap->setValue(iX0, iY1, iZ0, pMap->getValue(iX0, iY1, iZ0) + fA * (1 - fB) * fC             * m_aAtoms[i].getMass());
      fVarp += fA * (1.0 - fB) * fC * ((1.0 - fA) * (1.0 - fA) + fB * fB + (1.0 - fC) * (1.0 - fC));
      pMap->setValue(iX1, iY0, iZ0, pMap->getValue(iX1, iY0, iZ0) + (1.0 - fA) * fB * fC           * m_aAtoms[i].getMass());
      fVarp += (1.0 - fA) * fB * fC * (fA * fA + (1.0 - fB) * (1.0 - fB) + (1.0 - fC) * (1.0 - fC));
      pMap->setValue(iX0, iY1, iZ1, pMap->getValue(iX0, iY1, iZ1) + fA * (1 - fB) * (1.0 - fC)       * m_aAtoms[i].getMass());
      fVarp += fA * (1.0 - fB) * (1.0 - fC) * ((1.0 - fA) * (1.0 - fA) + fB * fB + fC * fC);
      pMap->setValue(iX1, iY1, iZ0, pMap->getValue(iX1, iY1, iZ0) + (1.0 - fA) * (1 - fB) * fC       * m_aAtoms[i].getMass());
      fVarp += (1.0 - fA) * (1.0 - fB) * fC * (fA * fA + fB * fB + (1.0 - fC) * (1.0 - fC));
      pMap->setValue(iX1, iY0, iZ1, pMap->getValue(iX1, iY0, iZ1) + (1.0 - fA) * fB * (1.0 - fC)     * m_aAtoms[i].getMass());
      fVarp += (1.0 - fA) * fB * (1.0 - fC) * (fA * fA + (1.0 - fB) * (1.0 - fB) + fC * fC);
      pMap->setValue(iX1, iY1, iZ1, pMap->getValue(iX1, iY1, iZ1) + (1.0 - fA) * (1 - fB) * (1.0 - fC) * m_aAtoms[i].getMass());
      fVarp += (1.0 - fA) * (1.0 - fB) * (1.0 - fC) * (fA * fA + fB * fB + fC * fC);
    }
  }
  fVarp /= (double)(iAtomsNum);

  //SVTLBBO << "blur: variability: " << fVarp << endl;
  //SVTLBBO << "blur: lattice smoothing (sigma = atom rmsd): " << fWidth * sqrt(fVarp) << " Angstroem" << endl;
  //SVTLBBO << "blur: convolve with kernel" << endl;
  //SVTLBBO << "blur: done." << endl;

  // set correct position of map relative to pdb
  pMap->setGrid(fMinx - ((Real64)(iMargin)*fWidth), fMiny - ((Real64)(iMargin)*fWidth), fMinz - ((Real64)(iMargin)*fWidth));

  // return
  return pMap;
};

/**
 * Project mass
 * \param pMap pointer to map that will hold projected structure
 * \param pTransformation pointer to transformation to apply
 */
template<class T>
void svt_point_cloud_pdb<T>::projectMass(svt_volume<Real64> *pMap, svt_matrix4<Real64> oTransformation, bool bCA, int iAtoms)
{

  // get properties of map
  Real64 fWidth = pMap->getWidth();
  Real64 fMinx  = pMap->getGridX();
  Real64 fMiny  = pMap->getGridY();
  Real64 fMinz  = pMap->getGridZ();

  // interpolate structure to protein map and keep track of variability - i.e. create a volumetric map with peaks at the atomic positions...
  Real64 fGx, fGy, fGz, fA, fB, fC;
  int iX0, iY0, iZ0, iX1, iY1, iZ1;
  unsigned int i;
  unsigned int iAtomsNum;
  if (iAtoms == -1)// default: use all atoms
    iAtomsNum = this->size();
  else
    iAtomsNum = iAtoms;
  svt_vector4<Real64> oCurrPoint;
  Real64 fXa, fYa, fZa;

  for (i = 0; i < iAtomsNum; i++) {
    if (!bCA || this->getAtom(i)->isCA()) {
      // Apply transform
      oCurrPoint = this->getPoint(i);

      //replaced oCurrPoint *= oTransformation; to improve time efficiency with the following 6 lines
      fXa =  oTransformation[0][0] * oCurrPoint.x() + oTransformation[0][1] * oCurrPoint.y() + oTransformation[0][2] * oCurrPoint.z() + oTransformation[0][3] * oCurrPoint.w();
      fYa =  oTransformation[1][0] * oCurrPoint.x() + oTransformation[1][1] * oCurrPoint.y() + oTransformation[1][2] * oCurrPoint.z() + oTransformation[1][3] * oCurrPoint.w();
      fZa =  oTransformation[2][0] * oCurrPoint.x() + oTransformation[2][1] * oCurrPoint.y() + oTransformation[2][2] * oCurrPoint.z() + oTransformation[2][3] * oCurrPoint.w();

      oCurrPoint.x(fXa);
      oCurrPoint.y(fYa);
      oCurrPoint.z(fZa);

      // compute position within grid
      fGx = ((oCurrPoint.x() - fMinx) / fWidth);
      fGy = ((oCurrPoint.y() - fMiny) / fWidth);
      fGz = ((oCurrPoint.z() - fMinz) / fWidth);

      iX0 = (int)(floor(fGx));
      iY0 = (int)(floor(fGy));
      iZ0 = (int)(floor(fGz));
      iX1 = iX0 + 1;
      iY1 = iY0 + 1;
      iZ1 = iZ0 + 1;

      // interpolate
      fA = (Real64)(iX1) - fGx;
      fB = (Real64)(iY1) - fGy;
      fC = (Real64)(iZ1) - fGz;

      pMap->setValue(iX0, iY0, iZ0, pMap->getValue(iX0, iY0, iZ0) + fA * fB * fC                 * m_aAtoms[i].getMass());
      pMap->setValue(iX0, iY0, iZ1, pMap->getValue(iX0, iY0, iZ1) + fA * fB * (1.0 - fC)           * m_aAtoms[i].getMass());
      pMap->setValue(iX0, iY1, iZ0, pMap->getValue(iX0, iY1, iZ0) + fA * (1 - fB) * fC             * m_aAtoms[i].getMass());
      pMap->setValue(iX1, iY0, iZ0, pMap->getValue(iX1, iY0, iZ0) + (1.0 - fA) * fB * fC           * m_aAtoms[i].getMass());
      pMap->setValue(iX0, iY1, iZ1, pMap->getValue(iX0, iY1, iZ1) + fA * (1 - fB) * (1.0 - fC)       * m_aAtoms[i].getMass());
      pMap->setValue(iX1, iY1, iZ0, pMap->getValue(iX1, iY1, iZ0) + (1.0 - fA) * (1 - fB) * fC       * m_aAtoms[i].getMass());
      pMap->setValue(iX1, iY0, iZ1, pMap->getValue(iX1, iY0, iZ1) + (1.0 - fA) * fB * (1.0 - fC)     * m_aAtoms[i].getMass());
      pMap->setValue(iX1, iY1, iZ1, pMap->getValue(iX1, iY1, iZ1) + (1.0 - fA) * (1 - fB) * (1.0 - fC) * m_aAtoms[i].getMass());
    }
  }

  // return
  return;
};

/**
 * Project distance
 * \param pMap pointer to map that will hold projected distances
 * \param pTransformation pointer to transformation to apply
 */
template<class T>
void svt_point_cloud_pdb<T>::projectDistance(svt_volume<Real64> *pMap, svt_matrix4<Real64> oTransformation,
    bool bBackbone)
{

  // set the maximum distance we want to consider
  const Real64 cMaxExtent = 15.0;

  // get properties of map
  Real64 fWidth = pMap->getWidth();
  Real64 fMinx  = pMap->getGridX();
  Real64 fMiny  = pMap->getGridY();
  Real64 fMinz  = pMap->getGridZ();

  // prepare distance kernel
  int iKernelDim    = 2 * ((int) ceil(2.0 * cMaxExtent / fWidth)) / 2 + 1;
  int iKernelCenter = iKernelDim / 2;

  svt_volume<Real64> *pDistanceKernel = new svt_volume<Real64>(iKernelDim, iKernelDim, iKernelDim);

  {
    int    i = 0;
    Real64 x, y, z;
    Real64 distance;

    for (int m = 0; m < iKernelDim; m++)
      for (int l = 0; l < iKernelDim; l++)
        for (int k = 0; k < iKernelDim; k++) {

          x = ((Real64)(k - iKernelCenter)) * fWidth;
          y = ((Real64)(l - iKernelCenter)) * fWidth;
          z = ((Real64)(m - iKernelCenter)) * fWidth;

          distance = sqrt(x * x + y * y + z * z);

          // set the central voxel to -100 so we can distinguish voxels
          // which were already examined
          if ((m == iKernelCenter) && (l == iKernelCenter) && (k == iKernelCenter))
            distance = -100.0;

          pDistanceKernel->setAt(i, distance);
          i++;
        }
  }

  // Generate local distance information for all atoms
  Real64 fGx, fGy, fGz;
  int iX, iY, iZ;
  unsigned int i;
  unsigned int iAtomsNum = this->size();
  svt_vector4<Real64> oCurrPoint;

  for (i = 0; i < iAtomsNum; i++) {

    // Check for backbone
    if (bBackbone && !this->getAtom(i)->isBackbone()) continue;

    // Apply transform
    oCurrPoint = this->getPoint(i);
    oCurrPoint *= oTransformation;

    // compute position within grid
    fGx = ((oCurrPoint.x() - fMinx) / fWidth);
    fGy = ((oCurrPoint.y() - fMiny) / fWidth);
    fGz = ((oCurrPoint.z() - fMinz) / fWidth);

    // ROUND to nearest index
    iX = (int)(nearbyint(fGx));
    iY = (int)(nearbyint(fGy));
    iZ = (int)(nearbyint(fGz));

    // Check if we already did this voxel
    if (pMap->getValue(iX, iY, iZ) < 0.0) continue;

    // Loop over distance kernel
    {
      int i = 0;
      int iKernelX, iKernelY, iKernelZ;

      for (int m = 0; m < iKernelDim; m++)
        for (int l = 0; l < iKernelDim; l++)
          for (int k = 0; k < iKernelDim; k++) {

            iKernelX = (k - iKernelCenter);
            iKernelY = (l - iKernelCenter);
            iKernelZ = (m - iKernelCenter);

            // Only replace if value is zero or larger than distance from kernel
            if ((pMap->getValue(iX + iKernelX, iY + iKernelY, iZ + iKernelZ) == 0.0) ||
                (pMap->getValue(iX + iKernelX, iY + iKernelY, iZ + iKernelZ) > pDistanceKernel->getValue(i)))
              pMap->setValue(iX + iKernelX, iY + iKernelY, iZ + iKernelZ, pDistanceKernel->getValue(i));

            i++;
          }
    }

  }

  // cleanup
  delete pDistanceKernel;

  // return
  return;
};

/**
 * Project mass and correlation - no new volume will get created, the correlation is calculated on the fly
 * \param pMap pointer to map that holds the target map
 * \param pTransformation pointer to transformation to apply
 * \param bCA if true, only the CA atoms are evaluated
 */
template<class T>
Real64 svt_point_cloud_pdb<T>::projectMassCorr(svt_volume<Real64> *pMap, svt_matrix4<Real64> oTransformation, bool bCA)
{
  // get properties of map
  Real64 fWidth = pMap->getWidth();
  Real64 fMinx  = pMap->getGridX();
  Real64 fMiny  = pMap->getGridY();
  Real64 fMinz  = pMap->getGridZ();

  Real64 fCorr = 0.0;
  Real64 fGx, fGy, fGz, fA, fB, fC;
  int iX0, iY0, iZ0, iX1, iY1, iZ1;
  unsigned int i;
  unsigned int iAtomsNum = this->size();
  svt_vector4<Real64> oCurrPoint;
  Real64 fMass;
  Real64 fNormB = 0.0;

  Real64 *pMapData = pMap->getData();
  unsigned int iSizeX = pMap->getSizeX();
  unsigned int iSizeY = pMap->getSizeY();
  unsigned int iSizeZ = pMap->getSizeZ();
  unsigned int iSizeXY = iSizeX * iSizeY;
  int iY0S, iY1S, iZ0S, iZ1S;
  Real64 fXa, fYa, fZa;

  for (i = 0; i < iAtomsNum; i++) {
    if (!bCA || this->getAtom(i)->isCA()) {
      // Apply transform
      oCurrPoint = this->getPoint(i);
      //replaced oCurrPoint *= oTransformation; to improve time efficiency with the following 6 lines
      fXa =  oTransformation[0][0] * oCurrPoint.x() + oTransformation[0][1] * oCurrPoint.y() + oTransformation[0][2] * oCurrPoint.z() + oTransformation[0][3] * oCurrPoint.w();
      fYa =  oTransformation[1][0] * oCurrPoint.x() + oTransformation[1][1] * oCurrPoint.y() + oTransformation[1][2] * oCurrPoint.z() + oTransformation[1][3] * oCurrPoint.w();
      fZa =  oTransformation[2][0] * oCurrPoint.x() + oTransformation[2][1] * oCurrPoint.y() + oTransformation[2][2] * oCurrPoint.z() + oTransformation[2][3] * oCurrPoint.w();

      oCurrPoint.x(fXa);
      oCurrPoint.y(fYa);
      oCurrPoint.z(fZa);
      //oCurrPoint *= oTransformation;

      // compute position within grid
      fGx = ((oCurrPoint.x() - fMinx) / fWidth);
      fGy = ((oCurrPoint.y() - fMiny) / fWidth);
      fGz = ((oCurrPoint.z() - fMinz) / fWidth);

      iX0 = (int)(floor(fGx));
      iY0 = (int)(floor(fGy));
      iZ0 = (int)(floor(fGz));
      iX1 = iX0 + 1;
      iY1 = iY0 + 1;
      iZ1 = iZ0 + 1;

      fMass = m_aAtoms[i].getMass();

      if (iX0 >= 0 && iX0 < (int)iSizeX && iX1 >= 0 && iX1 < (int)iSizeX  &&
          iY0 >= 0 && iY0 < (int)iSizeY && iY1 >= 0 && iY1 < (int)iSizeY  &&
          iZ0 >= 0 && iZ0 < (int)iSizeZ && iZ1 >= 0 && iZ1 < (int)iSizeZ) {

        // interpolate
        fA = (Real64)(iX1) - fGx;
        fB = (Real64)(iY1) - fGy;
        fC = (Real64)(iZ1) - fGz;

        iY0S = iY0 * iSizeX;
        iY1S = iY1 * iSizeX;
        iZ0S = iZ0 * iSizeXY;
        iZ1S = iZ1 * iSizeXY;


        fCorr += pMapData[ iX0 + iY0S + iZ0S ] * (fA        * fB        * fC        * fMass);
        fCorr += pMapData[ iX0 + iY0S + iZ1S ] * (fA        * fB        * (1.0 - fC)  * fMass);
        fCorr += pMapData[ iX0 + iY1S + iZ0S ] * (fA        * (1.0 - fB)  * fC        * fMass);
        fCorr += pMapData[ iX1 + iY0S + iZ0S ] * ((1.0 - fA)  * fB        * fC        * fMass);
        fCorr += pMapData[ iX0 + iY1S + iZ1S ] * (fA        * (1.0 - fB)  * (1.0 - fC)  * fMass);
        fCorr += pMapData[ iX1 + iY1S + iZ0S ] * ((1.0 - fA)  * (1.0 - fB)  * fC        * fMass);
        fCorr += pMapData[ iX1 + iY0S + iZ1S ] * ((1.0 - fA)  * fB        * (1.0 - fC)  * fMass);
        fCorr += pMapData[ iX1 + iY1S + iZ1S ] * ((1.0 - fA)  * (1.0 - fB)  * (1.0 - fC)  * fMass);

        fNormB += fMass * fMass;
      } else
        fNormB += fMass * fMass;
    }
  }

  Real64 fNormA = pMap->getNorm();
  fCorr /= (fNormA * fNormB);

  // return
  return fCorr;
};

/**
 * Blur and correlation - no new volume will get created, the correlation is calculated on the fly
 * \param pMap pointer to map that holds the target map
 * \param pTransformation pointer to transformation to apply
 * \param bCA if true, only the CA atoms are evaluated
 */
template<class T>
Real64 svt_point_cloud_pdb<T>::blurCorr(svt_volume<Real64> *pKernel, svt_volume<Real64> *pMap, svt_matrix4<Real64> oTransformation, bool bCA)
{
  // get properties of map
  Real64 fWidth = pMap->getWidth();
  Real64 fMinx  = pMap->getGridX();
  Real64 fMiny  = pMap->getGridY();
  Real64 fMinz  = pMap->getGridZ();

  Real64 fCorr = 0.0;
  unsigned int i;
  unsigned int iAtomsNum = this->size();
  T oCurrPoint;

  unsigned int iSize = pKernel->getSizeX();
  unsigned int iSizeS = iSize * iSize;
  Real64 fMass;
  int iKX, iKY, iKZ;
  int iDim = (int)((Real32)(pKernel->getSizeX()) * 0.5f);
  int iStart = -iDim;
  int iEnd = iDim + 1;
  Real64 fGx, fGy, fGz, fA, fB, fC;
  int iX0, iY0, iZ0, iX1, iY1, iZ1;
  Real64 fKernel;
  Real64 fMassABC, fMassAB1C, fMassA1BC, fMass1ABC, fMassA1B1C, fMass1A1BC, fMass1AB1C, fMass1A1B1C;
  int iY0KY, iY1KY, iZ0KZ, iZ1KZ, iRest;

  Real64 *pKernelData = pKernel->getData();
  Real64 *pMapData = pMap->getData();

  int iSizeX = pMap->getSizeX();
  int iSizeY = pMap->getSizeY();
  int iSizeZ = pMap->getSizeZ();
  int iSizeXY = iSizeX * iSizeY;

  int iLBd = iDim;
  int iUBdX = (int)(iSizeX) - iDim - 2;
  int iUBdY = (int)(iSizeY) - iDim - 2;
  int iUBdZ = (int)(iSizeZ) - iDim - 2;

  vector< T > &aPoints = this->getPoints();

  for (i = 0; i < iAtomsNum; i++) {
    if (!bCA || this->getAtom(i)->isCA()) {
      // Apply transform
      oCurrPoint = aPoints[i];
      oCurrPoint *= oTransformation;

      // compute position within grid
      fGx = ((oCurrPoint.x() - fMinx) / fWidth);
      fGy = ((oCurrPoint.y() - fMiny) / fWidth);
      fGz = ((oCurrPoint.z() - fMinz) / fWidth);

      iX0 = (int)(fGx);
      iY0 = (int)(fGy);
      iZ0 = (int)(fGz);
      iX1 = iX0 + 1;
      iY1 = iY0 + 1;
      iZ1 = iZ0 + 1;

      // interpolate
      fA = (Real64)(iX1) - fGx;
      fB = (Real64)(iY1) - fGy;
      fC = (Real64)(iZ1) - fGz;

      fMass = m_aAtoms[i].getMass();

      fMassABC    = fA * fB * fC                   * fMass;
      fMassAB1C   = fA * fB * (1.0 - fC)             * fMass;
      fMassA1BC   = fA * (1.0 - fB) * fC             * fMass;
      fMass1ABC   = (1.0 - fA) * fB * fC             * fMass;
      fMassA1B1C  = fA * (1.0 - fB) * (1.0 - fC)       * fMass;
      fMass1A1BC  = (1.0 - fA) * (1.0 - fB) * fC       * fMass;
      fMass1AB1C  = (1.0 - fA) * fB * (1.0 - fC)       * fMass;
      fMass1A1B1C = (1.0 - fA) * (1.0 - fB) * (1.0 - fC) * fMass;

      if (!(iX0 < iLBd || iX0 > iUBdX || iY0 < iLBd || iY0 > iUBdY || iZ0 < iLBd || iZ0 > iUBdZ)) {
        for (iKZ = iStart; iKZ < iEnd; iKZ++)
          for (iKY = iStart; iKY < iEnd; iKY++) {
            iY0KY = (iY0 + iKY) * iSizeX;
            iY1KY = (iY1 + iKY) * iSizeX;
            iZ0KZ = (iZ0 + iKZ) * iSizeXY;
            iZ1KZ = (iZ1 + iKZ) * iSizeXY;

            iRest = ((iKY - iStart) * iSize) + iKZ - iStart;

            for (iKX = iStart; iKX < iEnd; iKX++) {
              fKernel = pKernelData[((iKX - iStart) * iSizeS) + iRest];

              fCorr += pMapData[(iX0 + iKX) + iY0KY + iZ0KZ] * fKernel * fMassABC;
              fCorr += pMapData[(iX0 + iKX) + iY0KY + iZ1KZ] * fKernel * fMassAB1C;
              fCorr += pMapData[(iX0 + iKX) + iY1KY + iZ0KZ] * fKernel * fMassA1BC;
              fCorr += pMapData[(iX1 + iKX) + iY0KY + iZ0KZ] * fKernel * fMass1ABC;
              fCorr += pMapData[(iX0 + iKX) + iY1KY + iZ1KZ] * fKernel * fMassA1B1C;
              fCorr += pMapData[(iX1 + iKX) + iY1KY + iZ0KZ] * fKernel * fMass1A1BC;
              fCorr += pMapData[(iX1 + iKX) + iY0KY + iZ1KZ] * fKernel * fMass1AB1C;
              fCorr += pMapData[(iX1 + iKX) + iY1KY + iZ1KZ] * fKernel * fMass1A1B1C;
            }
          }

      } else {

        for (iKZ = iStart; iKZ < iEnd; iKZ++)
          for (iKY = iStart; iKY < iEnd; iKY++) {
            iY0KY = (iY0 + iKY) * iSizeX;
            iY1KY = (iY1 + iKY) * iSizeX;
            iZ0KZ = (iZ0 + iKZ) * iSizeXY;
            iZ1KZ = (iZ1 + iKZ) * iSizeXY;

            iRest = ((iKY - iStart) * iSize) + iKZ - iStart;

            for (iKX = iStart; iKX < iEnd; iKX++) {
              fKernel = pKernelData[((iKX - iStart) * iSizeS) + iRest];

              fCorr += pMap->getValue((iX0 + iKX) + iY0KY + iZ0KZ) * fKernel * fMassABC;
              fCorr += pMap->getValue((iX0 + iKX) + iY0KY + iZ1KZ) * fKernel * fMassAB1C;
              fCorr += pMap->getValue((iX0 + iKX) + iY1KY + iZ0KZ) * fKernel * fMassA1BC;
              fCorr += pMap->getValue((iX1 + iKX) + iY0KY + iZ0KZ) * fKernel * fMass1ABC;
              fCorr += pMap->getValue((iX0 + iKX) + iY1KY + iZ1KZ) * fKernel * fMassA1B1C;
              fCorr += pMap->getValue((iX1 + iKX) + iY1KY + iZ0KZ) * fKernel * fMass1A1BC;
              fCorr += pMap->getValue((iX1 + iKX) + iY0KY + iZ1KZ) * fKernel * fMass1AB1C;
              fCorr += pMap->getValue((iX1 + iKX) + iY1KY + iZ1KZ) * fKernel * fMass1A1B1C;
            }
          }

      }
    }
  }

  Real64 fNormA = pMap->getNorm();
  Real64 fNormB = pKernel->getNorm() * iAtomsNum;

  fCorr /= fNormA * fNormB;

  // return
  return fCorr;
};

/**
 * blur the pdb structure and thereby create an artificial low-resolution map
 * \param fWidth voxel width of the target map
 * \param fResolution resolution of the target map
 * \param fAdjX adjust x coordinates by this value (e.g. to uncenter) default: 0
 * \param fAdjY adjust y coordinates by this value (e.g. to uncenter) default: 0
 * \param fAdjZ adjust z coordinates by this value (e.g. to uncenter) default: 0
 */
template<class T>
svt_volume<Real64> *svt_point_cloud_pdb<T>::blur(Real64 fWidth, Real64 fResolution, Real64 fAdjX, Real64 fAdjY, Real64 fAdjZ, bool bProgress)
{
  // create gaussian kernel
  svt_volume<Real64> oKernel;
  oKernel.createSitusBlurringKernel(fWidth, fResolution);

  // bring lattice into register with origin
  T oMin = this->getMinCoord();
  T oMax = this->getMaxCoord();

  Real64 fMinx = (fWidth * floor((oMin.x() + fAdjX) / fWidth) - fAdjX);
  Real64 fMaxx = (fWidth * ceil((oMax.x() + fAdjX) / fWidth) - fAdjX);
  Real64 fMiny = (fWidth * floor((oMin.y() + fAdjY) / fWidth) - fAdjY);
  Real64 fMaxy = (fWidth * ceil((oMax.y() + fAdjY) / fWidth) - fAdjY);
  Real64 fMinz = (fWidth * floor((oMin.z() + fAdjZ) / fWidth) - fAdjZ);
  Real64 fMaxz = (fWidth * ceil((oMax.z() + fAdjZ) / fWidth) - fAdjZ);

  // allocate protein density map
  int iMargin = (int) ceil((Real64)(oKernel.getSizeX()) / 2.0);
  int iExtx = (int)(ceil((fMaxx - fMinx) / fWidth)) + (2 * iMargin) + 1;
  int iExty = (int)(ceil((fMaxy - fMiny) / fWidth)) + (2 * iMargin) + 1;
  int iExtz = (int)(ceil((fMaxz - fMinz) / fWidth)) + (2 * iMargin) + 1;

  //SVTLBBO << "blur: target map size: " << iExtx << " x " << iExty << " x " << iExtz << endl;

  svt_volume<Real64> *pMap = new svt_volume<Real64>(iExtx, iExty, iExtz);
  //pMap->setValue( 0.0 );
  pMap->setWidth(fWidth);

  //SVTLBBO << "blur: project mass" << endl;

  // interpolate structure to protein map and keep track of variability - i.e. create a volumetric map with peaks at the atomic positions...
  Real64 fVarp = 0.0;
  Real64 fGx, fGy, fGz, fA, fB, fC;
  int iX0, iY0, iZ0, iX1, iY1, iZ1;
  unsigned int i;
  unsigned int iAtomsNum = this->size();

  try {
    for (i = 0; i < iAtomsNum; i++) {
      // compute position within grid
      fGx = iMargin + ((this->getPoint(i).x() - fMinx) / fWidth);
      fGy = iMargin + ((this->getPoint(i).y() - fMiny) / fWidth);
      fGz = iMargin + ((this->getPoint(i).z() - fMinz) / fWidth);

      iX0 = (int)(floor(fGx));
      iY0 = (int)(floor(fGy));
      iZ0 = (int)(floor(fGz));
      iX1 = iX0 + 1;
      iY1 = iY0 + 1;
      iZ1 = iZ0 + 1;

      // interpolate
      fA = (Real64)(iX1) - fGx;
      fB = (Real64)(iY1) - fGy;
      fC = (Real64)(iZ1) - fGz;

      pMap->setValue(iX0, iY0, iZ0, pMap->getValue(iX0, iY0, iZ0) + fA * fB * fC                 * m_aAtoms[i].getMass());
      fVarp += fA * fB * fC * ((1.0 - fA) * (1.0 - fA) + (1.0 - fB) * (1.0 - fB) + (1.0 - fC) * (1.0 - fC));
      pMap->setValue(iX0, iY0, iZ1, pMap->getValue(iX0, iY0, iZ1) + fA * fB * (1.0 - fC)           * m_aAtoms[i].getMass());
      fVarp += fA * fB * (1.0 - fC) * ((1.0 - fA) * (1.0 - fA) + (1.0 - fB) * (1.0 - fB) + fC * fC);
      pMap->setValue(iX0, iY1, iZ0, pMap->getValue(iX0, iY1, iZ0) + fA * (1 - fB) * fC             * m_aAtoms[i].getMass());
      fVarp += fA * (1.0 - fB) * fC * ((1.0 - fA) * (1.0 - fA) + fB * fB + (1.0 - fC) * (1.0 - fC));
      pMap->setValue(iX1, iY0, iZ0, pMap->getValue(iX1, iY0, iZ0) + (1.0 - fA) * fB * fC           * m_aAtoms[i].getMass());
      fVarp += (1.0 - fA) * fB * fC * (fA * fA + (1.0 - fB) * (1.0 - fB) + (1.0 - fC) * (1.0 - fC));
      pMap->setValue(iX0, iY1, iZ1, pMap->getValue(iX0, iY1, iZ1) + fA * (1 - fB) * (1.0 - fC)       * m_aAtoms[i].getMass());
      fVarp += fA * (1.0 - fB) * (1.0 - fC) * ((1.0 - fA) * (1.0 - fA) + fB * fB + fC * fC);
      pMap->setValue(iX1, iY1, iZ0, pMap->getValue(iX1, iY1, iZ0) + (1.0 - fA) * (1 - fB) * fC       * m_aAtoms[i].getMass());
      fVarp += (1.0 - fA) * (1.0 - fB) * fC * (fA * fA + fB * fB + (1.0 - fC) * (1.0 - fC));
      pMap->setValue(iX1, iY0, iZ1, pMap->getValue(iX1, iY0, iZ1) + (1.0 - fA) * fB * (1.0 - fC)     * m_aAtoms[i].getMass());
      fVarp += (1.0 - fA) * fB * (1.0 - fC) * (fA * fA + (1.0 - fB) * (1.0 - fB) + fC * fC);
      pMap->setValue(iX1, iY1, iZ1, pMap->getValue(iX1, iY1, iZ1) + (1.0 - fA) * (1 - fB) * (1.0 - fC) * m_aAtoms[i].getMass());
      fVarp += (1.0 - fA) * (1.0 - fB) * (1.0 - fC) * (fA * fA + fB * fB + fC * fC);

    }
    fVarp /= (double)(iAtomsNum);


  } catch (int e) {
    delete pMap;
    return NULL;
  }

  //SVTLBBO << "blur: variability: " << fVarp << endl;
  //SVTLBBO << "blur: lattice smoothing (sigma = atom rmsd): " << fWidth * sqrt(fVarp) << " Angstroem" << endl;
  //SVTLBBO << "blur: convolve with kernel" << endl;

  // convolve
  pMap->convolve(oKernel, bProgress);

  //SVTLBBO << "blur: done." << endl;

  // set correct position of map relative to pdb
  pMap->setGrid(fMinx - ((Real64)(iMargin)*fWidth) + fAdjX, fMiny - ((Real64)(iMargin)*fWidth) + fAdjY, fMinz - ((Real64)(iMargin)*fWidth) + fAdjZ);

  // return
  return pMap;
};

/**
 * blur the pdb structure and thereby create an artificial low-resolution map
 * \param fWidth voxel width of the target map
 * \param fResolution resolution of the target map
 * \param fAdjX adjust x coordinates by this value (e.g. to uncenter) default: 0
 * \param fAdjY adjust y coordinates by this value (e.g. to uncenter) default: 0
 * \param fAdjZ adjust z coordinates by this value (e.g. to uncenter) default: 0
 */
template<class T>
svt_volume<Real64> *svt_point_cloud_pdb<T>::blur1D(Real64 fWidth, Real64 fResolution, Real64 fAdjX, Real64 fAdjY, Real64 fAdjZ, bool bProgress)
{
  // create gaussian kernel
  svt_volume<Real64> oKernel;
  oKernel.create1DBlurringKernel(fWidth, fResolution);

  // bring lattice into register with origin
  T oMin = this->getMinCoord();
  T oMax = this->getMaxCoord();

  Real64 fMinx = (fWidth * floor((oMin.x() + fAdjX) / fWidth) - fAdjX);
  Real64 fMaxx = (fWidth * ceil((oMax.x() + fAdjX) / fWidth) - fAdjX);
  Real64 fMiny = (fWidth * floor((oMin.y() + fAdjY) / fWidth) - fAdjY);
  Real64 fMaxy = (fWidth * ceil((oMax.y() + fAdjY) / fWidth) - fAdjY);
  Real64 fMinz = (fWidth * floor((oMin.z() + fAdjZ) / fWidth) - fAdjZ);
  Real64 fMaxz = (fWidth * ceil((oMax.z() + fAdjZ) / fWidth) - fAdjZ);

  // allocate protein density map
  int iMargin = (int) ceil((Real64)(oKernel.getSizeX()) / 2.0);
  int iExtx = (int)(ceil((fMaxx - fMinx) / fWidth)) + (2 * iMargin) + 1;
  int iExty = (int)(ceil((fMaxy - fMiny) / fWidth)) + (2 * iMargin) + 1;
  int iExtz = (int)(ceil((fMaxz - fMinz) / fWidth)) + (2 * iMargin) + 1;

  //SVTLBBO << "blur: target map size: " << iExtx << " x " << iExty << " x " << iExtz << endl;

  svt_volume<Real64> *pMap = new svt_volume<Real64>(iExtx, iExty, iExtz);
  //pMap->setValue( 0.0 );
  pMap->setWidth(fWidth);

  //SVTLBBO << "blur: project mass" << endl;

  // interpolate structure to protein map and keep track of variability - i.e. create a volumetric map with peaks at the atomic positions...
  Real64 fVarp = 0.0;
  Real64 fGx, fGy, fGz, fA, fB, fC;
  int iX0, iY0, iZ0, iX1, iY1, iZ1;
  unsigned int i;
  unsigned int iAtomsNum = this->size();

  try {

    for (i = 0; i < iAtomsNum; i++) {
      // compute position within grid
      fGx = iMargin + ((this->getPoint(i).x() - fMinx) / fWidth);
      fGy = iMargin + ((this->getPoint(i).y() - fMiny) / fWidth);
      fGz = iMargin + ((this->getPoint(i).z() - fMinz) / fWidth);

      iX0 = (int)(floor(fGx));
      iY0 = (int)(floor(fGy));
      iZ0 = (int)(floor(fGz));
      iX1 = iX0 + 1;
      iY1 = iY0 + 1;
      iZ1 = iZ0 + 1;

      // interpolate
      fA = (Real64)(iX1) - fGx;
      fB = (Real64)(iY1) - fGy;
      fC = (Real64)(iZ1) - fGz;

      pMap->setValue(iX0, iY0, iZ0, pMap->getValue(iX0, iY0, iZ0) + fA * fB * fC                 * m_aAtoms[i].getMass());
      fVarp += fA * fB * fC * ((1.0 - fA) * (1.0 - fA) + (1.0 - fB) * (1.0 - fB) + (1.0 - fC) * (1.0 - fC));
      pMap->setValue(iX0, iY0, iZ1, pMap->getValue(iX0, iY0, iZ1) + fA * fB * (1.0 - fC)           * m_aAtoms[i].getMass());
      fVarp += fA * fB * (1.0 - fC) * ((1.0 - fA) * (1.0 - fA) + (1.0 - fB) * (1.0 - fB) + fC * fC);
      pMap->setValue(iX0, iY1, iZ0, pMap->getValue(iX0, iY1, iZ0) + fA * (1 - fB) * fC             * m_aAtoms[i].getMass());
      fVarp += fA * (1.0 - fB) * fC * ((1.0 - fA) * (1.0 - fA) + fB * fB + (1.0 - fC) * (1.0 - fC));
      pMap->setValue(iX1, iY0, iZ0, pMap->getValue(iX1, iY0, iZ0) + (1.0 - fA) * fB * fC           * m_aAtoms[i].getMass());
      fVarp += (1.0 - fA) * fB * fC * (fA * fA + (1.0 - fB) * (1.0 - fB) + (1.0 - fC) * (1.0 - fC));
      pMap->setValue(iX0, iY1, iZ1, pMap->getValue(iX0, iY1, iZ1) + fA * (1 - fB) * (1.0 - fC)       * m_aAtoms[i].getMass());
      fVarp += fA * (1.0 - fB) * (1.0 - fC) * ((1.0 - fA) * (1.0 - fA) + fB * fB + fC * fC);
      pMap->setValue(iX1, iY1, iZ0, pMap->getValue(iX1, iY1, iZ0) + (1.0 - fA) * (1 - fB) * fC       * m_aAtoms[i].getMass());
      fVarp += (1.0 - fA) * (1.0 - fB) * fC * (fA * fA + fB * fB + (1.0 - fC) * (1.0 - fC));
      pMap->setValue(iX1, iY0, iZ1, pMap->getValue(iX1, iY0, iZ1) + (1.0 - fA) * fB * (1.0 - fC)     * m_aAtoms[i].getMass());
      fVarp += (1.0 - fA) * fB * (1.0 - fC) * (fA * fA + (1.0 - fB) * (1.0 - fB) + fC * fC);
      pMap->setValue(iX1, iY1, iZ1, pMap->getValue(iX1, iY1, iZ1) + (1.0 - fA) * (1 - fB) * (1.0 - fC) * m_aAtoms[i].getMass());
      fVarp += (1.0 - fA) * (1.0 - fB) * (1.0 - fC) * (fA * fA + fB * fB + fC * fC);

    }
    fVarp /= (double)(iAtomsNum);


  } catch (int e) {
    delete pMap;
    return NULL;
  }

  //SVTLBBO << "blur: variability: " << fVarp << endl;
  //SVTLBBO << "blur: lattice smoothing (sigma = atom rmsd): " << fWidth * sqrt(fVarp) << " Angstroem" << endl;
  //SVTLBBO << "blur: convolve with kernel" << endl;

  // convolve
  pMap->convolve1D3D(oKernel, bProgress);

  //SVTLBBO << "blur: done." << endl;

  // set correct position of map relative to pdb
  pMap->setGrid(fMinx - ((Real64)(iMargin)*fWidth) + fAdjX, fMiny - ((Real64)(iMargin)*fWidth) + fAdjY, fMinz - ((Real64)(iMargin)*fWidth) + fAdjZ);

  // return
  return pMap;
};

/**
 * Calculate rmsd between this and another structure.
 * The other structure should have equal or larger size. If larger, an oligomer is assumed,
 * and the overlap with a chain is determined. The most overlapping chain is then used as reference and the rmsd is computed.
 * \param rOlig reference to other structure
 * \param bAlign align the two structures first (default: false) - valid only for the case if the two structures have the same size, ignored otherwise!
 */
template<class T>
Real64 svt_point_cloud_pdb<T>::rmsd(svt_point_cloud_pdb<T> &rOlig, bool bAlign, Select iSelection, bool bShowProgress)
{
  Real64 fRMSD = 0.0;
  unsigned int j;

  unsigned int iChain = 0;
  bool bKearsley = bAlign;

  unsigned int iSizeThis = this->getAtomsNumber(iSelection);
  unsigned int iSizeOlig = rOlig.getAtomsNumber(iSelection);

  // Multimolecule docking or full molecule
  if (iSizeThis < iSizeOlig) {
    // and calculate the rmsd to the oligomeric subcomponent
    int aChains[256];
    for (j = 0; j < 256; j++)
      aChains[j] = 0;

    // calculate the nearest neighbor of each atom in the monomer regarding the oligomer and record the chain. The chain that occurs most is the chain that the rmsd is computed relative to.
    try {
      for (j = 0; j < this->size(); j++) {
        unsigned int iIndex = rOlig.nearestNeighbor((*this)[j]);
        if (rOlig.atom(iIndex).getChainID() - 65 >= 0 && rOlig.atom(iIndex).getChainID() - 65 < 255)
          aChains[(unsigned int)(rOlig.atom(iIndex).getChainID() - 65)]++;
        else {
          if (rOlig.atom(iIndex).getChainID() == '-')
            aChains[0]++;
          else
            SVTLBBO << "Warning: Unknown ChainID: " << rOlig.atom(iIndex).getChainID() << endl;
        }

      }
    } catch (int e) {
      return 1.0E10;
    }


    int iMax = 0;
    unsigned int iMaxIndex = 0;
    for (j = 0; j < 256; j++) {
      if (aChains[j] > iMax) {
        iMax = aChains[j];
        iMaxIndex = j;
      }
    }

    iChain = iMaxIndex * this->size();
    if (bShowProgress)
      SVTLBBO << "Compare to chain: " << iMaxIndex << endl;

    bKearsley = false;
  }
  if (iSizeThis > iSizeOlig) {
    return 1.0E10;
  }

  // least-square fit
  svt_matrix4<Real64> oMat;
  if (bKearsley)
    oMat = this->kearsley(*this, rOlig);

  // now calculate the actual rmsd
  unsigned int iSizeSelection = 0;
  try {
    unsigned int iSize = this->size();
    iSizeOlig = rOlig.size();

    switch (iSelection) {
      case ALL:
        for (j = 0; j < iSize; j++) {
          if (j + iChain < iSizeOlig) {
            fRMSD += (oMat * (*this)[j]).distanceSq(rOlig[j + iChain]);
            iSizeSelection++;
          } else
            return 1.0E10;
        }
        break;
      case BACKBONE:
        for (j = 0; j < iSize; j++) {
          if (j + iChain <= iSizeOlig) {
            if (m_aAtoms[j].isBackbone()) {
              fRMSD += (oMat * (*this)[j]).distanceSq(rOlig[j + iChain]);
              iSizeSelection++;
            }
          } else {
            return 1.0E10;
          }

        }
        break;
      case TRACE:
        for (j = 0; j < iSize; j++) {
          if (j + iChain <= iSizeOlig) {
            if (m_aAtoms[j].isCA()) {
              fRMSD += (oMat * (*this)[j]).distanceSq(rOlig[j + iChain]);
              iSizeSelection++;
            }
          } else
            return 1.0E10;

        }
        break;
    }

  } catch (int e) {
    return 1.0E10;
  }

  fRMSD /= iSizeSelection;
  fRMSD = sqrt(fRMSD);

  return fRMSD;
};

/**
 * Calculate dRMSD (distance RMSD - intramolecular distances) between this and another structure.
 * \param rOlig reference to other structure
 */
template<class T>
Real64 svt_point_cloud_pdb<T>::drmsd(svt_point_cloud_pdb<T> &rOlig)
{
  Real64 fdRMSD = 0.0;
  unsigned int iIndexAtom1, iIndexAtom2;
  Real64 fDist1, fDist2;
  unsigned int iCount = this->size() * (this->size() - 1);

  // now calculate the actual rmsd
  try {
    if (this->size() != rOlig.size()) {
      return 1.0E10;
    }


    for (iIndexAtom1 = 0; iIndexAtom1 < this->size(); iIndexAtom1++) {
      for (iIndexAtom2 = 0; iIndexAtom2 < iIndexAtom1; iIndexAtom2++) {
        fDist1 = ((*this)[iIndexAtom1]).distance(((*this)[iIndexAtom2]));
        fDist2 = (rOlig[iIndexAtom1]).distance(rOlig[iIndexAtom2]);
        fdRMSD += (fDist2 - fDist1) * (fDist2 - fDist1);
      }

    }

  } catch (int e) {
    return 1.0E10;
  }


  fdRMSD /=  iCount;
  fdRMSD = sqrt(fdRMSD);

  return fdRMSD;
};

/**
 * Get a chain
 * \param cChainID the chain ID
 * \return returns the atoms of the chain as pointer to another svt_point_cloud_pdb object
 */
template<class T>
svt_point_cloud_pdb<T> svt_point_cloud_pdb<T>::getChain(char cChainID)
{
  svt_point_cloud_pdb<T> pPDB;

  for (unsigned int i = 0; i < (*this).size(); i++)
    if ((*this).atom(i).getChainID() == cChainID)
      pPDB.addAtom((*this).atom(i), (*this)[i]);

  return pPDB;
};

/**
 * Get the trace
 * \return returns the CA as svt_point_cloud_pdb object
 */
template<class T>
svt_point_cloud_pdb<T> svt_point_cloud_pdb<T>::getTrace(unsigned int iSkip)
{
  unsigned int iSize = (*this).size();

  if (iSkip > iSize) {
    error_sba(85010, "You choose to omit more atoms than the structure has.\nPlease decrease the number of omitted atoms!");
    return *this;
  }

  svt_point_cloud_pdb<T> oPDB;
  unsigned int iSkipped = iSkip; // initialize iSkipped so that selected CA starts form the first CA

  for (unsigned int i = 0; i < iSize; i++) {
    // if this is a nucleic acid, the resname should consist of two spaces and a letter
    //
    if (
      // Atom is C alpha
      (m_aAtoms[i].isCA())

      || // or

      // Atom belongs to a nucleic acid (using P instead of CA)
      ((m_aAtoms[i].getResname()[0] == ' ' && m_aAtoms[i].getResname()[1] == ' ')
       && (m_aAtoms[i].getName()[0] == 'P')
      )
    ) {
      if (iSkipped == iSkip) {
        oPDB.addAtom(m_aAtoms[i], (*this)[i]);
        iSkipped = 0;
      } else
        iSkipped++;
    }
  }

  SVTLBBO << oPDB.size() << " atoms in the trace." << endl;

  return oPDB;
};


/**
 * gets the backbone of the structure
 * \return svt_point_cloud_pdb with the backbone
 */
template<class T>
svt_point_cloud_pdb<T> svt_point_cloud_pdb<T>::getBackbone()
{
  svt_point_cloud_pdb<T> oPDB;
  unsigned int iSize = (*this).size();
  for (unsigned int i = 0; i < iSize; i++) {
    if (m_aAtoms[i].isBackbone())
      oPDB.addAtom(m_aAtoms[i], (*this)[i]);
  }

  SVTLBBO << oPDB.size() << " atoms in the backbone." << endl;

  return oPDB;
};

/**
 * Get the number of atoms in
 * \param iSelection can be ALL, BACKBONE, or CA; with ALL- function equivalent with this->size()
 */
template<class T>
unsigned int svt_point_cloud_pdb<T>::getAtomsNumber(Select iSelection)
{
  unsigned int iNumber = 0;
  unsigned int iSize  = this->size();

  switch (iSelection) {
    case ALL:
      iNumber = iSize;
      break;
    case BACKBONE:
      for (unsigned int j = 0; j < iSize; j++) {
        if (m_aAtoms[j].isBackbone())
          iNumber++;
      }
      break;
    case TRACE:
      for (unsigned int j = 0; j < iSize; j++) {
        if (m_aAtoms[j].isCA())
          iNumber++;
      }
      break;
  }

  return iNumber;
};


///////////////////////////////////////////////////////////////////////////////
// Secondary Structure Information
///////////////////////////////////////////////////////////////////////////////

/**
 * is the secondary structure information already specified?
 * \param bSecStruct true if the atoms know their secondary structure membership
 */
template<class T>
void svt_point_cloud_pdb<T>::setSecStructAvailable(bool bSecStructAvailable)
{
  m_bSecStructAvailable = bSecStructAvailable;

  if (!m_bSecStructAvailable)
    m_eSecStructSource = SSE_NOTAVAILABLE;
  else
    m_eSecStructSource = SSE_OTHER;

};

/**
 * is the secondary structure information already specified?
 * return true if the atoms know their secondary structure membership
 */
template<class T>
bool svt_point_cloud_pdb<T>::getSecStructAvailable()
{
  return m_bSecStructAvailable;
};

/**
 * set the secondary structure source
 * \param 0 - not available, 1 - pdb, 2 - stride
 */
template<class T>
void svt_point_cloud_pdb<T>::setSecStructSource(int eSecStructSource)
{
  m_eSecStructSource = (svt_secStructSource)eSecStructSource;

  if (eSecStructSource != SSE_NOTAVAILABLE)
    m_bSecStructAvailable = true;
  else
    m_bSecStructAvailable = false;

};

/**
 * set the secondary structure source
 * \return 0 - not available, 1 - pdb, 2 - stride
 */
template<class T>
int svt_point_cloud_pdb<T>::getSecStructSource()
{
  return m_eSecStructSource;
};

/**
 * get the compressed list of secondary structure- the one used to hold the information wrote in pdb in HELIX and SHEET entry
 */
template<class T>
vector <svt_sse_pdb> &svt_point_cloud_pdb<T>::getSecStructCompressedList()
{
  return m_aSsePdb;
};

/**
 * set the compressed list of secondary structure- the one used to hold the information wrote in pdb in HELIX and SHEET entry
 */
template<class T>
void svt_point_cloud_pdb<T>::setSecStructCompressedList(vector <svt_sse_pdb> aSse)
{
  m_aSsePdb = aSse;

  //reassign the sse to the atoms
  for (unsigned int iIndex = 0; iIndex < m_aAtoms.size(); iIndex++) {
    char cSse = getSecStructFromCompressedList(m_aAtoms[iIndex]);
    m_aAtoms[iIndex].setSecStruct(cSse);
  }
};

/**
 * get the sse from the list for this atom
 * \param oAtom the atom
 */
template<class T>
char svt_point_cloud_pdb<T>::getSecStructFromCompressedList(svt_point_cloud_atom oAtom)
{
  char cSse = 'C';
  for (unsigned int iSse = 0; iSse < m_aSsePdb.size(); iSse++) {
    if (oAtom.getChainID()      == m_aSsePdb[iSse].m_cInitialResChainID &&
        oAtom.getChainID()      == m_aSsePdb[iSse].m_cTerminalResChainID &&
        oAtom.getResidueSeq()   >= m_aSsePdb[iSse].m_iInitialResSeq &&
        oAtom.getResidueSeq()   <= m_aSsePdb[iSse].m_iTerminalResSeq) {
      if ((m_aSsePdb[iSse].m_aType)[0] == 'H') {
        //consider the different types of helices
        switch (m_aSsePdb[iSse].m_iClass) {
          case 3: //Right-handed pi
            cSse = 'I';
            break;
          case 5: //Right-handed pi
            cSse = 'G';
            break;
          default:
            cSse = 'H';
        }
      } else { // is beta sheet
        cSse = 'E';
      }
    }
  }
  return cSse;
}

/**
 * calculate the secondary structure information
 */
template<class T>
void svt_point_cloud_pdb<T>::calcSecStruct()
{
  if (this->getSecStructAvailable() == true)
    return;

  if (this->m_aAtoms.size() == 0)
    return;

  //holds the current model or chain
  svt_point_cloud_pdb<T> oPDB;

  unsigned int iCurrentModel;
  char cCurrentChain;
  unsigned int iOrdResidueSeq = 0;
  int iPrevResSeq = -1;
  int iNoCA = 0;

  unsigned int iSize = this->m_aAtoms.size();

  if (iSize > 0) {
    iCurrentModel = m_aAtoms[0].getModel();
    cCurrentChain = m_aAtoms[0].getChainID();
  } else
    return;

  try {

    for (unsigned int iAtom = 0; iAtom < iSize; iAtom++) {
      oPDB.addAtom(*(*this).getAtom(iAtom), (*this)[iAtom]);

      if (oPDB.getAtom(oPDB.size() - 1)->getResidueSeq() != iPrevResSeq)
        iOrdResidueSeq++;
      iPrevResSeq = oPDB.getAtom(oPDB.size() - 1)->getResidueSeq();
      oPDB.getAtom(oPDB.size() - 1)->setOrdResidueSeq(iOrdResidueSeq);

      if (oPDB.getAtom(oPDB.size() - 1)->isCA())
        iNoCA++;

      if ((iAtom + 1 < iSize && (m_aAtoms[iAtom + 1].getModel() != iCurrentModel || m_aAtoms[iAtom + 1].getChainID() != cCurrentChain)) || iAtom == iSize - 1) { // is this the last atom in the chain,model or last atom in the pdb
        if (iNoCA >= 5) // run stride only if there are at least 5 residuals
          oPDB.runStride();

        for (unsigned int i = 0; i < oPDB.size(); i++) {
          m_aAtoms[iAtom - i].setSecStruct(oPDB.getAtom(oPDB.size() - i - 1)->getSecStruct());
          m_aAtoms[iAtom - i].setSecStructNumResidues(oPDB.getAtom(oPDB.size() - i - 1)->getSecStructNumResidues());
        }

        oPDB.deleteAllAtoms();
        oPDB.setSecStructAvailable(false);
        iNoCA = 0;

        if (iAtom + 1 < iSize) {
          iCurrentModel = m_aAtoms[iAtom + 1].getModel();
          cCurrentChain = m_aAtoms[iAtom + 1].getChainID();
          iOrdResidueSeq = 0;
        }
      }
    }


  } catch (int e) {
  }

  calcOrdinalChainIDs();
  calcSecStructLengths();
  //also set the secStructureAvailable(true);

  setSecStructSource(SSE_STRIDE);
  compressSecStruct();
};

/**
 * calculate the secondary structure information
 */
template<class T>
void svt_point_cloud_pdb<T>::runStride()
{
  /*   if ( this->m_aAtoms.size() == 0 )
         return;

     char *pBuffer;
     char aRecordName[7];

     // creates svt_stride object
     svt_stride oStride;

     // alocate memory for buffors
     vector< char *> aBuffer;
     vector< char *> aVectorRaport;

     unsigned int iAtomNumber= m_aAtoms.size();
     unsigned int i;

     char cChainID;

     for(i=0; i<iAtomNumber; i++)
     {
         cChainID = m_aAtoms[i].getChainID();
         if (cChainID == '-')
             cChainID = ' ';

         if (m_aAtoms[i].getHetAtm())
             sprintf(aRecordName,"%s","HETATM\0");
         else
             sprintf(aRecordName,"%s","ATOM  \0");

         pBuffer = new char[1024];

         sprintf(pBuffer, "%s%5i %2s%c%c%c%-2s %c%4i%c   %8.3f%8.3f%8.3f%6.2f%6.2f %3s  %4s%2s%2s\n",
                 aRecordName,
                 m_aAtoms[i].getPDBIndex(),
                 m_aAtoms[i].getName(),
                 m_aAtoms[i].getRemoteness(),
                 m_aAtoms[i].getBranch(),
                 m_aAtoms[i].getAltLoc(),
                 m_aAtoms[i].getResname(),
                 cChainID,
                 m_aAtoms[i].getResidueSeq(),
                 m_aAtoms[i].getICode(),
                 this->getPoint(i).x(),
                 this->getPoint(i).y(),
                 this->getPoint(i).z(),
                 m_aAtoms[i].getOccupancy(),
                 m_aAtoms[i].getTempFact(),
                 m_aAtoms[i].getNote(),
                 m_aAtoms[i].getSegmentID(),
                 m_aAtoms[i].getElement(),
                 m_aAtoms[i].getCharge());

         aBuffer.push_back(pBuffer);

     }

     // run stride as a class
     if(oStride.sseindexer(aBuffer, aVectorRaport))
     {
         // read in the stride output: aVectorRaport
         char pResId[6];
         int iResId;
         char cChain;
         char cSec;
         char *pCBuffer;
         map< long, char > aSecStruct;

         unsigned int iVRDim = aVectorRaport.size();
         for(i=0; i<iVRDim; i++)
         {
             pCBuffer=aVectorRaport[i];
             pResId[0] = pCBuffer[16];
             pResId[1] = pCBuffer[17];
             pResId[2] = pCBuffer[18];
             pResId[3] = pCBuffer[19];
             pResId[4] = 0;
             iResId = atoi(pResId);
             // get the chain information
             cChain = pCBuffer[9];
             // get the secondary structure information
             cSec = toupper(pCBuffer[24]);

             aSecStruct[ iResId + (cChain * 65537) ] = cSec;
         }

         // apply this secondary structure information to all atoms of the residue
         int iOrdResidueSeq=1;
         unsigned int j;

         for(j=0; j<iAtomNumber; j++)
         {
             if (j>0 && m_aAtoms[j].getOrdResidueSeq() != m_aAtoms[j-1].getOrdResidueSeq() && !m_aAtoms[j].getHetAtm())
                 iOrdResidueSeq++;

             if (m_aAtoms[j].getHetAtm())
                 m_aAtoms[j].setSecStruct('C');
             else
                 m_aAtoms[j].setSecStruct( aSecStruct[ iOrdResidueSeq + (m_aAtoms[j].getChainID() * 65537) ] );
         }

         for(i=0; i<iAtomNumber; i++)
             delete[] aBuffer[i];
         for(i=0; i<aVectorRaport.size(); i++)
             delete[] aVectorRaport[i];
     }
     else
         SVTLBBO << "Secondary structure assignment failed!" << endl;
         */
}

/**
 * Compress the secondary structure information into the ssePdb object
 * \ATTENTION: erases the existant compressed information and recompresses it from the pdb atom entry
 */
template<class T>
void svt_point_cloud_pdb<T>::compressSecStruct()
{
  unsigned int iNum = m_aAtoms.size();

  //delete the existant information
  m_aSsePdb.clear();

  char cLastSse = m_aAtoms[0].getSecStruct();
  char cCurrentSse;
  int iLastCA = 0;
  svt_sse_pdb oSse;
  sprintf(oSse.m_aType, " ");

  bool bIsSseOpen = false; //if true - a Sse was started, but not yet finished;
  int iCountHelices = 0, iCountSheets = 0;
  for (unsigned int iAtom = 1; iAtom < iNum; iAtom++) {
    if (m_aAtoms[iAtom].isCA()) {
      cCurrentSse = m_aAtoms[iAtom].getSecStruct();
      if (cCurrentSse != cLastSse) { // the secondary structure element just changed
        //close previous sse element
        if (strlen(oSse.m_aType) > 1 && bIsSseOpen) {
          bIsSseOpen = false;

          sprintf(oSse.m_aTerminalResname, "%s", m_aAtoms[iLastCA].getResname());
          oSse.m_cTerminalResChainID  = m_aAtoms[iLastCA].getChainID() ;
          oSse.m_iTerminalResSeq      = m_aAtoms[iLastCA].getResidueSeq() ;
          oSse.m_cTerminalICode       = m_aAtoms[iLastCA].getICode() ;

          if (oSse.m_aType[0] == 'H') {
            oSse.m_iLen = oSse.m_iTerminalResSeq - oSse.m_iInitialResSeq + 1;
            iCountHelices++;
            oSse.m_iNum = iCountHelices;
            sprintf(oSse.m_aID, "%3d", iCountHelices % 1000);
          }

          if (oSse.m_aType[0] == 'S') {
            iCountSheets++;
            oSse.m_iNum = iCountSheets;
            sprintf(oSse.m_aID, "%3d", iCountSheets % 1000);

            oSse.m_iNumStrands = 1;

            oSse.m_iSense = 0;

            strcpy(oSse.m_aCurAtom,  "");
            strcpy(oSse.m_aCurResname, "");
            oSse.m_cCurResChainID = ' ';
            strcpy(oSse.m_aCurResSeq, "");
            oSse.m_cCurICode = ' ';

            strcpy(oSse.m_aPrevAtom,  "");
            strcpy(oSse.m_aPrevResname, "");
            oSse.m_cPrevResChainID = ' ';
            strcpy(oSse.m_aPrevResSeq, "");
            oSse.m_cPrevICode = ' ';

          }


          m_aSsePdb.push_back(oSse);
        }

        if (cCurrentSse == 'H' || cCurrentSse == 'G' || cCurrentSse == 'I') { // got a new helix
          bIsSseOpen = true;

          sprintf(oSse.m_aType, "HELIX ");
          sprintf(oSse.m_aComment, " ");

          switch (cCurrentSse) {
            case 'H':
              oSse.m_iClass = 1;
              break;
            case 'G':
              oSse.m_iClass = 5;
              break;
            case 'I':
              oSse.m_iClass = 3;
              break;
          }

          sprintf(oSse.m_aInitialResname, "%s", m_aAtoms[iAtom].getResname());
          oSse.m_cInitialResChainID =   m_aAtoms[iAtom].getChainID() ;
          oSse.m_iInitialResSeq = m_aAtoms[iAtom].getResidueSeq() ;
          oSse.m_cInitialICode =  m_aAtoms[iAtom].getICode() ;
        }

        if (cCurrentSse == 'b' || cCurrentSse == 'B' || cCurrentSse == 'E') {
          bIsSseOpen = true;

          sprintf(oSse.m_aType, "SHEET ");

          sprintf(oSse.m_aInitialResname, "%s", m_aAtoms[iAtom].getResname());
          oSse.m_cInitialResChainID   = m_aAtoms[iAtom].getChainID() ;
          oSse.m_iInitialResSeq       = m_aAtoms[iAtom].getResidueSeq() ;
          oSse.m_cInitialICode        = m_aAtoms[iAtom].getICode() ;
        }
      }
      cLastSse = cCurrentSse;
      iLastCA = iAtom;
    }
  }
}


/**
 * set the distance cutoff value that determines the maximum distance for which consecutive residues are connected
 * \param fTraceCutoff the distance cutoff value
 */
template<class T>
void svt_point_cloud_pdb<T>::setTraceCutoff(float fTraceCutoff)
{
  m_fTraceCutoff = fTraceCutoff;
}
/**
 * get the distance cutoff value that determines the maximum distance for which consecutive residues are connected
 * \return the distance cutoff value
 */
template<class T>
float svt_point_cloud_pdb<T>::getTraceCutoff()
{
  return m_fTraceCutoff;
}
/**
 * compute the ordinal chainIDs. Those identify chains in the pdb based only on a C-alpha to C-alpha
 * distance criterion and are supposed to be calculated independently from the pdb chainIDs
 */
template<class T>
void svt_point_cloud_pdb<T>::calcOrdinalChainIDs()
{
  int i = 0, j, iOrdChainID = 1, iNumAtoms = m_aAtoms.size();
  int iAtom_first = 0, iAtom_last = 0;
  T oVecA, oVecB;


  if (iNumAtoms == 0)
    return;

  while (true) {
    if (m_aAtoms[i].isCA() || (i + 1) == iNumAtoms)
      break;
    ++i;
  }

  oVecB = T((*this)[i]);
  ++i;

  while (true) {
    if (i == iNumAtoms)
      break;

    if (m_aAtoms[i].isCA()) {
      oVecA = oVecB;
      oVecB = T((*this)[i]);

      if (oVecA.distance(oVecB) > m_fTraceCutoff) {

        // find last atom of previous residue,
        // and set ord chain id from iAtom_first to iAtom_last

        iAtom_last = i;


        // now, find the index of the last atom in the previous residue. until this one, we
        // have to set the ordinal residue sequence number
        //
        // do a lot of tests, because screwed pdb files may assign the same residue sequence
        // number to consequent residues, so it is not sufficient to test for resudue
        // sequence number alone

        while (m_aAtoms[iAtom_last].getResidueSeq() == m_aAtoms[i].getResidueSeq()
               && strcmp(m_aAtoms[iAtom_last].getResname(), m_aAtoms[i].getResname()) == 0
               && m_aAtoms[iAtom_last].getModel()      == m_aAtoms[i].getModel()
               && strcmp(m_aAtoms[iAtom_last].getSegmentID(), m_aAtoms[i].getSegmentID()) == 0) {
          --iAtom_last;
        }

        for (j = iAtom_first; j <= iAtom_last; ++j)
          m_aAtoms[j].setOrdChainID(iOrdChainID);

        // store the index of the first atom. from this one on, we will set the incremented
        // ordinal chain id for the next ordinal chain
        iAtom_first = iAtom_last + 1;
        ++iOrdChainID;

        // now go on from the last atom in the previous residue. note that with the ++i below,
        // we actually continue the iteration from the first atom of the current residue
        i = iAtom_last;
      }
    }
    ++i;
  }

  // set ordinal chainID of remaining atoms
  for (j = iAtom_first; j < iNumAtoms; ++j)
    m_aAtoms[j].setOrdChainID(iOrdChainID);

  // set ordinal chainID of all HETATMs to -1
  for (j = 0; j < iNumAtoms; ++j) {
    if (m_aAtoms[j].getHetAtm())
      m_aAtoms[j].setOrdChainID(-1);

    //printf("atom %6i ordResSeq %5i ordchainID %3i chainID %c\n", j, m_aAtoms[j].getOrdResidueSeq(), m_aAtoms[j].getOrdChainID(), m_aAtoms[j].getChainID());
  }
}
/**
 * compute and set the length of the secondary structure elements. Secondary structure must have
 * been computed (obviously...), and calcOrdinalChainIDs must have been called before for
 * calcSecStructLengths to make sense
 */
template<class T>
void svt_point_cloud_pdb<T>::calcSecStructLengths()
{
  //if PDB file is corrupt or empty meaning no atoms, do nothing, return.
  if (m_aAtoms.size() == 0)
    return;
  // compute and set the length of the secondary structure elements
  //
  char cLastStruct = m_aAtoms[0].getSecStruct(), cCurrentStruct;
  int iNumAtoms = m_aAtoms.size();
  int iFirstAtom = 0, iNumResidues = 0;
  int i, j;
  bool bSetLength = false;


  // if the first atom in the array is a C_alpha, count it
  //
  if (m_aAtoms[0].isCA())
    ++iNumResidues;


  // only one atom? nothing much to do
  //
  if (iNumAtoms == 1) {
    if (m_aAtoms[0].getSecStruct() != ' ')
      m_aAtoms[0].setSecStructNumResidues(0);
    else
      m_aAtoms[0].setSecStructNumResidues(iNumResidues);

    return;
  }


  // assign the length of secondary structure (i.e. the number of residues in it) to all atoms
  // in a secondary structure section. this length is defined as the number of residues, which
  // equals the number of C_alphas.
  //
  // a sec. structure part ends, if
  //   1. the sec. structure type changes,
  //   2. no sec. structure is assigned to the next atom (like for HETATMs)
  //   3. a new ordinal chain begins (calcOrdinalChainIDs must have been called before)
  //

  i = 1;

  while (true) {

    // at every C_alpha, increase the residue count (since there is one C_alpha per residue)
    //if (m_aAtoms[i].isCA())
    if (m_aAtoms[i].isCA())
      ++iNumResidues;


    cCurrentStruct = m_aAtoms[i].getSecStruct();


    // 1. does the sec. structure type change, or
    // 2. is it empty?
    switch (cCurrentStruct) {
      case 'H':
      case 'G':
      case 'I':
        if (cLastStruct != 'H' && cLastStruct != 'G' && cLastStruct != 'I')
          bSetLength = true;
        break;

      case 'E':
      case 'B':
        if (cLastStruct != 'E' && cLastStruct != 'B')
          bSetLength = true;
        break;

      case 'T':
      case 'C':
        if (cLastStruct != 'T' && cLastStruct != 'C')
          bSetLength = true;
        break;

      default://case ' ':
        bSetLength = true;
        break;
    }


    // 3. does a new (ordinal) chain begin?
    if (m_aAtoms[i - 1].getOrdChainID() != m_aAtoms[i].getOrdChainID())
      bSetLength = true;


    // if it was determined that the secondary structure element ends, or we arrived at the last
    // atom, setSecStructNumResidues(iNumResidues) must be called now for all atoms in the current
    // structure element (i.e., since m_aAtoms[iFirstAtom])
    if (bSetLength || i == iNumAtoms - 1) {
      if (i == iNumAtoms - 1)
        ++i;
      // the following else if is because we are already at the C_alpha of the next part, and
      // it was already added to iNumResidues, which is not correct
      else if (m_aAtoms[i].isCA())
        --iNumResidues;


      for (j = iFirstAtom; j < i; ++j)
        m_aAtoms[j].setSecStructNumResidues(iNumResidues);

      bSetLength = false;
      iFirstAtom = i;

      // if the current atom is a C_alpha, it needs to be counted,too
      if (m_aAtoms[i].isCA())
        iNumResidues = 1;
      else
        iNumResidues = 0;
    }


    // for debugging...
    //    printf("atom index: % 6i", i);
    //    printf("  sec struct: %c", m_aAtoms[i].getSecStruct());
    //    printf("  ord res: % 4i", m_aAtoms[i].getOrdResidueSeq());
    //    printf("  ord chainID: % 3i", m_aAtoms[i].getOrdChainID());
    //    printf("  iNumResidues: % 3i ", iNumResidues);
    //  printf("\n");


    cLastStruct = cCurrentStruct;

    if (i == iNumAtoms)
      break;
    else
      ++i;
  }

  // for debugging...
  //     for (i=0;i<iNumAtoms; ++i)
  //     {
  //    printf("atom index: % 6i", i);
  //    printf("  sec struct: %c", m_aAtoms[i].getSecStruct());
  //    printf("  ord res: % 4i", m_aAtoms[i].getOrdResidueSeq());
  //    printf("  ord chainID: % 3i", m_aAtoms[i].getOrdChainID());
  //    printf("  iNumResidues: % 3i", m_aAtoms[i].getSecStructNumResidues());
  //  printf("\n");
  //     }
}
/**
 * calculate the bonds between the atoms (according to their distance the bonds are guessed)
 */
template<class T>
void svt_point_cloud_pdb<T>::calcBonds(bool bShowProgress)
{
  if (m_aAtoms.size() == 0)
    return;

  SVTLBBO << "No information about covalent bonds found, therefore the bonds are guessed based on a distance criterion..." << endl;

  // cutoff is maximum vdw radius * 1.2
  Real32 fCutoff = 1.90f * 2.5f;
  Real32 fCutoffSq = fCutoff * fCutoff;

  deleteAllBonds();

  SVTLBBO << "  Space Partitioning..." << endl;

  // calculate the minimum and maximum atom coords
  Real32 fMinX = this->getMinXCoord();
  Real32 fMaxX = this->getMaxXCoord();
  Real32 fMinY = this->getMinYCoord();
  Real32 fMaxY = this->getMaxYCoord();
  Real32 fMinZ = this->getMinZCoord();
  Real32 fMaxZ = this->getMaxZCoord();

  SVTLBBO << "  Max: ( " << fMaxX << ", " << fMaxY << ", " << fMaxZ << " ) - Min: ( " << fMinX << ", " << fMinY << ", " << fMinZ << " )" << endl;

  // pre-sorting all atoms into boxes
  // how many boxes do we need?
  unsigned int iBoxesX = (unsigned int)(((fMaxX - fMinX) / fCutoff) + 0.5f) + 1;
  unsigned int iBoxesY = (unsigned int)(((fMaxY - fMinY) / fCutoff) + 0.5f) + 1;
  unsigned int iBoxesZ = (unsigned int)(((fMaxZ - fMinZ) / fCutoff) + 0.5f) + 1;
  unsigned int iBoxesXY = iBoxesX * iBoxesY;

  SVTLBBO << "  " << iBoxesX << " x " << iBoxesY << " x " << iBoxesZ << endl;

  // allocate memory
  vector< vector<unsigned int> > aBox((iBoxesX + 1) * (iBoxesY + 1) * (iBoxesZ + 1));
  int iBoxX, iBoxY, iBoxZ, iBoxIndex;
  // now sort all atoms into these boxes
  unsigned int i;
  for (i = 0; i < m_aAtoms.size(); i++) {
    iBoxX = (int)((((*this)[i].x() - fMinX) / fCutoff) + 0.5f);
    iBoxY = (int)((((*this)[i].y() - fMinY) / fCutoff) + 0.5f);
    iBoxZ = (int)((((*this)[i].z() - fMinZ) / fCutoff) + 0.5f);
    iBoxIndex = iBoxX + (iBoxY * iBoxesX) + (iBoxZ * iBoxesXY);
    aBox[iBoxIndex].push_back(i);
  }

  // now look for neighbors in these cubes
  vector< vector<unsigned int > > aNeighbor(m_aAtoms.size());
  unsigned int iIndex, iIndexNext, iAtom, iNumElements;
  unsigned int j, x, y, z;

  try {
    for (z = 0; z < iBoxesZ; z++) {
      for (y = 0; y < iBoxesY; y++)
        for (x = 0; x < iBoxesX; x++) {
          // calculate the index of the current block
          iIndex = x + (y * iBoxesX) + (z * iBoxesXY);
          // loop over all atoms in that block
          for (i = 0; i < aBox[iIndex].size(); i++) {
            iAtom = aBox[iIndex][i];

            // search in same block for neighbors
            for (j = i + 1; j < aBox[iIndex].size(); j++) {
              if ((*this)[iAtom].distanceSq((*this)[aBox[iIndex][j]]) > fCutoffSq)
                continue;
              aNeighbor[iAtom].push_back(aBox[iIndex][j]);
            }
            // search for neighbors in block x+1
            if (x < iBoxesX - 1) {
              iIndexNext = (x + 1) + (y * iBoxesX) + (z * iBoxesXY);
              iNumElements = aBox[iIndexNext].size();
              for (j = 0; j < iNumElements; j++) {
                if ((*this)[iAtom].distanceSq((*this)[aBox[iIndexNext][j]]) > fCutoffSq)
                  continue;
                aNeighbor[iAtom].push_back(aBox[iIndexNext][j]);
              }
            }
            // search for neighbors in block y+1
            if (y < iBoxesY - 1) {
              iIndexNext = x + ((y + 1) * iBoxesX) + (z * iBoxesXY);
              iNumElements = aBox[iIndexNext].size();
              for (j = 0; j < iNumElements; j++) {
                if ((*this)[iAtom].distanceSq((*this)[aBox[iIndexNext][j]]) > fCutoffSq)
                  continue;
                aNeighbor[iAtom].push_back(aBox[iIndexNext][j]);
              }
            }
            // search for neighbors in block z+1 block
            if (z < iBoxesZ - 1) {
              iIndexNext = x + (y * iBoxesX) + ((z + 1) * iBoxesXY);
              iNumElements = aBox[iIndexNext].size();
              for (j = 0; j < iNumElements; j++) {
                if ((*this)[iAtom].distanceSq((*this)[aBox[iIndexNext][j]]) > fCutoffSq)
                  continue;
                aNeighbor[iAtom].push_back(aBox[iIndexNext][j]);
              }
            }
            // search for neighbors in block x+1 and y+1
            if (x < iBoxesX - 1 && y < iBoxesY - 1) {
              iIndexNext = (x + 1) + ((y + 1) * iBoxesX) + (z * iBoxesXY);
              iNumElements = aBox[iIndexNext].size();
              for (j = 0; j < iNumElements; j++) {
                if ((*this)[iAtom].distanceSq((*this)[aBox[iIndexNext][j]]) > fCutoffSq)
                  continue;
                aNeighbor[iAtom].push_back(aBox[iIndexNext][j]);
              }
            }
            // search for neighbors in block x+1 and z+1
            if (x < iBoxesX - 1 && z < iBoxesZ - 1) {
              iIndexNext = (x + 1) + (y * iBoxesX) + ((z + 1) * iBoxesXY);
              iNumElements = aBox[iIndexNext].size();
              for (j = 0; j < iNumElements; j++) {
                if ((*this)[iAtom].distanceSq((*this)[aBox[iIndexNext][j]]) > fCutoffSq)
                  continue;
                aNeighbor[iAtom].push_back(aBox[iIndexNext][j]);
              }
            }
            // search for neighbors in block y+1 and z+1
            if (y < iBoxesY - 1 && z < iBoxesZ - 1) {
              iIndexNext = x + ((y + 1) * iBoxesX) + ((z + 1) * iBoxesXY);
              iNumElements = aBox[iIndexNext].size();
              for (j = 0; j < iNumElements; j++) {
                if ((*this)[iAtom].distanceSq((*this)[aBox[iIndexNext][j]]) > fCutoffSq)
                  continue;
                aNeighbor[iAtom].push_back(aBox[iIndexNext][j]);
              }
            }
            // search for neighbors in block x+1, y-1
            if (x < iBoxesX - 1 && y > 0) {
              iIndexNext = (x + 1) + ((y - 1) * iBoxesX) + (z * iBoxesXY);
              iNumElements = aBox[iIndexNext].size();
              for (j = 0; j < iNumElements; j++) {
                if ((*this)[iAtom].distanceSq((*this)[aBox[iIndexNext][j]]) > fCutoffSq)
                  continue;
                aNeighbor[iAtom].push_back(aBox[iIndexNext][j]);
              }
            }

            // search for neighbors in block x-1, z+1
            if (x > 0 && z < iBoxesZ - 1) {
              iIndexNext = (x - 1) + (y * iBoxesX) + ((z + 1) * iBoxesXY);
              iNumElements = aBox[iIndexNext].size();
              for (j = 0; j < iNumElements; j++) {
                if ((*this)[iAtom].distanceSq((*this)[aBox[iIndexNext][j]]) > fCutoffSq)
                  continue;
                aNeighbor[iAtom].push_back(aBox[iIndexNext][j]);
              }
            }
            // search for neighbors in block z+1, y-1
            if (y > 0 && z < iBoxesZ - 1) {
              iIndexNext = x + ((y - 1) * iBoxesX) + ((z + 1) * iBoxesXY);
              iNumElements = aBox[iIndexNext].size();
              for (j = 0; j < iNumElements; j++) {
                if ((*this)[iAtom].distanceSq((*this)[aBox[iIndexNext][j]]) > fCutoffSq)
                  continue;
                aNeighbor[iAtom].push_back(aBox[iIndexNext][j]);
              }
            }
            // search for neighbors in block z+1, y+1, x+1 block
            if (x < iBoxesX - 1 && y < iBoxesY - 1 && z < iBoxesZ - 1) {
              iIndexNext = (x + 1) + ((y + 1) * iBoxesX) + ((z + 1) * iBoxesXY);
              iNumElements = aBox[iIndexNext].size();
              for (j = 0; j < iNumElements; j++) {
                if ((*this)[iAtom].distanceSq((*this)[aBox[iIndexNext][j]]) > fCutoffSq)
                  continue;
                aNeighbor[iAtom].push_back(aBox[iIndexNext][j]);
              }
            }
            // search for neighbors in block x-1, y+1, z+1
            if (x > 0 && y < iBoxesY - 1 && z < iBoxesZ - 1) {
              iIndexNext = (x - 1) + ((y + 1) * iBoxesX) + ((z + 1) * iBoxesXY);
              iNumElements = aBox[iIndexNext].size();
              for (j = 0; j < iNumElements; j++) {
                if ((*this)[iAtom].distanceSq((*this)[aBox[iIndexNext][j]]) > fCutoffSq)
                  continue;
                aNeighbor[iAtom].push_back(aBox[iIndexNext][j]);
              }
            }
            // search for neighbors in block x+1, y-1, z+1
            if (x < iBoxesX - 1 && y > 0 && z < iBoxesZ - 1) {
              iIndexNext = (x + 1) + ((y - 1) * iBoxesX) + ((z + 1) * iBoxesXY);
              iNumElements = aBox[iIndexNext].size();
              for (j = 0; j < iNumElements; j++) {
                if ((*this)[iAtom].distanceSq((*this)[aBox[iIndexNext][j]]) > fCutoffSq)
                  continue;
                aNeighbor[iAtom].push_back(aBox[iIndexNext][j]);
              }
            }
            // search for neighbors in block x-1, y-1, z+1
            if (x > 0 && y > 0 && z < iBoxesZ - 1) {
              iIndexNext = (x - 1) + ((y - 1) * iBoxesX) + ((z + 1) * iBoxesXY);
              iNumElements = aBox[iIndexNext].size();
              for (j = 0; j < iNumElements; j++) {
                if ((*this)[iAtom].distanceSq((*this)[aBox[iIndexNext][j]]) > fCutoffSq)
                  continue;
                aNeighbor[iAtom].push_back(aBox[iIndexNext][j]);
              }
            }
          }

        }

    }

  } catch (int e) {
  }

  try {
    // now we have an neighbor array with all corresponding atom pairs within a certain distance. Now we must check which of these pairs really are bonds.
    Real32 fDist, fCut;
    int k, iNeighborNum;
    for (i = 0; i < m_aAtoms.size(); i++) {
      iNeighborNum = aNeighbor[i].size();

      for (k = 0; k < iNeighborNum; k++) {
        j = aNeighbor[i][k];

        if (m_aAtoms[i].getModel() == m_aAtoms[j].getModel()) {
          fDist = (*this)[i].distance((*this)[j]);
          fCut = (m_aAtoms[j].getVDWRadius() + m_aAtoms[i].getVDWRadius()) * 0.6f;

          if ((!m_aAtoms[i].isHydrogen() || !m_aAtoms[j].isHydrogen()) && (!m_aAtoms[i].isQPDB() && !m_aAtoms[j].isQPDB())) {
            if (fDist < fCut)
              addBond(i, j);
          }
        }
      }

    }

  } catch (int e) { }


  SVTLBBO << "Bond calculation found " << m_aBonds.size() << " covalent bonds" << endl;
}


/**
 * Calculates and returns the Carbon Alpha (CA) of the atom's residue
 * \param iIndexAtom an integer representing the atom's index in the pdb structure
 * \return an svt_point_cloud_atom - the Carbon Alpha
 */
template<class T>
unsigned int svt_point_cloud_pdb<T>::getCA(unsigned int iIndexAtom)
{
  /*
  cout << "Atom " << iIndexAtom << " " << m_aAtoms[iIndexAtom].getName();
  cout << " :Res " << m_aAtoms[iIndexAtom].getResname() << " " << m_aAtoms[iIndexAtom].getResidueSeq();
  cout << " IsCa:" << m_aAtoms[iIndexAtom].isCA() <<" ";
  */

  if (m_aAtoms[iIndexAtom].isCA()) return iIndexAtom;

  int iIndex = iIndexAtom - 1;
  while ((iIndex >= 0) && (m_aAtoms[iIndex].getResidueSeq() == m_aAtoms[iIndexAtom].getResidueSeq()) &&
         (m_aAtoms[iIndex].isCA() == 0)) {
    iIndex--;
  }
  if (iIndex >= 0 && m_aAtoms[iIndex].getResidueSeq() == m_aAtoms[iIndexAtom].getResidueSeq() && m_aAtoms[iIndex].isCA() == 1) return iIndex;

  iIndex = iIndexAtom;
  while (iIndex < (int)(m_aAtoms.size()) &&
         m_aAtoms[iIndex].getResidueSeq() == m_aAtoms[iIndexAtom].getResidueSeq() &&
         m_aAtoms[iIndex].isCA() == 0) {
    iIndex++;
  }

  if (iIndex < (int)(m_aAtoms.size()) && m_aAtoms[iIndex].getResidueSeq() == m_aAtoms[iIndexAtom].getResidueSeq() && m_aAtoms[iIndex].isCA() == 1) return iIndex;

  error_sba(85010, "Error: Residue without CA");
  return iIndexAtom;
}

/**
 * Determine it iIndexAtom is on the backbone
 * \param iIndexAtom an integer representing the atom's index in the pdb structure
 * \return unsigned int - iIndexAtom if the atom is on the backbone or the CA if not
 */
template<class T>
unsigned int svt_point_cloud_pdb<T>::getBackbone(unsigned int iIndexAtom)
{
  if (m_aAtoms[iIndexAtom].isBackbone())
    return iIndexAtom;
  else
    return getCA(iIndexAtom);

}

/**
 * Compute the Distances between any atom of the structure using the topology - graph connections given in the PSF file
 */
template<class T>
bool svt_point_cloud_pdb<T>::computePSFGraphDistMat()
{
  // check if the psf file has been read
  if (m_aBonds.size() == 0) {
    error_sba(85010, "Can not compute bonds graph matrix because no bonds/connections have been created/loaded!");
    return false;
  }

  //scale matrix
  m_oGraphDists.resize(m_aAtoms.size(), m_aAtoms.size());

  // init the matrix
  for (unsigned int iIndexAtom1 = 0; iIndexAtom1 < m_aAtoms.size(); iIndexAtom1++)
    for (unsigned int iIndexAtom2 = 0; iIndexAtom2 < m_aAtoms.size(); iIndexAtom2++)
      m_oGraphDists[iIndexAtom1][iIndexAtom2] = (iIndexAtom1 == iIndexAtom2) ? 0.0 : 1000000000;//numeric_limits<Real64>::max();

  Real64 fDistTmp;
  for (unsigned int iIndexBond = 0; iIndexBond < m_aBonds.size(); iIndexBond++) {
    // in psf the atoms are numbered from 1 therefor getIndex -1
    fDistTmp = (*this)[m_aBonds[iIndexBond].getIndexA()].distance((*this)[m_aBonds[iIndexBond].getIndexB()]);

    m_oGraphDists[m_aBonds[iIndexBond].getIndexA()][m_aBonds[iIndexBond].getIndexB()] = fDistTmp;
    m_oGraphDists[m_aBonds[iIndexBond].getIndexB()][m_aBonds[iIndexBond].getIndexA()] = fDistTmp;
  }

  //Floyd's Algorithm
  for (unsigned int iIndexAtom3 = 0; iIndexAtom3 < m_aAtoms.size(); iIndexAtom3++)
    for (unsigned int iIndexAtom1 = 0; iIndexAtom1 < m_aAtoms.size(); iIndexAtom1++)
      for (unsigned int iIndexAtom2 = 0; iIndexAtom2 < m_aAtoms.size(); iIndexAtom2++)
        if (m_oGraphDists[iIndexAtom1][iIndexAtom2] > m_oGraphDists[iIndexAtom1][iIndexAtom3] + m_oGraphDists[iIndexAtom3][iIndexAtom2])  m_oGraphDists[iIndexAtom1][iIndexAtom2]  =  m_oGraphDists[iIndexAtom1][iIndexAtom3] + m_oGraphDists[iIndexAtom3][iIndexAtom2];

  return true;
}

/**
 * \param bPSFRead a bool indicating if the psf files has been read;
 */
template<class T>
void svt_point_cloud_pdb<T>::setPSFRead(bool bPSFRead)
{
  m_bPSFRead = bPSFRead;
}

template <class T>
void svt_point_cloud_pdb<T>::selectAtomResidueSeq(int iResidueSeq)
{

  (*this).resetAtomEnum();
  int iAtom = (*this).enumAtomResidueSeq(iResidueSeq);

  while (iAtom != -1) {
    (*this).getAtom(iAtom)->setSelected(true);
    iAtom = (*this).enumAtomResidueSeq(iResidueSeq);
  }

};

/**
 * get the index of the selected/deselected atoms
 * \param bIsSelected the selection status = true - is selected; false  - is not selected
 * \return the index of the atoms with selected status indicated by bIsSelected
 */
template <class T>
vector<unsigned int> svt_point_cloud_pdb<T>::getSelection(bool bIsSelected)
{
  vector<unsigned int> aVec;
  for (unsigned int iIndex = 0; iIndex < (*this).size(); iIndex++) {
    if ((*this).getAtom(iIndex)->getSelected() == bIsSelected)
      aVec.push_back(iIndex);
  }
  return aVec;
};


/**
 * create symmetric oligomers
 * \param symmetry type
 * \param symmetry axis
 * \param order is equivalent with the number of monomers in the symmetric unit
 * \param fOffsetAxis1 is the offset form the first axis
 * \param fOffsetAxis2 is the offset form the second axis
 * \return pdb containing the symmetric oligomere
 */
template<class T>
svt_point_cloud_pdb<T> svt_point_cloud_pdb<T>::applySymmetry(unsigned int iOrder, const SymmetryType eType, char cAxis, Real64 fOffsetAxis1, Real64 fOffsetAxis2)
{
  svt_point_cloud_pdb<T> oPDB, oPDBMono;
  svt_point_cloud_atom *pAtom;
  Real64 fTheta, fCos, fSin;
  svt_vector4<Real64> oVec;
  svt_matrix4<Real64> oMat, oRot;
  char cChain[2];


  switch (eType) {
    case SYMMETRY_C:

      fTheta = 360.0 / (Real64)iOrder;
      for (unsigned int iIndex = 0; iIndex < iOrder; iIndex++) {
        oPDBMono.deleteAllAtoms();
        oPDBMono.deleteAllBonds();

        fCos = cos((Real64)iIndex * PI * fTheta / 180.0);
        fSin = sin((Real64)iIndex * PI * fTheta / 180.0);

        switch (cAxis) {
          case 'x':
            for (unsigned int iIndexAtom = 0; iIndexAtom < (*this).size(); iIndexAtom++) {
              oVec.x((*this)[iIndexAtom].x());
              oVec.y(fOffsetAxis1 + fCos * ((*this)[iIndexAtom].y() - fOffsetAxis1) - fSin * ((*this)[iIndexAtom].z() - fOffsetAxis2));
              oVec.z(fOffsetAxis2 + fSin * ((*this)[iIndexAtom].y() - fOffsetAxis1) + fCos * ((*this)[iIndexAtom].z() - fOffsetAxis2));

              oPDBMono.addAtom(*(*this).getAtom(iIndexAtom), oVec);

              pAtom = oPDBMono.getAtom(oPDBMono.size() - 1);
              sprintf(cChain, "%d", iIndex);
              pAtom->setChainID(cChain[0]);
              //pAtom->setPDBIndex(iIndex*(*this).size()+ iIndexAtom);
            }
            break;
          case 'y':
            for (unsigned int iIndexAtom = 0; iIndexAtom < (*this).size(); iIndexAtom++) {
              oVec.x(fOffsetAxis1 + fCos * ((*this)[iIndexAtom].x() - fOffsetAxis1) + fSin * ((*this)[iIndexAtom].z() - fOffsetAxis2));
              oVec.y((*this)[iIndexAtom].y());
              oVec.z(fOffsetAxis2 - fSin * ((*this)[iIndexAtom].x() - fOffsetAxis1) + fCos * ((*this)[iIndexAtom].z() - fOffsetAxis2));

              oPDBMono.addAtom(*(*this).getAtom(iIndexAtom), oVec);

              pAtom = oPDBMono.getAtom(oPDBMono.size() - 1);
              sprintf(cChain, "%d", iIndex);
              pAtom->setChainID(cChain[0]);
              //pAtom->setPDBIndex(iIndex*(*this).size()+ iIndexAtom);
            }
            break;
          case 'z':
            for (unsigned int iIndexAtom = 0; iIndexAtom < (*this).size(); iIndexAtom++) {
              oVec.x(fOffsetAxis1 + fCos * ((*this)[iIndexAtom].x() - fOffsetAxis1) - fSin * ((*this)[iIndexAtom].y() - fOffsetAxis2));
              oVec.y(fOffsetAxis2 + fSin * ((*this)[iIndexAtom].x() - fOffsetAxis1) + fCos * ((*this)[iIndexAtom].y() - fOffsetAxis2));
              oVec.z((*this)[iIndexAtom].z());

              oPDBMono.addAtom(*(*this).getAtom(iIndexAtom), oVec);

              pAtom = oPDBMono.getAtom(oPDBMono.size() - 1);
              sprintf(cChain, "%d", iIndex);
              pAtom->setChainID(cChain[0]);
              //pAtom->setPDBIndex(iIndex*(*this).size()+ iIndexAtom);
            }
            break;
        }
        oPDB.append(oPDBMono);
      }

      break;
    case SYMMETRY_D:
    case SYMMETRY_H:
    default:
      SVTLBBO << "Symmetry type not yet available!" << endl;
      return (*this);
  }

  //recompute the models
  oPDB.calcAtomModels();

  return oPDB;
};

///////////////////////////////////////////////////////////////////////////////
// svt_vector3
///////////////////////////////////////////////////////////////////////////////

#if !defined(NOLIMITS)
#include <limits>
#else
#include <float.h>
#endif


typedef Real32 Point3f[3];
typedef Real32 Vector3f[3];


//
// binary operators
//

// Sum / Difference / Product
template<class T>
inline svt_vector3<T> operator+(const svt_vector3<T> &p1, const svt_vector3<T> &p2);

template<class T>
inline svt_vector3<T> operator+(const svt_vector3<T> &p, const T &f);

template<class T>
inline svt_vector3<T> operator+(const T &f, const svt_vector3<T> &p);


template<class T>
inline svt_vector3<T> operator-(const svt_vector3<T> &p1, const svt_vector3<T> &p2);

template<class T>
inline svt_vector3<T> operator-(const svt_vector3<T> &p, const T &f);

template<class T>
inline svt_vector3<T> operator-(const T &f, const svt_vector3<T> &p);

template<class T>
inline svt_vector3<T> operator-(const svt_vector3<T> &p);


template<class T>
inline svt_vector3<T> operator*(const svt_vector3<T> &p, const T &f);

template<class T>
inline svt_vector3<T> operator*(const T &f, const svt_vector3<T> &p);

template<class T>
inline svt_vector3<T> operator/(const svt_vector3<T> &p1, const T &f);

// Scalar Product
template<class T>
inline T operator*(const svt_vector3<T> &p1, const svt_vector3<T> &p2);


// Vector Product
template<class T>
inline svt_vector3<T>
vectorProduct(const svt_vector3<T> &p1, const svt_vector3<T> &p2);

template<class T>
inline svt_vector3<T>
crossProduct(const svt_vector3<T> &p1, const svt_vector3<T> &p2);


///////////////////////////////////////////////////////////////////////////////
// declaration
///////////////////////////////////////////////////////////////////////////////



/** A 3 value template vector
  *@author Maik Boltes, Stefan Birmanns, Frank Delonge
  */
template<class T> class svt_vector3 : public svt_matrix<T>
{
  private:

    T stack_data[3];

  public:

    /**
     * Constructor
     * \param fX initial x coordinate
     * \param fY initial y coordinate
     * \param fZ initial z coordinate
     * \param fW initial w coordinate
     */
    svt_vector3(T fX, T fY, T fZ)
      : svt_matrix<T>(3, 1, stack_data)
    {
      x(fX);
      y(fY);
      z(fZ);
    }

    svt_vector3(T fValue = T(0))
      : svt_matrix<T>(3, 1, stack_data)
    {
      x(fValue);
      y(fValue);
      z(fValue);
    }

    svt_vector3(const Point3f &rVec)
      : svt_matrix<T>(3, 1, stack_data)
    {
      x(rVec[0]);
      y(rVec[1]);
      z(rVec[2]);
    }

    svt_vector3(svt_vector4<T> oVec4)
      : svt_matrix<T>(3, 1, stack_data)
    {
      x(oVec4[0]);
      y(oVec4[1]);
      z(oVec4[2]);
    }

    svt_vector3(const svt_matrix<T> &that) : svt_matrix<T>(that, stack_data) {}

    virtual ~svt_vector3() {}

    svt_vector3<T> &operator=(const svt_vector3<T> &that);

    svt_vector3<T> &operator=(const Point3f &that);

    T &operator[](unsigned i)
    {
      return svt_matrix<T>::m_pData[i];
    }

    const T &operator[](unsigned i) const
    {
      return svt_matrix<T>::m_pData[i];
    }

    //
    // arithmetic operators
    // These operators need to be redefined because
    // only the first 3 coordinates shall be considered
    // all operators return *this to allow daisy chaining
    //
    svt_vector3<T> &operator+=(const svt_vector3<T> &p)
    {
      (*this)[0] += p[0];
      (*this)[1] += p[1];
      (*this)[2] += p[2];
      return *this;
    }

    svt_vector3<T> &operator+=(const T &f)
    {
      (*this)[0] += f;
      (*this)[1] += f;
      (*this)[2] += f;
      return *this;
    }

    svt_vector3<T> &operator-=(const svt_vector3<T> &p)
    {
      (*this)[0] -= p[0];
      (*this)[1] -= p[1];
      (*this)[2] -= p[2];
      return *this;
    }

    svt_vector3<T> &operator-=(const T &f)
    {
      (*this)[0] -= f;
      (*this)[1] -= f;
      (*this)[2] -= f;
      return *this;
    }

    svt_vector3<T> &operator*=(const T &f)
    {
      (*this)[0] *= f;
      (*this)[1] *= f;
      (*this)[2] *= f;
      return *this;
    }

    svt_vector3<T> &operator/=(const T &f)
    {
      (*this)[0] /= f;
      (*this)[1] /= f;
      (*this)[2] /= f;
      return *this;
    }

    /**
     * get/set methods
     */
    T x() const
    {
      return (*this)[0];
    }
    void x(T value)
    {
      (*this)[0] = value;
    }

    T y() const
    {
      return (*this)[1];
    }
    void y(T value)
    {
      (*this)[1] = value;
    }

    T z() const
    {
      return (*this)[2];
    }
    void z(T value)
    {
      (*this)[2] = value;
    }

    /**
     * set all three coords of the vector at once
     * \param fX x coord
     * \param fY y coord
     * \param fZ z coord
     */
    void set(T fX, T fY, T fZ)
    {
      x(fX);
      y(fY);
      z(fZ);
    }

    void set(T value)
    {
      *this = value;
    }

    void set(const T *p)
    {
      memcpy(svt_matrix<T>::m_pData, p, 3 * sizeof(T));
    }

    /**
     * set from a Vector4f
     * \param aNewVec new Vector4f to set the svt_vector4 to
     */
    void set(const Point3f &aNewVec)
    {
      x(aNewVec[0]);
      y(aNewVec[1]);
      z(aNewVec[2]);
    }


    /**
     * get the squared length of the vector
     * \return length^2
     */
    T lengthSq() const
    {
      return x() * x() + y() * y() + z() * z();
    }

    /**
     * get the length of the vector
     * \return length
     */
    T length() const
    {
      return sqrt(lengthSq());
    }

    /**
     * calculate the distance between two vectors
     * \param oVec the other vector
     * \return distance
     */
    T distance(const svt_vector3<T> &oVec) const
    {
      return ((*this - oVec)).length();
    }

    /**
     * calculate the squared distance between two vectors
     * \param oVec the other vector
     * \return squared distance
     */
    T distanceSq(const svt_vector3<T> &oVec) const
    {
      return (*this - oVec).lengthSq();
    }

    /**
     * calculate the squared distance between vectors and a vector from type Vector3f
     * \param oVec the other vector
     * \return squared distance
     */
    T distanceSq(const Vector3f &oVec) const
    {
      return (x() - oVec[0]) * (x() - oVec[0]) + (y() - oVec[1]) * (y() - oVec[1]) + (z() - oVec[2]) * (z() - oVec[2]);
    }

    /**
     * check, if vectors are nearly equal using the fast maximum norm
     * \param oVec the other vector
     * \return true if nearly equal otherwise false
     */
    bool nearEqual(const Vector3f &oVec)
    {
      if ((fabs(x() - oVec[0]) < EQ_EPS) && (fabs(y() - oVec[1]) < EQ_EPS) && (fabs(z() - oVec[2]) < EQ_EPS))
        return true;
      else
        return false;
    }

    /**
     * check, if vectors are nearly equal using the fast maximum norm
     * \param oVec the other vector
     * \return true if nearly equal otherwise false
     */
    bool nearEqual(const svt_vector3<T> &oVec)
    {
      if ((fabs(x() - oVec[0]) < EQ_EPS) && (fabs(y() - oVec[1]) < EQ_EPS) && (fabs(z() - oVec[2]) < EQ_EPS))
        return true;
      else
        return false;
    }

    /**
     * normalize the vector, return *this to allow daisy chaining
     */
    svt_vector3<T> &normalize()
    {
      T fLength = length();
      // 0 only for (0,0,0), so we have not to divide
      if (isPositive(fLength))
        (*this) /= fLength;

      return *this;
    }

    /**
     * returns data pointer for old function calls
     */
    T *cast()
    {
      return svt_matrix<T>::m_pData;//why not stack_data?
    }


};


///////////////////////////////////////////////////////////////////////////////
// definition
///////////////////////////////////////////////////////////////////////////////

template<class T>
svt_vector3<T> &svt_vector3<T>::operator=(const svt_vector3<T> &that)
{
  set(that.m_pData);
  return *this;
}

template<class T>
svt_vector3<T> &svt_vector3<T>::operator=(const Point3f &that)
{
  set(that);
  return *this;
}



template<class T>
inline svt_vector3<T> operator-(const svt_vector3<T> &p)
{
  svt_vector3<T> v;
  v[0] = -p[0];
  v[1] = -p[1];
  v[2] = -p[2];
  return v;
}


//
// arithmetic operators
// (need to be redefined because only first 3 components are taken into consideration)
//

// add 2 points
template<class T>
inline svt_vector3<T> operator+(const svt_vector3<T> &p1, const svt_vector3<T> &p2)
{
  return svt_vector3<T>(p1.x() + p2.x(), p1.y() + p2.y(), p1.z() + p2.z());
}


// add point and scalar
template<class T>
inline svt_vector3<T> operator+(const svt_vector3<T> &p, const T &f)
{
  return svt_vector3<T>(p.x() + f, p.y() + f, p.z() + f);
}


// add scalar and point
template<class T>
inline svt_vector3<T> operator+(const T &f, const svt_vector3<T> &p)
{
  return svt_vector3<T>(p.x() + f, p.y() + f, p.z() + f);
}


// substract 2 points
template<class T>
inline svt_vector3<T> operator-(const svt_vector3<T> &p1, const svt_vector3<T> &p2)
{
  return svt_vector3<T>(p1.x() - p2.x(), p1.y() - p2.y(), p1.z() - p2.z());
}


// substract point and scalar
template<class T>
inline svt_vector3<T> operator-(const svt_vector3<T> &p, const T &f)
{
  return svt_vector3<T>(p.x() - f, p.y() - f, p.z() - f);
}


// substract scalar and point
template<class T>
inline svt_vector3<T> operator-(const T &f, const svt_vector3<T> &p)
{
  return svt_vector3<T>(f - p.x(), f - p.y(), f - p.z());
}


// multiply/devide point with scalar
template<class T>
inline svt_vector3<T> operator*(const svt_vector3<T> &p, const T &f)
{
  return svt_vector3<T>(f * p.x(), f * p.y(), f * p.z());
}


template<class T>
inline svt_vector3<T> operator*(const T &f, const svt_vector3<T> &p)
{
  return svt_vector3<T>(f * p.x(), f * p.y(), f * p.z());
}


template<class T>
inline svt_vector3<T> operator/(const svt_vector3<T> &p, const T &f)
{
  return svt_vector3<T>(p.x() / f, p.y() / f, p.z() / f);
}


// Scalar Product
template<class T>
inline T operator*(const svt_vector3<T> &p1, const svt_vector3<T> &p2)
{
  return p1.x() * p2.x() + p1.y() * p2.y() + p1.z() * p2.z();
}


// Vector Product
template<class T>
inline svt_vector3<T> vectorProduct(const svt_vector3<T> &p1, const svt_vector3<T> &p2)
{
  return svt_vector3<T>(p1.y() * p2.z() - p1.z() * p2.y(),
                        p1.z() * p2.x() - p1.x() * p2.z(),
                        p1.x() * p2.y() - p1.y() * p2.x());
}
template<class T>
inline svt_vector3<T> crossProduct(const svt_vector3<T> &p1, const svt_vector3<T> &p2)
{
  return vectorProduct(p1, p2);
}


#endif


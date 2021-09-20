/*********************************************************************
*                           V O L F L T R                            *
**********************************************************************
* Program is part of the Situs package URL: situs.biomachina.org     *
* (c) Zbigniew Starosolski and Willy Wriggers, 2011                  *
**********************************************************************
*                                                                    *
* Filter tool for denoising maps and image stacks                    *
*                                                                    *
**********************************************************************
* See legal statement for terms of distribution                      *
*********************************************************************/


#include "lib_svt.h"
#ifdef _OPENMP
#include <omp.h>
#endif

#define FCOUT cout << "volfltr> "


// type of Kernel:
typedef enum {
  BETA_THETA,
  FISHER_VARIANCE_MEAN,
  ORG_VARIANCE,
  MAX_VARIANCE,
  ORG_FISHER,
  MEDIAN,
  FISHER_MEDIAN,
  FISHER_ORG_MAX_INC,
  VARIANCE,
  FISHER_VARIANCE,
  FISHER_MAX_INC_WITH_CENTER_MAX,
  FISHER_MAX_INC_WITH_CENTER_MEAN,
} kernel_type;

// type of neighborhood model:
typedef enum {
  TWO_DIMENSIONAL_4_NEIGHBORHOOD,
  TWO_DIMENSIONAL_8_NEIGHBORHOOD,
  THREE_DIMENSIONAL_6_NEIGHBORHOOD,
  THREE_DIMENSIONAL_26_NEIGHBORHOOD,
} path_type;

// note there are some differences between Situs and Sculptor implementations:
// - some of the functions below are member functions of classes in Sculptor's SVT
// - the interger distance bug fix (below) was handled differently
// - differences in program output / progress bar
// - all original DPSV filters available here in filtrGeodesicPath


void stepGeodesic(vector < vector<svt_vector3< int> > > &oPathVector, vector<svt_vector3< int> > &oPathSoFar, svt_vector3< int> oToPos, int iStepNumber, int iFullLength, int iCubeSize, int iPath)
// recursive function, finds all possible paths which fulfils requested condition: Mask size( iCubeSize), Paths length(iFullLength), and type of neighborhood model (iPath)
// i.e. next stepGeodesic at path to oToPos
{
  if ((oToPos.x() < 0 || oToPos.x() >= iCubeSize || oToPos.y() < 0 || oToPos.y() >= iCubeSize || oToPos.z() < 0 || oToPos.z() >= iCubeSize))
    return;

  svt_vector3< int> oOrigin(0, 0, 0);
  for (int i = 0; i < iStepNumber; i++)
    if ((oPathSoFar[i].distanceSq(oToPos) == 0)) { //  || (i > 0 &&  (oOrigin.distance(oToPos) <= oOrigin.distance(oPathSoFar[i]))))
      return;
    }

  oPathSoFar[iStepNumber] = oToPos; // remember current stepGeodesic
  svt_vector3< int>  oMargin((int)(iCubeSize * 0.5), (int)(iCubeSize * 0.5), (int)(iCubeSize * 0.5));
  if (iStepNumber == iFullLength) { // if this path has full length
    vector <svt_vector3< int> > newpath(iFullLength + 1); //temporary vector for the found path to be added to oPathVector
    for (int i = 0; i <= iFullLength; i++)
      newpath[i] = (oPathSoFar[i] - oMargin);

    oPathVector.push_back(newpath);
    return;
  }

  svt_vector3< int> oTempPosX(1, 0, 0);
  svt_vector3< int> oTempPosY(0, 1, 0);
  svt_vector3< int> oTempPosZ(0, 0, 1);
  svt_vector3< int> oTempPosXY(1, 1, 0);
  svt_vector3< int> oTempPosXZ(1, 0, 1);
  svt_vector3< int> oTempPoYZ(0, 1, 1);
  svt_vector3< int> oTempPosXmY(1, -1, 0);
  svt_vector3< int> oTempPosXmZ(1, 0, -1);
  svt_vector3< int> oTempPoYmZ(0, 1, -1);

  switch (iPath) {
    case TWO_DIMENSIONAL_4_NEIGHBORHOOD:
      stepGeodesic(oPathVector, oPathSoFar, oToPos + oTempPosY, iStepNumber + 1, iFullLength, iCubeSize, iPath); // next stepGeodesic: x, y+1, z
      stepGeodesic(oPathVector, oPathSoFar, oToPos + oTempPosX, iStepNumber + 1, iFullLength, iCubeSize, iPath); // next stepGeodesic: x+1, y, z
      stepGeodesic(oPathVector, oPathSoFar, oToPos - oTempPosY, iStepNumber + 1, iFullLength, iCubeSize, iPath); // next stepGeodesic: x, y-1, z
      stepGeodesic(oPathVector, oPathSoFar, oToPos - oTempPosX, iStepNumber + 1, iFullLength, iCubeSize, iPath); // next stepGeodesic: x-1, y, z
      break;
    case TWO_DIMENSIONAL_8_NEIGHBORHOOD:
      stepGeodesic(oPathVector, oPathSoFar, oToPos + oTempPosY, iStepNumber + 1, iFullLength, iCubeSize, iPath); // next stepGeodesic: x, y+1, z
      stepGeodesic(oPathVector, oPathSoFar, oToPos + oTempPosX, iStepNumber + 1, iFullLength, iCubeSize, iPath); // next stepGeodesic: x+1, y, z
      stepGeodesic(oPathVector, oPathSoFar, oToPos - oTempPosY, iStepNumber + 1, iFullLength, iCubeSize, iPath); // next stepGeodesic: x, y-1, z
      stepGeodesic(oPathVector, oPathSoFar, oToPos - oTempPosX, iStepNumber + 1, iFullLength, iCubeSize, iPath); // next stepGeodesic: x-1, y, z
      stepGeodesic(oPathVector, oPathSoFar, oToPos + oTempPosXY, iStepNumber + 1, iFullLength, iCubeSize, iPath); // next stepGeodesic: x+1, y+1, z
      stepGeodesic(oPathVector, oPathSoFar, oToPos - oTempPosXY, iStepNumber + 1, iFullLength, iCubeSize, iPath); // next stepGeodesic: x-1, y-1, z
      stepGeodesic(oPathVector, oPathSoFar, oToPos + oTempPosXmY, iStepNumber + 1, iFullLength, iCubeSize, iPath); // next stepGeodesic: x+1, y-1, z
      stepGeodesic(oPathVector, oPathSoFar, oToPos - oTempPosXmY, iStepNumber + 1, iFullLength, iCubeSize, iPath); // next stepGeodesic: x-1, y+1, z
      break;
    case THREE_DIMENSIONAL_6_NEIGHBORHOOD:
      stepGeodesic(oPathVector, oPathSoFar, oToPos + oTempPosZ, iStepNumber + 1, iFullLength, iCubeSize, iPath); // next stepGeodesic: x, y, z+1
      stepGeodesic(oPathVector, oPathSoFar, oToPos + oTempPosY, iStepNumber + 1, iFullLength, iCubeSize, iPath); // next stepGeodesic: x, y+1, z
      stepGeodesic(oPathVector, oPathSoFar, oToPos + oTempPosX, iStepNumber + 1, iFullLength, iCubeSize, iPath); // next stepGeodesic: x+1, y, z
      stepGeodesic(oPathVector, oPathSoFar, oToPos - oTempPosZ, iStepNumber + 1, iFullLength, iCubeSize, iPath); // next stepGeodesic: x, y, z-1
      stepGeodesic(oPathVector, oPathSoFar, oToPos - oTempPosY, iStepNumber + 1, iFullLength, iCubeSize, iPath); // next stepGeodesic: x, y-1, z
      stepGeodesic(oPathVector, oPathSoFar, oToPos - oTempPosX, iStepNumber + 1, iFullLength, iCubeSize, iPath); // next stepGeodesic: x-1, y, z
      break;
    case THREE_DIMENSIONAL_26_NEIGHBORHOOD:
      svt_vector3 <int> oTempPosLoop;
      for (int tx = (oToPos - oTempPosX).x(); tx <= (oToPos + oTempPosX).x(); tx++)
        for (int ty = (oToPos - oTempPosY).y(); ty <= (oToPos + oTempPosY).y(); ty++)
          for (int tz = (oToPos - oTempPosZ).z(); tz <= (oToPos + oTempPosZ).z(); tz++) {
            oTempPosLoop.x(tx);
            oTempPosLoop.y(ty);
            oTempPosLoop.z(tz);
            stepGeodesic(oPathVector, oPathSoFar, oTempPosLoop, iStepNumber + 1, iFullLength, iCubeSize, iPath); // next stepGeodesic: x, y, z-1
          }
      break;
  }
}

void findGeodesicPaths(vector < vector<svt_vector3< int> > > &oPathVector, svt_vector3< int> oToPos, int iFullLength, int iCubeSize, int iPath)
// function fills vector oPathSoFar with path into a vector < vector<svt_vector3< int> > > , generated by recursive function stepGeodesic
{
  vector< svt_vector3< int> > oPathSoFar(iFullLength + 1);
  oPathVector.empty();
  stepGeodesic(oPathVector, oPathSoFar, oToPos, 0, iFullLength , iCubeSize, iPath);
}

void FisherDiscriminant(vector <Real64> &fVectorVal_VarianceOnePath, vector <unsigned int > &iVectorVal_VarianceOnePath_Order, unsigned int &iK_star)
{
  unsigned int iN = fVectorVal_VarianceOnePath.size();
  vector<Real64> fVectorVal_VarianceOnePath_Sorted;

  vector< pair<Real64, int> > fVectorVal_VarianceOnePath_Sorted_Pair;

  for (unsigned int i = 0; i < fVectorVal_VarianceOnePath.size(); i++)
    fVectorVal_VarianceOnePath_Sorted_Pair.push_back(pair<Real64, int>(fVectorVal_VarianceOnePath[i] , i));

  sort(fVectorVal_VarianceOnePath_Sorted_Pair.begin(), fVectorVal_VarianceOnePath_Sorted_Pair.end());

  for (unsigned int i = 0; i < fVectorVal_VarianceOnePath.size(); i++) {
    iVectorVal_VarianceOnePath_Order[i] = fVectorVal_VarianceOnePath_Sorted_Pair[i].second;
    fVectorVal_VarianceOnePath_Sorted.push_back(fVectorVal_VarianceOnePath_Sorted_Pair[i].first);
  }

  Real64 fM1, fM2, fV1, fV2;

  vector<Real64> fVectorF;
  fVectorF.resize(fVectorVal_VarianceOnePath_Sorted.size());

  for (unsigned int iK = 0; iK < iN; iK++) {
    fM1 = 0.0;
    fM2 = 0.0;
    fV1 = 0.0;
    fV2 = 0.0;
    for (unsigned int iKtmpM1 = 0; iKtmpM1 < iK; iKtmpM1++)
      fM1 = fM1 + fVectorVal_VarianceOnePath_Sorted[iKtmpM1];
    fM1 = (1 / double(iK + 1)) * fM1;
    for (unsigned int iKtmpM2 = iK + 1; iKtmpM2 < iN; iKtmpM2++)
      fM2 = fM2 + fVectorVal_VarianceOnePath_Sorted[iKtmpM2];
    fM2 = (1 / double(iN - iK + 1)) * fM2;
    for (unsigned int iKtmpV1 = 0; iKtmpV1 < iK; iKtmpV1++)
      fV1 = fV1 + (fVectorVal_VarianceOnePath_Sorted[iKtmpV1] - fM1) * (fVectorVal_VarianceOnePath_Sorted[iKtmpV1] - fM1);
    for (unsigned int iKtmpV2 = iK + 1; iKtmpV2 < iN; iKtmpV2++)
      fV2 = fV2 + (fVectorVal_VarianceOnePath_Sorted[iKtmpV2] - fM2) * (fVectorVal_VarianceOnePath_Sorted[iKtmpV2] - fM2);
    // Fisher discriminant ratio
    fVectorF[iK] = ((fM1 - fM2) * (fM1 - fM2)) / (fV1 + fV2);
  }
  //max element of vector fVectorF
  iK_star = max_element(fVectorF.begin(), fVectorF.end()) - fVectorF.begin();

}

void medianOfVector(vector <Real64> &fVectorfVal_MeanOnePath, Real64 &fMedianOfVector)
{
  vector<Real64> fVectorVal_MeanOnePath_Sorted = fVectorfVal_MeanOnePath;

  stable_sort(fVectorVal_MeanOnePath_Sorted.begin(), fVectorVal_MeanOnePath_Sorted.end());

  unsigned int iN = fVectorVal_MeanOnePath_Sorted.size();
  unsigned int iIndexOfMedian;
  if ((iN % 2) != 0) {  //is odd
    iIndexOfMedian  = int(floor(double(iN) * 0.5));
    fMedianOfVector = fVectorVal_MeanOnePath_Sorted[iIndexOfMedian];
  } else {
    iIndexOfMedian  = iN / 2;
    fMedianOfVector =   fVectorVal_MeanOnePath_Sorted[iIndexOfMedian];
  }
  fVectorVal_MeanOnePath_Sorted.clear();
}


double intDistance(svt_vector3< int> iVec1, svt_vector3< int> iVec2)
// computes Euklidean distance for voxels on paths, this is safe for integers, unlike lib_sba.h: T distance(const svt_vector3<T>& oVec)
// note in Sculptor this is done differently, see Zbigniew e-mail 12/05/11
{
  double xv1, xv2, yv1, yv2, zv1, zv2;

  xv1 = (double)iVec1.x();
  yv1 = (double)iVec1.y();
  zv1 = (double)iVec1.z();
  xv2 = (double)iVec2.x();
  yv2 = (double)iVec2.y();
  zv2 = (double)iVec2.z();

  return sqrt((xv1 - xv2) * (xv1 - xv2) + (yv1 - yv2) * (yv1 - yv2) + (zv1 - zv2) * (zv1 - zv2));
}


void filtrGeodesicPath(svt_volume<Real64> &oVol, Real64 fBeta, Real64 fTheta, unsigned int iMaskDim,  unsigned int iPathLength, bool bProgress, int iKernel, int iPath)
//function performs filtration of volume. Parmeters:  fBeta - beta, fTheta - we do not use in default Kernel, iMaskDim - Mask Size ,  iPathLength - path length,  bool bProgress - progress info, int iKernel - type of kernel - set to 10 for DPSV filter, iPath - type of neighborhood model (values explained in the code)
{
  int iX, iY, iZ, iMargin;
  unsigned int iXTempP1, iYTempP1, iZTempP1, iXTempPi_1, iYTempPi_1, iZTempPi_1;
  svt_volume<Real64> oWalk3d;

  // create second map locally for output
  oWalk3d.allocate((unsigned int)oVol.getSizeX(), (unsigned int)oVol.getSizeY(), (unsigned int)oVol.getSizeZ(), 0.0);
  oWalk3d.setGridX(oVol.getGridX());
  oWalk3d.setGridY(oVol.getGridY());
  oWalk3d.setGridZ(oVol.getGridZ());
  oWalk3d.setWidth(oVol.getWidth());

  if (bProgress) FCOUT << "Applying DPSV filter..." <<  endl;

  vector < vector<svt_vector3< int> > > oPathVector;  // paths filled in recursive stepGeodesic()
  vector <double> oSum2RoValPo;          // for given path stores precomputed sum of distances from element [1] to [2], [3] and all following ones
  vector <int> oPathFirstIndexVector;        // indexes in oPathVector of first path for a neighbor
  vector <int> oPathLastIndexVector;         // indexes in oPathVector of last path for a neighbor
  double fDist_P0_Pi_0;
  iMargin = (int)(iMaskDim * 0.5);

  svt_vector3<int> oFirstStep1(iMargin, iMargin, iMargin);

  // fill oPathVector - Vector containing set of paths
  findGeodesicPaths(oPathVector, oFirstStep1, iPathLength, iMaskDim, iPath);

  int iNoPath = oPathVector.size();
  FCOUT << "Number of paths in mask:  " <<  iNoPath << endl;

  //fill oSum2RoValPo - Vector containing set of distances
  oSum2RoValPo.resize(iNoPath);
  for (int iSetPath = 0; iSetPath < iNoPath; iSetPath++) {
    double fSumDistance = 0.0;
    for (unsigned int iV = 2; iV <= iPathLength; iV++) fSumDistance += intDistance(oPathVector[iSetPath][0], oPathVector[iSetPath][iV]);
    oSum2RoValPo[iSetPath] = fSumDistance;
  }

  //fill oPathFirstIndexVector and oPathLastIndexVector
  for (int iSetPath = 0; iSetPath < iNoPath; iSetPath++) {
    if (iSetPath == 0 || !(oPathVector[iSetPath - 1][1] == oPathVector[iSetPath][1]))
      oPathFirstIndexVector.push_back(iSetPath);
    if (iSetPath == iNoPath - 1 || !(oPathVector[iSetPath][1] == oPathVector[iSetPath + 1][1])) // last path through this neighbor
      oPathLastIndexVector.push_back(iSetPath);
  }

  int iNoNeighborPixel = oPathFirstIndexVector.size();
  FCOUT << "Model of neighborhood : " << iNoNeighborPixel << endl;
  vector <Real64> fVectorfSetLambdaTemp;
  fVectorfSetLambdaTemp.resize(oPathVector.size());
  vector <Real64> fVectorfVal_MeanOnePath;
  vector <unsigned int> iVectorfVal_MeanOnePath_Order;
  vector <Real64> fVectorVal_MaxABSDiffInPath;
  vector <Real64> fVectorVal_MaxABSDiffInPath_Sorted;
  vector <unsigned int> iVectorVal_MaxABSDiffInPath_Order;
  vector <Real64> fVectorVal_VarianceOnePath_Sorted;
  vector <Real64> fVectorVal_VarianceOnePath;
  vector <unsigned int> iVectorVal_VarianceOnePath_Order;

  unsigned int iNoOfChangedtoAverage = 0;
  unsigned int iNoOfAllPath = 0;

  int iXmax, iYmax, iZmax;
  iXmax = (int) oVol.getSizeX();
  iYmax = (int) oVol.getSizeY();
  iZmax = (int) oVol.getSizeZ();

  // switching between different kernels (set to 10: DPSV_filter)
  switch (iKernel) {
  /*----*/case BETA_THETA:
      for (iZ = 0; iZ < iZmax; iZ++) { // Z size of volume
        for (iX = 0; iX < iXmax; iX++) { // X size of volume
          for (iY = 0; iY < iYmax; iY++) { // Y size of volume
            double  fSetLambda = 0.0, fHatP0Temp = 0.0;
            for (int iNeighbor = 0; iNeighbor < iNoNeighborPixel; iNeighbor++) { // for all neighbors
              double fSetLambdaTemp = 0.0;
              iXTempP1 = oPathVector[oPathFirstIndexVector[iNeighbor]][1].x() + iX;
              iYTempP1 = oPathVector[oPathFirstIndexVector[iNeighbor]][1].y() + iY;
              iZTempP1 = oPathVector[oPathFirstIndexVector[iNeighbor]][1].z() + iZ;
              double fValP1   = oVol.getValue(iXTempP1, iYTempP1, iZTempP1);

              for (int iSetPath = oPathFirstIndexVector[iNeighbor]; iSetPath <= oPathLastIndexVector[iNeighbor]; iSetPath++) { //for all paths through iNeighbor
                double fLambda = oSum2RoValPo[iSetPath] * fTheta;
                for (unsigned int iV = 2; iV <= iPathLength; iV++) { // for path elements 2
                  iXTempPi_1 = oPathVector[iSetPath][iV].x() + iX;
                  iYTempPi_1 = oPathVector[iSetPath][iV].y() + iY;
                  iZTempPi_1 = oPathVector[iSetPath][iV].z() + iZ;
                  double fValPi_1 = oVol.getValue(iXTempPi_1, iYTempPi_1, iZTempPi_1);
                  fLambda  += fabs(fValP1 - fValPi_1);
                }
                fSetLambdaTemp += exp((-fBeta) * fLambda);
              }
              fHatP0Temp += fSetLambdaTemp * fValP1;
              fSetLambda += fSetLambdaTemp;
            }
            fHatP0Temp = fHatP0Temp / fSetLambda;
            double fHatP0 = (fHatP0Temp);
            oWalk3d.setValue(iX, iY, iZ, fHatP0);
          }
        }
        if (bProgress) cout << ".";
      }
      break;

  /*----*/case ORG_FISHER:
      for (iZ = 0; iZ < iZmax; iZ++) { // Z size of volume
        for (iX = 0; iX < iXmax; iX++) { // X size of volume
          for (iY = 0; iY < iYmax; iY++) { // Y size of volume
            double fSetLambda = 0.0, fHatP0Temp = 0.0;
            //double fValP0   = oVol.getValue(iX,iY,iZ);
            for (int iNeighbor = 0; iNeighbor < iNoNeighborPixel; iNeighbor++) { // for all neighbors
              double fSetLambdaTemp = 0.0;
              double fValP1   = oVol. getValue(oPathVector[oPathFirstIndexVector[iNeighbor]][1].x() + iX,
                                               oPathVector[oPathFirstIndexVector[iNeighbor]][1].y() + iY,
                                               oPathVector[oPathFirstIndexVector[iNeighbor]][1].z() + iZ);
              fVectorfVal_MeanOnePath.clear();
              iVectorfVal_MeanOnePath_Order.clear();
              for (int iSetPath = oPathFirstIndexVector[iNeighbor]; iSetPath <= oPathLastIndexVector[iNeighbor]; iSetPath++) { //for all paths through iNeighbor
                double fVal_MeanOnePathTmp = 0.0; //oSum2RoValPo[iSetPath]*fTheta; // function of spatial position * fTheta;
                for (unsigned int iV = 2; iV <= iPathLength; iV++) {  // caclulation Lambda and mean for path (elements 2 up to iPathLength)
                  double fValPi_0 = oVol.getValue(oPathVector[iSetPath][iV].x() + iX,    oPathVector[iSetPath][iV].y() + iY,  oPathVector[iSetPath][iV].z() + iZ);
                  fVal_MeanOnePathTmp += fabs(fValP1 - fValPi_0); // temporary value for calculation mean of One Path
                }
                fVectorfVal_MeanOnePath.push_back(fVal_MeanOnePathTmp);
              }
              unsigned int iK_star = 0;                    //###################### FISHER DISCRIMINANT FOR PATH VARIANCE CHOISE
              iVectorfVal_MeanOnePath_Order.resize(fVectorfVal_MeanOnePath.size(), 0);

              FisherDiscriminant(fVectorfVal_MeanOnePath, iVectorfVal_MeanOnePath_Order, iK_star);
              if ((iK_star == 0) || (iK_star ==  fVectorfVal_MeanOnePath.size() - 1)) // condition for filtering ( if is FIRST or LAST thats means there is no Partition)
                fSetLambdaTemp += exp((-fBeta) * fVectorfVal_MeanOnePath[iVectorfVal_MeanOnePath_Order[0]]);
              if ((iK_star > 0) && (iK_star <  fVectorfVal_MeanOnePath.size() - 1)) { // condition for filtering ( if it is beteen FIRST and LAST ther is patrition, and we took only that is below iK_star )
                fSetLambdaTemp += exp((-fBeta) * fVectorfVal_MeanOnePath[iVectorfVal_MeanOnePath_Order[iK_star]]);
                iNoOfChangedtoAverage += 1;
              }
              iNoOfAllPath += 1;
              fHatP0Temp += fSetLambdaTemp * fValP1;
              fSetLambda += fSetLambdaTemp;
            }
            double fVoxelValTmp = fHatP0Temp / fSetLambda;
            if ((fSetLambda < 1.0E-200) && (fSetLambda > -1.0E-200)) {
              fVoxelValTmp = 0.0;
            }
            oWalk3d.setValue(iX, iY, iZ, fVoxelValTmp);
          }
        }
        if (bProgress) cout << ".";
      }
      FCOUT << "ORG_FISHER case, parameter beta: " << fBeta << " fTheta: " << fTheta << endl;
      FCOUT << "Fraction of paths partitioned by Fisher disciminant analysis: "  << iNoOfChangedtoAverage << " / " << iNoOfAllPath << " = " << (((double)iNoOfChangedtoAverage) / ((double)iNoOfAllPath)) << endl;
      break;

  /*----*/case FISHER_MAX_INC_WITH_CENTER_MEAN:
      for (iZ = 0; iZ < iZmax; iZ++) { // Z size of volume
        for (iX = 0; iX < iXmax; iX++) { // X size of volume
          for (iY = 0; iY < iYmax; iY++) { // Y size of volume
            double fSetLambda = 0.0, fHatP0Temp = 0.0;
            double fValP0   = oVol.getValue(iX, iY, iZ);
            for (int iNeighbor = 0; iNeighbor < iNoNeighborPixel; iNeighbor++) { // for all neighbors
              double fSetLambdaTemp = 0.0;
              double fValP1   = oVol. getValue(oPathVector[oPathFirstIndexVector[iNeighbor]][1].x() + iX,
                                               oPathVector[oPathFirstIndexVector[iNeighbor]][1].y() + iY,
                                               oPathVector[oPathFirstIndexVector[iNeighbor]][1].z() + iZ);
              fVectorVal_MaxABSDiffInPath.clear();
              fVectorVal_MaxABSDiffInPath_Sorted.clear();
              fVectorfVal_MeanOnePath.clear();
              iVectorVal_MaxABSDiffInPath_Order.clear();

              for (int iSetPath = oPathFirstIndexVector[iNeighbor]; iSetPath <= oPathLastIndexVector[iNeighbor]; iSetPath++) { //for all paths through iNeighbor
                double fVal_MeanOnePathTmp = 0.0; //fabs(fValP0-fValP1) + oSum2RoValPo[iSetPath]*fTheta; // function of spatial position * fTheta;
                double fVal_MaxABSDiffInPath = 0.0;
                fDist_P0_Pi_0 = intDistance(oPathVector[iSetPath][0], oPathVector[iSetPath][1]);
                double fValCurrMAX = fabs(fValP0 - fValP1) / fDist_P0_Pi_0;
                for (unsigned int iV = 2; iV <= iPathLength; iV++) {  // caclulation Lambda and mean for path (elements 2 up to iPathLength)
                  double fValPi_0 = oVol.getValue(oPathVector[iSetPath][iV].x() + iX,    oPathVector[iSetPath][iV].y() + iY,  oPathVector[iSetPath][iV].z() + iZ);
                  double fValPi_1   = oVol.getValue(oPathVector[iSetPath][iV - 1].x() + iX, oPathVector[iSetPath][iV - 1].y() + iY, oPathVector[iSetPath][iV - 1].z() + iZ);
                  fVal_MeanOnePathTmp += fabs(fValPi_1 - fValPi_0); // temporary value for calculation mean of One Path
                  fDist_P0_Pi_0 = intDistance(oPathVector[iSetPath][0], oPathVector[iSetPath][iV]);
                  if (fValCurrMAX < fabs(fValPi_1 - fValPi_0) / fDist_P0_Pi_0)     // choosing Current Maximal increment for path.
                    fValCurrMAX = fabs(fValPi_1 - fValPi_0) / fDist_P0_Pi_0;
                  if (fValCurrMAX < fabs(fValP0 - fValPi_0) / fDist_P0_Pi_0) // comparing with increment to center voxel and take maximal.
                    fValCurrMAX = fabs(fValP0 - fValPi_0) / fDist_P0_Pi_0;
                  fVal_MaxABSDiffInPath = fValCurrMAX;
                }
                fVectorfVal_MeanOnePath.push_back(fVal_MeanOnePathTmp);
                fVectorVal_MaxABSDiffInPath.push_back(fVal_MaxABSDiffInPath);   // vector with variances for all path through iNeighbor
              }
              unsigned int iK_star = 0;                    //###################### FISHER DISCRIMINANT FOR PATH VARIANCE CHOISE
              iVectorVal_MaxABSDiffInPath_Order.resize(fVectorVal_MaxABSDiffInPath.size(), 0);
              FisherDiscriminant(fVectorVal_MaxABSDiffInPath, iVectorVal_MaxABSDiffInPath_Order, iK_star);
              if ((iK_star == 0) || (iK_star == fVectorVal_MaxABSDiffInPath.size() - 1)) // condition for filtering ( if is FIRST or LAST thats means there is no Partition)
                fSetLambdaTemp += exp((-fBeta) * fVectorfVal_MeanOnePath[iVectorVal_MaxABSDiffInPath_Order[0]]);
              if ((iK_star > 0) && (iK_star <  fVectorVal_MaxABSDiffInPath.size() - 1)) { // condition for filtering ( if it is beteen FIRST and LAST ther is patrition, and we took only that is below iK_star )
                fSetLambdaTemp += exp((-fBeta) * fVectorfVal_MeanOnePath[iVectorVal_MaxABSDiffInPath_Order[iK_star]]);
                iNoOfChangedtoAverage += 1;
              }
              iNoOfAllPath += 1;
              fHatP0Temp += fSetLambdaTemp * fValP1;
              fSetLambda += fSetLambdaTemp;
            }
            oWalk3d.setValue(iX, iY, iZ, ((fHatP0Temp / fSetLambda)));
          }
        }
        if (bProgress) cout << ".";
      }
      FCOUT << "FISHER_MAX_INC_WITH_CENTER_MEAN case, parameter beta: " << fBeta << " fTheta: " << fTheta << endl;
      FCOUT << "Fraction of paths partitioned by Fisher disciminant analysis: "  << iNoOfChangedtoAverage << " / " << iNoOfAllPath << " = " << (((double)iNoOfChangedtoAverage) / ((double)iNoOfAllPath)) << endl;
      break;

  /*----*/case FISHER_MAX_INC_WITH_CENTER_MAX:


      FCOUT << "Processing " << iZmax << " Z-sections..." << endl;

      double fSetLambda, fHatP0Temp, fValP0, fSetLambdaTemp, fValP1, fVal_MaxABSDiffInPath, fValCurrMAX, fValPi_0, fValPi_1, fVoxelValTmp;
      int iNeighbor, iSetPath;
      unsigned int iV, iK_star;


#ifdef _OPENMP
      #pragma omp parallel for \
      shared(oWalk3d,iNoOfChangedtoAverage,iNoOfAllPath,oPathFirstIndexVector,oPathLastIndexVector,oPathVector,iPathLength) \
      private(iX,iY,fSetLambda,fHatP0Temp,fValP0,iNeighbor,fSetLambdaTemp,fValP1,iSetPath,fVal_MaxABSDiffInPath,fValCurrMAX,iV,fValPi_0,fValPi_1,fDist_P0_Pi_0,iK_star,fVoxelValTmp,fVectorVal_MaxABSDiffInPath,iVectorVal_MaxABSDiffInPath_Order) \
      schedule(dynamic, 1)
#endif
      for (iZ = 0; iZ < iZmax; iZ++) { // Z size of volume
        for (iX = 0; iX < iXmax; iX++) { // X size of volume
          for (iY = 0; iY < iYmax; iY++) { // Y size of volume
            fSetLambda = 0.0;
            fHatP0Temp = 0.0;
            fValP0   = oVol.getValue(iX, iY, iZ);
            for (iNeighbor = 0; iNeighbor < iNoNeighborPixel; iNeighbor++) { // for all neighbors
              fSetLambdaTemp = 0.0;
              fValP1   = oVol. getValue(oPathVector[oPathFirstIndexVector[iNeighbor]][1].x() + iX,
                                        oPathVector[oPathFirstIndexVector[iNeighbor]][1].y() + iY,
                                        oPathVector[oPathFirstIndexVector[iNeighbor]][1].z() + iZ);
              // Loop: among path throught iNeighbor - typically there are 6 or 26 neighbors.

              fVectorVal_MaxABSDiffInPath.clear();

              for (iSetPath = oPathFirstIndexVector[iNeighbor]; iSetPath <= oPathLastIndexVector[iNeighbor]; iSetPath++) { //for all paths through iNeighbor
                fVal_MaxABSDiffInPath = 0.0;
                fDist_P0_Pi_0 = intDistance(oPathVector[iSetPath][0], oPathVector[iSetPath][1]);
                fValCurrMAX = fabs(fValP0 - fValP1) / fDist_P0_Pi_0;
                for (iV = 2; iV <= iPathLength; iV++) {   // calculation Lambda and mean for path (elements 2 up to iPathLength)
                  fValPi_0 = oVol.getValue(oPathVector[iSetPath][iV].x() + iX,    oPathVector[iSetPath][iV].y() + iY,  oPathVector[iSetPath][iV].z() + iZ);
                  fValPi_1   = oVol.getValue(oPathVector[iSetPath][iV - 1].x() + iX, oPathVector[iSetPath][iV - 1].y() + iY, oPathVector[iSetPath][iV - 1].z() + iZ);
                  fDist_P0_Pi_0 = intDistance(oPathVector[iSetPath][0], oPathVector[iSetPath][iV]);
                  if (fValCurrMAX < fabs(fValPi_1 - fValPi_0) / (fDist_P0_Pi_0))  // choosing Current Maximal increment for path.
                    fValCurrMAX = fabs(fValPi_1 - fValPi_0) / (fDist_P0_Pi_0);
                  if (fValCurrMAX < fabs(fValP0 - fValPi_0) / (fDist_P0_Pi_0)) // comparing with increment to center voxel and take maximal.
                    fValCurrMAX = fabs(fValP0 - fValPi_0) / (fDist_P0_Pi_0);
                  fVal_MaxABSDiffInPath = fValCurrMAX;
                }
                fVectorVal_MaxABSDiffInPath.push_back(fVal_MaxABSDiffInPath);   // vector with variances for all path through iNeighbor
              }
              iK_star = 0;                   //###################### FISHER DISCRIMINANT FOR PATH VARIANCE CHOICE
              iVectorVal_MaxABSDiffInPath_Order.resize(fVectorVal_MaxABSDiffInPath.size(), 0);
              FisherDiscriminant(fVectorVal_MaxABSDiffInPath, iVectorVal_MaxABSDiffInPath_Order, iK_star);
              if ((iK_star == 0) || (iK_star == fVectorVal_MaxABSDiffInPath.size() - 1)) // condition for filtering (if is FIRST or LAST thats means there is no partition -> we take first one)
                fSetLambdaTemp += exp((-fBeta) * fVectorVal_MaxABSDiffInPath[iVectorVal_MaxABSDiffInPath_Order[0]]);
              if ((iK_star > 0) && (iK_star <  fVectorVal_MaxABSDiffInPath.size() - 1)) { // condition for filtering (if it is beteen FIRST and LAST there is FDA partition -> we take iK_star )
                fSetLambdaTemp += exp((-fBeta) * fVectorVal_MaxABSDiffInPath[iVectorVal_MaxABSDiffInPath_Order[iK_star]]);
#ifdef _OPENMP
                #pragma omp critical
#endif
                iNoOfChangedtoAverage += 1;
              }
#ifdef _OPENMP
              #pragma omp critical
#endif
              iNoOfAllPath += 1;
              fHatP0Temp += fSetLambdaTemp * fValP1;
              fSetLambda += fSetLambdaTemp;
            }
            if (abs(fSetLambda) < 1.0E-200) fVoxelValTmp = 0.0;   //  avoid division by zero
            else fVoxelValTmp = fHatP0Temp / fSetLambda;
#ifdef _OPENMP
            #pragma omp critical
#endif
            oWalk3d.setValue(iX, iY, iZ, fVoxelValTmp);
          }
        }
#ifdef _OPENMP
        if (bProgress) {
          FCOUT << "Z-section " << iZ + 1 << " processed by thread nr. " <<  omp_get_thread_num() + 1 << " of " << omp_get_num_threads() << "." << endl;
        }
#else
        if (bProgress) FCOUT << "Z-section " << iZ + 1 << " processed serially." << endl;
#endif
      }
      FCOUT << "DPSV filter applied with mask size " << iMaskDim << ", path length " << iPathLength <<  ", beta "  << fBeta << ", " << iNoNeighborPixel << "-neighborhood model. " <<  endl;
      FCOUT << "Fraction of paths partitioned by Fisher disciminant analysis: "  << iNoOfChangedtoAverage << " / " << iNoOfAllPath << " = " << (((double)iNoOfChangedtoAverage) / ((double)iNoOfAllPath)) << endl;
      break;

  /*----*/case FISHER_VARIANCE_MEAN:
      for (iZ = 0; iZ < iZmax; iZ++) { // Z size of volume
        for (iX = 0; iX < iXmax; iX++) { // X size of volume
          for (iY = 0; iY < iYmax; iY++) { // Y size of volume
            //double fHatP0 =0.0;
            double fSetLambda = 0.0, fHatP0Temp = 0.0;
            double fVal_MeanOnePath = 0.0, fVal_MeanOnePathTmp = 0.0;
            // double fLambda=0.0;
            double fValPi_1 = 0.0;
            //  double fValP1 = 0.0;
            double fSetLambdaTemp = 0.0;
            double fVal_VarianceOnePathTmp = 0.0, fVal_VarianceOnePath = 0.0;

            for (int iNeighbor = 0; iNeighbor < iNoNeighborPixel; iNeighbor++) { // for all neighbors
              fSetLambdaTemp = 0.0;
              fVal_VarianceOnePathTmp = 0.0;
              double fValP1   = oVol. getValue(oPathVector[oPathFirstIndexVector[iNeighbor]][1].x() + iX,
                                               oPathVector[oPathFirstIndexVector[iNeighbor]][1].y() + iY,
                                               oPathVector[oPathFirstIndexVector[iNeighbor]][1].z() + iZ);
              fVectorVal_VarianceOnePath.clear();
              fVectorfVal_MeanOnePath.clear();
              for (int iSetPath = oPathFirstIndexVector[iNeighbor]; iSetPath <= oPathLastIndexVector[iNeighbor]; iSetPath++) { //for all paths through iNeighbor
                fVal_MeanOnePathTmp = fValP1;
                for (unsigned int iV = 2; iV <= iPathLength; iV++) {  // caclulation Lambda and mean for path (elements 2 up to iPathLength)
                  fValPi_1 = oVol.getValue(oPathVector[iSetPath][iV].x() + iX, oPathVector[iSetPath][iV].y() + iY, oPathVector[iSetPath][iV].z() + iZ);
                  fVal_MeanOnePathTmp = fVal_MeanOnePathTmp + fValPi_1;   // temporary value for calculation mean of One Path
                }
                fVal_MeanOnePath =  fVal_MeanOnePathTmp / iPathLength;    // mean for One Path in mask
                fVectorfVal_MeanOnePath.push_back(fVal_MeanOnePath);  // comment for the monent
                for (unsigned int iVV = 2; iVV <= iPathLength; iVV++) {   // caclulation variance  for path (elements 2 up to iPathLength)
                  fValPi_1 = oVol.getValue(oPathVector[iSetPath][iVV].x() + iX, oPathVector[iSetPath][iVV].y() + iY, oPathVector[iSetPath][iVV].z() + iZ);
                  fVal_VarianceOnePathTmp =  fVal_VarianceOnePathTmp + (fValPi_1 - fVal_MeanOnePath) * (fValPi_1 - fVal_MeanOnePath); // temporary value for calculation variance of One Path
                }
                fVal_VarianceOnePath = fVal_VarianceOnePathTmp / iPathLength;   // variance for one path
                fVectorVal_VarianceOnePath.push_back(fVal_VarianceOnePath);     // vector with variances for all path through iNeighbor
              }
              //###################### FISHER DISCRIMINANT FOR PATH VARIANCE CHOISE
              unsigned int iK_star = 0;
              iVectorVal_VarianceOnePath_Order.resize(fVectorVal_VarianceOnePath.size(), 0);
              FisherDiscriminant(fVectorVal_VarianceOnePath, iVectorVal_VarianceOnePath_Order, iK_star);
              //###################### FISHER DISCRIMINANT FOR PATH VARIANCE CHOISE
              if ((iK_star == 0) || (iK_star == fVectorVal_VarianceOnePath.size() - 1)) // condition for filtering ( if is FIRST or LAST thats means there is no Partition)
                fSetLambdaTemp += exp((-fBeta) * fVectorfVal_MeanOnePath[iVectorVal_VarianceOnePath_Order[0]]);
              if ((iK_star > 0) && (iK_star <  fVectorVal_VarianceOnePath.size() - 1)) { // condition for filtering ( if it is beteen FIRST and LAST ther is patrition, and we took only that is below iK_star )
                fSetLambdaTemp += exp((-fBeta) * fVectorfVal_MeanOnePath[iVectorVal_VarianceOnePath_Order[iK_star]]);
                iNoOfChangedtoAverage += 1;
              }
              iNoOfAllPath += 1;
              fHatP0Temp += fSetLambdaTemp * fValP1;
              fSetLambda += fSetLambdaTemp;
            }
            oWalk3d.setValue(iX, iY, iZ, fHatP0Temp / fSetLambda);

          }
        }
        if (bProgress) cout << ".";
      }
      FCOUT << "FISHER_VARIANCE_MEAN case, parameter beta: " << fBeta << " fTheta: " << fTheta << endl;
      FCOUT << "Fraction of paths partitioned by Fisher disciminant analysis: "  << iNoOfChangedtoAverage << " / " << iNoOfAllPath << " = " << (((double)iNoOfChangedtoAverage) / ((double)iNoOfAllPath)) << endl;
      break;

  /*----*/case FISHER_VARIANCE:
      for (iZ = 0; iZ < iZmax; iZ++) { // Z size of volume
        for (iX = 0; iX < iXmax; iX++) { // X size of volume
          for (iY = 0; iY < iYmax; iY++) { // Y size of volume
            double fSetLambda = 0.0, fHatP0Temp = 0.0;
            double fVal_MeanOnePath = 0.0, fVal_MeanOnePathTmp = 0.0;
            double fValPi_1 = 0.0;
            double fSetLambdaTemp = 0.0;
            double fVal_VarianceOnePathTmp = 0.0, fVal_VarianceOnePath = 0.0;

            for (int iNeighbor = 0; iNeighbor < iNoNeighborPixel; iNeighbor++) { // for all neighbors
              fSetLambdaTemp = 0.0;
              fVal_VarianceOnePathTmp = 0.0;
              double fValP1   = oVol. getValue(oPathVector[oPathFirstIndexVector[iNeighbor]][1].x() + iX,
                                               oPathVector[oPathFirstIndexVector[iNeighbor]][1].y() + iY,
                                               oPathVector[oPathFirstIndexVector[iNeighbor]][1].z() + iZ);
              fVectorfVal_MeanOnePath.clear();
              fVectorVal_VarianceOnePath.clear();
              for (int iSetPath = oPathFirstIndexVector[iNeighbor]; iSetPath <= oPathLastIndexVector[iNeighbor]; iSetPath++) { //for all paths through iNeighbor
                fVal_MeanOnePathTmp = fValP1;
                for (unsigned int iV = 2; iV <= iPathLength; iV++) {  // caclulation Lambda and mean for path (elements 2 up to iPathLength)
                  fValPi_1 = oVol.getValue(oPathVector[iSetPath][iV].x() + iX, oPathVector[iSetPath][iV].y() + iY, oPathVector[iSetPath][iV].z() + iZ);
                  fVal_MeanOnePathTmp = fVal_MeanOnePathTmp + fValPi_1;   // temporary value for calculation mean of One Path
                }
                fVal_MeanOnePath =  fVal_MeanOnePathTmp / iPathLength;    // mean for One Path in mask
                fVectorfVal_MeanOnePath.push_back(fVal_MeanOnePath);  // comment for the monent
                for (unsigned int iVV = 2; iVV <= iPathLength; iVV++) {   // caclulation variance  for path (elements 2 up to iPathLength)
                  fValPi_1 = oVol.getValue(oPathVector[iSetPath][iVV].x() + iX, oPathVector[iSetPath][iVV].y() + iY, oPathVector[iSetPath][iVV].z() + iZ);
                  fVal_VarianceOnePathTmp =  fVal_VarianceOnePathTmp + (fValPi_1 - fVal_MeanOnePath) * (fValPi_1 - fVal_MeanOnePath); // temporary value for calculation variance of One Path
                }
                fVal_VarianceOnePath = fVal_VarianceOnePathTmp / iPathLength;   // variance for one path
                fVectorVal_VarianceOnePath.push_back(fVal_VarianceOnePath);     // vector with variances for all path through iNeighbor
              }
              //###################### FISHER DISCRIMINANT FOR PATH VARIANCE CHOISE
              unsigned int iK_star = 0;
              iVectorVal_VarianceOnePath_Order.resize(fVectorVal_VarianceOnePath.size(), 0);
              FisherDiscriminant(fVectorVal_VarianceOnePath, iVectorVal_VarianceOnePath_Order, iK_star);
              //###################### FISHER DISCRIMINANT FOR PATH VARIANCE CHOISE

              if ((iK_star == 0) || (iK_star == fVectorVal_VarianceOnePath.size() - 1)) // condition for filtering ( if is FIRST or LAST thats means there is no Partition)
                fSetLambdaTemp += exp((-fBeta) * fVectorVal_VarianceOnePath[iVectorVal_VarianceOnePath_Order[0]]);
              if ((iK_star > 0) && (iK_star <  fVectorVal_VarianceOnePath.size() - 1)) { // condition for filtering ( if it is beteen FIRST and LAST ther is patrition, and we took only that is below iK_star )
                fSetLambdaTemp += exp((-fBeta) * fVectorVal_VarianceOnePath[iVectorVal_VarianceOnePath_Order[iK_star]]);
                iNoOfChangedtoAverage += 1;
              }
              iNoOfAllPath += 1;
              fHatP0Temp += fSetLambdaTemp * fValP1;
              fSetLambda += fSetLambdaTemp;
            }
            oWalk3d.setValue(iX, iY, iZ, fHatP0Temp / fSetLambda);

          }
        }
        if (bProgress) cout << ".";
      }
      FCOUT << "FISHER_VARIANCE case, parameter beta: " << fBeta << " fTheta: " << fTheta << endl;
      FCOUT << "Fraction of paths partitioned by Fisher disciminant analysis: "  << iNoOfChangedtoAverage << " / " << iNoOfAllPath << " = " << (((double)iNoOfChangedtoAverage) / ((double)iNoOfAllPath)) << endl;
      break;

  /*----*/case VARIANCE:
      for (iZ = 0; iZ < iZmax; iZ++) { // Z size of volume
        for (iX = 0; iX < iXmax; iX++) { // X size of volume
          for (iY = 0; iY < iYmax; iY++) { // Y size of volume
            double fSetLambda = 0.0, fHatP0Temp = 0.0;
            double fVal_MeanOnePath = 0.0, fVal_MeanOnePathTmp = 0.0;
            double fValPi_1 = 0.0;
            double fSetLambdaTemp = 0.0;
            double fVal_VarianceOnePathTmp = 0.0, fVal_VarianceOnePath = 0.0;

            for (int iNeighbor = 0; iNeighbor < iNoNeighborPixel; iNeighbor++) { // for all neighbors
              fSetLambdaTemp = 0.0;
              fVal_VarianceOnePathTmp = 0.0;
              double fValP1   = oVol. getValue(oPathVector[oPathFirstIndexVector[iNeighbor]][1].x() + iX,
                                               oPathVector[oPathFirstIndexVector[iNeighbor]][1].y() + iY,
                                               oPathVector[oPathFirstIndexVector[iNeighbor]][1].z() + iZ);
              for (int iSetPath = oPathFirstIndexVector[iNeighbor]; iSetPath <= oPathLastIndexVector[iNeighbor]; iSetPath++) { //for all paths through iNeighbor
                fVal_MeanOnePathTmp = fValP1;
                for (unsigned int iV = 2; iV <= iPathLength; iV++) {  // caclulation Lambda and mean for path (elements 2 up to iPathLength)
                  fValPi_1 = oVol.getValue(oPathVector[iSetPath][iV].x() + iX, oPathVector[iSetPath][iV].y() + iY, oPathVector[iSetPath][iV].z() + iZ);
                  fVal_MeanOnePathTmp = fVal_MeanOnePathTmp + fValPi_1;   // temporary value for calculation mean of One Path
                }
                fVal_MeanOnePath =  fVal_MeanOnePathTmp / iPathLength;    // mean for One Path in mask
                for (unsigned int iVV = 2; iVV <= iPathLength; iVV++) {   // caclulation variance  for path (elements 2 up to iPathLength)
                  fValPi_1 = oVol.getValue(oPathVector[iSetPath][iVV].x() + iX, oPathVector[iSetPath][iVV].y() + iY, oPathVector[iSetPath][iVV].z() + iZ);
                  fVal_VarianceOnePathTmp =  fVal_VarianceOnePathTmp + (fValPi_1 - fVal_MeanOnePath) * (fValPi_1 - fVal_MeanOnePath); // temporary value for calculation variance of One Path
                }
                fVal_VarianceOnePath = fVal_VarianceOnePathTmp / iPathLength;   // variance for one path
                fSetLambdaTemp += exp((-fBeta) * fVal_VarianceOnePath);
              }
              fHatP0Temp += fSetLambdaTemp * fValP1;
              fSetLambda += fSetLambdaTemp;
            }
            oWalk3d.setValue(iX, iY, iZ, fHatP0Temp / fSetLambda);

          }
        }
        if (bProgress) cout << ".";
      }
      FCOUT << "VARIANCE case, parameter beta: " << fBeta << " fTheta: " << fTheta << endl;
      FCOUT << "Fraction of paths partitioned by Fisher disciminant analysis: "  << iNoOfChangedtoAverage << " / " << iNoOfAllPath << " = " << (((double)iNoOfChangedtoAverage) / ((double)iNoOfAllPath)) << endl;
      break;

  /*----*/case MEDIAN:
      vector <Real64> fVectorVal_MedianOnePath;
      vector <Real64> fVectorfVal_MeanAllPath;
      for (iZ = 0; iZ < iZmax; iZ++) { // Z size of volume
        for (iX = 0; iX < iXmax; iX++) { // X size of volume
          for (iY = 0; iY < iYmax; iY++) { // Y size of volume
            double fValPi_1 = 0.0;
            double fMedianOnePath = 0.0;
            double fMedianOfVector2 = 0.0;
            fVectorfVal_MeanAllPath.clear();
            fValP0   = oVol.getValue(iX, iY, iZ); // value of first neighboor.
            for (int iNeighbor = 0; iNeighbor < iNoNeighborPixel; iNeighbor++) { // for all neighbors
              for (int iSetPath = oPathFirstIndexVector[iNeighbor]; iSetPath <= oPathLastIndexVector[iNeighbor]; iSetPath++) { //for all paths through iNeighbor
                fVectorVal_MedianOnePath.clear();
                fVectorfVal_MeanOnePath.clear();
                for (unsigned int iV = 1; iV <= iPathLength; iV++) {  // caclulation Lambda and mean for path (elements 2 up to iPathLength)
                  fValPi_1 = oVol.getValue(oPathVector[iSetPath][iV].x() + iX, oPathVector[iSetPath][iV].y() + iY, oPathVector[iSetPath][iV].z() + iZ);
                  fVectorVal_MedianOnePath.push_back(fValPi_1);
                }
                //###################### MEDIAN OF ONE PATH
                medianOfVector(fVectorVal_MedianOnePath, fMedianOnePath);
                fVectorfVal_MeanOnePath.push_back(fMedianOnePath); // vector of meian of each path
              }
              //###################### MEDIAN OF PATHS THROUGHT ONE NEIGHBOOR
              double fMedianOfVector = 0.0;
              medianOfVector(fVectorfVal_MeanOnePath, fMedianOfVector);
              fVectorfVal_MeanAllPath.push_back(fMedianOfVector);        // vector of median of each neighboor
            }
            //###################### MEDIAN OF ALL NEIGHBOORS

            medianOfVector(fVectorfVal_MeanAllPath, fMedianOfVector2);
            oWalk3d.setValue(iX, iY, iZ, (fMedianOfVector2));
          }
        }
        if (bProgress) cout << ".";
      }
      break;
  }; //end swith

  oVol = oWalk3d;
}


int main(int argc, char *argv[])
{

  if (argc < 7 || argc > 8) {
    FCOUT << "Usage: volfltr <input_volume> <output_volume> <mask_width> <path_length> <path_type> <beta> [<nprocs>]" << endl;
    FCOUT <<  endl;
    FCOUT << "Parameters: " << endl;
    FCOUT <<  endl;
    FCOUT << "    mask_width:  2*n+1 (odd integer), where integer n > 0 is the minimum path length" << endl;
    FCOUT << "    path_length: n (emphasize straight paths) or >n (emphasize curved paths)" << endl;
    FCOUT << "    path_type:   0 for 4-neighborhood 2D model (image stack)" << endl;
    FCOUT << "                 1 for 8-neighborhood 2D model (image stack)" << endl;
    FCOUT << "                 2 for 6-neighborhood 3D model (3D map)" << endl;
    FCOUT << "                 3 for 26-neighborhood 3D model (3D map)" << endl;
    FCOUT << "    beta:        weighting exponent (double precision), see user guide" << endl;
    FCOUT << "    nprocs:      (optional) number of parallel threads (default: number of cores of the CPU)" << endl;
    FCOUT << "Example of use: volftr infile.situs outfile.mrc 5 2 3 0.0001" << endl;
    return 0;
  }

  int iMaskSize = atoi(argv[3]);
  int iPathLength = atoi(argv[4]);
  int iPath = atoi(argv[5]);
  Real64 fBeta = atof(argv[6]);
  if (argc == 8) {
    unsigned int iThreads;
    iThreads = atoi(argv[7]);
#ifdef _OPENMP
    omp_set_num_threads(iThreads);
    FCOUT << "Using " << iThreads << " parallel threads." << endl;
#endif
  }

  int  iKernel = 10; //atoi( argv[x] );

  //   iKernel Types:
  //   BETA_THETA,  (0)
  //   FISHER_VARIANCE_MEAN, (1)
  //   ORG_VARIANCE, (2)
  //   MAX_VARIANCE, (3)
  //   ORG_FISHER,   (4)
  //   MEDIAN,       (5)
  //   FISHER_MEDIAN, (6)
  //   FISHER_ORG_MAX_INC, (7)
  //   VARIANCE,  (8)
  //   FISHER_VARIANCE, (9)
  //   FISHER_MAX_INC_WITH_CENTER_MAX, (10) <-- Default Filter Choice.
  //   FISHER_MAX_INC_WITH_CENTER_MEAN,(11)


  FCOUT << "Loading file " << argv[1] << endl;
  svt_volume<Real64> oVol1;               // input volume
  oVol1.load(argv[1]);

  vector < vector<svt_vector3< int> > > oPathVector;    // paths filled in recursive stepGeodesic()
  vector <int> oPathFirstIndexVector;             // indexes in oPathVector of first path for a neighboor
  vector <int> oPathLastIndexVector;              // indexes in oPathVector of last path for a neighboor

  int iMargin = (int)(iMaskSize * 0.5);

  svt_vector3<int> oFirstStep1(iMargin, iMargin, iMargin);

  // fill oPathVector
  findGeodesicPaths(oPathVector,  oFirstStep1, iPathLength,  iMaskSize, iPath);

  // run filtration algorithms
  filtrGeodesicPath(oVol1, fBeta, 0.0 , iMaskSize, iPathLength,  1, iKernel, iPath);

  // save output volume
  oVol1.save(argv[2]);

  return 0;
};

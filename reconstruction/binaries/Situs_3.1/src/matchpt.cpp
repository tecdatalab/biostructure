/*********************************************************************
*                          M A T C H P T                             *
**********************************************************************
* Program is part of the Situs package URL: situs.biomachina.org     *
* (c) Stefan Birmanns and Willy Wriggers, 2009 - 2016                *
**********************************************************************
*                                                                    *
* "matchpoint": Next generation point cloud matching tool            *
*                                                                    *
**********************************************************************
* See legal statement for terms of distribution                      *
*********************************************************************/

#include "lib_mpt.h"
#ifdef _OPENMP
#include <omp.h>
#endif
#include <string.h>
#include <iostream>
#include <math.h>


/**
 * Cleanup the result stack, remove ambiguous results.
 */
void cleanResults(pointcloud &oCV_PDB, vector< matchresult > &rResults, double fThreshold)
{
  vector< matchresult > aFinal;
  vector< vec4 > aComs;
  unsigned int iNum, i, j;
  vec4 oCOA = oCV_PDB.coa();
  vec4 oTransCOA;

  iNum = rResults.size();

  for (i = 0; i < iNum; i++) {
    oTransCOA = rResults[i].getMatrix() * oCOA;

    bool bFound = false;
    for (j = 0; j < aComs.size(); j++) {
      if (aFinal[j].compareMatch(rResults[i]) == 0 || aComs[j].distance(oTransCOA) < fThreshold) {
        bFound = true;
        break;
      }
    }

    if (!bFound) {
      aComs.push_back(oTransCOA);
      matchresult oTest = rResults[i];
      aFinal.push_back(oTest);
    }
  }
  rResults = aFinal;
};


/**
 * Main routine
 */
int main(int argc, char *argv[])
{
  ///////////////////////////////////////////////////////////////////////////
  // File I/O
  ///////////////////////////////////////////////////////////////////////////

  // check command line arguments
  if (argc < 5) {
    cout << "_______________________________________________________________________________" << endl;
    MPTO << "USAGE:   matchpt file1 file2 file3 file4 [options]" << endl;
    MPTO << endl;
    MPTO << "  file1: inputfile 1, Codebook vectors from quanvol in PDB format." << endl;
    MPTO << "         Use NONE if the codebook vectors should be calculated within matchpt." << endl;
    MPTO << "         In that case matchpt will compute and match a series of vector sets and will return" << endl;
    MPTO << "         the result with the smallest RMSD. File 3 also has to be set to NONE." << endl;
    MPTO << "  file2: inputfile 2, Density map." << endl;
    MPTO << "         Use NONE if no correlation calculation desired." << endl;
    MPTO << "  file3: inputfile 3, Codebook vectors from quanpdb in PDB format." << endl;
    MPTO << "         Use NONE if the codebook vectors should be calculated within matchpt." << endl;
    MPTO << "  file4: inputfile 4, High-resolution structure in PDB format."  << endl;
    MPTO << "         Use NONE if only the codebook vectors should be matched." << endl;
    MPTO << endl;
    MPTO << "OPTIONS:" << endl;
    MPTO << endl;
    MPTO << "  -explor <int>    - Number of solutions that are 'explored' and written to disk (default: 10)." << endl;
    MPTO << "  -anchor <float>  - Radius of initial anchor point triangle search " << endl;
    MPTO << "                     space in Angstrom (default: 12A, the larger the slower)." << endl;
    MPTO << "  -radius <float>  - Radius of the neighbor-search in Angstrom " << endl;
    MPTO << "                     (default: 10A, the larger the slower)." << endl;
    MPTO << "  -wildcards <int> - Wildcards: How many unmatched points are allowed." << endl;
    MPTO << "                     To avoid false positives, it should not be larger" << endl;
    MPTO << "                     than 10% of the number of points " << endl;
    MPTO << "                     (default: 0, the larger the slower)." << endl;
    MPTO << "  -penalty <float> - Wildcard penalty in Angstrom: How much should the solutions be" << endl;
    MPTO << "                     penalized if they include unmatched points (default 1A)." << endl;
    MPTO << "  -runs <int>      - Number of runs. The algorithm will try different" << endl;
    MPTO << "                     anchor point triangles, if set to > 1 (default: 3)." << endl;
    MPTO << "  -cluster <int>   - Number of statistically independent runs used in the clustering  " << endl;
    MPTO << "                     of the points and in the determination of their variabilities (default: 8)." << endl;
    MPTO << "  -ident <float>   - Distance threshold in Anstrom for removing identical solutions." << endl;
    MPTO << "                     Useful only for oligomeric systems" << endl;
    MPTO << "                     (default 0A, the higher, the more results are filtered out)." << endl;
    MPTO << "  -res <float>     - Resolution of file2 in Angstrom (default: 15A)" << endl;
    MPTO << "  -mincv <int>     - Minimum number of vectors per structure (file4) unit (default: 4)" << endl;
    MPTO << "  -maxcv <int>     - Maximum number of vectors per structure (file4) unit (default: 9)" << endl;
    MPTO << "  -ranking <int>   - Ranking mode used for vector series / solutions:" << endl;
    MPTO << "                        0: Min RMSD              /   RMSD (default)" << endl;
    MPTO << "                        1: Min RMSD              /   cross-correlation" << endl;
    MPTO << "                        2: Min variabilities     /   RMSD" << endl;
    MPTO << "                        3: Min variabilities     /   cross-correlation" << endl;
    MPTO << "                        4: Max cross-correlation /   RMSD" << endl;
    MPTO << "                        5: Max cross-correlation /   cross-correlation" << endl;
    MPTO << "  -units <int>     - Number of structure units contained in target volume (default: 1.0)" << endl;
    MPTO << "  -nprocs <int>    - Number of parallel threads (default: the number of cores of the CPU)" << endl;
    cout << "_______________________________________________________________________________" << endl;
    return 1;
  }

  // load map cv
  pointcloud oVOL_CV;
  if (strcmp("NONE", argv[1]) != 0 && strcmp("none", argv[1]) != 0) {
    MPTO << "Loading quanvol codebook vector file: " << argv[1] << endl;
    oVOL_CV.loadPDB(argv[1]);
    if (oVOL_CV.size() == 0) {
      MPTO << "Error: Failed to load file " << argv[1] << endl;
      exit(1);
    }
  } else MPTO << "File1 == NONE" << endl;

  // load map volume
  volume oVOL;
  if (strcmp("NONE", argv[2]) != 0 && strcmp("none", argv[2]) != 0) {
    oVOL.loadVolume(argv[2]);
  } else {
    MPTO << "File2 == NONE" << endl;
    MPTO << "No low-resolution map available, will perform matching of codebook vectors, " << endl;
    MPTO << "but unable to compute correlation coefficient." << endl;
  }

  // load atomic cv
  pointcloud oPDB_CV;
  if (strcmp("NONE", argv[3]) != 0 && strcmp("none", argv[3]) != 0) {
    MPTO << "Loading quanpdb codebook vector file: " << argv[3] << endl;
    oPDB_CV.loadPDB(argv[3]);
    if (oPDB_CV.size() == 0) {
      MPTO << "Error: Failed to load file " << argv[3] << endl;
      exit(1);
    }
  } else MPTO << "File3 == NONE" << endl;

  // load atomic structure
  pointcloud oPDB;
  if (strcmp("NONE", argv[4]) != 0 && strcmp("none", argv[4]) != 0) {
    MPTO << "Loading high-resolution structure PDB file: " << argv[4] << endl;
    oPDB.loadPDB(argv[4]);
    if (oPDB.size() == 0)
      MPTO << "Error: Failed to load " << argv[4] << endl;
  } else MPTO << "File4 == NONE" << endl;
  if (oPDB.size() == 0) {
    MPTO << "No high-resolution structure available, will only perform matching of " << endl;
    MPTO << "codebook vectors." << endl;
  }

  // were codebook vectors loaded?
  bool bCV = false;
  if (oPDB_CV.size() == 0 && oVOL_CV.size() == 0) {
    if (oPDB.size() == 0 || oVOL.size() == 0) {
      MPTO << "Neither codebook vectors, nor data sets were loaded - nothing to do, exit." << endl;
      exit(1);
    }
    bCV = false;
    MPTO << "No codebook vectors available, will compute a series of vector" << endl;
    MPTO << "sets of different sizes and rank them."  << endl;
  } else if (oPDB_CV.size() != 0 && oVOL_CV.size() != 0) {
    bCV = true;
  } else {
    MPTO << "Codebook vectors are only available for one of the two data sets. Please either " << endl;
    MPTO << "specify vectors for both, the atomic model and the volumetric map, or for none of them." << endl;
    exit(1);
  }


  ///////////////////////////////////////////////////////////////////////////
  // Options
  ///////////////////////////////////////////////////////////////////////////

  unsigned int iSolutions = 10;
  double fMatchLambda = 10.0;
  double fMatchGamma = 12.0;
  unsigned int iWildcards = 0;
  double fSkipPenalty = 1.0;
  unsigned int iMatchRuns = 3;
  double fThreshold = 0.0;
  double fResolution = 15.0f;
  unsigned int iMinCV = 4;
  unsigned int iMaxCV = 9;
  double fUnits = 1.0;
  int iranking = 0;
  int icluster = 8;
  unsigned int iThreads = 1;
  
  if (argc > 5) MPTO << endl;

  // parse command-line
  for (int i = 4; i < argc; i++) {
    // Solutions
    if (strcmp(argv[i], "-explor") == 0) {
      if (i + 1 < argc) {
        iSolutions = atoi(argv[i + 1]);
        i++;
        MPTO << "Number of explored solutions set to " << iSolutions << endl;
      }
    }

    // Gamma
    if (strcmp(argv[i], "-radius") == 0) {
      if (i + 1 < argc) {
        fMatchGamma = atof(argv[i + 1]);
        i++;
        MPTO << "Neighbor-search radius set to " << fMatchGamma << " Angstrom " << endl;
      }
    }

    // Lambda
    if (strcmp(argv[i], "-anchor") == 0) {
      if (i + 1 < argc) {
        fMatchLambda = atof(argv[i + 1]);
        i++;
        MPTO << "Anchor-search radius set to " << fMatchLambda << " Angstrom " << endl;
      }
    }

    // Wildcards
    if (strcmp(argv[i], "-wildcards") == 0) {
      if (i + 1 < argc) {
        iWildcards = atoi(argv[i + 1]);
        i++;
        MPTO << "Number of wildcards set to " << iWildcards << endl;
      }
    }

    // SkipPenalty
    if (strcmp(argv[i], "-penalty") == 0) {
      if (i + 1 < argc) {
        fSkipPenalty = atof(argv[i + 1]);
        i++;
        MPTO << "Skip-Penalty set to " << fSkipPenalty << " Angstrom " << endl;
      }
    }

    // Anchor Runs
    if (strcmp(argv[i], "-runs") == 0) {
      if (i + 1 < argc) {
        iMatchRuns = atoi(argv[i + 1]);
        i++;
        MPTO << "Number of anchor point runs set to " << iMatchRuns << endl;
      }
    }

    // Cluster Runs
    if (strcmp(argv[i], "-cluster") == 0) {
      if (i + 1 < argc) {
        icluster = atoi(argv[i + 1]);
        i++;
        MPTO << "Number of codebook cluster runs set to " << icluster << endl;
      }
    }

    // Ident
    if (strcmp(argv[i], "-ident") == 0) {
      if (i + 1 < argc) {
        fThreshold = atof(argv[i + 1]);
        i++;
        MPTO << "Distance threshold set to " << fThreshold << " Angstrom " << endl;
      }
    }

    // Resolution
    if (strcmp(argv[i], "-res") == 0) {
      if (i + 1 < argc) {
        fResolution = atof(argv[i + 1]);
        i++;
        MPTO << "Resolution of target map set to " << fResolution << " Angstrom " << endl;
      }
    }

    // minimum number of vectors
    if (strcmp(argv[i], "-mincv") == 0) {
      if (i + 1 < argc) {
        iMinCV = atoi(argv[i + 1]);
        i++;
        if (iMinCV < 4) iMinCV = 4;
        MPTO << "Minimum number of vectors " << iMinCV << endl;
      }
    }

    // maximum number of vectors
    if (strcmp(argv[i], "-maxcv") == 0) {
      if (i + 1 < argc) {
        iMaxCV = atoi(argv[i + 1]);
        i++;
        if (iMaxCV < iMinCV) iMaxCV = iMinCV;
        MPTO << "Maximum number of vectors " << iMaxCV << endl;
      }
    }

    // ranking options, see above
    if (strcmp(argv[i], "-ranking") == 0) {
      if (i + 1 < argc) {
        iranking = atoi(argv[i + 1]);
        if (iranking < 0) {
          iranking = 0;
          MPTO << "Warning: -ranking must be in (0,...,5), selecting 0" << endl;
        }
        if (iranking > 5) {
          iranking = 5;
          MPTO << "Warning: -ranking must be in (0,...,5), selecting 5" << endl;
        }
        if (oVOL.getSizeX() == 0 && (iranking == 1 || iranking > 2)) {
          iranking = 0;
          MPTO << "Warning: No map loaded, unable to rank by CC, using default RMSD criterion" << endl;
        }
        if (oPDB.size() == 0 && (iranking == 1 || iranking > 2)) {
          iranking = 0;
          MPTO << "Warning: No structure loaded, unable to rank by CC, using default RMSD criterion" << endl;
        }
        i++;

        if (bCV == false) {
          switch (iranking) {
            case 0:
              MPTO << "Ranking option set to " << iranking << " (min RMSD / RMSD)" << endl;
              break;
            case 1:
              MPTO << "Ranking option set to " << iranking << " (min RMSD / cross correlation)" << endl;
              break;
            case 2:
              MPTO << "Ranking option set to " << iranking << " (min combined variabilities / RMSD)" << endl;
              break;
            case 3:
              MPTO << "Ranking option set to " << iranking << " (min combined variabilities / cross correlation)" << endl;
              break;
            case 4:
              MPTO << "Ranking option set to " << iranking << " (max cross correlation / RMSD)" << endl;
              break;
            case 5:
              MPTO << "Ranking option set to " << iranking << " (max cross correlation / cross correlation)" << endl;
              break;
          }
        } else {
          switch (iranking) {
            case 0:
            case 2:
            case 4:
              MPTO << "Ranking option set to " << iranking << " (solutions ranked by RMSD)" << endl;
              break;
            case 1:
            case 3:
            case 5:
              MPTO << "Ranking option set to " << iranking << " (solutions ranked by cross correlation)" << endl;
              break;
          }
        }
      }
    }

    // number of subunits
    if (strcmp(argv[i], "-units") == 0) {
      if (i + 1 < argc) {
        fUnits = atof(argv[i + 1]);
        i++;
        MPTO << "Number of structure units contained in target volume set to " << fUnits << endl;
      }
    }

    // number of desired threads
    if (strcmp(argv[i], "-nprocs") == 0) {
      if (i + 1 < argc) {
        iThreads = atoi(argv[i + 1]);
      } 
    }  
  }

#ifdef _OPENMP
  omp_set_num_threads(iThreads);
  if (iThreads > 1) {
    MPTO << "Using " << iThreads << " parallel threads. If multithreading is slow, or results below are unstable, use serial mode." << endl;
  } else {
    MPTO << "Using one thread (serial mode)." << endl;    	
  }
#endif


  ///////////////////////////////////////////////////////////////////////////
  // Matching
  ///////////////////////////////////////////////////////////////////////////

  char pStr[2048];
  vector< matchresult > aMatches;
  vector< double > aVOL_var;
  vector< double > aPDB_var;
  vector< pointcloud > aVOL_CV;
  vector< pointcloud > aPDB_CV;
  double CC = 0.0;

  // shall we compute a series of codebook vector sets or just use the ones the user has provided?
  if (bCV == false) {

    vec4 oCOA = oPDB.coa();
    for (unsigned int i = 0; i < oPDB.size(); i++)
      oPDB[i] -= oCOA;

    // lowest ranking score variables
    double fscore = 1.0E10;
    int iscore = -1;

    MPTO << endl;
    MPTO << "Computation and clustering of sets of (high res./low res.) codebook vectors." << endl;

#ifdef _OPENMP
    #pragma omp parallel for schedule(dynamic)
#endif
    for (int iNumCV = iMinCV; iNumCV <= (int) iMaxCV; iNumCV++) {
      // codebook vector calculation

      clustering oCluster;
      oCluster.cluster((unsigned int)(floor(iNumCV * fUnits + 0.5)), icluster, oVOL);
      pointcloud oVOL_CV = oCluster.getCodebook();
      double oVOL_var = oCluster.getVariability();
      aVOL_var.push_back(oVOL_var);
      aVOL_CV.push_back(oVOL_CV);

      oCluster.cluster(iNumCV, icluster, oPDB);
      pointcloud oPDB_CV = oCluster.getCodebook();
      double oPDB_var = oCluster.getVariability();
      aPDB_var.push_back(oPDB_var);

      // move origin to origin of the codebook vectors as the transformations will be calculated for them...
      vec4 oCOA = oPDB_CV.coa();
      for (unsigned int i = 0; i < oPDB_CV.size(); i++)
        oPDB_CV[i] -= oCOA;
      aPDB_CV.push_back(oPDB_CV);
    }

    MPTO << endl;
    MPTO << "Point cloud matching for (high res./low res.) vectors: " << endl;

#ifdef _OPENMP
    #pragma omp parallel for schedule(dynamic)
#endif
    for (int i = 0; i <= (int)(iMaxCV - iMinCV); i++) {
      // matching for range of CV numbers

      // set parameters
      aPDB_CV[i].setNextPointCOA(true);
      aPDB_CV[i].setWildcards(iWildcards);
      aPDB_CV[i].setSkipPenalty(fSkipPenalty);
      aPDB_CV[i].setRuns(iMatchRuns);
      aPDB_CV[i].setZoneSize(3);
      aPDB_CV[i].setLambda(fMatchLambda);
      aPDB_CV[i].setGamma(fMatchGamma);

      // do the matching and sort by ranking criterion
      vector< matchresult > aCurrMatches;
      aPDB_CV[i].match(aVOL_CV[i], aCurrMatches);
      cleanResults(aPDB_CV[i], aCurrMatches, fThreshold);

      double fCC_max = -1.0E10;
      for (unsigned int j = 0; j < aCurrMatches.size(); j++) {
        pointcloud currTrans = aCurrMatches[j].getMatrix() * oPDB;
        double fcurr_CC = oVOL.correlation(currTrans, fResolution);
        if (fcurr_CC > fCC_max) fCC_max = fcurr_CC;
      }

      if (aCurrMatches.size() > 0) {

        switch (iranking) { /* selection criterion for number of CV used for matching */
          case 0:
          case 1:
            if (fscore > aCurrMatches[0].getScore()) {
              fscore = aCurrMatches[0].getScore();
              iscore    = i + iMinCV;
              aMatches = aCurrMatches;
            }
            break;
          case 2:
          case 3:
            if (fscore > aPDB_var[i] + aVOL_var[i]) {
              fscore = aPDB_var[i] + aVOL_var[i];
              iscore    = i + iMinCV;
              aMatches = aCurrMatches;
            }
            break;
          case 4:
          case 5:
            if (oVOL.getSizeX() != 0) {
              if (fscore > -fCC_max) {
                fscore = -fCC_max;
                iscore  = i + iMinCV;
                aMatches = aCurrMatches;
              }
              break;
            }
          default:
            MPTO << "Error: ranking option switch statement - if() in cases 4 5 or default" << endl;
            exit(1);
        }

        if (aCurrMatches.size() == 1) sprintf(pStr, "1 match found with RMSD %5.3f A, variabilities %5.3f/%5.3f A, CC: %5.3f", aCurrMatches[0].getScore(), aPDB_var[i], aVOL_var[i], fCC_max);
        else sprintf(pStr, "%d matches found with min RMSD %5.3f A, variabilities %5.3f/%5.3f A, max CC %5.3f", (int) aCurrMatches.size(), aCurrMatches[0].getScore(), aPDB_var[i], aVOL_var[i], fCC_max);

        MPTO << i + iMinCV << "/" << (floor((i + iMinCV) * fUnits + 0.5)) << " vectors: " << pStr << endl;
      } else {
        MPTO << i + iMinCV << "/" << (floor((i + iMinCV) * fUnits + 0.5)) << " vectors: No matches found. " << endl;
      }
    }

    if (aMatches.size() > 0) {
      oPDB_CV = aPDB_CV[ iscore - iMinCV ];
      oVOL_CV = aVOL_CV[ iscore - iMinCV ];
      // move origin to origin of the codebook vectors as the transformations will be calculated for them...
      oCOA = oPDB_CV.coa();
      for (unsigned int i = 0; i < oPDB.size(); i++) oPDB[i] -= oCOA;
    } else {
      MPTO << "Warning: No matches found for full range of vectors, check input files and parameters!" << endl;
    }

  } else {

    // move origin to origin of the codebook vectors as the transformations will be calculated for them...
    vec4 oCOA = oPDB_CV.coa();
    for (unsigned int i = 0; i < oPDB.size(); i++)
      oPDB[i]    -= oCOA;
    for (unsigned int i = 0; i < oPDB_CV.size(); i++)
      oPDB_CV[i] -= oCOA;

    // set parameters
    oPDB_CV.setNextPointCOA(true);
    oPDB_CV.setWildcards(iWildcards);
    oPDB_CV.setSkipPenalty(fSkipPenalty);
    oPDB_CV.setRuns(iMatchRuns);
    oPDB_CV.setZoneSize(3);
    oPDB_CV.setLambda(fMatchLambda);
    oPDB_CV.setGamma(fMatchGamma);

    // do the matching
    oPDB_CV.match(oVOL_CV, aMatches);
    cleanResults(oPDB_CV, aMatches, fThreshold);
  }

  ///////////////////////////////////////////////////////////////////////////
  // Output Results
  ///////////////////////////////////////////////////////////////////////////

  int iRealSolutions;
  char pStrM[256], pStr1[256], pStr2[256], pStr3[256], pStr4[256];

  if (aMatches.size() > 0) {

    if (iSolutions < aMatches.size()) iRealSolutions = iSolutions;
    else iRealSolutions = aMatches.size();
    if (iRealSolutions == 1) sprintf(pStr, " solution. ");
    else sprintf(pStr, " solutions. ");
    if (aMatches.size() == 1) sprintf(pStrM, " match");
    else sprintf(pStrM, " matches");

    MPTO << endl;

    if (bCV == false) {
      sprintf(pStr1, "Based on ");
      switch (iranking) {
        case 0:
        case 1:
          sprintf(pStr2, "min RMSD");
          break;
        case 2:
        case 3:
          sprintf(pStr2, "min combined variabilities");
          break;
        case 4:
        case 5:
          sprintf(pStr2, "max CC");
          break;
      }
      sprintf(pStr3, ", selecting ");
    } else {
      sprintf(pStr1, "De");
      sprintf(pStr2, "ri");
      sprintf(pStr3, "ved ");
    }

    switch (iranking) {
      case 0:
      case 2:
      case 4:
        sprintf(pStr4, " RMSD");
        break;
      case 1:
      case 3:
      case 5:
        sprintf(pStr4, " CC");
        break;
    }
    
    MPTO << pStr1 << pStr2 << pStr3 << aMatches.size() << pStrM << " from " << oPDB_CV.size() << "/" << oVOL_CV.size() << " vectors. Exploring top " << iRealSolutions << pStr4 << pStr << endl;
    MPTO << endl;
    MPTO << "Solution filename, codebook vector RMSD in Angstrom," << endl;

    if (oPDB.size() > 0 && oVOL.getSizeX() != 0) MPTO << "cross-correlation coefficient, and permutation" << endl;
    else MPTO << "and permutation" << endl;
    MPTO << "(order of low res fitted to high res vectors):" << endl;
    MPTO << endl;

    // precompute index list for possible sorting
    vector<unsigned int> lindex;
    for (unsigned int i = 0; i < aMatches.size(); i++) lindex.push_back(i);

    // precompute CC list and partial sort index list by CC (in case of CC ranking)
    if (oPDB.size() > 0 && oVOL.getSizeX() != 0) {
      vector<double> aCClist;
      for (unsigned int i = 0; i < aMatches.size(); i++) {
        pointcloud oiTrans = aMatches[i].getMatrix() * oPDB;
        aCClist.push_back(oVOL.correlation(oiTrans, fResolution));
      }
      if (iranking == 1 || iranking == 3 || iranking == 5) {
        for (unsigned int i = 0; i < iRealSolutions; i++) {
          double CCI = aCClist[lindex[i]];
          for (unsigned int j = i + 1; j < aMatches.size(); j++) {
            double CCJ = aCClist[lindex[j]];
            if (CCI < CCJ) {
              unsigned int itmp = lindex[j];
              lindex[j] = lindex[i];
              lindex[i] = itmp;
              CCI = CCJ;
            }
          }
        }
      }
    }

    // write output following (possibly partially sorted) lindex
    for (unsigned int i = 0; i < iRealSolutions; i++) {

      // output file
      if (oPDB.size() > 0) {
        sprintf(pStr, "mpt_%03i.pdb", i + 1);
        pointcloud oTrans = aMatches[lindex[i]].getMatrix() * oPDB;
        oTrans.replacePDB(argv[4], pStr);
        if (oVOL.getSizeX() != 0) CC = oVOL.correlation(oTrans, fResolution);
      } else {
        sprintf(pStr, "mpt_%03i.log", i + 1);
        FILE *pFile = fopen(pStr, "w");
        if (pFile != NULL) {
          mat4 oMat = aMatches[lindex[i]].getMatrix();
          for (unsigned int y = 0; y < 4; y++)
            fprintf(pFile, "%05.2f %05.2f %05.2f %05.2f\n", oMat[y][0], oMat[y][1], oMat[y][2], oMat[y][3]);
          fclose(pFile);
        }
      }

      // output some information
      if (oPDB.size() > 0) {
        if (oVOL.getSizeX() != 0) sprintf(pStr, "mpt_%03i.pdb - RMSD: %6.3f  CC: %6.3f - ", i + 1, aMatches[lindex[i]].getScore(), CC);
        else sprintf(pStr, "mpt_%03i.pdb - RMSD: %6.3f - ", i + 1, aMatches[lindex[i]].getScore());
      } else {
        sprintf(pStr, "mpt_%03i.log - RMSD: %6.3f - ", i + 1, aMatches[lindex[i]].getScore());
      }
      MPTO << pStr;
      aMatches[lindex[i]].printMatch();
      cout << endl;

      // output file
      sprintf(pStr, "mpt_CV_%03i.pdb", i + 1);
      pointcloud oTrans_CV = aMatches[lindex[i]].getMatrix() * oPDB_CV;
      oTrans_CV.writePDB(pStr);
    }
    if (bCV == false) oVOL_CV.writePDB("mpt_CV_map.pdb");

  } else {
    MPTO << "No solutions found, check input files and parameters." << endl;
  }
  MPTO << endl;
  MPTO << "All done." << endl;
}

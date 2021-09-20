/*********************************************************************
*                           V O L T R A C                            *
**********************************************************************
* Program is part of the Situs package URL: situs.biomachina.org     *
* (c) Mirabela Rusu and Willy Wriggers, 2011                         *
**********************************************************************
*                                                                    *
* Volume tracing tool for helices, filaments, etc                    *
*                                                                    *
**********************************************************************
* See legal statement for terms of distribution                      *
*********************************************************************/

#include "lib_svt.h"
#include <fstream>
#include <iostream>

#define VTOUT cout << "voltrac> "

#ifdef __GAFIT_FORK
#include <sys/wait.h>
#endif

#ifdef _OPENMP
#include <omp.h> /* openmp used here only to find number of default SMP threads, since POSIX does not support a portable solution*/
#endif

enum {
  NONE,       //0
  LOCAL       //1
};



inline const char *const BoolToString(bool b)
{
  return b ? "true" : "false";
}


//
// Main routine
//
int main(int argc, char *argv[])
{

  svt_sgenrand(svt_getToD());

  //
  // Check command line arguments
  //
  if (argc < 2) {
    cout << "_______________________________________________________________________________" << endl;
    VTOUT << "USAGE:" << endl;
    VTOUT << " voltrac <density map> -<options> " << endl;
    VTOUT << endl;
    VTOUT << "OPTIONS:" << endl;
    VTOUT << endl;
    VTOUT << "\t-res <float>       Target map resolution in A [default: 8]" << endl;
    VTOUT << endl;
    VTOUT << "\t-ntraces <int>     Number of traced objects [default: 20]" << endl;
    VTOUT << endl;
    VTOUT << "\t-expth <float>     Expansion threshold as percent (values between [0 100]) " << endl;
    VTOUT << "\t                   Decrease value for noisy maps [default: 70]" << endl;
    VTOUT << endl;
    VTOUT << "\t-nprocs <int>      Number of parallel threads [default: usually number of cores of CPU]" << endl;
    VTOUT << endl;
    VTOUT << "\t-ani <float>       Resolution anisotropy factor (Z vs XY)" << endl;
    VTOUT << "\t                   ani > 1 compresses Z [default: 1]" << endl;
    VTOUT << endl;
    VTOUT << "\t-lambda <float>    Orientation-dependent density attenuation factor (Z vs XY)" << endl;
    VTOUT << "\t                   lambda < 1 attenuates Z [default: 1 ]" << endl;
    VTOUT << endl;
    VTOUT << "\t-locnorm <float>   Apply local normalization using sigma in voxel units [default: 2.5]" << endl;
    VTOUT << "\t                   0 - no local normalization is applied" << endl;
    VTOUT << endl;
    VTOUT << "\t-postgauss <float> Gaussian smoothing after any local normalization, sigma in voxels [default: 1.5]" << endl;
    VTOUT << "\t                   0 - no Gaussian smoothing is applied after any local normalization" << endl;
    VTOUT << endl;
    VTOUT << "\t-popsize <int>     Genetic algorithm population size [default: 100]" << endl;
    VTOUT << "\t                   (increase value for large maps)" << endl;
    VTOUT << endl;
    VTOUT << "\t-maxgen <int>      Maximum limit of genetic algorithm generations [default: 10000]" << endl;
    VTOUT << endl;
    VTOUT << "\t-syncgen <int>     Generations before parallel population is synchronized [default: 100]" << endl;
    VTOUT << endl;
    VTOUT << "\t-garadius <float>  Radius of the search template in A [default: 2.0]" << endl;
    VTOUT << endl;
    VTOUT << "\t-galength <int>    Length of the search template in A [default: 20]" << endl;
    VTOUT << endl;
    VTOUT << "\t-expradius <float> Radius of the expansion template in A [default: 1.0]" << endl;
    VTOUT << endl;
    VTOUT << "\t-explength <int>   Length of the expansion template in A [default: 8]" << endl;
    VTOUT << endl;
    VTOUT << "\t-distseg <float>   Distance between segments within the search or expansion templates in A [default: 1]" << endl;
    VTOUT << endl;
    VTOUT << "\t-taburad <float>   Tabu region radius in A [default: 6]" << endl;
    VTOUT << endl;
    VTOUT << "\t-expstep <float>   Translation step of the template during the expansion in A" << endl;
    VTOUT << "\t                   [default: 1.4]" <<  endl;
    VTOUT << endl;
    VTOUT << "\t-outtempl          Output the search and expansion templates [default: none] " << endl;
    VTOUT << endl;
    VTOUT << "\t-expert            Show hidden usage info for expert users" << endl;
    cout << "_______________________________________________________________________________" << endl;

    return 10;
  }

  //
  // Options
  //
  string      oMapFname; // Map Filename
  Real64      fRes; // Resolution
  // Genetic Algorithm related parameters
  int         iPopSize; // Population size
  int         eReinsertion; // Reinsertion scheme
  Real64      fMP; // Mutation probability
  Real64      fMO; // Mutation offset
  Real64      fCP; // Crossover probability
  Real64      fSP; // Selective pressure
  string      oPath; // Output path
  int         iMaxGen; // maximum number of generations (the algorithm can terminate early based on stop criteria)
  int         iMaxThread; //  maximum or total thread number (the number of running threads drops off at the end)
  int         iMaxRun; // total number of runs, not really a maximum
  Real64      fAngularStepSize;// angular step size
  Real64      fPsiFrom; // angular search range
  Real64      fPsiTo;
  Real64      fThetaFrom;
  Real64      fThetaTo;
  Real64      fPhiFrom;
  Real64      fPhiTo;
  Real64      fTranspositionProb; // transposition probability
  Real64      fDistanceThreshold; // cutoff distance
  Real64      fDistancePenalty; // cutoff distance penalty
  bool        bMutateAll;  // all mutate on / off
  Real64      fMutateAllProportion; // all mutate Proportion
  string      oLogFname; // log file
  string      oConFname; // config file
  int         iWriteModelInterval; // how often should the top score should be outputed,
  // every generation (=1) or every 10 (=10) or never (=0)
  unsigned int iTabuWindowSize; // size of tabu window - see svt_ga.h for details
  Real64      fTabuThreshold; // tabu threshold - see svt_ga.h for details
  Real64      fTabuRegionSize; // tabu region size - see svt_ga.h for details
  int         iSyncGen; // sync after how many generations?
  int         iRefinementMaxMutPerGene;// RefinementMaxMutPerGene
  //voltrac
  int         iNoOfTubes;// Number of traced objects (e.g. alpha helices) to determine
  // general template parametes
  Real64      fDistBetweenRepeats;// distance between repeats
  unsigned int iTemplatePointCount;   // the number of points in the template:
  // n-1 are in the circle and one in the center
  //expansion template parameters
  Real64      fTemplateRadius; // the radius of the template used in the expansion
  int         iTemplateRepeats; // number of times the circle is repeated in the search template
  // search template parameters
  int         iSearchTemplateRepeats; // number of times the circle is repeated
  Real64      fSearchTemplateRadius; // the radius of the template used in the ga Search
  bool        bOutputTemplates; //should the templates be output
  Real64      fCrawlingStep; // distance covered in one crawling step
  Real64      fAcceptMoveRatio; // how much of the original score is allowed to accept a move
  // ~ 0.7*OrigScore; Expansion Threshold
  unsigned int iMaxFailedCrawls; // the number of failed crawls that are accepted before stopping
  Real64      fLoNormSigma; //local normalization sigma
  Real64      fPostGaussSigma; //local normalization smoothing sigma
  Real64      fAni; // anisotropic correction
  Real64      fLambda; // orientation-dependent correction factor


  // catch -expert flag
  for (int i = 1; i < argc; i++) if (strcmp(argv[i], "-expert") == 0) {
      cout << "_______________________________________________________________________________" << endl;
      VTOUT << "EXPERT USAGE:" << endl;
      VTOUT << " voltrac -config <filename>" << endl;
      VTOUT << " or" << endl;
      VTOUT << " voltrac <density map> -<standard options> -csave <config file> " << endl;
      VTOUT << endl;
      VTOUT << "\t-config <filename> Read config file [default: none]" << endl;
      VTOUT << endl;
      VTOUT << "\t-csave <filename>  Write map name and selected options to config file [default: none]" << endl;
      cout << "_______________________________________________________________________________" << endl;
      return 11;
    }

  if (argc == 3 && (strcmp(argv[1], "-config") == 0)) {     // have config file
    svt_config oConf(argv[2]);
    VTOUT << "Reading parameters from config file: " << argv[2] << endl;
    oMapFname            = oConf.getValue("MapFilename", "");
    if (oMapFname.size() == 0) {
      VTOUT << "Error: Did not detect MapFilename entry in config file." << endl;
      exit(1);
    }
    fRes                 = oConf.getValue("Resolution", 8.0);
    iPopSize             = oConf.getValue("PopulationSize", 100);
    eReinsertion         = oConf.getValue("ReinsertionScheme", 1);
    fMP                  = oConf.getValue("MutationProbability", 0.05);
    fMO                  = oConf.getValue("MutationOffset", 0.05);
    fCP                  = oConf.getValue("CrossoverProbability", 0.95);
    fSP                  = oConf.getValue("SelectivePressure", 1.3);
    oPath                = oConf.getValue("OutputPath", "./voltrac_results");
    iMaxGen              = oConf.getValue("MaxGenerations", 10000);
    iMaxRun              = oConf.getValue("MaxRuns", 1);
    iMaxThread           = oConf.getValue("MaxThreads", 0); // set initially to zero to catch default situation
    fAngularStepSize     = oConf.getValue("AngularStepSize", 1.0);
    fPsiFrom             = oConf.getValue("PsiFrom"  , 0.0);
    fPsiTo               = oConf.getValue("PsiTo"    , 360.0);
    fThetaFrom           = oConf.getValue("ThetaFrom", 0.0);
    fThetaTo             = oConf.getValue("ThetaTo"  , 180.0);
    fPhiFrom             = oConf.getValue("PhiFrom"  , 0.0);
    fPhiTo               = oConf.getValue("PhiTo"    , 0.0);
    fTranspositionProb   = oConf.getValue("TranspositionProbability", 0.0);
    fDistanceThreshold   = oConf.getValue("DistanceThreshold", 5.0);
    fDistancePenalty     = oConf.getValue("DistancePenalty", 0.90);
    bMutateAll           = oConf.getValue("MutateAll", true);
    fMutateAllProportion = oConf.getValue("MutateAllProportion", 1.00);
    oLogFname            = oConf.getValue("LogFile", "");
    iWriteModelInterval  = oConf.getValue("WriteModelInterval", 0);
    iTabuWindowSize      = oConf.getValue("TabuWindowSize", 35);
    fTabuThreshold       = oConf.getValue("TabuThreshold", 4.0);
    fTabuRegionSize      = oConf.getValue("TabuRegionSize", 6.0);
    iSyncGen             = oConf.getValue("ThreadSyncGenerations", 100);
    iRefinementMaxMutPerGene = oConf.getValue("RefinementMaxMutPerGene", 5);
    iNoOfTubes           = oConf.getValue("NumberOfTraces", 20);
    fTemplateRadius      = oConf.getValue("TemplateRadius", 1.0);
    fSearchTemplateRadius = oConf.getValue("SearchTemplateRadius", 2.0);
    iTemplatePointCount  = oConf.getValue("TemplatePointCount", 31);
    iTemplateRepeats     = oConf.getValue("TemplateRepeats", 8);
    iSearchTemplateRepeats  = oConf.getValue("SearchTemplateRepeats", 20);
    fDistBetweenRepeats  = oConf.getValue("DistBetweenRepeats", 1.0);
    bOutputTemplates     = oConf.getValue("OutputTemplates", false);
    fCrawlingStep        = oConf.getValue("CrawlingStepSize", 1.4); // experimental value
    fAcceptMoveRatio     = oConf.getValue("AcceptMoveRatio", 0.70);
    iMaxFailedCrawls     = oConf.getValue("MaxFailedCrawls", 2);
    fLoNormSigma         = oConf.getValue("LocalNormSigma", 2.5);
    fPostGaussSigma      = oConf.getValue("PostGaussSigma", 1.5);
    fAni                 = oConf.getValue("AnisotropicCorrection", 1);
    fLambda              = oConf.getValue("Lambda", 1);

  } else { // read options from command line

    // initial default parameters
    fRes                = 8.0;
    iPopSize            = 100;
    eReinsertion        = 1;
    fMP                 = 0.05;
    fMO                 = 0.05;
    fCP                 = 0.95;
    fSP                 = 1.3;
    oPath               = "./voltrac_results";
    iMaxGen             = 10000;
    iMaxRun             = 1;
    iMaxThread          = 0;  // set initially to zero to catch default situation
    fAngularStepSize    = 1.0;
    fPsiFrom            = 0.0;
    fPsiTo              = 360.0;
    fThetaFrom          = 0.0;
    fThetaTo            = 180.0;
    fPhiFrom            = 0.0;
    fPhiTo              = 0.0;
    fTranspositionProb  = 0.0;
    fDistanceThreshold  = 5.0;
    fDistancePenalty    = 0.90;
    bMutateAll          = true;
    fMutateAllProportion = 1.00;
    oLogFname            = "";
    oConFname            = "";
    iWriteModelInterval  = 0;
    iTabuWindowSize      = 35;
    fTabuThreshold       = 4.0;
    fTabuRegionSize      = 6.0;
    iSyncGen             = 100;
    iRefinementMaxMutPerGene = 5;
    iNoOfTubes           = 20;
    fTemplateRadius      = 1.0;
    fSearchTemplateRadius = 2.0;
    iTemplatePointCount  = 31;
    iTemplateRepeats     = 8;
    iSearchTemplateRepeats  = 20;
    fDistBetweenRepeats  = 1.0;
    bOutputTemplates     = false;
    fCrawlingStep        = 1.4; // experimental value
    fAcceptMoveRatio     = 0.70;
    iMaxFailedCrawls     = 2;
    fLoNormSigma         = 2.5; // sigma used for local norm
    fPostGaussSigma      = 1.5; // sigma used for any Gaussian smoothing after local norm
    fAni                 = 1.0;
    fLambda              = 1.0;

    if (strcmp("NONE", argv[1]) != 0 && strcmp("none", argv[1]) != 0)
      oMapFname =  argv[1];

    // parse command-line
    for (int i = 2; i < argc; i++) {
      if (strcmp(argv[i], "-res") == 0)
        if (i + 1 < argc) {
          fRes = atof(argv[i + 1]);
          i++;
        }

      if (strcmp(argv[i], "-ani") == 0)
        if (i + 1 < argc) {
          fAni = atof(argv[i + 1]);
          i++;
        }

      if (strcmp(argv[i], "-lambda") == 0)
        if (i + 1 < argc) {
          fLambda = atof(argv[i + 1]);
          i++;
        }

      if (strcmp(argv[i], "-ntraces") == 0)
        if (i + 1 < argc) {
          iNoOfTubes = atoi(argv[i + 1]);
          i++;
        }

      if (strcmp(argv[i], "-expth") == 0)
        if (i + 1 < argc) {
          fAcceptMoveRatio = atof(argv[i + 1]);
          i++;
        }

      if (strcmp(argv[i], "-locnorm") == 0)
        if (i + 1 < argc) {
          fLoNormSigma = atof(argv[i + 1]);
          i++;
        }

      if (strcmp(argv[i], "-postgauss") == 0)
        if (i + 1 < argc) {
          fPostGaussSigma = atof(argv[i + 1]);
          i++;
        }

      if (strcmp(argv[i], "-garadius") == 0)
        if (i + 1 < argc) {
          fSearchTemplateRadius = atof(argv[i + 1]);
          i++;
        }

      if (strcmp(argv[i], "-galength") == 0)
        if (i + 1 < argc) {
          iSearchTemplateRepeats = atoi(argv[i + 1]);
          i++;
        }

      if (strcmp(argv[i], "-expradius") == 0)
        if (i + 1 < argc) {
          fTemplateRadius = atof(argv[i + 1]);
          i++;
        }

      if (strcmp(argv[i], "-explength") == 0)
        if (i + 1 < argc) {
          iTemplateRepeats = atoi(argv[i + 1]);
          i++;
        }

      // distance between segments
      if (strcmp(argv[i], "-distseg") == 0)
        if (i + 1 < argc) {
          fDistBetweenRepeats = atof(argv[i + 1]);
          i++;
        }

      if (strcmp(argv[i], "-expstep") == 0)
        if (i + 1 < argc) {
          fCrawlingStep = atof(argv[i + 1]);
          i++;
        }

      if (strcmp(argv[i], "-outtempl") == 0) {
        bOutputTemplates = true;
      }

      if (strcmp(argv[i], "-popsize") == 0)
        if (i + 1 < argc) {
          iPopSize = atoi(argv[i + 1]);
          i++;
        }

      if (strcmp(argv[i], "-taburad") == 0)
        if (i + 1 < argc) {
          fTabuRegionSize = atof(argv[i + 1]);
          i++;
        }

      if (strcmp(argv[i], "-maxgen") == 0)
        if (i + 1 < argc) {
          iMaxGen = atoi(argv[i + 1]);
          i++;
        }

      if (strcmp(argv[i], "-syncgen") == 0)
        if (i + 1 < argc) {
          iSyncGen = atoi(argv[i + 1]);
          i++;
        }

      if (strcmp(argv[i], "-nprocs") == 0)
        if (i + 1 < argc) {
          iMaxThread = atoi(argv[i + 1]);
          i++;
        }

      if (strcmp(argv[i], "-csave") == 0)
        if (i + 1 < argc) {
          oConFname = argv[i + 1];
          i++;
        }
    }
  }

  // load the map
  svt_ga_vol oMap;
  VTOUT << "Loading volume from: " << oMapFname.c_str() << endl;
  oMap.load(oMapFname.c_str());
  if (oMap.size() == 0) {
    VTOUT << "Error: Failed to load file " << oMapFname.c_str() << endl;
    exit(1);
  }

  Real64 fOrigZWoAni = oMap.getGridZ();

  if (fAni != 1) // some anisotropic correction needed
    oMap.interpolate_map(oMap.getWidth(), oMap.getWidth(), oMap.getWidth()*fAni);

  // create voltrac_results folder, any existing will be deleted
  char pCmd[1256];
  if (oPath.size() != 0) {
    sprintf(pCmd, "if [ -d %s ]\n then\n rm -r %s\n fi\n mkdir -p %s\n", oPath.c_str(), oPath.c_str(), oPath.c_str());
    system(pCmd);
  }
  VTOUT << "Results folder created: " << oPath << endl;

  // if specified, redirect program output to optional log file in results folder
  ofstream oFile;
  streambuf *pBuffer = NULL;
  if (oLogFname.size() != 0) {
    char pFnameOut[2560];
    sprintf(pFnameOut, "%s/%s", oPath.c_str(), oLogFname.c_str());
    VTOUT << "Redirecting output to log file in results folder: " << pFnameOut << endl;
    oFile.open(pFnameOut);
    pBuffer = cout.rdbuf();
    cout.rdbuf(oFile.rdbuf());
  }

  // set iMaxThread to default number of cores if possible, using OPENMP for portability
#ifdef _SMP_
  if (iMaxThread != 0) {
    VTOUT << "The requested number of parallel threads (populations) is: " << iMaxThread << endl;
  } else {
#ifdef _OPENMP
    int nthreads_omp, tid_omp;
    #pragma omp parallel shared(nthreads_omp) private(tid_omp)
    {
      tid_omp = omp_get_thread_num();
      if (tid_omp == 0) nthreads_omp = omp_get_num_threads();
    }
    iMaxThread = nthreads_omp;
#else
    VTOUT << "Warning: Unable to determine number of cores of CPU. Set the number of threads for parallelism." << endl;
    iMaxThread = 1;
#endif
    VTOUT << "Number of parallel threads automatically assigned by program: " << iMaxThread << endl;
  }
  if (iMaxThread <= 0) {
    VTOUT << "Error: Number of threads not valid." << endl;
    exit(1);
  }
#else
  VTOUT << "Warning: Parallel (SMP) thread support not compiled in, running in serial mode." << endl;
  iMaxThread = 1;
#endif

  VTOUT << "Resolution of target map is " << fRes << " Angstrom " << endl;
  VTOUT << "Will rank and return " << iNoOfTubes << " traces " << endl;
  VTOUT << "The expansion threshold was set to " << fAcceptMoveRatio << endl;

  if (fLoNormSigma > 0) {
    VTOUT << "A local normalization will be applied with sigma = " << fLoNormSigma << " voxels" << endl;
    if (fPostGaussSigma > 0) VTOUT << "Next, a Gaussian smoothing will be applied with sigma = " << fPostGaussSigma << " voxels" << endl;
    VTOUT << "Next, densities will be thresholded at zero, negative densities will be ignored" << endl;
  } else VTOUT << "Note: map densities will be thresholded at zero, negative densities will be ignored" << endl;

  VTOUT << "The radius of the search template was set to " << fSearchTemplateRadius <<  " Angstrom " << endl;
  VTOUT << "The length of the search template was set to " << iSearchTemplateRepeats <<  " Angstrom " << endl;
  VTOUT << "The radius of the expansion template was set to " << fTemplateRadius <<  " Angstrom " << endl;
  VTOUT << "The length of the expansion template was set to " << iTemplateRepeats <<  " Angstrom " << endl;
  VTOUT << "Output genetic algorithm and expansion templates: " << BoolToString(bOutputTemplates) << endl;
  VTOUT << "The population size was set to " << iPopSize << endl;
  VTOUT << "The radius of the tabu region was set to " << fTabuRegionSize <<  " Angstrom" << endl;
  VTOUT << "Number of independent parallel generations before synchronization: " << iSyncGen << endl;
  VTOUT << "The expansion step was set to " << fCrawlingStep <<  " Angstrom " << endl;
  VTOUT << "Anisotropic correction was set to " << fAni << endl;
  VTOUT << "Orientation dependent correction factor was set to " << fLambda << endl;

  if (fAcceptMoveRatio > 1)
    fAcceptMoveRatio /= 100.0; // if given in percent, we need the actual value that is used for multiplications

  // optional local normalization and Gaussian post processing
  if (fLoNormSigma > 0) {
    oMap.locallyNormalize(fLoNormSigma, false);
    svt_volume<Real64> oGaussian;
    if (fPostGaussSigma > 0) {
      oGaussian.createGaussian(fPostGaussSigma, 1.0);
      oMap.convolve(oGaussian, false);
    }
  }

  // always discard negative values
  oMap.threshold(0.0 , oMap.getMaxDensity());

  // create GA
  unsigned int iGenes = 4; // number of links + x,y,z coordinates for the first point
  svt_gacylinder<svt_gacylinder_ind> oGA(iGenes);
  oGA.setOutputPath(oPath.c_str());   // set output path first
  oGA.setReinsertionScheme(eReinsertion);
  oGA.setMaxGen(iMaxGen);
  oGA.setTarget(oMap);
  oGA.setNoOfCylinder2Detect(iNoOfTubes);
  oGA.setResolution(fRes);
  oGA.setMutationProb(fMP);
  oGA.setCrossoverProb(fCP);
  oGA.setMutationOffset(fMO);
  oGA.setPopSize(iPopSize);
  oGA.setSelectivePressure(fSP);
  oGA.setAngularStepSize(fAngularStepSize);
  oGA.setAngularSearchRange(fPsiFrom, fPsiTo, fThetaFrom, fThetaTo, fPhiFrom, fPhiTo);
  oGA.setTranspositionProb(fTranspositionProb);
  oGA.setDistanceThreshold(fDistanceThreshold);
  oGA.setDistanceThresholdPenalty(fDistancePenalty);
  oGA.setWriteModelInterval(iWriteModelInterval);
  oGA.setMutateAll(bMutateAll);
  oGA.setMutateAllProportion(fMutateAllProportion);
  oGA.setTabuWindowSize(iTabuWindowSize);
  oGA.setTabuThreshold(fTabuThreshold);
  oGA.setTabuRegionSize(fTabuRegionSize);
  oGA.setRefinementMaxMutPerGene(iRefinementMaxMutPerGene);
  oGA.setMaxThread(iMaxThread);
  oGA.setSyncGen(iSyncGen);
  oGA.setApplyBlurring2Model(false);
  oGA.setTemplateRadius(fTemplateRadius);
  oGA.setSearchTemplateRadius(fSearchTemplateRadius);
  oGA.setTemplatePointCount(iTemplatePointCount);
  oGA.setTemplateRepeats(iTemplateRepeats);
  oGA.setSearchTemplateRepeats(iSearchTemplateRepeats);
  oGA.setOutputTemplates(bOutputTemplates);
  oGA.setDistBetweenRepeats(fDistBetweenRepeats);
  oGA.setCrawlingStepSize(fCrawlingStep);
  oGA.setAcceptMoveRatio(fAcceptMoveRatio);
  oGA.setMaxFailedCrawls(iMaxFailedCrawls);
  oGA.setAniCorr(fAni, fOrigZWoAni);
  oGA.setLambda(fLambda);


  // if specified, output the settings in a configuration file in the results folder
  if (oConFname.size() != 0) {
    char pFnameParam[2560];
    sprintf(pFnameParam, "%s/%s", oPath.c_str(), oConFname.c_str());
    VTOUT << "Writing map filename and parameters to config file " << pFnameParam << endl;
    FILE *pFileParam = fopen(pFnameParam, "w");
    fprintf(pFileParam, "MapFilename = %s\n", oMapFname.c_str());
    fprintf(pFileParam, "LogFile = %s\n", oLogFname.c_str());
    fclose(pFileParam);
    oGA.writeConfiguration(pFnameParam);
    pFileParam = fopen(pFnameParam, "a");
    fprintf(pFileParam, "MaxRuns = %d\n", iMaxRun);
    fprintf(pFileParam, "LocalNormSigma = %f\n", fLoNormSigma);
    fprintf(pFileParam, "PostGaussSigma = %f\n", fPostGaussSigma);
    fclose(pFileParam);
  }

  // run the GA
  VTOUT << endl;
  svt_population<svt_gacylinder_ind> oPop;
  for (int iRun = 0; iRun < iMaxRun; iRun++) {
    oGA.setRun(iRun);
    oGA.setMaxThread(iMaxThread);
    oGA.setMaxGen(iMaxGen);
    oGA.setSyncGen(iSyncGen);
    oGA.setDone(false);
    if (iMaxRun != 1) VTOUT << "Statistically independent run number " << iRun + 1 << " out of " << iMaxRun << endl;
    VTOUT << "Max. generations: " << iMaxGen << "; sync generations: " << iSyncGen << "; Max. number of synchronizations: " << (int)floor(iMaxGen / iSyncGen) << endl;
    oPop = oGA.execute();
  }

  // close output redirection
  if (oLogFname.size() != 0) {
    oFile.close();
    cout.rdbuf(pBuffer);
  }

  VTOUT << "All done." << endl;
  return 0;
}

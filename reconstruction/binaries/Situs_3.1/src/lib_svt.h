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

#ifndef __SITUS_LIB_SVT
#define __SITUS_LIB_SVT

#include <vector>
#include <math.h>
#include <stdio.h>
#include <iostream>
#include <string.h>
#include <float.h>
#include <pthread.h>
#include "lib_sba.h"

#ifdef __GAFIT_FORK
#include <sys/wait.h>
#endif

#ifndef WIN32
#include <pthread.h>
#endif

#ifdef WIN32
#include <windows.h>
#include <process.h>
#endif

#define EPS 1.0E-7
/* #define PI 3.14159265358979323846 */
#define MAXPARAM 100

using namespace std;

#define svt_ga_vol svt_volume<Real64>
#define svt_ga_vec svt_vector4<Real64>
#define svt_ga_mat svt_matrix4<Real64>

#define SVTLBO cout << "lib_svt> "



//forward declarations


typedef vector< Real64 >    svt_array_real64;
typedef enum {
  UNKNOWN,        //0
  RANDOM,         //1
  CROSSOVER,      //2
  MUTATION,       //3
  MUTATIONALL,    //4
  GENEREFINEMENT, //5
  TRANSPOSITION,  //6
  TABU,           //7
  AROUND_CENTER
} creation_method;

enum {
  REINSERTION_ELITIST,
  REINSERTION_ELITIST_UNIQUE,
  REINSERTION_GLOBALRANKING,
  REINSERTION_GLOBALRANKING_UNIQUE,
  REINSERTION_SHARING
};

enum {
  SVT_THREAD_PRIORITY_NORMAL,
  SVT_THREAD_PRIORITY_HIGH,
  SVT_THREAD_PRIORITY_LOW
};

///////////////////////////////////////////////////////////////////////////////
// SVT_CONFIG
///////////////////////////////////////////////////////////////////////////////

/**
 * This class stores the application and library configuration.
 * The basic idea is that the class stores a list of keywords and values. Each keyword can only have just one value.
 * The application just ask with the functions getValue for one keyword. If no value is stored for this keyword, the class
 * returns the default value (provided by the application itself).

 * @author Stefan Birmanns
 */
class svt_config
{
  protected:

    struct parameter {
      char name[64];
      char wert[256];
    } m_aParameter[MAXPARAM];

    int m_iItems;

    void parse(const char *pFilename);
    int findItem(const char *pName);

  public:
    /**
     * Constructor
     */
    svt_config(const char *pFname = NULL);
    ~svt_config(void);

    /**
     * get boolean value
     */
    bool  getValue(const char *para, bool  default_value);
    /**
     * get int value
     */
    int   getValue(const char *para, int   default_value);
    /**
     * get float value
     */
    float getValue(const char *para, float default_value);
    /**
     * get Real64 value
     */
    Real64 getValue(const char *para, Real64 default_value);
    /**
     * get a string
     */
    const char *getValue(const char *para, const char *default_value);

    /**
     * set boolean value
     */
    void setValue(const char *para, bool new_value);
    /**
     * set int value
     */
    void setValue(const char *para, int new_value);
    /**
     * set float value
     */
    void setValue(const char *para, float new_value);
    /**
     * set Real64 value
     */
    void setValue(const char *para, Real64 new_value);
    /**
     * set a string
     */
    void setValue(const char *para, const char *new_value);
};


///////////////////////////////////////////////////////////////////////////////
// SVT_SEMAPHORE
///////////////////////////////////////////////////////////////////////////////

/** a class implementing semaphore. Semaphores are kind of mutexes, but they solve the problem of the lost-signals.
 *@author Stefan Birmanns
 */
class svt_semaphore
{
  private:
    int init;
    int s;
    int del;
#ifdef WIN32
    HANDLE mux;
#endif

#ifndef WIN32
    pthread_mutex_t mux;
    pthread_cond_t pos;
#endif
  public:
    svt_semaphore(int i = 1);
    ~svt_semaphore();

  public:
    void P();
    void V();
    bool tryLock();
};


///////////////////////////////////////////////////////////////////////////////
// SVT_CREATETHREAD
///////////////////////////////////////////////////////////////////////////////

/**
 * Create a seperate thread with the function func
 * \param pFunc pointer to the function which should be executed in the seperate thread
 * \param pArg  pointer to the arguments for the function
 * \param iPriority priority of the thread (e.g. SVT_THREAD_PRIORITY_HIGH)
 */
void svt_createThread(void *(*pFunc)(void *), void *pArg, int iPriority = SVT_THREAD_PRIORITY_NORMAL);

/**
 * Terminate thread. Has to be called by the thread itself!
 */
void svt_terminateThread();




///////////////////////////////////////////////////////////////////////////////
// SVT_GA_IND
///////////////////////////////////////////////////////////////////////////////

class svt_ga_ind;
//
/**
 * individual class for the GA algorithm
 * \author Stefan Birmanns
 */
class svt_ga_ind
{
    // the genes of the individual (always from 0.0 - 1.0, should be remapped in the fitness function)
    svt_array_real64 m_oGenes;

    // fitness or objective function
    Real64 m_fFitness;
    // probability or sometimes also called fitness function
    Real64 m_fProbability;
    // accumulated probability
    Real64 m_fAccProbability;
    //age
    unsigned int m_iAge;

    //some other property value - can be a penalty, an reward etc
    Real64 m_fProp;

    //how was the individual generated: random
    creation_method m_eOrigin;

    // the niche ID of the individual
    int m_iNicheID;

  public:

    /**
     * Constructor
     */
    svt_ga_ind() :
      m_fFitness(0.0),
      m_fProbability(0.0),
      m_fAccProbability(0.0),
      m_iAge(0),
      m_fProp(0.0),
      m_eOrigin(UNKNOWN),
      m_iNicheID(-1)
    {
    };

    /**
     * Destructor
     */
    virtual ~svt_ga_ind()
    {
    };

    /**
     * get the number of genes in this ind
     * \return number of genes
     */
    inline int getGeneCount() const
    {
      return m_oGenes.size();
    };

    /**
     * add a gene
     * \param fGene the new gene
     */
    inline void addGene(Real64 fGene)
    {
      m_oGenes.push_back(fGene);
    };

    /**
     * get the genes of this ind
     * \return vector with genes
     */
    inline svt_array_real64 getGenes() const
    {
      return m_oGenes;
    };

    /**
     * set the genes of this ind
     * \param vector with genes
     */
    inline void setGenes(svt_array_real64 oGenes)
    {
      m_oGenes = oGenes;
    };

    /**
     * get a single gene of this ind
     * \param iIndex index of the gene
     * \return int with the gene
     */
    inline Real64 getGene(int iIndex) const
    {
      return m_oGenes[iIndex];
    };

    /**
     * set a single gene of this ind
     * \param iIndex index of gene
     * \param fValue new value
     */
    inline void setGene(int iIndex, Real64 fValue)
    {
      m_oGenes[iIndex] = fValue;
    };

    /**
     * get fitness of this ind
     * \return fitness value
     */
    inline Real64 getFitness() const
    {
      return m_fFitness;
    };

    /**
     * set fitness of this ind
     * \param fFitness new fitness value
     */
    inline void setFitness(Real64 fFitness)
    {
      m_fFitness = fFitness;
    };

    /**
     * compare to inds according to their fitness
     */
    inline bool operator<(const svt_ga_ind &oInd) const
    {
      return m_fFitness < oInd.getFitness();
    };

    /**
     * compare to inds according to their fitness
     */
    inline bool operator==(const svt_ga_ind &oInd) const
    {
      return m_fFitness == oInd.getFitness();
    };

    /**
     * compare to inds according to their genes
     * \return true if both are equal, false otherwise
     */
    inline bool compareGenes(const svt_ga_ind &oInd) const
    {
      svt_array_real64 oGenes = oInd.getGenes();

      for (unsigned int i = 0; i < m_oGenes.size(); i++)
        if (m_oGenes[i] != oGenes[i])
          return false;

      return true;
    };

    /**
     * get probability of the ind to get selected during the roulette wheel process
     * \return probability
     */
    inline Real64 getProbability() const
    {
      return m_fProbability;
    };

    /**
     * set probability of the ind to get selected during the roulette wheel process
     * \param fProbability probability
     */
    inline void setProbability(Real64 fProbability)
    {
      m_fProbability = fProbability;
    };

    /**
     * get accumulated probability of the ind to get selected during the roulette wheel process
     * \return probability
     */
    inline Real64 getAccProbability() const
    {
      return m_fAccProbability;
    };

    /**
     * set accumulated probability of the ind to get selected during the roulette wheel process
     * \param fProbability probability
     */
    inline void setAccProbability(Real64 fAccProbability)
    {
      m_fAccProbability = fAccProbability;
    };

    /**
     * print genes to cout
     */
    virtual void printGenes();

    /**
     * print genes to stdout
     */
    void printGenesPf();

    /**
     * print genes to file
     */
    void printGenes(FILE *file);


    /**
     * interpret genes of individual as single number
     * \return number
     */
    inline Real64 getValue()
    {
      Real64 fValue = 0.0f;

      for (unsigned int i = 0; i < m_oGenes.size(); i++)
        fValue += m_oGenes[i];

      return fValue;

    }

    /**
     * Get the Property
     * \return the property value
     */
    inline Real64 getProp() const
    {
      return m_fProp;
    }

    /**
    * set Property
    * \param the property value
    */
    inline void setProp(Real64 fProp)
    {
      m_fProp = fProp;
    }

    /**
     * get Origin
     * \return method used to generate the individual
     */
    inline creation_method getOrigin() const
    {
      return m_eOrigin;
    }

    /**
     * set Origin
     * \return method used to generate the individual
     */
    inline void setOrigin(creation_method eOrigin)
    {
      m_eOrigin = eOrigin;
    }


    /**
     * Get the age
     * \return the age of the individue  = how many generations did it lived
     */
    inline unsigned int getAge() const
    {
      return m_iAge;
    }

    /**
     * Get the age
     * \return the age of the individue  = how many generations did it lived
     */
    inline void setAge(unsigned int iAge)
    {
      m_iAge = iAge;
    }




    /**
     * increase the age with one year;
     * \ the age of the individue  = how many generations did it lived
     */
    inline void incAge()
    {
      m_iAge++;
    }

    /**
     * make age 0 - new born individue
     *
     */
    inline void resetAge()
    {
      m_iAge = 0;
    }


    /**
     * cutoff - cuts all genes off that are outside the [0.0 .. 1.0] interval.
     */
    inline void cutoff()
    {
      for (unsigned int i = 0; i < m_oGenes.size(); i++) {
        m_oGenes[i] = fabs(m_oGenes[i]);
        if (m_oGenes[i] > 1.0f)
          m_oGenes[i] = 1.0f;
      }
    }

    /**
     * calculate the distance between two individuals
     * \param rOther reference to the other individual
     * \return vector distance between the two gene-vectors
     */
    virtual Real64 distance(svt_ga_ind &rOther);

    //
    // share fitness
    //

    /**
     * Set niche ID
     */
    inline void setNiche(int iNicheID)
    {
      m_iNicheID = iNicheID;
    };

    /**
     * Get niche ID
     */
    inline int getNiche()
    {
      return m_iNicheID;
    };

};

///////////////////////////////////////////////////////////////////////////////
// SVT_GA
///////////////////////////////////////////////////////////////////////////////

template<class T> class svt_population : public vector<T> {};
/**
 * Genetic Algorithm
 * \author Stefan Birmanns
 */
template<class T> class svt_ga
{
  protected:

    // output debug messages
    bool m_bVerbose;

    // encoding: how many genes
    int m_iGenes;

    // the current population
    svt_population<T> m_oPop;

    // the next population
    svt_population<T> m_oNextPop;

    // a temporary population
    svt_population<T> m_oTempPop;

    // the best population - holds the best m_iBestPopSize individuals in the population
    svt_population<T> m_oBestPop;

    // the population size
    int m_iPopSize;

    // the number of best individuals to rememeber at each generation
    unsigned int m_iBestPopSize;

    // how many generations so far?
    int m_iGenerations;

    // selective pressure parameter [1.0 ... 2.0]
    Real64 m_fSP;
    // crossover possibility [0.0 .. 1.0]
    Real64 m_fCrossProb;
    // mutation possibility [0.0 .. 1.0]
    Real64 m_fMutationProb;

    // fixed mutation prob, for the dynamic updating
    Real64 m_fFixedMutationProb;
    Real64 m_fFixedMutationOffset;

    //mutation offset
    Real64 m_fMutationOffset;

    // type of reinsertion scheme used
    int m_iReinsertionScheme;

    // maximum number of generations
    int m_iMaxGen;

    //number of generations before synchronize
    int m_iSyncGen;

    //some statistics about the population's fitness;
    Real64 m_fAvgFitness, m_fMinFitness,  m_fMaxFitness;
    unsigned int m_iNoUniqueInd;

    //how many time did the algo find the same best individual
    int m_iIdentBestIsSame;

    //number of interations when best individual does not change - used to stop run
    int m_iIdentBestIsSameMax;

    // minimal distance cutoff towards the top individual
    Real64 m_fCutoffDistance;
    // minimal distance cutoff penalty factor
    Real64 m_fCutoffDistancePenalty;
    // shall we apply the mutation to all the individuals?
    bool m_bMutateAll;

    // how many should be mutated
    Real64 m_fMutateAllProportion;


    // probability for the transposition to happen
    Real64 m_fTranspositionProb;

    //is the GA done
    bool m_bDone;

    // indicates whether the thread is still running - needed once m_bDone was set true from true from exterior program
    bool m_bIsThreadRunning;

    // semaphore of the thread
    svt_semaphore m_oThreadSema;

    // tabu search array with the tabu region centers
    unsigned int m_iTabuWindowSize;       // size of the tabu window over which we average the top individual distances
    svt_population<T> m_oTabus;              // the array with the tabu regions
    svt_array_real64 m_oTabuWindow;       // the tabu window with the distances of the top individuals
    Real64 m_fTabuThreshold;              // when do we consider an average distance as being too small?
    Real64 m_fTabuRegionSize;             // if the distance between an individual and a stored tabu region is smaller than this value, the individual is discarded

    // stopping criterion
    Real64 m_fStopScore;

    // current run
    unsigned int m_iRun;

    //current parallel run
    unsigned int m_iParallelRun;

    // current thread
    unsigned int m_iThread;

    //max thread
    unsigned int m_iMaxThread;

    //The size of a niche
    Real64 m_fNicheSize;

    //The maximal allowed population per niche - expressed as a proportion
    Real64 m_fMaxPopPerNiche;

    //penalty for individuals in the same niche
    Real64 m_fSameNichePenalty;

    int m_iRefinementMaxMutPerGene;

    //how many times was the fitness updated
    unsigned int m_iFitnessUpdateCount;

    //time required for the update of one generation: in seconds
    Real64 m_fTimeGen;

    //array of GA that run in parallel
    vector< svt_ga * > m_oGA_Array;

    // a parent ga - null if this is the main thread; the thread that started the parallel threads
    svt_ga *m_pParentGA;

  public:

    /**
     * Constructor
     * \param iGenes how many genes?
     */
    svt_ga(int iGenes);

    /**
     * Destructor
     */
    virtual ~svt_ga();

    /**
     * run the genetic algorithm
     * \return vector with fitness values
     */
    svt_array_real64 run();

    //
    // Population
    //

    /**
     * generate initial population
     * \param iNum number of inds in this population
     */
    virtual void initPopulation(int iNum);
    /**
     * create a new individual
     */
    virtual T initIndividual();

    /**
     * generate new population (selection, recombination, mutation, reinsertion)
     */
    virtual void updatePopulation();

    /**
     * Get the current population
     */
    svt_population<T> getPopulation();

    /**
     * get best individual
     */
    T getBest();

    /**
     * Get the tabu regions
     */
    svt_population<T> getTabuRegions();

    /**
     *  get best tabu individual
     */
    T getBestTabu();

    /**
    * Set the tabu regions
    */
    void setTabuRegions(svt_population<T> &rTabuPop);
    /**
     * Delete all tabu Regions
     */
    void delTabuRegions();
    /**
     * Set the current population
     */
    void setPopulation(svt_population<T> &rPop);

    /**
     * Update the best population
     */
    Real64 updateBestPopulation();

    /**
    * get the fitness of the best population
    */
    svt_array_real64 getBestPopFitness();

    /**
     * update fitness
     */
    void updateFitness();

    /**
     * get how many times the fitness was updated
     */
    unsigned int getFitnessUpdateCount();

    /**
     * get highest fitness
     */
    Real64 getHighestFitness();

    /**
     * print individual with highest fitness
     */
    virtual void printHighestFitness();

    /**
     * insert an individual
     */
    void insertInd(T &rInd);

    //
    // Parameters
    //

    /**
     * set population size
     * \param iPopSize number of individuals in the population
     */
    void setPopSize(int iPopSize);

    /**
     * get population size
     */
    int getPopSize();

    /**
     * set the size of the best population
     * \param iBestPopSize number of individuals in the population
     */
    void setBestPopSize(int iBestPopSize);

    /**
     * get size of the best population - population of best individuals
     */
    int getBestPopSize();

    /**
     * set maximum number of generations
     * \param iMaxGen maximum number of generations
     */
    void setMaxGen(int iMaxGen);

    /**
     * get maximum number of generations
     * \return maximum number of generations
     */
    int getMaxGen() const;

    /**
     * set the number of generations before update
     * \param iSyncGen number of generations
     */
    void setSyncGen(int iSyncGen);

    /**
    * get the number of generations before update
    * \return iSyncGen number of generations
    */
    int getSyncGen();

    /**
     * get the current generation
     * \return m_iGenerations
     */
    int getCurrGen() const;
    /**
     * set the current generation
     * \param iGenerations the new current generation index
     */
    void setCurrGen(int iGenerations);

    /**
     * set the threshold score value after which the ga will stop, even if the maximum number of generations was not reached
     * \param fStopScore ga will stop once that score is exceeded
     */
    void setStopScore(Real64 fStopScore);
    /**
     * get the threshold score value after which the ga will stop, even if the maximum number of generations was not reached
     * \return ga will stop once that score is exceeded
     */
    Real64 getStopScore();

    /**
     * set run number
     * \param iRun run number
     */
    void setRun(int iRun);

    /**
     * get the run number
     * \return run number
     */
    int getRun() const;

    /**
     * set parallel run number
     * \param iRun run number
     */
    void setParallelRun(int iParallelRun);

    /**
     * get the run number
     * \return run number
     */
    int getParallelRun() const;


    /**
     * set thread number
     * \param iThread thread number
     */
    void setThread(int iThread);

    /**
     * get the thread number
     * \return thread number
     */
    int getThread() const;

    /**
     * set max thread number
     * \param iMaxThread thread number
     */
    void setMaxThread(int iMaxThread);

    /**
     * get the max thread number
     * \return max thread number
     */
    int getMaxThread() const;

    /**
     * set the selective pressure
     * \param fSP selective pressure factor [1.0 ... 2.0] (default: 1.2)
     */
    void setSelectivePressure(Real64 fSP);
    /**
     * get the selective pressure
     * \return selective pressure factor [1.0 ... 2.0] (default: 1.2)
     */
    Real64 getSelectivePressure() const;

    /**
     * set the crossover probability
     * \param fCrossProb crossover probability [0.0 ... 1.0] (default: 0.95)
     */
    void setCrossoverProb(Real64 fCrossProb);
    /**
     * get the crossover probability
     * \return crossover probability [0.0 ... 1.0] (default: 0.95)
     */
    Real64 getCrossoverProb() const;

    /**
     * set the mutation probability
     * \param fMutationProb mutation probability [0.0 ... 1.0] (default: 0.05)
     */
    void setMutationProb(Real64 fMutationProb);
    /**
     * get the mutation probability
     * \return mutation probability [0.0 ... 1.0] (default: 0.05)
     */
    Real64 getMutationProb();

    /**
     * set the mutation offset
     * \param fMutationOffset
     */
    void setMutationOffset(Real64 fMutationOffset);
    /**
     * get the mutation offset
     * \return the mutation offset
     */
    Real64 getMutationOffset();

    /**
     * set the number of mutation that is applied on one gene durring local refinement
     */
    void setRefinementMaxMutPerGene(int iRefinementMaxMutPerGene);
    /**
    * get the number of mutation that is applied on one gene durring local refinement
    */
    int getRefinementMaxMutPerGene();


    /**
     * set the transposition probability
     * \param fTranspositionProb transposition probability [0.0 ... 1.0] (default: 0.05)
     */
    void setTranspositionProb(Real64 fTranspositionProb);
    /**
     * get the transposition probability
     * \return transposition probability [0.0 ... 1.0] (default: 0.05)
     */
    Real64 getTranspositionProb();

    /**
     * Set the cutoff distance parameter. All individuals with a distance lower than the one set here, will get penalized.
     * \param fCutoffDistance the cutoff distance
     */
    void setDistanceThreshold(Real64 fCutoffDistance);
    /**
     * Get the cutoff distance parameter. All individuals with a distance lower than the one set here, will get penalized.
     * \return the cutoff distance
     */
    Real64 getDistanceThreshold();

    /**
     * Set the cutoff distance penalty parameter. All individuals with a distance lower than the one set here, will get multiplied with the factor set here.
     * \param fCutoffDistancePenalty the cutoff distance penalty
     */
    void setDistanceThresholdPenalty(Real64 fCutoffDistancePenalty);
    /**
     * Get the cutoff distance penalty parameter. All individuals with a distance lower than the one set here, will get multiplied with the factor set here.
     * \return the cutoff distance penalty
     */
    Real64 getDistanceThresholdPenalty();

    /**
     * Shall the entire population be mutated? Normally, the GA mutates only few individuals of the new population, based on the mutationprobability. These individuals mostly are
     * new gene sets, created through crossover (crossover probability is typically very high).
     * If the parameter gets set to true, all old individuals get mutated, which can speed up the convergence of the algorithm.
     * \param bMutateAll if set to true all individuals get mutated.
     */
    void setMutateAll(bool bMutateAll);
    /**
     * Shall the entire population be mutated? Normally, the GA mutates only few individuals of the new population, based on the mutationprobability. These individuals mostly are
     * new gene sets, created through crossover (crossover probability is typically very high).
     * If the parameter gets set to true, all old individuals get mutated, which can speed up the convergence of the algorithm.
     * \return if true all individuals get mutated.
     */
    bool getMutateAll();

    /**
     * What proportion of the entire population should be mutated
     * \param bMutateAllProportion 1 means all individuals, 0 means none
     */
    void setMutateAllProportion(Real64 fMutateAllProportion);

    /**
     * What proportion of the entire population should be mutated
     * \return bMutateAllProportion 1 means all individuals, 0 means none
     */
    Real64 getMutateAllProportion();

    /**
     * get the time to compute one generation (as computed during the last generation)
     */
    Real64 getTimeGen();

    //
    // Tabu Search
    //

    /**
     * Set the tabu search window size. The tabu search computes the gene distances of the top-individual over time, with a moving window. It averages all those distance values.
     * If the distances vary a lot, because constantly completely new solutions get to the top, everything is considered fine. If the average drops, and only small differences
     * can be seen, this probably means premature convergence. The size of the window can be adjusted with this function.
     * \param iTabuWindowSize new size of the tabu-search window
     */
    void setTabuWindowSize(unsigned int iTabuWindowSize);
    /**
     * Get the tabu search window size. The tabu search computes the gene distances of the top-individual over time, with a moving window. It averages all those distance values.
     * If the distances vary a lot, because constantly completely new solutions get to the top, everything is considered fine. If the average drops, and only small differences
     * can be seen, this probably means premature convergence. The size of the window can be accessed with this function.
     * \return size of the tabu-search window
     */
    unsigned int getTabuWindowSize();

    /**
     * At some point the distances of the top individuals get really small and we consider this as stagnation of the GA. With this function one can set the threshold, if the
     * average distance is lower, we store the top individual in an array and remove all individuals from this region.
     * \fTabuThreshold the new threshold below which we say the GA stagnates
     */
    void setTabuThreshold(Real64 fTabuThreshold);
    /**
     * At some point the distances of the top individuals get really small and we consider this as stagnation of the GA. With this function one can access the threshold, if the
     * average distance is lower, we store the top individual in an array and remove all individuals from this region.
     * \return threshold below which we say the GA stagnates
     */
    Real64 getTabuThreshold();

    /**
     * If the distance between an individual and a stored tabu region is smaller than this value, the individual is discarded.
     *\param fTabuRegionSize the new size of the tabu regions
     */
    void setTabuRegionSize(Real64 fTabuRegionSize);
    /**
     * If the distance between an individual and a stored tabu region is smaller than this value, the individual is discarded.
     *\return the size of the tabu regions
     */
    Real64 getTabuRegionSize();

    /**
     * check whether the ind is in one of the tabu regions
     * \param pInd
     **/
    bool isInTabuReg(T *pInd);

    /**
     * set the parent GA
     * \param pParentGa - the ga that started this thread ; NULL if the main thread
     */
    void setParentGA(svt_ga *pParentGA);

    /**
     * get the parent GA
     * \return pParentGa - the ga that started this thread ; NULL if the main thread
     */
    svt_ga *getParentGA();

    /**
     * refine an individual;
     * \param the individual that will be refined
     */
    virtual void refineInd(T *pInd);

    //
    //  Sharing
    //

    /**
     * set the Niche size
     * \param the new nicheSize
     */
    void setNicheSize(Real64 fNicheSize);

    /**
     * get the Niche size
     * \return nicheSize
     */
    Real64 getNicheSize();

    /**
     * set the maximum allowed population per Niche - expressed as a proportion of the original population
     * Remarks: beyong this value (5-10%) - individuals have their fitness set to 0
     * \param proportion of individuals allowed in one niche
     */
    void setMaxPopPerNiche(Real64 fMaxPopPerNiche);

    /**
     * get the maximum allowed population per Niche - expressed as a proportion of the original population
     * \param proportion of individuals allowed in one niche
     */
    Real64 getMaxPopPerNiche();

    /**
     * set the Niche distance penalty - penalize individuals in the same niche according to their rank to the top individual
     * \param how much will individuals be penalized
     */
    void setSameNichePenalty(Real64 fSameNichePenalty);

    /**
     * get the Niche distance penalty - penalize individuals in the same niche according to their rank to the top individual
     * \param how much will the
     */
    Real64 getSameNichePenalty();



    //
    // Threads
    //

    /**
    * set running state of the thread;
    */
    void setIsThreadRunning(bool bIsThreadRunning);

    /**
     * \return whether the thread is running;
     */
    bool getIsThreadRunning();

    /**
     * set the variable m_bDone - (should the run of the GA stop cause it reached finish condition)
     * \param bDone - the state
     */
    void setDone(bool bDone);

    /**
     * get the variable m_bDone - (should the run of the GA stop cause it reached the finish condition?)
     * \return m_bDone - the state
     */
    bool getDone() const;

  protected:

    //
    // Selection
    //

    /**
     * selection of a new generation
     */
    void selection();

    /**
     * sort population according to fitness
     */
    void sortPopulation();
    /**
     * sort next population according to fitness
     */
    void sortNextPopulation();

    //
    // Recombination
    //

    /**
     * recombination of the selected members of the old population to form a new generation
     */
    void recombination();

    /**
     * uniform (coin-flipping) crossover operator
     * \param rParentA reference to first parent object
     * \param rParentB reference to first parent object
     * \param pNewIndA pointer to new ind
     * \param pNewIndB pointer to second new ind
     */
    virtual void crossover(T &rParentA, T &rParentB, T *pNewIndA, T *pNewIndB);

    /**
     * 1 point crossover operator
     * \param rParentA reference to first parent object
     * \param rParentB reference to first parent object
     * \param pNewIndA pointer to new ind
     * \param pNewIndB pointer to second new ind
     */
    void crossover1Point(T &rParentA, T &rParentB, T *pNewIndA, T *pNewIndB);
    /**
     * 2 point crossover operator
     * \param rParentA reference to first parent object
     * \param rParentB reference to first parent object
     * \param pNewIndA pointer to new ind
     * \param pNewIndB pointer to second new ind
     */
    void crossover2Point(T &rParentA, T &rParentB, T *pNewIndA, T *pNewIndB);
    /**
     * 2 point shuffle crossover operator
     * \param rParentA reference to first parent object
     * \param rParentB reference to first parent object
     * \param pNewIndA pointer to new ind
     * \param pNewIndB pointer to second new ind
     */
    void crossover2PointShuffle(T &rParentA, T &rParentB, T *pNewIndA, T *pNewIndB);

    /**
     * uniform (coin-flipping) crossover operator
     * \param rParentA reference to first parent object
     * \param rParentB reference to first parent object
     * \param pNewIndA pointer to new ind
     * \param pNewIndB pointer to second new ind
     */
    void crossoverUniform(T &rParentA, T &rParentB, T *pNewIndA, T *pNewIndB);
    /**
     * intermediate crossover operator
     * \param rParentA reference to first parent object
     * \param rParentB reference to first parent object
     * \param pNewIndA pointer to new ind
     * \param pNewIndB pointer to second new ind
     */
    void crossoverIntermediate(T &rParentA, T &rParentB, T *pNewIndA, T *pNewIndB);
    /**
     * subgroup crossover operator
     * \param rParentA reference to first parent object
     * \param rParentB reference to first parent object
     * \param pNewInd pointer to new ind
     */
    void crossoverSubgroup(T &rParentA, T &rParentB, T *pNewInd);
    /**
     * arithmethic crossover operator - xoff = (alpha)*x1+(1-alpha)*x2 - alpha is always random
     * \param rParentA reference to first parent object
     * \param rParentB reference to first parent object
     * \param pNewIndA pointer to new ind
     * \param pNewIndB pointer to second new ind
     */
    void crossoverArithmetic(T &rParentA, T &rParentB, T *pNewIndA, T *pNewIndB);

    //
    // Mutation
    //

    /**
     * mutation
     */
    void mutation();

    /**
     * uniform mutation
     * \param iInd index of individual
     */
    void mutationRandom(int iInd);

    /**
     * uniform mutation
     * \param iInd index of individual
     */
    void mutationUniform(int iInd);
    /**
     * moving window mutation
     * \param iInd index of individual
     */
    void mutationMovingWindow(int iInd);

    /**
     * mutation with a normal distribution
     * \param iInd index of individual
     */
    void mutationNorm(int iInd);

    /**
     * mutation with a cauchy distribution
     * \param iInd index of individual
     */
    void mutationCauchy(int iInd);

    /**
     * mutation with a cauchy distribution
     * \param oInd the individual
     */
    void mutationCauchy(T *oInd);


    /**
    * mutation with a cauchy distribution - 3 mutation per individual
    * \param iInd index of individual
    */
    void mutationMultiCauchy(int iInd);

    /**
     * uniform mutation for 1 in 7 genes
     * \param iInd index of individual
     */
    void mutationMultipoint(int iInd);

    /**
     * random mutation for 1 in 7 genes
     * \param iInd index of individual
     */
    void mutationMultiRandom(int iInd);

    /**
     * Mutate the entire population (by doubling its size).
     */
    void mutationAllPop();

    /**
     * custom mutation (can be changed by derived class, default implementation just calls mutationBGA)
     * \param iInd index of individual
     */
    virtual void mutationCustom(int iInd);

    //
    // Add new random individuals
    //

    /**
     * adds new random individuals in the population
     */
    void addRandomIndividuals();

    //
    // Transposition
    //

    /**
     * transposition
     */
    virtual void transposition();

    /**
     * flips two genes between the same individual
     * \param rParentA reference to object
     * \param pNewIndA pointer to new ind
     */
    virtual void transpositionUniform(T &rParentA, T *pNewIndA);

    //
    // Reinsertion
    //

  public:

    /**
     * set Reinsertion scheme
     * \param iReinsertionScheme the reinsertion scheme to used
     */
    void setReinsertionScheme(unsigned int iReinsertionScheme);

  protected:

    /**
     * reinsertion
     */
    void reinsertion();
    /**
     * elitist reinsertion - replace the 50% worst parents with 50% best new individuals
     */
    void reinsertion_elitist();

    /**
     * elitist reinsertion - make duplicates fitness 0; replace the 50% worst parents with 50% best new individuals
     */
    void reinsertion_elitist_unique();

    /**
     * share fitness among multiple nishes
     */
    void reinsertion_sharing();

    /**
     * share the fitness between the individuals of the population: a %percent of individuals are allowed in one niche - the rest are just killed and need to populate other niches
     * \param oPop - the individuals come form these population
     */
    void shareFitness(vector< T> &oPop);

    /**
     * reinsertion - global reinsertion based on the global ranking.
     * Both, the old and the new population are ranked together and only the best individuals are inserted into the next gen.
     */
    void reinsertion_globalranking();
    /**
     * reinsertion - global reinsertion based on the global ranking but only unique individuals are allowed
     */
    void reinsertion_globalranking_unique();

    //
    // Fitness
    //

    /**
     * update fitness - this function has to get overloaded!
     * \param pInd pointer to individual that should get updated
     */
    virtual void updateFitness(T *pInd) = 0;

    /**
     * Is the current individual a valid individual
     */
    virtual bool isValid(T *pInd);

    /**
     * Function to check verify whether the integrity of the genes is maintained
     * \param pInd the individual for which to check genes
     * \return the if individual correct
     */
    virtual void makeValid(T *pInd);


    /**
     * Create an object
     */
    virtual svt_ga *createObject() = 0;

  public:

    /**
     * Penalize individuals that are similar to allow a more diverse population
     * \param the population
     * \param fCutoffDistance the gene distance between which they get penalized
     * \param fCufoffDistancePenalty how much do they get penalized
     */
    void penalizeSimilar(svt_population<T> &oPop, Real64 fCutoffDistance, Real64 fCutoffDistancePenalty);

    /**
     * Discard invalid(fitness value=0) individuals
     * \param oPop the population
     */
    void discardNullInd(svt_population<T> &oPop);

    //
    // Print out diagnostic information
    //

    /**
     * print results (to cout)
     */
    virtual void printResults();

    /**
     * print population
     * \param the population
     */
    void printPop(svt_population<T> &oPop);

    /**
     * print results (to cout)
     */
    virtual void printNextPop();

    /**
     * Print the Min fitness, the avergage fitness and the Max fitness
     */
    void printStatistics();

    /**
     * print the fitness of each individual of the population
     */
    void printPopFitness(char *pFname);

    /**
     * output results to files
     */
    virtual void outputResult(bool bTabuAdded = false) = 0;

    /**
     * output the configuration of the program
     */
    virtual void writeConfiguration(char *pFilename);

    /**
     * Write the top scoring solutions to the disk
     * \param oPop the population of solutions
     * \param iWriteSolutions how many solutions to write
     */
    virtual void writeSolutions(svt_population<T> &oPop, unsigned int iWriteSolutions, char *pFilename) = 0;


    ///////////////////////////////////////////////////////////////////////////////
    // run ga in thread
    ///////////////////////////////////////////////////////////////////////////////

    /**
     * function to create the thread
     */
    void initThread();

    /**
     * function to create the thread
     */
    void initThreads();


    /**
     * function to create the thread
     * \return a population
     */
    virtual svt_population<T> execute();


    /**
     * Clear the content of the threads
     */
    virtual void clearThreads();

    /**
     * Refine population
     * \param oPop what population
     * \param iNoInd4Refinement how many individuals are to be refined
     */
    void refine(svt_population<T> &oPop, unsigned int iNoInd4Refinement);



};

///////////////////////////////////////////////////////////////////////////////
// Implementation
///////////////////////////////////////////////////////////////////////////////

/**
 * Constructor
 * \param iGenes how many genes?
 */
template<class T>
svt_ga<T>::svt_ga(int iGenes) :
  m_bVerbose(false),
  m_iBestPopSize(1),
  m_fSP(1.3),
  m_fCrossProb(0.95),
  m_fMutationProb(0.05),
  m_fMutationOffset(0.05f),
  m_iReinsertionScheme(REINSERTION_ELITIST_UNIQUE),
  m_iMaxGen(1000),
  m_iSyncGen(100),
  m_fCutoffDistance(0.01875),
  m_fCutoffDistancePenalty(0.90),
  m_bMutateAll(true),
  m_fMutateAllProportion(1.00),
  m_fTranspositionProb(0.05),
  m_bDone(false),
  m_iTabuWindowSize(50),
  m_fTabuThreshold(0.0015),
  m_fTabuRegionSize(0.000625),
  m_fStopScore(0.0),
  m_iRun(0),
  m_iParallelRun(0),
  m_fNicheSize(0.035),
  m_fMaxPopPerNiche(0.05),
  m_fSameNichePenalty(0.99),
  m_iRefinementMaxMutPerGene(20),
  m_iFitnessUpdateCount(0),
  m_fTimeGen(0.0f),
  m_pParentGA(NULL)
{
  m_iGenerations = 0;
  m_iGenes = iGenes;
  m_iIdentBestIsSameMax = 40;
}

/**
 * Destructor
 */
template<class T>
svt_ga<T>::~svt_ga()
{
  for (unsigned int iThread = 0; iThread < m_iMaxThread && iThread < m_oGA_Array.size() ; iThread++)
    if (m_oGA_Array[iThread] != NULL)
      delete(m_oGA_Array[iThread]);
  m_oGA_Array.clear();
}

/**
 * run the genetic algorithm
 * \return vector with fitness values
 */
template<class T>
svt_array_real64 svt_ga<T>::run()
{
  if (m_bVerbose)
    printf("svt_ga function: run() \n");

  m_bDone = false;
  svt_array_real64 oFitVec;

  // init first population (randomly)
  //initPopulation( m_iPopSize );

  m_oTempPop.reserve(m_oPop.size() * 3);
  m_oNextPop.reserve(m_oPop.size() * 2);

  long int iTime ;
  // main loop
  while (m_bDone == false) {
    //get time
    iTime = svt_getToD();

    updatePopulation();

    //stop ga if m_bDone was set true form exterior program during thread run
    if (m_bDone) return oFitVec ;

    updateBestPopulation();
    outputResult();

    // store rmsd
    oFitVec.push_back(((-1.0) * (getHighestFitness() - 1.0E10)));

    m_iGenerations++;
    if (m_iGenerations > m_iMaxGen - 1)
      m_bDone = true;

    m_fTimeGen = (Real64)(svt_getToD() - iTime) / 1000.0f;
    //SVTLBO << "Gen:" << m_iGenerations << " time: " << m_fTimeGen << endl;
  }

  return oFitVec;
}

///////////////////////////////////////////////////////////////////////////////
// population
///////////////////////////////////////////////////////////////////////////////

/**
 * generate initial population
 * \param iNum number of inds in this population
 */
template<class T>
void svt_ga<T>::initPopulation(int iNum)
{
  m_oPop.clear();
  m_oTabuWindow.clear();

  for (int i = 0; i < iNum; i++) {
    // create new object
    T oInd = initIndividual();
    m_oPop.push_back(oInd);
  }

  m_iPopSize = iNum;
  m_iGenerations = 0;
}

/**
 * create a new individual
 */
template<class T>
T svt_ga<T>::initIndividual()
{
  // create new object
  T oInd;

  // random assignments
  for (int j = 0; j < m_iGenes; j++)
    oInd.addGene(svt_genrand());

  oInd.setOrigin(RANDOM);
  oInd.resetAge();

  makeValid(&oInd);

  return oInd;
}

/**
* generate new population (selection, recombination, mutation, reinsertion)
*/
template<class T>
void svt_ga<T>::updatePopulation()
{
  m_oThreadSema.P();
  for (unsigned int i = 0; i < m_oPop.size(); i++)
    m_oPop[i].incAge();


  // only update the fitness of all individuals for the first generation - the rest is done in reinsertion
  if (m_iGenerations == 0)
    svt_ga<T>::updateFitness();

  // selection
  selection();

  // recombination
  recombination();

  // mutation
  mutation();

  // shall we mutate the entire population?
  if (m_bMutateAll)
    mutationAllPop();

  //add random individuals
  //addRandomIndividuals();

  //flip genes between the same individual
  transposition();

  // reinsertion
  reinsertion();

  // bring all genes back to the [0.0 .. 1.0] interval
  for (unsigned int i = 0; i < m_oPop.size(); i++)
    makeValid(&m_oPop[i]);

  m_oThreadSema.V();
}

/**
 * Get the current population
 */
template<class T>
svt_population<T> svt_ga<T>::getPopulation()
{
  return m_oPop;
};

/**
 * get best individual
 */
template<class T>
T svt_ga<T>::getBest()
{
  T oInd;
  if (m_oBestPop.size() > 0)
    oInd = m_oBestPop[0];

  return oInd;
};

/**
 * Get the tabu regions
 */
template<class T>
svt_population<T> svt_ga<T>::getTabuRegions()
{
  return m_oTabus;
};

/**
 *  get best tabu individual
 */
template<class T>
T  svt_ga<T>::getBestTabu()
{
  T oInd;
  if (m_oTabus.size() > 0) {
    //sort(m_oTabus.begin(), m_oTabus.end());
    oInd =  m_oTabus[ m_oTabus.size() - 1 ]; // tabus should be sorted at all time
  }

  return oInd;
};



/**
 * Set the tabu regions
 */
template<class T>
void svt_ga<T>::setTabuRegions(svt_population<T> &rTabuPop)
{
  m_oTabus = rTabuPop;
};


/**
 * Delete all tabu Regions
 */
template<class T>
void svt_ga<T>::delTabuRegions()
{
  m_oTabus.clear();
};

/**
 * Set the current population
 */
template<class T>
void svt_ga<T>::setPopulation(svt_population<T> &rPop)
{
  m_oTabuWindow.clear();
  m_oPop.clear();
  m_oPop = rPop;
  m_iPopSize = rPop.size();
};

/**
 * Update the best population
 */
template<class T>
Real64 svt_ga<T>::updateBestPopulation()
{
  m_oThreadSema.P();

  sortPopulation();

  m_oBestPop.clear();
  for (unsigned int iIndex = 0; iIndex < m_iBestPopSize; iIndex++)
    m_oBestPop.push_back(m_oPop[ m_oPop.size() - iIndex - 1 ]);

  m_oThreadSema.V();

  return m_oPop[ m_oPop.size() - 1 ].getFitness();
}


/**
 * get the fitness of the best population
 */
template<class T>
svt_array_real64 svt_ga<T>::getBestPopFitness()
{
  svt_array_real64 oVec;
  if (m_oBestPop.size() == 0) // population was not yet been initialized or no generation has been created
    return oVec;

  sort(m_oBestPop.rbegin(), m_oBestPop.rend());

  for (unsigned int iIndex = 0; iIndex < m_iBestPopSize; iIndex++)
    oVec.push_back(m_oBestPop[iIndex].getFitness());

  return oVec;
};

/**
 * get how many times the fitness was updated
 */
template<class T>
unsigned int svt_ga<T>::getFitnessUpdateCount()
{
  return m_iFitnessUpdateCount;
};


/**
 * get highest fitness
 */
template<class T>
Real64 svt_ga<T>::getHighestFitness()
{
  sortPopulation();
  return m_oPop[ m_oPop.size() - 1 ].getFitness();
}

/**
 * print individual with highest fitness
 */
template<class T>
void svt_ga<T>::printHighestFitness()
{
  sortPopulation();
  printf("%3i [%2li] = %3.5f ", m_iGenerations, m_oPop.size() - 1, (-1.0) * (getHighestFitness() - 1.0E10));
  m_oPop[m_oPop.size() - 1].printGenes();
}

///////////////////////////////////////////////////////////////////////////////
// Selection
///////////////////////////////////////////////////////////////////////////////

/**
 * selection
 * rank-based selection
 */
template<class T>
void svt_ga<T>::selection()
{
  if (m_bVerbose)
    printf("svt_ga function: selection() - in\n");

  // first step: sort the current population
  sortPopulation();

  m_fAvgFitness = 0.0f;
  m_fMaxFitness = 0.0f;
  m_fMinFitness = 1e10;
  Real64 fFitness;

  // second step: calculate the probability, based on the selective pressure parameter
  Real64 fAcc = 0.0;
  unsigned iPopSize = m_oPop.size();
  for (unsigned int i = 0; i < iPopSize; i++) {
    m_oPop[i].setProbability((2.0 - m_fSP) + 2.0 * (m_fSP  - 1.0) * i / (m_oPop.size() - 1));
    fAcc += m_oPop[i].getProbability();
    m_oPop[i].setAccProbability(fAcc);

    //compute some statistics
    fFitness = m_oPop[i].getFitness();
    m_fAvgFitness += fFitness;

    if (fFitness < m_fMinFitness && fFitness > 0.0)
      m_fMinFitness = fFitness;

    if (fFitness > m_fMaxFitness)
      m_fMaxFitness = fFitness;
  }
  m_fAvgFitness /= m_oPop.size();

  // third step: choose the individuals for the next generation according to the accumulated probabilities
  m_oNextPop.clear();
  Real64 fRand;
  unsigned int iIndex;
  for (unsigned int i = 0; i < m_oPop.size(); i++) {
    fRand = svt_genrand();
    fRand *= m_oPop.size();

    iIndex = 0;
    while (iIndex < m_oPop.size() && m_oPop[iIndex++].getAccProbability() < fRand) { };
    m_oNextPop.push_back(m_oPop[iIndex - 1]);
    m_oNextPop[ m_oNextPop.size() - 1 ].setOrigin(UNKNOWN);
  }

  //consider the TABU regions to support the building block hypothesis
  if (m_oTabus.size() > 10) {
    sort(m_oTabus.begin(), m_oTabus.end());
    //SVTLBO << m_oTabus[0].getFitness() << " " << m_oTabus[m_oTabus.size() -1 ].getFitness() << endl;
  }
  for (unsigned int i = 0; i < m_oTabus.size() && i < 10; i++) {
    m_oNextPop.push_back(m_oTabus[m_oTabus.size() - 1 - i]);
    m_oNextPop[ m_oNextPop.size() - 1 ].setOrigin(TABU);
  }

  if (m_bVerbose)
    printf("svt_ga function: selection() - out\n");
}

/**
 * sort population according to fitness
 */
template<class T>
void svt_ga<T>::sortPopulation()
{
  sort(m_oPop.begin(), m_oPop.end());
}

/**
 * sort next population according to fitness
 */
template<class T>
void svt_ga<T>::sortNextPopulation()
{
  sort(m_oNextPop.begin(), m_oNextPop.end());
}

/**
 * insert an individual
 */
template<class T>
void svt_ga<T>::insertInd(T &rInd)
{
  m_oPop.push_back(rInd);
}

///////////////////////////////////////////////////////////////////////////////
// Parameters
///////////////////////////////////////////////////////////////////////////////

/**
 * set population size
 * \param iPopSize number of individuals in the population
 */
template<class T>
void svt_ga<T>::setPopSize(int iPopSize)
{
  m_iPopSize = iPopSize;
};

/**
 * get population size
 */
template<class T>
int svt_ga<T>::getPopSize()
{
  return m_iPopSize;
};


/**
 * set the size of the best population
 * \param iBestPopSize number of individuals in the population
 */
template<class T>
void svt_ga<T>::setBestPopSize(int iBestPopSize)
{
  m_iBestPopSize = iBestPopSize;
};

/**
 * get size of the best population - population of best individuals
 */
template<class T>
int svt_ga<T>::getBestPopSize()
{
  return m_iBestPopSize;
};


/**
 * set maximum number of generations
 * \param iMaxGen maximum number of generations
 */
template<class T>
void svt_ga<T>::setMaxGen(int iMaxGen)
{
  m_iMaxGen = iMaxGen;
};
/**
 * get maximum number of generations
 * \return maximum number of generations
 */
template<class T>
int svt_ga<T>::getMaxGen() const
{
  return m_iMaxGen;
};


/**
 * set the number of generations before update
 * \param iSyncGen number of generations
 */
template<class T>
void svt_ga<T>::setSyncGen(int iSyncGen)
{
  m_iSyncGen = iSyncGen;
};

/**
* get the number of generations before update
* \return iSyncGen number of generations
*/
template<class T>
int svt_ga<T>::getSyncGen()
{
  return m_iSyncGen;
};


/**
 * get the current generation
 * \return m_iGenerations
 */
template<class T>
int svt_ga<T>::getCurrGen() const
{
  return m_iGenerations;
};
/**
 * set the current generation
 * \param iGenerations the new current generation index
 */
template<class T>
void svt_ga<T>::setCurrGen(int iGenerations)
{
  m_iGenerations = iGenerations;
};

/**
 * set the threshold score value after which the ga will stop, even if the maximum number of generations was not reached
 * \param fStopScore ga will stop once that score is exceeded
 */
template<class T>
void svt_ga<T>::setStopScore(Real64 fStopScore)
{
  m_fStopScore = fStopScore;
};
/**
 * get the threshold score value after which the ga will stop, even if the maximum number of generations was not reached
 * \return ga will stop once that score is exceeded
 */
template<class T>
Real64 svt_ga<T>::getStopScore()
{
  return m_fStopScore;
};

/**
 * set run number
 * \param iRun run number
 */
template<class T>
void svt_ga<T>::setRun(int iRun)
{
  m_iRun = iRun;
};

/**
 * get the run number
 * \return run number
 */
template<class T>
int svt_ga<T>::getRun() const
{
  return m_iRun;
};

/**
 * set run number
 * \param iRun run number
 */
template<class T>
void svt_ga<T>::setParallelRun(int iParallelRun)
{
  m_iParallelRun = iParallelRun;
};

/**
 * get the run number
 * \return run number
 */
template<class T>
int svt_ga<T>::getParallelRun() const
{
  return m_iParallelRun;
};

/**
 * set thread number
 * \param iThread thread number
 */
template<class T>
void svt_ga<T>::setThread(int iThread)
{
  m_iThread = iThread;
};

/**
 * get the thread number
 * \return thread number
 */
template<class T>
int svt_ga<T>::getThread() const
{
  return m_iThread;
};

/**
 * set max thread number
 * \param iMaxThread thread number
 */
template<class T>
void svt_ga<T>::setMaxThread(int iMaxThread)
{
  m_iMaxThread = iMaxThread;
};

/**
 * get the max thread number
 * \return max thread number
 */
template<class T>
int svt_ga<T>::getMaxThread() const
{
  return m_iMaxThread;
};

/**
 * set the selective pressure
 * \param fSP selective pressure factor [1.0 ... 2.0] (default: 1.2)
 */
template<class T>
void svt_ga<T>::setSelectivePressure(Real64 fSP)
{
  m_fSP = fSP;
};
/**
 * get the selective pressure
 * \return selective pressure factor [1.0 ... 2.0] (default: 1.2)
 */
template<class T>
Real64 svt_ga<T>::getSelectivePressure() const
{
  return m_fSP;
};

/**
 * set the crossover probability
 * \param fCrossProb crossover probability [0.0 ... 1.0] (default: 0.95)
 */
template<class T>
void svt_ga<T>::setCrossoverProb(Real64 fCrossProb)
{
  m_fCrossProb = fCrossProb;
};
/**
 * get the crossover probability
 * \return crossover probability [0.0 ... 1.0] (default: 0.95)
 */
template<class T>
Real64 svt_ga<T>::getCrossoverProb() const
{
  return m_fCrossProb;
};

/**
 * set the mutation probability
 * \param fMutationProb mutation probability [0.0 ... 1.0] (default: 0.05)
 */
template<class T>
void svt_ga<T>::setMutationProb(Real64 fMutationProb)
{
  m_fMutationProb = fMutationProb;
  m_fFixedMutationProb = fMutationProb;
};
/**
 * get the mutation probability
 * \return mutation probability [0.0 ... 1.0] (default: 0.05)
 */
template<class T>
Real64 svt_ga<T>::getMutationProb()
{
  return m_fMutationProb;
};

/**
 * set the mutation offset
 * \param fMutationOffset
 */
template<class T>
void svt_ga<T>::setMutationOffset(Real64 fMutationOffset)
{
  m_fMutationOffset = fMutationOffset;
  m_fFixedMutationOffset = fMutationOffset;
};
/**
 * get the mutation offset
 * \return the mutation offset
 */
template<class T>
Real64 svt_ga<T>::getMutationOffset()
{
  return m_fMutationOffset;
};

/**
 * set the number of mutation that is applied on one gene durring local refinement
 */
template<class T>
void svt_ga<T>::setRefinementMaxMutPerGene(int iRefinementMaxMutPerGene)
{
  m_iRefinementMaxMutPerGene = iRefinementMaxMutPerGene;
};
/**
 * get the number of mutation that is applied on one gene durring local refinement
 */
template<class T>
int svt_ga<T>::getRefinementMaxMutPerGene()
{
  return m_iRefinementMaxMutPerGene;
};

/**
 * set the transposition probability
 * \param fTranspositionProb transposition probability [0.0 ... 1.0] (default: 0.05)
 */
template<class T>
void svt_ga<T>::setTranspositionProb(Real64 fTranspositionProb)
{
  m_fTranspositionProb = fTranspositionProb;
};
/**
 * get the transposition probability
 * \return transposition probability [0.0 ... 1.0] (default: 0.05)
 */
template<class T>
Real64 svt_ga<T>::getTranspositionProb()
{
  return m_fTranspositionProb;
};

/**
 * Set the cutoff distance parameter. All individuals with a distance lower than the one set here, will get penalized.
 * \param fCutoffDistance the cutoff distance
 */
template<class T>
void svt_ga<T>::setDistanceThreshold(Real64 fCutoffDistance)
{
  m_fCutoffDistance = fCutoffDistance;
};
/**
 * Get the cutoff distance parameter. All individuals with a distance lower than the one set here, will get penalized.
 * \return the cutoff distance
 */
template<class T>
Real64 svt_ga<T>::getDistanceThreshold()
{
  return m_fCutoffDistance;
};
/**
 * Set the cutoff distance penalty parameter. All individuals with a distance lower than the one set here, will get multiplied with the factor set here.
 * \param fCutoffDistancePenalty the cutoff distance penalty
 */
template<class T>
void svt_ga<T>::setDistanceThresholdPenalty(Real64 fCutoffDistancePenalty)
{
  m_fCutoffDistancePenalty = fCutoffDistancePenalty;
};
/**
 * Get the cutoff distance penalty parameter. All individuals with a distance lower than the one set here, will get multiplied with the factor set here.
 * \return the cutoff distance penalty
 */
template<class T>
Real64 svt_ga<T>::getDistanceThresholdPenalty()
{
  return m_fCutoffDistancePenalty;
};

/**
 * Shall the entire population be mutated? Normally, the GA mutates only few individuals of the new population, based on the mutationprobability. These individuals mostly are
 * new gene sets, created through crossover (crossover probability is typically very high).
 * If the parameter gets set to true, all old individuals get mutated, which can speed up the convergence of the algorithm.
 * \param bMutateAll if set to true all individuals get mutated.
 */
template<class T>
void svt_ga<T>::setMutateAll(bool bMutateAll)
{
  m_bMutateAll = bMutateAll;
};
/**
 * Shall the entire population be mutated? Normally, the GA mutates only few individuals of the new population, based on the mutationprobability. These individuals mostly are
 * new gene sets, created through crossover (crossover probability is typically very high).
 * If the parameter gets set to true, all old individuals get mutated, which can speed up the convergence of the algorithm.
 * \return if true all individuals get mutated.
 */
template<class T>
bool svt_ga<T>::getMutateAll()
{
  return m_bMutateAll;
};

/**
 * What proportion of the entire population should be mutated
 * \param bMutateAllProportion 1 means all individuals, 0 means none
 */
template<class T>
void svt_ga<T>::setMutateAllProportion(Real64 fMutateAllProportion)
{
  m_fMutateAllProportion = fMutateAllProportion;
};

/**
 * What proportion of the entire population should be mutated
 * \return bMutateAllProportion 1 means all individuals, 0 means none
 */
template<class T>
Real64 svt_ga<T>::getMutateAllProportion()
{
  return m_fMutateAllProportion;
};

/**
 * get the time to compute one generation (as computed during the last generation)
 */
template<class T>
Real64 svt_ga<T>::getTimeGen()
{
  return m_fTimeGen;
};

/**
 * set the variable m_bDone - (should the run of the GA stop cause it reached finish condition)
 * \param bDone - the state
 */
template<class T>
void svt_ga<T>::setDone(bool bDone)
{
  m_bDone = bDone;

  for (unsigned int iThread = 0; iThread < m_oGA_Array.size(); iThread++)
    m_oGA_Array[iThread]->setDone(bDone);
}

/**
 * get the variable m_bDone - (should the run of the GA stop cause it reached the finish condition?)
 * \return m_bDone - the state
 */
template<class T>
bool svt_ga<T>::getDone() const
{
  return m_bDone;
}

/**
 * set running state of the thread;
 */
template<class T>
void svt_ga<T>::setIsThreadRunning(bool bIsThreadRunning)
{
  m_bIsThreadRunning = bIsThreadRunning;

  for (unsigned int iThread = 0; iThread < m_oGA_Array.size(); iThread++)
    m_oGA_Array[iThread]->setIsThreadRunning(bIsThreadRunning);
};


/**
 * \return whether the thread is running;
 */
template<class T>
bool svt_ga<T>::getIsThreadRunning()
{
  return m_bIsThreadRunning;
};


///////////////////////////////////////////////////////////////////////////////
// Recombination
///////////////////////////////////////////////////////////////////////////////

/**
 * recombination
 */
template<class T>
void svt_ga<T>::recombination()
{
  if (m_bVerbose)
    printf("svt_ga function: recombination() - in\n");

  // select some individuals from the new population for recombination
  Real64 fRand;
  for (unsigned int i = 0; i < m_oNextPop.size(); i++) {
    fRand = svt_genrand();

    if (fRand < m_fCrossProb) {
      int iInd1, iInd2;
      T oInd1, oInd2;

      fRand = svt_genrand();
      iInd1 = (int)(fRand * (Real64)(m_oNextPop.size()));

      do {
        fRand = svt_genrand();
        iInd2 = (int)(fRand * (Real64)(m_oNextPop.size()));
      } while (iInd1 == iInd2);

      // choose randomly which crossover operator we want to apply
      fRand = svt_genrand();

      if (fRand < 0.16)
        crossover1Point(m_oNextPop[iInd1], m_oNextPop[iInd2], &oInd1, &oInd2);
      else if (fRand < 0.32)
        crossover2Point(m_oNextPop[iInd1], m_oNextPop[iInd2], &oInd1, &oInd2);
      else if (fRand < 0.48)
        crossover2PointShuffle(m_oNextPop[iInd1], m_oNextPop[iInd2], &oInd1, &oInd2);
      else if (fRand < 0.64)
        crossoverArithmetic(m_oNextPop[iInd1], m_oNextPop[iInd2], &oInd1, &oInd2);
      else if (fRand < 0.80)
        crossoverIntermediate(m_oNextPop[iInd1], m_oNextPop[iInd2], &oInd1, &oInd2);
      else if (fRand < 1.0)
        crossoverUniform(m_oNextPop[iInd1], m_oNextPop[iInd2], &oInd1, &oInd2);

      m_oNextPop[iInd1] = oInd1;
      m_oNextPop[iInd2] = oInd2;

      m_oNextPop[iInd1].setOrigin(CROSSOVER);
      m_oNextPop[iInd1].resetAge();
      makeValid(&m_oNextPop[iInd1]);
      m_oNextPop[iInd2].setOrigin(CROSSOVER);
      m_oNextPop[iInd2].resetAge();
      makeValid(&m_oNextPop[iInd2]);
    }


  }

  int iDel = 0;
  for (unsigned int i = 0; i < m_oNextPop.size(); i++) {
    if (m_oNextPop[i].getOrigin() == TABU) {
      m_oNextPop.erase(m_oNextPop.begin() + i);
      i--;
      iDel ++;
    }
  }

  for (unsigned int i = 0; i < m_oNextPop.size(); i++) {
    if (m_oNextPop[i].getOrigin() == TABU) {
      SVTLBO << "Found tabu, but I shouldn't " << endl;
    }
  }



  if (m_bVerbose)
    printf("svt_ga function: recombination() - out\n");

}

/**
 * uniform (coin-flipping) crossover operator
 * \param rParentA reference to first parent object
 * \param rParentB reference to first parent object
 * \param pNewIndA pointer to new ind
 * \param pNewIndB pointer to second new ind
 */
template<class T>
void svt_ga<T>::crossover(T &rParentA, T &rParentB, T *pNewIndA, T *pNewIndB)
{
  crossoverUniform(rParentA, rParentB, pNewIndA, pNewIndB);
};

/**
 * 2 point crossover operator
 * \param rParentA reference to first parent object
 * \param rParentB reference to first parent object
 * \param pNewIndA pointer to new ind
 * \param pNewIndB pointer to second new ind
 */
template<class T>
void svt_ga<T>::crossover2Point(T &rParentA, T &rParentB, T *pNewIndA, T *pNewIndB)
{
  svt_array_real64 oGenesA = rParentA.getGenes();
  svt_array_real64 oGenesB = rParentB.getGenes();
  svt_array_real64 oNewGenA;
  svt_array_real64 oNewGenB;
  int i, iPoint1, iPoint2;
  Real64 fRand;

  fRand = svt_genrand();
  iPoint1 = (int)(fRand * (Real64)(m_iGenes));
  do {
    fRand = svt_genrand();
    iPoint2 = (int)(fRand * (Real64)(m_iGenes));
  } while (iPoint1 == iPoint2);

  if (iPoint1 > iPoint2) {
    int iTemp = iPoint1;
    iPoint1 = iPoint2;
    iPoint2 = iTemp;
  }

  for (i = 0; i < iPoint1; i++) {
    oNewGenA.push_back(oGenesA[i]);
    oNewGenB.push_back(oGenesB[i]);
  }

  for (i = iPoint1; i < iPoint2; i++) {
    oNewGenA.push_back(oGenesB[i]);
    oNewGenB.push_back(oGenesA[i]);
  }

  for (i = iPoint2; i < m_iGenes; i++) {
    oNewGenA.push_back(oGenesA[i]);
    oNewGenB.push_back(oGenesB[i]);
  }

  pNewIndA->setGenes(oNewGenA);
  pNewIndB->setGenes(oNewGenB);
};

/**
 * 2 point shuffle crossover operator
 * \param rParentA reference to first parent object
 * \param rParentB reference to first parent object
 * \param pNewIndA pointer to new ind
 * \param pNewIndB pointer to second new ind
 */
template<class T>
void svt_ga<T>::crossover2PointShuffle(T &rParentA, T &rParentB, T *pNewIndA, T *pNewIndB)
{
  svt_array_real64 oGenesA = rParentA.getGenes();
  svt_array_real64 oGenesB = rParentB.getGenes();
  svt_array_real64 oNewGenA;
  svt_array_real64 oNewGenB;
  Real64 fRand;

  // left/right shuffle of parent a
  fRand = svt_genrand();
  if (fRand < 0.5) {
    for (int i = 1; i < m_iGenes; i++)
      oNewGenA.push_back(rParentA.getGene(i));
    oNewGenA.push_back(rParentA.getGene(0));
  } else {
    oNewGenA.push_back(rParentA.getGene(m_iGenes - 1));
    for (int i = 0; i < m_iGenes - 1; i++)
      oNewGenA.push_back(rParentA.getGene(i));
  }

  // left/right shuffle of parent b
  fRand = svt_genrand();
  if (fRand < 0.5) {
    for (int i = 1; i < m_iGenes; i++)
      oNewGenB.push_back(rParentB.getGene(i));
    oNewGenB.push_back(rParentB.getGene(0));
  } else {
    oNewGenB.push_back(rParentB.getGene(m_iGenes - 1));
    for (int i = 0; i < m_iGenes - 1; i++)
      oNewGenB.push_back(rParentB.getGene(i));
  }

  // now let us do the crossover
  T oIndA;
  oIndA.setGenes(oNewGenA);
  T oIndB;
  oIndB.setGenes(oNewGenB);
  crossover2Point(oIndA, oIndB, pNewIndA, pNewIndB);
};

/**
 * 1 point crossover operator
 * \param rParentA reference to first parent object
 * \param rParentB reference to first parent object
 * \param pNewIndA pointer to new ind
 * \param pNewIndB pointer to second new ind
 */
template<class T>
void svt_ga<T>::crossover1Point(T &rParentA, T &rParentB, T *pNewIndA, T *pNewIndB)
{
  svt_array_real64 oGenesA = rParentA.getGenes();
  svt_array_real64 oGenesB = rParentB.getGenes();
  svt_array_real64 oNewGenA;
  svt_array_real64 oNewGenB;

  Real64 fRand = svt_genrand();
  int iPoint1 = (int)(fRand * (Real64)(m_iGenes));

  if (svt_genrand() > 0.5) {
    for (int i = 0; i < iPoint1; i++) {
      oNewGenA.push_back(oGenesA[i]);
      oNewGenB.push_back(oGenesB[i]);
    }

    for (int i = iPoint1; i < m_iGenes; i++) {
      oNewGenA.push_back(oGenesB[i]);
      oNewGenB.push_back(oGenesA[i]);
    }
  } else {
    for (int i = 0; i < iPoint1; i++) {
      oNewGenA.push_back(oGenesB[i]);
      oNewGenB.push_back(oGenesA[i]);
    }

    for (int i = iPoint1; i < m_iGenes; i++) {
      oNewGenA.push_back(oGenesA[i]);
      oNewGenB.push_back(oGenesB[i]);
    }
  }

  pNewIndA->setGenes(oNewGenA);
  pNewIndB->setGenes(oNewGenB);
};

/**
 * uniform (coin-flipping) crossover operator
 * \param rParentA reference to first parent object
 * \param rParentB reference to first parent object
 * \param pNewIndA pointer to new ind
 * \param pNewIndB pointer to second new ind
 */
template<class T>
void svt_ga<T>::crossoverUniform(T &rParentA, T &rParentB, T *pNewIndA, T *pNewIndB)
{
  for (int j = 0; j < m_iGenes; j++) {
    Real64 fRand = svt_genrand();

    if (fRand < 0.5) {
      pNewIndA->addGene(rParentA.getGene(j));
      pNewIndB->addGene(rParentB.getGene(j));
    } else {
      pNewIndA->addGene(rParentB.getGene(j));
      pNewIndB->addGene(rParentA.getGene(j));
    }
  }
}

/**
 * intermediate crossover operator
 * \param rParentA reference to first parent object
 * \param rParentB reference to first parent object
 * \param pNewIndA pointer to new ind
 * \param pNewIndB pointer to second new ind
 */
template<class T>
void svt_ga<T>::crossoverIntermediate(T &rParentA, T &rParentB, T *pNewIndA, T *pNewIndB)
{
  for (int j = 0; j < m_iGenes; j++) {
    Real64 fRand = svt_genrand();
    fRand *= 1.5;
    fRand -= 0.25;

    pNewIndA->addGene(rParentA.getGene(j) + (fRand * (rParentB.getGene(j) - rParentA.getGene(j))));
    pNewIndB->addGene(rParentB.getGene(j) + (fRand * (rParentA.getGene(j) - rParentB.getGene(j))));
  }
}

/**
 * arithmethic crossover operator - xoff = (alpha)*x1+(1-alpha)*x2 - alpha is always random
 * \param rParentA reference to first parent object
 * \param rParentB reference to first parent object
 * \param pNewIndA pointer to new ind
 * \param pNewIndB pointer to second new ind
 */
template<class T>
void svt_ga<T>::crossoverArithmetic(T &rParentA, T &rParentB, T *pNewIndA, T *pNewIndB)
{
  for (int j = 0; j < m_iGenes; j++) {
    Real64 fRand = svt_genrand();

    pNewIndA->addGene(fRand * rParentA.getGene(j) + (1.0f - fRand)*rParentB.getGene(j));
    pNewIndB->addGene((1 - fRand)*rParentA.getGene(j) + (fRand)*rParentB.getGene(j));

  }
};

/**
 * subgroup crossover operator
 * \param rParentA reference to first parent object
 * \param rParentB reference to first parent object
 * \param pNewInd pointer to new ind
 */
template<class T>
void svt_ga<T>::crossoverSubgroup(T &rParentA, T &rParentB, T *pNewInd)
{
  svt_array_real64 oGenesA = rParentA.getGenes();
  svt_array_real64 oGenesB = rParentB.getGenes();
  svt_array_real64 oNewGen;
  int i, iPoint1, iPoint2, iPoint3;
  Real64 fRand;

  // first step: cut out a part from ind A
  fRand = svt_genrand();
  iPoint1 = (int)(fRand * (Real64)(m_iGenes));
  do {
    fRand = svt_genrand();
    iPoint2 = (int)(fRand * (Real64)(m_iGenes));
  } while (iPoint1 == iPoint2);

  if (iPoint1 > iPoint2) {
    int iTemp = iPoint1;
    iPoint1 = iPoint2;
    iPoint2 = iTemp;
  }

  for (i = iPoint1; i < iPoint2; i++)
    oNewGen.push_back(oGenesA[i]);

  // now cut out a part from ind B
  do {
    fRand = svt_genrand();
    iPoint3 = (int)(fRand * (Real64)(m_iGenes));
  } while (iPoint3 > iPoint2 - iPoint1);

  int iEnd = iPoint3 + (m_iGenes - (iPoint2 - iPoint1));

  //cout << "iPoint1: " << iPoint1 << " iPoint2: " << iPoint2 << " iPoint3: " << iPoint3 << " iEnd: " << iEnd << endl;
  for (i = iPoint3; i < iEnd; i++)
    oNewGen.push_back(oGenesB[i]);

  pNewInd->setGenes(oNewGen);
};

///////////////////////////////////////////////////////////////////////////////
// Mutation
///////////////////////////////////////////////////////////////////////////////

/**
 * mutation
 */
template<class T>
void svt_ga<T>::mutation()
{
  if (m_bVerbose)
    printf("svt_ga function: mutation() - in \n");

  Real64 fRand;
  for (unsigned int i = 0; i < m_oNextPop.size(); i++) {
    fRand = svt_genrand();
    if (fRand < m_fMutationProb || m_oNextPop[i].getOrigin() != CROSSOVER) { // mutate if individual was not modified through cross over
      mutationCustom(i);
      m_oNextPop[i].setOrigin(MUTATION);
      m_oNextPop[i].resetAge();
      makeValid(&m_oNextPop[i]);
    }
  }

  if (m_bVerbose)
    printf("svt_ga function: mutation() - out \n");

}

/**
 * uniform mutation
 * \param iInd index of individual
 */
template<class T>
void svt_ga<T>::mutationRandom(int iInd)
{
  // select gene for mutation svt_genrandomly
  Real64 fRand = svt_genrand();
  int iRandIndex = (int)(fRand * (Real64)(m_iGenes));

  m_oNextPop[iInd].setGene(iRandIndex, svt_genrand());
}

/**
 * uniform mutation
 * \param iInd index of individual
 */
template<class T>
void svt_ga<T>::mutationUniform(int iInd)
{
  // select gene for mutation svt_genrandomly
  Real64 fRand = svt_genrand();
  int iRandIndex = (int)(fRand * (Real64)(m_iGenes));

  // generate offset
  fRand = svt_genrand() * 0.05;

  if (svt_genrand() > 0.5)
    fRand *= -1.0;

  m_oNextPop[iInd].setGene(iRandIndex, m_oNextPop[iInd].getGene(iRandIndex) + fRand);
}

/**
 * mutation with a normal distribution
 * \param iInd index of individual
 */
template<class T>
void svt_ga<T>::mutationNorm(int iInd)
{
  // select gene for mutation svt_genrandomly
  Real64 fRand = svt_genrand();
  int iRandIndex = (int)(fRand * (Real64)(m_iGenes));

  // generate offset
  fRand = svt_ranNormal(0.0, m_fMutationOffset);

  m_oNextPop[iInd].setGene(iRandIndex, m_oNextPop[iInd].getGene(iRandIndex) + fRand);
};

/**
 * mutation with a cauchy distribution
 * \param iInd index of individual
 */
template<class T>
void svt_ga<T>::mutationCauchy(int iInd)
{
  // select gene for mutation svt_genrandomly
  Real64 fRand = svt_genrand();
  int iRandIndex = (int)(fRand * (Real64)(m_iGenes));
  Real64 fIntPart;

  // generate offset
  fRand = svt_ranCauchy(0.0, m_fMutationOffset);
  fRand = modf(fRand, &fIntPart); // keep only the fractional part of the Random number

  m_oNextPop[iInd].setGene(iRandIndex, m_oNextPop[iInd].getGene(iRandIndex) + fRand);
};

/**
 * mutation with a cauchy distribution
 * \param oInd the individual
 */
template<class T>
void svt_ga<T>::mutationCauchy(T *oInd)
{
  // select gene for mutation svt_genrandomly
  Real64 fRand = svt_genrand();
  int iRandIndex = (int)(fRand * (Real64)(m_iGenes));
  Real64 fIntPart;

  // generate offset
  fRand = svt_ranCauchy(0.0, m_fMutationOffset);
  fRand = modf(fRand, &fIntPart); // keep only the fractional part of the Random number

  oInd->setGene(iRandIndex, oInd->getGene(iRandIndex) + fRand);
};

/**
 * mutation with a cauchy distribution - 3 mutation per individual
 * \param iInd index of individual
 */
template<class T>
void svt_ga<T>::mutationMultiCauchy(int iInd)
{
  Real64 fRand;
  Real64 fMutRatio = 4.0f; // 1 in iMutRatio will be mutated
  unsigned int iGenes = m_oNextPop[iInd].getGeneCount();
  Real64 fMutCount = iGenes / fMutRatio; // how many mutations per individual

  Real64 fIntPart;

  for (unsigned int iIndex = 0; iIndex < iGenes; iIndex++) {
    fRand = svt_genrand();

    if (fRand < 1.0 / fMutRatio) { // do mutate - mutate 1 genes in 7
      fRand = svt_ranCauchy(0.0, m_fMutationOffset / fMutCount);
      fRand = modf(fRand, &fIntPart); // keep only the fractional part of the Random number

      m_oNextPop[iInd].setGene(iIndex, m_oNextPop[iInd].getGene(iIndex) + fRand);
    }
  }

};

/**
 * uniform mutation for 1 in 7 genes
 * \param iInd index of individual
 */
template<class T>
void svt_ga<T>::mutationMultipoint(int iInd)
{
  Real64 fRand;
  Real64 fMutRatio = 4.0f; // 1 in iMutRatio will be mutated
  unsigned int iGenes = m_oNextPop[iInd].getGeneCount();
  Real64 fMutCount = iGenes / fMutRatio; // how many mutations per individual

  for (unsigned int iIndex = 0; iIndex < iGenes; iIndex++) {
    fRand = svt_genrand();

    if (fRand < 1.0 / fMutRatio) { // do mutate - mutate 1 genes in 7
      fRand = svt_genrand() *  m_fMutationOffset / fMutCount;
      m_oNextPop[iInd].setGene(iIndex, m_oNextPop[iInd].getGene(iIndex) + fRand);
    }
  }

};

/**
 * random mutation for 1 in 7 genes
 * \param iInd index of individual
 */
template<class T>
void svt_ga<T>::mutationMultiRandom(int iInd)
{
  Real64 fRand;
  Real64 fMutRatio = 4.0f; // 1 in iMutRatio will be mutated
  unsigned int iGenes = m_oNextPop[iInd].getGeneCount();

  for (unsigned int iIndex = 0; iIndex < iGenes; iIndex++) {
    fRand = svt_genrand();

    if (fRand < 1.0 / fMutRatio) { // do mutate - mutate 1 genes in 7
      fRand = svt_genrand();
      m_oNextPop[iInd].setGene(iIndex, fRand);
    }
  }

};

/**
 * moving window mutation
 * \param iInd index of individual
 */
template<class T>
void svt_ga<T>::mutationMovingWindow(int iInd)
{
  // select gene for mutation svt_genrandomly
  Real64 fRand = svt_genrand();
  int iRandIndex = (int)(fRand * (Real64)(m_iGenes));

  // generate offset
  fRand = svt_genrand();
  fRand /= 2.0;
  fRand += (1.0f - (m_iGenerations / m_iMaxGen)) * 0.5;

  if (svt_genrand() > 0.5)
    fRand *= -1.0;

  m_oNextPop[iInd].setGene(iRandIndex, m_oNextPop[iInd].getGene(iRandIndex) + fRand);
}

/**
 * custom mutation (can be changed by derived class, default implementation just calls mutationUniform)
 * \param iInd index of individual
 */
template<class T>
void svt_ga<T>::mutationCustom(int iInd)
{
  mutationCauchy(iInd);
}

/**
 * Mutate the entire population (by doubling its size).
 */
template<class T>
void svt_ga<T>::mutationAllPop()
{
  if (m_bVerbose)
    printf("svt_ga function: mutationAllPop() - in\n");

  sortPopulation();

  unsigned int iSize = m_oPop.size();
  Real64 fDist;
  //SVTLBO << " MutateAll between" << iSize -1 << " : " << iSize-1-iSize*m_fMutateAllProportion+1 << endl;
  for (int i = iSize - 1; i > iSize - 1 - iSize * m_fMutateAllProportion + 1 && i >= 0; i--) {
    m_oNextPop.push_back(m_oPop[i]);
    mutationCustom(m_oNextPop.size() - 1);
    m_oNextPop[ m_oNextPop.size() - 1].setOrigin(MUTATIONALL);
    m_oNextPop[ m_oNextPop.size() - 1].resetAge();

    //check if the new individuals is "identical" with its parent  and discard then
    fDist = ((svt_ga_ind)m_oPop[i]).distance((svt_ga_ind &)(m_oNextPop[ m_oNextPop.size() - 1 ]));
    if (fDist < EPS)
      m_oNextPop.pop_back();
  }

  if (m_bVerbose)
    printf("svt_ga function: mutationAllPop() - out\n");

}

///////////////////////////////////////////////////////////////////////////////
// Add random individuals
///////////////////////////////////////////////////////////////////////////////

/**
 * adds new random individuals in the population
 */
template<class T>
void svt_ga<T>::addRandomIndividuals()
{
  unsigned int iPopSize =  m_oNextPop.size();
  unsigned int iCount = 0;
  for (unsigned int iIndex = 0; iIndex < iPopSize; iIndex++) {
    if (svt_genrand() < m_fMutationProb) {
      T oInd = initIndividual();
      updateFitness(&oInd);

      m_oNextPop.push_back(oInd);
      m_oNextPop[ m_oNextPop.size() - 1].setOrigin(MUTATION);
      m_oNextPop[ m_oNextPop.size() - 1].resetAge();

      iCount++;
    }
  }
}

///////////////////////////////////////////////////////////////////////////////
// Transposition
///////////////////////////////////////////////////////////////////////////////

/**
 * transposition
 */
template<class T>
void svt_ga<T>::transposition()
{
  if (m_bVerbose)
    printf("svt_ga function: transposition() - in\n");

  Real64 fRand;
  for (unsigned int i = 0; i < m_oNextPop.size(); i++) {
    fRand = svt_genrand();

    if (fRand < m_fTranspositionProb) {
      T oInd;
      transpositionUniform(m_oNextPop[i], &oInd);

      oInd.setOrigin(TRANSPOSITION);
      oInd.resetAge();

      m_oNextPop.push_back(oInd);
    }
  }

  if (m_bVerbose)
    printf("svt_ga function: transposition() - out\n");
}

/**
 * flips two genes between the same individual
 * \param rParentA reference to object
 * \param pNewIndA pointer to new ind
 */
template<class T>
void svt_ga<T>::transpositionUniform(T &rParentA, T *pNewIndA)
{
  if (m_bVerbose)
    printf("svt_ga function: transpositionUniform() - in\n");

  Real64 fRand;
  int iPoint1, iPoint2;

  fRand = svt_genrand();
  iPoint1 = (int)(fRand * (Real64)(m_iGenes));

  do {
    fRand = svt_genrand();
    iPoint2 = (int)(fRand * (Real64)(m_iGenes));
  } while (iPoint1 == iPoint2);

  (*pNewIndA) = rParentA;

  pNewIndA->setGene(iPoint1, rParentA.getGene(iPoint2));
  pNewIndA->setGene(iPoint2, rParentA.getGene(iPoint1));

  if (m_bVerbose)
    printf("svt_ga function: transpositionUniform() - out\n");

};

///////////////////////////////////////////////////////////////////////////////
// Reinsertion
///////////////////////////////////////////////////////////////////////////////

/**
 * set Reinsertion scheme
 * \param iReinsertionScheme the reinsertion scheme to used
 */
template<class T>
void svt_ga<T>::setReinsertionScheme(unsigned int iReinsertionScheme)
{
  m_iReinsertionScheme = iReinsertionScheme;
};

/**
 * reinsertion
 */
template<class T>
void svt_ga<T>::reinsertion()
{
  if (m_bVerbose)
    printf("svt_ga function: reinsertion() - in\n");

  for (unsigned int i = 0; i < m_oNextPop.size(); i++)
    updateFitness(&m_oNextPop[i]);

  // which method should we use
  switch (m_iReinsertionScheme) {
    default:
    case REINSERTION_ELITIST:
      reinsertion_elitist();
      break;
    case REINSERTION_ELITIST_UNIQUE:
      reinsertion_elitist_unique();
      break;
    case REINSERTION_GLOBALRANKING:
      reinsertion_globalranking();
      break;
    case REINSERTION_GLOBALRANKING_UNIQUE:
      reinsertion_globalranking_unique();
      break;
    case REINSERTION_SHARING:
      reinsertion_sharing();
      break;
  }

  if (m_bVerbose)
    printf("svt_ga function: reinsertion() - out\n");

}

/**
 * elitist reinsertion - replace the 50% worst parents with 50% best new individuals
 */
template<class T>
void svt_ga<T>::reinsertion_elitist()
{
  sortPopulation();
  sortNextPopulation();

  // temp population
  svt_population<T> oTempPop;

  // copy the best of the old population...
  for (unsigned int i = 0; i < m_oPop.size(); i++)
    oTempPop.push_back(m_oPop[ m_oPop.size() - i - 1]);
  // and the best of the new population
  for (unsigned int i = 0; i < m_oNextPop.size(); i++)
    oTempPop.push_back(m_oNextPop[ m_oNextPop.size() - i - 1]);

  sort(oTempPop.begin(), oTempPop.end());

  // do we actually have enough?
  while ((int)oTempPop.size() < m_iPopSize) oTempPop.push_back(m_oPop[ m_oPop.size() - 1]);

  // copy everything over
  m_oPop.clear();
  for (int i = 0; i < m_iPopSize; i++)
    m_oPop.push_back(oTempPop[ oTempPop.size() - i - 1]);
}

/**
 * elitist reinsertion. Only the best individuals from both, the old and the new generations, are copied over to the next population. Some enhancements: Distance penalty and insertion of a few, fully random individuals.
 */
template<class T>
void svt_ga<T>::reinsertion_elitist_unique()
{
  if (m_bVerbose)
    printf("svt_ga function: reinsertion_elitist_unique() - in\n");

  sortPopulation();

  // create a couple of best individual mutations
  T oInd;
  T oIndAdd = m_oPop[ m_oPop.size() - 1];
  bool bAdd = false;
  Real64 fRand, fIntPart;
  for (unsigned int i = 0; i < (unsigned int)m_iGenes; i++) {
    for (int j = 0; j < m_iRefinementMaxMutPerGene; j++) {
      oInd = m_oPop[ m_oPop.size() - 1];

      fRand = svt_ranCauchy(0.0, m_fMutationOffset);
      fRand = modf(fRand, &fIntPart);

      oInd.setGene(i, oInd.getGene(i) + fRand);
      makeValid(&oInd);
      updateFitness(&oInd);
      if (oInd.getFitness() > oIndAdd.getFitness()) {
        oIndAdd = oInd;
        bAdd = true;
      }
    }
  }

  if (bAdd) {
    //oIndAdd.resetAge();
    oIndAdd.setOrigin(GENEREFINEMENT);
    m_oNextPop.push_back(oIndAdd);
  }
  sortNextPopulation();

  // temp population
  m_oTempPop.clear();
  // copy the best of the old population...
  m_oTempPop.insert(m_oTempPop.begin(), m_oPop.end() - m_iPopSize, m_oPop.end());
  // ...and the next population
  m_oTempPop.insert(m_oTempPop.begin(), m_oNextPop.end() - m_oNextPop.size(), m_oNextPop.end());
  // do we actually have enough?
  while ((int)m_oTempPop.size() < m_iPopSize) m_oTempPop.push_back(m_oPop[ m_oPop.size() - 1]);
  // now sort everything
  sort(m_oTempPop.begin(), m_oTempPop.end());

  //
  //Check which individuals are in tabu regions
  //
  for (unsigned int i = 0; i < m_oTempPop.size(); i++) {
    if (isInTabuReg(&m_oTempPop[i]))
      m_oTempPop[i].setFitness(0.0);
  }
  sort(m_oTempPop.begin(), m_oTempPop.end());
  penalizeSimilar(m_oTempPop, m_fCutoffDistance, m_fCutoffDistancePenalty);
  sort(m_oTempPop.begin(), m_oTempPop.end());

  //
  // add tabu region
  //
  m_oTabuWindow.push_back(m_oTempPop[m_oTempPop.size() - 1].distance(m_oPop[m_oPop.size() - 1]));

  if (m_oTabuWindow.size() > m_iTabuWindowSize)
    m_oTabuWindow.erase(m_oTabuWindow.begin());

  Real64 fAvg = 0.0;
  for (unsigned int i = 0; i < m_oTabuWindow.size(); i++)
    fAvg += m_oTabuWindow[i];
  fAvg = fAvg / m_oTabuWindow.size();

  int i = m_oTempPop.size() - 1;
  if ((m_oTabuWindow.size() >= m_iTabuWindowSize - 2 && fAvg < m_fTabuThreshold)) {
    bool bFound = false;
    for (unsigned int j = 0; j < m_oTabus.size(); j++)
      if (m_oTempPop[i].distance(m_oTabus[j]) < m_fTabuRegionSize)
        bFound = true;

    //add to tabu and refine
    if (!bFound)
      refineInd(&m_oTempPop[i]);
  }

  //
  // Sort the newly created population and copy it over into the main array...
  //
  sort(m_oTempPop.begin(), m_oTempPop.end());
  m_oPop.clear();
  m_oPop.insert(m_oPop.begin(), m_oTempPop.end() - m_iPopSize, m_oTempPop.end());

  if (m_bVerbose)
    printf("svt_ga function: reinsertion_elitist_unique() - out\n");
}

/**
 * share fitness among multiple nishes
 *\FIXME add more here
 */
template<class T>
void svt_ga<T>::reinsertion_sharing()
{
  vector< T > oTempPop;

  // copy the best of the old population...
  oTempPop.insert(oTempPop.begin(), m_oPop.end() - m_iPopSize, m_oPop.end());
  // ...and the next population
  oTempPop.insert(oTempPop.begin(), m_oNextPop.end() - m_iPopSize, m_oNextPop.end());

  // do we actually have enough?
  while ((int)oTempPop.size() < m_iPopSize) oTempPop.push_back(m_oPop[ m_oPop.size() - 1]);

  sort(oTempPop.begin(), oTempPop.end());

  shareFitness(oTempPop);

  sort(oTempPop.begin(), oTempPop.end());


  // add a certain amount of random individuals
  Real64 fOrigPop = 0.9;
  Real64 fRandPop = 0.1;
  m_oPop.clear();
  for (int i = 0; i < m_iPopSize * fOrigPop; i++)
    m_oPop.push_back(oTempPop[ oTempPop.size() - i - 1]);

  for (int i = 0; i < m_iPopSize * fRandPop; i++) {
    T oInd = initIndividual();
    updateFitness(&oInd);
    m_oPop.push_back(oInd);
  }

  sortPopulation();
}

/**
 * share the fitness between the individuals of the population: a %percent of individuals are allowed in one niche - the rest are just killed and need to populate other niches
 * \param oPop - the individuals come form these population
 */
template<class T>
void svt_ga<T>::shareFitness(vector< T> &oPop)
{
  vector<int> oNiche;
  vector<int> oRank;
  vector<int> oCountsPerNiche;
  int iPopSize = oPop.size();
  Real64 fDist;
  svt_matrix<Real64> oDistMat(iPopSize, iPopSize);

  for (int i = (int)oPop.size() - 1; i > 0; i--) {
    for (int j = i - 1; j > 0; j--) {
      fDist = oPop[i].distance(oPop[j]);

      oDistMat[i][j] = fDist;
      oDistMat[j][i] = fDist;

      if (fDist < EPS)
        oPop[j].setFitness(0.0);
    }
  }

  // set Niche to all individuals to 0
  for (int i = (int)iPopSize - 1; i > 0; i--) {
    oNiche.push_back(0);
    oRank.push_back(0);
  }

  // the best individual is in niche 1
  oNiche[ iPopSize - 1 ] = 1;
  oRank [ iPopSize - 1 ] = 0;

  //current niche
  int iNiche = 1;
  unsigned int iCount = 1;
  int i = (int)iPopSize - 1;

  // compute distances - get niches
  while (i >= 0) {
    for (int j = i - 1; j > 0; j--) {
      if (oNiche[j] == 0) { // j does not have a niche
        // get distance to the top of the niche
        fDist = oDistMat[i][j];

        if (fDist < m_fNicheSize && oPop[j].getFitness() > 0) {
          oNiche[ j ] = iNiche;
          oRank[ j ] = iCount; // the rank
          iCount++;
        }
      }
    }

    //done counting how many are in this niche - so add them
    oCountsPerNiche.push_back(iCount);

    i--;

    // don't make niche heads if individuals are already placed in a niche or its fitness is 0
    while (i >= 0 && (oNiche[i] != 0 || oPop[i].getFitness() == 0))
      i--;

    // next niche
    if (i >= 0) {
      iNiche++;
      oNiche[i] = iNiche;
      oRank[i] = 0;
      iCount = 1;
    }
  }

  for (i = 0; i < (int)oCountsPerNiche.size(); i++)
    oCountsPerNiche[ i ] = 0;

  //here share fitness
  for (i = iPopSize - 1; i >= 0; i--) {
    // if ind is not identical with other
    if (oPop[i].getFitness() > 0) {
      oCountsPerNiche[ oNiche[i] ]++ ;

      //decrease fitness inside niche and discard individuals if to many in the niche
      //SVTLBO << i <<  " " << oNiche[i] << " " << oCountsPerNiche[ oNiche[i] ] <<  " " <<  m_iPopSize * m_fMaxPopPerNiche << " " << oPop[i].getFitness() << " ";
      if (oCountsPerNiche[ oNiche[i] ] < m_iPopSize * m_fMaxPopPerNiche) {
        oPop[i].setFitness(oPop[i].getFitness()*pow(m_fSameNichePenalty, oRank[i]));
        oPop[i].setNiche(oNiche[i]);
      } else {
        oPop[i].setFitness(0);
        oPop[i].setNiche(0);
      }
      //cout << oPop[i].getFitness() << " " << oRank[i] << " " <<  m_fSameNichePenalty << endl;
    } else {
      oPop[i].setNiche(0);
    }
  }
}


/**
 * reinsertion - global reinsertion based on the global ranking.
 * Both, the old and the new population are ranked together and only the best individuals are inserted into the next gen.
 */
template<class T>
void svt_ga<T>::reinsertion_globalranking()
{
  for (unsigned int i = 0; i < m_oNextPop.size(); i++)
    m_oPop.push_back(m_oNextPop[i]);

  sortPopulation();

  // erase the worst individuals
  m_oNextPop.clear();
  for (int i = 0; i < m_iPopSize; i++) {
    m_oNextPop.push_back(m_oPop[ m_oPop.size() - i - 1]);
  }

  m_oPop = m_oNextPop;
}

/**
 * reinsertion - global reinsertion based on the global ranking but only unique individuals are allowed
 */
template<class T>
void svt_ga<T>::reinsertion_globalranking_unique()
{
  for (unsigned int i = 0; i < m_oNextPop.size(); i++)
    m_oPop.push_back(m_oNextPop[i]);

  sortPopulation();

  // remove identical individuals
  m_oNextPop.clear();
  svt_array_real64 oUsed;
  Real64 fItem;
  for (int i = m_oPop.size() - 1; (i >= 0 && (int)m_oNextPop.size() < m_iPopSize); i--) {
    fItem = m_oPop[i].getValue();

    if (oUsed.size() == 0 || find(oUsed.begin(), oUsed.end(), fItem) == oUsed.end()) {
      m_oNextPop.push_back(m_oPop[i]);
      oUsed.push_back(fItem);
    }
  }

  m_oPop = m_oNextPop;
}

///////////////////////////////////////////////////////////////////////////////
// Fitness
///////////////////////////////////////////////////////////////////////////////

/**
 * update fitness
 */
template<class T>
void svt_ga<T>::updateFitness()
{
  if (m_bVerbose)
    printf("svt_ga function: updateFitness() - out\n");

  for (unsigned int i = 0; i < m_oPop.size(); i++) {
    m_oPop[i].incAge();
    updateFitness(&m_oPop[i]);
    //cout << "[" << i << " - " << m_iGenerations << " - " << m_iThread <<  "] " ;
    //m_oPop[i].printGenes();
  }

  if (m_bVerbose)
    printf("svt_ga function: updateFitness() - out\n");

}

/**
 * check if the ind fullfills the validity requirements
 */
template<class T>
bool svt_ga<T>::isValid(T *pInd)
{
  for (int j = 0; j < m_iGenes; j++)
    if (pInd->getGene(j) > 1.0 || pInd->getGene(j) < 0.0)
      return false;

  return true;
};

/**
 * Function to check whether the integrity of the genes is maintained
 * \param pInd the individual for which to check genes
 * \return the if individual correct
 */
template<class T>
void svt_ga<T>::makeValid(T *pInd)
{
  //     if (isValid(pInd))
  //         return;

  svt_vector4<Real64> oVec;
  Real64 fNewGene;

  for (int j = 0; j < m_iGenes; j++) {
    fNewGene = pInd->getGene(j);

    if (fabs(fNewGene) > 1.0f)
      fNewGene -= floor(fNewGene);

    if (fNewGene < 0.0f)
      fNewGene = fabs(fNewGene);

    if (fabs(fNewGene) >= 1.0f)
      SVTLBO << "Error: Gene Value: " << j % 4 << ": " << fNewGene << " Origin: " << pInd->getOrigin() << endl;

    pInd->setGene(j, fNewGene);
  }
}

/**
 * Penalize individuals that are similar to allow a more diverse population
 * \param the population
 * \param fCutoffDistance the gene distance between which they get penalized
 * \param fCufoffDistancePenalty how much do they get penalized
 */
template<class T>
void svt_ga<T>::penalizeSimilar(svt_population<T> &oPop, Real64 fCutoffDistance, Real64 fCutoffDistancePenalty)
{
  Real64 fDist = 0.0;
  unsigned int iReduced = 0;
  unsigned int iZero = 0;
  int iSize = oPop.size();

  iReduced = 0;
  for (int i = iSize - 1; i > 0; i--) {
    for (int j = i - 1; j >= 0; j--) {
      fDist = oPop[i].distance(oPop[j]);

      // distance small?
      if (fDist < fCutoffDistance) {

        iReduced++;
        if (iReduced < iSize * 0.7)
          oPop[j].setFitness(oPop[j].getFitness() * fCutoffDistancePenalty);
        else
          oPop[j].setFitness(0.0);
      }

      // distance zero?
      if (fDist < EPS) {
        oPop[j].setFitness(0.0);
        iZero++;
      }
    }

  }
  //SVTLBO << "iReduced: "<< iReduced << " iZero:" << iZero << endl;
}

/**
 * Discard invalid(fitness value=0) individuals
 * \param oPop the population
 */
template<class T>
void svt_ga<T>::discardNullInd(svt_population<T> &oPop)
{
  //sort such the 0 are at the end of the list
  sort(oPop.rbegin(), oPop.rend());

  //discard 0
  int iIndex = oPop.size() - 1;
  while (iIndex >= 0 && oPop[iIndex].getFitness() == 0.0) {
    oPop.pop_back();
    iIndex = oPop.size() - 1;
  }
}

///////////////////////////////////////////////////////////////////////////////
// Print/output diagnostic information
///////////////////////////////////////////////////////////////////////////////

/**
 * print results (to cout)
 */
template<class T>
void svt_ga<T>::printResults()
{
  char pOut[1024];
  for (unsigned int i = 0; i < m_oPop.size(); i++) {
    sprintf(pOut, "[%2i] = ", i);
    SVTLBO << pOut ;
    m_oPop[i].printGenes();
  }
}

/**
 * print population
 * \param the population
 */
template<class T>
void svt_ga<T>::printPop(svt_population<T> &oPop)
{
  for (unsigned int i = 0; i < oPop.size(); i++) {
    printf("[%2i] = ", i);
    oPop[i].printGenesPf();
  }
}

/**
 * print results (to cout)
 */
template<class T>
void svt_ga<T>::printNextPop()
{
  char pOut[1024];
  unsigned int iSize = m_oNextPop.size();
  for (int i = iSize - 1; i >= 0; i--) {
    sprintf(pOut, "[%3i] = %1d %1d %8.3f", i, m_oNextPop[i].getOrigin(), m_oNextPop[i].getAge(), m_oNextPop[i].getProp());
    SVTLBO << pOut ;
    m_oNextPop[i].printGenes();
  }
}

/**
 * print the fitness of each individual of the population
 */
template<class T>
void svt_ga<T>::printPopFitness(char *pFname)
{

  FILE *file;
  if (m_iGenerations == 0)
    file = fopen(pFname, "w");
  else
    file = fopen(pFname, "a");

  for (unsigned int i = 0; i < m_oPop.size(); i++)
    fprintf(file, " %10.8f", -1.0 * (-1.0E10 + m_oPop[i].getFitness()));
  fprintf(file, "\n");

  fclose(file);
}

/**
 * Print the Min fitness, the avergage fitness and the Max fitness
 */
template<class T>
void svt_ga<T>::printStatistics()
{
  printf("%d\t%10.8f\t%10.8f\t%10.8f\n", m_iNoUniqueInd, -1.0 * (-1.0E10 + m_fMinFitness), -1.0 * (-1.0E10 + m_fAvgFitness), -1.0 * (-1.0E10 + m_fMaxFitness));
};

/**
 * output the configuration of the program
 */
template<class T>
void svt_ga<T>::writeConfiguration(char *pFnameParam)
{
  FILE *pFileParam = fopen(pFnameParam, "a");

  fprintf(pFileParam, "PopulationSize = %i\n",           getPopSize());
  fprintf(pFileParam, "ReinsertionScheme = %i\n",        m_iReinsertionScheme);
  fprintf(pFileParam, "MutationProbability = %f\n",      getMutationProb());
  fprintf(pFileParam, "MutationOffset = %f\n",           getMutationOffset());
  fprintf(pFileParam, "CrossoverProbability = %f\n",     getCrossoverProb());
  fprintf(pFileParam, "SelectivePressure = %f\n",        getSelectivePressure());
  fprintf(pFileParam, "MaxGenerations = %i\n",           getMaxGen());
  fprintf(pFileParam, "MaxThreads = %i\n",               getMaxThread());
  fprintf(pFileParam, "TranspositionProbability = %f\n", getTranspositionProb());
  fprintf(pFileParam, "DistanceThreshold = %f\n",        getDistanceThreshold());
  fprintf(pFileParam, "DistancePenalty = %f\n",          getDistanceThresholdPenalty());
  if (m_bMutateAll)
    fprintf(pFileParam, "MutateAll = true\n");
  else
    fprintf(pFileParam, "MutateAll = false\n");
  fprintf(pFileParam, "MutateAllProportion = %f\n",      getMutateAllProportion());
  fprintf(pFileParam, "StopScore = %f\n",                getStopScore());
  fprintf(pFileParam, "TabuWindowSize = %i\n",           getTabuWindowSize());
  fprintf(pFileParam, "TabuThreshold = %f\n",            getTabuThreshold());
  fprintf(pFileParam, "TabuRegionSize = %f\n",           getTabuRegionSize());
  fprintf(pFileParam, "NicheSize = %f\n",                getNicheSize());
  fprintf(pFileParam, "MapPopPerNiche = %f\n",          getMaxPopPerNiche());
  fprintf(pFileParam, "SameNichePenalty = %f\n",         getSameNichePenalty());
  fprintf(pFileParam, "RefinementMaxMutPerGene = %d\n",  getRefinementMaxMutPerGene());

  fclose(pFileParam);
};

///////////////////////////////////////////////////////////////////////////////
// run ga in thread
///////////////////////////////////////////////////////////////////////////////

/**
 * Thread function - it starts the ga in a thread
 */
template<class T>
void *runThread(void *pData)
{
  if (!pData)
    return NULL;

  svt_ga<T> *pGA = (svt_ga<T> *)pData;
  pGA->setIsThreadRunning(true);
  pGA->run();
  pGA->setIsThreadRunning(false);

  return NULL;
}


/**
 * function to create the thread
 */
template<class T>
void svt_ga<T>::initThread()
{
  svt_createThread(&runThread<T>, (void *)this, SVT_THREAD_PRIORITY_NORMAL);
}

///////////////////////////////////////////////////////////////////////////////
// tabu search functions
///////////////////////////////////////////////////////////////////////////////

/**
 * Set the tabu search window size. The tabu search computes the gene distances of the top-individual over time, with a moving window. It averages all those distance values.
 * If the distances vary a lot, because constantly completely new solutions get to the top, everything is considered fine. If the average drops, and only small differences
 * can be seen, this probably means premature convergence. The size of the window can be adjusted with this function.
 * \param iTabuWindowSize new size of the tabu-search window
 */
template<class T>
void svt_ga<T>::setTabuWindowSize(unsigned int iTabuWindowSize)
{
  m_iTabuWindowSize = iTabuWindowSize;
};
/**
 * Get the tabu search window size. The tabu search computes the gene distances of the top-individual over time, with a moving window. It averages all those distance values.
 * If the distances vary a lot, because constantly completely new solutions get to the top, everything is considered fine. If the average drops, and only small differences
 * can be seen, this probably means premature convergence. The size of the window can be accessed with this function.
 * \return size of the tabu-search window
 */
template<class T>
unsigned int svt_ga<T>::getTabuWindowSize()
{
  return m_iTabuWindowSize;
};

/**
 * At some point the distances of the top individuals get really small and we consider this as stagnation of the GA. With this function one can set the threshold, if the
 * average distance is lower, we store the top individual in an array and remove all individuals from this region.
 * \fTabuThreshold the new threshold below which we say the GA stagnates
 */
template<class T>
void svt_ga<T>::setTabuThreshold(Real64 fTabuThreshold)
{
  m_fTabuThreshold = fTabuThreshold;
};
/**
 * At some point the distances of the top individuals get really small and we consider this as stagnation of the GA. With this function one can access the threshold, if the
 * average distance is lower, we store the top individual in an array and remove all individuals from this region.
 * \return threshold below which we say the GA stagnates
 */
template<class T>
Real64 svt_ga<T>::getTabuThreshold()
{
  return m_fTabuThreshold;
};

/**
 * If the distance between an individual and a stored tabu region is smaller than this value, the individual is discarded.
 *\param fTabuRegionSize the new size of the tabu regions
 */
template<class T>
void svt_ga<T>::setTabuRegionSize(Real64 fTabuRegionSize)
{
  m_fTabuRegionSize = fTabuRegionSize;
};

/**
 * If the distance between an individual and a stored tabu region is smaller than this value, the individual is discarded.
 *\return the size of the tabu regions
 */
template<class T>
Real64 svt_ga<T>::getTabuRegionSize()
{
  return m_fTabuRegionSize;
};

/**
 * check whether the ind is in one of the tabu regions
 * \param pInd
 **/
template<class T>
bool svt_ga<T>::isInTabuReg(T *pInd)
{
  bool bInTabuReg = false;
  for (unsigned int j = 0; j < this->m_oTabus.size(); j++) {
    if (pInd->distance(this->m_oTabus[j]) < this->m_fTabuRegionSize)
      bInTabuReg = true;
  }
  return bInTabuReg;
}


/**
 * set the parent GA
 * \param pParentGa - the ga that started this thread ; NULL if the main thread
 */
template<class T>
void svt_ga<T>::setParentGA(svt_ga *pParentGA)
{
  m_pParentGA = pParentGA;
};

/**
 * get the parent GA
 * \return pParentGa - the ga that started this thread ; NULL if the main thread
 */
template<class T>
svt_ga<T> *svt_ga<T>::getParentGA()
{
  return m_pParentGA;
};


/**
 * refine an individual;
 * \param the individual that will be refined
 */
template<class T>
void svt_ga<T>::refineInd(T *pInd)
{
  char pOut[1024];
  sprintf(pOut, "[%02d-%04d] Added new tabu region (now %d): ", m_iThread, m_iGenerations, (int)m_oTabus.size());
  SVTLBO << pOut << endl;

  pInd->setAge(m_oTabus.size() + 1);

  m_oTabus.push_back(*pInd);
  sort(m_oTabus.begin(), m_oTabus.end());
  m_oTabuWindow.clear();

  //any of the individuals are in the new added tabu region
  for (unsigned int i = 0; i < m_oTempPop.size(); i++) {
    if (m_oTempPop[i].distance(m_oTabus[m_oTabus.size() - 1]) < m_fTabuRegionSize)
      m_oTempPop[i].setFitness(0.0);
  }
  //outputResult(true);

  //if added tabu just before merge , then let the ga run for a little longer to refine this solution
  if (m_iMaxGen - m_iGenerations < 15)
    m_iGenerations = m_iMaxGen - 15;

  return;
};

///////////////////////////////////////////////////////////////////////////////
//  Sharing
///////////////////////////////////////////////////////////////////////////////

/**
 * set the Niche size
 * \param the new nicheSize
 */
template<class T>
void svt_ga<T>::setNicheSize(Real64 fNicheSize)
{
  m_fNicheSize = fNicheSize;
};

/**
 * get the Niche size
 * \return nicheSize
 */
template<class T>
Real64 svt_ga<T>::getNicheSize()
{
  return m_fNicheSize;
};

/**
 * set the maximum allowed population per Niche - expressed as a proportion of the original population
 * Remarks: beyong this value (5-10%) - individuals have their fitness set to 0
 * \param proportion of individuals allowed in one niche
 */
template<class T>
void svt_ga<T>::setMaxPopPerNiche(Real64 fMaxPopPerNiche)
{
  m_fMaxPopPerNiche = fMaxPopPerNiche;
};

/**
 * get the maximum allowed population per Niche - expressed as a proportion of the original population
 * \param proportion of individuals allowed in one niche
 */
template<class T>
Real64 svt_ga<T>::getMaxPopPerNiche()
{
  return m_fMaxPopPerNiche;
};

/**
 * set the Niche distance penalty - penalize individuals in the same niche according to their rank to the top individual
 * \param how much will individuals be penalized
 */
template<class T>
void svt_ga<T>::setSameNichePenalty(Real64 fSameNichePenalty)
{
  m_fSameNichePenalty = fSameNichePenalty;
};

/**
 * get the Niche distance penalty - penalize individuals in the same niche according to their rank to the top individual
 * \param how much will the
 */
template<class T>
Real64 svt_ga<T>::getSameNichePenalty()
{

  return m_fSameNichePenalty;
};

/**
 * function to create the thread
 * it assumes that the object that calls it already has the pop initialized
 */
template<class T>
svt_population<T> svt_ga<T>::execute()
{
  static svt_semaphore oSema;

  if (!oSema.tryLock())
    return m_oPop;

  //clear arrays/lists if runs were already executed
  for (unsigned int iThread = 0; iThread < m_oGA_Array.size() && iThread < m_iMaxThread; iThread++) {
    if (m_oGA_Array[iThread] != NULL)
      delete(m_oGA_Array[iThread]);
  }
  m_oGA_Array.clear();

  delTabuRegions();
  setDone(false);

  vector< svt_ga * > oGA_Array;
  for (unsigned int iThread = 0; iThread < m_iMaxThread; iThread++) {
    svt_ga *pGA_Tmp = (*this).createObject() ;
    pGA_Tmp->setThread(iThread);
    pGA_Tmp->setRun(m_iRun);

    if (m_iMaxThread == 1) {
      pGA_Tmp->setMaxGen(m_iMaxGen);
      m_iSyncGen = m_iMaxGen;
    } else
      pGA_Tmp->setMaxGen(m_iSyncGen);

    pGA_Tmp->setParentGA(this);
    oGA_Array.push_back(pGA_Tmp);
  }
  m_oGA_Array = oGA_Array;

  long iStartTimeMerge, iStartTime = svt_getElapsedTime();
  unsigned int iThreadGenerations = 0;

#ifdef __GAFIT_FORK
  vector<pid_t> oPid;
  pid_t iPid;
#endif
  bool bDone = false;
  svt_population<T> oFinalPop, oTotalPop, oTotalTabuPop;

  for (m_iParallelRun = 0; (int)m_iParallelRun <  m_iMaxGen / m_iSyncGen && bDone == false; m_iParallelRun++) {
    iStartTimeMerge = svt_getElapsedTime();

    SVTLBO << "Synchronization " << (int)(m_iParallelRun) + 1 << " out of max. " << (int)floor(m_iMaxGen / m_iSyncGen) << ".";
    if (m_iMaxThread != 1) cout << " Starting " << m_iMaxThread << " threads..." << endl;
    else cout << " Starting 1 serial thread..." << endl;

    // start the threads
    for (unsigned int iThread = 0; iThread < m_iMaxThread; iThread++) {
      m_oGA_Array[iThread]->setCurrGen(0);
      m_oGA_Array[iThread]->setDone(false);
#ifndef __GAFIT_FORK
      m_oGA_Array[iThread]->initThread();
#else
      if ((iPid = fork()) == 0) {
        oPid.push_back(iPid);
        m_oGA_Array[iThread]->run();
        exit(0);
      }
#endif
    }

#ifdef __GAFIT_FORK
    while (oPid.size() < (unsigned int)m_iMaxThread);
#endif

    // make sure threads are really running
    bool bSomeNotRunning = true;
    while (bSomeNotRunning) {
      bSomeNotRunning = false;
      for (unsigned int iThread = 0; iThread < m_iMaxThread; iThread++)
        if (!m_oGA_Array[iThread]->getIsThreadRunning())
          bSomeNotRunning = true;
    }

    // wait until they are finished
    bool bAllFinished = false;
    while (!bAllFinished) {
      bAllFinished = true;
      for (unsigned int iThread = 0; iThread < m_iMaxThread && bDone == false; iThread++) {
        if (!m_oGA_Array[iThread]->getDone())
          bAllFinished = false;
        else if (m_oGA_Array[iThread]->getCurrGen() < m_iSyncGen) {
          //SVTLBO << "Thread " << iThread << " finished with stopping criterion (score higher than)!" << endl;
          bDone = true;
          for (unsigned int i = 0; i < m_iMaxThread; i++)
            if (i != iThread)
              m_oGA_Array[i]->setDone(true);
          iThreadGenerations += m_oGA_Array[iThread]->getCurrGen();
        }
      }
      if (!bAllFinished)
        svt_sleep(50);
    }


    // make sure threads are really not running anymore
    bool bSomeRunning = true;
    while (bSomeRunning) {
      bSomeRunning = false;
      for (unsigned int iThread = 0; iThread < m_iMaxThread; iThread++)
        if (m_oGA_Array[iThread]->getIsThreadRunning())
          bSomeRunning = true;
    }

    if (m_bDone) { // the user stop it most probably and by now all threads are done
      oSema.V();
      return m_oPop;
    }

    //delete the content of the populations and start again
    oTotalPop.clear();
    oTotalTabuPop.clear();

    // combine populations
    for (unsigned int iThread = 0; iThread < m_iMaxThread; iThread++) {
      svt_population<T> oTempPop = m_oGA_Array[iThread]->getPopulation();
      svt_population<T> oTabuTempPop = m_oGA_Array[iThread]->getTabuRegions();

      oTotalPop.insert(oTotalPop.begin(), oTempPop.begin(), oTempPop.end());
      oTotalTabuPop.insert(oTotalTabuPop.begin(), oTabuTempPop.begin(), oTabuTempPop.end());
    }

    // keep track of how many generations we have so far
    if (!bDone)
      iThreadGenerations += m_iSyncGen;

    sort(oTotalPop.begin(), oTotalPop.end());
    m_oGA_Array[0]->penalizeSimilar(oTotalPop, 2.0 * getDistanceThreshold(), (1 - 2.0 * (1 - getDistanceThresholdPenalty()))); //penalty twice larger than usual
    sort(oTotalPop.begin(), oTotalPop.end());

    sort(oTotalTabuPop.begin(), oTotalTabuPop.end());
    m_oGA_Array[0]->penalizeSimilar(oTotalTabuPop, 0.0, 1.0);
    m_oGA_Array[0]->discardNullInd(oTotalTabuPop);
    sort(oTotalTabuPop.begin(), oTotalTabuPop.end());

    m_oTabus = oTotalTabuPop;

    svt_sleep(1000);

    SVTLBO << "Recombination of the populations:" << endl;
    SVTLBO << "  Highest fitness: " <<  oTotalPop[oTotalPop.size() - 1].getFitness() << endl;
    SVTLBO << "  Number of generations (per thread): " <<  iThreadGenerations << endl;
    SVTLBO << "  Number of tabu regions  " << (int)oTotalTabuPop.size() << endl;
    SVTLBO << "  Partial Elapsed time : " << (((svt_getElapsedTime() - iStartTimeMerge) / (1000.0f)) / 60.0f) << " min" << endl;

    if (!bDone) {
      svt_population<T> oNewPop;
      oNewPop.insert(oNewPop.begin(), oTotalPop.end() - m_iPopSize, oTotalPop.end());
      for (unsigned int iThread = 0; iThread < m_iMaxThread; iThread++) {
        m_oGA_Array[iThread]->setPopulation(oNewPop);
        m_oGA_Array[iThread]->setTabuRegions(oTotalTabuPop);
        if (iThread == 0)
          m_oGA_Array[iThread]->outputResult(true);
      }
    }
  }

  //Output runtime stats
  long iETime = svt_getElapsedTime() - iStartTime;
  SVTLBO << "Number of generations (per thread): " << iThreadGenerations <<  endl;
  SVTLBO << "Elapsed time: " << ((iETime / (1000.0f)) / 60.0f) << " min" <<  endl;

  //clear the content of the threads
  clearThreads();

  //Now combine the final population with the tabu regions
  oFinalPop = oTotalPop;
  if (oTotalTabuPop.size() > 0)
    oFinalPop.insert(oFinalPop.begin(), oTotalTabuPop.begin(), oTotalTabuPop.end());

  sort(oFinalPop.rbegin(), oFinalPop.rend());

  char pFilename[1024];
  strcpy(pFilename, "GlobalSearchSolution");
  writeSolutions(oFinalPop, 1, pFilename);

  oSema.V();

  return oFinalPop;
}

/**
 * Clear the content of the threads
 */
template<class T>
void svt_ga<T>::clearThreads()
{
  unsigned int iFitnessUpdateCount = 0;
  for (unsigned int iThread = 0; iThread < m_iMaxThread; iThread++) {
    iFitnessUpdateCount += m_oGA_Array[iThread]->getFitnessUpdateCount();

    //clear for the next run
    //m_oGA_Array[iThread]->delTabuRegions();
    //delete(m_oGA_Array[iThread]);
  }
  //delete clear all
  //m_oGA_Array.clear();

  SVTLBO << "Fitness computed  " << iFitnessUpdateCount << " times" << endl;
  SVTLBO << endl;
}

/**
 * Refine population
 * \param oPop what population
 * \param iNoInd4Refinement how many individuals are to be refined
 */
template<class T>
void svt_ga<T>::refine(svt_population<T> &oPop, unsigned int iNoInd4Refinement)
{
  if (oPop.size() > 0) {
    delTabuRegions();

    unsigned int iEffectivePopSize = iNoInd4Refinement < oPop.size() ? iNoInd4Refinement : oPop.size();
    SVTLBO << endl;
    SVTLBO << "Restart GA only with " << iEffectivePopSize << " tabus" << endl;
    SVTLBO << endl;

    long iStartTime = svt_getElapsedTime();

    //store settings
    unsigned int iPopSize    = m_iPopSize;
    unsigned int iMaxGen       = m_iMaxGen;
    unsigned int iMaxThread    = m_iMaxThread;
    Real64 fTabuThreshold      = m_fTabuThreshold;



    //set population here to assess the new scores of this populaiton
    setPopulation(oPop);
    updateFitness();
    oPop = getPopulation();
    sort(oPop.rbegin(), oPop.rend());
    oPop.erase(oPop.begin() + iEffectivePopSize, oPop.end());

    setPopulation(oPop);

    setMaxGen(101);
    setCurrGen(0);
    setThread(0);
    setMaxThread(0);
    setTabuThreshold(0.0);

    run();

    long iETime = svt_getElapsedTime() - iStartTime;
    SVTLBO << "Elapsed time: " << ((iETime / (1000.0f)) / 60.0f) << " min" <<  endl;

    //set settings used before refinement
    setPopSize(iPopSize);
    setMaxGen(iMaxGen);
    setMaxThread(iMaxThread);
    setTabuThreshold(fTabuThreshold);

    SVTLBO << "Number of tabu regions at the end of refine " << m_oTabus.size() << endl;

  }
}



/**
 * Thread function - it starts the ga in a thread
 */
template<class T>
void *runThreads(void *pData)
{
  if (!pData)
    return NULL;

  svt_ga<T> *pGA = (svt_ga<T> *)pData;

  pGA->setIsThreadRunning(true);
  pGA->execute();
  pGA->setIsThreadRunning(false);

  return NULL;
}


/**
 * function to create the thread
 */
template<class T>
void svt_ga<T>::initThreads()
{
  svt_createThread(&runThreads<T>, (void *)this, SVT_THREAD_PRIORITY_NORMAL);
}


///////////////////////////////////////////////////////////////////////////////
// svt_eulerAngles
///////////////////////////////////////////////////////////////////////////////

/**
 * This class generates a nondegenerate set of Euler angles
 * in the specified scan range. For the full sphere, the number
 * of points is almost identical to the spiral method but without
 * the helical slant and the weirdness at the poles. Angle
 * generation for subintervals is far superior and even cases
 * where the range can't be evenly devided by delta are handled
 * gracefully.
 * The method works by dividing the sphere into slices and
 * determining how many psi angles need to be in a slice to re-
 * produce the same spacing as on the equator.
*/
class svt_eulerAngles
{
  protected:

    float *m_pAngles;
    unsigned long m_iAngles;

  public:

    /**
     * Constructor
     */
    svt_eulerAngles();

    /**
     * Destructor
     */
    virtual ~svt_eulerAngles();

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
    void initTable(double fPsiFrom, double fPsiTo, double fThetaFrom, double fThetaTo, double fPhiFrom, double fPhiTo, double fDelta);

    /**
     * remove angles from table
     * \param fPsiFrom   lower boundary of the psi angles (in degrees)
     * \param fPsiTo     upper boundary of the psi angles (in degrees)
     * \param fThetaFrom lower boundary of the theta angles (in degrees)
     * \param fThetaTo   upper boundary of the theta angles (in degrees)
     * \param fPhiFrom   lower boundary of the phi angles (in degrees)
     * \param fPhiTo     upper boundary of the phi angles (in degrees)
     * WARNING: not fully tested
     */
    void removeAngles(double fPsiFrom, double fPsiTo, double fThetaFrom, double fThetaTo, double fPhiFrom, double fPhiTo);

    /**
     * add the oppsite angles: if angle = (psi, theta, phi)
     * add (-psi, theta, phi) (psi, theta, -phi) (-psi, theta, -phi)
     * eq with (2pi-psi, theta, phi) (psi, theta, 2pi-phi) (2pi-psi, theta, 2pi-phi)
     */
    void addOppositeAngles();


    /**
     * How many angles do we have stored in total?
     * \return unsigned long with the number of angles
     */
    unsigned long getAngleCount();

    /**
     * The angles follow the common PTP (Goldstein) convention. This function returns the Psi angle.
     * \param iIndex index into the table of angles
     * \return psi angle
     */
    float getPsi(unsigned long iIndex);
    /**
     * The angles follow the common PTP (Goldstein) convention. This function returns the Theta angle.
     * \param iIndex index into the table of angles
     * \return theta angle
     */
    float getTheta(unsigned long iIndex);
    /**
     * The angles follow the common PTP (Goldstein) convention. This function returns the Phi angle.
     * \param iIndex index into the table of angles
     * \return phi angle
     */
    float getPhi(unsigned long iIndex);

    /**
     * searches the angles that are within fAngleRange from the angle indicated by iIndex
     * \param iIndex the reference angle around which to search
     * \param fRange how far away from the reference
     * \return a list of indexes that indicates the angles in the angle list that are close to the angle idicated by iIndex
     */
    vector<long unsigned int> getNeighborAngles(unsigned long int iIndex, Real64 fAngleRange);

  protected:

    /**
     * This function precomputes the angle table. It is called automatically in the constructor.
     */
    unsigned long proportionalEulerAngles(double fPsiFrom, double fPsiTo, double fThetaFrom, double fThetaTo, double fPhiFrom, double fPhiTo, double fDelta);

    /**
     * This function precomputes the angle table. It is called automatically in the constructor.
     */
    unsigned long proportionalEulerAnglesThetaOrdering(double fPsiFrom, double fPsiTo, double fThetaFrom, double fThetaTo, double fPhiFrom, double fPhiTo, double fDelta);
};

///////////////////////////////////////////////////////////////////////////////
// SVT_GACYLINDER_IND
///////////////////////////////////////////////////////////////////////////////

/**
 * GA multifit individual
 * \author Mirabela Rusu
 */
class svt_gacylinder_ind : public svt_ga_ind
{
  protected:
    // represents the encoded solution in a very simplified reprentation
    // here we use 4 atoms for each unit (no genes/4) that represent the
    // Cartesian coordinate systems ( (0,0,0), (1,0,0), (0,1,0), (0,0,1) )
    // to which the rotation and translations were applied
    vector< svt_vector4<Real64> > m_oCoarsePhenotype;

    //wrote on disk
    bool m_bWrote;

    // the number of turns in the cylinder - the cylinder has 2*m_iTurns+1 circles inside the cylinder
    unsigned int m_iTurns;

    // the height of a turn
    Real64 m_fHeightTurn;

    // the transformation matrix that was applied to original vector
    svt_ga_mat m_oTrans;

    //fitness of the top part
    Real64 m_fFitnessTop;

    //fitness of the bottom part
    Real64 m_fFitnessBot;


  public:

    /**
     * Constructor
     */
    svt_gacylinder_ind();

    /**
     * destructor
     */
    ~svt_gacylinder_ind();

    /**
     * create the coarse phenotype (equivalent simple pdb) but don't fill yet atomic coordinates
     * \param number of points
     */
    void buildCoarsePhenotype();

    /**
     * update the coarse phenotype for the given unit
     * \param transformation matrix for that unit
     * \param number of units
     */
    void updateCoarsePhenotype(svt_ga_mat oMat);


    /**
     * get the coarse phenotype
     */
    svt_point_cloud_pdb<svt_ga_vec> getCoarsePhenotype();

    /**
     * calculate the distance between two individuals
     * \param rOther reference to the other individual
     * \return vector distance between the two gene-vectors
     */
    virtual Real64 distance(svt_gacylinder_ind &rOther);

    /**
     * set Wrote on disk
     * \param bWrote whether it was already wrote on disk
     */
    void setWrote(bool bWrote);

    /**
     * get Wrote on disk
     * \return bWrote whether it was already wrote on disk
     */
    bool getWrote();

    /**
     * get the number of turns
     * \return the number of turns
     */
    unsigned int getTurns();

    /**
     * get the number of turns
     * \param the number of turns
     */
    void setTurns(unsigned int iTurns);

    /**
     * get the height of a turn
     * \return the height
     */
    Real64 getHeightTurn();

    /**
     * set the height of a turn
     * \param the height
     */
    void setHeightTurn(Real64 fHeightTurn);

    /**
     * get the Transformation
     * \return the transformation
     */
    svt_ga_mat getTrans();

    /**
     * set the Transformation
     * \param the transformation
     */
    void setTrans(svt_ga_mat &rTrans);


    /**
     * set Fitness Top
     * \param the new fitness
     */
    void setFitnessTop(Real64 fFitness);

    /**
     * get Fitness Top
     * \return the new fitness
     */
    Real64 getFitnessTop();

    /**
     * set Fitness Bot
     * \param the new fitness
     */
    void setFitnessBot(Real64 fFitness);

    /**
     * get Fitness Top
     * \return the new fitness
     */
    Real64 getFitnessBot();

};

///////////////////////////////////////////////////////////////////////////////
// SVT_TUBE
///////////////////////////////////////////////////////////////////////////////

/***************************************************************************
                          svt_tube
                          a tube
                          ---------
    begin                : 04/26/2010
    author               : Mirabela Rusu
    email                : Mirabela.Rusu@uth.tmc.edu
 ***************************************************************************/

/**
 *  a cylinder
 * \author Mirabela Rusu
 */
class svt_tube
{
  protected:
    // the gacylinder ind that create the cylinder
    vector<svt_gacylinder_ind > m_oElements;

    //another score
    Real64 m_fAvgFitness;

    //another score
    Real64 m_fSumFitness;

    //the length of the tube - it represents the number of step/turns that were added to the tube
    unsigned int m_iSize;

    //the length of the tube in A
    Real64 m_fLength; //

    // the score based on the map
    Real64 m_fMapScore;

    //default NO: once the tube is computed then true
    bool m_bWasTubeComputed;

    //default NO: once the tube is computed then true
    bool m_bWasFineTubeComputed;

    // the points on the axis
    svt_point_cloud_pdb<svt_ga_vec> m_oTube;

    // the points on the axis
    svt_point_cloud_pdb<svt_ga_vec> m_oFineTube;

    // the hr tube
    svt_point_cloud_pdb<svt_ga_vec> m_oHRTube;

    //
    svt_vector4<Real64> m_oFirstAddedElem;

    //penatly - score that is substracted from the original score to penalize for characteristics (e.g. axes crossing)
    Real64 m_fPenalty;

    // curvature of the axis
    Real64 m_fCurvature;

    //should the score be recomputed
    bool m_bRecomputeScore;

    //scores
    vector<Real64> m_oScores;

    //
    svt_ga_vec m_oDirection;
  public:

    /**
     * Constructor
     */
    svt_tube();

    /**
     * destructor
     */
    virtual ~svt_tube() {};

    /**
     * add a new element to tube
     * \param oElem a svt_gacylinder_ind that in considered in this cylinder
     */
    void add2Tube(svt_gacylinder_ind oElem, bool bFlip = false);

    /**
     * remove the last element of the tube
     */
    void pop();

    /**
     * get the list of individuals added to the tube
     * \return the individuals
     */
    vector<svt_gacylinder_ind > getElements();

    /**
     * set the map score
     * \param fScore the score assigned
     */
    void setMapScore(Real64 fScore);

    /**
     * get the map score
     * \return the score of the tube
     */
    Real64 getMapScore();

    /**
     * set the score
     * \param fScore the score assigned
     */
    void setAvgScore(Real64 fScore);

    /**
     * get the score
     * \return the score of the tube
     */
    Real64 getAvgScore();

    /**
     * Compute avg score
     */
    void computeAvgScore();

    /**
    * set the score
    * \param iIndex - which element in vector score
    * \param fScore - the value
    */
    void setScore(int iIndex, Real64 fScore);

    /**
    * returns the list of scores
    * \param a vector with the different scores
    */
    vector<Real64> getScores() ;

    /**
       * get the number of points
       * \return the number of points
       */
    unsigned int size();

    /**
       * get the number of turns
       * \return the number of turns
       */
    Real64 getTurns(Real64 fRatio = 1.0);

    /**
     * Compute the length in A0
     */
    void computeLength();

    /**
     * get the length = sum distances between points
     * assumes that the points are in order = point 0 is closest 1 and point 2 follows point 1... etc
     * \return the length of the tube
     */
    Real64 getLength();

    /**
     * get the first element
     */
    svt_vector4<Real64> getFirstAddedElem();

    /**
     * set penalty
     * \param the penalty
     */
    void setPenalty(Real64 fPenalty);

    /**
     * increase penalty
     * \param the penalty
     */
    void addPenalty(Real64 fPenalty);

    /**
     * get penalty
     * \return the penalty
     */
    Real64 getPenalty();

    /**
     * get curvature
     */
    Real64 getCurvature();

    /**
     * overload < operator
     * \param that another svt_tube element
     */
    bool operator<(const svt_tube &that) const;

    /**
     * overload < operator using the max density of the map
     * \param that another svt_tube element
     */
    static bool lt_mapScore(svt_tube first, svt_tube second);

    /**
     * overload < operator using the max density of the map
     * \param that another svt_tube element
     */
    static bool lt_length(svt_tube first, svt_tube second);


    /**
     * overload < operator using the max density of the map
     * \param that another svt_tube element
     */
    static bool lt_wPenatly(svt_tube first, svt_tube second);

    /**
     * overload < operator using the max density of the map
     * \param that another svt_tube element
     */
    static bool lt_score(svt_tube first, svt_tube second);

    /**
     * get tube as
     * \param a template, e.g. circle
     * \param should return the fine version
     * \return the pdb tube
     */
    svt_point_cloud_pdb<svt_ga_vec> getTube(svt_point_cloud_pdb<svt_ga_vec> *pTemplate = NULL, bool bFine = false);

    /**
     * get the direction of the tube
     */
    svt_ga_vec getDirection();

    /**
     * print the tube
     */
    void print();

    /**
     * compute the volume underneath the tube
     */
    void fillExplored(svt_ga_vol &rVol, svt_ga_vol *pVol);

    /**
     *get the high resolution version of the tube
     */
    svt_point_cloud_pdb<svt_ga_vec> getHRTube();

    /**
     * creates a highresolution version of the helix
     */
    void createHighResTube(svt_ga_vol &rVol, svt_point_cloud_pdb<svt_ga_vec> *pTemplate);

    /**
     * Estimate curvature of the axis
     */
    void estimate_curvature();

    /**
     * discard the points that are at the ends of the tube is their score is < fScore
     */
    void discardPointsAtEnds(Real64 fScore);
};


///////////////////////////////////////////////////////////////////////////////
// SVT_GACYLINDER
///////////////////////////////////////////////////////////////////////////////

/***************************************************************************
                          svt_gacylinder
                          fit multiple domains into map
                          ---------
    begin                : 01/21/2009
    author               : Mirabela Rusu
    email                : Mirabela.Rusu@uth.tmc.edu
 ***************************************************************************/


/**
 * GA-based fitting algorithm
 * \author Mirabela Rusu
 */
template<class T> class svt_gacylinder : public svt_ga<T>
{
  protected:
    //the template of the cylinder
    svt_point_cloud_pdb<svt_ga_vec> m_oTemplate, m_oBigTemplate, m_oTurn, m_oBigTurn, m_oCircle1, m_oCircle2, m_oCircle3;

    //target volume
    svt_ga_vol m_oTar;

    //the target structure for validation purposes only
    svt_point_cloud_pdb<svt_ga_vec> m_oTarStr;

    // the axis of the target structure
    svt_point_cloud_pdb<svt_ga_vec> m_oTarStrAxes;

    //the target structure for validation purposes only
    svt_point_cloud_pdb<svt_ga_vec> m_oTarStrCaInHelix;

    // a pdb that contains only the coa of the atoms in each axis
    svt_point_cloud_pdb<svt_ga_vec> m_oTarStrAxesCenters;

    // holds information about the helices; i.e. length
    vector< vector< Real64 > > m_oTarStrAxesInfo;

    // a volume representing the blurred axes
    svt_ga_vol m_oTarStrAxesVol;

    //the target structure for validation purposes only
    svt_point_cloud_pdb<svt_ga_vec> m_oCoarseTarget;

    //the model
    svt_point_cloud_pdb<svt_ga_vec> m_oModel;

    //model volume
    svt_ga_vol m_oModelVol;

    //the coordinate of the 0,0,0 voxel
    Real64 m_fOrigX, m_fOrigY, m_fOrigZ;

    //no of voxels on x,y, and z
    unsigned int m_iSizeX, m_iSizeY, m_iSizeZ;

    // voxel width
    Real64 m_fWidth;

    //the resolution
    Real64 m_fRes;

    // the size of the search spacetransClear()
    svt_ga_vec m_oSearchSpace;

    // the scoring function to be used
    //score m_eScore;

    //the kernel for the blurring of the model
    svt_ga_vol m_oKernel;

    // set of non-degenerate euler angles
    static svt_eulerAngles m_oAngles;

    // the angular step size
    Real64 m_fDelta;

    // the angular search range
    Real64 m_fPsiFrom;
    Real64 m_fPsiTo;
    Real64 m_fThetaFrom;
    Real64 m_fThetaTo;
    Real64 m_fPhiFrom;
    Real64 m_fPhiTo;


    bool m_bSetTranslSearchRange;
    // the translational search range
    Real64 m_fXFrom;
    Real64 m_fXTo;
    Real64 m_fYFrom;
    Real64 m_fYTo;
    Real64 m_fZFrom;
    Real64 m_fZTo;


    // result output
    char m_pPath[1256];

    // how many generations should gacylinder wait until it outputs a new model file (0 turns output off)
    unsigned int m_iWriteModelInterval;

    //SpringPotential
    Real64 m_fSpringPotential;

    //the max distance deviation between two points: dist(p2-p1) in state1 - dist(p2-p1) in state2
    Real64 m_fMaxDistDev;

    //how much should be disregarded on each side
    Real64 m_fBorder;

    //mask durring the correlation
    bool m_bMask;

    //value between 0-1; default 0.95 that indicates that a crawl is accepted if the score is larger that m_fAcceptMovePercent * original score
    Real64 m_fAcceptMoveRatio;

    //Number of "failed" (score < m_fAcceptMovePercent * original score) crawl steps allowed before stopping the crawl
    unsigned int m_iMaxFailedCrawls;

    //a cylinder
    vector <svt_tube> m_oCylinders;

    //a cylinder
    vector <  svt_point_cloud_pdb<svt_ga_vec> > m_oCylindersPdb;

    //pdb of the axis with the anisotromic decompression applied
    vector <  svt_point_cloud_pdb<svt_ga_vec> > m_oCylindersPdbAniDecorr;


    //any cylinders were already identified
    bool m_bCanOutputResults;

    //is refining now
    bool m_bRefining;

    //created durring the crawling; added only when the crawling was completely done
    svt_population<T> m_oTabusInCrawl;

    //was a "valid" cylinder found in this run?
    bool m_bFoundCylinder;

    // number of helices to detect
    unsigned int m_iNoOfCylinder2Detect;

    // apply blurring to the model
    bool m_bApplyBlurring2Model;

    // the radius of the expansion template in A
    Real64 m_fTemplateRadius;

    // the radius of the ga search template in A
    Real64 m_fSearchTemplateRadius;

    // the number of points in the template ( m_iTemplatePointCount - 1 points in the circle and 1 at the center  )
    unsigned int m_iTemplatePointCount;

    // number of times the circle is repeated in the expansion template
    int m_iTemplateRepeats;

    // number of times the circle is repeated in the search template
    int m_iSearchTemplateRepeats;

    // distance between repeats
    Real64 m_fDistBetweenRepeats;

    // distance covered in one crawling step
    Real64 m_fCrawlingStepSize;

    // the number of steps per turn
    Real64 m_fStepsPerTurn;

    //explored volume
    svt_ga_vol m_oExploredVol;

    //the amount ot the explored volume
    Real64 m_fExploredPercent;

    //an idealized hardcodded ala helix - created using modeller
    svt_point_cloud_pdb<svt_ga_vec> m_oIdealHelix;

    //should a helix be fitted on the axis
    bool m_bFitHelix;

    //
    svt_point_cloud_pdb<svt_ga_vec> m_oHelices;

    //the maximum distanced to consider 2 atoms / tubes as corresponding/mapped to each other
    Real64 m_fMaxDist4Map;

    //the rmsd output to print in outputResults; is computed in discardCylindres
    char m_pRmsdOut[1024];

    // the number of individuals created in the initInd; helps figuring if the space is full with tabu regions and the search should be stopped
    int m_iNoIndInInit;

    //the average score of all solutions
    Real64 m_fAvgScoreAll;

    //should the templates be updated
    bool m_bOutputTemplates;

    // Anisotric correction related terms:
    Real64 m_fAni;
    Real64 m_fOrigZWoAni;

    // the correction Lamdba - see update fitness
    Real64 m_fLambda;
  public:

    /**
     * Constructor
     */
    svt_gacylinder(unsigned int iGenes);

    /**
     * Destructor
     */
    virtual ~svt_gacylinder();

    /**
     * create a new individual
     */
    virtual T initIndividual();

    /**
     * generate initial population
     * \param iNum number of inds in this population
     */
    virtual void initPopulation(int iNum, bool bOutputInput = false);

    /**
     * initialize the angles
     */
    void initAngles();

    /**
     * Create an object
     */
    virtual svt_gacylinder<T> *createObject();


    /**
     * create the volume template
     */
    void initTemplate();

    /**
     * add helices found in threads to the global list of helices; discard the one that overlap
     * \param bFromThreads should only the helices from threads should be looked at? or all helices
     * \return the average length of the added helices
     */
    Real64 discardCylinders(bool bFromThreads = true);

    /**
     * discard points from the cylinders according to their score // use here the average score of all accepted points in crawling
     */
    void discardPointsoOnCylinders();

    /**
     * delete all helices
     */
    void  clearCylinders();

    /**
     * post process tubes
     */
    void postProcessTubes();
    /**
     * sort the cylinders and fills the Cylinder pdb
     */
    void sortCylinders();

    /**
     * function to create the thread
     * \return a population
     */
    virtual svt_population<T> execute();


    //
    // Parameters
    //

    /**
     * Set resolution of the target map
     * \param fRes the resolution of the target map
     */
    void setResolution(Real64 fRes);

    /**
     * Get resolution of the target map
     * \return the resolution of the target map
     */
    Real64 getResolution();

    /**
     * Set the masking when computing the correlation
     * \param bMask bool indicating whether to mask
     */
    void setMask(bool bMask);

    /**
     * Set the angular step size. In general: The smaller the value, the better the accuracy. In contrast to an exhaustive search, the runtime will also not be longer, if a finer step size
     * is chosen. The only limitation is the memory - as the table needs to be stored, very low numbers might result in an excessive use of memory. Recommended: ~0.5 to 1.0.
     * \param fDelta the angular step size (default 0.5)
     */
    void setAngularStepSize(Real64 fDelta);
    /**
     * Get the angular step size. In general: The smaller the value, the better the accuracy. In contrast to an exhaustive search, the runtime will also not be longer, if a finer step size
     * is chosen. The only limitation is the memory - as the table needs to be stored, very low numbers might result in an excessive use of memory. Recommended: ~0.5 to 1.0.
     * \return the angular step size (default 0.5)
     */
    Real64 getAngularStepSize();

    /**
     * Returns the number of angles
     */
    Real64 getAnglesCount();

    /**
     * Get angles
     */
    svt_vector4<Real64> &getAngle(long unsigned int iIndex);

    /**
     * Set the ranges for the angular search.
     * \param fPsiFrom   lower limit of the psi angles
     * \param fPsiTo     upper limit of the psi angles
     * \param fThetaFrom lower limit of the theta angles
     * \param fThetaTo   upper limit of the theta angles
     * \param fPhiFrom   lower limit of the phi angles
     * \param fPhiTo     upper limit of the phi angles
     */
    void setAngularSearchRange(Real64 fPsiFrom, Real64 fPsiTo, Real64 fThetaFrom, Real64 fThetaTo, Real64 fPhiFrom, Real64 fPhiTo);

    /**
     * Set the ranges for the translational search relative to the centers of the units
     * (if fXFrom = -20 and fXTo = 20 the units moves -20 and 20 A from the current position)
     * \param fXFrom     upper limit of the X
     * \param fXTo       lower limit of the X
     * \param fYFrom     upper limit of the Y
     * \param fYTo       lower limit of the Y
     * \param fZFrom     upper limit of the Z
     * \param fZTo       lower limit of the Z
     */
    void setRelativeTranslSearchRange(Real64 fXFrom, Real64 fXTo, Real64 fYFrom, Real64 fYTo, Real64 fZFrom, Real64 fZTo);

    //
    // Units, Maps, PDB files...
    //
    /**
     * set target
     * \param oTar the map
     */
    void setTarget(svt_ga_vol &oTar);

    /**
     * set target str
     * only for validation purposes
     * \param oTarStr the target structure
     */
    void setTarStr(svt_point_cloud_pdb<svt_ga_vec>  &oTarStr, bool bComputeAxis = true);

    /**
     * set the tempate of the cylinder
     * \param oPdb the target structure
     */
    void setTemplate(svt_point_cloud_pdb<svt_ga_vec>  &oPdb);

    /**
     * set coarse target
     * \param oCoarseTarget the coarse version of the target structure
     */
    void setCoarseTarget(svt_point_cloud_pdb<svt_ga_vec>  &oCoarseTarget);

    /**
     * Get model
     * \return the model that was last generated by calcTransformation
     */
    inline svt_point_cloud_pdb<svt_ga_vec>  &getModel()
    {
      return m_oModel;
    }

    /**
     * Get model volume
     * \return the volume of the model that was last generated by calcTransformation
     */
    inline svt_ga_vol  &getModelVol()
    {
      return m_oModelVol;
    }

    /**
     * set the accept Move score percentage:e.g. 0.95 means that moves are accepted if withing 0.95 of the original score
     * \param the score percentage
     */
    void setAcceptMoveRatio(Real64 fAcceptMoveRatio);

    /**
     * get the accept Move score percentage: e.g. 0.95 means that moves are accepted if withing 0.95 of the original score
     * \return the score percentage
     */
    Real64 getAcceptMoveRatio();

    /**
     * set the max number of failed Crawls before stoping the search
     * \param iMaxFailedCrawls e.g. 2 indicates 2 times tried before stoping the crawl
     */
    void setMaxFailedCrawls(unsigned int iMaxFailedCrawls);

    /**
     * get the max number of failed Crawls before stoping the search
     * \return iMaxFailedCrawls e.g. 2 indicates 2 times tried before stoping the crawl
     */
    unsigned int getMaxFailedCrawls();

    /**
     * how many Cylinders shouls be searched for
     */
    void setNoOfCylinder2Detect(int iNoOfCylinder2Detect);

    /**
     * how many Cylinders shouls be searched for
     */
    int getNoOfCylinder2Detect();

    /**
     * was "valid" cylinder found
     */
    void setFoundCylinder(bool bFoundCylinder);

    /**
     *  was "valid" cylinder found
     */
    bool wasFoundCylinder();

    /**
     *  set apply Blurring to model
     * \param state
     */
    void setApplyBlurring2Model(bool bBlurring);

    /**
     * set the radius of the template
     * \param the radius of the template
     */
    void setTemplateRadius(Real64 fTemplateRadius);

    /**
     * get the radius of the template
     * \param the radius of the template
     */
    Real64 getTemplateRadius();

    /**
     * set the radius of the template used for the ga Search
     * \param the radius of the template
     */
    void setSearchTemplateRadius(Real64 fSearchTemplateRadius);

    /**
     * get the radius of the template used for the ga Search
     * \param the radius of the template
     */
    Real64 getSearchTemplateRadius();

    /**
     * set the number of points in one of the circles of the template
     * \param  number of points in one of the circles of the template
     */
    void setTemplatePointCount(unsigned int iTemplatePointCount);

    /**
     * get the number of points in one of the circles of the template
     * \param  number of points in one of the circles of the template
     */
    unsigned int getTemplatePointCount();

    /**
     * set the number of circles copied in the expansion template
     * \param the number of repeats to be set
     */
    void setTemplateRepeats(unsigned int iTemplateRepeats);

    /**
     * get the number of circles copied in the expansion template
     * \return the number of repeats
     */
    unsigned int getTemplateRepeats();

    /**
     * set the number of circles copied in the search template
     * \param the number of repeats to be set
     */
    void setSearchTemplateRepeats(unsigned int iSearchTemplateRepeats);

    /**
     * get the number of circles copied in the search template
     * \return the number of repeats
     */
    unsigned int getSearchTemplateRepeats();

    /**
     * set the distance between two repeats
     * \param the distance between repeats
     */
    void setDistBetweenRepeats(Real64 fDistBetweenRepeats);

    /**
     * get the distance between two repeats
     * \return the distance between repeats
     */
    Real64 getDistBetweenRepeats();

    /**
     * set the size of a crawling step
     * \param the size of a crawling step
     */
    void setCrawlingStepSize(Real64 fCrawlingStepSize);

    /**
    * get the size of a crawling step
    * \return the size of a crawling step
    */
    Real64 getCrawlingStepSize();

    /**
     * set fit high resolution helix
     */
    void setFitHelix(bool bFitHelix);
    /**
     * get fit high resolution helix
     */
    bool getFitHelix();

    /**
     * are results available and can be outputed
     */
    bool canOutputResults();

    /**
     * set the outputTemplate option
     * \param the bOutputTemplates
     */
    void setOutputTemplates(bool bOutputTemplates);

    /**
     * get the outputTemplate option
     * \return OutputTemplates
     */
    bool getOutputTemplates();

    /**
     * set anisotropic correction related parameters
     * \param fAni the correction, e.g. fAni=2 will compress the map 2 times
     * \param fOrigZWoAni the origin of the map before ani correction - use for decompression
     */
    void setAniCorr(Real64 fAni, Real64 fOrigZWoAni);

    /**
     * get anisotropic correction
     * \param fAni the correction, e.g. fAni=2 will compress the map 2 times
     */
    Real64 getAniCorr();

    /**
     * get Origin on Z axis of the map without anisotropic correction
     */
    Real64 getOrigZWoAni();

    /**
     * get Origin on Z axis of the map with anisotropic correction
     */
    Real64 getOrigZ();

    /**
     * Set the correction factor lambda
     */
    void setLambda(Real64 fLambda);

    /**
     * Get the correction factor lambda
     */
    Real64 getLambda();


    /**
     * init kernels
     */
    void initKernels();


    //
    // Fitness function
    //

    /**
     * returns the transformation matrix
     */
    svt_ga_mat getTrans(T *pInd);

    /**
     * Updates the coordinates of the model based on the genes of the individual
     * attention: it does not update the volume of the model (see update volume)
     */
    void updateModel(T *pInd);

    /**
     * Updates the volume of the model
     * Attention: it does not update the model - the pdb remains the same
     * \param pInd the individual
     * \param bCoarse should the coarse model be used
     */
    void updateVolume(T *pInd);

    /**
     * update fitness
     * \param pInd pointer to individual that gets updated
     */
    virtual void updateFitness(T *pInd);

    /**
     * refine an individual;
     * \param the individual that will be refined
     */
    virtual void refineInd(T *pInd);

    /**
     * refine an individual; by default it does nothing - overload in the classes that inherit
     * \param the individual that will be refined
     * \param iteration
     */
    void refine(T *pInd, int iIter);

    /**
     * refine the translation and rotation of an individual - calls the refinetransl and refineRot a few times
     * \param the individual that will be refined
     */
    void refineGenes(T *pInd, svt_vector4<Real64> *pCenter = NULL);

    /**
     * refine the translation of an individual
     * \param the individual that will be refined
     */
    void refineTransl(T *pInd);
    /**
     * refine the translation of an individual, only allow movements on Spheres
     * \param the individual that will be refined
     */
    void refineTranslOnSphere(T *pInd, svt_vector4<Real64> *pCenter);

    /**
     * refine the translation of an individual
     * \param the individual that will be refined
     */
    void refineRot(T *pInd, Real64 fMaxCorr = 0.9925);

    /**
     * refine the rotation of an individual using random angles
     * \param the individual that will be refined
     */
    void refineRotRandom(T *pInd);

    /**
     * crawl on the tube and search for similar scoring cylinder placements
     * \param the individual that will be refined
     * \param iDirection indicates the direction  +1 for forward and -1 for backwards
     */
    void crawl(T *pInd, int iDirection);

    /**
     * refine the length of an individual
     * \param the individual that will be refined
     */
    void refineInLength(T *pInd);

    /**
     * add 2 Cylinder
     */
    void add2Tube(T *pInd, bool bAddNew = false, bool bFlip = false);

    /**
     * set the cylinderes
     * \param the vector containing the helices
     */
    void setCylinders(vector<svt_tube> &oCylinders);
    /**
     * get the cylinderes
     * \return the vector containing the helices
     */
    vector<svt_tube> &getCylinders();
    /**
     * set the cylinderes
     * \param the vector containing the helices
     */
    void setCylindersPdb(vector< svt_point_cloud_pdb<svt_ga_vec> > &oCylinders);
    /**
     * get the cylinderes
     * \return the vector containing the helices
     */
    vector< svt_point_cloud_pdb<svt_ga_vec> > &getCylindersPdb();
    /**
     * Calculate the full correlation corresponding to the current individual (with blur)
     */
    Real64 getCorrelation();

    /**
     * Calculate the rmsd of the individual towards the target model
     * \return the rmsd
     */
    Real64 getRMSD();

    //
    // Mutation
    //

    /**
     * custom mutation (can be changed by derived class, default implementation just calls mutationBGA)
     * \param iInd index of individual
     */
    virtual void mutationCustom(int iInd);

    /**
     * mutation with a cauchy distribution
     * \param iInd index of individual
     */
    void mutationCauchy(int iInd, int iRandIndex, Real64 fRatio = 1.0);

    /**
     * mutation all the genes with cauchy of large standard deviation
     */
    void mutationAllCauchy(int iInd);

    /**
     * Compute Average Score of all solutions
     */
    void computeAvgScoreAll();

    //
    // Transposition
    //

    //
    // Output statistics, result files, etc
    //

    /**
     * output result
     */
    virtual void outputResult(bool bTabuAdded = false);

    /**
     * update result
     * \param the number of tubes to investigate
     * \param an int used for validation only
     */
    virtual void updateResults(unsigned int iNoOfTubes = 0, int iNum = 0);

    /**
      * output the best model
      */
    void outputBest();

    /**
     * print results (to cout)
     */
    void printResults();

    /**
     * writes the genes and scores of all individuals into file
     */
    void writePop(char *pFname_target);

    /**
     * Set the output path path
     * \param pPath pointer to array of char
     */
    void setOutputPath(const char *pPath);

    /**
     * Get the output path path
     * \param pPath pointer to array of char
     */
    const char *getOutputPath();

    /**
     * How many generations should gacylinder wait until it outputs a new model file (0 turns output off)
     * \param iWriteModelInterval number of generations
     */
    void setWriteModelInterval(unsigned int iWriteModelInterval);
    /**
     * How many generations should gacylinder wait until it outputs a new model file (0 turns output off)
     * \return number of generations
     */
    unsigned int getWriteModelInterval();

    /**
     * output the configuration of the program
     */
    virtual void writeConfiguration(char *pFilename);


    //
    // Run in thread
    //

    /**
     * Write the top scoring solutions to the disk
     * \param oPop the population of solutions
     * \param iWriteSolutions how many solutions to write
     */
    void writeSolutions(svt_population<T> &oPop, unsigned int iWriteSolutions, char *pFileName);

};

///////////////////////////////////////////////////////////////////////////////
// Constructor
///////////////////////////////////////////////////////////////////////////////
template<class T>
svt_eulerAngles svt_gacylinder<T>::m_oAngles;

/**
 * Constructor
 */
template<class T>
svt_gacylinder<T>::svt_gacylinder(unsigned int iGenes) : svt_ga<T>(iGenes),
  m_fWidth(0.0f),
  m_fRes(8.0f),
  m_fDelta(0.5f),
  m_fPsiFrom(0.0f),
  m_fPsiTo(360.0f),
  m_fThetaFrom(0.0f),
  m_fThetaTo(180.0f),
  m_fPhiFrom(0.0f),
  m_fPhiTo(360.0f),
  m_fXFrom(0.0f),
  m_fXTo(0.0f),
  m_fYFrom(0.0f),
  m_fYTo(0.0f),
  m_fZFrom(0.0f),
  m_fZTo(0.0f),
  m_fBorder(5.0f),
  m_bMask(true),
  m_fAcceptMoveRatio(0.90),
  m_iMaxFailedCrawls(2),
  m_bCanOutputResults(false),
  m_bRefining(false),
  m_iNoOfCylinder2Detect(0),
  m_bApplyBlurring2Model(false),
  m_fTemplateRadius(1.0),
  m_fSearchTemplateRadius(2.0),
  m_iTemplatePointCount(11),
  m_iTemplateRepeats(8),
  m_iSearchTemplateRepeats(20),
  m_fDistBetweenRepeats(5.1 / 4.0),
  m_fCrawlingStepSize(1.4/*5.4/3.6*/),    //experimental value
  m_fExploredPercent(0.0),
  m_bFitHelix(false),
  m_fMaxDist4Map(4.0),
  m_fAni(1.0),
  m_fLambda(1.0)
{
  strcpy(m_pPath, "");
  m_fStepsPerTurn = 5.4 / m_fCrawlingStepSize;
  sprintf(m_pRmsdOut, " ");
};


/**
 * Destructor
 */
template<class T>
svt_gacylinder<T>::~svt_gacylinder() {};

/**
 * create a new individual
 */
template<class T>
T svt_gacylinder<T>::initIndividual()
{
  T oInd;
  svt_array_real64 oGenes;
  // create new object - that overlap a bit with the map
  do {
    unsigned int iX, iY, iZ;
    do {
      oGenes.clear();
      for (int j = 0; j < 3; j++)
        oGenes.push_back(svt_genrand());
      iX = int (oGenes[0] * m_oTar.getSizeX());
      iY = int (oGenes[1] * m_oTar.getSizeY());
      iZ = int (oGenes[2] * m_oTar.getSizeZ());
    } while (m_oTar.getValue(iX, iY, iZ) <= 0.0);

    for (int j = 3; j < this->m_iGenes; j++)
      oGenes.push_back(svt_genrand());
    oInd.setGenes(oGenes);
    oInd.setTurns(m_iTemplateRepeats);
    oInd.setHeightTurn(m_fDistBetweenRepeats);

    this->makeValid(&oInd);
    updateFitness(&oInd);
    m_iNoIndInInit++;
  } while (oInd.getFitness() < 1e-6 || this->isInTabuReg(&oInd));

  oInd.setOrigin(RANDOM);
  oInd.resetAge();

  return oInd;
}

/**
 * initialize the angles
 */
template<class T>
void svt_gacylinder<T>::initAngles()
{
  // create the angle table
  if (m_oAngles.getAngleCount() <= 0) {
    SVTLBO << "Create angular search table (P/T/P, Delta, #Angles, memory): ( " << m_fPsiFrom << " - " << m_fPsiTo << " / " << m_fThetaFrom << " - " << m_fThetaTo << " / " << m_fPhiFrom << " - " << m_fPhiTo << ", " << m_fDelta << ", ";

    m_oAngles.initTable(m_fPsiFrom, m_fPsiTo, m_fThetaFrom, m_fThetaTo, m_fPhiFrom, m_fPhiTo, m_fDelta);

    cout << m_oAngles.getAngleCount() << ", " << floor((m_oAngles.getAngleCount() * sizeof(float) * 3) / 1048576) << "mb )" << endl;
  }

}

/**
 * generate initial population
 * \param iNum number of inds in this population
 */
template<class T>
void svt_gacylinder<T>::initPopulation(int iNum, bool bOutputInput)
{
  m_iNoIndInInit = 0;

  svt_ga<T>::initPopulation(iNum);

  if (strlen(m_pPath) != 0 && bOutputInput) {
    // output the files
    char pFname[1256];
    if (m_fWidth != 0.0) {
      sprintf(pFname, "%s/target.mrc", m_pPath);
      //m_oTar.save( pFname );
    }

    if (m_oTarStr.size() > 0) {
      sprintf(pFname, "%s/target.pdb", m_pPath);
      m_oTarStr.writePDB(pFname);
    }

  }
}

/**
 * Create an object
 */
template<class T>
svt_gacylinder<T> *svt_gacylinder<T>::createObject()
{
  return new svt_gacylinder<T>(*this);
};

/**
 * create the volume template
 */
template<class T>
void svt_gacylinder<T>::initTemplate()
{
  svt_ga_vec oVec;
  if (m_oTurn.size() == 0) { // create template
    m_oTurn.delAllPoints();
    m_oBigTurn.delAllPoints();

    m_oCircle1.delAllPoints();
    m_oCircle2.delAllPoints();
    m_oCircle3.delAllPoints();

    oVec.x(0.0);
    oVec.y(0.0);
    oVec.z(0.0);

    svt_point_cloud_atom oAtom;
    oAtom.setName("C");
    oAtom.setRemoteness('A');
    oAtom.setResidueSeq(1);

    m_oTurn.addAtom(oAtom, oVec);
    m_oBigTurn.addAtom(oAtom, oVec);

    Real64 fDelta = 360.0 / Real64(m_iTemplatePointCount - 1);
    for (Real64 fAngle = 0; fAngle < 360; fAngle += fDelta) {
      //create the circle used by the expansion template
      oVec.x(m_fTemplateRadius * cos(deg2rad(fAngle)));
      oVec.y(m_fTemplateRadius * sin(deg2rad(fAngle)));

      oAtom.setMass(1.0);
      oAtom.setResidueSeq(int(fAngle / fDelta) + 2);
      m_oTurn.addAtom(oAtom, oVec);

      //create the circle used by the search template
      oVec.x(m_fSearchTemplateRadius * cos(deg2rad(fAngle)));
      oVec.y(m_fSearchTemplateRadius * sin(deg2rad(fAngle)));

      oAtom.setMass(1.0);
      m_oBigTurn.addAtom(oAtom, oVec);

      oVec.x(2 * m_fTemplateRadius * cos(deg2rad(fAngle)));
      oVec.y(2 * m_fTemplateRadius * sin(deg2rad(fAngle)));
      m_oCircle1.addAtom(oAtom, oVec);

      oVec.x(4 * m_fTemplateRadius * cos(deg2rad(fAngle)));
      oVec.y(4 * m_fTemplateRadius * sin(deg2rad(fAngle)));
      m_oCircle2.addAtom(oAtom, oVec);

      oVec.x(5 * m_fTemplateRadius * cos(deg2rad(fAngle)));
      oVec.y(5 * m_fTemplateRadius * sin(deg2rad(fAngle)));
      m_oCircle3.addAtom(oAtom, oVec);
    }
  }

  /*
   if (m_oBigTurn.size() == 0) // create template
   {
       m_oBigTurn.delAllPoints();

       oVec.x( 0.0 );
       oVec.y( 0.0 );
       oVec.z( 0.0 );

       svt_point_cloud_atom oAtom;
       oAtom.setName( "C");
       oAtom.setRemoteness('A');

       oAtom.setMass(1.0);
       m_oBigTurn.addAtom( oAtom, oVec);

       Real64 fDelta = 360.0/Real64(m_iTemplatePointCount-1);
       for (Real64 fAngle = 0; fAngle < 360; fAngle += fDelta)
       {
           oVec.x ( m_fSearchTemplateRadius * cos (deg2rad( fAngle ) ) );
           oVec.y ( m_fSearchTemplateRadius * sin (deg2rad( fAngle ) ) );

           oAtom.setMass(1.0);
           m_oBigTurn.addAtom( oAtom, oVec );
       }
   }
  */
  if (m_oTemplate.size() == 0) { //no template yet
    m_oTemplate = m_oTurn;

    oVec.x(0.0);
    oVec.y(0.0);

    svt_ga_mat oMat;
    svt_point_cloud_pdb<svt_ga_vec> oPdb;
    for (int iTurn = 0; iTurn < m_iTemplateRepeats; iTurn++) {
      oVec.z(m_fDistBetweenRepeats * iTurn);
      oMat.setTranslation(oVec);

      oPdb = oMat * m_oTurn;
      m_oTemplate.append(oPdb);
    }

    //center
    svt_ga_vec oCoa = m_oTemplate.coa();

    oMat.loadIdentity();
    oMat.setTranslation(-oCoa);
    m_oTemplate = oMat * m_oTemplate;

    //m_oModel = m_oTemplate;
  }

  if (m_oBigTemplate.size() == 0) { //no template yet
    m_oBigTemplate = m_oBigTurn;

    oVec.x(0.0);
    oVec.y(0.0);

    svt_ga_mat oMat;
    svt_point_cloud_pdb<svt_ga_vec> oPdb;
    //for helices used 2.5 ; for filopodia - 2
    for (int iTurn = 0; iTurn < m_iSearchTemplateRepeats; iTurn++) {
      oVec.z(m_fDistBetweenRepeats * iTurn);
      oMat.setTranslation(oVec);

      oPdb = oMat * m_oBigTurn;
      m_oBigTemplate.append(oPdb);
    }

    //center
    svt_ga_vec oCoa = m_oBigTemplate.coa();

    oMat.loadIdentity();
    oMat.setTranslation(-oCoa);
    m_oBigTemplate = oMat * m_oBigTemplate;
  }

  if (strlen(m_pPath)) {
    // output the files
    char pFname[1256];
    if (m_oTemplate.size() > 0 && m_bOutputTemplates) {
      sprintf(pFname, "%s/ExpansionTemplate.pdb", m_pPath);
      m_oTemplate.writePDB(pFname);
    }

    if (m_oBigTemplate.size() > 0 && m_bOutputTemplates) {
      sprintf(pFname, "%s/GATemplate.pdb", m_pPath);
      m_oBigTemplate.writePDB(pFname);
    }

  }

};

/**
 * add helices found in threads to the global list of helices; discard the one that overlap
 * \param bFromThreads should only the helices from threads should be looked at? or all helices
 * \return the average length of the added helices
 */
template<class T>
Real64 svt_gacylinder<T>::discardCylinders(bool bFromThreads)
{
  //
  // get the cylinders to be considered : all vs the one detected at this thread
  //
  vector<svt_tube> oCylinders, oTmpList;
  if (bFromThreads) {
    for (unsigned int iThread = 0; iThread < this->m_oGA_Array.size(); iThread++) {
      oCylinders = ((svt_gacylinder<T> *)this->m_oGA_Array[ iThread ])->getCylinders();
      oTmpList.insert(oTmpList.begin(), oCylinders.begin(), oCylinders.end());

      ((svt_gacylinder<T> *)this->m_oGA_Array[ iThread ])->clearCylinders();
    }
  } else
    oTmpList = m_oCylinders;


  //
  //check how they intersect - how much
  //
  svt_point_cloud_pdb<svt_ga_vec> oPdb1, oPdb2;
  svt_ga_vol oVol1, oVol2;

  oVol2 = m_oTar;
  oVol1 = m_oTar;

  //assume all helices should be added
  vector<bool> oAdd;
  for (unsigned int iIndex = 0; iIndex < oTmpList.size(); iIndex++)
    oAdd.push_back(true);

  for (unsigned int iIndex = 0; iIndex < oTmpList.size(); iIndex++) {
    oPdb1 = oTmpList[iIndex].getTube();
    oAdd.push_back(true);
    //oTmpList[iIndex].estimate_curvature();

    oVol1.setValue(0.0);
    oPdb1.projectMass(&oVol1);
    oVol1.convolve1D3D(m_oKernel, false);  // don't normalize

    for (unsigned int iIndex1 = iIndex + 1; iIndex1 < oTmpList.size(); iIndex1++) {
      oPdb2 = oTmpList[iIndex1].getTube();

      oVol2.setValue(0.0);
      oPdb2.projectMass(&oVol2);
      oVol2.convolve1D3D(m_oKernel, false);
      Real64 fCorr = oVol2.correlation(oVol1, false);

      //they overlap alot?
      if (fCorr > 0.5) {
        if (oTmpList[iIndex].getAvgScore()*oTmpList[iIndex].getAvgScore()*oTmpList[iIndex].getLength() <
            oTmpList[iIndex1].getAvgScore()*oTmpList[iIndex1].getAvgScore()*oTmpList[iIndex1].getLength())
          oAdd[ iIndex ] = false;
        else
          oAdd[ iIndex1 ] = false;
      }
    }
  }


  //delete all before putting back the one that don't overlap;
  if (!bFromThreads)
    m_oCylinders.clear();

  // for the one that are to be added; compute the other scores
  // also create high-resolution helix
  Real64 fAvg = 0;
  unsigned int iCount = 0;
  for (unsigned int iIndex = 0; iIndex < oTmpList.size(); iIndex++) {
    if (oAdd[iIndex]) {
      fAvg += oTmpList[iIndex].getTurns(m_fStepsPerTurn);
      iCount++;

      if (m_bFitHelix)
        oTmpList[iIndex].createHighResTube(m_oTar, &m_oIdealHelix);

      m_oCylinders.push_back(oTmpList[iIndex]);
      oTmpList[iIndex].fillExplored(m_oTar, &m_oExploredVol);

      oPdb1 = m_oCylinders[m_oCylinders.size() - 1].getTube();
      oPdb2 = m_oCylinders[m_oCylinders.size() - 1].getTube(&m_oCircle3);

      svt_ga_mat oMat;
      Real64 fScore;
      fScore = oPdb1.projectMassCorr(&m_oTar, oMat, false);
      m_oCylinders[m_oCylinders.size() - 1].setScore(0, fScore);
      fScore = oPdb2.projectMassCorr(&m_oTar, oMat, false);
      m_oCylinders[m_oCylinders.size() - 1].setScore(1, fScore);

      Real64 fAvgCcExpansion = m_oCylinders[m_oCylinders.size() - 1].getAvgScore();
      Real64 fLen = m_oCylinders[m_oCylinders.size() - 1].getLength();
      Real64 fCCInterior = m_oCylinders[m_oCylinders.size() - 1].getScores()[0];
      Real64 fCCExterior = m_oCylinders[m_oCylinders.size() - 1].getScores()[1];
      m_oCylinders[m_oCylinders.size() - 1].setScore(2,
          fAvgCcExpansion * fAvgCcExpansion * fCCInterior * fCCInterior / (fCCExterior * fCCExterior)*fLen);
      if (bFromThreads)
        m_oCylinders[m_oCylinders.size() - 1].setScore(3, this->m_iParallelRun);
    }
  }

  if (iCount > 0)
    fAvg /= (Real64)iCount;

  //
  //at the end and when a target is given compute some statistics - RMSD
  //
  if (!bFromThreads && m_oTarStr.size() > 0) {
    svt_point_cloud_pdb<svt_ga_vec> oPdb, oPdbHel;
    vector <int> oAxesMapped, oCAMapped;
    vector <int> oPredMapped; // the predicted points that were already
    vector <Real64> oAxesPointsMinDist, oCAMinDist, oPredMinDis;
    vector <Real64> oPredPointsMinDist; // the min distance from the predicted points to the actual axes


    for (unsigned int iA2 = 0 ; iA2 < m_oTarStrAxes.size(); iA2++) {
      oAxesMapped.push_back(0);
      oAxesPointsMinDist.push_back(1e10);
      m_oTarStrAxes.getAtom(iA2)->setTempFact(0.0);
    }

    for (unsigned int iA2 = 0 ; iA2 < m_oTarStrCaInHelix.size(); iA2++) {
      oCAMapped.push_back(0);
      oCAMinDist.push_back(1e10);
      m_oTarStrCaInHelix.getAtom(iA2)->setTempFact(0.0);
    }

    svt_point_cloud_pdb<svt_ga_vec> oPredPoints;
    for (unsigned int iIndex = 0; iIndex < m_oCylinders.size() && iIndex < m_iNoOfCylinder2Detect; iIndex++) {
      //rmsd on the axis
      oPdb = m_oCylinders[iIndex].getTube();
      for (unsigned int iA1 = 0 ; iA1 < oPdb.size(); iA1++) {
        oPredPoints.addAtom(*oPdb.getAtom(iA1), oPdb[iA1]);
        oPredPoints.getAtom(oPredPoints.size() - 1)->setTempFact(0.0);
        oPredPointsMinDist.push_back(1e10);
      }
    }

    Real64 fDist, fDistMin, fRmsd = 0.0;

    Real64 fMaxDist4MapSq = m_fMaxDist4Map * m_fMaxDist4Map;
    int iCountPredPointsMatched = 0;

    for (unsigned int iA1 = 0 ; iA1 < oPredPoints.size(); iA1++) {
      fDistMin =  oPredPoints[iA1].distanceSq(m_oTarStrAxes[0]);
      for (unsigned int iA2 = 1 ; iA2 < m_oTarStrAxes.size(); iA2++) {
        fDist = oPredPoints[iA1].distanceSq(m_oTarStrAxes[iA2]);
        if (fDist < fDistMin)
          fDistMin = fDist;
      }

      oPredPointsMinDist[iA1] = fDistMin;
      if (fDistMin < fMaxDist4MapSq) {
        iCountPredPointsMatched++;
        oPredPoints.getAtom(iA1)->setTempFact(7.0);
      } else
        oPredPoints.getAtom(iA1)->setTempFact(1.0);

    }

    int iCountAxesPointsMatched = 0;
    for (unsigned int iA2 = 0 ; iA2 < m_oTarStrAxes.size(); iA2++) {
      fDistMin = oPredPoints[0].distanceSq(m_oTarStrAxes[iA2]);
      for (unsigned int iA1 = 1 ; iA1 < oPredPoints.size(); iA1++) {
        fDist = oPredPoints[iA1].distanceSq(m_oTarStrAxes[iA2]);
        if (fDist < fDistMin)
          fDistMin = fDist;
      }

      oAxesPointsMinDist[iA2] = fDistMin;
      if (fDistMin < fMaxDist4MapSq) {
        iCountAxesPointsMatched++;
        m_oTarStrAxes.getAtom(iA2)->setTempFact(7.0);
      } else
        m_oTarStrAxes.getAtom(iA2)->setTempFact(1.0);
    }

    iCount = 0;
    fRmsd = 0;
    for (unsigned int iA2 = 0 ; iA2 < m_oTarStrAxes.size(); iA2++)
      if (oAxesPointsMinDist[iA2] < fMaxDist4MapSq) {
        fRmsd += oAxesPointsMinDist[iA2];
        iCount++;
      }

    if (iCount != 0)
      fRmsd /= iCount;

    fRmsd = sqrt(fRmsd);

    Real64 fSe  = iCountAxesPointsMatched / (Real64)m_oTarStrAxes.size();
    Real64 fPpv = iCountPredPointsMatched / (Real64)oPredPoints.size();

    sprintf(m_pRmsdOut, " RMSD: %6.3f %6.2f %6.2f %4d %4d %4d %4d %4d", fRmsd, fSe * 100.0, fPpv * 100.0, iCountAxesPointsMatched, m_oTarStrAxes.size() , iCountPredPointsMatched, oPredPoints.size(), m_oTarStrAxes.size());
  }

  return fAvg;
}

/**
 * discard points from the cylinders according to their score // use here the average score of all accepted points in crawling
 */
template<class T>
void  svt_gacylinder<T>::discardPointsoOnCylinders()
{
  char pOut[1024];

  // compute the score
  computeAvgScoreAll();

  SVTLBO << "Average score " << m_fAvgScoreAll << endl;

  for (unsigned int iIndex = 0; iIndex < m_oCylinders.size(); iIndex++) {
    sprintf(pOut, " %2d %3d %8.3f ", iIndex, (int)m_oCylinders[iIndex].size(),  m_oCylinders[iIndex].getAvgScore());
    m_oCylinders[iIndex].discardPointsAtEnds(m_fAvgScoreAll * 0.9);
    sprintf(pOut, "%s %2d %3d %8.3f ", pOut, iIndex, (int)m_oCylinders[iIndex].size(),  m_oCylinders[iIndex].getAvgScore());
    SVTLBO << pOut << endl;
  }

  sortCylinders();

  while (m_oCylinders.size() > 0 && m_oCylinders[m_oCylinders.size() - 1].getAvgScore() == 0)
    m_oCylinders.pop_back();


};

/**
 * delete all helices
 */
template<class T>
void svt_gacylinder<T>::clearCylinders()
{
  m_oCylinders.clear();
};

/**
 * post process tubes
 */
template<class T>
void svt_gacylinder<T>::postProcessTubes()
{
  svt_ga_vol oVol1, oVol2;
  svt_point_cloud_pdb<svt_ga_vec> oPdb1, oPdb2;
  oVol2 = m_oTar;
  oVol1 = m_oTar;

  //assume all helices should be added
  vector<bool> oAdd;
  for (unsigned int iIndex = 0; iIndex < m_oCylinders.size(); iIndex++) {
    oAdd.push_back(true);
    m_oCylinders[iIndex].setPenalty(0.0);
  }

  bool bFoundMultipleTimes;
  for (unsigned int iTube1 = 0; iTube1 < m_oCylinders.size(); iTube1++) {
    oPdb1 = m_oCylinders[iTube1].getTube();

    oVol1.setValue(0.0);
    oPdb1.projectMass(&oVol1);
    oVol1.convolve1D3D(m_oKernel, false);  // don't normalize

    bFoundMultipleTimes = false;

    if (oAdd[iTube1]) // not yet decided to discard
      for (unsigned int iTube2 = 0; iTube2 < m_oCylinders.size() ; iTube2++) {
        if (iTube2 != iTube1) {
          oPdb2 = m_oCylinders[iTube2].getTube();

          oVol2.setValue(0.0);
          oPdb2.projectMass(&oVol2);
          oVol2.convolve1D3D(m_oKernel, false);
          Real64 fCorr = oVol2.correlation(oVol1, false);

          //they overlap?
          if (fCorr > 0.5) {
            bFoundMultipleTimes = true;
            if (m_oCylinders[iTube1].getAvgScore() < m_oCylinders[iTube2].getAvgScore())
              oAdd[ iTube1 ] = false;
            else
              oAdd[ iTube2 ] = false;
          }
        }
      }

    if (!bFoundMultipleTimes)
      oAdd[ iTube1 ] = false;
  }

  vector<svt_tube> oCylinders;
  for (unsigned int iIndex = 0; iIndex < m_oCylinders.size(); iIndex++)
    if (oAdd[iIndex])
      oCylinders.push_back(m_oCylinders[iIndex]);

  clearCylinders();
  m_oCylinders = oCylinders;

  sortCylinders();
};

/**
 * sort the cylinders and fills the Cylinder pdb
 */
template<class T>
void svt_gacylinder<T>::sortCylinders()
{
  sort(m_oCylinders.rbegin(), m_oCylinders.rend(), svt_tube::lt_score);

  m_oCylindersPdb.clear();

  for (unsigned int iIndex = 0; iIndex < m_oCylinders.size(); iIndex++)
    m_oCylindersPdb.push_back(m_oCylinders[iIndex].getTube(NULL, true));
};

/**
 * function to create the thread
 * it assumes that the object that calls it already has the pop initialized
 */
template<class T>
svt_population<T> svt_gacylinder<T>::execute()
{
  //clear arrays/lists if runs were already executed
  (*this).m_oGA_Array.clear();
  svt_ga<T>::delTabuRegions();
  (*this).setDone(false);
  m_oCylinders.clear();
  m_oCylindersPdb.clear();
  m_oExploredVol.setValue(0.0);

  sprintf(m_pRmsdOut, " ");

  initAngles();
  initTemplate();
  initKernels();

  vector< svt_gacylinder * > oGA_Array;
  for (unsigned int iThread = 0; iThread < (*this).m_iMaxThread; iThread++) {
    svt_gacylinder *pGA_Tmp = (*this).createObject() ;
    pGA_Tmp->setThread(iThread);
    pGA_Tmp->setRun((*this).m_iRun);
    pGA_Tmp->setMaxGen((*this).m_iSyncGen);
    oGA_Array.push_back(pGA_Tmp);
  }

  for (unsigned int iThread = 0; iThread < (*this).m_iMaxThread; iThread++)
    (*this).m_oGA_Array.push_back(oGA_Array[iThread]);

  long iStartTime = svt_getElapsedTime();
  unsigned int iThreadGenerations = 0;

#ifdef __GAFIT_FORK
  vector<pid_t> oPid;
  pid_t iPid;
#endif
  bool bDone = false;
  svt_population<T> oFinalPop, oTotalPop, oTotalTabuPop;
  int iCount2Short = 0; //counts how many times too short cylinders were found
  Real64 fAvgLength;

  for (int iParallelRun = 0; iParallelRun < (*this).m_iMaxGen / (*this).m_iSyncGen && bDone == false; iParallelRun++) {
    this->m_iParallelRun = iParallelRun;

    SVTLBO << "Synchronization " << iParallelRun + 1 << " out of max. " << (int)floor((*this).m_iMaxGen / (*this).m_iSyncGen) << ".";
    if ((*this).m_iMaxThread != 1) cout << " Starting " << (*this).m_iMaxThread << " threads..." << endl;
    else cout << " Starting 1 serial thread..." << endl;

    (*this).initPopulation((*this).m_iPopSize);

    // start the threads
    for (unsigned int iThread = 0; iThread < (*this).m_iMaxThread; iThread++) {
      (*this).m_oGA_Array[iThread]->setCurrGen(0);
      ((svt_gacylinder<T> *)(*this).m_oGA_Array[iThread])->setFoundCylinder(true);
      (*this).m_oGA_Array[iThread]->setPopulation((*this).m_oPop);

      (*this).m_oGA_Array[iThread]->setDone(false);
#ifndef __GAFIT_FORK
      (*this). m_oGA_Array[iThread]->initThread();
#else
      if ((iPid = fork()) == 0) {
        oPid.push_back(iPid);
        (*this).m_oGA_Array[iThread]->run();
        exit(0);
      }
#endif
    }

#ifdef __GAFIT_FORK
    while (oPid.size() < (unsigned int)m_iMaxThread);
#endif

    // make sure threads are really running
    bool bSomeNotRunning = true;
    while (bSomeNotRunning) {
      bSomeNotRunning = false;
      for (unsigned int iThread = 0; iThread < (*this).m_iMaxThread; iThread++)
        if (!(*this).m_oGA_Array[iThread]->getIsThreadRunning())
          bSomeNotRunning = true;
    }

    // wait until they are finished
    bool bAllFinished = false;
    while (!bAllFinished) {
      bAllFinished = true;
      for (unsigned int iThread = 0; iThread < (*this).m_iMaxThread && bDone == false; iThread++) {
        if (!(*this).m_oGA_Array[iThread]->getDone())
          bAllFinished = false;
        else if ((*this).m_oGA_Array[iThread]->getCurrGen() < (*this).m_iSyncGen) {
          bDone = true;
          for (unsigned int i = 0; i < (*this).m_iMaxThread; i++)
            if (i != iThread)
              (*this).m_oGA_Array[i]->setDone(true);
          iThreadGenerations += (*this).m_oGA_Array[iThread]->getCurrGen();
        }
      }

      if (!bAllFinished)
        svt_sleep(500);
    }

    // make sure threads are really not running anymore
    bool bSomeRunning = true;
    while (bSomeRunning) {
      bSomeRunning = false;
      for (unsigned int iThread = 0; iThread < (*this).m_iMaxThread; iThread++)
        if ((*this).m_oGA_Array[iThread]->getIsThreadRunning())
          bSomeRunning = true;
    }

    //delete the content of the populations and start again
    oTotalPop.clear();
    oTotalTabuPop.clear();

    // combine populations
    for (unsigned int iThread = 0; iThread < (*this).m_iMaxThread; iThread++) {
      svt_population<T> oTempPop = (*this).m_oGA_Array[iThread]->getPopulation();
      svt_population<T> oTabuTempPop = (*this).m_oGA_Array[iThread]->getTabuRegions();

      oTotalPop.insert(oTotalPop.begin(), oTempPop.begin(), oTempPop.end());
      oTotalTabuPop.insert(oTotalTabuPop.begin(), oTabuTempPop.begin(), oTabuTempPop.end());
    }

    // keep track of how many generations we have so far
    if (!bDone)
      iThreadGenerations += (*this).m_iSyncGen;

    sort(oTotalPop.begin(), oTotalPop.end());
    (*this).m_oGA_Array[0]->penalizeSimilar(oTotalPop, 2.0 * (*this).getDistanceThreshold(), (1 - 2.0 * (1 - (*this).getDistanceThresholdPenalty()))); //penalty twice larger than usual
    sort(oTotalPop.begin(), oTotalPop.end());

    sort(oTotalTabuPop.begin(), oTotalTabuPop.end());
    (*this).m_oGA_Array[0]->penalizeSimilar(oTotalTabuPop, 0.0, 1.0);
    (*this).m_oGA_Array[0]->discardNullInd(oTotalTabuPop);
    sort(oTotalTabuPop.begin(), oTotalTabuPop.end());

    (*this).m_oTabus = oTotalTabuPop;

    if (!bDone) {
      for (unsigned int iThread = 0; iThread < (*this).m_iMaxThread; iThread++) {
        (*this).m_oGA_Array[iThread]->setTabuRegions(oTotalTabuPop);
        if (iThread == 0)
          (*this). m_oGA_Array[iThread]->outputResult(true);
      }
    }

    // discard the cylinders that were found in this parallel run
    fAvgLength = discardCylinders();
    if (fAvgLength < 2)
      iCount2Short++;
    else
      iCount2Short = 0;

    if (iCount2Short >= 4) {
      bDone = true; // only short helices were found in the last runs
      SVTLBO << "Stopping: only short traces were found!" << endl;
    }

    m_bCanOutputResults = true;

    sortCylinders();
    updateResults();

    if (m_fExploredPercent >= 85) {
      bDone = true;
      SVTLBO << "Stopping: more than 85 percent of the map was covered!" << endl;
    }

    if (m_iNoIndInInit / (Real64)this->m_iPopSize  >= 40) { // space is full with tabu regions
      bDone = true;
      SVTLBO << "Stopping: search templates can't be placed within map as space is full with tabu regions!" << endl;
    }

    // found 3 times more than requested
    if (m_iNoOfCylinder2Detect > 0 && m_oCylinders.size() > m_iNoOfCylinder2Detect * 3.0) { // found 3 times more than requested
      bDone = true;
      SVTLBO << "Stopping: found more than three times the requested " << m_iNoOfCylinder2Detect << " traces!" << endl;
    }
  }

  sortCylinders();
  discardCylinders(false);
  sortCylinders();
  SVTLBO << "A total of " << m_oCylinders.size() <<  " unique traces were found" << endl;
  updateResults(m_iNoOfCylinder2Detect);
  updateResults(m_oCylinders.size(), 1);

  //Output runtime stats
  long iETime = svt_getElapsedTime() - iStartTime;
  SVTLBO << "Number of generations (per thread): " << iThreadGenerations <<  endl;
  SVTLBO << "Elapsed time: " << ((iETime / (1000.0f)) / 60.0f) << " min" <<  endl;

  //clear the content of the threads
  this->clearThreads();

  return oFinalPop;
}


///////////////////////////////////////////////////////////////////////////////
// Parameters
///////////////////////////////////////////////////////////////////////////////

/**
 * Set resolution of the target map
 * \param fRes the resolution of the target map
 */
template<class T>
void svt_gacylinder<T>::setResolution(Real64 fRes)
{
  m_fRes = fRes;
};
/**
 * Get resolution of the target map
 * \return the resolution of the target map
 */
template<class T>
Real64 svt_gacylinder<T>::getResolution()
{
  return m_fRes;
};

/**
 * Set the masking when computing the correlation
 * \param bMask bool indicating whether to mask
 */
template<class T>
void svt_gacylinder<T>::setMask(bool bMask)
{
  m_bMask = bMask;
};

/**
 * Set the angular step size. In general: The smaller the value, the better the accuracy. In contrast to an exhaustive search, the runtime will also not be longer, if a finer step size
 * is chosen. The only limitation is the memory - as the table needs to be stored, very low numbers might result in an excessive use of memory. Recommended: ~0.5 to 1.0.
 * \param fDelta the angular step size (default 0.5)
 */
template<class T>
void svt_gacylinder<T>::setAngularStepSize(Real64 fDelta)
{
  m_fDelta = fDelta;
};
/**
 * Get the angular step size. In general: The smaller the value, the better the accuracy. In contrast to an exhaustive search, the runtime will also not be longer, if a finer step size
 * is chosen. The only limitation is the memory - as the table needs to be stored, very low numbers might result in an excessive use of memory. Recommended: ~0.5 to 1.0.
 * \return the angular step size (default 0.5)
 */
template<class T>
Real64 svt_gacylinder<T>::getAngularStepSize()
{
  return m_fDelta;
};

/**
 * Returns the number of angles
 */
template<class T>
Real64 svt_gacylinder<T>::getAnglesCount()
{
  return m_oAngles.getAngleCount();
};

/**
 * Get angles
 */
template<class T>
svt_vector4<Real64> &svt_gacylinder<T>::getAngle(long unsigned int iIndex)
{
  svt_vector4<Real64> oAngle;

  if (iIndex >= 0 && iIndex < m_oAngles.getAngleCount()) {
    oAngle.x(m_oAngles.getPsi(iIndex));
    oAngle.y(m_oAngles.getTheta(iIndex));
    oAngle.z(m_oAngles.getPhi(iIndex));
  } else {
    SVTLBO << "Angle out of boundaries! Exiting.." << endl;
  }

  return oAngle;
};


/**
 * Set the ranges for the angular search.
 * \param fPsiFrom   lower limit of the psi angles
 * \param fPsiTo     upper limit of the psi angles
 * \param fThetaFrom lower limit of the theta angles
 * \param fThetaTo   upper limit of the theta angles
 * \param fPhiFrom   lower limit of the phi angles
 * \param fPhiTo     upper limit of the phi angles
 */
template<class T>
void svt_gacylinder<T>::setAngularSearchRange(Real64 fPsiFrom, Real64 fPsiTo, Real64 fThetaFrom, Real64 fThetaTo, Real64 fPhiFrom, Real64 fPhiTo)
{
  m_fPsiFrom   = fPsiFrom;
  m_fPsiTo     = fPsiTo;
  m_fThetaFrom = fThetaFrom;
  m_fThetaTo   = fThetaTo;
  m_fPhiFrom   = fPhiFrom;
  m_fPhiTo     = fPhiTo;
};


/**
 * Set the ranges for the translational search relative to the centers of the units
 * (if fXFrom = -20 and fXTo = 20 the units moves -20 and 20 A from the current position)
 * \param fXFrom     upper limit of the X
 * \param fXTo       lower limit of the X
 * \param fYFrom     upper limit of the Y
 * \param fYTo       lower limit of the Y
 * \param fZFrom     upper limit of the Z
 * \param fZTo       lower limit of the Z
 */
template<class T>
void svt_gacylinder<T>::setRelativeTranslSearchRange(Real64 fXFrom, Real64 fXTo, Real64 fYFrom, Real64 fYTo, Real64 fZFrom, Real64 fZTo)
{


  if ((fXFrom != fXTo || fYFrom != fYTo || fZFrom != m_fZTo)) { // the values were already initialized
    m_fXFrom    = fXFrom + m_fBorder;
    m_fXTo      = fXTo - m_fBorder;
    m_fYFrom    = fYFrom + m_fBorder;
    m_fYTo      = fYTo - m_fBorder;
    m_fZFrom    = fZFrom + m_fBorder;
    m_fZTo      = fZTo - m_fBorder;
  }
};

///////////////////////////////////////////////////////////////////////////////
// Units, Map, PDB files...
///////////////////////////////////////////////////////////////////////////////
/**
 * set the tempate of the cylinder
 * \param oPdb the target structure
 */
template<class T>
void svt_gacylinder<T>::setTemplate(svt_point_cloud_pdb<svt_ga_vec>  &oPdb)
{
  m_oTurn = oPdb;

  for (unsigned int iIndex = 0; iIndex < m_oTemplate.size(); iIndex++) {
    m_oTurn.getAtom(iIndex)->setName("C");
    m_oTurn.getAtom(iIndex)->setRemoteness('A');
    m_oTurn.getAtom(iIndex)->adjustMass();
  }


  svt_ga_vec oCoa = m_oTurn.coa();

  svt_matrix4<Real64> oMat;
  oMat.loadIdentity();
  oMat.setTranslation(-oCoa);

  m_oTurn = oMat * m_oTurn;
}


/**
 * set target
 * \param oTar the map
 */
template<class T>
void svt_gacylinder<T>::setTarget(svt_ga_vol &oTar)
{
  m_oTar = oTar;

  m_iSizeX = m_oTar.getSizeX();
  m_iSizeY = m_oTar.getSizeY();
  m_iSizeZ = m_oTar.getSizeZ();

  m_oSearchSpace.x(m_iSizeX);
  m_oSearchSpace.y(m_iSizeY);
  m_oSearchSpace.z(m_iSizeZ);

  m_fOrigX = m_oTar.getGridX();
  m_fOrigY = m_oTar.getGridY();
  m_fOrigZ = m_oTar.getGridZ();

  m_fWidth = m_oTar.getWidth();

  m_oModelVol.allocate(m_iSizeX, m_iSizeY, m_iSizeZ, 0.0f);
  m_oModelVol.setGrid(m_fOrigX, m_fOrigY, m_fOrigZ);
  m_oModelVol.setWidth(m_fWidth);

  m_fXFrom    = m_fOrigX + m_fBorder;
  m_fXTo  = m_fOrigX + m_iSizeX * m_fWidth - m_fBorder;
  m_fYFrom    = m_fOrigY + m_fBorder;
  m_fYTo  = m_fOrigY + m_iSizeY * m_fWidth - m_fBorder;
  m_fZFrom    = m_fOrigZ + m_fBorder;
  m_fZTo  = m_fOrigZ + m_iSizeZ * m_fWidth - m_fBorder;

  if (m_fRes != 0)
    initKernels();

  //set the exlored volume to null
  m_oExploredVol = m_oTar;
  m_oExploredVol.setValue(0.0);

  //
  // Output input
  //
  char pFname[256];
  if (strlen(m_pPath) != 0) {
    sprintf(pFname, "%s/aaa.situs", m_pPath);
    //m_oTar.save( pFname );
  }
};

/**
 * init blurring kernels
 */
template<class T>
void svt_gacylinder<T>::initKernels()
{
  m_oKernel.create1DBlurringKernel(m_fWidth, m_fRes / 2.0);
};

/**
 * set the accept Move score percentage:e.g. 0.95 means that moves are accepted if withing 0.95 of the original score
 * \param the score percentage
 */
template<class T>
void svt_gacylinder<T>::setAcceptMoveRatio(Real64 fAcceptMoveRatio)
{
  m_fAcceptMoveRatio = fAcceptMoveRatio;
};

/**
 * get the accept Move score percentage: e.g. 0.95 means that moves are accepted if withing 0.95 of the original score
 * \param the score percentage
 */
template<class T>
Real64 svt_gacylinder<T>::getAcceptMoveRatio()
{
  return m_fAcceptMoveRatio;
};


/**
 * set the max number of failed Crawls before stoping the search
 * \param iMaxFailedCrawls e.g. 2 indicates 2 times tried before stoping the crawl
 */
template<class T>
void svt_gacylinder<T>::setMaxFailedCrawls(unsigned int iMaxFailedCrawls)
{
  m_iMaxFailedCrawls = iMaxFailedCrawls;
};

/**
 * get the max number of failed Crawls before stoping the search
 * \return iMaxFailedCrawls e.g. 2 indicates 2 times tried before stoping the crawl
 */
template<class T>
unsigned int svt_gacylinder<T>::getMaxFailedCrawls()
{
  return m_iMaxFailedCrawls;
};

/**
 * how many Cylinders shouls be searched for
 */
template<class T>
void svt_gacylinder<T>::setNoOfCylinder2Detect(int iNoOfCylinder2Detect)
{
  m_iNoOfCylinder2Detect = iNoOfCylinder2Detect;
};

/**
 * how many Cylinders shouls be searched for
 */
template<class T>
int svt_gacylinder<T>::getNoOfCylinder2Detect()
{
  return m_iNoOfCylinder2Detect;
};

/**
 * was "valid" cylinder found
 */
template<class T>
void svt_gacylinder<T>::setFoundCylinder(bool bFoundCylinder)
{
  m_bFoundCylinder = bFoundCylinder;
}

/**
 *  was "valid" cylinder found
 */
template<class T>
bool svt_gacylinder<T>::wasFoundCylinder()
{
  return m_bFoundCylinder;
};

/**
 *  set apply Blurring to model
 * \param state
 */
template<class T>
void svt_gacylinder<T>::setApplyBlurring2Model(bool bBlurring)
{
  m_bApplyBlurring2Model = bBlurring;
};

/**
 * set the radius of the template
 * \param the radius of the template
 */
template<class T>
void svt_gacylinder<T>::setTemplateRadius(Real64 fTemplateRadius)
{
  m_fTemplateRadius = fTemplateRadius;
};

/**
 * get the radius of the template
 * \param the radius of the template
 */
template<class T>
Real64 svt_gacylinder<T>::getTemplateRadius()
{
  return m_fTemplateRadius;
};

/**
 * set the radius of the template
 * \param the radius of the template
 */
template<class T>
void svt_gacylinder<T>::setSearchTemplateRadius(Real64 fSearchTemplateRadius)
{
  m_fSearchTemplateRadius = fSearchTemplateRadius;
};

/**
 * get the radius of the template
 * \param the radius of the template
 */
template<class T>
Real64 svt_gacylinder<T>::getSearchTemplateRadius()
{
  return m_fSearchTemplateRadius;
};

/**
 * set the number of points in one of the circles of the template
 * \param  number of points in one of the circles of the template
 */
template<class T>
void svt_gacylinder<T>::setTemplatePointCount(unsigned int iTemplatePointCount)
{
  m_iTemplatePointCount = iTemplatePointCount;
};

/**
 * get the number of points in one of the circles of the template
 * \param  number of points in one of the circles of the template
 */
template<class T>
unsigned int svt_gacylinder<T>::getTemplatePointCount()
{
  return m_iTemplatePointCount;
};

/**
 * set the number of copies of the template
 * \param the number of repeats to be set
 */
template<class T>
void svt_gacylinder<T>::setTemplateRepeats(unsigned int iTemplateRepeats)
{
  m_iTemplateRepeats = iTemplateRepeats;
};

/**
 * get the number of copies of the circle in the search template
 * \return the number of repeats
 */
template<class T>
unsigned int svt_gacylinder<T>::getSearchTemplateRepeats()
{
  return m_iSearchTemplateRepeats;
};

/**
 * set the number of copies of the circle in the search template
 * \param the number of repeats to be set
 */
template<class T>
void svt_gacylinder<T>::setSearchTemplateRepeats(unsigned int iSearchTemplateRepeats)
{
  m_iSearchTemplateRepeats = iSearchTemplateRepeats;
};

/**
 * get the number of copies of the template
 * \return the number of repeats
 */
template<class T>
unsigned int svt_gacylinder<T>::getTemplateRepeats()
{
  return m_iTemplateRepeats;
};

/**
 * set the distance between two repeats
 * \param the distance between repeats
 */
template<class T>
void svt_gacylinder<T>::setDistBetweenRepeats(Real64 fDistBetweenRepeats)
{
  m_fDistBetweenRepeats = fDistBetweenRepeats;
};

/**
 * get the distance between two repeats
 * \return the distance between repeats
 */
template<class T>
Real64 svt_gacylinder<T>::getDistBetweenRepeats()
{
  return m_fDistBetweenRepeats;
};

/**
 * set the size of a crawling step
 * \param the size of a crawling step
 */
template<class T>
void svt_gacylinder<T>::setCrawlingStepSize(Real64 fCrawlingStepSize)
{
  m_fCrawlingStepSize = fCrawlingStepSize;
  m_fStepsPerTurn = 5.4 / m_fCrawlingStepSize;
};

/**
* get the size of a crawling step
* \return the size of a crawling step
*/
template<class T>
Real64 svt_gacylinder<T>::getCrawlingStepSize()
{
  return m_fCrawlingStepSize;
};

/**
 * set fit high resolution helix
 */
template<class T>
void svt_gacylinder<T>::setFitHelix(bool bFitHelix)
{
  m_bFitHelix = bFitHelix;
};

/**
 * get fit high resolution helix
 */
template<class T>
bool svt_gacylinder<T>::getFitHelix()
{
  return m_bFitHelix;
};


/**
 * are results available and can be outputed
 */
template<class T>
bool svt_gacylinder<T>::canOutputResults()
{
  return m_bCanOutputResults;
};

/**
 * set the outputTemplate option
 * \param the bOutputTemplates
 */
template<class T>
void svt_gacylinder<T>::setOutputTemplates(bool bOutputTemplates)
{
  m_bOutputTemplates = bOutputTemplates;
};

/**
 * get the outputTemplate option
 * \return OutputTemplates
 */
template<class T>
bool svt_gacylinder<T>::getOutputTemplates()
{
  return m_bOutputTemplates;
};

/**
 * set anisotropic correction related parameters
 * \param fAni the correction, e.g. fAni=2 will compress the map 2 times
 * \param fOrigZWoAni the origin of the map before ani correction - use for decompression
 */
template<class T>
void svt_gacylinder<T>::setAniCorr(Real64 fAni, Real64 fOrigZWoAni)
{
  m_fAni = fAni;
  m_fOrigZWoAni = fOrigZWoAni;
};

/**
 * get anisotropic correction
 * \param fAni the correction, e.g. fAni=2 will compress the map 2 times
 */
template<class T>
Real64 svt_gacylinder<T>::getAniCorr()
{
  return m_fAni;
};

/**
 * get Origin on Z axis of the map without anisotropic correction
 */
template<class T>
Real64 svt_gacylinder<T>::getOrigZWoAni()
{
  return m_fOrigZWoAni;
};

/**
 * get Origin on Z axis of the map with anisotropic correction
 */
template<class T>
Real64 svt_gacylinder<T>::getOrigZ()
{
  return m_fOrigZ;
};

/**
 * Set the correction factor lambda
 */
template<class T>
void svt_gacylinder<T>::setLambda(Real64 fLambda)
{
  m_fLambda = fLambda;
};

/**
 * Get the correction factor lambda
 */
template<class T>
Real64 svt_gacylinder<T>::getLambda()
{
  return m_fLambda;
};

///////////////////////////////////////////////////////////////////////////////
// Routines to calculate the fitness
///////////////////////////////////////////////////////////////////////////////
/**
 * returns the transformation matrix
 */
template<class T>
svt_ga_mat svt_gacylinder<T>::getTrans(T *pInd)
{
  this->makeValid(pInd);

  //translation
  svt_ga_vec oTransl;
  oTransl.x(m_fXFrom + pInd->getGene(0) * (m_fXTo - m_fXFrom));
  oTransl.y(m_fYFrom + pInd->getGene(1) * (m_fYTo - m_fYFrom));
  oTransl.z(m_fZFrom + pInd->getGene(2) * (m_fZTo - m_fZFrom));

  //Rotation
  Real64 fPsi     = m_oAngles.getPsi((long unsigned int)(pInd->getGene(3) * m_oAngles.getAngleCount()));
  Real64 fTheta   = m_oAngles.getTheta((long unsigned int)(pInd->getGene(3) * m_oAngles.getAngleCount()));
  Real64 fPhi     = m_oAngles.getPhi((long unsigned int)(pInd->getGene(3) * m_oAngles.getAngleCount()));

  svt_matrix4<Real64> oMat;
  oMat.loadIdentity();
  oMat.rotatePTP(fPhi, fTheta, fPsi);
  oMat.setTranslation(oTransl);

  pInd->setTrans(oMat);

  return oMat;
}


/**
 * Updates the coordinates of the model based on the genes of the individual
 * attention: it does not update the volume of the model (see update volume)
 */
template<class T>
void svt_gacylinder<T>::updateModel(T *pInd)
{
  svt_matrix4<Real64> oMat;
  oMat.loadIdentity();

  oMat = getTrans(pInd);

  if (m_bRefining) {
    if (m_oModel.size() != m_oTemplate.size())
      m_oModel = m_oTemplate;

    m_oModel = oMat * m_oTemplate;
  } else {
    if (m_oModel.size() != m_oBigTemplate.size())
      m_oModel = m_oBigTemplate;

    m_oModel = oMat * m_oBigTemplate;
  }

  pInd->updateCoarsePhenotype(oMat);
};

/**
 * Updates the volume of the model
 * Attention: it does not update the model - the pdb remains the same
 */
template<class T>
void svt_gacylinder<T>::updateVolume(T *pInd)
{
  updateModel(pInd);

  //svt_ga_mat oMat = getTrans(pInd);
  m_oModelVol.setValue(0.0f);
  m_oModel.projectMass(&m_oModelVol);

  if (m_bApplyBlurring2Model)
    m_oModelVol.convolve1D3D(m_oKernel, false);
};

/**
 * update fitness
 * \param pInd pointer to individual that gets updated
 */
template<class T>
void svt_gacylinder<T>::updateFitness(T *pInd)
{
  svt_ga_mat oMat;
  updateModel(pInd);

  Real64 fFitness, fSin;
  fFitness = m_oModel.projectMassCorr(&m_oTar, oMat, false);
  if (m_fLambda > 1 || m_fLambda < 0) // the set lambda is not between 0 and 1
    pInd->setFitness(fFitness);
  else {
    // code uses Goldstein convention for Euler angles, this is different from tomo paper where Theta is measured from (X,Y) plane
    // hence we replace cos^2 in tomo paper with sin^2 to give lambda correction
    fSin = sin(m_oAngles.getTheta((long unsigned int)(pInd->getGene(3) * m_oAngles.getAngleCount())));
    pInd->setFitness((m_fLambda + (1 - m_fLambda)*fSin * fSin)*fFitness);
  }


  this->m_iFitnessUpdateCount++;
}

/**
 * add 2 Cylinder
 */
template<class T>
void svt_gacylinder<T>::add2Tube(T *pInd, bool bAddNew, bool bFlip)
{
  if (bAddNew) {
    svt_tube oCylinder;
    m_oCylinders.push_back(oCylinder);

    //delete all the tabus in the temporary tabu list
    m_oTabusInCrawl.clear();
  }

  if (m_oCylinders.size() > 0) {
    svt_ga_mat oMat = pInd->getTrans();
    m_oCylinders[m_oCylinders.size() - 1].add2Tube(*pInd, bFlip);

  }

  //add to tabu list tabu
  m_oTabusInCrawl.push_back(*pInd);
}


/**
 * get the cylinderes
 * \return the vector containing the tubes
 */
template<class T>
vector<svt_tube> &svt_gacylinder<T>::getCylinders()
{
  return m_oCylinders;
};

/**
 * set the cylinderes
 * \param the vector containing the tubes
 */
template<class T>
void svt_gacylinder<T>::setCylinders(vector<svt_tube>   &oCylinders)
{
  m_oCylinders = oCylinders;
};
/**
 * set the cylinderes
 * \param the vector containing the helices
 */
template<class T>
void svt_gacylinder<T>::setCylindersPdb(vector< svt_point_cloud_pdb<svt_ga_vec> > &oCylinders)
{
  m_oCylindersPdb = oCylinders;
};

/**
 * get the cylinderes
 * \return the vector containing the helices
 */
template<class T>
vector< svt_point_cloud_pdb<svt_ga_vec> > &svt_gacylinder<T>::getCylindersPdb()
{
  return m_oCylindersPdb;
};

/**
 * refine an individual; by default it does nothing - overload in the classes that inherit
 * \param the individual that will be refined
 */
template<class T>
void svt_gacylinder<T>::refineInd(T *pInd)
{
  refine(pInd, 0);
};

/**
 * refine an individual; by default it does nothing - overload in the classes that inherit
 * \param the individual that will be refined
 * \param iteration
 */
template<class T>
void svt_gacylinder<T>::refine(T *pInd, int iIter)
{
  m_bRefining = true;
  updateFitness(pInd);

  Real64 fFitness = pInd->getFitness();

  //refine the rotations randomly
  refineRotRandom(pInd);

  //refine the rotation and translation a few times
  refineGenes(pInd);

  //refine the rotations randomly
  refineRotRandom(pInd);

  //refine the rotation and translation a few times
  refineGenes(pInd);

  add2Tube(pInd, true);

  //refine the length
  refineInLength(pInd);

  //m_bMask = false;
  m_bRefining = false;
  updateFitness(pInd);

  if (m_oTabusInCrawl.size() > 1) {
    // found better fitness while crawling and didn't refine already 2 times => keep refining
    if (m_oTabusInCrawl[0].getFitness() > fFitness && iIter < 2 && !this->isInTabuReg(&m_oTabusInCrawl[0])) {
      //get better Ind
      pInd->setGenes(m_oTabusInCrawl[0].getGenes());

      //delete current cylinder
      m_oCylinders.pop_back();
      m_oTabusInCrawl.clear();

      refine(pInd, iIter + 1);
      return;
    }

    sort(m_oTabusInCrawl.rbegin(), m_oTabusInCrawl.rend());
    char pFname[265];
    sprintf(pFname, "%s/Tabu%02d.pdb", m_pPath, this->m_iRun);
    svt_point_cloud_pdb<svt_ga_vec> oPdb;

    for (unsigned int iIndex = 0; iIndex < m_oTabusInCrawl.size(); iIndex++) {
      this->m_oTabus.push_back(m_oTabusInCrawl[iIndex]);
      /*oPdb = m_oTabusInCrawl[iIndex].getCoarsePhenotype();
      if (strlen(m_pPath)!=0)
          oPdb.writePDB( pFname, true );*/
    }
    // char pOut[256];
    // sprintf(pOut, "[%02d-%04d] Added %d tabu regions. ", this->m_iThread, this->m_iGenerations, (int)m_oTabusInCrawl.size()); SVTLBO << pOut << endl;

    m_bFoundCylinder = true;
  } else {
    m_bFoundCylinder = false;
    m_oCylinders.pop_back();
  }

  m_oTabusInCrawl.clear();
  sort(this->m_oTabus.begin(), this->m_oTabus.end());
  this->m_oTabuWindow.clear();
  //done in this generatation
  (*this).m_iGenerations = (*this).m_iMaxGen + 1; // stop run at the end of the refinement

};

/**
 * refine the translation and rotation of an individual - calls the refinetransl and refineRot a few times
 * \param the individual that will be refined
 */
template<class T>
void svt_gacylinder<T>::refineGenes(T *pInd, svt_vector4<Real64> *pCenter)
{
  //refine the translation
  if (pCenter == NULL)
    refineTransl(pInd);
  else
    refineTranslOnSphere(pInd, pCenter);

  //refine rotation
  refineRot(pInd);
};


/**
 * refine the translation of an individual
 * \param the individual that will be refined
 */
template<class T>
void svt_gacylinder<T>::refineTransl(T *pInd)
{
  char    pFname[1024];
  bool    bFound      = false;
  Real64  fX          = 0;
  Real64  fY          = 0;
  Real64  fZ          = 0;
  Real64  fAddGeneX   = 1.0 / (8.0 * Real64(m_iSizeX) * m_fWidth); //half angstron
  Real64  fAddGeneY   = 1.0 / (8.0 * Real64(m_iSizeY) * m_fWidth);
  Real64  fAddGeneZ   = 1.0 / (8.0 * Real64(m_iSizeZ) * m_fWidth);

  T      *pNewInd = new T(*pInd);
  updateFitness(pNewInd);
  Real64  fMax    = pNewInd->getFitness();

  for (int iIndexX = -this->m_iRefinementMaxMutPerGene * 2; iIndexX <= this->m_iRefinementMaxMutPerGene * 2; iIndexX++) {
    pNewInd->setGene(0, pInd->getGene(0) + fAddGeneX * iIndexX);
    for (int iIndexY = -this->m_iRefinementMaxMutPerGene * 2; iIndexY <= this->m_iRefinementMaxMutPerGene * 2; iIndexY++) {
      pNewInd->setGene(1, pInd->getGene(1) + fAddGeneY * iIndexY);
      for (int iIndexZ = -this->m_iRefinementMaxMutPerGene * 2; iIndexZ <= this->m_iRefinementMaxMutPerGene * 2; iIndexZ++) {
        pNewInd->setGene(2, pInd->getGene(2) + fAddGeneZ * iIndexZ);

        if (pNewInd->getGene(0) <= 1 && pNewInd->getGene(0) >= 0 &&
            pNewInd->getGene(1) <= 1 && pNewInd->getGene(1) >= 0 &&
            pNewInd->getGene(2) <= 1 && pNewInd->getGene(2) >= 0) {
          this->makeValid(pNewInd);
          updateFitness(pNewInd);

          if (fMax < pNewInd->getFitness()) {
            bFound = true;
            fMax = pNewInd->getFitness();
            fX = pNewInd->getGene(0);
            fY = pNewInd->getGene(1);
            fZ = pNewInd->getGene(2);

            if (strlen(m_pPath) != 0) {
              sprintf(pFname, "%s/ModelTabuRef%02d%02d.pdb", m_pPath, this->m_iRun, this->m_iThread);
              //m_oModel.writePDB( pFname, true );
            }
          }
        }
      }
    }
  }

  if (bFound) { // did found a better score
    pInd->setGene(0, fX);
    pInd->setGene(1, fY);
    pInd->setGene(2, fZ);
    updateFitness(pInd);
  }

  delete pNewInd;
}

/**
 * refine the translation of an individual, only allow movements on Spheres
 * \param the individual that will be refined
 */
template<class T>
void svt_gacylinder<T>::refineTranslOnSphere(T *pInd, svt_vector4<Real64> *pCenter)
{
  char    pFname[1024];
  bool    bFound  = false;
  Real64  fX      = 0;
  Real64  fY      = 0;
  Real64  fZ      = 0;

  T      *pNewInd = new T(*pInd);
  updateFitness(pNewInd);
  Real64  fMax    = pNewInd->getFitness();

  svt_ga_vec oTransl;
  oTransl.x(m_fXFrom + pNewInd->getGene(0) * (m_fXTo - m_fXFrom));
  oTransl.y(m_fYFrom + pNewInd->getGene(1) * (m_fYTo - m_fYFrom));
  oTransl.z(m_fZFrom + pNewInd->getGene(2) * (m_fZTo - m_fZFrom));

  svt_ga_vec oAxis    = oTransl - (*pCenter);
  Real64 fR           = oAxis.length();
  Real64 fNewR, fNewTheta, fNewPhi;
  svt_ga_vec oNewCenter;

  svt_point_cloud_pdb<svt_ga_vec> oPdb, oCenters;
  svt_point_cloud_atom oAtom;
  oPdb = pNewInd->getCoarsePhenotype();
  oCenters.addAtom(oAtom, *pCenter);
  oCenters.addAtom(oAtom, oPdb[0]);

  Real64 fPsi         = m_oAngles.getPsi((long unsigned int)(pNewInd->getGene(3) * m_oAngles.getAngleCount()));
  Real64 fTheta       = m_oAngles.getTheta((long unsigned int)(pNewInd->getGene(3) * m_oAngles.getAngleCount()));
  Real64 fPhi         = m_oAngles.getPhi((long unsigned int)(pNewInd->getGene(3) * m_oAngles.getAngleCount()));

  svt_ga_mat oMat;
  oMat.loadIdentity();
  oMat.rotatePTP(fPhi, fTheta, fPsi);

  fNewR = fR; // add 0.5 angstrom
  for (fNewTheta = deg2rad(5.00); fNewTheta <= deg2rad(10.0); fNewTheta += deg2rad(5.00)) {
    //fNewTheta = deg2rad(20.0)*iIndexY; // 20 deg
    for (fNewPhi = 0; fNewPhi < deg2rad(360); fNewPhi += deg2rad(45.0)) {
      //fNewPhi = deg2rad(45.0)*iIndexZ; // 45 deg

      oNewCenter.x(fNewR * cos(fNewPhi) * sin(fNewTheta));
      oNewCenter.y(fNewR * sin(fNewPhi) * sin(fNewTheta));
      oNewCenter.z(fNewR * cos(fNewTheta));

      oNewCenter = *pCenter + oMat * oNewCenter;

      pNewInd->setGene(0, (oNewCenter.x() - m_fXFrom) / (m_fXTo - m_fXFrom));
      pNewInd->setGene(1, (oNewCenter.y() - m_fYFrom) / (m_fYTo - m_fYFrom));
      pNewInd->setGene(2, (oNewCenter.z() - m_fZFrom) / (m_fZTo - m_fZFrom));

      if (pNewInd->getGene(0) <= 1 && pNewInd->getGene(0) >= 0 && pNewInd->getGene(1) <= 1 && pNewInd->getGene(1) >= 0 && pNewInd->getGene(2) <= 1 && pNewInd->getGene(2) >= 0) {

        this->makeValid(pNewInd);

        updateFitness(pNewInd);

        oPdb = pNewInd->getCoarsePhenotype();
        oCenters.addAtom(oAtom, oPdb[0]);

        if (fMax < pNewInd->getFitness()) {
          bFound = true;
          fMax = pNewInd->getFitness();
          fX = pNewInd->getGene(0);
          fY = pNewInd->getGene(1);
          fZ = pNewInd->getGene(2);
          if (strlen(m_pPath) != 0) {
            sprintf(pFname, "%s/ModelTabuRefCylinder%02d%02d.pdb", m_pPath, this->m_iRun, this->m_iThread);
            //m_oModel.writePDB( pFname, true );
          }
        }
      }
    }
  }

  if (bFound) { // did found a better score
    pInd->setGene(0, fX);
    pInd->setGene(1, fY);
    pInd->setGene(2, fZ);
    updateFitness(pInd);
  }

  delete pNewInd;
}

/**
 * refine the translation of an individual
 * \param the individual that will be refined
 */
template<class T>
void svt_gacylinder<T>::refineRot(T *pInd, Real64 fMaxCorr)
{
  char pFname[1024];
  bool bFound = false;
  long int iMaxIndex  = -1;

  T      *pNewInd = new T(*pInd);
  updateFitness(pNewInd);
  Real64  fMax    = pNewInd->getFitness();

  svt_ga_vec oVec, oNull;
  oVec.x(0.0f);
  oVec.y(0.0f);
  oVec.z(1.0f);
  oNull.x(0.0f);
  oNull.y(0.0f);
  oNull.z(0.0f);

  svt_ga_mat oMat;
  oMat = pInd->getTrans();
  oMat.setTranslation(oNull);

  svt_ga_vec oInitAxis = oVec * oMat;
  svt_ga_vec oCurrAxis;
  Real64 fDotProd;
  int iCount = 0;

  for (long unsigned int iAngle = 0; iAngle < m_oAngles.getAngleCount(); iAngle++) {
    Real64 fPsi     = m_oAngles.getPsi(iAngle);
    Real64 fTheta   = m_oAngles.getTheta(iAngle);
    Real64 fPhi     = m_oAngles.getPhi(iAngle);

    oMat.loadIdentity();
    oMat.rotatePTP(fPhi, fTheta, fPsi);
    oCurrAxis = oVec * oMat;

    fDotProd = oInitAxis.x() * oCurrAxis.x() + oInitAxis.y() * oCurrAxis.y() + oInitAxis.z() * oCurrAxis.z();

    if (fDotProd > fMaxCorr) {
      iCount++;
      if (iCount % 2 == 0) { // compute fitness every 5 angles

        pNewInd->setGene(3, iAngle / (Real64)m_oAngles.getAngleCount());
        updateFitness(pNewInd);

        if (fMax < pNewInd->getFitness()) {
          bFound = true;
          fMax = pNewInd->getFitness();
          iMaxIndex = iAngle;
          if (strlen(m_pPath) != 0) {
            sprintf(pFname, "%s/ModelTabuRotBetter%02d%02d.pdb", m_pPath, this->m_iRun, this->m_iThread);
            //m_oModel.writePDB( pFname, true );
          }
        }
      }
    }
  }

  if (bFound) {
    pInd->setGene(3, iMaxIndex / (Real64)m_oAngles.getAngleCount());
    updateFitness(pInd);
  }

  delete pNewInd;
};

/**
 * refine the rotation of an individual using random angles
 * \param the individual that will be refined
 */
template<class T>
void svt_gacylinder<T>::refineRotRandom(T *pInd)
{
  char pFname[1024];
  bool bFound = false;
  Real64 fGene = 0;

  T      *pNewInd = new T(*pInd);
  updateFitness(pNewInd);
  Real64  fMax    = pNewInd->getFitness();

  if (strlen(m_pPath) != 0) {
    sprintf(pFname, "%s/ModelTabuRotBetter%02d%02d.pdb", m_pPath, this->m_iRun, this->m_iThread);
    //m_oModel.writePDB( pFname, true );
  }

  //randomly sample 10% of the angles
  for (long unsigned int iAngle = 0; iAngle < m_oAngles.getAngleCount()/*0.10*/; iAngle++) {
    //pNewInd->setGene( 3, svt_genrand() );
    pNewInd->setGene(3, iAngle / Real64(m_oAngles.getAngleCount()));
    updateFitness(pNewInd);

    if (fMax < pNewInd->getFitness()) {
      bFound = true;
      fMax = pNewInd->getFitness();
      fGene = pNewInd->getGene(3);

      if (strlen(m_pPath) != 0) {
        sprintf(pFname, "%s/ModelTabuRotBetter%02d%02d.pdb", m_pPath, this->m_iRun, this->m_iThread);
        //m_oModel.writePDB( pFname, true );
      }
    }
  }

  if (bFound) {
    pInd->setGene(3, fGene);
    updateFitness(pInd);

    if (strlen(m_pPath) != 0) {
      sprintf(pFname, "%s/ModelTabuRotBetter%02d%02d.pdb", m_pPath, this->m_iRun, this->m_iThread);
      //m_oModel.writePDB( pFname, true );
    }


  }

  delete pNewInd;
};

/**
 * crawl on the tube and search for similar scoring cylinder placements
 * \param the individual that will be refined
 * \param iDirection indicates the direction  +1 for forward and -1 for backwards
 */
template<class T>
void svt_gacylinder<T>::crawl(T *pInd, int iDirection)
{
  //should the elements in the list of accepted step should be flipped
  //by default no, except when changing direction
  //char pOut[1024];
  bool bFlip = false;

  if (iDirection == -1)
    bFlip = true;

  //char pOut[256];
  if (iDirection != 1 && iDirection != -1) {
    SVTLBO << "The direction can be either +1 (forward) or -1 (backwards)" << endl;
    return;
  }

  char pFname[1024];
  svt_point_cloud_pdb<svt_ga_vec> oPdb;
  T      *pNewInd         = new T(*pInd);
  updateFitness(pNewInd);
  Real64  fStartScore         = pNewInd->getFitness();

  //
  // Move on the cylinder in the direction indicated by the axis of the start cylinder
  //
  //center of the start cylinder
  svt_ga_vec oTransl;
  oTransl.x(m_fXFrom + pNewInd->getGene(0) * (m_fXTo - m_fXFrom));
  oTransl.y(m_fYFrom + pNewInd->getGene(1) * (m_fYTo - m_fYFrom));
  oTransl.z(m_fZFrom + pNewInd->getGene(2) * (m_fZTo - m_fZFrom));

  svt_vector4<Real64> oOldTrans;

  int     iNoFailureToImprove = 0; // number of steps that failed to improve; counted since last improvement
  while (iNoFailureToImprove < (int)m_iMaxFailedCrawls) {

    //translate in the direction of the axes of the previous cylinder
    oPdb = pNewInd->getCoarsePhenotype();
    if (iDirection == 1)
      oOldTrans = oTransl;
    else
      oOldTrans = oTransl + 2 * iDirection * m_fCrawlingStepSize * (oPdb[1] - oPdb[0]);

    //oTransl - is the translation of the previous cylinder (before refinement)
    oTransl = oTransl + iDirection * m_fCrawlingStepSize * (oPdb[1] - oPdb[0]);

    pNewInd->setGene(0, (oTransl.x() - m_fXFrom) / (m_fXTo - m_fXFrom));
    pNewInd->setGene(1, (oTransl.y() - m_fYFrom) / (m_fYTo - m_fYFrom));
    pNewInd->setGene(2, (oTransl.z() - m_fZFrom) / (m_fZTo - m_fZFrom));

    if (pNewInd->getGene(0) <= 1 && pNewInd->getGene(0) >= 0 &&
        pNewInd->getGene(1) <= 1 && pNewInd->getGene(1) >= 0 &&
        pNewInd->getGene(2) <= 1 && pNewInd->getGene(2) >= 0) {
      //refine rotation and translation
      refineGenes(pNewInd, &oOldTrans);
      updateFitness(pNewInd);

      // is the score good enough - then make tabu
      if (fStartScore * m_fAcceptMoveRatio < pNewInd->getFitness()) //restart the stop counter
        iNoFailureToImprove = 0;
      else //inc the no of failure
        iNoFailureToImprove++;

      //a good step or the first failure then add to cylinder
      if (iNoFailureToImprove <= 1) {
        if (!bFlip)
          add2Tube(pNewInd);
        else {
          add2Tube(pNewInd, false, bFlip);
          bFlip = false;
        }

        //print
        if (strlen(m_pPath) != 0) {
          sprintf(pFname, "%s/Model%02d%02d%02d.pdb", m_pPath, this->m_iRun, this->m_iThread, (int)this->m_oPop.size() - 1);
          //m_oModel.writePDB( pFname, true );
        }
      } else {
        //fail- remove the last two points that were added as they constitue failures
        if (m_oCylinders.size() > 1) {
          m_oCylinders[m_oCylinders.size() - 1].pop();
          m_oCylinders[m_oCylinders.size() - 1].pop();
        }
      }
    } else { // i'M outside map and should stop here
      iNoFailureToImprove = m_iMaxFailedCrawls;
    }
  }
  delete pNewInd;
}

/**
 * refine the length of an individual
 * \param the individual that will be refined
 */
template<class T>
void svt_gacylinder<T>::refineInLength(T *pInd)
{
  if (strlen(m_pPath) != 0) {
    char pFname[1256];
    updateFitness(pInd);
    sprintf(pFname, "%s/Model%02d%02d99.pdb", m_pPath, this->m_iRun, this->m_iThread);
    //m_oModel.writePDB( pFname, true );
  }

  // crawl forward
  T *pNewInd = new T(*pInd);
  crawl(pNewInd, +1);
  delete pNewInd;

  //crawl backworks
  pNewInd = new T(*pInd);
  crawl(pNewInd, -1);
  delete pNewInd;

}

/**
 * Calculate the full correlation corresponding to a certain individual (with blur)
 */
template<class T>
Real64 svt_gacylinder<T>::getCorrelation()
{
  m_oModelVol.convolve1D3D(m_oKernel, false);  // don't normalize
  return m_oModelVol.correlation(m_oTar, false);
}

/**
 * Calculate the rmsd of the individual towards the target model
 * \return the rmsd
 */
template<class T>
Real64 svt_gacylinder<T>::getRMSD()
{
  if (m_oTarStr.size() > 0) {
    Real64 fRMSD = m_oModel.rmsd(m_oTarStr, false, ALL, false);
    return fRMSD;
  } else
    return 0;

}
///////////////////////////////////////////////////////////////////////////////
// Mutation
///////////////////////////////////////////////////////////////////////////////

/**
 * custom mutation (can be changed by derived class, default implementation just calls mutationBGA)
 * \param iInd index of individual
 */
template<class T>
void svt_gacylinder<T>::mutationCustom(int iInd)
{
  Real64 fRand = svt_genrand();
  if (fRand < 0.30)
    svt_ga<T>::mutationCauchy(iInd);
  else if (fRand < 0.60)
    svt_ga<T>::mutationMultiCauchy(iInd);
  else if (fRand < 0.90)
    svt_ga<T>::mutationRandom(iInd);
  else
    mutationAllCauchy(iInd);
};

/**
 * mutation with a cauchy distribution
 * \param iInd index of individual
 */
template<class T>
void svt_gacylinder<T>::mutationCauchy(int iInd, int iRandIndex, Real64 fRatio)
{

  Real64 fNewGene, fIntPart, fRand;

  fRand = svt_ranCauchy(0.0, this->m_fMutationOffset * fRatio);

  fRand = modf(fRand, &fIntPart);
  fNewGene = this->m_oNextPop[iInd].getGene(iRandIndex) + fRand;

  //bring the gene back into the 0 - 1 range
  this->m_oNextPop[iInd].setGene(iRandIndex, fNewGene);

  this-> makeValid(&this->m_oNextPop[iInd]);
};

/**
 * mutation
 */
template<class T>
void svt_gacylinder<T>::mutationAllCauchy(int iInd)
{
  for (unsigned int iIndex = 0; iIndex < (unsigned int)this->m_iGenes; iIndex++)
    mutationCauchy(iInd, iIndex);
}

/**
 * Compute Average Score of all solutions
 */
template<class T>
void svt_gacylinder<T>::computeAvgScoreAll()
{
  unsigned int iCount = 0;
  m_fAvgScoreAll = 0;

  for (unsigned int iIndex = 0; iIndex < m_oCylinders.size(); iIndex++) {
    m_fAvgScoreAll += m_oCylinders[iIndex].getAvgScore() * m_oCylinders[iIndex].size();
    iCount += m_oCylinders[iIndex].size();
  }

  if (iCount > 0)
    m_fAvgScoreAll /= (Real64)iCount;

}
///////////////////////////////////////////////////////////////////////////////
// Output statistics, result files, etc
///////////////////////////////////////////////////////////////////////////////

/**
 * output result
 */
template<class T>
void svt_gacylinder<T>::outputResult(bool bTabuAdded)
{
  char pFname[256], pOut[256];
  if (m_iWriteModelInterval != 0 && this->m_iGenerations % m_iWriteModelInterval == 0) {
    if ((int)this->m_oPop.size() == this->m_iPopSize) {
      svt_ga<T>::sortPopulation();
      int iSize = 1;// this->m_oPop.size();
      for (int i = 0; i < iSize; i++) {
        updateModel(&this->m_oPop[ this->m_oPop.size() - i - 1 ]);

        //sprintf(pOut, "[%02d-%02d-%04d] %8.6f  - ", this->m_iRun, this->m_iThread, this->m_iGenerations, this->m_oPop[ this->m_oPop.size() - i - 1 ].getFitness());
        SVTLBO << pOut;
        this->m_oPop[ this->m_oPop.size() - i - 1 ].printGenes() ;

        if (strlen(m_pPath) != 0 && this->m_iWriteModelInterval != 0 && this->m_iGenerations % this->m_iWriteModelInterval == 0) {
          sprintf(pFname, "%s/Model%02d%02d%02d.pdb", m_pPath, this->m_iRun, this->m_iThread, (int)this->m_oPop.size() - i - 1);
          //m_oModel.writePDB( pFname, true );
        }
      }
    }
  }



  svt_ga<T>::sortPopulation();
  int iSize = this->m_oPop.size();
  svt_point_cloud_pdb<svt_ga_vec> oPdb;
  svt_point_cloud_atom oAtom;
  for (int i = 0; i < iSize; i++) {
    if (this->m_oPop[ this->m_oPop.size() - i - 1 ].getOrigin() == TABU) {
      SVTLBO << "fitness is 0" << endl;
      oPdb.addAtom(oAtom, this->m_oPop[ this->m_oPop.size() - i - 1 ].getCoarsePhenotype()[0]);
    }
  }

  if (strlen(m_pPath) != 0 && oPdb.size() > 0) {
    sprintf(pFname, "%s/ModelDiscard%02d%02d%02d.pdb", m_pPath, this->m_iRun, this->m_iThread, this->m_iGenerations);
    //  oPdb.writePDB( pFname, true );
  }

}

/**
 * update result
 */
template<class T>
void svt_gacylinder<T>::updateResults(unsigned int iNoOfTubes, int iNum)
{
  int iIndex = 0;
  svt_tube oCylinder;

  svt_point_cloud_pdb<svt_ga_vec>  oPdb, oAllCylinder, oAllCylinderAniDecorr, oTPCylinder, oRestCylinder, oPdbHel;
  //iNoOfTubes was not set; then show all
  if (iNoOfTubes == 0)
    iNoOfTubes = m_oCylindersPdb.size();

  Real64 fTurnsOff = 0;
  Real64 fPercentCovered = 0;
  unsigned int iFalsePositive = 0, iFalseNegative = 0, iTruePositive = 0, iTrueNegative = 0, iDetectedMoreThanOnce = 0;
  unsigned int iCountCylinder = 0;
  char pOut[1256];
  vector<svt_gacylinder_ind> oInds;
  if (iNum == 0) {
    sprintf(pOut, "Index  Score  Points  Length (A)  Turns (for alpha-helices)\n");
    SVTLBO << pOut ;
  }
  for (iIndex = 0; iIndex < (int)iNoOfTubes && iIndex < (int)m_oCylindersPdb.size(); iIndex++) {
    if (iNum == 0) {
      sprintf(pOut, "%4d: %8.5f %5d %7.3f    %7.3f", iIndex + 1, m_oCylinders[iIndex].getScores()[2], m_oCylinders[iIndex].size(), m_oCylinders[iIndex].getLength() , m_oCylindersPdb[iIndex].size() / m_fStepsPerTurn);
      SVTLBO << pOut ;
    }
    oPdb = m_oCylindersPdb[iIndex];
    //oInds = m_oCylinders[iIndex].getElements();
    if (oPdb.size() > 0) {
      svt_ga_vec oVec;
      vector<Real64> oVecMinDist;
      vector<int> oVecModel;
      for (unsigned int i = 0; i < oPdb.size(); i++) {
        oVecMinDist.push_back(1e10);
        oVecModel.push_back(-1);

        oVec = oPdb[i];
        Real64 fDist;
        for (unsigned int iAtom = 0; iAtom < m_oTarStrAxes.size(); iAtom++) {
          fDist = oVec.distance(m_oTarStrAxes[iAtom]);
          if (fDist < oVecMinDist[i]) {
            oVecMinDist[i] = fDist;
            oVecModel[i] = m_oTarStrAxes.getAtom(iAtom)->getModel();
          }
        }

        if (oVecMinDist[i] > m_fMaxDist4Map)
          oVecModel[i] = -1;
      }

      sort(oVecModel.rbegin(), oVecModel.rend());
      int i = oVecModel.size() - 1;
      while (i >= 0 && oVecModel[i] == -1) {
        oVecModel.pop_back();
        i--;
      }

      unsigned int iNum1, iMaxNum = 0;
      int iModel = -1;
      for (unsigned int i = 0; i < oVecModel.size(); i++) {
        if ((i == 0) || (i > 0 && oVecModel[i] != oVecModel[i - 1])) {
          iNum1 = 0;
          for (unsigned int j = 0; j < oVecModel.size(); j++)
            if (oVecModel[j] == oVecModel[i])
              iNum1++;

          if (iNum1 > iMaxNum) {
            iMaxNum = iNum1;
            iModel  = oVecModel[i];
          }
        }
      }
      if (iMaxNum < 4) // not enough points correspond to each other
        iModel = -1;

      Real64 fMinDist = iMaxNum;

      if (iModel == -1) { //matching helix not found
        iFalsePositive ++;
        if (iNum == 0)
          cout << endl;
      } else {
        int iCylinderInTar = -1;
        for (int i = 0; i < (int)m_oTarStrAxesInfo.size(); i++) {
          if (m_oTarStrAxesInfo[i][0] == iModel) {
            m_oTarStrAxesInfo[i][2] =  m_oTarStrAxesInfo[i][2] + 1;
            iCylinderInTar = i;
          }
        }

        if (iCylinderInTar != -1) { // found
          fPercentCovered += m_oCylindersPdb[iIndex].size() / m_fStepsPerTurn > m_oTarStrAxesInfo[iCylinderInTar][1] ? 100 : (100.0 * m_oCylindersPdb[iIndex].size() / m_fStepsPerTurn) / m_oTarStrAxesInfo[iCylinderInTar][1] ;
          fTurnsOff +=  abs(m_oTarStrAxesInfo[iCylinderInTar][1] - m_oCylindersPdb[iIndex].size() / m_fStepsPerTurn);

          sprintf(pOut, " --- %8.3f %8.3f %8.3f %8.3f %4.0f %4.2f ", fMinDist, m_oTarStrAxesInfo[iCylinderInTar][0], m_oTarStrAxesInfo[iCylinderInTar][1], (100.0 * m_oCylindersPdb[iIndex].size() / m_fStepsPerTurn) / m_oTarStrAxesInfo[iCylinderInTar][1], m_oTarStrAxesInfo[iCylinderInTar][3], abs(m_oTarStrAxesInfo[iCylinderInTar][1] - m_oCylindersPdb[iIndex].size() / m_fStepsPerTurn));
          cout << pOut << endl;
          iCountCylinder ++;
          for (unsigned int iAtom = 0; iAtom < oPdb.size(); iAtom++) {
            oPdb.getAtom(iAtom)->setModel(iIndex);
            oTPCylinder.addAtom(*oPdb.getAtom(iAtom), oPdb[iAtom]);
          }
        } else if (iNum == 0) cout << endl;
      }
    } else if (iNum == 0) cout << endl;

    for (unsigned int iAtom = 0; iAtom < oPdb.size(); iAtom++) {
      oPdb.getAtom(iAtom)->setModel(iIndex);
      oAllCylinder.addAtom(*oPdb.getAtom(iAtom), oPdb[iAtom]);
      //ani correction : decompress the helices
      if (m_fAni != 1) {
        oPdb[iAtom].z(m_fOrigZWoAni + m_fAni * (oPdb[iAtom].z() - m_fOrigZ));
        oAllCylinderAniDecorr.addAtom(*oPdb.getAtom(iAtom), oPdb[iAtom]);
      }
    }
  }
  fPercentCovered /= (Real64)iCountCylinder; //average by the no of helices
  fTurnsOff /= (Real64)iCountCylinder;

  //compute rates - false positive, false negatives
  for (unsigned int i = 0; i < m_oTarStrAxesInfo.size(); i++) {
    if (m_oTarStrAxesInfo[i][2] == 0) { // cylinder i was not detected
      SVTLBO << "Cylinder " << i << " was not detected !" << endl;
      //iFalseNegative ++;

      bool bFound = false;
      //search again but this time from the center of the target axis
      unsigned int iCyl;
      for (iCyl = 0; iCyl < m_oCylindersPdb.size() && iCyl < iNoOfTubes; iCyl++) {
        oPdb = m_oCylindersPdb[iCyl];
        for (unsigned int iAtom = 0; iAtom < oPdb.size(); iAtom++) {
          if (m_oTarStrAxesCenters[i].distance(oPdb[iAtom]) <  m_fMaxDist4Map) {
            sprintf(pOut, " again---- %3d %6d %4.2f ", iCyl, m_oCylindersPdb[iCyl].size(), m_oCylindersPdb[iCyl].size() / m_fStepsPerTurn);
            cout << pOut << endl;

            bFound = true;

            iAtom = oPdb.size(); // exit with this condition from looking at this helix
            iCyl = m_oCylindersPdb.size(); // completely exit out of the loop; just interested in the first top scoring helix
          }
        }
      }

      if (bFound && iCyl < iNoOfTubes) {
        //iFalseNegative --;
        iTruePositive ++;
      }
    }

    if (m_oTarStrAxesInfo[i][2] >= 1)  // cylinder i was detected two times
      iTruePositive ++;

    if (m_oTarStrAxesInfo[i][2] > 1) { // cylinder i was detected two times
      SVTLBO << "Cylinder " << i << " was detected  two times!" << endl;
      iDetectedMoreThanOnce++;
    }
  }

  //the remaning Tubes
  for (iIndex = iNoOfTubes; iIndex < (int)m_oCylindersPdb.size(); iIndex++) {
    if (iNum == 0) {
      sprintf(pOut, "%4d: %8.5f %5d %7.3f    %7.3f", iIndex + 1, m_oCylinders[iIndex].getScores()[2], m_oCylinders[iIndex].size(), m_oCylinders[iIndex].getLength() , m_oCylindersPdb[iIndex].size() / m_fStepsPerTurn);
      SVTLBO << pOut ;
    }

    oPdb = m_oCylindersPdb[iIndex];
    if (oPdb.size() > 0) {
      svt_ga_vec oVec;
      vector<Real64> oVecMinDist;
      vector<int> oVecModel;
      for (unsigned int i = 0; i < oPdb.size(); i++) {
        oVecMinDist.push_back(1e10);
        oVecModel.push_back(-1);

        oVec = oPdb[i];
        Real64 fDist;
        for (unsigned int iAtom = 0; iAtom < m_oTarStrAxes.size(); iAtom++) {
          fDist = oVec.distance(m_oTarStrAxes[iAtom]);
          if (fDist < oVecMinDist[i]) {
            oVecMinDist[i] = fDist;
            oVecModel[i] = m_oTarStrAxes.getAtom(iAtom)->getModel();
          }
        }

        if (oVecMinDist[i] > m_fMaxDist4Map)
          oVecModel[i] = -1;
      }

      sort(oVecModel.rbegin(), oVecModel.rend());
      int i = oVecModel.size() - 1;
      while (i >= 0 && oVecModel[i] == -1) {
        oVecModel.pop_back();
        i--;
      }

      unsigned int iNum = 0, iMaxNum = 0;
      int iModel = -1;
      for (unsigned int i = 0; i < oVecModel.size(); i++) {
        if ((i == 0) || (i > 0 && oVecModel[i] != oVecModel[i - 1])) {
          iNum = 0;
          for (unsigned int j = 0; j < oVecModel.size(); j++)
            if (oVecModel[j] == oVecModel[i])
              iNum++;

          if (iNum > iMaxNum) {
            iMaxNum = iNum;
            iModel  = oVecModel[i];
          }
        }
      }
      if (iMaxNum < 4) // not enough points correspond to each other
        iModel = -1;

      Real64 fMinDist = iMaxNum;

      if (iModel == -1) { //too far to be considered an cylinder
        iTrueNegative++;
        if (iNum == 0) cout << endl;
      } else {
        int iCylinderInTar = -1;
        for (int i = 0; i < (int)m_oTarStrAxesInfo.size(); i++) {
          if (m_oTarStrAxesInfo[i][0] == iModel) {
            m_oTarStrAxesInfo[i][2] =  m_oTarStrAxesInfo[i][2] + 1;
            iCylinderInTar = i;
          }
        }

        if (iCylinderInTar != -1) { // found
          sprintf(pOut, " --- %8.3f %8.3f %8.3f %8.3f %4.0f %4.2f", fMinDist, m_oTarStrAxesInfo[iCylinderInTar][0], m_oTarStrAxesInfo[iCylinderInTar][1], (100.0 * m_oCylindersPdb[iIndex].size() / m_fStepsPerTurn) / m_oTarStrAxesInfo[iCylinderInTar][1], m_oTarStrAxesInfo[iCylinderInTar][3], abs(m_oTarStrAxesInfo[iCylinderInTar][1] - m_oCylindersPdb[iIndex].size() / m_fStepsPerTurn));
          cout << pOut << endl;
          iFalseNegative++;
        } else if (iNum == 0) cout << endl;
      }
    } else if (iNum == 0) cout << endl;


    for (unsigned int iAtom = 0; iAtom < oPdb.size(); iAtom++) {
      oPdb.getAtom(iAtom)->setModel(iIndex);
      oRestCylinder.addAtom(*oPdb.getAtom(iAtom), oPdb[iAtom]);
    }
  }

  oAllCylinder.calcAtomModels();
  oRestCylinder.calcAtomModels();
  oTPCylinder.calcAtomModels();

  char pFnameHelices[256];
  if (iNum == 0)
    sprintf(pFnameHelices, "Best%dHelices", m_iNoOfCylinder2Detect);
  else
    strcpy(pFnameHelices, "AllHelices");

  if (strlen(m_pPath) != 0) {
    char pFname[256];
    sprintf(pFname, "%s/%s_%02d.pdb", m_pPath, pFnameHelices, this->m_iRun);
    if (m_fAni != 1) {
      if (oAllCylinderAniDecorr.size() > 0)
        oAllCylinderAniDecorr.writePDB(pFname);
    } else {
      if (oAllCylinder.size() > 0)
        oAllCylinder.writePDB(pFname);
    }


  }

  svt_ga_vol oMapDetectedCylinder;
  if (oTPCylinder.size() > 0)
    oMapDetectedCylinder = *oTPCylinder.blur(1.0, 4.0);

  m_fExploredPercent = m_oTar.correlation(m_oExploredVol) * 100 ;

  for (unsigned int iIndex = 0; iIndex < m_oTarStrAxesInfo.size(); iIndex++)
    m_oTarStrAxesInfo[iIndex][2] = 0.0;
};

/**
 * output the best model
 */
template<class T>
void svt_gacylinder<T>::outputBest()
{
  svt_ga<T>::sortPopulation();

  updateFitness(&this->m_oPop[ this->m_oPop.size() - 1]);

  Real64 fRmsd =  m_oModel.rmsd(m_oTarStr, false, ALL,  false);
  char pOut[1024];
  sprintf(pOut,  "Best individual at generation %04d has CC %15.8f - RMSD: %10.8f CC_after_blurring: ", this->m_iGenerations - 1, this->m_oPop[ this->m_oPop.size() - 1].getFitness(), fRmsd);
  SVTLBO << pOut;

  m_oModelVol.convolve1D3D(m_oKernel, false);  // don't normalize
  Real64 fCC = m_oModelVol.correlation(m_oTar, false);

  sprintf(pOut , "%6.5f\n", fCC);
  cout << pOut;
  if (strlen(m_pPath) != 0) {
    char pFname[256];
    sprintf(pFname, "%s/BestModel%02d%02d_%04d.pdb", m_pPath, this->m_iRun, this->m_iThread, this->m_iGenerations - 1);
    //m_oModel.writePDB( pFname );
  }
}

/**
 * print results (to cout)
 */
template<class T>
void svt_gacylinder<T>::printResults()
{
  for (int i = this->m_oPop.size() - 1; i >= 0; i--) {
    printf("[%3i] = %1d %1d %8.3f", i, this->m_oPop[i].getOrigin(), this->m_oPop[i].getAge(), this->m_oPop[i].getProp());
    this->m_oPop[i].printGenes();
  }
}

/**
 * Set the output path path
 * \param pPath pointer to array of char
 */
template<class T>
void svt_gacylinder<T>::setOutputPath(const char *pPath)
{
  strcpy(m_pPath, pPath);
};

/**
 * Get the output path path
 * \param pPath pointer to array of char
 */
template<class T>
const char *svt_gacylinder<T>::getOutputPath()
{
  return m_pPath;
};

/**
 * How many generations should gacylinder wait until it outputs a new model file (0 turns output off)
 * \param iWriteModelInterval number of generations
 */
template<class T>
void svt_gacylinder<T>::setWriteModelInterval(unsigned int iWriteModelInterval)
{
  m_iWriteModelInterval = iWriteModelInterval;
};
/**
 * How many generations should gacylinder wait until it outputs a new model file (0 turns output off)
 * \return number of generations
 */
template<class T>
unsigned int svt_gacylinder<T>::getWriteModelInterval()
{
  return m_iWriteModelInterval;
};

/**
 * output the configuration of the program:w
 *
 */
template<class T>
void svt_gacylinder<T>::writeConfiguration(char *pFnameParam)
{
  svt_ga<T>::writeConfiguration(pFnameParam);

  FILE *pFileParam = fopen(pFnameParam, "a");

  fprintf(pFileParam, "Resolution = %f\n",                   m_fRes);
  fprintf(pFileParam, "VoxelWidth = %f\n",                   m_oTar.getWidth());
  fprintf(pFileParam, "OutputPath = %s\n",                   m_pPath);
  fprintf(pFileParam, "AngularStepSize = %f\n",              getAngularStepSize());
  fprintf(pFileParam, "PsiFrom = %f\n"  ,                    m_fPsiFrom);
  fprintf(pFileParam, "PsiTo = %f\n"    ,                    m_fPsiTo);
  fprintf(pFileParam, "ThetaFrom = %f\n",                    m_fThetaFrom);
  fprintf(pFileParam, "ThetaTo = %f\n"  ,                    m_fThetaTo);
  fprintf(pFileParam, "PhiFrom = %f\n"  ,                    m_fPhiFrom);
  fprintf(pFileParam, "PhiTo = %f\n"    ,                    m_fPhiTo);
  fprintf(pFileParam, "XFrom = %f\n"    ,                    m_fXFrom);
  fprintf(pFileParam, "XTo = %f\n"      ,                    m_fXTo);
  fprintf(pFileParam, "YFrom = %f\n"    ,                    m_fYFrom);
  fprintf(pFileParam, "YTo = %f\n"      ,                    m_fYTo);
  fprintf(pFileParam, "ZFrom = %f\n"    ,                    m_fZFrom);
  fprintf(pFileParam, "ZTo = %f\n"      ,                    m_fZTo);
  fprintf(pFileParam, "WriteModelInterval = %i\n",           getWriteModelInterval());
  fprintf(pFileParam, "NumberOfTraces = %i\n",               m_iNoOfCylinder2Detect);
  if (m_bApplyBlurring2Model)
    fprintf(pFileParam, "ApplyBlurring2Model = true\n");
  else
    fprintf(pFileParam, "ApplyBlurring2Model = false\n");
  fprintf(pFileParam, "TemplateRadius = %f\n",           getTemplateRadius());
  fprintf(pFileParam, "SearchTemplateRadius = %f\n",     getSearchTemplateRadius());
  fprintf(pFileParam, "TemplatePointCount = %i\n",       getTemplatePointCount());
  fprintf(pFileParam, "TemplateRepeats = %i\n",          getTemplateRepeats());
  fprintf(pFileParam, "SearchTemplateRepeats = %i\n",    getSearchTemplateRepeats());
  fprintf(pFileParam, "DistBetweenRepeats = %f\n",       getDistBetweenRepeats());
  fprintf(pFileParam, "CrawlingStepSize = %f\n",         getCrawlingStepSize());
  fprintf(pFileParam, "AcceptMoveRatio = %f\n",          getAcceptMoveRatio() * 100);
  fprintf(pFileParam, "MaxFailedCrawls = %i\n",          getMaxFailedCrawls());

  if (m_bFitHelix)
    fprintf(pFileParam, "FitHelix = true\n");
  else
    fprintf(pFileParam, "FitHelix = false\n");
  fprintf(pFileParam, "Lambda = %f\n",                   getLambda());

  if (m_bOutputTemplates)
    fprintf(pFileParam, "OutputTemplate = true\n");
  else
    fprintf(pFileParam, "OutputTemplate = false\n");


  fclose(pFileParam);
};



///////////////////////////////////////////////////////////////////////////////
// run ga in thread
///////////////////////////////////////////////////////////////////////////////

/**
 * Write the top scoring solutions to the disk
 * \param oPop the population of solutions
 * \param iWriteSolutions how many solutions to write
 */
template<class T>
void svt_gacylinder<T>::writeSolutions(svt_population<T> &oPop, unsigned int iWriteSolutions, char *pFilename)
{
  int iSize = oPop.size();
  if (iSize > 0 && strlen(m_pPath) != 0) {
    unsigned int iEffectiveWroteSolutions = (int)iWriteSolutions < iSize ? iWriteSolutions : iSize;
    SVTLBO << "Output of the " << iEffectiveWroteSolutions << " best individuals:" << endl;

    sort(oPop.rbegin(), oPop.rend());
    oPop.erase(oPop.begin() + iEffectiveWroteSolutions, oPop.end());

    //initPopulation( m_iPopSize, false );

    char pStr[2560];
    for (unsigned int i = 0; i < iEffectiveWroteSolutions; i++) {
      sprintf(pStr, "%s/%s_%02i%03i.pdb", m_pPath, pFilename, this->m_iRun, i + 1);

      updateModel(&oPop[i]);
      updateVolume(&oPop[i]);
      m_oModel.writePDB(pStr);

      if (m_oTarStr.size() > 0)
        sprintf(pStr, "  [%02ld] %s/%s_%02i%02li.pdb - Score: %1.5f CC: %1.6f RMSD: %1.6f \n", oPop.size() - i, m_pPath, pFilename, this->m_iRun, oPop.size() - i, oPop[i].getFitness(), getCorrelation(), getRMSD());
      else
        sprintf(pStr, "  [%02ld] %s/%s_%02i%02li.pdb - Score: %1.5f \n", oPop.size() - i, m_pPath, pFilename, this->m_iRun, oPop.size() - i, oPop[i].getFitness());
      SVTLBO << pStr;
    }
    SVTLBO << endl;


  }
}

#endif



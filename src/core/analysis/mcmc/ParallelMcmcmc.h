#ifndef ParallelMcmcmc_H
#define ParallelMcmcmc_H

#include "Cloneable.h"
#include "Mcmc.h"
#include "Model.h"
#include "Monitor.h"
#include "Move.h"
#include <vector>
#include "../../utils/Options.h" // included for omp.h ... don't quite understand why it doesn't propagate through from main.cpp...

/**
 * @brief Parallel Metropolis-Coupled Markov chain Monte Carlo (MCMCMC) algorithm class.
 * 
 * This file contains the declaration of the Markov chain Monte Carlo (MCMC) algorithm class. 
 * An MCMC object manages the MCMC analysis by setting up the chain, calling the moves, the monitors and etc.
 *
 *
 *
 * @copyright Copyright 2009-
 * @author The RevBayes Development Core Team (Michael Landis & Sebastian Hoehna)
 * @since Version 1.0, 2013-05-20
 *
 */
class ParallelMcmcmc : Cloneable {
    
public:
    ParallelMcmcmc(Model* m, const std::vector<Move*> &moves, const std::vector<Monitor*> &mons, std::string fn = "", std::string sT="random", int ev = 1,
                int nc=4, int np=4, int si=1000, double dt=0.1, double st=1.0, double sh=1.0, bool saveall = false, size_t steppingStones = 0);
    ParallelMcmcmc(const ParallelMcmcmc &m);
    virtual                                            ~ParallelMcmcmc(void);                                                          //!< Virtual destructor
    
    // public methods
    void                                                burnin(int g, int ti);
    ParallelMcmcmc*                                     clone(void) const;
    void                                                printOperatorSummary(void) const;
    void                                                run(int g);
    void                                                readStream(size_t b = 0);
    bool                                                restore(void);
    unsigned int                                        getCurrentGeneration(void);

private:
    void                                                initialize(void);
    void                                                swapChains(void);
    double                                              computeBeta(double d, double s, size_t i);   // incremental temperature schedule

    void												fromStream(std::istream& is, bool keep = true, bool keepCold = false);
    void												toStream(std::ostream& os);
    void                                                monitorSteppingStone(size_t gen);

    size_t                                              numChains;
    size_t                                              numProcesses;
    std::vector<size_t>                                 chainIdxByHeat;
    std::vector<std::vector<size_t> >                   chainsPerProcess;
    std::vector<Mcmc*>                                  chains;
    std::string                                         scheduleType;
    unsigned int                                        currentGeneration;
    unsigned int                                        swapInterval;
    
    bool												saveall;

    std::string											filename;
    std::fstream										stream;
    int 												every;
    size_t												numNodes;

    size_t                                              activeIndex;  // index of coldest chain, i.e. which one samples the posterior
    double                                              delta;        // delta-T, temperature increment for computeBeta
    double                                              sigma;        // scales power in heating schedule
    double                                              startingHeat; // default 1.0
    
    size_t                                              numSteppingStones;
    std::vector<std::vector<double> >                   steppingStones;

};


#endif /* defined(__rb_mlandis__ParallelMcmcmc__) */

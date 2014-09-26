#ifndef AbstractBirthDeathProcess_H
#define AbstractBirthDeathProcess_H

#include "Taxon.h"
#include "TimeTree.h"
#include "TypedDagNode.h"
#include "TypedDistribution.h"

namespace RevBayesCore {
    
    class Clade;
    
    /**
     * @file
     * This file contains the declaration of the random variable class for constant rate birth-death process.
     * This class is derived from the stochastic node and each instance will represent a random variable.
     *
     * @brief Declaration of the constant rate Birth-Death process class.
     *
     * @copyright Copyright 2009-
     * @author The RevBayes Development Core Team (Sebastian Hoehna)
     * @since 2014-01-17, version 1.0
     *
     */
    class AbstractBirthDeathProcess : public TypedDistribution<TimeTree> {
        
    public:
        AbstractBirthDeathProcess(const TypedDagNode<double> *o, const std::string &cdt, 
                                  const std::vector<Taxon> &tn, const std::vector<Clade> &c);        
        
        // pure virtual member functions
        virtual AbstractBirthDeathProcess*                  clone(void) const = 0;                                                                              //!< Create an independent clone
        
                
        // public member functions you may want to overwrite
        double                                              computeLnProbability(void);                                                                         //!< Compute the log-transformed probability of the current value.
        virtual void                                        swapParameter(const DagNode *oldP, const DagNode *newP);                                            //!< Implementation of swaping parameters
        virtual void                                        redrawValue(void);                                              //!< Draw a new random value from the distribution

    protected:  
        // pure virtual helper functions
        virtual double                                      computeLnProbabilityTimes(void) const = 0;                                                                         //!< Compute the log-transformed probability of the current value.
        virtual std::vector<double>*                        simSpeciations(size_t n, double origin) const = 0;                                        //!< Simulate n speciation events.
        virtual double                                      pSurvival(double start, double end) const = 0;                                                      //!< Compute the probability of survival of the process (without incomplete taxon sampling).
        virtual void                                        prepareProbComputation(void);
        
        // helper functions
        void                                                attachTimes(TimeTree *psi, std::vector<TopologyNode *> &tips, size_t index, 
                                                                        const std::vector<double> *times, double T);
        void                                                buildRandomBinaryTree(std::vector<TopologyNode *> &tips);
        std::vector<double>*                                divergenceTimesSinceOrigin(void) const;                                                             //!< Extract the divergence times from the tree.
        int                                                 diversity(double t) const;                                                                          //!< Diversity at time t.
        std::vector<double>*                                getAgesOfInternalNodesFromMostRecentSample(void) const;                                             //!< Get the ages of all internal nodes since the time of the most recent tip age.
        std::vector<double>*                                getAgesOfTipsFromMostRecentSample(void) const;                                                      //!< Get the ages of all tip nodes since the time of the most recent tip age.
        bool                                                matchesConstraints(void);
        void                                                simulateTree(void);

        // members
        std::string                                         condition;                                                                                          //!< The condition of the process (none/survival/#taxa).
        std::vector<Clade>                                  constraints;                                                                                        //!< Topological constrains.
        const TypedDagNode<double>*                         origin;                                                                                             //!< Time since the origin.
        size_t                                              numTaxa;                                                                                            //!< Number of taxa (needed for correct initialization).
        std::vector<Taxon>                                  taxa;                                                                                         //!< Taxon names that will be attached to new simulated trees.
        double                                              logTreeTopologyProb;                                                                                //!< Log-transformed tree topology probability (combinatorial constant).
        
    };
    
}

#endif

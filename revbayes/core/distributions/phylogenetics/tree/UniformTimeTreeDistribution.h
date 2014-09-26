/**
 * @file
 * This file contains the declaration of the uniform time tree distribution.
 *
 * @brief Declaration of the uniform time tree distribution class.
 *
 * @author Fredrik Ronquist
 *
 */

#ifndef UniformTimeTreeDistribution_H
#define UniformTimeTreeDistribution_H

#include "TimeTree.h"
#include "TypedDagNode.h"
#include "TypedDistribution.h"

namespace RevBayesCore {
    
    class Clade;
    
    class UniformTimeTreeDistribution : public TypedDistribution<TimeTree> {
        
    public:
        UniformTimeTreeDistribution(
                                        const TypedDagNode<double>*                 originT,
                                        const std::vector<std::string>&             taxaNames
                                    );                                                                                  //!< Constructor
        UniformTimeTreeDistribution(const UniformTimeTreeDistribution &x);                                              //!< Copy constructor

        virtual                                            ~UniformTimeTreeDistribution(void);                          //!< Virtual destructor
        
        // public member functions
        UniformTimeTreeDistribution*                        clone(void) const;                                          //!< Create an independent clone
        double                                              computeLnProbability(void);                                 //!< Compute ln prob of current value
        void                                                redrawValue(void);                                          //!< Draw a new random value from distribution
        void                                                swapParameter(const DagNode *oldP, const DagNode *newP);    //!< Swap distribution parameters
        
    private:

        // helper functions
        void                                                attachTimes(TimeTree *psi, std::vector<TopologyNode *> &tips, size_t index, const std::vector<double> &times, double T);
        void                                                buildRandomBinaryHistory(std::vector<TopologyNode *> &tips);
        void                                                simulateTree(void);
        
        // members
        const TypedDagNode<double>*                         originTime;
        size_t                                              numTaxa;
        std::vector<std::string>                            taxonNames;
    };
    
}

#endif

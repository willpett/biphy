/**
 * @file
 * This file contains the declaration of the random variable class for constant rate birth-death process.
 * This class is derived from the stochastic node and each instance will represent a random variable.
 *
 * @brief Declaration of the constant rate Birth-Death process class.
 *
 * (c) Copyright 2009- under GPL version 3
 * @date Last modified: $Date:$
 * @author The RevBayes Development Core Team
 * @license GPL version 3
 * @version 1.0
 * @since 2012-06-17, version 1.0
 * @interface TypedDagNode
 *
 * $Id:$
 */

#ifndef UniformConstrainedTimeTreeDistribution_H
#define UniformConstrainedTimeTreeDistribution_H

#include "TimeTree.h"
#include "TypedDagNode.h"
#include "TypedDistribution.h"

namespace RevBayesCore {
    
    class UniformConstrainedTimeTreeDistribution : public TypedDistribution<TimeTree> {
        
    public:
    	UniformConstrainedTimeTreeDistribution(const TypedDagNode<double>* originT, const std::vector<std::string> &tn,
    										const std::vector<Clade> &c = std::vector<Clade>(),
    										const Clade &o = Clade(std::vector<std::string>(),0));
    	UniformConstrainedTimeTreeDistribution(const UniformConstrainedTimeTreeDistribution &n);                                                                                          //!< Copy constructor
        virtual                                            ~UniformConstrainedTimeTreeDistribution(void);                                                                    //!< Virtual destructor
        
        // public member functions
        UniformConstrainedTimeTreeDistribution*             clone(void) const;                                                                                  //!< Create an independent clone
        double                                              computeLnProbability(void);
        void                                                redrawValue(void);
        void                                                swapParameter(const DagNode *oldP, const DagNode *newP);                                            //!< Implementation of swaping parameters
        
    private:
        
        // helper functions
        void                                                attachTimes(TimeTree *psi, std::vector<TopologyNode *> &tips, size_t index, const std::vector<double> &times, double T);
	void                                                buildRandomBinaryTree(std::vector<TopologyNode *> &tips, unsigned int size);
        void                                                rearrangeRandomBinaryTree(std::vector<TopologyNode*> &tips, std::vector<TopologyNode *> &children);
        void                                                simulateTree(void);
        void                                                rearrangeTree(void);
        bool                                                matchesConstraints(void);
        bool                                                hasOutgroup(void);
        
        // members
        const TypedDagNode<double>*                         originTime;
        unsigned int                                        numTaxa;
        std::vector<std::string>                            taxonNames;
        std::vector<Clade>                                  constraints;
        Clade						    outgroup;
    };
    
}

#endif

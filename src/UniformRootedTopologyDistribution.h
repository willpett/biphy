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

#ifndef UniformRootedTopologyDistribution_H
#define UniformRootedTopologyDistribution_H

#include "Topology.h"
#include "TypedDagNode.h"
#include "TypedDistribution.h"

namespace RevBayesCore {
    
    class UniformRootedTopologyDistribution : public TypedDistribution<Topology> {
        
    public:
    	UniformRootedTopologyDistribution(const std::vector<std::string> &tn,
    														const std::vector<Clade> &c = std::vector<Clade>(),
    														const Clade &o = Clade(std::vector<std::string>(),0));
    	UniformRootedTopologyDistribution(const UniformRootedTopologyDistribution &n);                                                                                          //!< Copy constructor
        virtual                                            ~UniformRootedTopologyDistribution(void);                                                                    //!< Virtual destructor
        
        // public member functions
        UniformRootedTopologyDistribution*                  clone(void) const;                                                                                  //!< Create an independent clone
        double                                              computeLnProbability(void);
        void                                                redrawValue(void);
        void                                                swapParameter(const DagNode *oldP, const DagNode *newP);                                            //!< Implementation of swaping parameters
        
    private:
        
        // helper functions
        void                                                buildRandomBinaryTree(std::vector<TopologyNode *> &tips, unsigned int size);
        //void                                                rearrangeRandomBinaryTree(std::vector<TopologyNode*> &tips, std::vector<TopologyNode *> &children);
        void                                                simulateTree(void);
        //void                                                rearrangeTree(void);
        bool                                                matchesConstraints(void);
        bool                                                hasOutgroup(void);
        
        // members
        unsigned int                                        numTaxa;
        std::vector<std::string>                            taxonNames;
        std::vector<Clade>                                  constraints;
        Clade						    					outgroup;
        double                                              logTreeTopologyProb;

        
    };
    
}

#endif

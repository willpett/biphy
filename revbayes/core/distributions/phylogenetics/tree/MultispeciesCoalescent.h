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

#ifndef MultispeciesCoalescent_H
#define MultispeciesCoalescent_H

#include "TimeTree.h"
#include "TypedDagNode.h"
#include "TypedDistribution.h"

namespace RevBayesCore {
    
    class Clade;
    
    class MultispeciesCoalescent : public TypedDistribution<TimeTree> {
        
    public:
        MultispeciesCoalescent(const TypedDagNode<TimeTree> *st, const std::map<std::string, std::string> &g2S);        
        MultispeciesCoalescent(const MultispeciesCoalescent &n);                                                                                          //!< Copy constructor
        virtual                                            ~MultispeciesCoalescent(void);                                                                    //!< Virtual destructor
        
        // public member functions
        MultispeciesCoalescent*                             clone(void) const;                                                                                  //!< Create an independent clone
        double                                              computeLnProbability(void);
        void                                                redrawValue(void);
        void                                                swapParameter(const DagNode *oldP, const DagNode *newP);                                            //!< Implementation of swaping parameters
        void                                                setNes(TypedDagNode<std::vector<double> >* inputNes);
        void                                                setNe(TypedDagNode<double>* inputNe);

    private:
        
        double                                              getNe(size_t index) const;
        
        // helper functions
        void                                                attachTimes(TimeTree *psi, std::vector<TopologyNode *> &tips, size_t index, const std::vector<double> &times);
        void                                                buildRandomBinaryTree(std::vector<TopologyNode *> &tips);
        void                                                simulateTree(void);
        
        // members
        std::map<std::string, std::string>                  gene2species;
        const TypedDagNode<TimeTree>*                       speciesTree;
        const TypedDagNode<std::vector<double> >*           Nes;
        const TypedDagNode<double >*                        Ne;
        size_t                                              numTaxa;
        double                                              logTreeTopologyProb;
        
    };
    
}

#endif

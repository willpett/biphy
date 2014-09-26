/**
 * @file
 * This file contains the declaration of the tree assembly function. 
 * The function assembles the tree from a topology and a set of branch length.
 * This is useful so that we can separate the distribution on branch lengths from
 * the distribution on the topology.
 *
 * @brief Declaration of the tree assembly function function.
 *
 * (c) Copyright 2009- under GPL version 3
 * @date Last modified: $Date$
 * @author The RevBayes Development Core Team
 * @license GPL version 3
 * @version 1.0
 * @since 2013-03-04, version 1.0
 * @interface Function
 *
 * $Id$
 */



#ifndef TreeAssemblyFunction_H
#define TreeAssemblyFunction_H

#include "BranchLengthTree.h"
#include "Topology.h"
#include "TypedDagNode.h"
#include "TypedFunction.h"

#include <vector>

namespace RevBayesCore {
    
    class TreeAssemblyFunction : public TypedFunction<BranchLengthTree> {
        
    public:
        TreeAssemblyFunction(const TypedDagNode<Topology> *t, const TypedDagNode<std::vector<double> > *b);
        TreeAssemblyFunction(const TreeAssemblyFunction &n);                                                                              //!< Copy constructor
        virtual                                            ~TreeAssemblyFunction(void);                                                    //!< Virtual destructor
        
        // public member functions
        TreeAssemblyFunction*                               clone(void) const;                                                              //!< Create an independent clone
        void                                                keep(DagNode* affecter);
        void                                                restore(DagNode *restorer);   
        void                                                touch(DagNode *toucher );
        void                                                update(void);
        
    protected:
        void                                                swapParameterInternal(const DagNode *oldP, const DagNode *newP);                        //!< Implementation of swaping parameters
        
    private:
        
        // members
        const TypedDagNode<Topology>*                       tau;
        const TypedDagNode<std::vector<double> >*           brlen;
        std::set<size_t>                                    touchedNodeIndices;
    };
    
}

#endif

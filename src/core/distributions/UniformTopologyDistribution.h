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

#ifndef UniformTopologyDistribution_H
#define UniformTopologyDistribution_H

#include "Tree.h"
#include "TypedDagNode.h"
#include "TypedDistribution.h"

class UniformTopologyDistribution : public TypedDistribution<Tree> {
    
public:
    UniformTopologyDistribution(const std::vector<std::string> &tn, const Clade &o = Clade(std::vector<std::string>()), bool rooted = false, bool rigid = false, const Tree* ct = NULL);
    UniformTopologyDistribution(const UniformTopologyDistribution &n);                                                                                          //!< Copy constructor
    virtual                                            ~UniformTopologyDistribution(void);                                                                    //!< Virtual destructor
    
    // public member functions
    UniformTopologyDistribution*                        clone(void) const;                                                                                  //!< Create an independent clone
    double                                              computeLnProbability(void);
    void                                                redrawValue(void);
    void                                                swapParameter(const DagNode *oldP, const DagNode *newP);                                            //!< Implementation of swaping parameters
    // virtual methods
	virtual void                    					setValue(Tree *v);
	virtual void                    					setValue(const Tree &v);

	std::vector<bool>                                   getConstrainedNodes() const;

private:
    
    // helper functions
    void                                                simulateClade(std::vector<TopologyNode *> &tips);
    void                                                rearrangeRandomBinaryTree(std::vector<TopologyNode*> &tips, std::vector<TopologyNode *> &children);
    void                                                simulateTree(void);
    void                                                rearrangeTree(void);
    bool                                                matchesConstraints(void);
    bool                                                hasOutgroup(void);
    void                                                resetNodeIndices(void);
    
    // members
    bool                                                rooted;
    bool												rigid;
    unsigned int                                        numTaxa;
    std::vector<std::string>                            taxonNames;
    Clade                                               outgroup;
    std::vector<Clade>                                  constraints;
    std::map<Clade, double>                             constraint_brlens;
    double                                              logTreeTopologyProb;

    
};

#endif

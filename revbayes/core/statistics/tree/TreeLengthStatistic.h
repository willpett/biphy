//
//  TreeLengthStatistic.h
//  revbayes_mlandis
//
//  Created by Michael Landis on 1/10/13.
//  Copyright (c) 2013 Michael Landis. All rights reserved.
//

#ifndef __revbayes_mlandis__TreeLengthStatistic__
#define __revbayes_mlandis__TreeLengthStatistic__

#include "Tree.h"
#include "TypedDagNode.h"
#include "TypedFunction.h"

#include <vector>
#include <string>

class TreeLengthStatistic : public TypedFunction<double> {
    
public:
    TreeLengthStatistic(const TypedDagNode<Tree> *t);                                                                                   //!< Default constructor
    TreeLengthStatistic(const TreeLengthStatistic& t);                                                                                      //!< Copy constructor
    virtual                                    ~TreeLengthStatistic(void);                                                                  //!< Destructor
    
    TreeLengthStatistic&                        operator=(const TreeLengthStatistic& t);
    
    // Basic utility functions
    TreeLengthStatistic*                        clone(void) const;                                                                          //!< Clone object
    void                                        update(void);                                                                               //!< Clone the function
    
protected:
    void                                        swapParameterInternal(const DagNode *oldP, const DagNode *newP);                            //!< Implementation of swaping parameters
    
private:
    // members
    const TypedDagNode<Tree>*               tree;
    
};

//#include "TreeLengthStatistic.h"

using namespace RevBayesCore;


TreeLengthStatistic::TreeLengthStatistic(const TypedDagNode<Tree> *t) : TypedFunction<double>( new double(0.0) ), tree( t ) {
    // add the tree parameter as a parent
    addParameter( tree );
    
    update();
}


TreeLengthStatistic::TreeLengthStatistic(const TreeLengthStatistic &n) : TypedFunction<double>( n ), tree( n.tree ) {
    // no need to add parameters, happens automatically
}


TreeLengthStatistic::~TreeLengthStatistic( void ) {
    // We don't delete the parameters, because they might be used somewhere else too. The model needs to do that!
}



TreeLengthStatistic* TreeLengthStatistic::clone( void ) const {
    return new TreeLengthStatistic( *this );
}


void TreeLengthStatistic::update( void ) {
    
    double treeHeight = tree->getValue().getRoot().getAge();
    
    std::vector<TopologyNode*> nodes = tree->getValue().getNodes();
    std::vector<TopologyNode*>::iterator it = nodes.begin();
    std::vector<TopologyNode*>::iterator it_end = nodes.end();
    double treeLength = 0.0;
    for ( ; it != it_end; it++)
        treeLength += (*it)->getBranchLength();
    
    *value = treeLength;
}


void TreeLengthStatistic::swapParameterInternal(const DagNode *oldP, const DagNode *newP) {
    if (oldP == tree) {
        tree = static_cast<const TypedDagNode<Tree>* >( newP );
    }
}

#endif /* defined(__revbayes_mlandis__TreeLengthStatistic__) */

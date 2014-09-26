/**
 * @file
 * This file contains the implementation of a time-calibrated tree. This version of a tree
 * contains nodes with times.
 *
 * @brief Implementation of a time calibrated tree
 *
 * (c) Copyright 2009- under GPL version 3
 * @date Last modified: $Date: 2012-07-05 16:47:08 +0200 (Thu, 05 Jul 2012) $
 * @author The RevBayes core development team
 * @license GPL version 3
 * @version 1.0
 * @since 2012-07-17, version 1.0
 *
 * $Id: BranchLengthTree.cpp 1651 2012-07-05 14:47:08Z hoehna $
 */

#include "BranchLengthTree.h"
#include "RbException.h"
#include "RbOptions.h"
#include "TreeChangeEventListener.h"
#include "Topology.h"
#include "TopologyNode.h"

#include <cmath>

using namespace RevBayesCore;

// Declarations
std::ostream& operator<<(std::ostream& o, const BranchLengthTree& x);


/* Default constructor */
BranchLengthTree::BranchLengthTree(void) : Tree() {
    
}


/* Copy constructor */
BranchLengthTree::BranchLengthTree(const BranchLengthTree& t) : Tree(t) {
    
    // just copy the branch lengths
    branchLengths = t.branchLengths;
}


/* Destructor */
BranchLengthTree::~BranchLengthTree(void) {
}


BranchLengthTree& BranchLengthTree::operator=(const BranchLengthTree &t) {
    
    if (this != &t) {
        Tree::operator=(t);
        
        // just copy the branch lengths
        branchLengths   = t.branchLengths;
    }
    
    return *this;
}


/* Clone function */
BranchLengthTree* BranchLengthTree::clone(void) const {
    
    return new BranchLengthTree(*this);
}



double BranchLengthTree::getAge(size_t idx) const {
    
    const TopologyNode &n = topology->getNode( idx );
    
    if ( n.isTip() ) {
        return 0;
    } else {
        double max = 0;
        
        for (size_t i = 0; i < n.getNumberOfChildren(); ++i) {
            const TopologyNode &c = n.getChild( i );
            max = fmax(max, c.getAge() + c.getBranchLength());
        }
        return max;
    }
}


double BranchLengthTree::getBranchLength(size_t idx) const {
    return branchLengths[idx];
}



double BranchLengthTree::getTime(size_t idx) const {
    throw RbException("A branch-length tree does not provide node times.");
    return -1;
}


double BranchLengthTree::getTreeLength(void) const {
    double treeLength = 0;
    const std::vector<TopologyNode*> &nodes = topology->getNodes();
    for (size_t i = 0; i < nodes.size(); i++)
        treeLength += nodes[i]->getBranchLength();
    return treeLength;
}


void BranchLengthTree::resizeElementVectors(size_t n) {
    // remove all elements
    branchLengths.clear();
    
    // resize to new dimension
    branchLengths.resize( n );
}


void BranchLengthTree::setBranchLength(size_t idx, double bl) {
    
    // fire a tree change event
    const std::set<TreeChangeEventListener*> &listeners = changeEventHandler.getListeners();
    for (std::set<TreeChangeEventListener*>::iterator it = listeners.begin(); it != listeners.end(); ++it) {
        (*it)->fireTreeChangeEvent( topology->getNode(idx) );
    }
    
    branchLengths[idx] = bl;
}



std::ostream& operator<<(std::ostream& o, const BranchLengthTree& x) {
    o << x.getNewickRepresentation();
    
    return o;
}

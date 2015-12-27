//
//  TreeConverter.cpp
//  RevBayesCore
//
//  Created by Sebastian Hoehna on 7/17/12.
//  Copyright 2012 __MyCompanyName__. All rights reserved.
//

#include "Tree.h"
#include "TreeUtilities.h"

#include <algorithm>
#include <iostream>

void TreeUtilities::constructTimeTreeRecursively(TopologyNode *tn, const TopologyNode &n, std::vector<TopologyNode*> &nodes, std::vector<double> &ages, double depth) {
    
    // set the name
    tn->setName( n.getName() );
    
    // remember the node
    nodes.push_back( tn );
    
    // set the age
    double a = depth - n.getBranchLength();
    if ( a < 1E-10 ) 
    {
        a = 0.0;
    }
    ages.push_back( a );
    
    // create children
    for (size_t i = 0; i < n.getNumberOfChildren(); ++i) {
        const TopologyNode& child = n.getChild( i );
        TopologyNode* newChild = new TopologyNode();
        
        // set parent child relationship
        newChild->setParent( tn );
        tn->addChild( newChild );
        
        // start recursive call
        constructTimeTreeRecursively(newChild, child, nodes, ages, a);
    }
}


void TreeUtilities::rescaleSubtree(TopologyNode *n, double factor) {
    // we only rescale internal nodes
    if ( !n->isTip() ) {
        // rescale the age of the node
        double newAge = n->getAge() * factor;
        n->setAge(newAge);
        
        // assertion that we have binary trees
#ifdef ASSERTIONS_TREE
        if ( n->getNumberOfChildren() != 2 ) {
            throw Exception("NNI is only implemented for binary trees!");
        }
#endif
        
        // rescale both children
        rescaleSubtree( &n->getChild(0), factor);
        rescaleSubtree( &n->getChild(1), factor);
    }
}


void TreeUtilities::rescaleTree(TopologyNode *n, double factor) {
    // rescale the time of the node
    double newAge = n->getAge() * factor;
    n->setAge( newAge);
    
    // recursive call for internal nodes
    if ( !n->isTip() ) {
        
        // assertion that we have binary trees
#ifdef ASSERTIONS_TREE
        if ( n->getNumberOfChildren() != 2 ) {
            throw Exception("NNI is only implemented for binary trees!");
        }
#endif
        
        // rescale both children
        rescaleTree( &n->getChild(0), factor);
        rescaleTree( &n->getChild(1), factor);
    }
}



std::string TreeUtilities::uniqueNewickTopology(const Tree &t) 
{
    return uniqueNewickTopologyRecursive( t.getRoot() );
}


std::string TreeUtilities::uniqueNewickTopologyRecursive(const TopologyNode &n) 
{
    // check whether this is an internal node
    if ( n.isTip() ) 
    {
        return n.getName();
    } 
    else 
    {
        std::string newick = "(";
        std::vector<std::string> children;
        for (size_t i = 0; i < n.getNumberOfChildren(); ++i) 
        {
            children.push_back( uniqueNewickTopologyRecursive(n.getChild( i ) ) );
        }
        sort(children.begin(), children.end());
        for (std::vector<std::string>::iterator it = children.begin(); it != children.end(); ++it) 
        {
            if ( it != children.begin() ) 
            {
                newick += ",";
            }
            newick += *it;
        }
        newick += ")";
        newick += n.getName();
        
        return newick;
    }
    
}

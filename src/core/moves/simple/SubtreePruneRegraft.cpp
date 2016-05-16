//
//  NearestNeighborInterchange.cpp
//  RevBayesCore
//
//  Created by Sebastian Hoehna on 7/12/12.
//  Copyright 2012 __MyCompanyName__. All rights reserved.
//

#include "SubtreePruneRegraft.h"
#include "RandomNumberFactory.h"
#include "RandomNumberGenerator.h"
#include "TopologyNode.h"
#include "TreeUtilities.h"

#include <cmath>

SubtreePruneRegraft::SubtreePruneRegraft(StochasticNode<Tree> *v, double w, bool outgrp) : SimpleMove( v, w), variable( v ), outgroup(outgrp) {
    
}



/* Clone object */
SubtreePruneRegraft* SubtreePruneRegraft::clone( void ) const {
    
    return new SubtreePruneRegraft( *this );
}



const std::string& SubtreePruneRegraft::getMoveName( void ) const {
    static std::string name = "SPR";
    
    return name;
}


bool SubtreePruneRegraft::isDescendant(const TopologyNode &n, const TopologyNode &p) {
    
    if ( n.isRoot() ) {
        return false;
    }
    
    if ( &n == &p ) {
        return true;
    }
    
    return isDescendant(n.getParent(), p);
}


/** Perform the move */
double SubtreePruneRegraft::performSimpleMove( void ) {
    
    // Get random number generator    
    RandomNumberGenerator* rng     = GLOBAL_RNG;
    
    Tree& tau = variable->getValue();
    
    // pick a random node which is not the root and neithor the direct descendant of the root
    TopologyNode* node = NULL;
    do {
        try{
            double u = rng->uniform01();
            size_t index = std::floor(tau.getNumberOfNodes() * u);
            node = &tau.getNode(index);
        }
        catch(...)
        {
            continue;
        }
    } while ( node == NULL || node->isRoot() || node->getParent().isRoot() || (node->getParent().getParent().isRoot() && outgroup));
    
    // now we store all necessary values
    storedChoosenNode   = node;
    TopologyNode &parent = node->getParent();
    TopologyNode &grandparent = parent.getParent();
    storedBrother = &parent.getChild( 0 );
    
    // check if we got the correct child
    if ( node == storedBrother ) {
        storedBrother = &parent.getChild( 1 );
    }
    
    // pick a random new parent node
    TopologyNode* newBrother = NULL;
    do {
        try
        {
            double u = rng->uniform01();
            size_t index = std::floor(tau.getNumberOfNodes() * u);
            newBrother = &tau.getNode(index);
        }
        catch(...)
        {
            continue;
        }
    } while ( newBrother == NULL ||
            newBrother->isRoot() ||
    		newBrother == &parent ||
			newBrother == storedBrother ||
			newBrother == storedChoosenNode ||
			&(newBrother->getParent()) == storedChoosenNode ||
			(newBrother->getParent().isRoot() && outgroup)
			//|| isDescendant(*newBrother,*node)
		);
    
    TopologyNode &newGrandparent = newBrother->getParent();

    // if the regrafting subtree contains the root
    if(!isDescendant(*newBrother,*node))
    {
    	prunedroot = false;
		// now prune
		grandparent.removeChild( &parent );
		parent.removeChild( storedBrother );
		grandparent.addChild( storedBrother );
		storedBrother->setParent( &grandparent );

		// re-attach
		newGrandparent.removeChild( newBrother );
		parent.addChild( newBrother );
		newGrandparent.addChild( &parent );
		parent.setParent( &newGrandparent );
		newBrother->setParent( &parent );
    }
    // if the pruned subtree contains the root
    else
    {
    	prunedroot = true;
    	// prune
		std::vector<TopologyNode*> storedchildren = storedChoosenNode->getChildren();
		std::vector<size_t> storedindices;
		std::vector<double> storedbrlens;
		for(size_t i = 0 ; i < storedchildren.size(); i++)
		{
			storedindices.push_back(storedchildren[i]->getIndex());
			storedbrlens.push_back(storedchildren[i]->getBranchLength());
			storedChoosenNode->removeChild(storedchildren[i]);
			storedchildren[i]->setParent(NULL);
		}
		TopologyNode* newParent = &(newBrother->getParent());
		TopologyNode* regraft = newParent->reverseParentChild();

		for(size_t i = 0 ; i < storedchildren.size(); i++)
		{
			if(regraft != storedchildren[i])
			{
				regraft->addChild(storedchildren[i]);
				storedchildren[i]->setParent(regraft);

				storedBrother = storedchildren[i];
			}
			else
			{
				newParent->setIndex(storedindices[i]);
				newParent->setBranchLength(storedbrlens[i]);
			}
		}

		newParent->removeChild(newBrother);

		// re-attach
		storedChoosenNode->addChild(newBrother);
		newBrother->setParent(storedChoosenNode);
		storedChoosenNode->addChild(newParent);
		newParent->setParent(storedChoosenNode);

		tau.orderNodesByIndex();
    }
    
    return 0.0;
}


void SubtreePruneRegraft::rejectSimpleMove( void ) {
    
    // undo the proposal
	if(!prunedroot)
	{
		TopologyNode &parent = storedChoosenNode->getParent();
		TopologyNode &grandparent = parent.getParent();
		TopologyNode* oldBrother = &parent.getChild( 0 );
		TopologyNode &newGrandparent = storedBrother->getParent();

		// check if we got the correct child
		if ( storedChoosenNode == oldBrother ) {
			oldBrother = &parent.getChild( 1 );
		}

		// now prune
		grandparent.removeChild( &parent );
		parent.removeChild( oldBrother );
		grandparent.addChild( oldBrother );
		oldBrother->setParent( &grandparent );

		// re-attach
		newGrandparent.removeChild( storedBrother );
		parent.addChild( storedBrother );
		newGrandparent.addChild( &parent );
		parent.setParent( &newGrandparent );
		storedBrother->setParent( &parent );
	}
	else
	{
		// prune
		std::vector<TopologyNode*> storedchildren = storedChoosenNode->getChildren();
		std::vector<size_t> storedindices;
		std::vector<double> storedbrlens;
		for(size_t i = 0 ; i < storedchildren.size(); i++)
		{
			storedindices.push_back(storedchildren[i]->getIndex());
			storedbrlens.push_back(storedchildren[i]->getBranchLength());
			storedChoosenNode->removeChild(storedchildren[i]);
			storedchildren[i]->setParent(NULL);
		}

		TopologyNode* newParent = &(storedBrother->getParent());
		TopologyNode* oldBrother = newParent->reverseParentChild();

		for(size_t i = 0 ; i < storedchildren.size(); i++)
		{
			if(oldBrother != storedchildren[i])
			{
				oldBrother->addChild(storedchildren[i]);
				storedchildren[i]->setParent(oldBrother);
			}
			else
			{
				newParent->setIndex(storedindices[i]);
				newParent->setBranchLength(storedbrlens[i]);
			}
		}

		newParent->removeChild(storedBrother);

		// re-attach
		storedChoosenNode->addChild(storedBrother);
		storedBrother->setParent(storedChoosenNode);
		storedChoosenNode->addChild(newParent);
		newParent->setParent(storedChoosenNode);

		variable->getValue().orderNodesByIndex();
	}
    
}


void SubtreePruneRegraft::swapNode(DagNode *oldN, DagNode *newN) {
    // call the parent method
    SimpleMove::swapNode(oldN, newN);
    
    variable = static_cast<StochasticNode<Tree>* >(newN) ;
}


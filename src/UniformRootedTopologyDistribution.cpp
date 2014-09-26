
#include "UniformRootedTopologyDistribution.h"
#include "RandomNumberFactory.h"
#include "RandomNumberGenerator.h"
#include "RbConstants.h"
#include "Topology.h"
#include "TopologyNode.h"

#include <algorithm>
#include <cmath>

using namespace RevBayesCore;

UniformRootedTopologyDistribution::UniformRootedTopologyDistribution(const std::vector<std::string> &tn, const std::vector<Clade> &c, const Clade &o) : TypedDistribution<Topology>( new Topology() ),
		numTaxa( tn.size() ), taxonNames( tn ), constraints( c ), outgroup( o ) {
    
    double lnFact = 0.0;
    for (size_t i = 2; i < numTaxa; i++) {
        lnFact += std::log(i);
    }

    logTreeTopologyProb = (numTaxa - 1) * RbConstants::LN2 - 2.0 * lnFact - std::log( numTaxa ) ;
    if(hasOutgroup()){
    	std::vector<std::string> ingroup;
		for (size_t i = 0; i < taxonNames.size(); i++) {
			if(std::find(outgroup.begin(),outgroup.end(),taxonNames[i]) == outgroup.end()){
				ingroup.push_back(taxonNames[i]);
			}
		}
    	constraints.push_back(Clade(ingroup,0));
    }
    simulateTree();
    
}



UniformRootedTopologyDistribution::UniformRootedTopologyDistribution(const UniformRootedTopologyDistribution &v) : TypedDistribution<Topology>( v ), numTaxa( v.numTaxa ), taxonNames( v.taxonNames ), constraints( v.constraints), outgroup( v.outgroup), logTreeTopologyProb( v.logTreeTopologyProb ) {
    //rearrangeTree();
}


UniformRootedTopologyDistribution::~UniformRootedTopologyDistribution() {
    
}


void UniformRootedTopologyDistribution::buildRandomBinaryTree(std::vector<TopologyNode*> &tips, unsigned int size) {
    
    if (tips.size() < size) {
        // Get the rng
        RandomNumberGenerator* rng = GLOBAL_RNG;
        
        // randomly draw one node from the list of tips
        size_t index = static_cast<size_t>( floor(rng->uniform01()*tips.size()) );
        
        // get the node from the list
        TopologyNode* parent = tips.at(index);
        
        // remove the randomly drawn node from the list
        tips.erase(tips.begin()+index);
        
        // add a left child
        TopologyNode* leftChild = new TopologyNode(0);
        parent->addChild(leftChild);
        leftChild->setParent(parent);
        tips.push_back(leftChild);
        
        // add a right child
        TopologyNode* rightChild = new TopologyNode(0);
        parent->addChild(rightChild);
        rightChild->setParent(parent);
        tips.push_back(rightChild);
        
        // recursive call to this function
        buildRandomBinaryTree(tips,size);
    }
}

void UniformRootedTopologyDistribution::rearrangeRandomBinaryTree(std::vector<TopologyNode*> &tips, std::vector<TopologyNode*> &children) {

	if (children.size() > 0) {
		// Get the rng
		RandomNumberGenerator* rng = GLOBAL_RNG;
		size_t index1 = static_cast<size_t>( floor(rng->uniform01()*tips.size()) );
		TopologyNode* parent = tips.at(index1);
		tips.erase(tips.begin()+index1);

		std::vector<TopologyNode*> ch = parent->getChildren();
		for (std::vector<TopologyNode*>::iterator it = ch.begin(); it != ch.end(); ++it){
			parent->removeChild(*it,true);
		}

		// randomly draw one node from the list of tips
		TopologyNode* leftChild;
		TopologyNode* rightChild;
		index1 = static_cast<size_t>( floor(rng->uniform01()*children.size()) );
        leftChild = children.at(index1);
        size_t index2 = static_cast<size_t>( floor(rng->uniform01()*children.size()) );
		rightChild = children.at(index2);

		while((leftChild->isTip() && rightChild->isTip() && tips.size() == 0 && children.size() > 2) || rightChild == leftChild){
			index1 = static_cast<size_t>( floor(rng->uniform01()*children.size()) );
			leftChild = children.at(index1);
			index2 = static_cast<size_t>( floor(rng->uniform01()*children.size()) );
			rightChild = children.at(index2);
		}

        parent->addChild(leftChild);
		parent->addChild(rightChild);
		leftChild->setParent(parent);
		rightChild->setParent(parent);
		children.erase(std::find(children.begin(),children.end(),leftChild));
		children.erase(std::find(children.begin(),children.end(),rightChild));

		if(!leftChild->isTip())
			tips.push_back(leftChild);
		if(!rightChild->isTip())
			tips.push_back(rightChild);
		rearrangeRandomBinaryTree(tips,children);
    }
}


UniformRootedTopologyDistribution* UniformRootedTopologyDistribution::clone( void ) const {
    return new UniformRootedTopologyDistribution( *this );
}


double UniformRootedTopologyDistribution::computeLnProbability( void ) {
    
	if ( !matchesConstraints() )
	{
		return RbConstants::Double::neginf;
	}

    return logTreeTopologyProb;
    
}

bool UniformRootedTopologyDistribution::matchesConstraints( void ) {

    const TopologyNode &root = value->getRoot();

    for (std::vector<Clade>::iterator it = constraints.begin(); it != constraints.end(); ++it)
    {
        if ( !root.containsClade( *it, true ) )
        {
            return false;
        }
    }

    return true;
}

void UniformRootedTopologyDistribution::redrawValue( void ) {
    rearrangeTree();
	//simulateTree();
}

bool UniformRootedTopologyDistribution::hasOutgroup( void ) {
    return outgroup.size() > 0;
}

void UniformRootedTopologyDistribution::simulateTree( void ) {

	RandomNumberGenerator* rng = GLOBAL_RNG;
    TopologyNode* root = new TopologyNode();
    std::vector<TopologyNode* > nodes;
//    nodes.push_back(root);


    // add a left child
    TopologyNode* leftChild = new TopologyNode(0);
    root->addChild(leftChild);
    leftChild->setParent(root);
    nodes.push_back(leftChild);

    std::vector<std::string> placed;

    // build the outgroup
	if(hasOutgroup()){
		buildRandomBinaryTree(nodes,outgroup.size());
		for (size_t i=0; i<outgroup.size(); i++) {
			size_t index = size_t( floor(rng->uniform01() * nodes.size()) );

			// get the node from the list
			TopologyNode* node = nodes.at(index);

			// remove the randomly drawn node from the list
			nodes.erase(nodes.begin()+index);

			// set name
			std::string name = outgroup.getTaxonName(i);
			node->setName(name);
			placed.push_back(name);
		}
	}
    // add a right child
    TopologyNode* rightChild = new TopologyNode(0);
    root->addChild(rightChild);
    rightChild->setParent(root);
    nodes.push_back(rightChild);

    size_t numTips = numTaxa-outgroup.size();
    if(constraints.size() > hasOutgroup()){
    	for (size_t i=0; i<constraints.size()-hasOutgroup(); i++) {
    		numTips -= constraints[i].size();
    		numTips++;
    	}
    }
    buildRandomBinaryTree(nodes,numTips);

	for (size_t i=0; i<constraints.size()-hasOutgroup(); i++) {
		size_t index = size_t( floor(rng->uniform01() * nodes.size()) );

		// get the node from the list
		TopologyNode* node = nodes.at(index);

		// remove the randomly drawn node from the list
		nodes.erase(nodes.begin()+index);

		// add node to empty vector
		std::vector<TopologyNode* > tips(1,node);
		// build binary subtree for this constraint
		buildRandomBinaryTree(tips,constraints[i].size());

		//label the tips of this subtree
		for (size_t j=0; j<numTaxa; j++) {
			if(std::find(constraints[i].begin(),constraints[i].end(),taxonNames[j]) != constraints[i].end()){
				size_t index = size_t( floor(rng->uniform01() * tips.size()) );

				// get the node from the list
				TopologyNode* node = tips.at(index);

				// remove the randomly drawn node from the list
				tips.erase(tips.begin()+index);

				// set name
				std::string& name = taxonNames[j];
				node->setName(name);
				placed.push_back(name);
			}
		}
	}
	//add the rest of the taxa
    for (size_t i=0; i<numTaxa; i++) {
    	if(std::find(placed.begin(),placed.end(),taxonNames[i]) == placed.end()){
			size_t index = size_t( floor(rng->uniform01() * nodes.size()) );

			// get the node from the list
			TopologyNode* node = nodes.at(index);

			// remove the randomly drawn node from the list
			nodes.erase(nodes.begin()+index);

			// set name
			std::string& name = taxonNames[i];
			node->setName(name);
    	}
    }

    Topology * tau = new Topology();
    tau->setRoot(root);
    tau->setRooted(true);

    value = tau;
}

void UniformRootedTopologyDistribution::rearrangeTree( void ) {

	RandomNumberGenerator* rng = GLOBAL_RNG;

	TopologyNode* root = &(value->getRoot());
    std::vector<TopologyNode* > nodes = value->getNodes();

    std::vector<TopologyNode* > leaves;
    std::vector<TopologyNode* > internal;

    leaves.insert(leaves.begin(),nodes.begin(),nodes.begin()+numTaxa);
    internal.insert(internal.begin(),nodes.begin()+numTaxa,nodes.end()-1);

    if(hasOutgroup()){
    	std::vector<TopologyNode* > outgrp;
    	for (size_t i=0; i<outgroup.size(); i++) {
    		for (std::vector<TopologyNode*>::iterator it = leaves.begin(); it != leaves.end(); ++it){
				if((*it)->getName() == outgroup.getTaxonName(i)){
					outgrp.push_back((*it));
					leaves.erase(std::find(leaves.begin(),leaves.end(),(*it)));
					break;
				}
			}
		}
    	if(outgrp.size() == 1){
    		std::vector<TopologyNode*> ch = root->getChildren();
			for (std::vector<TopologyNode*>::iterator it = ch.begin(); it != ch.end(); ++it){
				root->removeChild(*it,true);
			}
			root->addChild(outgrp[0]);
			outgrp[0]->setParent(root);
    	}else{
    		size_t index = static_cast<size_t>( floor(rng->uniform01()*internal.size()) );
    		TopologyNode * node = internal.at(index);
    		internal.erase(internal.begin()+index);

    		root->addChild(node);
    		node->setParent(root);

    		std::vector<TopologyNode* > t;
    		std::vector<TopologyNode* > ch;
    		for (size_t i=0; i<outgrp.size()-1; i++) {
    			size_t index = static_cast<size_t>( floor(rng->uniform01()*internal.size()) );
    			TopologyNode * node = internal.at(index);
    			internal.erase(internal.begin()+index);
    			ch.push_back(node);
    		}
    		ch.insert(ch.begin(),outgrp.begin(),outgrp.end());

			t.push_back(node);
			rearrangeRandomBinaryTree(t,ch);
    	}
    	size_t index = static_cast<size_t>( floor(rng->uniform01()*internal.size()) );
		TopologyNode * node = internal.at(index);
		internal.erase(internal.begin()+index);

		root->addChild(node);
		node->setParent(root);
		root = node;
    }
    internal.insert(internal.begin(),leaves.begin(),leaves.end());

    std::vector<TopologyNode* > tips;
    tips.push_back(root);
	rearrangeRandomBinaryTree(tips,internal);
}

void UniformRootedTopologyDistribution::swapParameter(const DagNode *oldP, const DagNode *newP) {
    
}

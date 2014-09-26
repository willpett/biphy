
#include "UniformConstrainedTimeTreeDistribution.h"
#include "RandomNumberFactory.h"
#include "RandomNumberGenerator.h"
#include "RbMathCombinatorialFunctions.h"
#include "RbConstants.h"
#include "Topology.h"
#include "TopologyNode.h"

#include <algorithm>
#include <cmath>

using namespace RevBayesCore;

UniformConstrainedTimeTreeDistribution::UniformConstrainedTimeTreeDistribution(	const TypedDagNode<double>* originT,
        									const std::vector<std::string> &tn,
										const std::vector<Clade> &c,
										const Clade &o) 
	: TypedDistribution<TimeTree>( new TimeTree() ),
        originTime( originT ),
	taxonNames( tn ),
	constraints( c ),
	outgroup( o ) {

	addParameter( originTime );
    numTaxa = taxonNames.size(); 
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



UniformConstrainedTimeTreeDistribution::UniformConstrainedTimeTreeDistribution(const UniformConstrainedTimeTreeDistribution &v) : 
	TypedDistribution<TimeTree>( v ), 
	originTime( v.originTime ),
	numTaxa( v.numTaxa ), 
	taxonNames( v.taxonNames ), 
	constraints( v.constraints), 
	outgroup( v.outgroup) 
{
    //rearrangeTree();
}


UniformConstrainedTimeTreeDistribution::~UniformConstrainedTimeTreeDistribution() {
    
}


void UniformConstrainedTimeTreeDistribution::buildRandomBinaryTree(std::vector<TopologyNode*> &tips, unsigned int size) {
    
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

void UniformConstrainedTimeTreeDistribution::rearrangeRandomBinaryTree(std::vector<TopologyNode*> &tips, std::vector<TopologyNode*> &children) {

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


UniformConstrainedTimeTreeDistribution* UniformConstrainedTimeTreeDistribution::clone( void ) const {
    return new UniformConstrainedTimeTreeDistribution( *this );
}


double UniformConstrainedTimeTreeDistribution::computeLnProbability( void ) {
    
    if ( !matchesConstraints() )
    {
	return RbConstants::Double::neginf;
    }

    // Variable declarations and initialization
    double lnProb = 0.0;
    double originT = originTime->getValue();

    // we need to check as well that all ages are smaller than the origin
    // this can simply be checked if the root age is smaller than the origin
    if ( originT < value->getRoot().getAge() )
    {
        return RbConstants::Double::neginf;
    }

    // Take the uniform draws into account
    lnProb = (numTaxa - 2) * log( 1.0 / originT );

    // Take the ordering effect into account
    lnProb += RbMath::lnFactorial( numTaxa - 2 );

    // We return now; apparently we are not responsible for the topology probability
    return lnProb;
}

bool UniformConstrainedTimeTreeDistribution::matchesConstraints( void ) {

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

void UniformConstrainedTimeTreeDistribution::redrawValue( void ) {
    rearrangeTree();
	//simulateTree();
}

bool UniformConstrainedTimeTreeDistribution::hasOutgroup( void ) {
    return outgroup.size() > 0;
}

void UniformConstrainedTimeTreeDistribution::simulateTree( void ) {
    std::cout << "simulating tree\n";
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
    std::cout << "built topology\n";  
    TimeTree *psi = new TimeTree(); 
    psi->setTopology( tau, true );
    std::cout << "set topology\n";
 
    // Now simulate the speciation times counted from originTime
    std::vector<double> intNodeTimes;
    double              t_or = originTime->getValue();
    intNodeTimes.push_back( 0.0 );  // For time of mrca
    for ( size_t i=0; i<numTaxa-2; i++ )
    {
        double t = rng->uniform01() * t_or;
        intNodeTimes.push_back( t );
    }

    // Sort the speciation times from 0.0 (root node) to the largest value
    std::sort( intNodeTimes.begin(), intNodeTimes.end() );

    // Attach times
    nodes.clear();
    nodes.push_back( root );
    attachTimes(psi, nodes, 0, intNodeTimes, t_or);
    for (size_t i = 0; i < numTaxa; ++i) {
        TopologyNode& node = tau->getTipNode(i);
        psi->setAge( node.getIndex(), 0.0 );
    }

    value = psi;
}

void UniformConstrainedTimeTreeDistribution::rearrangeTree( void ) {

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

    // Now simulate the speciation times counted from originTime
    std::vector<double> intNodeTimes;
    double              t_or = originTime->getValue();
    intNodeTimes.push_back( 0.0 );  // For time of mrca
    for ( size_t i=0; i<numTaxa-2; i++ )
    {
        double t = rng->uniform01() * t_or;
        intNodeTimes.push_back( t );
    }

    // Sort the speciation times from 0.0 (root node) to the largest value
    std::sort( intNodeTimes.begin(), intNodeTimes.end() );

    // Attach times
    nodes.clear();
    nodes.push_back( root );
    attachTimes(value, nodes, 0, intNodeTimes, t_or);
    for (size_t i = 0; i < numTaxa; ++i) {
        TopologyNode& node = value->getTipNode(i);
        value->setAge( node.getIndex(), 0.0 );
    }
}

void UniformConstrainedTimeTreeDistribution::swapParameter(const DagNode *oldP, const DagNode *newP) {
    if (oldP == originTime) {
        originTime = static_cast<const TypedDagNode<double>* >( newP );
    }    
}

void UniformConstrainedTimeTreeDistribution::attachTimes(
                                                TimeTree*                       psi,
                                                std::vector<TopologyNode *>&    nodes,
                                                size_t                          index,
                                                const std::vector<double>&      interiorNodeTimes,
                                                double                          originTime
                                              )
{

    if (index < numTaxa-1)
    {
        // Get the rng
        RandomNumberGenerator* rng = GLOBAL_RNG;

        // Randomly draw one node from the list of nodes
        size_t node_index = static_cast<size_t>( floor(rng->uniform01()*nodes.size()) );

        // Get the node from the list
        TopologyNode* parent = nodes.at(node_index);
        psi->setAge( parent->getIndex(), originTime - interiorNodeTimes[index] );

        // Remove the randomly drawn node from the list
        nodes.erase(nodes.begin()+node_index);

        // Add the left child if an interior node
        TopologyNode* leftChild = &parent->getChild(0);
        if ( !leftChild->isTip() )
        {
            nodes.push_back(leftChild);
        }

        // Add the right child if an interior node
        TopologyNode* rightChild = &parent->getChild(1);
        if ( !rightChild->isTip() )
        {
            nodes.push_back(rightChild);
        }

        // Recursive call to this function
        attachTimes(psi, nodes, index+1, interiorNodeTimes, originTime);
    }
}

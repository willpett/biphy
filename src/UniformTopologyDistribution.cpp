
#include "RandomNumberFactory.h"
#include "RandomNumberGenerator.h"
#include "TopologyNode.h"

#include <algorithm>
#include <cmath>
#include "UniformTopologyDistribution.h"
#include "Constants.h"

UniformTopologyDistribution::UniformTopologyDistribution(const std::vector<std::string> &tn, const Clade &o, bool r) : TypedDistribution<Tree>( new Tree() ),
		numTaxa( tn.size() ), taxonNames( tn ), rooted(r), outgroup( o ) {
    
    double branchLnFact = 0.0;
    double nodeLnFact = 0.0;
    for (size_t i = 2; i < 2*numTaxa-3; i++) {
    	branchLnFact += std::log(i);
    	if(i <= numTaxa - 2)
    		nodeLnFact += std::log(i);
    }

    logTreeTopologyProb = (numTaxa - 2) * Constants::LN2 + nodeLnFact - branchLnFact - (rooted ? log(2*numTaxa-3) : 0);
    
    simulateTree();
    
}



UniformTopologyDistribution::UniformTopologyDistribution(const UniformTopologyDistribution &v) : 
    TypedDistribution<Tree>( v ), 
    numTaxa( v.numTaxa ), 
    taxonNames( v.taxonNames ),
    outgroup( v.outgroup),
    logTreeTopologyProb( v.logTreeTopologyProb ),
    rooted(v.rooted)
{
    //rearrangeTree();
}


UniformTopologyDistribution::~UniformTopologyDistribution() {
    
}

void UniformTopologyDistribution::setValue( Tree *v ) {

	TypedDistribution<Tree>::setValue(v);

	resetNodeIndices();
}

void UniformTopologyDistribution::setValue( const Tree &v ) {

	TypedDistribution<Tree>::setValue(v);

	resetNodeIndices();
}

void UniformTopologyDistribution::buildRandomBinaryTree(std::vector<TopologyNode*> &tips, unsigned int size) {
    
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
        TopologyNode* leftChild = new TopologyNode();
        parent->addChild(leftChild);
        leftChild->setParent(parent);
        tips.push_back(leftChild);
        
        // add a right child
        TopologyNode* rightChild = new TopologyNode();
        parent->addChild(rightChild);
        rightChild->setParent(parent);
        tips.push_back(rightChild);
        
        // recursive call to this function
        buildRandomBinaryTree(tips,size);
    }
}


UniformTopologyDistribution* UniformTopologyDistribution::clone( void ) const {
    return new UniformTopologyDistribution( *this );
}


double UniformTopologyDistribution::computeLnProbability( void ) {
    
	if ( !matchesConstraints() )
		return Constants::Double::neginf;

    return logTreeTopologyProb;
    
}

bool UniformTopologyDistribution::matchesConstraints( void ) {


    const TopologyNode &root = value->getRoot();
    
    bool match = true;
    
    if(hasOutgroup())
    {
        match = false;
            
        const std::vector<TopologyNode*> children = root.getChildren();
        
        for(std::vector<TopologyNode*>::const_iterator it = children.begin(); it != children.end(); it++)
        {
            std::vector<std::string> taxa;
            (*it)->getTaxa(taxa);
            
            bool childFound = true;
            for (std::vector<std::string>::const_iterator y_it = outgroup.begin(); y_it != outgroup.end(); ++y_it)
            {
                bool found = false;
                for (std::vector<std::string>::const_iterator it = taxa.begin(); it != taxa.end(); ++it)
                {
                    if ( (*y_it) == (*it) )
                    {
                        found = true;
                        break;
                    }
                }
                
                if (!found)
                {
                    childFound = false;
                    break;
                }
            }
            
            if(childFound && taxa.size() == outgroup.size())
            {
                match = true;
                break;
            }
        }
    }
    
    return match;
}

void UniformTopologyDistribution::redrawValue( void ) {
    //rearrangeTree();
	simulateTree();
}

bool UniformTopologyDistribution::hasOutgroup( void ) {
    return outgroup.size() > 0;
}

void UniformTopologyDistribution::simulateTree( void ) {

	RandomNumberGenerator* rng = GLOBAL_RNG;
    TopologyNode* root = new TopologyNode();
    
    std::vector<TopologyNode* > nodes;

    // add a left child
    TopologyNode* leftChild = new TopologyNode();
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
    TopologyNode* rightChild = new TopologyNode();
    root->addChild(rightChild);
    rightChild->setParent(root);
    nodes.push_back(rightChild);
    
    // add a middle child
    if(!rooted)
    {
        TopologyNode* middleChild = new TopologyNode();
        root->addChild(middleChild);
        middleChild->setParent(root);
        nodes.push_back(middleChild);
    }

    size_t numTips = numTaxa-outgroup.size();
    buildRandomBinaryTree(nodes,numTips);
    
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
    
    value->setRoot(root);
    value->setRooted(rooted);

    resetNodeIndices();
}

void UniformTopologyDistribution::resetNodeIndices()
{
	TopologyNode& leftChild = value->getRoot().getChild(0);
	TopologyNode& rightChild = value->getRoot().getChild(1);

	size_t leftIndex = leftChild.getIndex();
	size_t rightIndex = rightChild.getIndex();

	if(leftIndex > numTaxa + 1)
	{
		value->getNode(numTaxa).setIndex(leftIndex);
		leftChild.setIndex(numTaxa);
	}
	if(rightIndex > numTaxa + 1)
	{
		value->getNode(numTaxa+1).setIndex(rightIndex);
		rightChild.setIndex(numTaxa+1);
	}

	for (size_t i=0; i<numTaxa; i++)
		value->getTipNodeWithName(taxonNames[i]).setIndex(i);

	value->orderNodesByIndex();
}

void UniformTopologyDistribution::swapParameter(const DagNode *oldP, const DagNode *newP) {
    
}

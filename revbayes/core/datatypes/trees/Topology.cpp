//
//  Topology.cpp
//  RevBayesCore
//
//  Created by Sebastian Hoehna on 3/4/13.
//  Copyright 2013 __MyCompanyName__. All rights reserved.
//

#include "RbException.h"
#include "Topology.h"
#include "TopologyNode.h"
#include "NewickConverter.h"
#include "TreeUtilities.h"

#include <ostream>


using namespace RevBayesCore;


/* Default constructor */
Topology::Topology(void) : 
    root( NULL ),
    binary( true ),
    rooted( false ),
    numTips( 0 ), 
    numNodes( 0 )
{
    
}


/* Copy constructor */
Topology::Topology(const Topology& t) : 
    root( NULL ) 
{
    
    // set the parameters
    binary      = t.binary;
    numTips     = t.numTips;
    numNodes    = t.numNodes;
    rooted      = t.rooted;
    
    // need to perform a deep copy of the BranchLengthTree nodes
    if (t.root != NULL) 
    {
        TopologyNode * newRoot = t.getRoot().clone();

        // set the root. This will also set the nodes vector.
        setRoot(newRoot);
    }
    
}


/* Destructor */
Topology::~Topology(void) 
{

    nodes.clear();
    
    delete root;
}


Topology& Topology::operator=(const Topology &t) 
{
    
    if (this != &t) 
    {

        nodes.clear();
        delete root;
        root = NULL;
        
        binary      = t.binary;
        numTips     = t.numTips;
        numNodes    = t.numNodes;
        rooted      = t.rooted;
        
        TopologyNode* newRoot = t.root->clone();

        // set the root. This will also set the nodes vector.
        setRoot(newRoot);
    }
    
    return *this;
}


/* Clone function */
Topology* Topology::clone(void) const 
{
    
    return new Topology(*this);
}


std::vector<std::string> Topology::getTipNames( void ) const 
{
    std::vector<std::string> names;
    for (size_t i = 0; i < getNumberOfTips(); ++i) 
    {
        const TopologyNode& n = getTipNode( i );
        names.push_back( n.getName() );
    }
    
    return names;
}

/* fill the nodes vector by a phylogenetic traversal recursively starting with this node. 
 * The tips fill the slots 0,...,n-1 followed by the internal nodes and then the root.
 */
void Topology::fillNodesByPhylogeneticTraversal(TopologyNode* node) {
    
    // now call this function recursively for all your children
    for (size_t i=0; i<node->getNumberOfChildren(); i++) {
        fillNodesByPhylogeneticTraversal(&node->getChild(i));
    }
    
    if (node->isTip()) {
        // all the tips go to the beginning
        nodes.insert(nodes.begin(), node);
    } 
    else {
        // this is phylogenetic ordering so the internal nodes come last
        nodes.push_back(node);
    }
}


const std::string& Topology::getNewickRepresentation( void ) const {
    
    return root->computeNewick();
}


std::vector<TopologyNode *> Topology::getNodes( void ) const {
    
    return nodes;
}


TopologyNode& Topology::getNode(size_t idx) {
    
    if ( idx >= nodes.size() ) 
    {
    	std::stringstream ss;
    	ss << "Index "<< idx << " out of bounds in getNode of Topology";
        throw RbException(ss.str());
    }
    
    return *nodes[idx];
}


const TopologyNode& Topology::getNode(size_t idx) const {
    
    if ( idx >= nodes.size() ) 
    {
    	std::stringstream ss;
    	ss << "Index "<< idx << " out of bounds in getNode of Topology";
        throw RbException(ss.str());
    }
    
    return *nodes[idx];
}


/** Calculate the number of interior nodes in the BranchLengthTree by deducing the number of
 tips from number of nodes, and then subtract 1 more if the BranchLengthTree is rooted. */
size_t Topology::getNumberOfInteriorNodes( void ) const {
    
    size_t preliminaryNumIntNodes = getNumberOfNodes() - getNumberOfTips();
    
    if ( isRooted() )
        return preliminaryNumIntNodes - 1;
    else
        return preliminaryNumIntNodes;
}


size_t Topology::getNumberOfNodes( void ) const {
    
    return nodes.size();
}


std::string Topology::getPlainNewickRepresentation( void ) const {
    
    return root->computePlainNewick();
}

size_t Topology::getMrca(size_t index1, size_t index2) const {

	std::vector<size_t> chain1;
	std::vector<size_t> chain2;

	const TopologyNode& node1 = getNode(index1);
	const TopologyNode& node2 = getNode(index2);

	if(node1.isRoot())
		return node1.getIndex();

	if(node2.isRoot())
		return node2.getIndex();

	size_t parentIndex = node1.getParent().getIndex();

	chain1.push_back(parentIndex);
	while(!getNode(parentIndex).isRoot()){
		parentIndex = getNode(parentIndex).getParent().getIndex();
		chain1.insert(chain1.begin(),parentIndex);
	}

	parentIndex = node2.getParent().getIndex();

	chain2.push_back(parentIndex);
	while(!getNode(parentIndex).isRoot()){
		parentIndex = getNode(parentIndex).getParent().getIndex();
		chain2.insert(chain2.begin(),parentIndex);
	}

	size_t mrca = getRoot().getIndex();

	for(size_t index = 0; index < std::min(chain1.size(), chain2.size()); index++){
		if(chain1[index] == chain2[index]){
			mrca = chain1[index];
		}else{
			break;
		}
	}

	return mrca;
}


/** Calculate and return the number of tips on the BranchLengthTree by going through the vector
 of nodes, querying each about its tip status. */
size_t Topology::getNumberOfTips( void ) const {
    
    size_t n = 0;
    for (size_t i=0; i<nodes.size(); i++) {
        if (nodes[i]->isTip() == true)
            n++;
    }
    return n;
}


/** We provide this function to allow a caller to randomly pick one of the interior nodes.
 This version assumes that the root is always the last and the tips the first in the nodes vector. */
const TopologyNode& Topology::getInteriorNode( size_t indx ) const {
    
    // \TODO: Bound checking, maybe draw from downpass array instead
    return *nodes[ indx + getNumberOfTips() ];
}


TopologyNode& Topology::getRoot( void ) {
    
    return *root;
}


const TopologyNode& Topology::getRoot( void ) const {
    
    return *root;
}


/** We provide this function to allow a caller to randomly pick one of the interior nodes.
 This version assumes that the tips are first in the nodes vector. */
TopologyNode& Topology::getTipNode( size_t indx ) {
    
#ifdef ASSERTIONS_BranchLengthTree
    if ( indx >= getNumberOfTips() ) {
        throw RbException("Index out of bounds in getTipNode() of BranchLengthTree!");
    }
    if (!nodes[ indx ]->isTip()) {
        throw RbException("Node at index is not a tip but should have been!");
    }
#endif
    
    // \TODO: Bound checking
    return *nodes[ indx ];
}


/** We provide this function to allow a caller to randomly pick one of the interior nodes.
 This version assumes that the tips are first in the nodes vector. */
const TopologyNode& Topology::getTipNode( size_t indx ) const {
    
    // \TODO: Bound checking
    return *nodes[ indx ];
}


bool Topology::isBinary( void ) const {
    return binary;
}


bool Topology::isRooted( void ) const {
    return rooted;
}


void Topology::setRooted(bool tf) {
    rooted = tf;
}

void Topology::reRoot(size_t index) {
	TopologyNode* node = nodes[index];

	if(node->isRoot())
		return;

	TopologyNode* parent = &node->getParent();

	if(rooted){
		if(parent->isRoot())
			return;

		TopologyNode* childA = &root->getChild(0);
		TopologyNode* childB = &root->getChild(1);

		// remove root's children
		root->removeChild(childA);
		root->removeChild(childB);

		//set the children as roots
		childA->setParent(NULL);
		childB->setParent(NULL);

		//reattach the root and reset parent-child relationships
		parent->removeChild(node);

		root->addChild(node);
		node->setParent(root);

		root->addChild(parent);
		parent->setParent(root);

		TopologyNode* r = setChild(&parent->getParent(), parent);

		//check which subtree contains the root
		//and reattach appropriately
		if(r == childA){
			childA->addChild(childB);
			childB->setParent(childA);
		}else{
			childB->addChild(childA);
			childA->setParent(childB);
		}
	}else{
		std::vector<TopologyNode*> rootChildren = root->getChildren();
		std::vector<TopologyNode*> nodeChildren = node->getChildren();

		for(size_t i = 0; i < nodeChildren.size(); i++){
			node->removeChild(nodeChildren[i]);
			root->addChild(nodeChildren[i]);
			nodeChildren[i]->setParent(root);
		}

		for(size_t i = 0; i < rootChildren.size(); i++){
			if(rootChildren[i] != node){
				root->removeChild(rootChildren[i]);
				node->addChild(rootChildren[i]);
				rootChildren[i]->setParent(node);
			}
		}

		if(parent->isRoot())
			return;

		root->addChild(parent);
		parent->setParent(root);

		parent->removeChild(node);
		node->setParent(NULL);

		setChild(&parent->getParent(), parent);
	}
}

TopologyNode* Topology::setChild(TopologyNode* oldParent, TopologyNode* oldChild) {
	oldParent->removeChild(oldChild);

	TopologyNode* r = oldParent;
	if(!oldParent->isRoot())
		r = setChild(&oldParent->getParent(), oldParent);

	oldChild->addChild(oldParent);
	oldParent->setParent(oldChild);

	return r;
}

void Topology::setRoot( TopologyNode* r) {

    // set the root
    root = r;
    
//    root->setTopology( this );
    
    rooted = (root->getNumberOfChildren() == 2);

    nodes.clear();
    
    // bootstrap all nodes from the root and add the in a pre-order traversal
    // fillNodesByPreorderTraversal(r);
    fillNodesByPhylogeneticTraversal(r);
    std::vector<TopologyNode*> newnodes(nodes.size());
    for (unsigned int i = 0; i < nodes.size(); ++i) {
    	if(nodes[i]->getIndex() == -1 || ((nodes.size() != numNodes) && (numNodes > 0)))
    		nodes[i]->setIndex(i);

    	newnodes[nodes[i]->getIndex()] = nodes[i];
    }
    nodes = newnodes;
    
    numNodes = nodes.size();
}



std::ostream& RevBayesCore::operator<<(std::ostream& o, const Topology& x) {
    o << x.getNewickRepresentation();
    
    return o;
}

std::istream& RevBayesCore::operator>>(std::istream& is, Topology& x) {
    std::string tmp;
    is >> tmp;

	NewickConverter c;
    BranchLengthTree * tree = c.convertFromNewick(tmp);
    x = tree->getTopology();
    return is;
}

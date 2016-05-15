#include "DagNode.h"
#include "Tree.h"
#include "TopologyNode.h"
#include "TypedDagNode.h"
#include "NewickConverter.h"

#include <cmath>
#include "Constants.h"
#include "Exception.h"
#include "Options.h"

/* Default constructor */
Tree::Tree(void) :
    changeEventHandler(),
    root( NULL ),
    binary( true ),
    rooted( false ),
    numTips( 0 ),
    numNodes( 0 )
{
    
}


/* Copy constructor */
Tree::Tree(const Tree& t) : 
    changeEventHandler( ),
    root( NULL ),
    binary( t.binary ),
    rooted( t.rooted ),
    numTips( t.numTips ),
    numNodes( t.numNodes )
{
        
    // need to perform a deep copy of the BranchLengthTree nodes
    if (t.root != NULL)
    {
        TopologyNode * newRoot = t.getRoot().clone();
        
        // set the root. This will also set the nodes vector.
        // do not reorder node indices when copying (WP)
        setRoot(newRoot, false);
    }
    
}


/* Destructor */
Tree::~Tree(void) 
{
    
    nodes.clear();
    
    delete root;
    
}


Tree& Tree::operator=(const Tree &t)
{
    
    if (this != &t) 
    {
        // nothing really to do here, should be done in the derived classes
        // @todo: Find a better solution - Sebastian
        // Problem: If we redraw the tree because the initial states are invalid, 
        // then we somehow need to remember the tree event change listeners.
        // But it is not nice if the tree distribution needs to remember this!!!
//        changeEventHandler = t.changeEventHandler;
        
        
        
        nodes.clear();
        delete root;
        root = NULL;

        binary      = t.binary;
        numTips     = t.numTips;
        numNodes    = t.numNodes;
        rooted      = t.rooted;
        
        TopologyNode* newRoot = t.root->clone();
        
        // set the root. This will also set the nodes vector
        // do not reorder node indices when copying (WP)
        setRoot(newRoot, false);
        
    }
    
    return *this;
}


bool Tree::operator==(const Tree &t) const
{
    
    return getNewickRepresentation() == t.getNewickRepresentation();
}


bool Tree::operator!=(const Tree &t) const
{
    
    return !operator==(t);
}


bool Tree::operator<(const Tree &t) const
{
    
    return getNewickRepresentation() < t.getNewickRepresentation();
}


bool Tree::operator<=(const Tree &t) const
{
    
    return operator<(t) || operator==(t);
}


void Tree::addBranchParameter(std::string const &name, const std::vector<double> &parameters, bool internalOnly)
{
    
    getRoot().addBranchParameters(name,parameters,internalOnly);
    
}


void Tree::addNodeParameter(std::string const &name, const std::vector<double> &parameters, bool internalOnly)
{
    
    getRoot().addNodeParameters(name,parameters,internalOnly);
    
}

void Tree::addNodeParameter(std::string const &name, const std::vector<std::string*> &parameters, bool internalOnly)
{
    
    getRoot().addNodeParameters(name,parameters,internalOnly);
    
}


void Tree::clearParameters( void )
{
    
    clearNodeParameters();
    clearBranchParameters();
    
}

void Tree::clearBranchParameters( void )
{
    
    getRoot().clearBranchParameters();

}


void Tree::clearNodeParameters( void )
{
    
    getRoot().clearNodeParameters();
	
}


/* Clone function */
Tree* Tree::clone(void) const
{

    return new Tree(*this);
}

/* fill the nodes vector by a phylogenetic traversal recursively starting with this node.
 * The tips fill the slots 0,...,n-1 followed by the internal nodes and then the root.
 */
void Tree::fillNodesByPhylogeneticTraversal(TopologyNode* node)
{

    // now call this function recursively for all your children
    for (size_t i=0; i<node->getNumberOfChildren(); i++)
    {
        fillNodesByPhylogeneticTraversal(&node->getChild(i));
    }

    if (node->isTip())
    {
        // all the tips go to the beginning
        nodes.insert(nodes.begin(), node);
    }
    else
    {
        // this is phylogenetic ordering so the internal nodes come last
        nodes.push_back(node);
    }
}


/** We provide this function to allow a caller to randomly pick one of the interior nodes.
 This version assumes that the root is always the last and the tips the first in the nodes vector. */
const TopologyNode& Tree::getInteriorNode( size_t indx ) const
{
    
    // \TODO: Bound checking, maybe draw from downpass array instead
    return *nodes[ indx + getNumberOfTips() ];
}


std::string Tree::getNewickRepresentation() const
{
    
    return root->computeNewick();
}



TopologyNode& Tree::getNode(size_t idx)
{
    
    if ( idx >= nodes.size() )
    {
        throw Exception("Index out of bounds in getNode.");
    }
    
    return *nodes[idx];
}


const TopologyNode& Tree::getNode(size_t idx) const
{
    
    if ( idx >= nodes.size() )
    {
        throw Exception("Index out of bounds in getNode.");
    }
    
    return *nodes[idx];
}


const std::vector<TopologyNode*>& Tree::getNodes(void) const
{
    
    return nodes;
}




/** 
 * Calculate the number of interior nodes in the BranchLengthTree by deducing the number of
 * tips from number of nodes, and then subtract 1 more if the BranchLengthTree is rooted. 
 */
size_t Tree::getNumberOfInteriorNodes( void ) const
{
    
    size_t preliminaryNumIntNodes = getNumberOfNodes() - getNumberOfTips();
    
    if ( isRooted() )
    {
        return preliminaryNumIntNodes - 1;
    }
    else
    {
        return preliminaryNumIntNodes;
    }
    
}


size_t Tree::getNumberOfNodes(void) const
{
    
    return numNodes;
}


/** 
 * return the number of tips.
 */
size_t Tree::getNumberOfTips( void ) const
{
    
    return numTips;
}


TopologyNode& Tree::getRoot(void)
{
    return *root;
}


const TopologyNode& Tree::getRoot(void) const
{
    return *root;
}


/**
 * Get all the species names for this topology.
 * This might include duplicates.
 */
std::vector<std::string> Tree::getSpeciesNames() const
{
    std::vector< std::string > snames;
    for (size_t i = 0; i < getNumberOfTips(); ++i)
    {
        const TopologyNode& n = getTipNode( i );
        snames.push_back( n.getTaxon() );
    }
    
    return snames;
}

std::vector<std::string> Tree::getTaxa() const
{
    std::vector< std::string > taxa;
    for (size_t i = 0; i < getNumberOfTips(); ++i)
    {
        const TopologyNode& n = getTipNode( i );
        taxa.push_back( n.getTaxon() );
    }
    
    return taxa;
}


/**
 * Get the tip index for this name.
 */
size_t Tree::getTipIndex( const std::string &name ) const
{
    for (size_t i = 0; i < getNumberOfTips(); ++i)
    {
        const TopologyNode& n = getTipNode( i );
        if ( name == n.getName() )
        {
            return n.getIndex();
        }
    }
    
    // if name not found
    throw Exception("Could not find tip node with name '" + name + "' in tree." );
}


std::vector<std::string> Tree::getTipNames() const
{
    
    std::vector<std::string> names;
    for (size_t i = 0; i < getNumberOfTips(); ++i)
    {
        const TopologyNode& n = getTipNode( i );
        names.push_back( n.getName() );
    }
    
    return names;
}


/**
 * We provide this function to allow a caller to randomly pick one of the tip nodes.
 * This version assumes that the tips are first in the nodes vector.
 */
TopologyNode& Tree::getTipNode( size_t index )
{
    
//    
//    if ( index >= getNumberOfTips() )
//    {
//        throw Exception("Index out of bounds in getTipNode()!");
//    }
//    if (!nodes[ index ]->isTip())
//    {
//        throw Exception("Node at index is not a tip but should have been!");
//    }
    
    return *nodes[ index ];
}


const TopologyNode& Tree::getTipNode(size_t index) const
{
    
    //
    //    if ( index >= getNumberOfTips() )
    //    {
    //        throw Exception("Index out of bounds in getTipNode()!");
    //    }
    //    if (!nodes[ index ]->isTip())
    //    {
    //        throw Exception("Node at index is not a tip but should have been!");
    //    }
    
    return *nodes[index];
}


/**
 * Get the tip node with the given name.
 * The name should correspond to the taxon name, not the species name.
 * This will throw an error if the name doesn't exist.
 */
TopologyNode& Tree::getTipNodeWithName( const std::string &n )
{
    // get the index of this name
    size_t index = getTipIndex( n );
    
    return *nodes[ index ];
}



/**
 * Get the tip node with the given name.
 * The name should correspond to the taxon name, not the species name.
 * This will throw an error if the name doesn't exist.
 */
const TopologyNode& Tree::getTipNodeWithName( const std::string &n ) const
{
    // get the index of this name
    size_t index = getTipIndex( n );
    
    return *nodes[ index ];
}


/**
 * Get all the tip nodes with this species name.
 * If there is none, then we return an empty vector.
 */
std::vector<TopologyNode*> Tree::getTipNodesWithSpeciesName( const std::string &name )
{
    // create the vector of the tip nodes with this species name
    std::vector<TopologyNode*> tipNodes;
    
    // loop over all tips
    for (size_t i = 0; i < getNumberOfTips(); ++i)
    {
        // get the i-th tip
        TopologyNode& n = getTipNode( i );
        
        // test if the species name matches
        if ( name == n.getSpeciesName() )
        {
            // add this tip node to our list
            tipNodes.push_back( &n );
        }
    }
    
    // return the vector
    return tipNodes;
}


double Tree::getTmrca(const TopologyNode &n)
{
    
    return root->getTmrca( n );
}


TreeChangeEventHandler& Tree::getTreeChangeEventHandler( void ) const
{
    
    return changeEventHandler;
}


//bool Tree::hasSameTopology(const Tree &t) const
//{
//    
//    return hasSameTopology( t.getTopology() );
//}
//
//
//bool Tree::hasSameTopology(const Topology &t) const
//{
//    
//    return topology->getPlainNewickRepresentation() == t.getPlainNewickRepresentation();
//}


bool Tree::isBinary(void) const 
{
    
    return binary;
}


bool Tree::isBroken( void ) const
{
    
    for (size_t i = 0; i < getNumberOfInteriorNodes(); ++i)
    {
        
        const TopologyNode &n = getInteriorNode( i );
        double age = n.getAge();
        
        for (size_t j = 0; j < n.getNumberOfChildren(); ++j)
        {
            const TopologyNode &child = n.getChild( j );
            
            double est_age = child.getAge() + child.getBranchLength();
            
            if ( std::fabs(age-est_age) > 1E-4 )
            {
                return true;
            }
            
        }
        
    }
    
    for (size_t i = 0; i < getNumberOfNodes(); ++i)
    {
        
        const TopologyNode &n = getNode( i );
        
        if ( n.isRoot() == false )
        {
            double my_age = n.getAge();
            double my_parents_age = n.getParent().getAge();
            
            if ( std::fabs( my_parents_age - my_age - n.getBranchLength() ) > 1E-4 )
            {
                return true;
            }
        
        }
        
    }
    
    
    return false;
}


bool Tree::isRooted(void) const 
{
    
    return rooted;
}


bool Tree::isUltrametric( void ) const
{
    
    double tip_age = getTipNode( 0 ).getAge();
    for (size_t i = 1; i < getNumberOfTips(); ++i)
    {
        
        if ( std::fabs(tip_age-getTipNode(i).getAge()) > 1E-4 )
        {
            return false;
        }
        
    }
    
    return true;
}


void Tree::makeInternalNodesBifurcating(void)
{
    
    getRoot().makeBifurcating();
    
    // we need to reset the root so that the vector of nodes get filled again with the new number of nodes
    setRoot( &getRoot() );

}


// method to order nodes by their existing index
// used when reading in tree with existing node indexes we need to keep
void Tree::orderNodesByIndex()
{

    std::vector<TopologyNode*> nodes_copy = std::vector<TopologyNode*>(numNodes);
    for (int i = 0; i < numNodes; i++)
    {
        for (int j = 0; j < numNodes; j++)
        {
            if (i == nodes[j]->getIndex())
            {
                nodes_copy[i] = nodes[j];
            }
        }
    }
    nodes = nodes_copy;

}


void Tree::reroot(const std::string &outgroup)
{
    std::vector<std::string> tipnames = getTipNames();
    size_t outgroupIndex = tipnames.size();
    for (size_t i=0; i<tipnames.size(); ++i)
    {
        if ( tipnames[i] == outgroup )
        {
            outgroupIndex = i;
            break;
        }
    }

    if ( outgroupIndex == tipnames.size() )
    {
        throw Exception("Cannot reroot the tree because we could not find an outgroup with name '" + outgroup + "'.");
    }

    // reset parent/child relationships
	TopologyNode& outgroupNode = getTipNode( outgroupIndex );
	outgroupNode.getParent().reverseParentChild( );
    outgroupNode.getParent().setParent( NULL );

	// set the new root
	setRoot( &outgroupNode.getParent() );

}

void Tree::reroot(TopologyNode &n)
{
	// reset parent/child relationships
	n.getParent().reverseParentChild();
    n.getParent().setParent( NULL );

	// set the new root
	setRoot( &n.getParent() );

}


void Tree::setRooted(bool tf)
{
    rooted = tf;
}


void Tree::setRoot( TopologyNode* r, bool resetIndex )
{

    // delete the old root
    if ( r != root )
    {
        delete root;
    }
    
    // set the root
    root = r;

    nodes.clear();

    // bootstrap all nodes from the root and add the in a pre-order traversal
    fillNodesByPhylogeneticTraversal(r);

    binary = true;

    std::vector<bool> used_indices(nodes.size(), false);
    for (unsigned int i = 0; i < nodes.size(); ++i)
    {
        if(nodes[i]->getIndex() != -1)
        {
            used_indices[nodes[i]->getIndex()] = true;
        }
    }

    std::vector<size_t> unused_indices;
    for (unsigned int i = 0; i < nodes.size(); ++i)
    {
        if(!used_indices[i])
        {
            unused_indices.push_back(i);
        }
    }

    for (unsigned int i = 0; i < nodes.size(); ++i)
    {
        if(nodes[i]->getIndex() == -1)
        {
            nodes[i]->setIndex(unused_indices.front());
            unused_indices.erase(unused_indices.begin());
        }

        if(nodes[i]->getNumberOfChildren() != 2 && !nodes[i]->isRoot() && !nodes[i]->isTip())
        {
            binary = false;
        }
    }

    numNodes = nodes.size();

    orderNodesByIndex();
    // count the number of tips
    numTips = 0;
    for (size_t i = 0; i < numNodes; ++i)
    {
        numTips += ( nodes[i]->isTip() ? 1 : 0);
    }
    
    
    root->setTree( this );

}


std::ostream& operator<<(std::ostream& o, const Tree& x)
{
    
    o << x.getNewickRepresentation();
        
    return o;
}

std::istream& operator>>(std::istream& o, Tree& x)
{
    std::string tmp;
    o >> tmp;
    
    NewickConverter c;
    x = *(c.convertFromNewick(tmp));
        
    return o;
}

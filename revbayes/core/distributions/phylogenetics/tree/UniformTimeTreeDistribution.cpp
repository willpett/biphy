/**
 * @file
 * This file contains the implementation of the uniform time tree distribution.
 *
 * @brief Implementation of the uniform time tree distribution class.
 *
 * @author Fredrik Ronquist
 *
 */

#include "Clade.h"
#include "RandomNumberFactory.h"
#include "RandomNumberGenerator.h"
#include "RbConstants.h"
#include "RbMathCombinatorialFunctions.h"
#include "TopologyNode.h"
#include "Topology.h"
#include "UniformTimeTreeDistribution.h"

#include <algorithm>
#include <cmath>

using namespace RevBayesCore;

UniformTimeTreeDistribution::UniformTimeTreeDistribution(
                                                            const TypedDagNode<double>*         originT,
                                                            const std::vector<std::string>&     taxaNames
                                                         )
    :   TypedDistribution<TimeTree>( new TimeTree() ),
        originTime( originT ),
        taxonNames( taxaNames )
{
    
    addParameter( originTime );
    
    numTaxa = (int)( taxonNames.size() );
    
    simulateTree();
    
}


UniformTimeTreeDistribution::UniformTimeTreeDistribution(const UniformTimeTreeDistribution &x)
    :   TypedDistribution<TimeTree>( x ),
        originTime( x.originTime ),
        numTaxa( x.numTaxa ),
        taxonNames( x.taxonNames )
{
    // parameters are automatically copied
}


UniformTimeTreeDistribution::~UniformTimeTreeDistribution()
{
    
}


/**
 * Recursive call to attach ordered interior node times to the time tree psi. Call it initially with the
 * root of the tree.
 */
void UniformTimeTreeDistribution::attachTimes(
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


/** Build random binary tree to size numTaxa. The result is a draw from the uniform distribution on histories. */
void UniformTimeTreeDistribution::buildRandomBinaryHistory(std::vector<TopologyNode*> &tips) {
    
    if (tips.size() < numTaxa)
    {
        // Get the rng
        RandomNumberGenerator* rng = GLOBAL_RNG;

        // Randomly draw one node from the list of tips
        size_t index = static_cast<size_t>( floor(rng->uniform01()*tips.size()) );

        // Get the node from the list
        TopologyNode* parent = tips.at(index);
        
        // Remove the randomly drawn node from the list
        tips.erase(tips.begin()+index);
        
        // Add a left child
        TopologyNode* leftChild = new TopologyNode(0);
        parent->addChild(leftChild);
        leftChild->setParent(parent);
        tips.push_back(leftChild);
        
        // Add a right child
        TopologyNode* rightChild = new TopologyNode(0);
        parent->addChild(rightChild);
        rightChild->setParent(parent);
        tips.push_back(rightChild);
        
        // Recursive call to this function
        buildRandomBinaryHistory(tips);
    }
}


/* Clone function */
UniformTimeTreeDistribution* UniformTimeTreeDistribution::clone( void ) const {
    
    return new UniformTimeTreeDistribution( *this );
}


/* Compute probability */
double UniformTimeTreeDistribution::computeLnProbability( void ) {
    
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


void UniformTimeTreeDistribution::redrawValue( void ) {
    simulateTree();
}


/** Simulate the tree conditioned on the time of origin */
void UniformTimeTreeDistribution::simulateTree( void ) {
    
    // Get the rng
    RandomNumberGenerator* rng = GLOBAL_RNG;

    // Create the time tree object (topology + times)
    TimeTree *psi = new TimeTree();

    // Create an empty topology
    Topology *tau = new Topology();

    // Root the topology by setting the appropriate flag
    tau->setRooted( true );
    
    // Create the root node and a vector of nodes
    TopologyNode* root = new TopologyNode();
    std::vector<TopologyNode* > nodes;
    nodes.push_back(root);
    
    // Draw a random tree history
    buildRandomBinaryHistory(nodes);
    
    // Set the tip names
    for (size_t i=0; i<numTaxa; i++) {
        size_t index = size_t( floor(rng->uniform01() * nodes.size()) );
        
        // Get the node from the list
        TopologyNode* node = nodes.at(index);
        
        // Remove the randomly drawn node from the list
        nodes.erase(nodes.begin()+index);
        
        // Set name
        std::string& name = taxonNames[i];
        node->setName(name);
    }
    
    // Initialize the topology by setting the root
    tau->setRoot(root);
    
    // Connect the tree with the topology
    psi->setTopology( tau, true );
    
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
    
    // Finally store the new value
    value = psi;
    
}


void UniformTimeTreeDistribution::swapParameter(const DagNode *oldP, const DagNode *newP) {
    if (oldP == originTime) {
        originTime = static_cast<const TypedDagNode<double>* >( newP );
    }
}

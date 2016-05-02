#include "RandomNumberFactory.h"
#include "RandomNumberGenerator.h"
#include "TypedDagNode.h"

#include <cmath>
#include <iostream>

#include "PartialVectorScaleMove.h"
#include "Exception.h"

/** 
 * Constructor
 *
 * Here we simply allocate and initialize the move object.
 */
PartialVectorScaleMove::PartialVectorScaleMove( StochasticNode<double>* p, std::vector<ContinuousStochasticNode*> br, double l, bool t, double w ) : Move( p, w, t ),
        lambda( l ),
        partial(p),
        storedpartial(0.0),
        variables(br),
        storedvariables(std::vector<double>(br.size(), 0.0))
{
    // we need to allocate memory for the stored value
    for(size_t i = 0; i < variables.size(); i++)
    {
        this->nodes.insert(variables[i]);
    }
}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'B'.
 *
 * \return A new copy of the model. 
 */
PartialVectorScaleMove* PartialVectorScaleMove::clone( void ) const
{

    return new PartialVectorScaleMove( *this );
}


/**
 * Get moves' name of object 
 *
 * \return The moves' name.
 */
const std::string& PartialVectorScaleMove::getMoveName( void ) const
{
    static std::string name = "VectorScale";
    
    return name;
}


/** 
 * Perform the move.
 *
 * A scaling move draws a random uniform number u ~ unif(-0.5,0.5)
 * and scales the current vale by a scaling factor
 * sf = exp( lambda * u )
 * where lambda is the tuning parameter of the move to influence the size of the proposals.
 *
 * \return The hastings ratio.
 */
double PartialVectorScaleMove::performMove( double& probRatio )
{
        
    // Get random number generator    
    RandomNumberGenerator* rng     = GLOBAL_RNG;

    double& val = partial->getValue();

    storedpartial = val;

    // Generate new value (no reflection, so we simply abort later if we propose value here outside of support)
    double u = rng->uniform01();
    double m = std::exp( lambda * ( u - 0.5 ) );

    double hr = log(m);


    u = rng->uniform01();

    double scaleFactor = 1.0;
    if(u < val)
    {
        // scale the upper fraction
        scaleFactor = val + m*(1.0 - val);
        val = val/scaleFactor;

        hr += log(val) - log(storedpartial);
    }
    else
    {
        // scale the lower fraction
        scaleFactor = m*val + (1.0 - val);
        val = m*val/scaleFactor;

        hr += log(1.0 - val) - log(1.0 - storedpartial);
    }

    hr += (variables.size() - 2)*log(scaleFactor);

    if ( hr != hr || hr == Constants::Double::inf )
    {
        return Constants::Double::neginf;
    }

    probRatio = 0.0;

    std::set<DagNode* > affectedNodes;

    // set the new branch lengths
    for(size_t i = 0; i < variables.size(); i++)
    {
        double& br = variables[i]->getValue();
        storedvariables[i] = br;
        br *= scaleFactor;

        variables[i]->touch();

        probRatio += variables[i]->getLnProbabilityRatio();

        if ( probRatio != Constants::Double::inf && probRatio != Constants::Double::neginf )
        {
            variables[i]->getAffectedNodes(affectedNodes);
        }
    }

    // touch the node
    partial->touch();

    // calculate the probability ratio for the node we just changed
    probRatio += partial->getLnProbabilityRatio();

    if ( probRatio != Constants::Double::inf && probRatio != Constants::Double::neginf )
    {
        partial->getAffectedNodes(affectedNodes);

        for (std::set<DagNode* >::iterator i=affectedNodes.begin(); i!=affectedNodes.end(); ++i)
        {
            DagNode* theAffectedNode = *i;
            //std::cout << theAffectedNode->getName() << "  " << theAffectedNode->getLnProbabilityRatio() << " " << theAffectedNode->getLnProbability() << "\n";
            probRatio += theAffectedNode->getLnProbabilityRatio();
        }
    }

    return hr;
}


/**
 * Print the summary of the move.
 *
 * The summary just contains the current value of the tuning parameter.
 * It is printed to the stream that it passed in.
 *
 * \param[in]     o     The stream to which we print the summary.
 */
void PartialVectorScaleMove::printParameterSummary(std::ostream &o) const
{
    
    o << "lambda = " << lambda;

}

void PartialVectorScaleMove::acceptMove( void ) {

}


/**
 * Reject the move.
 *
 * Since the move stores the previous value and it is the only place
 * where complex undo operations are known/implement, we need to revert
 * the value of the variable/DAG-node to its original value.
 */
void PartialVectorScaleMove::rejectMove( void )
{
    partial->setValue(storedpartial);
    partial->touch();

    for(size_t i = 0; i < variables.size(); i++)
    {
        variables[i]->setValue(storedvariables[i]);
        variables[i]->touch();
    }
}


/**
 * Swap the current variable for a new one.
 *
 * \param[in]     oldN     The old variable that needs to be replaced.
 * \param[in]     newN     The new variable.
 */
void PartialVectorScaleMove::swapNode(DagNode *oldN, DagNode *newN)
{
    Move::swapNode(oldN, newN);

    // call the parent method
    if(oldN == partial)
    {
        partial = static_cast<StochasticNode<double>* >(newN);
    }
    else
    {
        for(size_t i = 0; i < variables.size(); i++)
        {
            if(oldN == variables[i])
            {
                variables[i] = static_cast<ContinuousStochasticNode* >(newN);
            }
        }
    }
    
}

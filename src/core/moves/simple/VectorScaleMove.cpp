#include "RandomNumberFactory.h"
#include "RandomNumberGenerator.h"
#include "TypedDagNode.h"

#include <cmath>
#include <iostream>

#include "VectorScaleMove.h"
#include "Exception.h"

/** 
 * Constructor
 *
 * Here we simply allocate and initialize the move object.
 */
VectorScaleMove::VectorScaleMove( std::vector<ContinuousStochasticNode*> br, double l, bool t, double w ) : SimpleMove( br.front(), w, t ),
        lambda( l ),
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
VectorScaleMove* VectorScaleMove::clone( void ) const
{

    return new VectorScaleMove( *this );
}


/**
 * Get moves' name of object 
 *
 * \return The moves' name.
 */
const std::string& VectorScaleMove::getMoveName( void ) const
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
double VectorScaleMove::performSimpleMove( void )
{
        
    // Get random number generator    
    RandomNumberGenerator* rng     = GLOBAL_RNG;

    // Generate new value (no reflection, so we simply abort later if we propose value here outside of support)
    double u = rng->uniform01();
    double m = std::exp( lambda * ( u - 0.5 ) );

    // set the new branch lengths
    for(size_t i = 0; i < variables.size(); i++)
    {
        double& br = variables[i]->getValue();
        storedvariables[i] = br;
        br *= m;

        variables[i]->touch();
    }
    
    // compute the Hastings ratio
    double lnHastingsratio = log(m)*variables.size();

    return lnHastingsratio;
}


/**
 * Print the summary of the move.
 *
 * The summary just contains the current value of the tuning parameter.
 * It is printed to the stream that it passed in.
 *
 * \param[in]     o     The stream to which we print the summary.
 */
void VectorScaleMove::printParameterSummary(std::ostream &o) const
{
    
    o << "lambda = " << lambda;

}


/**
 * Reject the move.
 *
 * Since the move stores the previous value and it is the only place
 * where complex undo operations are known/implement, we need to revert
 * the value of the variable/DAG-node to its original value.
 */
void VectorScaleMove::rejectSimpleMove( void )
{
    for(size_t i = 0; i < variables.size(); i++)
    {
        variables[i]->setValue(new double(storedvariables[i]));
    }
}


/**
 * Swap the current variable for a new one.
 *
 * \param[in]     oldN     The old variable that needs to be replaced.
 * \param[in]     newN     The new variable.
 */
void VectorScaleMove::swapNode(DagNode *oldN, DagNode *newN)
{
    // call the parent method
    SimpleMove::swapNode(oldN, newN);
    
    for(size_t i = 0; i < variables.size(); i++)
    {
        if(oldN == variables[i])
        {
            variables[i] = static_cast<ContinuousStochasticNode* >(newN);
        }
    }
    
}


/**
 * Tune the move to accept the desired acceptance ratio.
 *
 * The acceptance ratio for this move should be around 0.44.
 * If it is too large, then we increase the proposal size,
 * and if it is too small, then we decrease the proposal size.
 */
void VectorScaleMove::tune( void ) {
    
    double rate = numAccepted / double(numTried);
    
    if ( rate > 0.44 ) 
    {
        lambda *= (1.0 + ((rate-0.44)/0.56) );
    }
    else 
    {
        lambda /= (2.0 - rate/0.44 );
    }
    
}


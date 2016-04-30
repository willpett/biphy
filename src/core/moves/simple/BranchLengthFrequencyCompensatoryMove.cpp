#include "RandomNumberFactory.h"
#include "RandomNumberGenerator.h"
#include "TypedDagNode.h"

#include <cmath>
#include <iostream>

#include "BranchLengthFrequencyCompensatoryMove.h"
#include "Exception.h"

/** 
 * Constructor
 *
 * Here we simply allocate and initialize the move object.
 */
BranchLengthFrequencyCompensatoryMove::BranchLengthFrequencyCompensatoryMove( StochasticNode<double> *pi, std::vector<ContinuousStochasticNode*> br, double l, bool t, double w ) : SimpleMove( pi, w, t ),
        lambda( l ),
        frequency( pi ),
        storedFrequency( 0.0 ),
        branchlengths(br),
        storedBranchlengths(std::vector<double>(br.size(), 0.0))
{
    // we need to allocate memory for the stored value
    for(size_t i = 0; i < branchlengths.size(); i++)
    {
        this->nodes.insert(branchlengths[i]);
    }
}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'B'.
 *
 * \return A new copy of the model. 
 */
BranchLengthFrequencyCompensatoryMove* BranchLengthFrequencyCompensatoryMove::clone( void ) const
{

    return new BranchLengthFrequencyCompensatoryMove( *this );
}


/**
 * Get moves' name of object 
 *
 * \return The moves' name.
 */
const std::string& BranchLengthFrequencyCompensatoryMove::getMoveName( void ) const
{
    static std::string name = "BranchFrequency";
    
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
double BranchLengthFrequencyCompensatoryMove::performSimpleMove( void )
{
        
    // Get random number generator    
    RandomNumberGenerator* rng     = GLOBAL_RNG;
    
    double &val = frequency->getValue();
    
    // copy value
    storedFrequency = val;

    // Generate new value (no reflection, so we simply abort later if we propose value here outside of support)
    double u = rng->uniform01();
    double m = std::exp( lambda * ( u - 0.5 ) );

    double scalingFactor = (1.0 + val*(m - 1.0));

    // compute the Hastings ratio
    double lnHastingsratio = 2.0*log(m) + log( m - val * (1.0 - val)*(m - 1.0)*(m - 1.0));
    lnHastingsratio -= (branchlengths.size() - 4)*log(scalingFactor);
    lnHastingsratio = abs(lnHastingsratio);

    // set the new frequency
    val *= m/scalingFactor;

    // set the new branch lengths
    for(size_t i = 0; i < branchlengths.size(); i++)
    {
        double& br = branchlengths[i]->getValue();
        storedBranchlengths[i] = br;
        br *= scalingFactor;

        branchlengths[i]->touch();
    }
    
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
void BranchLengthFrequencyCompensatoryMove::printParameterSummary(std::ostream &o) const
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
void BranchLengthFrequencyCompensatoryMove::rejectSimpleMove( void )
{
    // swap current value and stored value
    frequency->setValue( new double(storedFrequency) );
    
    for(size_t i = 0; i < branchlengths.size(); i++)
    {
        branchlengths[i]->setValue(new double(storedBranchlengths[i]));
    }
}


/**
 * Swap the current variable for a new one.
 *
 * \param[in]     oldN     The old variable that needs to be replaced.
 * \param[in]     newN     The new variable.
 */
void BranchLengthFrequencyCompensatoryMove::swapNode(DagNode *oldN, DagNode *newN)
{
    // call the parent method
    SimpleMove::swapNode(oldN, newN);
    
    if(oldN == frequency)
    {
        frequency = static_cast<StochasticNode<double>* >(newN) ;
    }
    else
    {
        for(size_t i = 0; i < branchlengths.size(); i++)
        {
            if(oldN == branchlengths[i])
            {
                branchlengths[i] = static_cast<ContinuousStochasticNode* >(newN);
            }
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
void BranchLengthFrequencyCompensatoryMove::tune( void ) {
    
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


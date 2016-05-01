/**
 * @file
 * This file contains the implementation of SlidingMove,
 * which moves a real value uniformly within a sliding window.
 *
 * @brief Implementation of SlidingMove
 *
 * (c) Copyright 2009- under GPL version 3
 * @date Last modified: $Date$
 * @author The RevBayes Development Core Team
 * @license GPL version 3
 * @version 1.0
 * @since 2012-06-21, version 1.0
 *
 * $Id$
 */

#include "SlidingMove.h"
#include "RandomNumberFactory.h"
#include "RandomNumberGenerator.h"
#include "TypedDagNode.h"
#include "ContinuousDistribution.h"
#include <cmath>
#include <cassert>
#include "Exception.h"

/* Constructor with value */
SlidingMove::SlidingMove( StochasticNode<double> *n, double d, bool t, double w ) : SimpleMove( n, w, t ), delta( d ), variable( n ),storedValue( 0.0 ) {
    
    // we need to allocate memory for the stored value
}


/* Clone object */
SlidingMove* SlidingMove::clone( void ) const {

    return new SlidingMove( *this );
}



const std::string& SlidingMove::getMoveName( void ) const {
    static std::string name = "Sliding";
    
    return name;
}


/** Perform the move */
double SlidingMove::performSimpleMove( void ) {

    // Get random number generator    
    RandomNumberGenerator* rng     = GLOBAL_RNG;
    
    double &val = variable->getValue();
    
    // copy all value
    storedValue = val;

    double min = static_cast<ContinuousDistribution *>( &(variable->getDistribution()) )->getMin();
    double max = static_cast<ContinuousDistribution *>( &(variable->getDistribution()) )->getMax();

    double u      = rng->uniform01();
    double newVal = val + ( delta * ( u - 0.5 ) );

    /* reflect the new value */
    do {
        if ( newVal < min )
            newVal = 2.0 * min - newVal;
        else if ( newVal > max )
            newVal = 2.0 * max - newVal;
    } while ( newVal < min || newVal > max );

    // FIXME: not the most efficient way of handling multiple reflections :-P

    val = newVal;
	
    return 0.0;
}


void SlidingMove::printParameterSummary(std::ostream &o) const {
    o << "delta = " << delta;
}


void SlidingMove::rejectSimpleMove( void ) {
    // swap current value and stored value
    variable->setValue( new double(storedValue) );
}


void SlidingMove::swapNode(DagNode *oldN, DagNode *newN) {
    // call the parent method
    SimpleMove::swapNode(oldN, newN);
    
    variable = static_cast<StochasticNode<double>* >(newN);
}


void SlidingMove::tune( void ) 
{
}



//
//  Add.cpp
//  RevBayesCore
//
//  Created by Sebastian Hoehna on 8/13/12.
//  Copyright 2012 __MyCompanyName__. All rights reserved.
//

#include "LnFunction.h"
#include <cmath>

using namespace RevBayesCore;

LnFunction::LnFunction(const TypedDagNode<double> *x) : ContinuousFunction( new double(0.0) ), a( x ) {
    addParameter( x );

}


LnFunction* LnFunction::clone( void ) const {
    return new LnFunction(*this);
}


void LnFunction::swapParameterInternal(const DagNode *oldP, const DagNode *newP) {
    if (oldP == a) {
        a = static_cast<const TypedDagNode<double>* >( newP );
    }
}

void LnFunction::update( void ) {
    *value = log( a->getValue() );
}



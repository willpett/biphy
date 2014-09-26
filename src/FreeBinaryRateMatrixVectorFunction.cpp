//
//  FreeBinaryRateMatrixFunction.cpp
//  rb_mlandis
//
//  Created by Michael Landis on 4/4/14.
//  Copyright (c) 2014 Michael Landis. All rights reserved.
//

#include "FreeBinaryRateMatrixVectorFunction.h"
#include "RbException.h"

using namespace RevBayesCore;

FreeBinaryRateMatrixVectorFunction::FreeBinaryRateMatrixVectorFunction(const TypedDagNode<std::vector<double> > *tr) : TypedFunction<RbVector<RateMatrix> >( new RbVector<RateMatrix>() ), transitionRates( tr ) {
    // add the lambda parameter as a parent
    addParameter( transitionRates );
    const std::vector<double>& r = transitionRates->getValue();
    for (size_t i = 0; i < r.size(); ++i) {
    	value->push_back(new RateMatrix_FreeBinary());
    }
    
    update();
}


FreeBinaryRateMatrixVectorFunction::FreeBinaryRateMatrixVectorFunction(const FreeBinaryRateMatrixVectorFunction &n) : TypedFunction<RbVector<RateMatrix> >( n ), transitionRates( n.transitionRates ) {
    // no need to add parameters, happens automatically
}


FreeBinaryRateMatrixVectorFunction::~FreeBinaryRateMatrixVectorFunction( void ) {
    // We don't delete the parameters, because they might be used somewhere else too. The model needs to do that!
}



FreeBinaryRateMatrixVectorFunction* FreeBinaryRateMatrixVectorFunction::clone( void ) const {
    return new FreeBinaryRateMatrixVectorFunction( *this );
}


void FreeBinaryRateMatrixVectorFunction::update( void ) {
    // get the information from the arguments for reading the file
    const std::vector<double>& r = transitionRates->getValue();
    for (size_t i = 0; i < value->size(); ++i) {
    	std::vector<double> pi;
    	pi.push_back(r[i]);
    	pi.push_back(1.0-r[i]);
    	((RateMatrix_FreeBinary&)(*value)[i]).setTransitionRates(r);
    	(*value)[i].updateMatrix();
    }
}



void FreeBinaryRateMatrixVectorFunction::swapParameterInternal(const DagNode *oldP, const DagNode *newP) {
    if (oldP == transitionRates) {
        transitionRates = static_cast<const TypedDagNode<std::vector<double> >* >( newP );
    }
}


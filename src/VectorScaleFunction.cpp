//
//  VectorScaleFunction.cpp
//  RevBayesCore
//
//  Created by Sebastian Hoehna on 11/15/12.
//  Copyright 2012 __MyCompanyName__. All rights reserved.
//

#include "VectorScaleFunction.h"
#include "TypedDagNode.h"

using namespace RevBayesCore;

VectorScaleFunction::VectorScaleFunction(const TypedDagNode<std::vector<double> > *vec, const TypedDagNode<double> *sca) : TypedFunction<std::vector<double> >( new std::vector<double>() ), vector( vec ), scalar(sca) {
    // add the lambda parameter as a parent
    addParameter( vector);
    addParameter( scalar);
    
    update();
}


VectorScaleFunction::VectorScaleFunction(const VectorScaleFunction &n) : TypedFunction<std::vector<double> >( n ), vector(n.vector), scalar(n.scalar) {
    // no need to add parameters, happens automatically
    
    update();
}


VectorScaleFunction::~VectorScaleFunction( void ) {
    // We don't delete the parameters, because they might be used somewhere else too. The model needs to do that!
}



VectorScaleFunction* VectorScaleFunction::clone( void ) const {
    return new VectorScaleFunction( *this );
}


void VectorScaleFunction::update( void ) {
    
    // empty current simplex
    value->clear();
    
    std::vector<double> vec_values = vector->getValue();
    double scale = scalar->getValue();
    for (size_t i = 0; i < vec_values.size(); i++) {
        value->push_back(vec_values[i]*scale);
    }
}

void VectorScaleFunction::swapParameterInternal(const DagNode *oldP, const DagNode *newP) {
    
	if (oldP == vector) {
		vector = static_cast<const TypedDagNode<std::vector<double> >* >( newP );
	}else if(oldP == scalar){
		scalar = static_cast<const TypedDagNode<double>* >( newP );
	}
    
}

#include "EntropyFunction.h"
#include <cmath>

EntropyFunction::EntropyFunction(const TypedDagNode<std::vector<double> > *v) : TypedFunction<double>( new double(0.0) ), vals( v ) {
    // add the parameters as parents
    this->addParameter( vals );
    
    update();
}


EntropyFunction::EntropyFunction(const EntropyFunction &n) : TypedFunction<double>( n ), vals( n.vals ) {
    // no need to add parameters, happens automatically
    
    update();
}


EntropyFunction::~EntropyFunction( void ) {
    // We don't delete the parameters, because they might be used somewhere else too. The model needs to do that!
}



EntropyFunction* EntropyFunction::clone( void ) const {
    return new EntropyFunction( *this );
}


void EntropyFunction::update( void ) {
    
    double m = 0;
    const std::vector<double> &v = vals->getValue();
    for ( std::vector<double>::const_iterator it = v.begin(); it != v.end(); ++it) {
        m -= (*it)*log(*it);
    }
    
    *this->value = m;
    
}



void EntropyFunction::swapParameterInternal(const DagNode *oldP, const DagNode *newP) {
    
    if ( oldP == vals ) {
        vals = static_cast<const TypedDagNode<std::vector<double> >* >( newP );
    }
    
}


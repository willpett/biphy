#include "TruncatedDistribution.h"
#include "RandomNumberFactory.h"
#include "RbConstants.h"

using namespace RevBayesCore;

TruncatedDistribution::TruncatedDistribution(TypedDistribution<double> *fin, const TypedDagNode<double> *mi, const TypedDagNode<double> *ma) : TypedDistribution<double>( new double( 0.0 ) ), min( mi ), max( ma ), f( fin ){
    // add the parameters to the parents set
    addParameter( min );
    addParameter( max );
    
    f->redrawValue();
	while(f->getValue() < min->getValue() || f->getValue() > max->getValue()){
		f->redrawValue();
	}
	*value = f->getValue();
}


TruncatedDistribution::TruncatedDistribution(const TruncatedDistribution &n) : TypedDistribution<double>( n ), min( n.min ), max( n.max ), f( n.f ) {
    // parameters are automatically copied
}


TruncatedDistribution::~TruncatedDistribution( void ) {
    // We don't delete the parameters, because they might be used somewhere else too. The model needs to do that!
}


TruncatedDistribution* TruncatedDistribution::clone( void ) const {
    
    return new TruncatedDistribution( *this );
}


double TruncatedDistribution::computeLnProbability( void ) {
    
	if(*value < min->getValue() || *value > max->getValue())
		return RbConstants::Double::neginf;
	f->setValue(*value);
    return f->computeLnProbability();
}


double TruncatedDistribution::getMax( void ) const {
    
    return max->getValue();
}


double TruncatedDistribution::getMin( void ) const {
    
    return min->getValue();
}


void TruncatedDistribution::redrawValue( void ) {
    
	f->redrawValue();
	while(f->getValue() < min->getValue() || f->getValue() > max->getValue()){
		f->redrawValue();
	}
    *value = f->getValue();
    
}


void TruncatedDistribution::swapParameter(const DagNode *oldP, const DagNode *newP) {
    
    if (oldP == min) 
    {
        min = static_cast<const TypedDagNode<double>* >( newP );
    }
    else if (oldP == max) 
    {
        max = static_cast<const TypedDagNode<double>* >( newP );
    }
    
}

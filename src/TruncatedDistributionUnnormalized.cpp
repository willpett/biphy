#include "TruncatedDistributionUnnormalized.h"
#include "RandomNumberFactory.h"
#include "RbConstants.h"

using namespace RevBayesCore;

TruncatedDistributionUnnormalized::TruncatedDistributionUnnormalized(TypedDistribution<double> *fin, const TypedDagNode<double> *mi, const TypedDagNode<double> *ma) : TypedDistribution<double>( new double( 0.0 ) ), min( mi ), max( ma ), f( fin ){
    // add the parameters to the parents set
    addParameter( min );
    addParameter( max );
    
    f->redrawValue();
	while(f->getValue() < min->getValue() || f->getValue() > max->getValue()){
		f->redrawValue();
	}
	*value = f->getValue();
}


TruncatedDistributionUnnormalized::TruncatedDistributionUnnormalized(const TruncatedDistributionUnnormalized &n) : TypedDistribution<double>( n ), min( n.min ), max( n.max ), f( n.f->clone() ) {
    // parameters are automatically copied
}


TruncatedDistributionUnnormalized::~TruncatedDistributionUnnormalized( void ) {
    // We don't delete the parameters, because they might be used somewhere else too. The model needs to do that!
}


TruncatedDistributionUnnormalized* TruncatedDistributionUnnormalized::clone( void ) const {
    
    return new TruncatedDistributionUnnormalized( *this );
}


double TruncatedDistributionUnnormalized::computeLnProbability( void ) {
    
	if(*value < min->getValue() || *value > max->getValue())
		return RbConstants::Double::neginf;
	f->setValue(*value);
    return f->computeLnProbability();
}


double TruncatedDistributionUnnormalized::getMax( void ) const {
    
    return max->getValue();
}


double TruncatedDistributionUnnormalized::getMin( void ) const {
    
    return min->getValue();
}


void TruncatedDistributionUnnormalized::redrawValue( void ) {
    
	f->redrawValue();
	while(f->getValue() < min->getValue() || f->getValue() > max->getValue()){
		f->redrawValue();
	}
    *value = f->getValue();
    
}


void TruncatedDistributionUnnormalized::swapParameter(const DagNode *oldP, const DagNode *newP) {
    
    if (oldP == min) 
    {
        min = static_cast<const TypedDagNode<double>* >( newP );
    }
    else if (oldP == max) 
    {
        max = static_cast<const TypedDagNode<double>* >( newP );
    }
    
}


#include "NormalDistribution.h"
#include "DistributionNormal.h"
#include "RandomNumberFactory.h"
#include "Constants.h"

NormalDistribution::NormalDistribution(const TypedDagNode<double> *m, const TypedDagNode<double> *s) : ContinuousDistribution( new double( 0.0 ) ),
    mean( m ),
    stDev( s )
{
    // add the parameters to our set (in the base class)
    // in that way other class can easily access the set of our parameters
    // this will also ensure that the parameters are not getting deleted before we do
    addParameter( mean );
    addParameter( stDev );
    
    *value = Statistics::Normal::rv(mean->getValue(), stDev->getValue(), *GLOBAL_RNG);
}


NormalDistribution::~NormalDistribution( void ) {
    // We don't delete the parameters, because they might be used somewhere else too. The model needs to do that!
}


double NormalDistribution::cdf( void ) const {
    return Statistics::Normal::cdf( mean->getValue(), stDev->getValue(), *value);
}


NormalDistribution* NormalDistribution::clone( void ) const {
    return new NormalDistribution( *this );
}


double NormalDistribution::computeLnProbability( void ) {
    return Statistics::Normal::lnPdf(mean->getValue(), stDev->getValue(), *value);
}


double NormalDistribution::getMax( void ) const {
    return Constants::Double::inf;
}


double NormalDistribution::getMin( void ) const {
    return Constants::Double::neginf;
}


double NormalDistribution::quantile(double p) const {
    return Statistics::Normal::quantile(mean->getValue(), stDev->getValue(), p);
}


void NormalDistribution::redrawValue( void ) {
    *value = Statistics::Normal::rv(mean->getValue(), stDev->getValue(), *GLOBAL_RNG);
}


/** Swap a parameter of the distribution */
void NormalDistribution::swapParameter(const DagNode *oldP, const DagNode *newP)
{
    if (oldP == mean)
    {
        mean = static_cast<const TypedDagNode<double>* >( newP );
    }
    if (oldP == stDev)
    {
        stDev = static_cast<const TypedDagNode<double>* >( newP );
    }
}

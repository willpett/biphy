#include "PoissonDistribution.h"

#include "Constants.h"
#include "DistributionPoisson.h"
#include "RandomNumberFactory.h"

PoissonDistribution::PoissonDistribution(const TypedDagNode<double> *l) : TypedDistribution<int>( new int( 1 ) ),
    lambda( l )
{
    // add the parameters to our set (in the base class)
    // in that way other class can easily access the set of our parameters
    // this will also ensure that the parameters are not getting deleted before we do
    addParameter( lambda );
    
    *value = Statistics::Poisson::rv(lambda->getValue(), *GLOBAL_RNG);
}


PoissonDistribution::~PoissonDistribution( void ) {
    // We don't delete the parameters, because they might be used somewhere else too. The model needs to do that!
}


PoissonDistribution* PoissonDistribution::clone( void ) const {
    return new PoissonDistribution( *this );
}


double PoissonDistribution::computeLnProbability( void )
{
    
    double v = *value;
    
    // check that the value is inside the boundaries
    if ( v < 0.0 )
    {
        return Constants::Double::neginf;
    }
    
    return Statistics::Poisson::lnPdf(lambda->getValue(), *value);
}


void PoissonDistribution::redrawValue( void ) {
    *value = Statistics::Poisson::rv(lambda->getValue(), *GLOBAL_RNG);
}


/** Swap a parameter of the distribution */
void PoissonDistribution::swapParameter(const DagNode *oldP, const DagNode *newP)
{
    if (oldP == lambda)
    {
        lambda = static_cast<const TypedDagNode<double>* >( newP );
    }
}

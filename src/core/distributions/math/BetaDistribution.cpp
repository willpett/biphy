#include "BetaDistribution.h"
#include "DistributionBeta.h"
#include "RandomNumberFactory.h"
#include <limits>
#include "MathLogic.h"

BetaDistribution::BetaDistribution(const TypedDagNode<double> *a, const TypedDagNode<double> *b) : ContinuousDistribution( new double( 0.5 ) ), alpha( a ), beta( b ) {
    // add the parameters to the parents set
    addParameter( alpha );
    addParameter( beta );
    
    redrawValue();
}


BetaDistribution::~BetaDistribution( void ) {
    // We don't delete the parameters, because they might be used somewhere else too. The model needs to do that!
}


double BetaDistribution::cdf( void ) const {
    return Statistics::Beta::cdf(alpha->getValue(), beta->getValue(), *value);
}


BetaDistribution* BetaDistribution::clone( void ) const {
    return new BetaDistribution( *this );
}


double BetaDistribution::computeLnProbability( void ) {
    return Statistics::Beta::lnPdf(alpha->getValue(), beta->getValue(), *value);
}


double BetaDistribution::getMax( void ) const {
    return 1.0;
}


double BetaDistribution::getMin( void ) const {
    return 0.0;
}


double BetaDistribution::quantile(double p) const {
    return Statistics::Beta::quantile(alpha->getValue(), beta->getValue(), p);
}


void BetaDistribution::redrawValue( void ) {
    *value = Statistics::Beta::rv(alpha->getValue(), beta->getValue(), *GLOBAL_RNG);
}


void BetaDistribution::swapParameter(const DagNode *oldP, const DagNode *newP) {
    if (oldP == alpha) {
        alpha = static_cast<const TypedDagNode<double>* >( newP );
    }
    if (oldP == beta) {
        beta = static_cast<const TypedDagNode<double>* >( newP );
    }
}

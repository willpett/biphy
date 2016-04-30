#include "GammaDistribution.h"

#include "MathLogic.h"
#include "Constants.h"
#include "DistributionGamma.h"
#include "RandomNumberFactory.h"

GammaDistribution::GammaDistribution(const TypedDagNode<double> *sh, const TypedDagNode<double> *r) : ContinuousDistribution( new double( 1.0 ) ), shape( sh ), rate( r ) {
    // add the parameters to the parents set
    addParameter( shape );
    addParameter( rate );
    
    redrawValue();
}


GammaDistribution::~GammaDistribution( void ) {
    // We don't delete the parameters, because they might be used somewhere else too. The model needs to do that!
}


double GammaDistribution::cdf( void ) const {
    return Statistics::Gamma::cdf(shape->getValue(), rate->getValue(), *value);
}


GammaDistribution* GammaDistribution::clone( void ) const {
    return new GammaDistribution( *this );
}


double GammaDistribution::computeLnProbability( void ) {
    return Statistics::Gamma::lnPdf(shape->getValue(), rate->getValue(), *value);
}


double GammaDistribution::getMax( void ) const {
    return Constants::Double::inf;
}


double GammaDistribution::getMin( void ) const {
    return 0.0;
}


double GammaDistribution::quantile(double p) const {
    return Statistics::Gamma::quantile(shape->getValue(), rate->getValue(), p);
}


void GammaDistribution::redrawValue( void ) {
    *value = Statistics::Gamma::rv(shape->getValue(), rate->getValue(), *GLOBAL_RNG);
    while(	Math::compEssentiallyEqual(*value,1.0,std::numeric_limits<double>::epsilon()) ||
    		Math::compEssentiallyEqual(*value,0.0,std::numeric_limits<double>::epsilon()) ||
    		Math::isNan(*value))
    	*value = Statistics::Gamma::rv(shape->getValue(), rate->getValue(), *GLOBAL_RNG);
}


void GammaDistribution::swapParameter(const DagNode *oldP, const DagNode *newP) {
    if (oldP == shape) {
        shape = static_cast<const TypedDagNode<double>* >( newP );
    }
    if (oldP == rate) {
        rate = static_cast<const TypedDagNode<double>* >( newP );
    }
}


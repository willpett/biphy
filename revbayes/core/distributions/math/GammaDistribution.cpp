#include "GammaDistribution.h"
#include "DistributionGamma.h"
#include "RandomNumberFactory.h"
#include "RbConstants.h"
#include "RbMathLogic.h"

using namespace RevBayesCore;

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
    return RbStatistics::Gamma::cdf(shape->getValue(), rate->getValue(), *value);
}


GammaDistribution* GammaDistribution::clone( void ) const {
    return new GammaDistribution( *this );
}


double GammaDistribution::computeLnProbability( void ) {
    return RbStatistics::Gamma::lnPdf(shape->getValue(), rate->getValue(), *value);
}


double GammaDistribution::getMax( void ) const {
    return RbConstants::Double::inf;
}


double GammaDistribution::getMin( void ) const {
    return 0.0;
}


double GammaDistribution::quantile(double p) const {
    return RbStatistics::Gamma::quantile(shape->getValue(), rate->getValue(), p);
}


void GammaDistribution::redrawValue( void ) {
    *value = RbStatistics::Gamma::rv(shape->getValue(), rate->getValue(), *GLOBAL_RNG);
    while(	RbMath::compEssentiallyEqual(*value,1.0,std::numeric_limits<double>::epsilon()) ||
    		RbMath::compEssentiallyEqual(*value,0.0,std::numeric_limits<double>::epsilon()) ||
    		RbMath::isNan(*value))
    	*value = RbStatistics::Gamma::rv(shape->getValue(), rate->getValue(), *GLOBAL_RNG);
}


void GammaDistribution::swapParameter(const DagNode *oldP, const DagNode *newP) {
    if (oldP == shape) {
        shape = static_cast<const TypedDagNode<double>* >( newP );
    }
    if (oldP == rate) {
        rate = static_cast<const TypedDagNode<double>* >( newP );
    }
}


#include "DistributionNegativeBinomial.h"
#include "NegativeBinomialDistribution.h"
#include "RandomNumberFactory.h"
#include "RbConstants.h"

using namespace RevBayesCore;

NegativeBinomialDistribution::NegativeBinomialDistribution(const TypedDagNode<int> *m, const TypedDagNode<double> *q) : TypedDistribution<int>( new int( 0 ) ),
    r( m ),
    p( q )
{
    
    // add the parameters to our set (in the base class)
    // in that way other class can easily access the set of our parameters
    // this will also ensure that the parameters are not getting deleted before we do
    addParameter( m );
    addParameter( q );
    
    if(p->getValue() == 1.0)
    	*value = RbConstants::Integer::inf;
    else
    	*value = RbStatistics::NegativeBinomial::rv(r->getValue(), p->getValue(), *GLOBAL_RNG);
}


NegativeBinomialDistribution::~NegativeBinomialDistribution(void) {

    // We don't delete the parameters, because they might be used somewhere else too. The model needs to do that!
}



NegativeBinomialDistribution* NegativeBinomialDistribution::clone( void ) const {

    return new NegativeBinomialDistribution( *this );
}


double NegativeBinomialDistribution::computeLnProbability( void )
{
    
    // check that the value is inside the boundaries
    if ( *value < 0 )
    {
        return RbConstants::Double::neginf;
    }
    
    if(p->getValue() == 1.0)
    	if(*value == RbConstants::Integer::inf){
    		return 0.0;
    	}else{
    		return RbConstants::Double::neginf;
    	}

    return RbStatistics::NegativeBinomial::lnPdf(r->getValue(), p->getValue(), *value);
}



void NegativeBinomialDistribution::redrawValue( void ) {

	if(p->getValue() == 1.0){
		*value = RbConstants::Integer::inf;
	}else{
		*value = RbStatistics::NegativeBinomial::rv(r->getValue(), p->getValue(), *GLOBAL_RNG);
	}
    //std::cerr << r->getValue() << "\t" << p->getValue() << "\t" << *value << std::endl;
}


/** Swap a parameter of the distribution */
void NegativeBinomialDistribution::swapParameter(const DagNode *oldP, const DagNode *newP) {

    if (oldP == p)
    {
        p = static_cast<const TypedDagNode<double>* >( newP );
    }
    else if (oldP == r)
    {
        r = static_cast<const TypedDagNode<int>* >( newP );
    }

}



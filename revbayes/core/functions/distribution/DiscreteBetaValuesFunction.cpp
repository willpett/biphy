#include "DiscreteBetaValuesFunction.h"
#include "DistributionBeta.h"

DiscreteBetaValuesFunction::DiscreteBetaValuesFunction(const TypedDagNode<double> *a, const TypedDagNode<double> *b, std::vector<double> inbreaks, bool inmean) :
	TypedFunction< std::vector<double> >( new std::vector<double>() ), alpha(a), beta(b), mean(inmean), breaks(inbreaks) {
    // add the lambda parameter as a parent
    this->addParameter( alpha );
    this->addParameter( beta );

    update();
}


DiscreteBetaValuesFunction::DiscreteBetaValuesFunction(const DiscreteBetaValuesFunction &n) : TypedFunction< std::vector<double> >( n ), alpha(n.alpha), beta(n.beta), breaks(n.breaks), mean(n.mean) {
    // no need to add parameters, happens automatically
    
    update();
}



DiscreteBetaValuesFunction::~DiscreteBetaValuesFunction( void ) {
    // We don't delete the parameters, because they might be used somewhere else too. The model needs to do that!
}




DiscreteBetaValuesFunction* DiscreteBetaValuesFunction::clone( void ) const {
    return new DiscreteBetaValuesFunction( *this );
}



void DiscreteBetaValuesFunction::update( void ) {
    
    // empty current vector
	this->value->clear();

	double a = alpha->getValue();
	double b = beta->getValue();

	size_t numVals = (breaks.size()+1);

	if(mean)
	{
		// use means

		std::vector<double> vals;

		for(size_t i = 0; i < breaks.size(); i++)
		{
			double quantile = Statistics::Beta::quantile(a, b, breaks[i]);

			vals.push_back(Statistics::Beta::cdf(a + 1, b, quantile));
			double interval = vals[i];
			if(i > 0)
			{
				interval -= vals[i-1];
			}
			this->value->push_back(numVals*interval*a/(a+b));
		}

		double interval = 1.0 - vals.back();
		this->value->push_back(numVals*interval*a/(a+b));
	}
	else
	{
		// use medians
		std::vector<double> quantiles;

		double mean = 0.0;
		for(size_t i = 0; i < breaks.size(); i++)
		{
			double last = i > 0 ? breaks[i-1] : 0.0;
			quantiles.push_back(Statistics::Beta::quantile(a, b, (breaks[i] - last)/2.0));
			mean += quantiles.back();
		}

		quantiles.push_back(Statistics::Beta::quantile(a, b, (1.0 - breaks.back())/2.0));
		mean += quantiles.back();

		mean /= double(numVals);

		for(size_t i = 0; i < quantiles.size(); i++)
			this->value->push_back(quantiles[i]/mean);
	}
    
}

void DiscreteBetaValuesFunction::swapParameterInternal(const DagNode *oldP, const DagNode *newP) {
    
    if (oldP == alpha)
    {
    	alpha = static_cast<const TypedDagNode<double>* >( newP );
	}
    else if (oldP == beta)
    {
		beta = static_cast<const TypedDagNode<double>* >( newP );
	}
    
}

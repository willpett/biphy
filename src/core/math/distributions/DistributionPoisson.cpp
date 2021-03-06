/**
 * @file DistributionPoisson
 * This file contains the functions of the poisson distribution.
 *
 * @brief Implementation of the poisson distribution.
 *
 * (c) Copyright 2009- under GPL version 3
 * @date Last modified: $Date$
 * @author The RevBayes core development team
 * @license GPL version 3
 * @version 1.0
 * @since 2011-03-17, version 1.0
 *
 * $Id$
 */

#include <cmath>
#include <iostream>

#include "DistributionPoisson.h"

#include "Exception.h"
#include "MathCombinatorialFunctions.h"
#include "StatisticsHelper.h"

double Statistics::Poisson::pdf(double lambda, int x) {
    
	return exp(x * std::log(lambda) - lambda - Math::lnFactorial(x));
}

/*!
 * This function calculates the natural log of the probability for a
 * Poisson distribution. 
 *
 * \brief Natural log of Poisson probability.
 * \param lambda is the rate parameter of the Poisson. 
 * \param x is the value of the random variable. 
 * \return Returns the natural log of the probability. 
 * \throws Does not throw an error.
 */
double Statistics::Poisson::lnPdf(double lambda, int x) {
    
    return ( x * std::log(lambda) - lambda - Math::lnFactorial(x) );
}

/*!
 * This function calculates the cumulative probability for a
 * Poisson distribution. 
 *
 * \brief Poisson cumulative probability.
 * \param lambda is the rate parameter of the Poisson. 
 * \param x is the value of the random variable. 
 * \return Returns the cumulative probability. 
 * \throws Does not throw an error.
 */
double Statistics::Poisson::cdf(double lambda, int x) {
    
	if ( x < 0 )
		return 0.0;
	double next = exp(-lambda);
	double cdf = next;
	for (int i=1; i<=x; i++)
        {
		double last = next;
		next = last * lambda / (double)i;
		cdf += next;
        }
	return cdf;
}

/*!
 * This function returns the quantile of a Poisson probability 
 * distribution.
 *
 * \brief Poisson(lambda) quantile.
 * \param lambda is the rate parameter of the Poisson. 
 * \param p is the probability up to the quantile. 
 * \return Returns the quantile.
 * \throws Does not throw an error.
 */
double Statistics::Poisson::quantile(double lambda, double p) {
    
	/* Starting with x = 0, find the first value for which
     CDF(X-1) <= CDF <= CDF(X). */
	double sum = 0.0;
	int xmax = 100;
	for (int i=0; i<=xmax; i++)
        {
		double sumOld = sum;
		double newVal = 0.0;
		if ( i == 0 )
            {
			newVal = exp(-lambda);
			sum = newVal;
            }
		else
            {
			double last = newVal;
			newVal = last * lambda / ( double ) ( i );
			sum += newVal;
            }
		if ( sumOld <= p && p <= sum )
			return i;
        }
    
	return xmax;
}

/*!
 * This function generates a Poisson-distributed random 
 * variable with parameter lambda.
 *
 * \brief Poisson(lambda) random variable.
 * \param lambda the rate parameter of the Poisson. 
 * \param rng is a pointer to a random number object. 
 * \return This function returns a Poisson-distributed integer.
 * \throws Does not throw an error.
 */
int Statistics::Poisson::rv(double lambda, RandomNumberGenerator& rng) {
    
	if (lambda < 17.0)
        {
		if (lambda < 1.0e-6)
            {
                if (lambda == 0.0) 
				return 0;
            
			if (lambda < 0.0)
                {
                std::ostringstream s;
                s << "Parameter negative in poisson function";
                throw (Exception(s));
                }
            
			/* For extremely small lambda we calculate the probabilities of x = 1
             and x = 2 (ignoring higher x). The reason for using this 
             method is to prevent numerical inaccuracies in other methods. */
            //			return Statistics::Helper::poissonLow(lambda, *rng);
                
            // MJL 071713: Poisson rv are supported on 0..inf, so why would a small rate return 1 instead of 0?
            // return 1;
            return 0;
            }
		else 
            {
			/* use the inversion method */
			return Statistics::Helper::poissonInver(lambda, rng);
            // MJL 071713: Same as above.
            return 0;
            //return 1;
            }
        }
	else 
        {
		if (lambda > 2.0e9) 
            {
			/* there should be an error here */
			throw Exception( "Parameter too big in poisson function" );
            }
		/* use the ratio-of-uniforms method */
        return Statistics::Helper::poissonRatioUniforms(lambda, rng);
        
        return 1;
        }
}


/**
 * @file DistributionGamma
 * This file contains the functions of the gamma distribution.
 *
 * @brief Implementation of the gamma distribution.
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

#include "DistributionGamma.h"

#include "Constants.h"
#include "MathFunctions.h"
#include "MathLogic.h"
#include "DistributionNormal.h"
#include "StatisticsHelper.h"

/*!
 * This function calculates the probability density 
 * for a gamma-distributed random variable.
 *
 * \brief Gamma probability density.
 * \param the shape parameter of the gamma. 
 * \param the rate parameter of the gamma. 
 * \param x is the gamma random variable. 
 * \return Returns the probability density.
 * \throws Does not throw an error.
 */
double Statistics::Gamma::pdf(double shape, double rate, double x) {

	return (pow(rate, shape) / Math::gamma(shape)) * pow(x, shape - 1.0) * exp(-x * rate);
}


/*!
 * This function calculates the probability density 
 * for a gamma-distributed random variable.
 *
 * \brief Gamma probability density.
 * \param the shape parameter of the gamma. 
 * \param the rate parameter of the gamma. 
 * \param x is the gamma random variable. 
 * \return Returns the probability density.
 * \throws Does not throw an error.
 */
double Statistics::Gamma::pdf(double shape, double rate, double x, bool isLog) {

    //    double pr;
    //    if (shape < 0 || scale <= 0) {
    //        std::ostringstream s;
    //        s << "Cannot compute the pdf for the gamma distribution for shape = " << shape << " and scale = " << scale;
    //        throw (Exception(s));
    //	    }
    // if (x < 0)
    //	        return 0.0;
    // 
    // if (shape == 0) /* point mass at 0 */
    //	        return (x == 0)? Constants::Double::inf : 0.0;
    //
    // if (x == 0) {
    //	        if (shape < 1) return Constants::Double::inf;
    //	        if (shape > 1) return 0.0;
    //	        /* else */
    //	        return isLog ? -log(scale) : 1 / scale;
    //	    }
    //
    // if (shape < 1) {
    //	        pr = Statistics::Poisson::pdf(shape, x/scale, isLog);
    //	        return isLog ?  pr + log(shape/x) : pr*shape/x;
    //	    }
    // /* else  shape >= 1 */
    // pr = Statistics::Poisson::pdf(shape-1, x/scale, isLog);
    // return isLog ? pr - log(scale) : pr/scale;
    
    return isLog ? pdf(shape, rate, exp(x)) : pdf(shape, rate, x);
}

/*!
 * This function calculates the natural log of the probability density 
 * for a gamma-distributed random variable.
 *
 * \brief Natural log of gamma probability density.
 * \param the shape parameter of the gamma. 
 * \param the rate parameter of the gamma. 
 * \param x is the gamma random variable. 
 * \return Returns the natural log of the probability density.
 * \throws Does not throw an error.
 */
double Statistics::Gamma::lnPdf(double shape, double rate, double x) {

	return shape * log(rate) - Math::lnGamma(shape) + (shape - 1.0) * log(x) - x * rate;
}

/*!
 * This function calculates the cumulative probability  
 * for a gamma-distributed random variable.
 *
 * \brief Gamma cumulative probability.
 * \param the shape parameter of the gamma. 
 * \param the rate parameter of the gamma. 
 * \param x is the gamma random variable. 
 * \return Returns the cumulative probability.
 * \throws Does not throw an error.
 */
double Statistics::Gamma::cdf(double shape, double rate, double x) {
    
	return Math::incompleteGamma( rate*x, shape, Math::lnGamma(shape) );
}

/*!
 * This function returns the quantile of a gamma probability 
 * distribution.
 *
 * \brief Gamma quantile.
 * \param the shape parameter. 
 * \param the rate parameter. 
 * \param p is the probability up to the quantile. 
 * \return Returns the quantile.
 * \throws Does not throw an error.
 */
double Statistics::Gamma::quantile(double shape, double rate, double p) {

	double ret = Constants::Double::nan;
	try{
		ret = Statistics::Helper::pointChi2(p, 2.0 * shape) / (2.0 * rate);
	}
	catch(...)
	{
		ret = Constants::Double::nan;
	}
	return ret <= 0 ? Constants::Double::nan : ret;
}


double Statistics::Gamma::rv(double shape, double rate, RandomNumberGenerator& rng) {
    
	return (Statistics::Helper::rndGamma(shape, rng) / rate);
}


/**
 * @file DistributionDirichlet
 * This file contains the functions of the dirichlet distribution.
 *
 * @brief Implementation of the dirichlet distribution.
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
#include "DistributionDirichlet.h"

#include "Exception.h"
#include "MathFunctions.h"
#include "StatisticsHelper.h"

/*!
 * This function calculates the probability density 
 * for a Dirichlet-distributed random variable.
 *
 * \brief Dirichlet probability density.
 * \param a is a reference to a vector of doubles containing the Dirichlet parameters. 
 * \param z is a reference to a vector of doubles containing the random variables. 
 * \return Returns the probability density.
 * \throws Throws an Exception::ERROR.
 */
double Statistics::Dirichlet::pdf(const std::vector<double> &a, const std::vector<double> &z) {
	
	size_t n = a.size();
	double zSum = 0.0;
	for (size_t i=0; i<n; i++)
		zSum += z[i];
    
	double tol = 0.0001;
	if ( tol < fabs( zSum - 1.0 ) )
        {
        std::ostringstream s;
        s << "Fatal error in Dirichlet PDF";
        throw Exception(s);
        }
    
	double aSum = 0.0;
	for (size_t i=0; i<n; i++)
		aSum += a[i];
    
	double aProd = 1.0;
	for (size_t i=0; i<n; i++)
		aProd *= Math::gamma(a[i]);
    
	double pdf = Math::gamma(aSum) / aProd;
    
	for (size_t i=0; i<n; i++)
		pdf = pdf * pow( z[i], a[i] - 1.0 );
    
	return pdf;
}



/*!
 * This function calculates the natural log of the probability density 
 * for a Dirichlet-distributed random variable.
 *
 * \brief Natural log of Dirichlet probability density.
 * \param a is a reference to a vector of doubles containing the Dirichlet parameters. 
 * \param z is a reference to a vector of doubles containing the random variables. 
 * \return Returns the natural log of the probability density.
 * \throws Does not throw an error.
 */
double Statistics::Dirichlet::lnPdf(const std::vector<double> &a, const std::vector<double> &z) {
    
	size_t n = a.size(); //!< we assume that a and z have the same size
	double alpha0 = 0.0;
	for (size_t i=0; i<n; i++)
		alpha0 += a[i];
	double lnP = Math::lnGamma(alpha0);
	for (size_t i=0; i<n; i++)
		lnP -= Math::lnGamma(a[i]);
	for (size_t i=0; i<n; i++)
		lnP += (a[i] - 1.0) * std::log(z[i]);	
	return lnP;
}



/*!
 * This function generates a Dirichlet-distributed random variable.
 *
 * \brief Dirichlet random variable.
 * \param a is a reference to a vector of doubles containing the Dirichlet parameters. 
 * \param rng is a pointer to a random number object. 
 * \return Returns a vector containing the Dirichlet random variable.
 * \throws Does not throw an error.
 */
std::vector<double> Statistics::Dirichlet::rv(const std::vector<double> &a, RandomNumberGenerator& rng) {
    
	size_t n = a.size();
    std::vector<double> z(n);
	double sum = 0.0;
	for(size_t i=0; i<n; i++)
        {
		z[i] = Statistics::Helper::rndGamma(a[i], rng);
		sum += z[i];
        }
	for(size_t i=0; i<n; i++)
		z[i] /= sum;
    
	return z;
}


/**
 * @file DistributionNegativeBinomial
 * This file contains the functions of the binomial distribution.
 *
 * @brief Implementation of the binomial distribution.
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

#include "DistributionBeta.h"
#include "DistributionBinomial.h"
#include "DistributionGamma.h"
#include "DistributionPoisson.h"
#include "DistributionNegativeBinomial.h"

#include "Constants.h"
#include "Exception.h"
#include "MathFunctions.h"
#include "MathHelper.h"
#include "MathLogic.h"

/*!
 * This function calculates the probability density 
 * for a binomially-distributed random variable.
 *
 * \brief NegativeBinomial probability density.
 * \param r is the number of failures.
 * \param p is the success probability. 
 * \param x is the number of successes. 
 * \return Returns the probability density.
 * \throws Does not throw an error.
 */
double Statistics::NegativeBinomial::cdf(double r, double p, double x) {

    if(Math::isInt(r))
        {
        std::ostringstream s;
        s << "Cannot compute cdf of the binomial distribution because r = " << r << " is not an interger";
        throw (Exception(s));
        }
    r = int(r);
    /* PR#8560: r=1 is a valid value */
    if(r < 1 || p < 0 || p > 1)
        {
        std::ostringstream s;
        s << "Cannot compute cdf of the binomial distribution for r = " << r << " and p = " << p;
        throw (Exception(s));
        }
    
    if (x < 0) 
        return 0.0;
    
    x = floor(x + 1e-7);
    
    return 1.0 - Statistics::Beta::cdf(p, x + 1, r);
}

/*!
 * This function calculates the probability density 
 * for a binomially-distributed random variable.
 *
 * \brief NegativeBinomial probability density.
 * \param r is the number of failures.
 * \param p is the success probability. 
 * \param x is the number of successes. 
 * \return Returns the probability density.
 * \throws Does not throw an error.
 */
double Statistics::NegativeBinomial::lnPdf(double n, double p, double x) {

    double q = 1.0 - p;
    return pdf(n, p, q, x, true);
}

/*!
 * This function calculates the probability density 
 * for a binomially-distributed random variable.
 *
 * \brief NegativeBinomial probability density.
 * \param r is the number of failures.
 * \param p is the success probability. 
 * \param x is the number of successes. 
 * \return Returns the probability density.
 * \throws Does not throw an error.
 */
double Statistics::NegativeBinomial::pdf(double n, double p, double x) {

    double q = 1.0 - p;
    return pdf(n,p,q,x,false);
}

/*!
 * This function calculates the probability density 
 * for a binomially-distributed random variable.
 *
 * From R:
 *
 *     (1) pdf() has both p and q arguments, when one may be represented
 *         more accurately than the other (in particular, in df()).
 *     (2) pdf() does NOT check that inputs x and n are integers. This
 *         should be done in the calling function, where necessary.
 *         -- but is not the case at all when called e.g., from df() or dbeta() !
 *     (3) Also does not check for 0 <= p <= 1 and 0 <= q <= 1 or NaN's.
 *         Do this in the calling function.
 *
 * \brief NegativeBinomial probability density.
 * \param r is the number of failures.
 * \param p is the success probability. 
 * \param x is the number of successes. 
 * \return Returns the probability density.
 * \throws Does not throw an error.
 */
double Statistics::NegativeBinomial::pdf(double r, double p, double q, double x, bool asLog) {

	double binom = Statistics::Binomial::pdf(x + r - 1, p , q , x, asLog);

    return (asLog ? binom + log(q) : binom*q);
}

/*!
 * This function calculates the probability density 
 * for a binomially-distributed random variable.
 *
 * From R:
 *
 *     (1) pdf() has both p and q arguments, when one may be represented
 *         more accurately than the other (in particular, in df()).
 *     (2) pdf() does NOT check that inputs x and n are integers. This
 *         should be done in the calling function, where necessary.
 *         -- but is not the case at all when called e.g., from df() or dbeta() !
 *     (3) Also does not check for 0 <= p <= 1 and 0 <= q <= 1 or NaN's.
 *         Do this in the calling function.
 *
 * \brief NegativeBinomial probability density.
 * \param r is the number of failures.
 * \param p is the success probability. 
 * \param x is the number of successes. 
 * \return Returns the probability density.
 * \throws Does not throw an error.
 */
#include "DistributionNormal.h"

double Statistics::NegativeBinomial::quantile(double p, double n, double pr) {
    throw Exception("Quantiles not implemented for NegBinomial");
}

/*!
 * This function draws a random variable  
 * from a binomial distribution.
 *
 * From R:
 *
 *     (1) pdf() has both p and q arguments, when one may be represented
 *         more accurately than the other (in particular, in df()).
 *     (2) pdf() does NOT check that inputs x and n are integers. This
 *         should be done in the calling function, where necessary.
 *         -- but is not the case at all when called e.g., from df() or dbeta() !
 *     (3) Also does not check for 0 <= p <= 1 and 0 <= q <= 1 or NaN's.
 *         Do this in the calling function.
 *
 *  REFERENCE
 *
 *	Kachitvichyanukul, V. and Schmeiser, B. W. (1988).
 *	NegativeBinomial random variate generation.
 *	Communications of the ACM 31, 216-222.
 *	(Algorithm BTPEC).
 *
 *
 * \brief NegativeBinomial probability density.
 * \param r is the number of failures.
 * \param p is the success probability. 
 * \param x is the number of successes. 
 * \return Returns the probability density.
 * \throws Does not throw an error.
 */

#include <stdlib.h>
#include <limits.h>
#include "RandomNumberGenerator.h"

#define repeat for(;;)

int Statistics::NegativeBinomial::rv(double nin, double pp, RandomNumberGenerator &rng)
{
	double y = Statistics::Gamma::rv(nin,(1-pp)/pp,rng);
	if(y > 2.0e9)
		return Constants::Integer::inf;
	return Statistics::Poisson::rv(y,rng);
}



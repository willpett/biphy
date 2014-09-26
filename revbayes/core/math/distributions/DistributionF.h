/**
 * @file DistributionF
 * This file contains the functions of the F distribution.
 *
 * @brief Implementation of the F distribution.
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


#ifndef DistributionF_H
#define DistributionF_H

namespace RevBayesCore {
    
    class RandomNumberGenerator;

    namespace RbStatistics {
    
        namespace F {
        
            double                      pdf(double m, double n, double x);                                    /*!< F(a,b) probability density */
            double                      pdf(double m, double n, double x, bool give_log);                     /*!< F(a,b) probability density */
            double                      lnPdf(double a, double b, double x);                                  /*!< F(a,b) log_e probability density */
            double                      cdf(double df1, double df2, double x);                                /*!< F(a,b) cumulative probability */
            double                      quantile(double a, double b, double p);                               /*!< F(a,b) quantile */
            double                      rv(double a, double b, RandomNumberGenerator& rng);                   /*!< F(a,b) random variable */
        }
    }
}

#endif

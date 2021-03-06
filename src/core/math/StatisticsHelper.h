/**
 * @file
 * This file contains commonly used statistics functions that are used
 * in RevBayes. The probability density (pdf), log of the probability
 * density (lnPdf), cumulative probability (cdf), and quantiles of
 * common probability distributions.
 *
 * @brief Namespace containing statistical functions
 *
 * (c) Copyright 2009- under GPL version 3
 * @date Last modified: $Date$
 * @author The RevBayes Development Core Team
 * @license GPL version 3
 *
 * $Id$
 */
#ifndef StatisticsHelper_H
#define StatisticsHelper_H

#include <vector>

#include "Exception.h"
#include "RandomNumberGenerator.h"

namespace Statistics {

    namespace Helper {
    
        double                      dppConcParamFromNumTables(double tables, double num);
        double                      dppExpectNumTableFromConcParam(double conp, double num);
        double                      pointChi2(double prob, double v);
        int                         poissonInver(double lambda, RandomNumberGenerator& rng);
        int                         poissonLow(double lambda, RandomNumberGenerator& rng);
        int                         poissonRatioUniforms(double lambda, RandomNumberGenerator& rng);
        double                      rndGamma(double s, RandomNumberGenerator& rng);
        double                      rndGamma1(double s, RandomNumberGenerator& rng);
        double                      rndGamma2(double s, RandomNumberGenerator& rng);
    
        template <class T> void		randomlySelectFromVectorWithReplacement(std::vector<T>& sourceV, std::vector<T>& destV, size_t k, RandomNumberGenerator& rng) {
        
            destV.clear();
            for (int i=0; i<k; i++)
                destV.push_back( sourceV[(int)(rng.uniform01()*(sourceV.size()))] );
        }
        template <class T> void		randomlySelectFromVectorWithoutReplacement(std::vector<T>& sourceV, std::vector<T>& destV, size_t k, RandomNumberGenerator& rng) {
        
            if ( (size_t)sourceV.size() < k ){
                std::cerr << "sourceV\n";
                throw (Exception("Attempting to sample too many elements from source vector"));
            }
            destV.clear();
            std::vector<T> tmpV = sourceV;
            size_t n = tmpV.size();
            for (size_t i=0; i<k; i++) {
                int whichElement = (int)(rng.uniform01()*(n-i));
                destV.push_back( tmpV[whichElement] );
                tmpV[whichElement] = tmpV[n-i-1];
            }
        }
        template <class T> void		permuteVector(std::vector<T>& v, RandomNumberGenerator* rng) {
        
            Statistics::Helper::randomlySelectFromVectorWithoutReplacement<T>(v, v, v.size(), rng);
        }
    }
    
}

#endif

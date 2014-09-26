//
//  TraceAnalysisContinuous.cpp
//  RevBayesGui
//
//  Created by Sebastian Hoehna on 4/4/11.
//  Copyright 2011 __MyCompanyName__. All rights reserved.
//

#include "TraceAnalysisContinuous.h"
#include "DistributionNormal.h"
#include "RbConstants.h"

#include <cmath>
#include <math.h>

using namespace RevBayesCore;


/**
 * 
 */
TraceAnalysisContinuous::TraceAnalysisContinuous() {
    // initialize data to invalid values
    invalidateTraceStatistics();
//    burnin               = new int;
    burnin              = RbConstants::Integer::nan;
    

}

void TraceAnalysisContinuous::analyseMean(const std::vector<double>& values) {
    // make sure the burnin is valid
    int b = 0;
    
    analyseMean(values, b);
}

void TraceAnalysisContinuous::analyseMean(const std::vector<double>& values, int b) {
    double m = 0;
    int size = (int)values.size();
    for (int i=b; i<size; i++) {
        m += values.at(i);
    }
    
    mean = m/(size-b);
}

void TraceAnalysisContinuous::analyseMean(const std::vector<std::vector<double> >& values, const std::vector<int>& burnins) {
    double m = 0;
    int sampleSize = 0;
    // get the number of chains
    int chains = (int)values.size();
    // iterate over all chains
    for (int j=0; j<chains; j++) {
        // get the chain
        const std::vector<double>& chain = values.at(j);
        int b = burnins.at(j);
        int chainSize = (int) chain.size() - b;
        // add this chain size to the total sample size
        sampleSize += chainSize;
        for (int i=burnin; i<chainSize; i++) {
            m += chain.at(i);
        }
    }
    
    mean = m/sampleSize;
}

void TraceAnalysisContinuous::analyseMean(const std::vector<double>& values, int begin, int end) {
    double m = 0;
    for (int i=begin; i<end; i++) {
        m += values.at(i);
    }
    
    mean = m/(end-begin);
}

/**
 * Analyze trace
 *
 */
void TraceAnalysisContinuous::analyseCorrelation(const std::vector<double>& values) {
    // make sure the burnin is valid
    int b = 0;
    
    analyseCorrelation(values,b);
}

/**
 * Analyze trace
 *
 * @param burnin   the number of sampes to discard
 *
 * @attention This method assumes that the mean was either made invalid before execution or is calculated apropriately for this burnin.
 */
void TraceAnalysisContinuous::analyseCorrelation(const std::vector<double>& values, int b) {
    // if we have not yet calculated the mean, do this now
    if (mean == RbConstants::Double::max) {
        analyseMean(values,b);
    }
    
    int samples = (int) values.size() - burnin;
    int maxLag = (samples - 1 < MAX_LAG ? samples - 1 : MAX_LAG);
    
    double* gammaStat = new double[maxLag];
    // setting values to 0
    for (int i=0; i<maxLag; i++) {
        gammaStat[i] = 0;
    }
    double varStat = 0.0;
    
    for (int lag = 0; lag < maxLag; lag++) {
        for (int j = 0; j < samples - lag; j++) {
            double del1 = values.at(burnin + j) - mean;
            double del2 = values.at(burnin + j + lag) - mean;
            gammaStat[lag] += (del1 * del2);
        }
        
        gammaStat[lag] /= ((double) (samples - lag));
        
        if (lag == 0) {
            varStat = gammaStat[0];
        } else if (lag % 2 == 0) {
            // fancy stopping criterion :)
            if (gammaStat[lag - 1] + gammaStat[lag] > 0) {
                varStat += 2.0 * (gammaStat[lag - 1] + gammaStat[lag]);
            }
            // stop
            else
                maxLag = lag;
        }
    }
    
    // standard error of mean
    sem = sqrt(varStat / samples);
    
    // auto correlation time
    act = varStat / gammaStat[0];
    
    // effective sample size
    ess = samples / act;
    
    // standard deviation of autocorrelation time
//    seAct = (2.0 * sqrt(2.0 * (2.0 * (double) (maxLag + 1)) / samples) * (varStat / gammaStat[0]) * stepSize);
    
    delete[] gammaStat;
}

/**
 * Analyze trace
 *
 * @param begin     begin index for analysis
 * @param end       end index for analysis
 *
 * @attention This method assumes that the mean was either made invalid before execution or is calculated apropriately for this burnin.
 */
void TraceAnalysisContinuous::analyseCorrelation(const std::vector<double>& values, int begin, int end) {
    // if we have not yet calculated the mean, do this now
    if (mean == RbConstants::Double::max) {
        analyseMean(values,burnin);
    }
    
    int samples = end - begin;
    int maxLag = (samples - 1 < MAX_LAG ? samples - 1 : MAX_LAG);
    
    double* gammaStat = new double[maxLag];
    // setting values to 0
    for (int i=0; i<maxLag; i++) {
        gammaStat[i] = 0;
    }
    double varStat = 0.0;
    
    for (int lag = 0; lag < maxLag; lag++) {
        for (int j = 0; j < samples - lag; j++) {
            double del1 = values.at(begin + j) - mean;
            double del2 = values.at(begin + j + lag) - mean;
            gammaStat[lag] += (del1 * del2);
        }
        
        gammaStat[lag] /= ((double) (samples - lag));
        
        if (lag == 0) {
            varStat = gammaStat[0];
        } else if (lag % 2 == 0) {
            // fancy stopping criterion :)
            if (gammaStat[lag - 1] + gammaStat[lag] > 0) {
                varStat += 2.0 * (gammaStat[lag - 1] + gammaStat[lag]);
            }
            // stop
            else
                maxLag = lag;
        }
    }
    
    // standard error of mean
    sem = sqrt(varStat / samples);
    
    // auto correlation time
    act = varStat / gammaStat[0];
    
    // effective sample size
    ess = samples / act;
    
    // standard deviation of autocorrelation time
    //    seAct = (2.0 * sqrt(2.0 * (2.0 * (double) (maxLag + 1)) / samples) * (varStat / gammaStat[0]) * stepSize);
    
    delete[] gammaStat;
}

void TraceAnalysisContinuous::invalidateTraceStatistics() {
    
    // initialize data to invalid values
    act                 = RbConstants::Double::max;
    ess                 = RbConstants::Double::max;
    mean                = RbConstants::Double::max;
    sem                 = RbConstants::Double::max;
}


/**
 * @return the mean
 */
double TraceAnalysisContinuous::getMean() {
    if (mean == RbConstants::Double::max) {
        // throw an error that the ACT needs recalculation
        //analyseMean(values);
    }
    
    return mean;
}

/**
 * @return the standard error of the mean
 */
double TraceAnalysisContinuous::getStdErrorOfMean() {
    if (sem == RbConstants::Double::max) {
        // throw an error that the ACT needs recalculation
        //analyseCorrelation(values);
    }
    
    return sem;
}

/**
 * @return the autocorrelation time
 */
double TraceAnalysisContinuous::getAct() {
    if (act == RbConstants::Double::max) {
        // throw an error that the ACT needs recalculation
        //analyseCorrelation(values);
    }
    
    return act;
}

/**
 * @return the eSS
 */
double TraceAnalysisContinuous::getEss() {
    if (ess == RbConstants::Double::max) {
        // throw an error that the ACT needs recalculation
//        analyseCorrelation(values);
    }
    
    return ess;
}


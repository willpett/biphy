//
//  GeneralRateMatrix.cpp
//  RevBayesCore
//
//  Created by Sebastian Hoehna on 4/25/13.
//  Copyright 2013 __MyCompanyName__. All rights reserved.
//

#include "GeneralRateMatrix.h"

using namespace RevBayesCore;


GeneralRateMatrix::GeneralRateMatrix(size_t n) : RateMatrix(n), stationaryFreqs( std::vector<double>(numStates,1.0/n) ), transitionRates( std::vector<double>(numStates*numStates-numStates, 1.0/n) ) {
    
}


GeneralRateMatrix::GeneralRateMatrix(const GeneralRateMatrix &rm) : RateMatrix(rm) {
    
    stationaryFreqs   = rm.stationaryFreqs;
    transitionRates   = rm.transitionRates;
}

GeneralRateMatrix::~GeneralRateMatrix(void) {
    // nothing to do
}


GeneralRateMatrix& GeneralRateMatrix::operator=(const GeneralRateMatrix &rm) {
    
    if (this != &rm)
    {
        RateMatrix::operator=( rm );
        
        stationaryFreqs   = rm.stationaryFreqs;
        transitionRates   = rm.transitionRates;
    }
    
    return *this;
}

/** Set the exchangeability rates directly. We assume that we know
 what the exchangeability rates are when this function is called. */
void GeneralRateMatrix::setTransitionRates(const std::vector<double>& tr) {
    
    transitionRates = tr;
    
    // set flags
    needsUpdate = true;
}


/** Set the stationary frequencies directly. We assume that we know
 what the stationary frequencies are when this function is called. */
void GeneralRateMatrix::setStationaryFrequencies(const std::vector<double>& f) {
    
    stationaryFreqs = f;
    
    // set flags
    needsUpdate = true;
}



const std::vector<double>& GeneralRateMatrix::getStationaryFrequencies( void ) const {
    return stationaryFreqs;
}


void GeneralRateMatrix::updateMatrix( void ) {
    
    if ( needsUpdate ) 
    {
        
        // rescale 
        rescaleToAverageRate( 1.0 );
        
        // now update the eigensystem
//        updateEigenSystem();
        
        // clean flags
        needsUpdate = false;
    }
}


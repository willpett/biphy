//
//  SingleRandomMoveSchedule.cpp
//  revbayes_mlandis
//
//  Created by Michael Landis on 2/6/13.
//  Copyright (c) 2013 Michael Landis. All rights reserved.
//

#include "SingleRandomMoveSchedule.h"
#include "RandomNumberFactory.h"
#include "RandomNumberGenerator.h"

using namespace RevBayesCore;

SingleRandomMoveSchedule::SingleRandomMoveSchedule(const std::vector<Move*> &s) : MoveSchedule( s ) {
    
    sumOfWeights = 0.0;
    for (std::vector<Move*>::const_iterator it = moves.begin(); it != moves.end(); ++it) {
        sumOfWeights+= (*it)->getUpdateWeight();
        weights.push_back( (*it)->getUpdateWeight() );
    }
}


SingleRandomMoveSchedule::~SingleRandomMoveSchedule() {
    // we own nothing
}


SingleRandomMoveSchedule* SingleRandomMoveSchedule::clone( void ) const {
    return new SingleRandomMoveSchedule(*this);
}

double SingleRandomMoveSchedule::getNumberMovesPerIteration( void ) const {
    return 1.0;
}


Move* SingleRandomMoveSchedule::nextMove( unsigned long gen ) {
    
    RandomNumberGenerator* rng = GLOBAL_RNG;
    double u = sumOfWeights * rng->uniform01();
    
    int index = 0;
    while ( weights[index] < u || !moves[index]->isActive( gen ) ) {
        u -= weights[index];
        ++index;
    }
    
    return moves[index];
}

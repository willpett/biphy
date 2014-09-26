//
//  SequenctialMoveSchedule.cpp
//  RevBayesCore
//
//  Created by Sebastian Hoehna on 8/17/12.
//  Copyright 2012 __MyCompanyName__. All rights reserved.
//

#include "SequenctialMoveSchedule.h"
#include "RandomNumberFactory.h"
#include "RandomNumberGenerator.h"

using namespace RevBayesCore;

SequentialMoveSchedule::SequentialMoveSchedule(const std::vector<Move*> &s) : MoveSchedule( s ), currentMove( 0 ), usedPropOfCurrentMove( 0 ) {
    
    movesPerIteration = 0.0;
    for (std::vector<Move*>::const_iterator it = moves.begin(); it != moves.end(); ++it) {
        movesPerIteration += (*it)->getUpdateWeight();
    }
}


SequentialMoveSchedule::~SequentialMoveSchedule() {
    // we own nothing
}


SequentialMoveSchedule* SequentialMoveSchedule::clone( void ) const {
    return new SequentialMoveSchedule(*this);
}

double SequentialMoveSchedule::getNumberMovesPerIteration( void ) const {
    return movesPerIteration;
}


Move* SequentialMoveSchedule::nextMove( unsigned long gen ) {
    
    bool found = false;
    do {
        double remainingWeight = moves[currentMove]->getUpdateWeight() - usedPropOfCurrentMove;
        if ( remainingWeight >= 1.0 ) {
            usedPropOfCurrentMove++;
            found = true;
        }
        else if ( remainingWeight > 0.0 ){
            usedPropOfCurrentMove++;
            RandomNumberGenerator* rng = GLOBAL_RNG;
            found = remainingWeight > rng->uniform01();
        } 
        else {
            usedPropOfCurrentMove = 0.0;
            do {
                currentMove++;
                if ( currentMove >= moves.size()) {
                    currentMove = 0;
                }
            } while ( moves[currentMove]->isActive( gen ) );
        }
    } while ( !found );
    
    return moves[currentMove];
}

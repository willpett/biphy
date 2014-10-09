/**
 * @file
 * This file contains the implementation of a move on the DPP based on Neal (2000) Algorithm 8
 * this move changes the value assigned to each table
 * works on doubles
 *
 * @brief Implementation of DPPBetaSimplexMove
 *
 * (c) Copyright 2009- under GPL version 3
 * @date Last modified: $Date: 2012-05-11 14:54:35 +0200 (Fri, 11 May 2012) $
 * @author The RevBayes Development Core Team
 * @license GPL version 3
 * @version 1.0
 * @since 2012-08-14, version 1.0
 *
 * $Id$
 */

#include "DPPBetaSimplexMove.h"
#include "DistributionBeta.h"
#include "RandomNumberFactory.h"
#include "RandomNumberGenerator.h"
#include "DirichletProcessPriorDistribution.h"
#include "RbStatisticsHelper.h"

#include <cmath>

RevBayesCore::DPPBetaSimplexMove::DPPBetaSimplexMove(StochasticNode<std::vector<double> > *v, double a, double w) : Move( v, w, false ), variable( v ) {
    
	// set isGibbs to true
	alpha = a;
}


/** Clone object */
RevBayesCore::DPPBetaSimplexMove* RevBayesCore::DPPBetaSimplexMove::clone( void ) const {
    
    return new DPPBetaSimplexMove( *this );
}



const std::string& RevBayesCore::DPPBetaSimplexMove::getMoveName( void ) const {
    static std::string name = "DPP Beta Simplex Move (double)";
    
    return name;
}

bool RevBayesCore::DPPBetaSimplexMove::isGibbs( void ) const {
    
    return true;
}


/** Perform the move */
void RevBayesCore::DPPBetaSimplexMove::performGibbsMove( void ) {
    
    // Get random number generator    
    RandomNumberGenerator* rng     = GLOBAL_RNG;
	
	DirichletProcessPriorDistribution<double>& dist = static_cast<DirichletProcessPriorDistribution<double> &>( variable->getDistribution() );
	
	dist.createRestaurantVectors();
	int numTables = dist.getNumberOfCategories();
	int numElements = dist.getNumberOfElements();
	std::vector<double> tableVals = dist.getTableParameters();
	std::vector<int> allocVec = dist.getElementAllocationVec();
	std::vector<double>& elementVals = variable->getValue();
	TypedDistribution<double>* g0 = dist.getBaseDistribution();
	
	double storedValue;
	variable->touch();
		
	// get old lnL
	double oldLnl = getCurrentLnProbabilityForMove();

	// randomly draw a new index
	double i = floor(rng->uniform01()*double(numTables));
	double currentValue = tableVals[i];

	// first, we get the parameters of the Beta for the forward move
	double af = alpha * currentValue + 1.0;
	double bf = alpha * (1.0-currentValue) + 1.0;

	// then, we propose new values
	double new_value = RbStatistics::Beta::rv( af, bf, *rng );

	tableVals[i] = new_value;

	// and calculate the Dirichlet parameters for the (imagined) reverse move
	double ar = alpha * new_value + 1.0;
	double br = alpha * (1.0-new_value) + 1.0;

	// finally, we calculate the log of the Hastings ratio
	double lnProposalRatio = RbStatistics::Beta::lnPdf(ar, br, currentValue) - RbStatistics::Beta::lnPdf(af, bf, new_value);

	// Assign new value to elements
	for(int j=0; j<numElements; j++){
		if(allocVec[j] == i)
			elementVals[j] = tableVals[i];
	}

	g0->getValue() = tableVals[i]; // new
	double priorRatio = g0->computeLnProbability();
	g0->getValue() = storedValue; // old
	priorRatio -= g0->computeLnProbability();

	variable->touch();
	double newLnl = getCurrentLnProbabilityForMove();
	double lnProbRatio = newLnl - oldLnl;

	double r = priorRatio + lnProbRatio + lnProposalRatio;
	double u = log(rng->uniform01());
	if ( u < r ) //accept
		variable->keep();
	else{ // reject
		for(int j=0; j<numElements; j++){
			if(allocVec[j] == i)
				elementVals[j] = storedValue;
		}
		tableVals[i] = storedValue;
		variable->touch();
		variable->keep();
	}
		
    dist.createRestaurantVectors();
}



void RevBayesCore::DPPBetaSimplexMove::swapNode(DagNode *oldN, DagNode *newN) {
    // call the parent method
    Move::swapNode(oldN, newN);
    
    variable = static_cast<StochasticNode<std::vector<double> >* >( newN );
}

double RevBayesCore::DPPBetaSimplexMove::getCurrentLnProbabilityForMove(void) {
	
	std::set<DagNode*> affected;
	variable->getAffectedNodes( affected );
	double lnProb = 0.0;
	for (std::set<DagNode*>::iterator it = affected.begin(); it != affected.end(); ++it) {
		double lp = (*it)->getLnProbability();
		lnProb += lp;
	}
	return lnProb;
}


double RevBayesCore::DPPBetaSimplexMove::safeExponentiation(double x) {
	
	if (x < -300.0)
		return 0.0;
	else if (x > 0.0)
		return 1.0;
	else
		return exp(x);
}

void RevBayesCore::DPPBetaSimplexMove::acceptMove( void ) {
    
}

double RevBayesCore::DPPBetaSimplexMove::performMove(double& probRatio) {
   return 0.0;
}

void RevBayesCore::DPPBetaSimplexMove::rejectMove( void ) {
    
}




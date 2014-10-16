/**
 * @file
 * This file contains the declaration of DolloBirthRateGibbsMove, which performs the DPP move based on Neal (2000) Algorithm 8
 * this move changes the value assigned to each table
 *
 * @brief Declaration of DolloBirthRateGibbsMove
 *
 * (c) Copyright 2009- under GPL version 3
 * @date Last modified: $Date: 2012-05-11 14:54:35 +0200 (Fri, 11 May 2012) $
 * @author The RevBayes Development Core Team
 * @license GPL version 3
 * @version 1.0
 * @since 2012-07-20, version 1.0
 *
 * $Id: DolloBirthRateGibbsMove.h $
 */

#ifndef DolloBirthRateGibbsMove_H
#define DolloBirthRateGibbsMove_H

#include <ostream>

#include "DolloBranchHeterogeneousCharEvoModel.h"

#include "DistributionGamma.h"
#include "Move.h"
#include "StochasticNode.h"
#include "RandomNumberFactory.h"
#include "RandomNumberGenerator.h"
#include "RbStatisticsHelper.h"

#include <cmath>

namespace RevBayesCore {

	template<class charType, class treeType>
    class DolloBirthRateGibbsMove : public Move {
    
    public:
        DolloBirthRateGibbsMove(StochasticNode<double>* l, StochasticNode<AbstractCharacterData>* m, double w);                                                                      //!< Internal constructor
    
        // Basic utility functions
        DolloBirthRateGibbsMove*								clone(void) const;                                                                  //!< Clone object
        void                                                    swapNode(DagNode *oldN, DagNode *newN);
        const std::string&                                      getMoveName(void) const;                                                            //!< Get the name of the move for summary printing
		bool													isGibbs(void) const;

    protected:
        void													performGibbsMove(void);                                                            //!< Perform move
        void													acceptMove(void);                                                                   //!< Accept the InferenceMoveSimple
        double													performMove(double& probRatio);                                                     //!< Perform the InferenceMoveSimple
        void													rejectMove(void);                                                                   //!< Reject the InferenceMoveSimple
    
    private:
        StochasticNode<double>*                   				lambda;
        StochasticNode<AbstractCharacterData>* 					model;
 
    };
    
}

template<class charType, class treeType>
RevBayesCore::DolloBirthRateGibbsMove<charType,treeType>::DolloBirthRateGibbsMove(StochasticNode<double > *l, StochasticNode<AbstractCharacterData>* m, double w) : Move( l, w, false ), lambda( l ), model(m) {

	// set isGibbs to true
	nodes.insert(model);
}


/** Clone object */
template<class charType, class treeType>
RevBayesCore::DolloBirthRateGibbsMove<charType,treeType>* RevBayesCore::DolloBirthRateGibbsMove<charType,treeType>::clone( void ) const {

    return new DolloBirthRateGibbsMove( *this );
}


template<class charType, class treeType>
const std::string& RevBayesCore::DolloBirthRateGibbsMove<charType,treeType>::getMoveName( void ) const {
    static std::string name = "Dollo Birth Rate Gibbs Move";

    return name;
}

template<class charType, class treeType>
bool RevBayesCore::DolloBirthRateGibbsMove<charType,treeType>::isGibbs( void ) const {

    return true;
}

template<class charType, class treeType>
void RevBayesCore::DolloBirthRateGibbsMove<charType,treeType>::performGibbsMove( void ) {

    // Get random number generator
    RandomNumberGenerator* rng     = GLOBAL_RNG;

    DolloBranchHeterogeneousCharEvoModel<charType,treeType>& dist = static_cast<DolloBranchHeterogeneousCharEvoModel<charType,treeType> &>( model->getDistribution() );

    double alpha = dist.getNumSites()+1;
    double beta = dist.getTotalMass();

    double newval = RbStatistics::Gamma::rv(alpha, beta, *rng);
    lambda->setValue(newval);
    lambda->keep();
}

template<class charType, class treeType>
void RevBayesCore::DolloBirthRateGibbsMove<charType,treeType>::swapNode(DagNode *oldN, DagNode *newN) {
    // call the parent method
    Move::swapNode(oldN, newN);

    if(oldN == model){
        std::cerr << "model" << std::endl;
    	model = static_cast<StochasticNode<AbstractCharacterData>* >( newN );
    }else if(oldN == lambda){
    	lambda = static_cast<StochasticNode<double>* >( newN );
    }
}

template<class charType, class treeType>
void RevBayesCore::DolloBirthRateGibbsMove<charType,treeType>::acceptMove( void ) {

}

template<class charType, class treeType>
double RevBayesCore::DolloBirthRateGibbsMove<charType,treeType>::performMove(double& probRatio) {
   return 0.0;
}

template<class charType, class treeType>
void RevBayesCore::DolloBirthRateGibbsMove<charType,treeType>::rejectMove( void ) {

}



#endif


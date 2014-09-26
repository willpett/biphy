/**
 * @file
 * This file contains the declaration of a complete path scaling move, which changes
 * a all values of the path of a random walk (Brownian motion) by scaling them.
 *
 * @brief Declaration of CompletePathScalingMove
 *
 * (c) Copyright 2009- under GPL version 3
 * @date Last modified: $Date: 2012-05-11 14:54:35 +0200 (Fri, 11 May 2012) $
 * @author The RevBayes Development Core Team
 * @license GPL version 3
 * @version 1.0
 * @since 2012-08-14, version 1.0
 *
 * $Id: SimplexMove.h 1535 2012-05-11 12:54:35Z hoehna $
 */

#ifndef CompletePathRescaleMove_H
#define CompletePathRescaleMove_H

#include <ostream>
#include <set>
#include <string>

#include "UnivariateFunction.h"
#include "StochasticNode.h"
#include "SimpleMove.h"

namespace RevBayesCore {
    
    class CompletePathRescaleMove : public SimpleMove {
        
    public:
        CompletePathRescaleMove(StochasticNode<UnivariateFunction>* node, double lambda, bool tuning, double weight);                                         //!< Internal constructor
        
        // Basic utility functions
        CompletePathRescaleMove*                clone(void) const;                                                                  //!< Clone object
        void                                    swapNode(DagNode *oldN, DagNode *newN);
        
    protected:
        const std::string&                      getMoveName(void) const;                                                            //!< Get the name of the move for summary printing
        double                                  performSimpleMove(void);                                                            //!< Perform move
        void                                    printParameterSummary(std::ostream &o) const;
        void                                    rejectSimpleMove(void);
        void                                    tune(void);
       
    private:
        
        StochasticNode<UnivariateFunction>*     variable;
        double                                  lambda;
        std::vector<double>                     storedValue;
        
    };
    
}

#endif


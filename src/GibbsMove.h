/**
 * @file
 * This file contains the declaration of Subtree-Prune-and-Regraft, 
 * which randomly draws a node in the tree and exchanges its two neighbors.
 *
 * @brief Declaration of GibbsMove
 *
 * (c) Copyright 2009- under GPL version 3
 * @date Last modified: $Date: 2012-07-12 16:14:14 +0200 (Thu, 12 Jul 2012) $
 * @author The RevBayes Development Core Team
 * @license GPL version 3
 * @version 1.0
 * @since 2012-07-12, version 1.0
 *
 * $Id: ScaleMove.h 1677 2012-07-12 14:14:14Z hoehna $
 */

#ifndef GibbsMove_H
#define GibbsMove_H

#include <ostream>
#include <set>
#include <string>

#include "Move.h"
#include "StochasticNode.h"
#include "Topology.h"

namespace RevBayesCore {
    
	template<class valueType>
    class GibbsMove : public Move {
        
    public:
        GibbsMove( StochasticNode<valueType > *n, double weight);                                            //!<  constructor
        
        // Basic utility functions
        GibbsMove*                  				clone(void) const;                                                                  //!< Clone object
        void                                        swapNode(DagNode *oldN, DagNode *newN);
        virtual bool                                isGibbs(void) const;
        const std::string&                          getMoveName(void) const;
        
    protected:
        virtual void                                performGibbsMove(void);                                                           //!< Get the name of the move for summary printing
        void										acceptMove(void);                                                                   //!< Accept the InferenceMoveSimple
		double										performMove(double& probRatio);                                                     //!< Perform the InferenceMoveSimple
		void										rejectMove(void);
        
    private:
        
        // member variables
        StochasticNode<valueType>*    				 variable;
        
    };

	template<class valueType>
	GibbsMove<valueType>::GibbsMove(StochasticNode<valueType> *v, double w) : Move( v, w), variable( v ) {

	}


	/* Clone object */
	template<class valueType>
	GibbsMove<valueType>* GibbsMove<valueType>::clone( void ) const {

	    return new GibbsMove( *this );
	}

	template<class valueType>
	bool GibbsMove<valueType>::isGibbs( void ) const {

	    return true;
	}

	template<class valueType>
	const std::string& GibbsMove<valueType>::getMoveName( void ) const {
	    static std::string name = "GibbsMove";

	    return name;
	}


	/** Perform the move */
	template<class valueType>
	void GibbsMove<valueType>::performGibbsMove( void ) {

	    variable->redraw();
	}

	template<class valueType>
	void GibbsMove<valueType>::acceptMove( void ) {

	}

	template<class valueType>
	double GibbsMove<valueType>::performMove(double& probRatio) {
		return 0.0;
	}


	template<class valueType>
	void GibbsMove<valueType>::rejectMove( void ) {

	}

	template<class valueType>
	void GibbsMove<valueType>::swapNode(DagNode *oldN, DagNode *newN) {
	    // call the parent method
	    Move::swapNode(oldN, newN);

	    variable = static_cast<StochasticNode<valueType >* >(newN) ;
	}
    
}

#endif


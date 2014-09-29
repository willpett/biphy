#ifndef MultiMove_H
#define MultiMove_H

#include <ostream>
#include <set>
#include <string>

#include "SimpleMove.h"
#include "StochasticNode.h"

namespace RevBayesCore {
    
    /**
     * The scaling operator. 
     *
     * A scaling move draws a random uniform number u ~ unif(-0.5,0.5)
     * and scales the current vale by a scaling factor
     * sf = exp( lambda * u )
     * where lambda is the tuning parameter of the move to influence the size of the proposals.
     *
     * @copyright Copyright 2009-
     * @author The RevBayes Development Core Team (Sebastian Hoehna)
     * @since 2009-09-08, version 1.0
     *
     */
    class MultiMove : public Move {

    public:
    	MultiMove(std::vector<Move*> n, DagNode* m, double w, bool autoTune = false);                                             //!<  constructor

        // Basic utility functions
    	MultiMove*              	clone(void) const;                                                                  //!< Clone object
    	void                        swapNode(DagNode *oldN, DagNode *newN);

    	std::vector<Move*>			getMoves(void) const;
    	const std::set<DagNode*>&   getDagNodes(void) const;

    protected:
        const std::string&              getMoveName(void) const;                                                            //!< Get the name of the move for summary printing
        double                          performMove(double &probRatio);                                                            //!< Perform move
        void                            acceptMove();
        void                            rejectMove(void);
        void                            tune(void);                                                                         //!< Tune the move to achieve a better acceptance/rejection ratio

    private:
        std::vector<Move*> 				moves;
        std::set<DagNode*> 		    	move_nodes;
    };
    
}

#endif


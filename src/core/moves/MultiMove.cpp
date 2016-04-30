#include "MultiMove.h"
#include "ConstantNode.h"
#include "RandomNumberFactory.h"
#include "RandomNumberGenerator.h"
#include <cmath>
#include <iostream>
#include "Exception.h"

/** 
 * Constructor
 *
 * Here we simply allocate and initialize the move object.
 */
MultiMove::MultiMove(std::vector<Move*> n, DagNode* m, double w, bool t) : Move(m, w, t), moves(n)
{
	for (std::vector<Move*>::iterator it = moves.begin(); it != moves.end(); it++){
			std::set<DagNode*> set = (*it)->getDagNodes();
			for(std::set<DagNode*>::iterator jt = set.begin(); jt != set.end(); jt++){
				move_nodes.insert( *jt );
				nodes.insert(*jt);
			}
	}
}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'B'.
 *
 * \return A new copy of the model. 
 */
MultiMove* MultiMove::clone( void ) const
{

    return new MultiMove( *this );
}


/**
 * Get moves' name of object 
 *
 * \return The moves' name.
 */
const std::string& MultiMove::getMoveName( void ) const
{
    static std::string name = "MultiMove";
    
    return name;
}

std::vector<Move*> MultiMove::getMoves( void ) const
{
    return moves;
}


/** 
 * Perform the move.
 *
 * A scaling move draws a random uniform number u ~ unif(-0.5,0.5)
 * and scales the current vale by a scaling factor
 * sf = exp( lambda * u )
 * where lambda is the tuning parameter of the move to influence the size of the proposals.
 *
 * \return The hastings ratio.
 */
double MultiMove::performMove( double &probRatio )
{
    double lnHastingsratio = 0.0;
    probRatio = 0.0;
    
    for(size_t i = 0; i < moves.size(); i++){
    	//std::cerr << moves[i]->getMoveName() << "\n";
    	double prob = 0.0;
		lnHastingsratio += moves[i]->perform(prob);
		//std::cerr << "\t" << prob << "\n";
		probRatio += prob;
    }
    //std::cerr << probRatio << "\t" << lnHastingsratio << "\n";
    return lnHastingsratio;
}


/**
 * Reject the move.
 *
 * Since the move stores the previous value and it is the only place
 * where complex undo operations are known/implement, we need to revert
 * the value of the variable/DAG-node to its original value.
 */
void MultiMove::rejectMove( void )
{
    // swap current value and stored value
	for(size_t i = 0; i < moves.size(); i++){
		moves[i]->reject();
	}
    
	//std::cerr << "rejected\n";
}

/**
 * Reject the move.
 *
 * Since the move stores the previous value and it is the only place
 * where complex undo operations are known/implement, we need to revert
 * the value of the variable/DAG-node to its original value.
 */
void MultiMove::acceptMove( void )
{
    // swap current value and stored value
	for(size_t i = 0; i < moves.size(); i++){
		moves[i]->accept();
	}
	//std::cerr << "accepted\n";
}


/**
 * Tune the move to accept the desired acceptance ratio.
 *
 * The acceptance ratio for this move should be around 0.44.
 * If it is too large, then we increase the proposal size,
 * and if it is too small, then we decrease the proposal size.
 */
void MultiMove::tune( void ) {
	for (std::vector<Move*>::iterator it = moves.begin(); it != moves.end(); it++)
	{
		(*it)->autoTune();
	}
}

void MultiMove::swapNode(DagNode *oldN, DagNode *newN) {
    // error catching
    if ( move_nodes.find(oldN) != move_nodes.end() ) {
    	for (std::vector<Move*>::iterator it = moves.begin(); it != moves.end(); it++){
    		std::set<DagNode*> set = (*it)->getDagNodes();
    		if(set.find(oldN) != set.end()){
    			(*it)->swapNode(oldN,newN);
    			//std::cerr << "swapped " << *it << "\n";
    		}
    	}
    }else{
    	Move::swapNode(oldN,newN);
    }
}


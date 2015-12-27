#ifndef MultiMove_H
#define MultiMove_H

#include <ostream>
#include <set>
#include <string>

#include "SimpleMove.h"
#include "StochasticNode.h"

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

#endif


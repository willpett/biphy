#ifndef ScaleMove_H
#define ScaleMove_H

#include <ostream>
#include <set>
#include <string>

#include "SimpleMove.h"
#include "StochasticNode.h"

class ScaleMove : public SimpleMove {

public:
    ScaleMove( StochasticNode<double> *n, double l, bool tuning, double w);                                             //!<  constructor

    // Basic utility functions
    ScaleMove*                      clone(void) const;                                                                  //!< Clone object
    void                            swapNode(DagNode *oldN, DagNode *newN);                                             //!< Swap the DAG nodes on which the move is working on

protected:
    const std::string&              getMoveName(void) const;                                                            //!< Get the name of the move for summary printing
    double                          performSimpleMove(void);                                                            //!< Perform move
    void                            printParameterSummary(std::ostream &o) const;                                       //!< Print the parameter summary
    void                            rejectSimpleMove(void);                                                             //!< Reject the move
    void                            tune(void);                                                                         //!< Tune the move to achieve a better acceptance/rejection ratio

private:
    // parameters
    double                          lambda;                                                                             //!< The scaling parameter of the move  

    StochasticNode<double>*         variable;                                                                           //!< The variable the move is working on
    double                          storedValue;                                                                        //!< The stored value of the move used for rejections.

};

#endif


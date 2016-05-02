#ifndef PartialVectorScaleMove_H
#define PartialVectorScaleMove_H

#include <ostream>
#include <set>
#include <string>

#include "Move.h"
#include "ContinuousStochasticNode.h"

class PartialVectorScaleMove : public Move {

public:
    PartialVectorScaleMove( StochasticNode<double>* p, std::vector<ContinuousStochasticNode *> br, double l, bool tuning, double w);                                             //!<  constructor

    // Basic utility functions
    PartialVectorScaleMove*                      clone(void) const;                                                                  //!< Clone object
    void                            swapNode(DagNode *oldN, DagNode *newN);                                             //!< Swap the DAG nodes on which the move is working on

protected:
    const std::string&              getMoveName(void) const;                                                            //!< Get the name of the move for summary printing
    double                          performMove(double& probRatio);                                                            //!< Perform move
    void                            printParameterSummary(std::ostream &o) const;                                       //!< Print the parameter summary
    void                            rejectMove(void);
    void                            acceptMove(void);

private:
    // parameters
    double                          lambda;                                                                             //!< The scaling parameter of the move  

    StochasticNode<double>*                   partial;
    double                                    storedpartial;
    std::vector<ContinuousStochasticNode* >   variables;
    std::vector<double>                       storedvariables;

};

#endif


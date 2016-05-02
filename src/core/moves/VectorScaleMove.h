#ifndef VectorScaleMove_H
#define VectorScaleMove_H

#include <ostream>
#include <set>
#include <string>

#include "Move.h"
#include "ContinuousStochasticNode.h"

class VectorScaleMove : public Move {

public:
    VectorScaleMove( std::vector<ContinuousStochasticNode *> br, double l, bool tuning, double w);                                             //!<  constructor

    // Basic utility functions
    VectorScaleMove*                      clone(void) const;                                                                  //!< Clone object
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

    std::vector<ContinuousStochasticNode* >   variables;
    std::vector<double>                     storedvariables;

};

#endif


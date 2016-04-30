#ifndef BranchLengthFrequencyCompensatoryMove_H
#define BranchLengthFrequencyCompensatoryMove_H

#include <ostream>
#include <set>
#include <string>

#include "SimpleMove.h"
#include "ContinuousStochasticNode.h"

class BranchLengthFrequencyCompensatoryMove : public SimpleMove {

public:
    BranchLengthFrequencyCompensatoryMove( StochasticNode<double> *pi, std::vector<ContinuousStochasticNode *> br, double l, bool tuning, double w);                                             //!<  constructor

    // Basic utility functions
    BranchLengthFrequencyCompensatoryMove*                      clone(void) const;                                                                  //!< Clone object
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

    StochasticNode<double>*                 frequency;
    std::vector<ContinuousStochasticNode* >   branchlengths;
    double                                  storedFrequency;
    std::vector<double>                     storedBranchlengths;

};

#endif


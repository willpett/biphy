#ifndef VectorSingleElementScaleMove_H
#define VectorSingleElementScaleMove_H

#include "SimpleMove.h"
#include "StochasticNode.h"

#include <ostream>
#include <vector>
#include <string>

class VectorSingleElementScaleMove : public SimpleMove {
    
public:
    VectorSingleElementScaleMove(StochasticNode<std::vector<double> >* n, double l, bool t, double w);                         //!< Constructor
    
    // Basic utility functions
    VectorSingleElementScaleMove*                     clone(void) const;                                                                  //!< Clone this object.
    const std::string&                          getMoveName(void) const;                                                            //!< Get the name of the move for summary printing.
    void                                        swapNode(DagNode *oldN, DagNode *newN);                                             //!< Swap the variable if it was replaced.
    
protected:
    
    double                                  performSimpleMove(void);                                                            //!< Perform move
    void                                    printParameterSummary(std::ostream &o) const;
    void                                    rejectSimpleMove(void);
    void                                    acceptSimpleMove(void);
    void                                    tune(void);
    void                                    touch( DagNode *toucher );
    
private:
    
    StochasticNode<std::vector<double> >*   variable;

    double                                      lambda;                                                                             //!< The scale parameter of the move (larger lambda -> larger proposals).
    size_t                                      index;                                                                              //!< The index of the last modified element.
    double                                      storedValue;                                                                        //!< The stored value of the last modified element.
    
};

#endif


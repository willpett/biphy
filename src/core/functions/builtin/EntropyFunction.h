/**
 * @file
 * This file contains the declaration of the deterministic variable class for Vectors.
 * This class is derived from the deterministic node Unot each instance will represent a deterministic variable
 * computing the Vector of its parameters.
 *
 * @brief Declaration of the deterministic variable for Vectors.
 *
 * (c) Copyright 2009- under GPL version 3
 * @date Last modified: $Date$
 * @author The RevBayes Development Core Team
 * @license GPL version 3
 * @version 1.0
 * @since 2012-07-06, version 1.0
 * @interface TypedDagNode
 *
 * $Id$
 */



#ifndef EntropyFunction_H
#define EntropyFunction_H

#include "TypedFunction.h"
#include "TypedDagNode.h"

#include <vector>

class EntropyFunction : public TypedFunction<double> {
    
public:
    EntropyFunction(const TypedDagNode<std::vector<double> > * v);
    EntropyFunction(const EntropyFunction &n);                                                                                        //!< Copy constructor
    virtual                                            ~EntropyFunction(void);                                                       //!< Virtual destructor
    
    // public member functions
    EntropyFunction*                                       clone(void) const;                                                          //!< Create an independent clone
    void                                                update(void);
    
protected:
    void                                                swapParameterInternal(const DagNode *oldP, const DagNode *newP);            //!< Implementation of swaping parameters
    
private:
    
    // members
    const TypedDagNode<std::vector<double> >*           vals;
    
};


#endif

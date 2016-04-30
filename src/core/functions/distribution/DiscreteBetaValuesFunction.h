/**
 * @file
 * This file contains the declaration of the deterministic variable class for Vectors.
 * This class is derived from the deterministic node and each instance will represent a deterministic variable
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



#ifndef DiscreteBetaValuesFunction_H
#define DiscreteBetaValuesFunction_H

#include "TypedFunction.h"
#include "TypedDagNode.h"

#include <vector>

class DiscreteBetaValuesFunction : public TypedFunction< std::vector<double> > {
    
public:
    DiscreteBetaValuesFunction(const TypedDagNode<double> *a, const TypedDagNode<double> *b, std::vector<double> breaks, bool mean = true);
    DiscreteBetaValuesFunction(const DiscreteBetaValuesFunction &n);                                                                                        //!< Copy constructor
    virtual                                            ~DiscreteBetaValuesFunction(void);                                                       //!< Virtual destructor
    
    // public member functions
    DiscreteBetaValuesFunction*                         clone(void) const;
    void                                                update(void);
    
protected:
    void                                                swapParameterInternal(const DagNode *oldP, const DagNode *newP);            //!< Implementation of swaping parameters
    
private:
    
    // members
    const TypedDagNode<double>* 						alpha;
    const TypedDagNode<double>* 						beta;

    std::vector<double>									breaks;

    bool												mean;
    
};

#endif

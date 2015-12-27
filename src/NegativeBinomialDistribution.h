/**
 * @file
 * This file contains the declaration of the binomially distribution class.
 *
 * @brief Declaration of the binomial distribution.
 *
 * (c) Copyright 2009- under GPL version 3
 * @date Last modified: $Date:$
 * @author The RevBayes Development Core Team
 * @license GPL version 3
 * @version 1.0
 * @since 2012-06-17, version 1.0
 * @interface TypedDagNode
 *
 * $Id:$
 */



#ifndef NegativeBinomialDistribution_H
#define NegativeBinomialDistribution_H

#include "TypedDagNode.h"
#include "TypedDistribution.h"

class NegativeBinomialDistribution : public TypedDistribution<int> {
    
public:
    NegativeBinomialDistribution(const TypedDagNode<int> *n, const TypedDagNode<double> *p);
    virtual                                            ~NegativeBinomialDistribution(void);                                             //!< Virtual destructor
    
    // public member functions
    NegativeBinomialDistribution*                               clone(void) const;                                                      //!< Create an independent clone
    double                                              computeLnProbability(void);
    void                                                redrawValue(void);

protected:
    // Parameter management functions
    void                                                swapParameter(const DagNode *oldP, const DagNode *newP);        //!< Swap a parameter
    
private:
    
    // members
    const TypedDagNode<int>*                            r;
    const TypedDagNode<double>*                         p;
};
#endif

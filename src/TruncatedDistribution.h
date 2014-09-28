/**
 * @file
 * This file contains the declaration of the uniformly distributed random variable class.
 * This class is derived from the stochastic node and each instance will represent a random variable
 * from a normal distribution in the model graph.
 *
 * @brief Declaration of the stochastic DAG node base class.
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



#ifndef TruncatedDistribution_H
#define TruncatedDistribution_H

#include "ContinuousDistribution.h"
#include "TypedDagNode.h"

namespace RevBayesCore {
    
	class TruncatedDistribution : public TypedDistribution<double> {
        
    public:
        TruncatedDistribution(TypedDistribution<double> *f, const TypedDagNode<double> *min, const TypedDagNode<double> *max);
        TruncatedDistribution(const TruncatedDistribution &n);                                                                              //!< Copy constructor
        virtual                                            ~TruncatedDistribution(void);                                                  //!< Virtual destructor
        
        // public member functions
        TruncatedDistribution*                              clone(void) const;                                                          //!< Create an independent clone
        double                                              computeLnProbability(void);
        double                                              getMax(void) const;
        double                                              getMin(void) const;                                                   //!< Qu
        void                                                redrawValue(void);
        void                                                swapParameter(const DagNode *oldP, const DagNode *newP);                    //!< Implementation of swaping parameters
        
    private:
        
        // members
        const TypedDagNode<double>*                          min;
        const TypedDagNode<double>*                          max;
        TypedDistribution<double>*				     		 f;
        
    };
    
}

#endif

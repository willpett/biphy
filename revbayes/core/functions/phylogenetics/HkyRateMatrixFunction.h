/**
 * @file
 * This file contains the declaration of the HKY rate matrix function class.
 * This class is derived from the function class and is used to
 * compute the rate matrix of a general time reversible (Hky) Markov chain.
 *
 * @brief Declaration of the HKY rate matrix function.
 *
 * (c) Copyright 2009- under GPL version 3
 * @date Last modified: $Date$
 * @author The RevBayes Development Core Team
 * @license GPL version 3
 * @version 1.0
 * @since 2012-07-06, version 1.0
 * @interface Function
 *
 * $Id$
 */



#ifndef HkyRateMatrixFunction_H
#define HkyRateMatrixFunction_H

#include "RateMatrix_HKY.h"
#include "TypedDagNode.h"
#include "TypedFunction.h"

#include <vector>

namespace RevBayesCore {
    
    class HkyRateMatrixFunction : public TypedFunction<RateMatrix> {
        
    public:
        HkyRateMatrixFunction(const TypedDagNode<double> *k, const TypedDagNode<std::vector<double> > *bf);
        HkyRateMatrixFunction(const HkyRateMatrixFunction &n);                                                                              //!< Copy constructor
        virtual                                            ~HkyRateMatrixFunction(void);                                                    //!< Virtual destructor
        
        // public member functions
        HkyRateMatrixFunction*                              clone(void) const;                                                              //!< Create an independent clone
        void                                                update(void);
        
    protected:
        void                                                swapParameterInternal(const DagNode *oldP, const DagNode *newP);                        //!< Implementation of swaping parameters
        
    private:
        
        // members
        const TypedDagNode<std::vector<double> >*           baseFrequencies;
        const TypedDagNode<double>*                         kappa;
        
    };
    
}

#endif

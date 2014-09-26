//
//  FreeBinaryRateMatrixFunction.h
//  rb_mlandis
//
//  Created by Michael Landis on 4/4/14.
//  Copyright (c) 2014 Michael Landis. All rights reserved.
//

#ifndef FreeBinaryRateMatrixVectorFunction_H
#define FreeBinaryRateMatrixVectorFunction_H

#include "RateMatrix_FreeBinary.h"
#include "RbVector.h"
#include "TypedDagNode.h"
#include "TypedFunction.h"

#include <vector>

namespace RevBayesCore {
    
    class FreeBinaryRateMatrixVectorFunction : public TypedFunction<RbVector<RateMatrix> > {
        
    public:
        FreeBinaryRateMatrixVectorFunction(const TypedDagNode<std::vector<double> > *bf);
        FreeBinaryRateMatrixVectorFunction(const FreeBinaryRateMatrixVectorFunction &n);                                                                              //!< Copy constructor
        virtual                                            ~FreeBinaryRateMatrixVectorFunction(void);                                                    //!< Virtual destructor
        
        // public member functions
        FreeBinaryRateMatrixVectorFunction*                 clone(void) const;                                                              //!< Create an independent clone
        void                                                update(void);
        
    protected:
        void                                                swapParameterInternal(const DagNode *oldP, const DagNode *newP);                        //!< Implementation of swaping parameters
        
    private:
        
        // members
        const TypedDagNode<std::vector<double> >*           transitionRates;
        
    };
    
}

#endif /* defined(FreeBinaryRateMatrixVectorFunction_H) */

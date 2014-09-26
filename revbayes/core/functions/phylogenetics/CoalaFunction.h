/**
 * @file
 * This file contains the declaration of the F81 rate matrix function class.
 * This class is derived from the function class and is used to
 * compute the rate matrix of a Felsenstein (F81) Markov chain.
 *
 * @brief Declaration of the F81 rate matrix function.
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



#ifndef CoalaFunction_H
#define CoalaFunction_H

#include "CorrespondenceAnalysis.h"
#include "TypedDagNode.h"
#include "TypedFunction.h"

#include <vector>

namespace RevBayesCore {
    
    class CoalaFunction : public TypedFunction< std::vector<double> > {
        
    public:
        CoalaFunction(const TypedDagNode<std::vector<double> > *coords, const MatrixReal &ca, const std::vector<double> & cw);
        CoalaFunction(const CoalaFunction &n);                                                                              //!< Copy constructor
        virtual                                            ~CoalaFunction(void);                                                    //!< Virtual destructor
        
        // public member functions
        CoalaFunction*                                      clone(void) const;                                                              //!< Create an independent clone
        void                                                update(void);
        
    protected:
        void                                                swapParameterInternal(const DagNode *oldP, const DagNode *newP);                        //!< Implementation of swaping parameters
        
    private:
        
        // members
        const TypedDagNode<std::vector<double> >*           coordinates;
        MatrixReal                                          coa;
        std::vector<double>                                 colWeights;
    };
    
}

#endif

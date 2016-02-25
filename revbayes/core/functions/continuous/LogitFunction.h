/**
 * @file
 * This file contains the declaration of the binary multiplication function, f(x) = a - b.
 *
 * @brief Declaration of functions.
 *
 * (c) Copyright 2009-
 * @date Last modified: $Date: 2012-06-20 22:57:09 +0200 (Wed, 20 Jun 2012) $
 * @author The RevBayes Development Core Team
 * @license GPL version 3
 * @version 1.0
 * @since 2012-08-13, version 1.0
 *
 * $Id: RandomNumberFactory.h 1643 2012-06-20 20:57:09Z hoehna $
 */

#ifndef LogitFunction_H
#define LogitFunction_H

#include "TypedFunction.h"
#include "TypedDagNode.h"

#include <cmath>

class LogitFunction : public TypedFunction<double> {
    
public:
    LogitFunction(const TypedDagNode<double> *p, bool inv = false);
    
    LogitFunction*                      	clone(void) const;                                                  //!< Create a clon.
    void                                    update(void);                                                       //!< Recompute the value
    
protected:
    void                                    swapParameterInternal(const DagNode *oldP, const DagNode *newP);    //!< Implementation of swaping parameters
    
private:
    const TypedDagNode<double>*     p;
    bool							inverse;
};

LogitFunction::LogitFunction(const TypedDagNode<double> *l, bool inv) :
		TypedFunction<double>( new double() ), p( l ), inverse( inv ) {
    this->addParameter( p );

}



LogitFunction* LogitFunction::clone( void ) const {
    return new LogitFunction(*this);
}


void LogitFunction::swapParameterInternal(const DagNode *oldP, const DagNode *newP) {
    if (oldP == p) {
        p = static_cast<const TypedDagNode<double>* >( newP );
    }
}



void LogitFunction::update( void ) {

	// p = exp(logit)/(1 + exp(logit))
	if(inverse)
	{
		double logit = p->getValue();

		*this->value = exp(logit)/(1 + exp(logit));
	}
	// logit = log(p / (1-p)
	else
	{
		double pval = p->getValue();

		*this->value = log(pval / (1 - pval));
	}
}

#endif

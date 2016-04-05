#ifndef LnCorrectionFunction_H
#define LnCorrectionFunction_H

#include "TypedFunction.h"
#include "StochasticNode.h"
#include "BinarySubstitutionModel.h"

#include <cmath>

class LnCorrectionFunction : public TypedFunction<std::vector<RealNumber> > {
    
public:
    LnCorrectionFunction(const StochasticNode<BinaryCharacterData >* m);
    
    // public member functions
    LnCorrectionFunction*                               clone(void) const;                                                          //!< Create an independent clone
    void                                                update(void);
    
protected:
    void                                                swapParameterInternal(const DagNode *oldP, const DagNode *newP);            //!< Implementation of swaping parameters
    
private:
    
    // members
    const StochasticNode<BinaryCharacterData >*        parameter;
    
};

#endif

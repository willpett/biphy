#include "LnCorrectionFunction.h"
#include "BinarySubstitutionModel.h"

 LnCorrectionFunction::LnCorrectionFunction(const StochasticNode<BinaryCharacterData> *p) : TypedFunction< RealNumber >( new RealNumber(0.0) ),
    parameter( p )
{
     if(dynamic_cast<const BinarySubstitutionModel*>(&(parameter->getDistribution())) == 0){
        throw Exception("LnCorrectionFunction requires BinarySubstitutionModel");
    }

    this->addParameter( p );
}


LnCorrectionFunction* LnCorrectionFunction::clone( void ) const
{
    return new LnCorrectionFunction( *this );
}


void LnCorrectionFunction::update( void ) {
    
    const BinarySubstitutionModel* model = dynamic_cast<const BinarySubstitutionModel* >(&(parameter->getDistribution()));

    *this->value = exp(model->getLnCorrection());
    
}


void LnCorrectionFunction::swapParameterInternal(const DagNode *oldP, const DagNode *newP) {
    
    if (oldP == parameter) {
        parameter = static_cast<const StochasticNode<BinaryCharacterData>* >( newP );
    }
    
}

#ifndef LnCorrectionFunction_H
#define LnCorrectionFunction_H

#include "TypedFunction.h"
#include "StochasticNode.h"

namespace RevBayesCore {
    
	template<class charType, class treeType>
    class LnCorrectionFunction : public TypedFunction<double> {
        
    public:
        LnCorrectionFunction(const StochasticNode<AbstractCharacterData>* m);
        
        // public member functions
        LnCorrectionFunction*                               clone(void) const;                                                          //!< Create an independent clone
        void                                                update(void);
        
    protected:
        void                                                swapParameterInternal(const DagNode *oldP, const DagNode *newP);            //!< Implementation of swaping parameters
        
    private:
        
        // members
        const StochasticNode<AbstractCharacterData>*        parameter;
        
    };
    
}

template <class charType, class treeType>
RevBayesCore::LnCorrectionFunction< charType, treeType>::LnCorrectionFunction(const StochasticNode<AbstractCharacterData> *p) : TypedFunction< double >( new double(0.0) ),
	parameter( p )
{
	 if(dynamic_cast<const AbstractSiteCorrectionModel<charType, treeType>*>(&(parameter->getDistribution())) == 0){
		throw RbException("LnCorrectionFunction requires AbstractSiteCorrectionModel");
	}

    this->addParameter( p );
}

template <class charType, class treeType>
RevBayesCore::LnCorrectionFunction<charType, treeType>* RevBayesCore::LnCorrectionFunction<charType, treeType>::clone( void ) const
{
	return new LnCorrectionFunction( *this );
}

template <class charType, class treeType>
void RevBayesCore::LnCorrectionFunction< charType, treeType>::update( void ) {
    
	const AbstractSiteCorrectionModel<charType,treeType>* model = dynamic_cast<const AbstractSiteCorrectionModel<charType, treeType>* >(&(parameter->getDistribution()));

	*this->value = exp(model->getLnCorrection());
    
}

template <class charType, class treeType>
void RevBayesCore::LnCorrectionFunction< charType, treeType>::swapParameterInternal(const DagNode *oldP, const DagNode *newP) {
    
    if (oldP == parameter) {
		parameter = static_cast<const StochasticNode<AbstractCharacterData>* >( newP );
    }
    
}

#endif

#ifndef VectorScaleFunction_H
#define VectorScaleFunction_H

#include "TypedFunction.h"
#include "TypedDagNode.h"

#include <vector>

namespace RevBayesCore {
    
	template <class valueType>
    class VectorScaleFunction : public TypedFunction<std::vector<valueType> > {
        
    public:
        VectorScaleFunction(const TypedDagNode<std::vector<valueType> > *vector, const TypedDagNode<valueType> *scalar);
        VectorScaleFunction(const VectorScaleFunction &n);                                                                                      //!< Copy constructor
        virtual                                            ~VectorScaleFunction(void);                                                      //!< Virtual destructor
        
        // public member functions
        VectorScaleFunction*                                clone(void) const;                                                          //!< Create an independent clone
        void                                                update(void);
        void                                                keep(DagNode* affecter);
        void                                                restore(DagNode *restorer);
        void                                                touch(DagNode *toucher );
        
    protected:
        void                                                swapParameterInternal(const DagNode *oldP, const DagNode *newP);            //!< Implementation of swaping parameters
        
    private:
        
        // members
        const TypedDagNode<std::vector<valueType> > *vector;
        const TypedDagNode<valueType> *scalar;
        
    };

template <class valueType>
void VectorScaleFunction<valueType>::keep( DagNode *toucher ) {
    this->dagNode->clearTouchedElementIndices();
}

template <class valueType>
void VectorScaleFunction<valueType>::restore( DagNode *toucher ) {

    this->update();

    this->dagNode->clearTouchedElementIndices();
}

template <class valueType>
void VectorScaleFunction<valueType>::touch( DagNode *toucher ) {

    if (toucher == vector)
	{
    	const std::set<size_t> &indices = vector->getTouchedElementIndices();

    	for(std::set<size_t>::const_iterator it = indices.begin(); it != indices.end(); it++)
    		this->dagNode->addTouchedElementIndex( *it );
	}
}

template <class valueType>
VectorScaleFunction<valueType>::VectorScaleFunction(const TypedDagNode<std::vector<valueType> > *vec, const TypedDagNode<valueType> *sca) : TypedFunction<std::vector<valueType> >( new std::vector<valueType>() ), vector( vec ), scalar(sca) {
    // add the lambda parameter as a parent
    this->addParameter( vector);
    this->addParameter( scalar);

    update();
}

template <class valueType>
VectorScaleFunction<valueType>::VectorScaleFunction(const VectorScaleFunction &n) : TypedFunction<std::vector<valueType> >( n ), vector(n.vector), scalar(n.scalar) {
    // no need to add parameters, happens automatically

    update();
}

template <class valueType>
VectorScaleFunction<valueType>::~VectorScaleFunction( void ) {
    // We don't delete the parameters, because they might be used somewhere else too. The model needs to do that!
}


template <class valueType>
VectorScaleFunction<valueType>* VectorScaleFunction<valueType>::clone( void ) const {
    return new VectorScaleFunction( *this );
}

template <class valueType>
void VectorScaleFunction<valueType>::update( void ) {

    // empty current simplex
    this->value->clear();

    std::vector<double> vec_values = vector->getValue();
    double scale = scalar->getValue();
    for (size_t i = 0; i < vec_values.size(); i++) {
    	this->value->push_back(vec_values[i]*scale);
    }
}

template <class valueType>
void VectorScaleFunction<valueType>::swapParameterInternal(const DagNode *oldP, const DagNode *newP) {

	if (oldP == vector) {
		vector = static_cast<const TypedDagNode<std::vector<valueType> >* >( newP );
	}else if(oldP == scalar){
		scalar = static_cast<const TypedDagNode<double>* >( newP );
	}

}

}
#endif

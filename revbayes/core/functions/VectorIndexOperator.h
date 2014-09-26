/**
 * @file
 * This file contains the declaration of the deterministic variable class for VectorIndexOperators.
 * This class is derived from the deterministic node and each instance will represent a deterministic variable
 * computing the VectorIndexOperator of its parameters.
 *
 * @brief Declaration of the deterministic variable for VectorIndexOperators.
 *
 * (c) Copyright 2009- under GPL version 3
 * @date Last modified: $Date$
 * @author The RevBayes Development Core Team
 * @license GPL version 3
 * @version 1.0
 * @since 2012-07-06, version 1.0
 * @interface TypedDagNode
 *
 * $Id$
 */



#ifndef VectorIndexOperator_H
#define VectorIndexOperator_H

#include "TypedFunction.h"

namespace RevBayesCore {
    
    template <class valueType>
    class VectorIndexOperator : public TypedFunction<valueType> {
        
    public:
        VectorIndexOperator(const TypedDagNode<std::vector<valueType> >* v, const TypedDagNode<int>* idx);
        VectorIndexOperator(const VectorIndexOperator &n);                                                                          //!< Copy constructor
        virtual                                            ~VectorIndexOperator(void);                                              //!< Virtual destructor
        
        // public member functions
        VectorIndexOperator*                                clone(void) const;                                                      //!< Create an independent clone
        void                                                update(void);
        
    protected:
        void                                                swapParameterInternal(const DagNode *oldP, const DagNode *newP);        //!< Implementation of swaping parameters
        
    private:
        
        // members
        const TypedDagNode<std::vector<valueType> >*        vector;
        const TypedDagNode<int>*                            index;
        
    };
    
}


template <class valueType>
RevBayesCore::VectorIndexOperator<valueType>::VectorIndexOperator( const TypedDagNode<std::vector<valueType> >* v, const TypedDagNode<int> *idx) : TypedFunction<valueType>( new valueType() ), vector( v ), index( idx ) {
    // add the vector parameter as a parent
    this->addParameter( vector );
    this->addParameter( index );
    
    update();
}


template <class valueType>
RevBayesCore::VectorIndexOperator<valueType>::VectorIndexOperator(const VectorIndexOperator<valueType> &n) : TypedFunction<valueType>( n ), vector( n.vector ), index( n.index ) {
    // no need to add parameters, happens automatically
    
    update();
}


template <class valueType>
RevBayesCore::VectorIndexOperator<valueType>::~VectorIndexOperator( void ) {
    // We don't delete the parameters, because they might be used somewhere else too. The model needs to do that!
}



template <class valueType>
RevBayesCore::VectorIndexOperator<valueType>* RevBayesCore::VectorIndexOperator<valueType>::clone( void ) const {
    return new VectorIndexOperator<valueType>( *this );
}


template <class valueType>
void RevBayesCore::VectorIndexOperator<valueType>::update( void ) {
    
    const std::vector<valueType> &v = vector->getValue();
    *(this->value) = v[index->getValue() - 1];
}



template <class valueType>
void RevBayesCore::VectorIndexOperator<valueType>::swapParameterInternal(const DagNode *oldP, const DagNode *newP) {
    
    if (oldP == vector) {
        vector = static_cast<const TypedDagNode<std::vector<valueType> >* >( newP );
    }
    else if (oldP == index) {
        index = static_cast<const TypedDagNode<int>* >( newP );
    }
    
}

#endif

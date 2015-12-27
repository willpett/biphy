#ifndef SumFunction_H
#define SumFunction_H

#include "TypedFunction.h"
#include "TypedDagNode.h"

template<class valueType>
class SumFunction : public TypedFunction<valueType> {
    
public:
    SumFunction(const TypedDagNode<std::vector<valueType> > * v);
    virtual                                            ~SumFunction(void);                                                         //!< Virtual destructor
    
    // public member functions
    SumFunction*                                        clone(void) const;                                                          //!< Create an independent clone
    void                                                update(void);
    
protected:
    void                                                swapParameterInternal(const DagNode *oldP, const DagNode *newP);            //!< Implementation of swaping parameters
    
private:
    
    // members
    const TypedDagNode<std::vector<valueType> >*              vals;
    
};

template<class valueType>
SumFunction<valueType>::SumFunction(const TypedDagNode<std::vector<valueType> > *v) : TypedFunction<valueType>( new valueType(0.0) ), vals( v )
{
    // add the parameters as parents
    this->addParameter( vals );

    update();
}

template<class valueType>
SumFunction<valueType>::~SumFunction( void )
{
    // We don't delete the parameters, because they might be used somewhere else too. The model needs to do that!
}


template<class valueType>
SumFunction<valueType>* SumFunction<valueType>::clone( void ) const
{
    return new SumFunction( *this );
}

template<class valueType>
void SumFunction<valueType>::update( void )
{

    double m = 0;
    const std::vector<valueType> &v = vals->getValue();
    for (size_t i = 0; i < v.size(); i++)
    {
        m += v[i];
    }

    *(this->value) = m ;

}


template<class valueType>
void SumFunction<valueType>::swapParameterInternal(const DagNode *oldP, const DagNode *newP)
{

    if ( oldP == vals )
    {
        vals = static_cast<const TypedDagNode<std::vector<valueType> >* >( newP );
    }

}


#endif

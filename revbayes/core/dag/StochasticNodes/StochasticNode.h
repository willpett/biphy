/**
 * @file
 * This file contains the declaration of the stochastic DAG node class, which is our base class for all stochastic DAG nodes with a specific type.
 * This class is used as the base class for all DAG nodes representing random variables. We derive each time from this class to implement
 * stochastic DAG nodes for the specific distributions.
 *
 * @brief Declaration of the stochastic DAG node base class.
 *
 * (c) Copyright 2009- under GPL version 3
 * @date Last modified: $Date$
 * @author The RevBayes Development Core Team
 * @license GPL version 3
 * @version 1.0
 * @since 2012-06-17, version 1.0
 * @interface TypedDagNode
 *
 * $Id$
 */



#ifndef StochasticNode_H
#define StochasticNode_H

#include "DynamicNode.h"
#include "BinaryCharacterData.h"

template <class valueType>
class TypedDistribution;

template <class valueType>
class StochasticNode : public DynamicNode<valueType> {

public:
 
    StochasticNode(const std::string &n, TypedDistribution<valueType> *d);
    StochasticNode(const StochasticNode<valueType> &n);                                                                             //!< Copy constructor
    virtual                                            ~StochasticNode(void);                                                       //!< Virtual destructor

    // pure virtual methods
    virtual StochasticNode<valueType>*                  clone(void) const;
    virtual void										extractValue(std::istream& is);

    // methods    
    void                                                clamp(valueType *val);                                                      //!< Clamp an observation to this random variable
    virtual TypedDistribution<valueType>&               getDistribution(void);
    virtual const TypedDistribution<valueType>&         getDistribution(void) const;
    virtual double                                      getLnProbability(void);
    virtual double                                      getLnProbabilityRatio(void);
    valueType&                                          getValue(void);
    const valueType&                                    getValue(void) const;
    bool                                                isClamped(void) const;                                                      //!< Is this DAG node clamped?
    bool                                                isStochastic(void) const;                                                   //!< Is this DAG node stochastic?
    void                                                printStructureInfo(std::ostream &o) const;                                  //!< Print the structural information (e.g. name, value-type, distribution/function, children, parents, etc.)
    void                                                redraw(void);                                                               //!< Redraw the current value of the node (applies only to stochastic nodes)
    virtual void                                        reInitialized(void);                                                        //!< The DAG was re-initialized so maybe you want to reset some stuff (delegate to distribution)
    virtual void                                        setValue(valueType *val, bool touch=true);                                  //!< Set the value of this node
    virtual void                                        setValue(const valueType &val, bool touch=true);                            //!< Set the value of this node
    void                                                setIgnoreRedraw(bool tf=true);

protected:    
    
    virtual void                                        getAffected(std::set<DagNode *>& affected, DagNode* affecter);              //!< Mark and get affected nodes
    virtual void                                        keepMe(DagNode* affecter);                                                  //!< Keep value of this and affected nodes
    virtual void                                        restoreMe(DagNode *restorer);                                               //!< Restore value of this nodes
    void                                                swapParameter(const DagNode *oldP, const DagNode *newP);                    //!< Swap the parameter of this node (needs overwriting in deterministic and stochastic nodes)
    virtual void                                        touchMe(DagNode *toucher);                                                  //!< Tell affected nodes value is reset

    // protected members
    bool                                                clamped;
    bool                                                ignoreRedraw;
    double                                              lnProb;                                                                     //!< Current log probability
    bool                                                needsProbabilityRecalculation;                                              //!< Do we need recalculation of the ln prob?
    double                                              storedLnProb;
    TypedDistribution<valueType>*                       distribution;

};

template<>
inline void StochasticNode<BinaryCharacterData >::extractValue(std::istream &is) {
    throw Exception("StochasticNode<BinaryCharacterData>::extractValue not implemented");
}


#include "../../utils/Constants.h"
#include "../../utils/Options.h"
#include "TypedDistribution.h"

template<class valueType>
void StochasticNode<valueType>::extractValue(std::istream &is) {
	if ( Utils::is_vector<valueType>::value )
	{
		size_t numElements = Utils::sub_vector<valueType>::size( getValue() );

		valueType tmp(getValue());
		is >> tmp;
		setValue(tmp);
	}
	else
	{
		valueType val;
		is >> val;
		setValue(val);
	}
}

template<class valueType>
StochasticNode<valueType>::StochasticNode( const std::string &n, TypedDistribution<valueType> *d ) : DynamicNode<valueType>( n ), clamped( false ), ignoreRedraw(false), lnProb( Constants::Double::neginf ), needsProbabilityRecalculation( true ), distribution( d ) {
    
    this->type = DagNode::STOCHASTIC;
    
    // set myself as the DAG node of the distribution
    distribution->setStochasticNode( this );
    
    // get the parameters from the distribution and add them as my parents in the DAG
    const std::set<const DagNode*> distParents = distribution->getParameters();
    
    for (std::set<const DagNode*>::iterator it = distParents.begin(); it != distParents.end(); ++it) 
    {
        this->addParent( *it );
    }
}



template<class valueType>
StochasticNode<valueType>::StochasticNode( const StochasticNode<valueType> &n ) : DynamicNode<valueType>( n ), clamped( n.clamped ), ignoreRedraw(n.ignoreRedraw), needsProbabilityRecalculation( true ), distribution( n.distribution->clone() ) {
    
    this->type = DagNode::STOCHASTIC;
    
    // set myself as the DAG node of the distribution
    distribution->setStochasticNode( this );
}


template<class valueType>
StochasticNode<valueType>::~StochasticNode( void ) {
    delete distribution;
}



template<class valueType>
void StochasticNode<valueType>::clamp(valueType *val) {
    // clamp the node with the value
    // we call set value because some derived classes might have special implementations for setting values (e.g. mixtures)
    setValue( val );
    
    clamped = true;

}



template<class valueType>
StochasticNode<valueType>* StochasticNode<valueType>::clone( void ) const {
    
    return new StochasticNode<valueType>( *this );
}



/** 
 * Get the affected nodes.
 * This call is started by the parent and since we don't have one this is a dummy implementation!
 */
template<class valueType>
void StochasticNode<valueType>::getAffected(std::set<DagNode *> &affected, DagNode* affecter) {
    
    // insert this node as one of the affected
    affected.insert( this );
    
    // call for potential specialized handling (e.g. internal flags)
    distribution->getAffected( affected, affecter );
    
}


template<class valueType>
TypedDistribution<valueType>& StochasticNode<valueType>::getDistribution( void ) {
    
    return *distribution;
}


template<class valueType>
const TypedDistribution<valueType>& StochasticNode<valueType>::getDistribution( void ) const {
    
    return *distribution;
}


template<class valueType>
double StochasticNode<valueType>::getLnProbability( void ) {
    
    if ( needsProbabilityRecalculation ) 
    {
        
        // compute and store log-probability
        lnProb = distribution->computeLnProbability();
        
        // reset flag
        needsProbabilityRecalculation = false;
    }
    
    return this->heat*lnProb;
}


template<class valueType>
double StochasticNode<valueType>::getLnProbabilityRatio( void ) {
    
	return this->heat * ( getLnProbability() - storedLnProb );
}


template<class valueType>
valueType& StochasticNode<valueType>::getValue( void ) {
    
    return distribution->getValue();
}


template<class valueType>
const valueType& StochasticNode<valueType>::getValue( void ) const {
    
    return distribution->getValue();
}


template<class valueType>
bool StochasticNode<valueType>::isClamped( void ) const {
    
    return clamped;
}


template<class valueType>
bool StochasticNode<valueType>::isStochastic( void ) const {
    
    return true;
}


/**
 * Keep the current value of the node. 
 * At this point, we also need to make sure we update the stored ln probability.
 */
template<class valueType>
void StochasticNode<valueType>::keepMe( DagNode* affecter ) {
    
    if ( this->touched ) 
    {
        
        storedLnProb = 1.0E9;       // An almost impossible value for the density
        if ( needsProbabilityRecalculation ) 
        {
            lnProb = distribution->computeLnProbability();
        }
        
        distribution->keep( affecter );
        
        // clear the list of touched element indices
        this->touchedElements.clear();
                
    }
    
    needsProbabilityRecalculation   = false;
    
    
    // delegate call
    DynamicNode<valueType>::keepMe( affecter );
    
}


template<class valueType>
void StochasticNode<valueType>::printStructureInfo( std::ostream &o ) const
{
    o << "_variableType = Stochastic DAG node" << std::endl;
    o << "_distribution = " << "<unnamed>" << std::endl;
    o << "_touched      = " << ( this->touched ? "TRUE" : "FALSE" ) << std::endl;
    o << "_clamped      = " << ( clamped ? "TRUE" : "FALSE" ) << std::endl;
    o << "_value        = " << this->distribution->getValue() << std::endl;
    o << "_lnProb       = " << lnProb << std::endl;
    if ( this->touched )
        o << "_storedLnProb = " << storedLnProb << std::endl;    
    
    o << "_parents      = ";
    this->printParents(o);
    o << std::endl;
    
    o << "_children     = ";
    this->printChildren(o);
    o << std::endl;
}


template<class valueType>
void StochasticNode<valueType>::redraw( void ) {

    // draw the value
    if (!ignoreRedraw)
        distribution->redrawValue();
    
    // touch this node for probability recalculation
    this->touch();
    
}


template<class valueType>
void StochasticNode<valueType>::reInitialized( void ) {
    
    distribution->reInitialized();
        
}


/** Restore the old value of the node and tell affected */
template<class valueType>
void StochasticNode<valueType>::restoreMe(DagNode *restorer) {
    
    if ( this->touched ) 
    {
        
        lnProb          = storedLnProb;
        storedLnProb    = 1.0E6;    // An almost impossible value for the density
                        
        // reset flags that recalculation is not needed
        needsProbabilityRecalculation = false;
        
        // call for potential specialized handling (e.g. internal flags)
        distribution->restore(restorer);
        
        // clear the list of touched element indices
        this->touchedElements.clear();
        
    }
    
    // delegate call
    DynamicNode<valueType>::restoreMe( restorer );
    
}


template<class valueType>
void StochasticNode<valueType>::setValue(valueType *val, bool forceTouch) {
    // set the value
    distribution->setValue( val );

    if ( forceTouch ) 
    {
        // touch this node for probability recalculation
        this->touch();
    }

}


template<class valueType>
void StochasticNode<valueType>::setValue(const valueType &val, bool forceTouch) {
    
    // set the value
    distribution->setValue( val );
    
    if ( forceTouch ) 
    {
        // touch this node for probability recalculation
        this->touch();
    }
    
}

template <class valueType>
void StochasticNode<valueType>::setIgnoreRedraw(bool tf)
{
    ignoreRedraw = tf;
}

template <class valueType>
void StochasticNode<valueType>::swapParameter(const DagNode *oldP, const DagNode *newP) {
    distribution->swapParameter(oldP, newP);
}



/** touch this node for recalculation */
template<class valueType>
void StochasticNode<valueType>::touchMe( DagNode *toucher ) {
    
    if (!this->touched) {
        
        storedLnProb = lnProb;
        
    }
        
    needsProbabilityRecalculation = true;
    
    // call for potential specialized handling (e.g. internal flags), we might have been touched already by someone else, so we need to delegate regardless
    distribution->touch( toucher );
    
    // delegate call
    DynamicNode<valueType>::touchMe( toucher );
}

#endif


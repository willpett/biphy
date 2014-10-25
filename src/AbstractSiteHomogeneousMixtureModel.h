#ifndef AbstractSiteHomogeneousMixtureModel_H
#define AbstractSiteHomogeneousMixtureModel_H

#include "AbstractCharEvoModel.h"
#include "DnaState.h"
#include "RateMatrix.h"
#include "RbVector.h"
#include "TopologyNode.h"
#include "TransitionProbabilityMatrix.h"
#include "TreeChangeEventListener.h"
#include "TypedDistribution.h"

namespace RevBayesCore {
    
    template<class charType, class treeType>
    class AbstractSiteHomogeneousMixtureModel : public virtual AbstractCharEvoModel<charType, treeType> {
        
    public:
        AbstractSiteHomogeneousMixtureModel(const TypedDagNode< treeType > *t, size_t nChars, bool c, size_t nSites, size_t nMix = 1);
        AbstractSiteHomogeneousMixtureModel(const AbstractSiteHomogeneousMixtureModel &n);                                                                                                //!< Copy constructor
        virtual                                            ~AbstractSiteHomogeneousMixtureModel(void);                                                                   //!< Virtual destructor
        

        void                                                setSiteRates(const TypedDagNode< std::vector< double > > *r);
        virtual void                                        swapParameter(const DagNode *oldP, const DagNode *newP);                                    //!< Implementation of swaping parameters
        
    protected:
        const std::vector< double > &						getSiteRates(void);
        virtual void                                        touchSpecialization(DagNode *toucher);

        size_t                                              mixtureOffset;
        size_t                                              numSiteRates;
        
        const TypedDagNode< std::vector< double > >*        siteRates;
        std::vector< double >								oneRate;
    protected:
        bool												rateVariationAcrossSites;



    };
    
}


#include "ConstantNode.h"
#include "DiscreteCharacterState.h"
#include "RandomNumberFactory.h"
#include "TopologyNode.h"
#include "TransitionProbabilityMatrix.h"

#include <cmath>
#include <cstring>

template<class charType, class treeType>
RevBayesCore::AbstractSiteHomogeneousMixtureModel<charType, treeType>::AbstractSiteHomogeneousMixtureModel(const TypedDagNode< treeType > *t, size_t nChars, bool c, size_t nSites, size_t nMix) :
	AbstractCharEvoModel<charType, treeType>(t, nChars, c, nSites),
	numSiteRates( nMix ), oneRate(std::vector<double>(1,1.0))
{
    
    // initialize with default parameters
    siteRates                   = NULL;

    rateVariationAcrossSites 	= false;

    mixtureOffset               =  this->numPatterns*this->numChars;
}


template<class charType, class treeType>
RevBayesCore::AbstractSiteHomogeneousMixtureModel<charType, treeType>::AbstractSiteHomogeneousMixtureModel(const AbstractSiteHomogeneousMixtureModel &d) : AbstractCharEvoModel<charType, treeType>( d ) {
    // parameters are automatically copied
    // initialize with default parameters
    siteRates                   = d.siteRates;
    mixtureOffset 				= d.mixtureOffset;
    numSiteRates				= d.numSiteRates;
    oneRate						= d.oneRate;
    
    
    // flags specifying which model variants we use
    rateVariationAcrossSites                    = d.rateVariationAcrossSites;
}


template<class charType, class treeType>
RevBayesCore::AbstractSiteHomogeneousMixtureModel<charType, treeType>::~AbstractSiteHomogeneousMixtureModel( void ) {
    // We don't delete the parameters, because they might be used somewhere else too. The model needs to do that!
    
}


template<class charType, class treeType>
void RevBayesCore::AbstractSiteHomogeneousMixtureModel<charType, treeType>::setSiteRates(const TypedDagNode< std::vector< double > > *r) {
    
    // remove the old parameter first
    if ( siteRates != NULL )
    {
        this->removeParameter( siteRates );
        siteRates = NULL;
    }
    
    if ( r != NULL )
    {
        // set the value
        rateVariationAcrossSites = true;
        siteRates = r;
        this->numSiteRates = r->getValue().size();
        this->resizeLikelihoodVectors();
    }
    else
    {
        // set the value
        rateVariationAcrossSites = false;
        siteRates = NULL;
        this->numSiteRates = 1;
        this->resizeLikelihoodVectors();
        
    }
    
    // add the parameter
    this->addParameter( siteRates );

    // redraw the current value
    if ( this->dagNode != NULL && !this->dagNode->isClamped() )
    {
        this->redrawValue();
    }
}

template<class charType, class treeType>
const std::vector<double> & RevBayesCore::AbstractSiteHomogeneousMixtureModel<charType, treeType>::getSiteRates(void) {

	if(rateVariationAcrossSites){
		return siteRates->getValue();
	}else{
		return oneRate;
	}
}




template<class charType, class treeType>
void RevBayesCore::AbstractSiteHomogeneousMixtureModel<charType, treeType>::swapParameter(const DagNode *oldP, const DagNode *newP) {
    
    if (oldP == siteRates)
    {
        siteRates = static_cast<const TypedDagNode< std::vector< double > >* >( newP );
    }
    else
    {
    	AbstractCharEvoModel<charType, treeType>::swapParameter(oldP,newP);
    }
    
}

template<class charType, class treeType>
void RevBayesCore::AbstractSiteHomogeneousMixtureModel<charType, treeType>::touchSpecialization( DagNode* affecter ) {

    AbstractCharEvoModel<charType, treeType>::touchSpecialization(affecter);

}


#endif

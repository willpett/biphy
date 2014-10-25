#ifndef AbstractBranchHeterogeneousModel_H
#define AbstractBranchHeterogeneousModel_H

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
    class AbstractBranchHeterogeneousModel : public virtual AbstractCharEvoModel<charType, treeType> {
        
    public:
        AbstractBranchHeterogeneousModel(const TypedDagNode< treeType > *t, size_t nChars, bool c, size_t nSites);
        AbstractBranchHeterogeneousModel(const AbstractBranchHeterogeneousModel &n);                                                                                                //!< Copy constructor
        virtual                                            ~AbstractBranchHeterogeneousModel(void);                                                                   //!< Virtual destructor
        
        // public member functions
        void                                                setClockRate(const TypedDagNode< double > *r);
        void                                                setClockRate(const TypedDagNode< std::vector< double > > *r);
        void                                                setRateMatrix(const TypedDagNode< RateMatrix > *rm);
        void                                                setRateMatrix(const TypedDagNode< RbVector< RateMatrix > > *rm);
        void                                                setRootFrequencies(const TypedDagNode< std::vector< double > > *f);
        virtual void                                        swapParameter(const DagNode *oldP, const DagNode *newP);                                    //!< Implementation of swaping parameters
        
    protected:
        
        const std::vector<double>&                          getRootFrequencies(void);
        const RateMatrix*									getRateMatrix(size_t index);
        double                                              getClockRate(size_t index);
        virtual void                                        touchSpecialization(DagNode *toucher);
        
        
    private:        
        // members
        const TypedDagNode< double >*                       homogeneousClockRate;
        const TypedDagNode< std::vector< double > >*        heterogeneousClockRates;
        const TypedDagNode< RateMatrix >*                   homogeneousRateMatrix;
        const TypedDagNode< RbVector< RateMatrix > >*       heterogeneousRateMatrices;
        const TypedDagNode< std::vector< double > >*        rootFrequencies;
        
        
        // flags specifying which model variants we use
        bool                                                branchHeterogeneousClockRates;
        bool                                                branchHeterogeneousSubstitutionMatrices;
    };
    
}


#include "ConstantNode.h"
#include "DiscreteCharacterState.h"
#include "RandomNumberFactory.h"
#include "TopologyNode.h"

#include <cmath>
#include <cstring>

template<class charType, class treeType>
RevBayesCore::AbstractBranchHeterogeneousModel<charType, treeType>::AbstractBranchHeterogeneousModel(const TypedDagNode< treeType > *t, size_t nChars, bool c, size_t nSites) :
	AbstractCharEvoModel<charType, treeType>(t, nChars, c, nSites)
{
    // initialize with default parameters
    homogeneousClockRate        = new ConstantNode<double>("clockRate", new double(1.0) );
    heterogeneousClockRates     = NULL;
    homogeneousRateMatrix       = NULL;
    heterogeneousRateMatrices   = NULL;
    rootFrequencies             = NULL;

    
    // flags specifying which model variants we use
    branchHeterogeneousClockRates               = false;
    branchHeterogeneousSubstitutionMatrices     = false;
    
    // add the parameters to the parents list
    this->addParameter( homogeneousClockRate );
    this->addParameter( homogeneousRateMatrix );
}


template<class charType, class treeType>
RevBayesCore::AbstractBranchHeterogeneousModel<charType, treeType>::AbstractBranchHeterogeneousModel(const AbstractBranchHeterogeneousModel &d) : AbstractCharEvoModel<charType, treeType>( d ) {
    // parameters are automatically copied
    // initialize with default parameters
    homogeneousClockRate        = d.homogeneousClockRate;
    heterogeneousClockRates     = d.heterogeneousClockRates;
    homogeneousRateMatrix       = d.homogeneousRateMatrix;
    heterogeneousRateMatrices   = d.heterogeneousRateMatrices;
    rootFrequencies             = d.rootFrequencies;
    
    
    // flags specifying which model variants we use
    branchHeterogeneousClockRates               = d.branchHeterogeneousClockRates;
    branchHeterogeneousSubstitutionMatrices     = d.branchHeterogeneousSubstitutionMatrices;
}


template<class charType, class treeType>
RevBayesCore::AbstractBranchHeterogeneousModel<charType, treeType>::~AbstractBranchHeterogeneousModel( void ) {
    // We don't delete the parameters, because they might be used somewhere else too. The model needs to do that!
    
}


template<class charType, class treeType>
const std::vector<double>& RevBayesCore::AbstractBranchHeterogeneousModel<charType, treeType>::getRootFrequencies( void ) {
    
    if ( branchHeterogeneousSubstitutionMatrices || rootFrequencies != NULL ) 
    {
    	return rootFrequencies->getValue();
    } 
    else 
    {
        return homogeneousRateMatrix->getValue().getStationaryFrequencies();
    }

}



template<class charType, class treeType>
void RevBayesCore::AbstractBranchHeterogeneousModel<charType, treeType>::setClockRate(const TypedDagNode< double > *r) {
    
    // remove the old parameter first
    if ( homogeneousClockRate != NULL ) 
    {
        this->removeParameter( homogeneousClockRate );
        homogeneousClockRate = NULL;
    }
    else // heterogeneousClockRate != NULL
    {
        this->removeParameter( heterogeneousClockRates );
        heterogeneousClockRates = NULL;
    }
    
    // set the value
    branchHeterogeneousClockRates = false;
    homogeneousClockRate = r;
    
    // add the parameter
    this->addParameter( homogeneousClockRate );
    
    // redraw the current value
    if ( this->dagNode != NULL && !this->dagNode->isClamped() ) 
    {
        this->redrawValue();
    }
    
}



template<class charType, class treeType>
void RevBayesCore::AbstractBranchHeterogeneousModel<charType, treeType>::setClockRate(const TypedDagNode< std::vector< double > > *r) {
    
    // remove the old parameter first
    if ( homogeneousClockRate != NULL ) 
    {
        this->removeParameter( homogeneousClockRate );
        homogeneousClockRate = NULL;
    }
    else // heterogeneousClockRate != NULL
    {
        this->removeParameter( heterogeneousClockRates );
        heterogeneousClockRates = NULL;
    }
    
    // set the value
    branchHeterogeneousClockRates = true;
    heterogeneousClockRates = r;
    
    // add the parameter
    this->addParameter( heterogeneousClockRates );
    
    // redraw the current value
    if ( this->dagNode != NULL && !this->dagNode->isClamped() ) 
    {
        this->redrawValue();
    }
    
}


template<class charType, class treeType>
void RevBayesCore::AbstractBranchHeterogeneousModel<charType, treeType>::setRateMatrix(const TypedDagNode< RateMatrix > *rm) {
    
    // remove the old parameter first
    if ( homogeneousRateMatrix != NULL ) 
    {
        this->removeParameter( homogeneousRateMatrix );
        homogeneousRateMatrix = NULL;
    }
    else // heterogeneousRateMatrix != NULL
    {
        this->removeParameter( heterogeneousRateMatrices );
        heterogeneousRateMatrices = NULL;
    }
    
    // set the value
    branchHeterogeneousSubstitutionMatrices = false;
    homogeneousRateMatrix = rm;
    
    // add the parameter
    this->addParameter( homogeneousRateMatrix );
    
    // redraw the current value
    if ( this->dagNode != NULL && !this->dagNode->isClamped() ) 
    {
        this->redrawValue();
    }
    
}


template<class charType, class treeType>
void RevBayesCore::AbstractBranchHeterogeneousModel<charType, treeType>::setRateMatrix(const TypedDagNode< RbVector< RateMatrix > > *rm) {
    
    // remove the old parameter first
    if ( homogeneousRateMatrix != NULL ) 
    {
        this->removeParameter( homogeneousRateMatrix );
        homogeneousRateMatrix = NULL;
    }
    else // heterogeneousRateMatrix != NULL
    {
        this->removeParameter( heterogeneousRateMatrices );
        heterogeneousRateMatrices = NULL;
    }
    
    // set the value
    branchHeterogeneousSubstitutionMatrices = true;
    heterogeneousRateMatrices = rm;
    
    // add the parameter
    this->addParameter( heterogeneousRateMatrices );
    
    // redraw the current value
    if ( this->dagNode != NULL && !this->dagNode->isClamped() ) 
    {
        this->redrawValue();
    }
    
}


template<class charType, class treeType>
void RevBayesCore::AbstractBranchHeterogeneousModel<charType, treeType>::setRootFrequencies(const TypedDagNode< std::vector< double > > *f) {
    
    // remove the old parameter first
    if ( rootFrequencies != NULL ) 
    {
        this->removeParameter( rootFrequencies );
        rootFrequencies = NULL;
    }
    
    if ( f != NULL )
    {
        // set the value
//        branchHeterogeneousSubstitutionMatrices = true;
        rootFrequencies = f;
    
        // add the parameter
        this->addParameter( rootFrequencies );
    }
    else
    {
        branchHeterogeneousSubstitutionMatrices = false;
    }
    
    // redraw the current value
    if ( this->dagNode != NULL && !this->dagNode->isClamped() ) 
    {
        this->redrawValue();
    }
}




template<class charType, class treeType>
void RevBayesCore::AbstractBranchHeterogeneousModel<charType, treeType>::swapParameter(const DagNode *oldP, const DagNode *newP) {
    
    if (oldP == homogeneousClockRate) 
    {
        homogeneousClockRate = static_cast<const TypedDagNode< double >* >( newP );
    }
    else if (oldP == heterogeneousClockRates) 
    {
        heterogeneousClockRates = static_cast<const TypedDagNode< std::vector< double > >* >( newP );
    }
    else if (oldP == homogeneousRateMatrix) 
    {
        homogeneousRateMatrix = static_cast<const TypedDagNode< RateMatrix >* >( newP );
    }
    else if (oldP == heterogeneousRateMatrices) 
    {
        heterogeneousRateMatrices = static_cast<const TypedDagNode< RbVector< RateMatrix > >* >( newP );
    }
    else if (oldP == rootFrequencies) 
    {
        rootFrequencies = static_cast<const TypedDagNode< std::vector< double > >* >( newP );
    }
    else 
    {
    	AbstractCharEvoModel<charType, treeType>::swapParameter(oldP,newP);
    }
    
}

template<class charType, class treeType>
void RevBayesCore::AbstractBranchHeterogeneousModel<charType, treeType>::touchSpecialization( DagNode* affecter ) {
    
    // if the topology wasn't the culprit for the touch, then we just flag everything as dirty
    if ( affecter == heterogeneousClockRates ) 
    {
        const std::set<size_t> &indices = heterogeneousClockRates->getTouchedElementIndices();
        
        // maybe all of them have been touched or the flags haven't been set properly
        if ( indices.size() == 0 ) 
        {
            // just delegate the call
        	AbstractCharEvoModel<charType, treeType>::touchSpecialization( affecter );
        } 
        else 
        {
            const std::vector<TopologyNode *> &nodes = this->tau->getValue().getNodes();
            // flag recomputation only for the nodes
            for (std::set<size_t>::iterator it = indices.begin(); it != indices.end(); ++it) 
            {
                this->recursivelyFlagNodeDirty( *nodes[*it] );
            } 
        }
    }
    else if ( affecter == heterogeneousRateMatrices )
    {
        
        const std::set<size_t> &indices = heterogeneousRateMatrices->getTouchedElementIndices();
        
        // maybe all of them have been touched or the flags haven't been set properly
        if ( indices.size() == 0 ) 
        {
            // just delegate the call
        	AbstractCharEvoModel<charType, treeType>::touchSpecialization( affecter );
        } 
        else 
        {
            const std::vector<TopologyNode *> &nodes = this->tau->getValue().getNodes();
            // flag recomputation only for the nodes
            for (std::set<size_t>::iterator it = indices.begin(); it != indices.end(); ++it) 
            {
                this->recursivelyFlagNodeDirty( *nodes[*it] );
            } 
        }
    }
    else if ( affecter == rootFrequencies )
    {
        
        const TopologyNode &root = this->tau->getValue().getRoot();
        this->recursivelyFlagNodeDirty( root );
    }
    else
    {
    	AbstractCharEvoModel<charType, treeType>::touchSpecialization( affecter );
    }
    
}

template<class charType, class treeType>
const RevBayesCore::RateMatrix* RevBayesCore::AbstractBranchHeterogeneousModel<charType, treeType>::getRateMatrix(size_t nodeIdx) {

    // first, get the rate matrix for this branch
    const RateMatrix *rm;
    if ( this->branchHeterogeneousSubstitutionMatrices == true )
    {
        rm = &this->heterogeneousRateMatrices->getValue()[nodeIdx];
    }
    else
    {
        rm = &this->homogeneousRateMatrix->getValue();
    }

    return rm;
}

template<class charType, class treeType>
double RevBayesCore::AbstractBranchHeterogeneousModel<charType, treeType>::getClockRate(size_t nodeIdx) {

    double branchTime;
    if ( this->branchHeterogeneousClockRates == true )
    {
        branchTime = this->heterogeneousClockRates->getValue()[nodeIdx];
    }
    else
    {
        branchTime = this->homogeneousClockRate->getValue();
    }

    return branchTime;

}


#endif

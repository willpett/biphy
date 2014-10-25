#ifndef BinaryCharEvoModel_H
#define BinaryCharEvoModel_H

#include "AbstractSiteCorrectionModel.h"
#include "StandardState.h"
#include "RateMatrix.h"
#include "RbVector.h"
#include "TopologyNode.h"
#include "TransitionProbabilityMatrix.h"
#include "TreeChangeEventListener.h"
#include "TypedDistribution.h"

namespace RevBayesCore {
    
	enum CorrectionType { NONE 					= 0x00,
						  NO_ABSENT_SITES		= 0x01,
						  NO_PRESENT_SITES		= 0x02,
						  NO_CONSTANT_SITES		= 0x03,
						  NO_SINGLETON_GAINS 	= 0x04,
						  NO_SINGLETON_LOSSES	= 0x08,
						  NO_SINGLETONS			= 0x0C,
						  NO_UNINFORMATIVE		= 0x0F
						};

    template<class treeType>
    class BinaryCharEvoModel : public AbstractSiteCorrectionModel<StandardState, treeType> {

    public:
    	BinaryCharEvoModel(const TypedDagNode< treeType > *p, bool c, size_t nSites, int type = 0);
        BinaryCharEvoModel(const BinaryCharEvoModel &n);                                                                                                //!< Copy constructor
        virtual                                            ~BinaryCharEvoModel(void);                                                                   //!< Virtual destructor
        
        // public member functions
        BinaryCharEvoModel*         						clone(void) const;
        virtual void                                        redrawValue(void);

    protected:
        
        virtual void                                        setCorrectionPatterns();
        size_t 												simulate( const TopologyNode &node, std::vector<StandardState> &taxa, size_t rateIndex);
        
        int													type;

    };
    
}


#include "ConstantNode.h"
#include "DiscreteCharacterState.h"
#include "RandomNumberFactory.h"
#include "TopologyNode.h"
#include "TransitionProbabilityMatrix.h"

#include <cmath>
#include <cstring>

template<class treeType>
RevBayesCore::BinaryCharEvoModel<treeType>::BinaryCharEvoModel(const TypedDagNode< treeType > *t, bool c, size_t nSites, int type) :
	AbstractCharEvoModel<StandardState, treeType>(t, 2, c, nSites),
	AbstractSiteCorrectionModel<StandardState, treeType>(t, 2, c , nSites), type(type)
{
	size_t numTips = this->tau->getValue().getNumberOfTips();

	this->numCorrectionSites = (bool)(type & NO_ABSENT_SITES)
							 + (bool)(type & NO_PRESENT_SITES)
							 + numTips*(bool)(type & NO_SINGLETON_GAINS)
							 + numTips*(bool)(type & NO_SINGLETON_LOSSES);

	this->resizeLikelihoodVectors();
}


template<class treeType>
RevBayesCore::BinaryCharEvoModel<treeType>::BinaryCharEvoModel(const BinaryCharEvoModel &d) :
	AbstractCharEvoModel<StandardState, treeType>( d ),
	AbstractSiteCorrectionModel<StandardState, treeType>(d), type(d.type)
{
}


template<class treeType>
RevBayesCore::BinaryCharEvoModel<treeType>::~BinaryCharEvoModel( void ) {
    // We don't delete the parameters, because they might be used somewhere else too. The model needs to do that!
    
}


template<class treeType>
RevBayesCore::BinaryCharEvoModel<treeType>* RevBayesCore::BinaryCharEvoModel<treeType>::clone( void ) const {
    
    return new BinaryCharEvoModel<treeType>( *this );
}

template<class treeType>
void RevBayesCore::BinaryCharEvoModel<treeType>::setCorrectionPatterns(){

	size_t numTips = this->tau->getValue().getNumberOfTips();
	std::vector<TopologyNode*> nodes = this->tau->getValue().getNodes();
	// allocate and fill the cells of the matrices

	size_t nonSingletons = (bool)(type & NO_ABSENT_SITES) + (bool)(type & NO_PRESENT_SITES);
	size_t i = 0;
	for (std::vector<TopologyNode*>::iterator it = nodes.begin(); it != nodes.end(); ++it)
	{
		if ( (*it)->isTip() )
		{
			std::string name = (*it)->getName();

			// resize the column
			this->correctionCharMatrix[name].resize(this->numCorrectionSites);
			this->correctionGapMatrix[name].resize(this->numCorrectionSites);

			for (size_t site = 0; site < this->numCorrectionSites; ++site)
			{
				this->correctionGapMatrix[name][site] = false;
				if(site < nonSingletons){
					this->correctionCharMatrix[name][site] = site + !(type & NO_ABSENT_SITES);
				}else{
					if(site - nonSingletons < numTips)
						if(site - nonSingletons == i){
							this->correctionCharMatrix[name][site] = !(type & NO_SINGLETON_LOSSES);
						}else{
							this->correctionCharMatrix[name][site] = (bool)(type & NO_SINGLETON_LOSSES);
						}
					else
						if(site - nonSingletons == numTips + i){
							this->correctionCharMatrix[name][site] = 1;
						}else{
							this->correctionCharMatrix[name][site] = 0;
						}
				}
			}
			i++;
		}
	}
}

template<class treeType>
void RevBayesCore::BinaryCharEvoModel<treeType>::redrawValue( void ) {

    // delete the old value first
    delete this->value;

    // create a new character data object
    this->value = new DiscreteCharacterData<StandardState>();

    size_t numTips = this->tau->getValue().getNumberOfTips();

    for (size_t t = 0; t < numTips; ++t)
	{
    	DiscreteTaxonData<StandardState> data;
    	data.setTaxonName( this->tau->getValue().getNode(t).getName() );
		this->value->addTaxonData(data);
	}

    RandomNumberGenerator* rng = GLOBAL_RNG;
    //std::cerr << omega << std::endl;
    size_t cap = 0;
    for ( size_t i = 0; i < this->numSites; i++ )
    {
    	double u = rng->uniform01();
        size_t rateIndex = (int)(u*this->numSiteRates);

        std::vector<StandardState> taxa(this->tau->getValue().getNumberOfNodes(), StandardState());

        const std::vector< double > &pi = this->getRootFrequencies();

        if(rng->uniform01() < pi[1])
        	taxa[this->tau->getValue().getRoot().getIndex()]++;

		// recursively simulate the sequences
		size_t numLeaves = simulate( this->tau->getValue().getRoot(), taxa, rateIndex );

		if((this->type & NO_ABSENT_SITES) && numLeaves == 0){
			i--;
			continue;
		}else if((this->type & NO_PRESENT_SITES) && numLeaves == numTips){
			i--;
			continue;
		}else if((this->type & NO_SINGLETON_GAINS) && numLeaves == 1){
			i--;
			continue;
		}else if((this->type & NO_SINGLETON_LOSSES) && numLeaves == numTips - 1){
			i--;
			continue;
		}

		// add the taxon data to the character data
		for (size_t t = 0; t < numTips; ++t)
		{
			if(t == 22 && taxa[t].getState() == 2)
				cap++;
			this->value->getTaxonData(t).addCharacter(taxa[t]);
		}
    }

    std::cerr << (double)cap/this->numSites << std::endl;
    // compress the data and initialize internal variables
    this->compress();

}

template<class treeType>
size_t RevBayesCore::BinaryCharEvoModel<treeType>::simulate( const TopologyNode &node, std::vector<StandardState> &taxa, size_t rateIndex) {

	size_t numLeaves = 0;
    // get the children of the node
    const std::vector<TopologyNode*>& children = node.getChildren();

    // get the sequence of this node
    size_t nodeIndex = node.getIndex();
    StandardState &parentState = taxa[ nodeIndex ];

    if(parentState.getState() == 2 && node.isTip())
    	numLeaves++;

    // simulate the sequence for each child
    RandomNumberGenerator* rng = GLOBAL_RNG;
    for (std::vector< TopologyNode* >::const_iterator it = children.begin(); it != children.end(); ++it)
    {
        const TopologyNode &child = *(*it);

        // update the transition probability matrix
        this->updateTransitionProbabilities( child.getIndex(), child.getBranchLength() );

        StandardState &childState = taxa[ child.getIndex() ];

        unsigned long cp = parentState.getState();

		double *pij = this->transitionProbMatrices[ rateIndex ][cp-1];

		if(rng->uniform01() < pij[1])
			childState++;

		numLeaves += simulate( child, taxa, rateIndex );
	}

    return numLeaves;

}


#endif

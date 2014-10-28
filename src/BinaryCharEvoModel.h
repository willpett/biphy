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
        DiscreteTaxonData<StandardState> *					getDolloCompatible(void);

    protected:
        
        virtual void                                        setCorrectionPatterns();
        size_t 												simulate( const TopologyNode &node, std::vector<StandardState> &taxa, size_t rateIndex);
        
        size_t												numBirths(size_t site, const std::vector<DiscreteTaxonData<StandardState> > &mapping, const TopologyNode &node);

        int													type;
        DiscreteTaxonData<StandardState> *					dolloCompatible;

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
	AbstractSiteCorrectionModel<StandardState, treeType>(t, 2, c , nSites), type(type),
	dolloCompatible(new DiscreteTaxonData<StandardState>())
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
	AbstractSiteCorrectionModel<StandardState, treeType>(d), type(d.type), dolloCompatible(new DiscreteTaxonData<StandardState>())
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

	this->correctionCharMatrix.clear();
	this->correctionGapMatrix.clear();

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
    size_t numNodes = this->tau->getValue().getNumberOfNodes();

    RandomNumberGenerator* rng = GLOBAL_RNG;

    const std::vector< double > &pi = this->getRootFrequencies();

    const TopologyNode &root = this->tau->getValue().getRoot();
    size_t rootIndex = this->tau->getValue().getRoot().getIndex();

    std::vector< DiscreteTaxonData<StandardState> > taxa = std::vector< DiscreteTaxonData<StandardState> >(numTips, DiscreteTaxonData<StandardState>() );

    for ( size_t i = 0; i < this->numSites; i++ )
    {
    	double u = rng->uniform01();
        size_t rateIndex = (int)(u*this->numSiteRates);

        std::vector<StandardState> siteData(numNodes, StandardState());

        if(rng->uniform01() < pi[1])
        	siteData[rootIndex]++;

		// recursively simulate the sequences
		size_t numLeaves = simulate( root, siteData, rateIndex );

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
			taxa[t].addCharacter(siteData[this->tau->getValue().getNode(t).getIndex()]);
		}
    }

    // add the taxon data to the character data
	for (size_t i = 0; i < this->tau->getValue().getNumberOfTips(); ++i)
	{
		taxa[i].setTaxonName(this->tau->getValue().getNode(i).getName());
		this->value->addTaxonData( taxa[i] );
	}

    // compress the data and initialize internal variables
    this->compress();

}

template<class treeType>
size_t RevBayesCore::BinaryCharEvoModel<treeType>::simulate( const TopologyNode &node, std::vector<StandardState> &data, size_t rateIndex) {

	size_t numLeaves = 0;
    // get the children of the node
    const std::vector<TopologyNode*>& children = node.getChildren();

    // get the sequence of this node
    size_t nodeIndex = node.getIndex();
    StandardState &parentState = data[ nodeIndex ];

    // simulate the sequence for each child
    RandomNumberGenerator* rng = GLOBAL_RNG;
    for (std::vector< TopologyNode* >::const_iterator it = children.begin(); it != children.end(); ++it)
    {
        const TopologyNode &child = *(*it);

        // update the transition probability matrix
        this->updateTransitionProbabilities( child.getIndex(), child.getBranchLength() );

        StandardState &childState = data[ child.getIndex() ];

        unsigned long cp = parentState.getState();

		double *pij = this->transitionProbMatrices[ rateIndex ][cp-1];

		if(rng->uniform01() < pij[1])
			childState++;

		if(child.isTip())
		    numLeaves = childState.getState() == 2;
		else
			numLeaves += simulate( child, data, rateIndex );
	}

    return numLeaves;

}

template<class treeType>
RevBayesCore::DiscreteTaxonData<RevBayesCore::StandardState>* RevBayesCore::BinaryCharEvoModel<treeType>::getDolloCompatible(void) {

	delete dolloCompatible;
	dolloCompatible = new DiscreteTaxonData<StandardState>();

    const std::vector<DiscreteTaxonData<StandardState> > &mapping = this->getMapping();

    const TopologyNode &root = this->tau->getValue().getRoot();
    const TopologyNode &left = root.getChild(0);
    const TopologyNode &right = root.getChild(1);

    for (size_t i = 0; i < this->numSites; ++i)
    {
    	StandardState rootState = mapping[root.getIndex()].getCharacter(i);
    	StandardState leftState = mapping[left.getIndex()].getCharacter(i);
    	StandardState rightState = mapping[right.getIndex()].getCharacter(i);

    	size_t births = 0;
    	if(rootState.getState() == 1 && leftState.getState() == 2)
			births++;
		if(rootState.getState() == 1 && rightState.getState() == 2)
			births++;

    	if(!left.isTip())
    		births += numBirths(i,mapping,left);
    	if(!right.isTip())
    		births += numBirths(i,mapping,right);

    	StandardState c;
    	c.setToFirstState();
    	if(births <= 1)
    		c++;

    	dolloCompatible->addCharacter(c);
    }

    return dolloCompatible;
}

template<class treeType>
size_t RevBayesCore::BinaryCharEvoModel<treeType>::numBirths(size_t site, const std::vector<DiscreteTaxonData<StandardState> > &mapping, const TopologyNode &node) {

	const TopologyNode &left = node.getChild(0);
	const TopologyNode &right = node.getChild(1);

	StandardState nodeState = mapping[node.getIndex()].getCharacter(site);
	StandardState leftState = mapping[left.getIndex()].getCharacter(site);
	StandardState rightState = mapping[right.getIndex()].getCharacter(site);

	size_t births = 0;
	if(nodeState.getState() == 1 && leftState.getState() == 2)
		births++;
	if(nodeState.getState() == 1 && rightState.getState() == 2)
		births++;

	if(!left.isTip())
		births += numBirths(site,mapping,left);
	if(!right.isTip())
		births += numBirths(site,mapping,right);

	return births;

}


#endif

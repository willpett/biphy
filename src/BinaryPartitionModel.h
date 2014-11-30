#ifndef BinaryPartitionModel_H
#define BinaryPartitionModel_H

#include "StandardState.h"
#include "RateMatrix.h"
#include "RbVector.h"
#include "RbStatisticsHelper.h"
#include "DistributionPoisson.h"
#include "DistributionGamma.h"
#include "PartitionedSiteCorrectionModel.h"
#include "TopologyNode.h"
#include "TransitionProbabilityMatrix.h"
#include "TreeChangeEventListener.h"
#include "TypedDistribution.h"

namespace RevBayesCore {

	enum BinaryCorrectionType {
						  NO_ABSENT_SITES		= 0x01,
						  NO_PRESENT_SITES		= 0x02,
						  NO_SINGLETON_GAINS 	= 0x04,
						  NO_SINGLETON_LOSSES	= 0x08,
						  NO_SINGLETONS			= 0x0C
						};

    template<class treeType>
    class BinaryPartitionModel : public PartitionedSiteCorrectionModel<StandardState, treeType> {

    public:
    	BinaryPartitionModel(const TypedDagNode< treeType > *p, bool c, size_t nSites, int type = 0);
        BinaryPartitionModel(const BinaryPartitionModel &n);                                                                                                //!< Copy constructor
        virtual                                            ~BinaryPartitionModel(void);                                                                   //!< Virtual destructor
        
        // public member functions
        BinaryPartitionModel*         						clone(void) const;

        virtual void                                        redrawValue(void);
        DiscreteTaxonData<StandardState> *					getDolloCompatible(void);

    protected:
        
		//void                                                computeRootCorrections(size_t root, size_t l, size_t r);
		//void                                                computeInternalNodeCorrections(const TopologyNode &n, size_t nIdx, size_t l, size_t r);
		//void                                                computeTipCorrections(const TopologyNode &node, size_t nIdx);

        void                                        		touchSpecialization(DagNode *toucher);
        virtual void                                        setCorrectionPatterns();
        
        size_t												numBirths(size_t site, const std::vector<DiscreteTaxonData<StandardState> > &mapping, const TopologyNode &node);

        DiscreteTaxonData<StandardState> *					dolloCompatible;

        std::string											presentState;
        std::string											absentState;

    private:

        size_t 												simulate( const TopologyNode &node, std::vector<StandardState> &taxa, size_t rateIndex);

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
RevBayesCore::BinaryPartitionModel<treeType>::BinaryPartitionModel(const TypedDagNode< treeType > *t, bool c, size_t nSites, int type) :
	AbstractCharEvoModel<StandardState, treeType>(t, 2, c, nSites),
	PartitionedSiteCorrectionModel<StandardState, treeType>(t, 2, c , nSites, type),
	dolloCompatible(new DiscreteTaxonData<StandardState>())
{
	presentState = "1";
	absentState = "0";

}


template<class treeType>
RevBayesCore::BinaryPartitionModel<treeType>::BinaryPartitionModel(const BinaryPartitionModel &d) :
	AbstractCharEvoModel<StandardState, treeType>( d ),
	PartitionedSiteCorrectionModel<StandardState, treeType>(d), dolloCompatible(new DiscreteTaxonData<StandardState>())
{
}


template<class treeType>
RevBayesCore::BinaryPartitionModel<treeType>::~BinaryPartitionModel( void ) {
    // We don't delete the parameters, because they might be used somewhere else too. The model needs to do that!
    
}


template<class treeType>
RevBayesCore::BinaryPartitionModel<treeType>* RevBayesCore::BinaryPartitionModel<treeType>::clone( void ) const {
    
    return new BinaryPartitionModel<treeType>( *this );
}

template<class treeType>
void RevBayesCore::BinaryPartitionModel<treeType>::setCorrectionPatterns(){

	this->numCorrectionPartitions = 0;
	this->numSitesPerPartition.clear();
	this->numGapsPerPartition.clear();
	this->numCorrectionSitesPerPartition.clear();
	this->perMixtureCorrections.clear();
	this->perPartitionLnCorrections.clear();
	this->patternsPerPartition.clear();

	if(this->type == NONE)
		return;

	std::map<std::string, size_t> gapPatterns;
	for(size_t site = 0; site < this->numPatterns; site++){
		std::stringstream gapstream;
		for(std::map<std::string, std::vector<bool> >::iterator it = this->gapMatrix.begin(); it != this->gapMatrix.end(); it++)
			gapstream << it->second[site];

		std::string gaps = gapstream.str();

		if(gapPatterns.find(gaps) == gapPatterns.end()){
			if(count(gaps.begin(), gaps.end(), '1') == 0)
				this->patternsPerPartition.insert(this->patternsPerPartition.begin(),site);
			else
				this->patternsPerPartition.push_back(site);
		}

		gapPatterns[gaps] += this->patternCounts[site];
	}

	size_t numTips = this->tau->getValue().getNumberOfTips();
	std::vector<TopologyNode*> nodes = this->tau->getValue().getNodes();

	size_t numCorrectionSites = 0;
	if(this->type & NO_ABSENT_SITES)
		numCorrectionSites++;
	if(this->type & NO_PRESENT_SITES)
		numCorrectionSites++;
	if(this->type & NO_SINGLETON_GAINS)
		numCorrectionSites += numTips;
	if(this->type & NO_SINGLETON_LOSSES)
		numCorrectionSites += numTips;

	this->numCorrectionPartitions = this->patternsPerPartition.size();

	for(std::map<std::string, size_t>::iterator it = gapPatterns.begin(); it != gapPatterns.end(); it++){
		size_t gaps = count(it->first.begin(),it->first.end(), '1');

		if(gaps > 0){
			this->numSitesPerPartition.push_back(it->second);
			this->numGapsPerPartition.push_back(gaps);

			size_t numSites = numCorrectionSites;
			if(this->type & NO_SINGLETON_GAINS)
				numSites -= gaps;
			if(this->type & NO_SINGLETON_LOSSES)
				numSites -= gaps;

			this->numCorrectionSitesPerPartition.push_back(numSites);
		}else{
			this->numSitesPerPartition.insert(this->numSitesPerPartition.begin(), it->second);
			this->numGapsPerPartition.insert(this->numGapsPerPartition.begin(), 0);
			this->numCorrectionSitesPerPartition.insert(this->numCorrectionSitesPerPartition.begin(), numCorrectionSites);
		}
	}

	this->perMixtureCorrections = std::vector<std::vector<double> >(this->numCorrectionPartitions,std::vector<double>(this->numSiteRates, 0.0));
	this->perPartitionLnCorrections = std::vector<double>(this->numCorrectionPartitions, 0.0);

	// allocate and fill the cells of the matrices

	this->correctionCharMatrix.clear();
	this->correctionGapMatrix.clear();

	size_t i = 0;
	for (std::vector<TopologyNode*>::iterator it = nodes.begin(); it != nodes.end(); ++it)
	{
		if ( (*it)->isTip() )
		{
			std::string name = (*it)->getName();

			this->correctionCharMatrix[name].resize(this->numCorrectionPartitions);
			this->correctionGapMatrix[name].resize(this->numCorrectionPartitions);

			for (size_t partition = 0; partition < this->numCorrectionPartitions; ++partition)
			{

				this->correctionCharMatrix[name][partition].resize(this->numCorrectionSitesPerPartition[partition]);
				this->correctionGapMatrix[name][partition].resize(this->numCorrectionSitesPerPartition[partition]);

				size_t present = 1 + this->usingAmbiguousCharacters;
				size_t absent = this->usingAmbiguousCharacters;

				size_t nonSingletons = bool(this->type & NO_ABSENT_SITES) + bool(this->type & NO_PRESENT_SITES);

				for (size_t site = 0; site < this->numCorrectionSitesPerPartition[partition]; site++)
				{
					this->correctionGapMatrix[name][partition][site] = this->gapMatrix[name][this->patternsPerPartition[partition]];

					if(this->correctionGapMatrix[name][partition][site]){
						//std::cerr << "gap" << std::endl;
						this->correctionCharMatrix[name][partition][site] = 0;
						continue;
					}

					if(site < nonSingletons){
						if(site == 0)
							this->correctionCharMatrix[name][partition][site] = (this->type & NO_ABSENT_SITES) ? absent : present;
						else
							this->correctionCharMatrix[name][partition][site] = present;
					}else{
						if(site - nonSingletons < numTips)
							if(site - nonSingletons == i){
								this->correctionCharMatrix[name][partition][site] = (this->type & NO_SINGLETON_LOSSES) ? absent : present;
							}else{
								this->correctionCharMatrix[name][partition][site] = (this->type & NO_SINGLETON_LOSSES) ? present : absent;
							}
						else
							if(site - nonSingletons == numTips + i){
								this->correctionCharMatrix[name][partition][site] = present;
							}else{
								this->correctionCharMatrix[name][partition][site] = absent;
							}
					}
				}
			}
			i++;
		}
	}

}

template<class treeType>
void RevBayesCore::BinaryPartitionModel<treeType>::redrawValue( void ) {

	delete this->value;

	// create a new character data object
	this->value = new DiscreteCharacterData<StandardState>();

	size_t numTips = this->tau->getValue().getNumberOfTips();
	size_t numNodes = this->tau->getValue().getNumberOfNodes();

	RandomNumberGenerator* rng = GLOBAL_RNG;

	const TopologyNode &root = this->tau->getValue().getRoot();
	size_t rootIndex = this->tau->getValue().getRoot().getIndex();

	std::vector< DiscreteTaxonData<StandardState> > taxa = std::vector< DiscreteTaxonData<StandardState> >(numTips, DiscreteTaxonData<StandardState>() );

	// first sample a mean number of characters (gamma)
	// from the marginal posterior gamma | N ~ Gamma(N, exp(lnCorrection) )
	// given gamma ~ 1/lambda
	if(this->type != NONE){
		double gamma = RbStatistics::Gamma::rv(this->numSites, exp(this->perPartitionLnCorrections.front()), *rng);

		// then resample numSites from Poisson( exp(lnCorrection)*gamma )
		this->numSites = RbStatistics::Poisson::rv( exp(this->perPartitionLnCorrections.front())*gamma, *rng);
	}

	// sample the rate categories conditioned the likelihood of the unobservable site-patterns
	std::vector<size_t> perSiteRates;
	double total = 0.0;
	if(this->type != NONE){
		for ( size_t i = 0; i < this->numSiteRates; ++i ){
			total += 1 - this->perMixtureCorrections.front()[i];
		}
	}else{
		total = this->numSiteRates;
	}

	for ( size_t i = 0; i < this->numSites; ++i )
	{
		// draw the state
		double u = rng->uniform01()*total;
		size_t rateIndex = 0;

		double tmp = 0.0;
		while(tmp < u){
			if(this->type != NONE)
				tmp += 1 - this->perMixtureCorrections.front()[rateIndex];
			else
				tmp += 1.0;

			if(tmp < u)
				rateIndex++;
		}
		perSiteRates.push_back( rateIndex );
	}

    const std::vector< double > &pi = this->getRootFrequencies();

    // then sample the site-pattern, rejecting those that match the unobservable ones
    for ( size_t i = 0; i < this->numSites; i++ )
    {
    	size_t rateIndex = perSiteRates[i];

        std::vector<StandardState> siteData(numNodes, StandardState());

        if(rng->uniform01() < pi[1])
        	siteData[rootIndex].setState(presentState);
        else
        	siteData[rootIndex].setState(absentState);

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
			taxa[t].addCharacter(siteData[t]);
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
size_t RevBayesCore::BinaryPartitionModel<treeType>::simulate( const TopologyNode &node, std::vector<StandardState> &data, size_t rateIndex) {

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
			childState.setState(presentState);
		else
			childState.setState(absentState);

		if(child.isTip())
		    numLeaves += childState.getStringValue() == this->presentState;
		else
			numLeaves += simulate( child, data, rateIndex );
	}

    return numLeaves;

}

template<class treeType>
RevBayesCore::DiscreteTaxonData<RevBayesCore::StandardState>* RevBayesCore::BinaryPartitionModel<treeType>::getDolloCompatible(void) {

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
		if(rootState.getStringValue() == this->absentState && leftState.getStringValue() == this->presentState)
			births++;
		if(rootState.getStringValue() == this->absentState && rightState.getStringValue() == this->presentState)
			births++;

    	if(!left.isTip())
    		births += numBirths(i,mapping,left);
    	if(!right.isTip())
    		births += numBirths(i,mapping,right);

    	StandardState c;
    	if(births <= 1)
    		c.setState("1");
    	else
    		c.setState("0");

    	dolloCompatible->addCharacter(c);
    }

    return dolloCompatible;
}

template<class treeType>
void RevBayesCore::BinaryPartitionModel<treeType>::touchSpecialization( DagNode* affecter ) {

    // if the topology wasn't the culprit for the touch, then we just flag everything as dirty
    if ( affecter == this->dagNode )
    {
    	if(((DiscreteCharacterData<StandardState> *)this->value)->getNumberOfCharacters() > 0){
			std::string symbols = ((DiscreteCharacterData<StandardState> *)this->value)->getCharacter(0,0).getStateLabels();

			presentState = symbols[1];
			absentState = symbols[0];
    	}

	}

    PartitionedSiteCorrectionModel<StandardState,treeType>::touchSpecialization(affecter);

}

template<class treeType>
size_t RevBayesCore::BinaryPartitionModel<treeType>::numBirths(size_t site, const std::vector<DiscreteTaxonData<StandardState> > &mapping, const TopologyNode &node) {

	const TopologyNode &left = node.getChild(0);
	const TopologyNode &right = node.getChild(1);

	StandardState nodeState = mapping[node.getIndex()].getCharacter(site);
	StandardState leftState = mapping[left.getIndex()].getCharacter(site);
	StandardState rightState = mapping[right.getIndex()].getCharacter(site);

	size_t births = 0;
	if(nodeState.getStringValue() == this->absentState && leftState.getStringValue() == this->presentState)
		births++;
	if(nodeState.getStringValue() == this->absentState && rightState.getStringValue() == this->presentState)
		births++;

	if(!left.isTip())
		births += numBirths(site,mapping,left);
	if(!right.isTip())
		births += numBirths(site,mapping,right);

	return births;

}


#endif

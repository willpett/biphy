#ifndef BinaryCharEvoModel_H
#define BinaryCharEvoModel_H

#include "AbstractSiteCorrectionModel.h"
#include "StandardState.h"
#include "RateMatrix.h"
#include "RbVector.h"
#include "RbStatisticsHelper.h"
#include "TopologyNode.h"
#include "TransitionProbabilityMatrix.h"
#include "TreeChangeEventListener.h"
#include "TypedDistribution.h"

#include "DistributionBinomial.h"
#include "DistributionNegativeBinomial.h"

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
        
        virtual void                        				computeRootCorrection(size_t root, size_t l, size_t r);
		virtual void                        				computeInternalNodeCorrection(const TopologyNode &n, size_t nIdx, size_t l, size_t r);
		virtual void                        				computeTipCorrection(const TopologyNode &node, size_t nIdx);

        void                                        		touchSpecialization(DagNode *toucher);
        virtual void                                        setCorrectionPatterns();
        size_t 												simulate( const TopologyNode &node, std::vector<StandardState> &taxa, size_t rateIndex);
        
        size_t												numBirths(size_t site, const std::vector<DiscreteTaxonData<StandardState> > &mapping, const TopologyNode &node);

        int													type;
        DiscreteTaxonData<StandardState> *					dolloCompatible;

        std::string											presentState;
        std::string											absentState;

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
	this->numCorrectionSites = 4*(type > 0);

	this->resizeLikelihoodVectors();

	presentState = "1";
	absentState = "0";

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
	this->correctionCharMatrix.clear();
	this->correctionGapMatrix.clear();
}

template<class treeType>
void RevBayesCore::BinaryCharEvoModel<treeType>::redrawValue( void ) {

	//this->computeLnProbability();

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


    // first sample a total number of characters (M)
    // from the marginal posterior M - N | N ~ NegBinomial(N+1, exp(lnCorrection) )
    if(this->numCorrectionSites > 0){
		//double gamma = RbStatistics::Gamma::rv(this->N,exp(this->perSiteCorrection), *rng);
    	double M = RbStatistics::NegativeBinomial::rv(this->N + 1,exp(this->perSiteCorrection), *rng);

		// resample numSites from Binomial( M + N, exp(lnCorrection) )
		this->numSites = RbStatistics::Binomial::rv( M + this->N, exp(this->perSiteCorrection), *rng);
    }

    // sample the rate categories conditioned the likelihood of the unobservable site-patterns
    std::vector<size_t> perSiteRates;
	double total = 0.0;
	for ( size_t i = 0; i < this->numSiteRates; ++i ){
		total += 1 - this->perMixtureCorrections[i];
	}

    for ( size_t i = 0; i < this->numSites; ++i )
	{
		// draw the state
		double u = rng->uniform01()*total;
		size_t rateIndex = 0;

		double tmp = 0.0;
		while(tmp < u){
			tmp += 1 - this->perMixtureCorrections[rateIndex];
			if(tmp < u)
				rateIndex++;
		}
		perSiteRates.push_back( rateIndex );
	}

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
void RevBayesCore::BinaryCharEvoModel<treeType>::computeTipCorrection(const TopologyNode &node, size_t nodeIndex)
{
	if(this->numCorrectionSites == 0)
		return;

	std::vector<double>::iterator p_node = this->partialLikelihoods.begin() + this->activeLikelihood[nodeIndex]*this->activeLikelihoodOffset + nodeIndex*this->nodeOffset;

	// iterate over all mixture categories

	for (size_t mixture = 0; mixture < this->numSiteRates; ++mixture)
	{

		size_t offset = mixture*this->mixtureOffset + this->numPatterns*this->siteOffset;

		std::vector<double>::iterator     	u      = p_node + offset;

		//Probability of no observed presence in all leaves descending from this node
		u[0] = 1.0;
		u[1] = 0.0;

		//Probability of observed presence in only one leaf descending from this node
		u[2] = 0.0;
		u[3] = 1.0;

		//Probability of observed presence in all leaves descending from this node
		u[4] = 0.0;
		u[5] = 1.0;

		//Probability of observed presence in all but one leaves descending from this node
		u[6] = 1.0;
		u[7] = 0.0;

	}
}

template<class treeType>
void RevBayesCore::BinaryCharEvoModel<treeType>::computeInternalNodeCorrection(const TopologyNode &node, size_t nodeIndex, size_t left, size_t right)
{

	if(this->numCorrectionSites == 0)
		return;

	// get the pointers to the partial likelihoods for this node and the two descendant subtrees
	std::vector<double>::const_iterator   p_left  = this->partialLikelihoods.begin() + this->activeLikelihood[left]*this->activeLikelihoodOffset + left*this->nodeOffset;
	std::vector<double>::const_iterator   p_right = this->partialLikelihoods.begin() + this->activeLikelihood[right]*this->activeLikelihoodOffset + right*this->nodeOffset;
	std::vector<double>::iterator         p_node  = this->partialLikelihoods.begin() + this->activeLikelihood[nodeIndex]*this->activeLikelihoodOffset + nodeIndex*this->nodeOffset;

	const TopologyNode &left_node = this->tau->getValue().getNode(left);
	const TopologyNode &right_node = this->tau->getValue().getNode(right);

	// iterate over all mixture categories

	for (size_t mixture = 0; mixture < this->numSiteRates; ++mixture)
	{
		this->updateTransitionProbabilities( left, left_node.getBranchLength() );
		TransitionProbabilityMatrix	pij    		= this->transitionProbMatrices[mixture];

		this->updateTransitionProbabilities( right, right_node.getBranchLength() );
		TransitionProbabilityMatrix	pik    		= this->transitionProbMatrices[mixture];

		size_t offset = mixture*this->mixtureOffset + this->numPatterns*this->siteOffset;

		std::vector<double>::iterator     			u       = p_node + offset;
		std::vector<double>::const_iterator     	u_j  	= p_left + offset;
		std::vector<double>::const_iterator     	u_k 	= p_right + offset;

		for(size_t site = 0; site < 8; site += 4){
			u[site] = 0;
			u[site + 1] = 0;
			u[site + 2] = 0;
			u[site + 3] = 0;
			for(size_t cj = 0; cj < 2; cj++){
				for(size_t ck = 0; ck < 2; ck++){
					u[site] 	+= pij[0][cj]*u_j[site + cj]*pik[0][ck]*u_k[site + ck];
					u[site + 1] += pij[1][cj]*u_j[site + cj]*pik[1][ck]*u_k[site + ck];

					u[site + 2] += pij[0][cj]*u_j[site + cj]*pik[0][ck]*u_k[site + 2 + ck] + pij[0][cj]*u_j[site + 2 + cj]*pik[0][ck]*u_k[site + ck];
					u[site + 3] += pij[1][cj]*u_j[site + cj]*pik[1][ck]*u_k[site + 2 + ck] + pij[1][cj]*u_j[site + 2 + cj]*pik[1][ck]*u_k[site + ck];
				}
			}
		}
	}

}

template<class treeType>
void RevBayesCore::BinaryCharEvoModel<treeType>::computeRootCorrection( size_t root, size_t left, size_t right)
{
	this->perSiteCorrection = 0.0;

	if(this->numCorrectionSites > 0){

		const std::vector<double> &f                    = this->getRootFrequencies();

		// get the pointers to the partial likelihoods for this node and the two descendant subtrees
		std::vector<double>::const_iterator   p_left  = this->partialLikelihoods.begin() + this->activeLikelihood[left]*this->activeLikelihoodOffset + left*this->nodeOffset;
		std::vector<double>::const_iterator   p_right = this->partialLikelihoods.begin() + this->activeLikelihood[right]*this->activeLikelihoodOffset + right*this->nodeOffset;
		std::vector<double>::iterator         p_node  = this->partialLikelihoods.begin() + this->activeLikelihood[root]*this->activeLikelihoodOffset + root*this->nodeOffset;

		const TopologyNode &left_node = this->tau->getValue().getNode(left);
		const TopologyNode &right_node = this->tau->getValue().getNode(right);

		// iterate over all mixture categories
		for (size_t mixture = 0; mixture < this->numSiteRates; ++mixture)
		{
			this->perMixtureCorrections[mixture] = 0.0;

			//get the 1->1 transition probabilities for each branch
			this->updateTransitionProbabilities( left, left_node.getBranchLength() );
			TransitionProbabilityMatrix	pij    		= this->transitionProbMatrices[mixture];

			this->updateTransitionProbabilities( right, right_node.getBranchLength() );
			TransitionProbabilityMatrix	pik    		= this->transitionProbMatrices[mixture];

			size_t offset = mixture*this->mixtureOffset + this->numPatterns*this->siteOffset;

			std::vector<double>::iterator     			u       = p_node + offset;
			std::vector<double>::const_iterator     	u_j  	= p_left + offset;
			std::vector<double>::const_iterator     	u_k 	= p_right + offset;

			for(size_t site = 0; site < 8; site += 4){
				u[site] = 0;
				u[site + 1] = 0;
				u[site + 2] = 0;
				u[site + 3] = 0;
				for(size_t cj = 0; cj < 2; cj++){
					for(size_t ck = 0; ck < 2; ck++){
						u[site] 	+= pij[0][cj]*u_j[site + cj]*pik[0][ck]*u_k[site + ck];
						u[site + 1] += pij[1][cj]*u_j[site + cj]*pik[1][ck]*u_k[site + ck];

						u[site + 2] += pij[0][cj]*u_j[site + cj]*pik[0][ck]*u_k[site + 2 + ck] + pij[0][cj]*u_j[site + 2 + cj]*pik[0][ck]*u_k[site + ck];
						u[site + 3] += pij[1][cj]*u_j[site + cj]*pik[1][ck]*u_k[site + 2 + ck] + pij[1][cj]*u_j[site + 2 + cj]*pik[1][ck]*u_k[site + ck];
					}
				}
			}

			bool jl = left_node.isTip();
			bool kl = right_node.isTip();


			//std::cerr << u[0] << "\t" << u[1] << std::endl;

			double prob = 1.0;
			if(this->type & NO_ABSENT_SITES)
				prob -= u[0]*f[0] + u[1]*f[1];
			if(this->type & NO_SINGLETON_GAINS)
				prob -= u[2]*f[0] + u[3]*f[1];
			if(this->type & NO_PRESENT_SITES)
				prob -= u[4]*f[0] + u[5]*f[1];
			//If both of this node's children are leaves, then u[1] = u[3]
			if((this->type & NO_SINGLETON_LOSSES) && !(jl && kl))
				prob -= u[6]*f[0] + u[7]*f[1];

			// correct rounding errors
			if(prob < 0)
				prob = 0;

			//std::cerr << prob << std::endl;

			this->perMixtureCorrections[mixture] = prob;
		}

		// sum the likelihoods for all mixtures together
		for (size_t mixture = 0; mixture < this->numSiteRates; ++mixture)
			this->perSiteCorrection += this->perMixtureCorrections[mixture];

		// normalize the log-probability
		this->perSiteCorrection = log(this->perSiteCorrection) - log(this->numSiteRates);
	}
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
void RevBayesCore::BinaryCharEvoModel<treeType>::touchSpecialization( DagNode* affecter ) {

    // if the topology wasn't the culprit for the touch, then we just flag everything as dirty
    if ( affecter == this->dagNode )
    {
    	if(((DiscreteCharacterData<StandardState> *)this->value)->getNumberOfCharacters() > 0){
			std::string symbols = ((DiscreteCharacterData<StandardState> *)this->value)->getCharacter(0,0).getStateLabels();

			presentState = symbols[1];
			absentState = symbols[0];
    	}

	}else
	{
		AbstractSiteCorrectionModel<StandardState,treeType>::touchSpecialization(affecter);
    }

}

template<class treeType>
size_t RevBayesCore::BinaryCharEvoModel<treeType>::numBirths(size_t site, const std::vector<DiscreteTaxonData<StandardState> > &mapping, const TopologyNode &node) {

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

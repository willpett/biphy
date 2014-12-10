#ifndef PartitionedSiteCorrectionModel_H
#define PartitionedSiteCorrectionModel_H

#include "AbstractCharacterData.h"
#include "DiscreteTaxonData.h"
#include "DistributionPoisson.h"
#include "DistributionGamma.h"
#include "GeneralCharEvoModel.h"
#include "RateMatrix.h"
#include "RbStatisticsHelper.h"
#include "TopologyNode.h"
#include "TransitionProbabilityMatrix.h"
#include "Tree.h"
#include "TreeChangeEventListener.h"
#include "TypedDistribution.h"

#include <memory.h>
#include <algorithm>

namespace RevBayesCore {

	enum CorrectionType { NONE 					= 0x00,
					  	  NO_CONSTANT_SITES		= 0x03,
						  NO_UNINFORMATIVE		= 0x0F
						};

    template<class charType, class treeType>
    class PartitionedSiteCorrectionModel : public GeneralCharEvoModel<charType, treeType> {
        
    public:
        // Note, we need the size of the alignment in the constructor to correctly simulate an initial state
        PartitionedSiteCorrectionModel(const TypedDagNode< treeType > *t, size_t nChars, bool c, size_t nSites, int type = 0);
        PartitionedSiteCorrectionModel(const PartitionedSiteCorrectionModel &n);                                                                                          //!< Copy constructor
        virtual                                                            ~PartitionedSiteCorrectionModel(void);
        
        virtual void                                       					redrawValue(void);

    protected:
        // helper method for this and derived classes
        virtual void                                        				computeRootLikelihood(size_t root, size_t l, size_t r);
		virtual void                                        				computeInternalNodeLikelihood(const TopologyNode &n, size_t nIdx, size_t l, size_t r);
		virtual void                                        				computeTipLikelihood(const TopologyNode &node, size_t nIdx);
		virtual void                                        				resizeLikelihoodVectors(void);
        
		virtual void                                        				touchSpecialization(DagNode *toucher);

        virtual void                                                        computeRootCorrections(size_t root, size_t l, size_t r);
        virtual void                                                        computeInternalNodeCorrections(const TopologyNode &n, size_t nIdx, size_t l, size_t r);
        virtual void                                                        computeTipCorrections(const TopologyNode &node, size_t nIdx);
        
        virtual void                                                        setCorrectionPatterns();

        std::map<std::string,std::vector<std::vector<unsigned long> > >     correctionCharMatrix;
        std::map<std::string,std::vector<std::vector<bool> > >              correctionGapMatrix;

        size_t																numCorrectionPartitions;
        std::vector<size_t>													numCorrectionSitesPerPartition;
        std::vector<size_t>													patternsPerPartition;
        std::vector<size_t>													numSitesPerPartition;
        std::vector<size_t>													numGapsPerPartition;

        std::vector<std::vector<double> >									perMixtureCorrections;
        std::vector<double>													perPartitionLnCorrections;

        int																	type;

    private:
        
        std::map<size_t, size_t> 											simulate( const TopologyNode &node, std::vector<charType> &taxa, size_t rateIndex);
    
    };
    
}


#include "DiscreteCharacterState.h"
#include "DiscreteCharacterData.h"
#include "RandomNumberFactory.h"
#include "RandomNumberGenerator.h"
#include "TopologyNode.h"
#include "TransitionProbabilityMatrix.h"

#include <cmath>

template<class charType, class treeType>
RevBayesCore::PartitionedSiteCorrectionModel<charType, treeType>::PartitionedSiteCorrectionModel(const TypedDagNode<treeType> *t, size_t nChars, bool c, size_t nSites, int intype) :
	GeneralCharEvoModel<charType, treeType>(t, nChars, c, nSites),
	AbstractCharEvoModel<charType, treeType>(t, nChars, c, nSites),
	numCorrectionPartitions(),
	numCorrectionSitesPerPartition(),
	numSitesPerPartition(),
	perMixtureCorrections(),
	perPartitionLnCorrections(),
	type(intype)
{
}


template<class charType, class treeType>
RevBayesCore::PartitionedSiteCorrectionModel<charType, treeType>::PartitionedSiteCorrectionModel(const PartitionedSiteCorrectionModel &n) :
	GeneralCharEvoModel<charType, treeType>( n ),
	AbstractCharEvoModel<charType, treeType>( n ),
	correctionCharMatrix(n.correctionCharMatrix),
	correctionGapMatrix(n.correctionGapMatrix),
	numCorrectionPartitions(n.numCorrectionPartitions),
	numCorrectionSitesPerPartition(n.numCorrectionSitesPerPartition),
	numSitesPerPartition(n.numSitesPerPartition),
	perMixtureCorrections(n.perMixtureCorrections),
	numGapsPerPartition(n.numGapsPerPartition),
	patternsPerPartition(n.patternsPerPartition),
	perPartitionLnCorrections(n.perPartitionLnCorrections),
	type(n.type)
{
	/*
	for(std::map<std::string,std::vector<std::vector<bool> > >::iterator it = correctionGapMatrix.begin(); it != correctionGapMatrix.end(); it++){
		std::cerr << it->first << std::endl;
		for(size_t p = 0; p < numCorrectionPartitions; p++){
			for(size_t site = 0; site < numCorrectionSitesPerPartition[p]; site++){
				std::cerr << (it->second)[p][site] << "\t";
			}
		}
		std::cerr << std::endl;
	}
	*/
}


template<class charType, class treeType>
RevBayesCore::PartitionedSiteCorrectionModel<charType, treeType>::~PartitionedSiteCorrectionModel( void ) {
    // We don't delete the paramoms, because they might be used somewhere else too. The model needs to do that!
}

template<class charType, class treeType>
void RevBayesCore::PartitionedSiteCorrectionModel<charType, treeType>::computeRootLikelihood(size_t root, size_t left, size_t right)
{

    GeneralCharEvoModel<charType, treeType>::computeRootLikelihood(root, left, right);

    computeRootCorrections(root, left, right);

    for(size_t partition = 0; partition < numCorrectionPartitions; partition++){
    	//std::cerr << numSitesPerPartition[partition] << "\t" << perPartitionLnCorrections[partition] << std::endl;
    	this->lnProb -= numSitesPerPartition[partition]*perPartitionLnCorrections[partition];
    }
}

template<class charType, class treeType>
void RevBayesCore::PartitionedSiteCorrectionModel<charType, treeType>::computeInternalNodeLikelihood(const TopologyNode &node, size_t nodeIndex, size_t left, size_t right)
{

    GeneralCharEvoModel<charType, treeType>::computeInternalNodeLikelihood(node, nodeIndex, left, right);

    computeInternalNodeCorrections(node, nodeIndex, left, right);

}

template<class charType, class treeType>
void RevBayesCore::PartitionedSiteCorrectionModel<charType, treeType>::computeTipLikelihood(const TopologyNode &node, size_t nIdx)
{

    GeneralCharEvoModel<charType, treeType>::computeTipLikelihood(node, nIdx);

    computeTipCorrections(node, nIdx);

}

template<class charType, class treeType>
void RevBayesCore::PartitionedSiteCorrectionModel<charType, treeType>::computeRootCorrections( size_t root, size_t left, size_t right)
{
	if(type == NONE)
		return;

    // get the root frequencies
    const std::vector<double> &f                    = this->getRootFrequencies();
    std::vector<double>::const_iterator f_end       = f.end();
    std::vector<double>::const_iterator f_begin     = f.begin();

    // get the pointers to the partial likelihoods of the left and right subtree
    std::vector<double>::const_iterator p_left   = this->partialLikelihoods.begin() + this->activeLikelihood[left]*this->activeLikelihoodOffset + left*this->nodeOffset;
    std::vector<double>::const_iterator p_right  = this->partialLikelihoods.begin() + this->activeLikelihood[right]*this->activeLikelihoodOffset + right*this->nodeOffset;

	// get pointers the likelihood for both subtrees
	std::vector<double>::const_iterator   p_partition_left     = p_left + this->numPatterns*this->siteOffset;
	std::vector<double>::const_iterator   p_partition_right    = p_right + this->numPatterns*this->siteOffset;

	for (size_t partition = 0; partition < numCorrectionPartitions; ++partition)
	{
		perPartitionLnCorrections[partition] = 0.0;

		std::vector<double>::const_iterator   p_mixture_partition_left     = p_partition_left;
		std::vector<double>::const_iterator   p_mixture_partition_right    = p_partition_right;

		for (size_t mixture = 0; mixture < this->numSiteRates; ++mixture)
		{
			perMixtureCorrections[partition][mixture] = 0.0;

			std::vector<double>::const_iterator   p_site_mixture_partition_left     = p_mixture_partition_left;
			std::vector<double>::const_iterator   p_site_mixture_partition_right    = p_mixture_partition_right;

			for (size_t site = 0; site < numCorrectionSitesPerPartition[partition]; ++site)
			{
				// temporary variable storing the likelihood
				double tmp = 0.0;
				// get the pointer to the stationary frequencies
				std::vector<double>::const_iterator f_j             = f_begin;
				// get the pointers to the likelihoods for this site and mixture category
				std::vector<double>::const_iterator p_site_left_j   = p_site_mixture_partition_left;
				std::vector<double>::const_iterator p_site_right_j  = p_site_mixture_partition_right;
				// iterate over all starting states
				for (; f_j != f_end; ++f_j)
				{
					// add the probability of starting from this state
					tmp += *p_site_left_j * *p_site_right_j * *f_j;

					// increment pointers
					++p_site_left_j; ++p_site_right_j;
				}
				// add the likelihood for this mixture category
				perMixtureCorrections[partition][mixture] += tmp;

				// increment the pointers to the next site
				p_site_mixture_partition_left+=this->siteOffset; p_site_mixture_partition_right+=this->siteOffset;

			} // end-for over all sites (=patterns)

			p_mixture_partition_left += this->mixtureOffset; p_mixture_partition_right += this->mixtureOffset;

			perPartitionLnCorrections[partition] += 1.0 - perMixtureCorrections[partition][mixture];
		}

		// increment the pointers to the next mixture category
		p_partition_left += numCorrectionSitesPerPartition[partition]*this->siteOffset;
		p_partition_right += numCorrectionSitesPerPartition[partition]*this->siteOffset;

		// normalize the log-probability
		perPartitionLnCorrections[partition] = log(perPartitionLnCorrections[partition]) - log(this->numSiteRates);

    }

}

template<class charType, class treeType>
void RevBayesCore::PartitionedSiteCorrectionModel<charType,treeType>::computeInternalNodeCorrections(const TopologyNode &node, size_t nodeIndex, size_t left, size_t right)
{
	if(type == NONE)
		return;

	// compute the transition probability matrix
	this->updateTransitionProbabilities( nodeIndex, node.getBranchLength() );

	// get the pointers to the partial likelihoods for this node and the two descendant subtrees
	std::vector<double>::const_iterator   p_left  = this->partialLikelihoods.begin() + this->activeLikelihood[left]*this->activeLikelihoodOffset + left*this->nodeOffset;
	std::vector<double>::const_iterator   p_right = this->partialLikelihoods.begin() + this->activeLikelihood[right]*this->activeLikelihoodOffset + right*this->nodeOffset;
	std::vector<double>::iterator         p_node  = this->partialLikelihoods.begin() + this->activeLikelihood[nodeIndex]*this->activeLikelihoodOffset + nodeIndex*this->nodeOffset;

	// iterate over all mixture categories
	for (size_t mixture = 0; mixture < this->numSiteRates; ++mixture)
	{
		// the transition probability matrix for this mixture category
		const double *                       tp_begin    = this->transitionProbMatrices[mixture].theMatrix;

		// get the pointers to the likelihood for this mixture category
		size_t offset 	= mixture*this->mixtureOffset + this->numPatterns*this->siteOffset;

		for (size_t partition = 0; partition < numCorrectionPartitions; ++partition)
		{
			std::vector<double>::iterator          p_site_mixture          = p_node + offset;
			std::vector<double>::const_iterator    p_site_mixture_left     = p_left + offset;
			std::vector<double>::const_iterator    p_site_mixture_right    = p_right + offset;
			// compute the per site probabilities

			for (size_t site = 0; site < numCorrectionSitesPerPartition[partition] ; ++site)
			{
				// get the pointers for this mixture category and this site
				const double*       tp_a    = tp_begin;
				// iterate over the possible starting states
				for (size_t c1 = 0; c1 < this->numChars; ++c1)
				{
					// temporary variable
					double sum = 0.0;

					// iterate over all possible terminal states
					for (size_t c2 = 0; c2 < this->numChars; ++c2 )
					{
						sum += p_site_mixture_left[c2] * p_site_mixture_right[c2] * tp_a[c2];

					} // end-for over all distination character

					// store the likelihood for this starting state
					p_site_mixture[c1] = sum;

					// increment the pointers to the next starting state
					tp_a+=this->numChars;

				} // end-for over all initial characters

				// increment the pointers to the next site
				p_site_mixture_left+=this->siteOffset; p_site_mixture_right+=this->siteOffset; p_site_mixture+=this->siteOffset;

			} // end-for over all sites (=patterns)

			offset += numCorrectionSitesPerPartition[partition]*this->siteOffset;
		}

	} // end-for over all mixtures (=rate-categories)

}

template<class charType, class treeType>
void RevBayesCore::PartitionedSiteCorrectionModel<charType, treeType>::computeTipCorrections(const TopologyNode &node, size_t nodeIndex)
{
	if(type == NONE)
		return;

	std::vector<double>::iterator p_node = this->partialLikelihoods.begin() + this->activeLikelihood[nodeIndex]*this->activeLikelihoodOffset + nodeIndex*this->nodeOffset;

    std::string name = node.getName();

    const std::vector<std::vector<bool> > &gap_node = correctionGapMatrix[name];
    const std::vector<std::vector<unsigned long> > &char_node = correctionCharMatrix[name];

    // compute the transition probabilities
    this->updateTransitionProbabilities( nodeIndex, node.getBranchLength() );

    std::vector<double>::iterator   p_mixture      = p_node;

    // iterate over all mixture categories
    for (size_t mixture = 0; mixture < this->numSiteRates; ++mixture)
    {
        // the transition probability matrix for this mixture category
    	const double *	tp_begin    		= this->transitionProbMatrices[mixture].theMatrix;

    	size_t offset = mixture*this->mixtureOffset + this->numPatterns*this->siteOffset;

    	for (size_t partition = 0; partition < numCorrectionPartitions; ++partition)
    	{
			// iterate over all sites
			for (size_t site = 0; site != numCorrectionSitesPerPartition[partition]; ++site)
			{

				std::vector<double>::iterator     	p_site_mixture      = p_mixture + offset + site*this->siteOffset;
				// is this site a gap?
				if ( gap_node[partition][site] )
				{
					// since this is a gap we need to assume that the actual state could have been any state
					// iterate over all initial states for the transitions
					for (size_t c1 = 0; c1 < this->numChars; ++c1)
					{

						// store the likelihood
						p_site_mixture[c1] = 1.0;

					}
				}
				else // we have observed a character
				{

					// get the original character
					unsigned long org_val = char_node[partition][site];

					// iterate over all possible initial states
					for (size_t c1 = 0; c1 < this->numChars; ++c1)
					{

						if ( this->usingAmbiguousCharacters )
						{
							// compute the likelihood that we had a transition from state c1 to the observed state org_val
							// note, the observed state could be ambiguous!
							unsigned long val = org_val;

							// get the pointer to the transition probabilities for the terminal states
							const double * d = tp_begin+(this->numChars*c1);

							double tmp = 0.0;

							while ( val != 0 ) // there are still observed states left
							{
								// check whether we observed this state
								if ( (val & 1) == 1 )
								{
									// add the probability
									tmp += *d;
								}

								// remove this state from the observed states
								val >>= 1;

								// increment the pointer to the next transition probability
								++d;
							} // end-while over all observed states for this character

							// store the likelihood
							p_site_mixture[c1] = tmp;

						}
						else // no ambiguous characters in use
						{

							// store the likelihood
							p_site_mixture[c1] = tp_begin[c1*this->numChars+org_val];
							//l+= org_val;

						}

					} // end-for over all possible initial character for the branch

				}

			} // end-for over sites

			offset += numCorrectionSitesPerPartition[partition]*this->siteOffset;

    	} // end-for over partitions

    } // end-for over mixtures
}

template<class charType, class treeType>
void RevBayesCore::PartitionedSiteCorrectionModel<charType, treeType>::touchSpecialization( DagNode* affecter ) {

	GeneralCharEvoModel<charType, treeType>::touchSpecialization(affecter);

	if(affecter == this->dagNode){
		setCorrectionPatterns();
		resizeLikelihoodVectors();
	}

	if(affecter == this->siteRates)
		perMixtureCorrections = std::vector<std::vector<double> >(this->numCorrectionPartitions,std::vector<double>(this->numSiteRates, 0.0));
}

template<class charType, class treeType>
void RevBayesCore::PartitionedSiteCorrectionModel<charType, treeType>::setCorrectionPatterns(){

	numCorrectionPartitions = 0;
	numSitesPerPartition.clear();
	numGapsPerPartition.clear();
	numCorrectionSitesPerPartition.clear();
	perMixtureCorrections.clear();
	perPartitionLnCorrections.clear();
	patternsPerPartition.clear();

	if(type == NONE)
		return;

	std::map<std::string, size_t> gapPatterns;
	for(size_t site = 0; site < this->numPatterns; site++){
		std::stringstream gapstream;
		for(std::map<std::string, std::vector<bool> >::iterator it = this->gapMatrix.begin(); it != this->gapMatrix.end(); it++)
			gapstream << it->second[site];

		std::string gaps = gapstream.str();

		if(gapPatterns.find(gaps) == gapPatterns.end()){
			if(std::count(gaps.begin(), gaps.end(), '1') == 0)
				patternsPerPartition.insert(patternsPerPartition.begin(),site);
			else
				patternsPerPartition.push_back(site);
		}

		gapPatterns[gaps] += this->patternCounts[site];
	}

	size_t numTips = this->tau->getValue().getNumberOfTips();
	std::vector<TopologyNode*> nodes = this->tau->getValue().getNodes();

	size_t numCorrectionSites = 0;
	if(type & NO_CONSTANT_SITES)
		numCorrectionSites += this->numChars;
	if(type & NO_UNINFORMATIVE)
		numCorrectionSites += numTips*this->numChars*(this->numChars - 1)/2;

	numCorrectionPartitions = patternsPerPartition.size();

	for(std::map<std::string, size_t>::iterator it = gapPatterns.begin(); it != gapPatterns.end(); it++){

		size_t gaps = count(it->first.begin(),it->first.end(), '1');

		if(gaps > 0){
			numSitesPerPartition.push_back(it->second);
			numGapsPerPartition.push_back(gaps);

			if(type & NO_UNINFORMATIVE)
				numCorrectionSitesPerPartition.push_back(numCorrectionSites - gaps*this->numChars*(this->numChars - 1)/2);
			else
				numCorrectionSitesPerPartition.push_back(numCorrectionSites);
		}else{
			numSitesPerPartition.insert(numSitesPerPartition.begin(), it->second);
			numGapsPerPartition.insert(numGapsPerPartition.begin(), 0);
			numCorrectionSitesPerPartition.insert(numCorrectionSitesPerPartition.begin(), numCorrectionSites);
		}
	}

	perMixtureCorrections = std::vector<std::vector<double> >(numCorrectionPartitions,std::vector<double>(this->numSiteRates, 0.0));
	perPartitionLnCorrections = std::vector<double>(numCorrectionPartitions, 0.0);

	// allocate and fill the cells of the matrices

	correctionCharMatrix.clear();
	correctionGapMatrix.clear();

	size_t i = 0;
	for (std::vector<TopologyNode*>::iterator it = nodes.begin(); it != nodes.end(); ++it)
	{
		if ( (*it)->isTip() )
		{
			std::string name = (*it)->getName();

			correctionCharMatrix[name].resize(numCorrectionPartitions);
			correctionGapMatrix[name].resize(numCorrectionPartitions);

			for (size_t partition = 0; partition < numCorrectionPartitions; ++partition)
			{
				correctionCharMatrix[name][partition].resize(numCorrectionSitesPerPartition[partition]);
				correctionGapMatrix[name][partition].resize(numCorrectionSitesPerPartition[partition]);

				size_t char1 = 0;
				size_t char2 = 1;

				size_t nonSingletons = bool(type & NO_CONSTANT_SITES)*this->numChars;
				for (size_t site = 0; site < numCorrectionSitesPerPartition[partition]; site++)
				{
					correctionGapMatrix[name][partition][site] = this->gapMatrix[name][patternsPerPartition[partition]];

					if(correctionGapMatrix[name][partition][site]){
						correctionCharMatrix[name][partition][site] = 0;
						continue;
					}

					if(site < nonSingletons){
						correctionCharMatrix[name][partition][site] = 1 + site + this->usingAmbiguousCharacters;
					}else{
						if((site - nonSingletons) % numGapsPerPartition[partition] == i){
							correctionCharMatrix[name][partition][site] = char1 + 1 + this->usingAmbiguousCharacters;
						}else{
							correctionCharMatrix[name][partition][site] = char2 + 1 + this->usingAmbiguousCharacters;
						}
					}

					char2++;

					if(char2 == char1){
						char2++;
					}

					if(char2 == this->numChars){
						char1++;
						char2 = 0;
					}
				}
			}
			i++;
		}
	}
}

template<class charType, class treeType>
void RevBayesCore::PartitionedSiteCorrectionModel<charType, treeType>::redrawValue( void ) {

	//this->computeLnProbability();

    // delete the old value first
    delete this->value;

    // create a new character data object
    this->value = new DiscreteCharacterData<charType>();

    size_t numTips = this->tau->getValue().getNumberOfTips();
    size_t numNodes = this->tau->getValue().getNumberOfNodes();

    RandomNumberGenerator* rng = GLOBAL_RNG;

    const TopologyNode &root = this->tau->getValue().getRoot();
    size_t rootIndex = this->tau->getValue().getRoot().getIndex();

    std::vector< DiscreteTaxonData<charType> > taxa = std::vector< DiscreteTaxonData<charType> >(numTips, DiscreteTaxonData<charType>() );

    // first sample a mean number of characters (gamma)
    // from the marginal posterior gamma | N ~ Gamma(N, exp(lnCorrection) )
    // given gamma ~ 1/lambda
    if(type != NONE){
		double gamma = RbStatistics::Gamma::rv(this->numSites, exp(perPartitionLnCorrections.front()), *rng);

		// then resample numSites from Poisson( exp(lnCorrection)*gamma )
		this->numSites = RbStatistics::Poisson::rv( exp(perPartitionLnCorrections.front())*gamma, *rng);
    }

    // sample the rate categories conditioned the likelihood of the unobservable site-patterns
    std::vector<size_t> perSiteRates;
	double total = 0.0;
	if(type != NONE){
		for ( size_t i = 0; i < this->numSiteRates; ++i ){
			total += 1 - perMixtureCorrections.front()[i];
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
			if(type != NONE)
				tmp += 1 - perMixtureCorrections.front()[rateIndex];
			else
				tmp += 1.0;

			if(tmp < u)
				rateIndex++;
		}
		perSiteRates.push_back( rateIndex );
	}

    const std::vector< double > &stationaryFreqs = this->getRootFrequencies();

    // then sample the site-pattern, rejecting those that match the unobservable ones
    for ( size_t i = 0; i < this->numSites; i++ )
    {
    	size_t rateIndex = perSiteRates[i];

        std::vector<charType> siteData(numNodes, charType());

        // create the character
		charType c = siteData[rootIndex];
		c.setToFirstState();
		// draw the state
		double u = rng->uniform01();
		std::vector< double >::const_iterator freq = stationaryFreqs.begin();
		while ( true )
		{
			u -= *freq;

			if ( u > 0.0 )
			{
				++c;
				++freq;
			}
			else
			{
				break;
			}

		}

		// recursively simulate the sequences
		std::map<size_t, size_t> numLeaves = simulate( root, siteData, rateIndex );

		if((type & NO_CONSTANT_SITES) && numLeaves.size() == 1){
			i--;
			continue;
		}else if((type & NO_UNINFORMATIVE) && numLeaves.size() == 2){
			for(std::map<size_t, size_t>::iterator it = numLeaves.begin(); it != numLeaves.end(); it++)
				if(it->second == 1){
					i--;
					continue;
				}
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

template<class charType, class treeType>
std::map<size_t, size_t> RevBayesCore::PartitionedSiteCorrectionModel<charType, treeType>::simulate( const TopologyNode &node, std::vector<charType> &data, size_t rateIndex) {

	std::map<size_t, size_t> numLeaves;
    // get the children of the node
    const std::vector<TopologyNode*>& children = node.getChildren();

    // get the sequence of this node
    size_t nodeIndex = node.getIndex();
    charType &parentState = data[ nodeIndex ];

    // simulate the sequence for each child
    RandomNumberGenerator* rng = GLOBAL_RNG;
    for (std::vector< TopologyNode* >::const_iterator it = children.begin(); it != children.end(); ++it)
    {
        const TopologyNode &child = *(*it);

        // update the transition probability matrix
        this->updateTransitionProbabilities( child.getIndex(), child.getBranchLength() );

        unsigned long pstate = parentState.getState();
		size_t p = 0;
		while ( pstate != 1 )
		{
			// shift to the next state
			pstate >>= 1;
			// increase the index
			++p;
		}

		double *pij = this->transitionProbMatrices[ rateIndex ][ p ];

        charType &childState = data[ child.getIndex() ];

		// draw the state
		double u = rng->uniform01();
		while ( true )
		{
			u -= *pij;

			if ( u > 0.0 )
			{
				++childState;
				++pij;
			}
			else
			{
				break;
			}
		}

		if(child.isTip()){
		    numLeaves[childState.getState()]++;
		}else{
			std::map<size_t, size_t> l = simulate( child, data, rateIndex );
			for(std::map<size_t, size_t>::iterator it = l.begin(); it != l.end(); it++)
				numLeaves[it->first] += it->second;
		}
	}

    return numLeaves;

}


template<class charType, class treeType>
void RevBayesCore::PartitionedSiteCorrectionModel<charType, treeType>::resizeLikelihoodVectors( void ) {

    // we resize the partial likelihood vectors to the new dimensions
    size_t numSiteRates = this->numSiteRates;
    size_t numPatterns = this->numPatterns;
    size_t numChars = this->numChars;
    size_t numNodes = this->tau->getValue().getNumberOfNodes();

    size_t numCorrectionSites = 0;
    for(size_t partition = 0; partition < numCorrectionPartitions; partition++)
    	numCorrectionSites += numCorrectionSitesPerPartition[partition];

    this->partialLikelihoods.clear();
    this->partialLikelihoods.resize(2*numNodes*numSiteRates*(numPatterns+numCorrectionSites)*numChars);

    this->transitionProbMatrices = std::vector<TransitionProbabilityMatrix>(this->numSiteRates, TransitionProbabilityMatrix(this->numChars) );

    // set the offsets for easier iteration through the likelihood vector
    this->activeLikelihoodOffset      =  numNodes*numSiteRates*(numPatterns+numCorrectionSites)*numChars;
    this->nodeOffset                  =  numSiteRates*(numPatterns+numCorrectionSites)*numChars;
    this->mixtureOffset               =  (numPatterns+numCorrectionSites)*numChars;
}


#endif

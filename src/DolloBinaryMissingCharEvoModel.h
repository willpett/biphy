#ifndef DolloBinaryMissingCharEvoModel_H
#define DolloBinaryMissingCharEvoModel_H

#include "DolloBinaryCharEvoModel.h"
#include "TopologyNode.h"

namespace RevBayesCore {
    
    template<class treeType>
    class DolloBinaryMissingCharEvoModel : public DolloBinaryCharEvoModel<treeType> {
        
    public:
    	DolloBinaryMissingCharEvoModel(const TypedDagNode< treeType > *p, const TypedDagNode<Topology> *t, bool c, size_t nSites, int type = 0);
        DolloBinaryMissingCharEvoModel(const DolloBinaryMissingCharEvoModel &n);                                                                                                //!< Copy constructor
        virtual                                            ~DolloBinaryMissingCharEvoModel(void);                                                                   //!< Virtual destructor
        
        // public member functions
        DolloBinaryMissingCharEvoModel*     clone(void) const;
		void                                swapParameter(const DagNode *oldP, const DagNode *newP);
		double 								getMissingRate(void);
		void								setMissingRate(const TypedDagNode< double > *r);
        
    protected:

	   virtual void                        computeTipCorrection(const TopologyNode &node, size_t nodeIndex);
	   virtual void                        computeTipLikelihood(const TopologyNode &node, size_t nodeIndex);

	   const TypedDagNode< double >*		fracMissing;

    };
    
}

#include <cmath>
#include <cstring>

template<class treeType>
RevBayesCore::DolloBinaryMissingCharEvoModel<treeType>::DolloBinaryMissingCharEvoModel(const TypedDagNode< treeType > *p, const TypedDagNode<Topology> *t, bool c, size_t nSites, int type) :
	DolloBinaryCharEvoModel<treeType>(p, t, c , nSites, type),
	AbstractCharEvoModel<StandardState, treeType>(p, 2, c , nSites),
	fracMissing(NULL)
{
}


template<class treeType>
RevBayesCore::DolloBinaryMissingCharEvoModel<treeType>::DolloBinaryMissingCharEvoModel(const DolloBinaryMissingCharEvoModel &d) :
	DolloBinaryCharEvoModel<treeType>(d),
	AbstractCharEvoModel<StandardState, treeType>(d),
	fracMissing(d.fracMissing)
{
}


template<class treeType>
RevBayesCore::DolloBinaryMissingCharEvoModel<treeType>::~DolloBinaryMissingCharEvoModel( void ) {
    // We don't delete the parameters, because they might be used somewhere else too. The model needs to do that!
    
}


template<class treeType>
RevBayesCore::DolloBinaryMissingCharEvoModel<treeType>* RevBayesCore::DolloBinaryMissingCharEvoModel<treeType>::clone( void ) const {
    
    return new DolloBinaryMissingCharEvoModel<treeType>( *this );
}

template<class treeType>
void RevBayesCore::DolloBinaryMissingCharEvoModel<treeType>::computeTipCorrection(const TopologyNode &node, size_t nodeIndex)
{
	if(this->numCorrectionSites == 0)
		return;

	std::vector<double>::iterator p_node = this->partialLikelihoods.begin() + this->activeLikelihood[nodeIndex]*this->activeLikelihoodOffset + nodeIndex*this->nodeOffset;

	// iterate over all mixture categories

	double miss = getMissingRate();

	for (size_t mixture = 0; mixture < this->numSiteRates; ++mixture)
	{

		size_t offset = mixture*this->mixtureOffset + this->numPatterns*this->siteOffset;

		std::vector<double>::iterator     	u      = p_node + offset;

		//Probability of no observed presence in all leaves descending from this node
		// = probability of no observation
		u[0] = miss;

		//Probability of observed presence in only one leaf descending from this node
		// = probability of observation
		u[1] = 1.0 - miss;

		//Probability of observed presence in all leaves descending from this node
		// = probability of observation
		u[2] = 1.0 - miss;

		//Probability of observed presence in all but one leaves descending from this node
		// = probability of no observation
		u[3] = miss;

		//get the 1->1 transition probabilities for each branch
		this->updateTransitionProbabilities( nodeIndex, node.getBranchLength() );
		double pr    		= this->transitionProbMatrices[mixture][1][1];

		double r = 1.0;
		if(this->rateVariationAcrossSites == true)
			r = this->siteRates->getValue()[mixture];

		//Probability of number of observations M descending from this node
		double prob = 1.0;

		if(this->type & NO_SINGLETON_GAINS)
			prob -= 1.0 - miss; //  P( M != 1) = 1 - u[1] = 1 - u[2] = miss
		else if(this->type & NO_ABSENT_SITES)
			prob -= miss; // P( M != 0) = 1 - u[0] = 1 - u[3] = 1.0 - miss

		if(prob < 0.0)
			prob = 0.0;

		this->totalmass[nodeIndex][mixture][0] = prob;
		this->totalmass[nodeIndex][mixture][1] = (1-pr)/r;

	}
}

template<class treeType>
void RevBayesCore::DolloBinaryMissingCharEvoModel<treeType>::computeTipLikelihood(const TopologyNode &node, size_t nodeIndex)
{

	std::vector<double>::iterator p_node = this->partialLikelihoods.begin() + this->activeLikelihood[nodeIndex]*this->activeLikelihoodOffset + nodeIndex*this->nodeOffset;

	std::string name = node.getName();
	const std::vector<bool> &gap_node = this->gapMatrix[name];
	const std::vector<unsigned long> &char_node = this->charMatrix[name];

	double miss = getMissingRate();

	// compute the transition probabilities
	this->updateTransitionProbabilities( nodeIndex, node.getBranchLength() );

	// iterate over all mixture categories
	for (size_t mixture = 0; mixture < this->numSiteRates; ++mixture)
	{
		// the transition probability matrix for this mixture category
		const double *	tp_begin    		= this->transitionProbMatrices[mixture].theMatrix;

		// iterate over all sites
		for (size_t site = 0; site != this->numPatterns; ++site)
		{

			size_t offset = mixture*this->mixtureOffset + site*this->siteOffset;

			std::vector<double>::iterator     	p_site_mixture      = p_node + offset;

			// is this site a gap?
			if ( gap_node[site] )
			{
				// since this is a gap we need to assume that the actual state could have been any state
				// iterate over all initial states for the transitions
				for (size_t c1 = 0; c1 < this->numChars; ++c1)
				{

					// store the likelihood
					p_site_mixture[c1] = miss;

				}
			}
			else // we have observed a character
			{

				// get the original character
				unsigned long org_val = char_node[site];

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

					}

					p_site_mixture[c1] *= 1.0 - miss;

				}

			}

		}

	}

    computeTipCorrection(node, nodeIndex);

}

template<class treeType>
void RevBayesCore::DolloBinaryMissingCharEvoModel<treeType>::setMissingRate(const TypedDagNode< double > *r) {

    // remove the old parameter first
    if ( fracMissing != NULL )
    {
        this->removeParameter( fracMissing );
        fracMissing = NULL;
    }

    fracMissing = r;

    // add the parameter
    this->addParameter( fracMissing );

    // redraw the current value
    if ( this->dagNode != NULL && !this->dagNode->isClamped() )
    {
        this->redrawValue();
    }

}

template<class treeType>
double RevBayesCore::DolloBinaryMissingCharEvoModel<treeType>::getMissingRate() {

    double missingRate;
    if ( fracMissing != NULL )
    {
    	missingRate = fracMissing->getValue();
    }
    else
    {
    	missingRate = 0.0;
    }

    return missingRate;

}

template<class treeType>
void RevBayesCore::DolloBinaryMissingCharEvoModel<treeType>::swapParameter(const DagNode *oldP, const DagNode *newP) {
    
    if (oldP == fracMissing)
    {
    	fracMissing = static_cast<const TypedDagNode<double >* >( newP );
    }
    else
    {
    	DolloBinaryCharEvoModel<treeType>::swapParameter(oldP,newP);
    }
    
}

#endif

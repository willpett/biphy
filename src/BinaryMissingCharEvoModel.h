#ifndef BinaryMissingCharEvoModel_H
#define BinaryMissingCharEvoModel_H

#include "BinaryCharEvoModel.h"
#include "TopologyNode.h"

namespace RevBayesCore {

    template<class treeType>
    class BinaryMissingCharEvoModel : public BinaryCharEvoModel<treeType> {

    public:
    	BinaryMissingCharEvoModel(const TypedDagNode< treeType > *p, bool c, size_t nSites, int type = 0);
        BinaryMissingCharEvoModel(const BinaryMissingCharEvoModel &n);
        virtual                             ~BinaryMissingCharEvoModel(void);
        
        // public member functions
        BinaryMissingCharEvoModel*         	clone(void) const;
        void                                swapParameter(const DagNode *oldP, const DagNode *newP);
        double 								getMissingRate(void);
        void 								setMissingRate(const TypedDagNode< double > *r);

    protected:
        
        virtual void                        computeTipLikelihood(const TopologyNode &node, size_t nIdx);
		virtual void                        computeTipCorrection(const TopologyNode &node, size_t nIdx);

        const TypedDagNode< double >*		fracMissing;
    };
    
}

#include <cmath>
#include <cstring>

template<class treeType>
RevBayesCore::BinaryMissingCharEvoModel<treeType>::BinaryMissingCharEvoModel(const TypedDagNode< treeType > *t, bool c, size_t nSites, int type) :
	BinaryCharEvoModel<treeType>(t, c , nSites, type),
	AbstractCharEvoModel<StandardState, treeType>(t, 2, c , nSites),
	fracMissing(NULL)
{
	this->numCorrectionSites = 4*(type > 0);

	this->resizeLikelihoodVectors();
}


template<class treeType>
RevBayesCore::BinaryMissingCharEvoModel<treeType>::BinaryMissingCharEvoModel(const BinaryMissingCharEvoModel &d) :
	BinaryCharEvoModel<treeType>(d),
	AbstractCharEvoModel<StandardState, treeType>(d),
	fracMissing(d.fracMissing)
{
}


template<class treeType>
RevBayesCore::BinaryMissingCharEvoModel<treeType>::~BinaryMissingCharEvoModel( void ) {
    // We don't delete the parameters, because they might be used somewhere else too. The model needs to do that!
    
}


template<class treeType>
RevBayesCore::BinaryMissingCharEvoModel<treeType>* RevBayesCore::BinaryMissingCharEvoModel<treeType>::clone( void ) const {
    
    return new BinaryMissingCharEvoModel<treeType>( *this );
}

template<class treeType>
void RevBayesCore::BinaryMissingCharEvoModel<treeType>::computeTipCorrection(const TopologyNode &node, size_t nodeIndex)
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
		u[0] = 1.0;
		u[1] = miss;

		//Probability of observed presence in only one leaf descending from this node
		u[2] = 0.0;
		u[3] = 1.0 - miss;

		//Probability of observed presence in all leaves descending from this node
		u[4] = 0.0;
		u[5] = 1.0 - miss;

		//Probability of observed presence in all but one leaves descending from this node
		u[6] = 1.0;
		u[7] = miss;

	}
}

template<class treeType>
void RevBayesCore::BinaryMissingCharEvoModel<treeType>::computeTipLikelihood(const TopologyNode &node, size_t nodeIndex)
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
void RevBayesCore::BinaryMissingCharEvoModel<treeType>::setMissingRate(const TypedDagNode< double > *r) {

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
double RevBayesCore::BinaryMissingCharEvoModel<treeType>::getMissingRate() {

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
void RevBayesCore::BinaryMissingCharEvoModel<treeType>::swapParameter(const DagNode *oldP, const DagNode *newP) {

    if (oldP == fracMissing)
    {
    	fracMissing = static_cast<const TypedDagNode<double >* >( newP );
    }
    else
    {
    	BinaryCharEvoModel<treeType>::swapParameter(oldP,newP);
    }

}


#endif

#ifndef AbstractSiteCorrectionModel_H
#define AbstractSiteCorrectionModel_H

#include "AbstractCharacterData.h"
#include "DiscreteTaxonData.h"
#include "DnaState.h"
#include "GeneralCharEvoModel.h"
#include "RateMatrix.h"
#include "StandardState.h"
#include "TopologyNode.h"
#include "TransitionProbabilityMatrix.h"
#include "Tree.h"
#include "TreeChangeEventListener.h"
#include "TypedDistribution.h"

#include <memory.h>

namespace RevBayesCore {
    
    /**
     * @brief Declaration of the character state evolution along a tree class.
     *
     * This file contains the distribution class for a character state evolving along a tree.
     * This abstract base class can be derived for any character evolution model with homogeneous mixture sites. A
     * homogeneous mixture model over sites is a model where all sites are drawn from the same distribution and the
     * specific instance of the per site parameter is integrated over. The per site parameter could be a rate scaler (e.g. the + gamma models)
     * or different rate matrices or anything else.
     *
     * The pruning algorithm is implemented in this base class and calles some few pure virtual methods. 
     * The important functions you have to override are:
     * - getRootFrequencies()
     * - updateTransitionProbabilities()
     *
     * The data is stored for convenience in this class in a matrix (std::vector<std::vector< unsigned > >) and can
     * be compressed.
     *
     * The current implementation assumes that all mixture categories have the same a priori probability.
     *
     *
     * @copyright Copyright 2009-
     * @author The RevBayes Development Core Team (Sebastian Hoehna)
     * @since 2012-06-17, version 1.0
     */
    template<class charType, class treeType>
    class AbstractSiteCorrectionModel : public GeneralCharEvoModel<charType, treeType> {
        
    public:
        // Note, we need the size of the alignment in the constructor to correctly simulate an initial state
        AbstractSiteCorrectionModel(const TypedDagNode< treeType > *t, size_t nChars, bool c, size_t nSites);
        AbstractSiteCorrectionModel(const AbstractSiteCorrectionModel &n);                                                                                          //!< Copy constructor
        virtual                                                            ~AbstractSiteCorrectionModel(void);
        
    protected:
        // helper method for this and derived classes
        virtual void                                        				computeRootLikelihood(size_t root, size_t l, size_t r);
		virtual void                                        				computeInternalNodeLikelihood(const TopologyNode &n, size_t nIdx, size_t l, size_t r);
		virtual void                                        				computeTipLikelihood(const TopologyNode &node, size_t nIdx);
		virtual void                                        				resizeLikelihoodVectors(void);
        
		virtual void                                        				touchSpecialization(DagNode *toucher);

        virtual void                                                        computeRootCorrection(size_t root, size_t l, size_t r);
        virtual void                                                        computeInternalNodeCorrection(const TopologyNode &n, size_t nIdx, size_t l, size_t r);
        virtual void                                                        computeTipCorrection(const TopologyNode &node, size_t nIdx);
        
        virtual void                                                        setCorrectionPatterns() = 0;


        size_t																N;
        double																perSiteCorrection;

        std::map<std::string,std::vector<unsigned long> >                   correctionCharMatrix;
        std::map<std::string,std::vector<bool> >                   			correctionGapMatrix;

        size_t																numCorrectionSites;

        std::vector<double>													perMixtureCorrections;

    private:
        
        
    
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
RevBayesCore::AbstractSiteCorrectionModel<charType, treeType>::AbstractSiteCorrectionModel(const TypedDagNode<treeType> *t, size_t nChars, bool c, size_t nSites) :
	GeneralCharEvoModel<charType, treeType>(t, nChars, c, nSites),
	AbstractCharEvoModel<charType, treeType>(t, nChars, c, nSites),
	perSiteCorrection(0.0), numCorrectionSites(0), perMixtureCorrections(this->numSiteRates,0.0), N(nSites)
{
}


template<class charType, class treeType>
RevBayesCore::AbstractSiteCorrectionModel<charType, treeType>::AbstractSiteCorrectionModel(const AbstractSiteCorrectionModel &n) :
	GeneralCharEvoModel<charType, treeType>( n ),
	AbstractCharEvoModel<charType, treeType>( n ),
	perSiteCorrection(n.perSiteCorrection),
	numCorrectionSites(n.numCorrectionSites),
	correctionCharMatrix(n.correctionCharMatrix),
	correctionGapMatrix(n.correctionGapMatrix),
	perMixtureCorrections(n.perMixtureCorrections),
	N(n.N)
{
}


template<class charType, class treeType>
RevBayesCore::AbstractSiteCorrectionModel<charType, treeType>::~AbstractSiteCorrectionModel( void ) {
    // We don't delete the paramoms, because they might be used somewhere else too. The model needs to do that!
}

template<class charType, class treeType>
void RevBayesCore::AbstractSiteCorrectionModel<charType, treeType>::computeRootLikelihood(size_t root, size_t left, size_t right)
{

    GeneralCharEvoModel<charType, treeType>::computeRootLikelihood(root, left, right);

    computeRootCorrection(root, left, right);

    this->lnProb -= log(this->numSites) + this->numSites*perSiteCorrection;

}

template<class charType, class treeType>
void RevBayesCore::AbstractSiteCorrectionModel<charType, treeType>::computeInternalNodeLikelihood(const TopologyNode &node, size_t nodeIndex, size_t left, size_t right)
{

    GeneralCharEvoModel<charType, treeType>::computeInternalNodeLikelihood(node, nodeIndex, left, right);

    computeInternalNodeCorrection(node, nodeIndex, left, right);

}

template<class charType, class treeType>
void RevBayesCore::AbstractSiteCorrectionModel<charType, treeType>::computeTipLikelihood(const TopologyNode &node, size_t nIdx)
{

    GeneralCharEvoModel<charType, treeType>::computeTipLikelihood(node, nIdx);

    computeTipCorrection(node, nIdx);

}

template<class charType, class treeType>
void RevBayesCore::AbstractSiteCorrectionModel<charType, treeType>::computeRootCorrection( size_t root, size_t left, size_t right)
{
    // reset the likelihood
    perSiteCorrection = 0.0;

    if(numCorrectionSites == 0)
    	return;

    // get the root frequencies
    const std::vector<double> &f                    = this->getRootFrequencies();
    std::vector<double>::const_iterator f_end       = f.end();
    std::vector<double>::const_iterator f_begin     = f.begin();

    // get the pointers to the partial likelihoods of the left and right subtree
    std::vector<double>::const_iterator p_left   = this->partialLikelihoods.begin() + this->activeLikelihood[left]*this->activeLikelihoodOffset + left*this->nodeOffset;
    std::vector<double>::const_iterator p_right  = this->partialLikelihoods.begin() + this->activeLikelihood[right]*this->activeLikelihoodOffset + right*this->nodeOffset;

    // get pointers the likelihood for both subtrees
    std::vector<double>::const_iterator   p_mixture_left     = p_left;
    std::vector<double>::const_iterator   p_mixture_right    = p_right;
    // iterate over all mixture categories
    for (size_t mixture = 0; mixture < this->numSiteRates; ++mixture)
    {
    	perMixtureCorrections[mixture] = 0.0;
        // get pointers to the likelihood for this mixture category
    	std::vector<double>::const_iterator   p_site_mixture_left     = p_mixture_left+this->numPatterns*this->siteOffset;
    	std::vector<double>::const_iterator   p_site_mixture_right    = p_mixture_right+this->numPatterns*this->siteOffset;
        // iterate over all sites
        for (size_t site = 0; site < numCorrectionSites; ++site)
        {
            // temporary variable storing the likelihood
            double tmp = 0.0;
            // get the pointer to the stationary frequencies
            std::vector<double>::const_iterator f_j             = f_begin;
            // get the pointers to the likelihoods for this site and mixture category
            std::vector<double>::const_iterator p_site_left_j   = p_site_mixture_left;
            std::vector<double>::const_iterator p_site_right_j  = p_site_mixture_right;
            // iterate over all starting states
            for (; f_j != f_end; ++f_j)
            {
                // add the probability of starting from this state
                tmp += *p_site_left_j * *p_site_right_j * *f_j;

                // increment pointers
                ++p_site_left_j; ++p_site_right_j;
            }
            // add the likelihood for this mixture category
            perMixtureCorrections[mixture] += tmp;

            // increment the pointers to the next site
            p_site_mixture_left+=this->siteOffset; p_site_mixture_right+=this->siteOffset;

        } // end-for over all sites (=patterns)

        // increment the pointers to the next mixture category
        p_mixture_left+=this->mixtureOffset; p_mixture_right+=this->mixtureOffset;

    } // end-for over all mixtures (=rate categories)

    // sum the log-likelihoods for all sites together
    for (size_t mixture = 0; mixture < this->numSiteRates; ++mixture)
    {
        perSiteCorrection += perMixtureCorrections[mixture];
    }
    // normalize the log-probability
    perSiteCorrection /= this->numSiteRates;
    perSiteCorrection = log(1-perSiteCorrection);

}

template<class charType, class treeType>
void RevBayesCore::AbstractSiteCorrectionModel<charType,treeType>::computeInternalNodeCorrection(const TopologyNode &node, size_t nodeIndex, size_t left, size_t right)
{
	if(numCorrectionSites == 0)
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
		size_t offset 	= mixture*this->mixtureOffset+this->numPatterns*this->siteOffset;
		std::vector<double>::iterator          p_site_mixture          = p_node + offset;
		std::vector<double>::const_iterator    p_site_mixture_left     = p_left + offset;
		std::vector<double>::const_iterator    p_site_mixture_right    = p_right + offset;
		// compute the per site probabilities

		for (size_t site = 0; site < numCorrectionSites ; ++site)
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

	} // end-for over all mixtures (=rate-categories)

}

template<class charType, class treeType>
void RevBayesCore::AbstractSiteCorrectionModel<charType, treeType>::computeTipCorrection(const TopologyNode &node, size_t nodeIndex)
{
	if(numCorrectionSites == 0)
		return;

	std::vector<double>::iterator p_node = this->partialLikelihoods.begin() + this->activeLikelihood[nodeIndex]*this->activeLikelihoodOffset + nodeIndex*this->nodeOffset;

    std::string name = node.getName();
    const std::vector<bool> &gap_node = correctionGapMatrix[name];
    const std::vector<unsigned long> &char_node = correctionCharMatrix[name];

    // compute the transition probabilities
    this->updateTransitionProbabilities( nodeIndex, node.getBranchLength() );

    std::vector<double>::iterator   p_mixture      = p_node;

    // iterate over all mixture categories
    for (size_t mixture = 0; mixture < this->numSiteRates; ++mixture)
    {
        // the transition probability matrix for this mixture category
    	const double *	tp_begin    		= this->transitionProbMatrices[mixture].theMatrix;

        // iterate over all sites
        for (size_t site = 0; site != numCorrectionSites; ++site)
        {

        	size_t offset = mixture*this->mixtureOffset + (this->numPatterns+site)*this->siteOffset;

        	std::vector<double>::iterator     	p_site_mixture      = p_mixture + offset;
            // is this site a gap?
            if ( gap_node[site] )
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
                        //l+= org_val;

                    }

                } // end-for over all possible initial character for the branch

            }

        }

    }
}

template<class charType, class treeType>
void RevBayesCore::AbstractSiteCorrectionModel<charType, treeType>::touchSpecialization( DagNode* affecter ) {

	if(affecter == this->siteRates)
		perMixtureCorrections = std::vector<double>(this->numSiteRates,0.0);

	GeneralCharEvoModel<charType, treeType>::touchSpecialization(affecter);

}


template<class charType, class treeType>
void RevBayesCore::AbstractSiteCorrectionModel<charType, treeType>::resizeLikelihoodVectors( void ) {
    // we resize the partial likelihood vectors to the new dimensions
    size_t numSiteRates = this->numSiteRates;
    size_t numPatterns = this->numPatterns;
    size_t numChars = this->numChars;
    size_t numNodes = this->tau->getValue().getNumberOfNodes();

    this->partialLikelihoods.clear();
    this->partialLikelihoods.resize(2*numNodes*numSiteRates*(numPatterns+numCorrectionSites)*numChars);

    this->transitionProbMatrices = std::vector<TransitionProbabilityMatrix>(this->numSiteRates, TransitionProbabilityMatrix(this->numChars) );

    // set the offsets for easier iteration through the likelihood vector
    this->activeLikelihoodOffset      =  numNodes*numSiteRates*(numPatterns+numCorrectionSites)*numChars;
    this->nodeOffset                  =  numSiteRates*(numPatterns+numCorrectionSites)*numChars;
    this->mixtureOffset               =  (numPatterns+numCorrectionSites)*numChars;

    this->setCorrectionPatterns();
}


#endif

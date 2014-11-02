#ifndef GeneralCharEvoModel_H
#define GeneralCharEvoModel_H

#include "AbstractSiteHomogeneousMixtureModel.h"
#include "AbstractBranchHeterogeneousModel.h"
#include "DnaState.h"
#include "RateMatrix.h"
#include "RbVector.h"
#include "TopologyNode.h"
#include "TransitionProbabilityMatrix.h"
#include "TreeChangeEventListener.h"
#include "TypedDistribution.h"

namespace RevBayesCore {
    
    template<class charType, class treeType>
    class GeneralCharEvoModel : public AbstractSiteHomogeneousMixtureModel<charType, treeType>, public AbstractBranchHeterogeneousModel<charType, treeType> {
        
    public:
        GeneralCharEvoModel(const TypedDagNode< treeType > *t, size_t nChars, bool c, size_t nSites);
        GeneralCharEvoModel(const GeneralCharEvoModel &n);                                                                                                //!< Copy constructor
        virtual                                            ~GeneralCharEvoModel(void);                                                                   //!< Virtual destructor
        
        // public member functions
        GeneralCharEvoModel*             					clone(void) const;
        virtual void                                        swapParameter(const DagNode *oldP, const DagNode *newP);
        virtual void                                        redrawValue(void);
        const std::vector< DiscreteTaxonData<charType> >&   getMapping(void);
        
    protected:
        
        virtual void                                        computeRootLikelihood(size_t root, size_t l, size_t r);
        virtual void                                        computeInternalNodeLikelihood(const TopologyNode &n, size_t nIdx, size_t l, size_t r);
        virtual void                                        computeTipLikelihood(const TopologyNode &node, size_t nIdx);
        virtual void                                        touchSpecialization(DagNode *toucher);
        virtual void                                        updateTransitionProbabilities(size_t nodeIdx, double brlen);
        virtual void                                        resizeLikelihoodVectors(void);

        virtual void                                        redrawMapping(void);
        virtual void                                        simulateMapping(const TopologyNode &node, size_t nodeIndex, std::vector<size_t> perSiteRates);

        
    private:        
        
        void                                                simulate(const TopologyNode& node, std::vector< DiscreteTaxonData< charType > > &t, const std::vector<size_t> &perSiteRates);
        std::vector< DiscreteTaxonData<charType> > 			mapping;

        std::vector<std::vector<double> > 					per_node_Likelihoods;

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
RevBayesCore::GeneralCharEvoModel<charType, treeType>::GeneralCharEvoModel(const TypedDagNode<treeType> *t, size_t nChars, bool c, size_t nSites) :
	AbstractSiteHomogeneousMixtureModel<charType, treeType>(t, nChars, c, nSites),
	AbstractBranchHeterogeneousModel<charType, treeType>(t, nChars, c, nSites),
	AbstractCharEvoModel<charType, treeType>(t, nChars, c, nSites)
{
}


template<class charType, class treeType>
RevBayesCore::GeneralCharEvoModel<charType, treeType>::GeneralCharEvoModel(const GeneralCharEvoModel &d) :
	AbstractBranchHeterogeneousModel<charType, treeType>( d ),
	AbstractSiteHomogeneousMixtureModel<charType, treeType>( d ),
	AbstractCharEvoModel<charType, treeType>(d)
{
}


template<class charType, class treeType>
RevBayesCore::GeneralCharEvoModel<charType, treeType>::~GeneralCharEvoModel( void ) {
    // We don't delete the parameters, because they might be used somewhere else too. The model needs to do that!
    
}


template<class charType, class treeType>
RevBayesCore::GeneralCharEvoModel<charType, treeType>* RevBayesCore::GeneralCharEvoModel<charType, treeType>::clone( void ) const {
    
    return new GeneralCharEvoModel<charType, treeType>( *this );
}

template<class charType, class treeType>
void RevBayesCore::GeneralCharEvoModel<charType, treeType>::computeRootLikelihood( size_t root, size_t left, size_t right)
{
    // reset the likelihood
    this->lnProb = 0.0;

    // get the root frequencies
    const std::vector<double> &f                    = this->getRootFrequencies();
    std::vector<double>::const_iterator f_end       = f.end();
    std::vector<double>::const_iterator f_begin     = f.begin();
    
    // get the pointers to the partial likelihoods of the left and right subtree
    std::vector<double>::const_iterator p_left   = this->partialLikelihoods.begin() + this->activeLikelihood[left]*this->activeLikelihoodOffset + left*this->nodeOffset;
    std::vector<double>::const_iterator p_right  = this->partialLikelihoods.begin() + this->activeLikelihood[right]*this->activeLikelihoodOffset + right*this->nodeOffset;
    
    // create a vector for the per mixture likelihoods
    // we need this vector to sum over the different mixture likelihoods
    std::vector<double> per_mixture_Likelihoods = std::vector<double>(this->numPatterns,0.0);
    
    // iterate over all mixture categories
    for (size_t mixture = 0; mixture < this->numSiteRates; ++mixture) 
    {
        // iterate over all sites
        for (size_t site = 0; site < this->numPatterns; ++site)
        {

            // temporary variable storing the likelihood
            double tmp = 0.0;

            // get the pointer to the stationary frequencies
            std::vector<double>::const_iterator f_j             = f_begin;

            size_t offset = mixture*this->mixtureOffset + site*this->siteOffset;

            // get the pointers to the likelihoods for this site and mixture category
            std::vector<double>::const_iterator p_site_left_j   = p_left + offset;
            std::vector<double>::const_iterator p_site_right_j  = p_right + offset;
            // iterate over all starting states
            for (; f_j != f_end; ++f_j) 
            {
                // add the probability of starting from this state
                tmp += *p_site_left_j * *p_site_right_j * *f_j;
                
                // increment pointers
                ++p_site_left_j; ++p_site_right_j;
            }
            // add the likelihood for this mixture category
            per_mixture_Likelihoods[site] += tmp;
            
        }
    }
    
    // sum the log-likelihoods for all sites together
    std::vector< size_t >::const_iterator patterns = this->patternCounts.begin();
    for (size_t site = 0; site < this->numPatterns; ++site, ++patterns)
    {
        this->lnProb += log( per_mixture_Likelihoods[site] ) * *patterns;
    }
    // normalize the log-probability
    this->lnProb -= log( this->numSiteRates ) * this->numSites;
    
}


template<class charType, class treeType>
void RevBayesCore::GeneralCharEvoModel<charType, treeType>::computeInternalNodeLikelihood(const TopologyNode &node, size_t nodeIndex, size_t left, size_t right)
{   
    // compute the transition probability matrix
    updateTransitionProbabilities( nodeIndex, node.getBranchLength() );

    // get the pointers to the partial likelihoods for this node and the two descendant subtrees
    std::vector<double>::const_iterator   p_left  = this->partialLikelihoods.begin() + this->activeLikelihood[left]*this->activeLikelihoodOffset + left*this->nodeOffset;
    std::vector<double>::const_iterator   p_right = this->partialLikelihoods.begin() + this->activeLikelihood[right]*this->activeLikelihoodOffset + right*this->nodeOffset;
    std::vector<double>::iterator         p_node  = this->partialLikelihoods.begin() + this->activeLikelihood[nodeIndex]*this->activeLikelihoodOffset + nodeIndex*this->nodeOffset;

    // iterate over all mixture categories
    for (size_t mixture = 0; mixture < this->numSiteRates; ++mixture)
    {
        // the transition probability matrix for this mixture category
    	const double *		tp_begin    = this->transitionProbMatrices[mixture].theMatrix;

        for (size_t site = 0; site < this->numPatterns ; ++site)
        {
            
        	// get the pointers to the likelihood for this mixture category
			size_t offset = mixture*this->mixtureOffset + site*this->siteOffset;

			std::vector<double>::iterator          	p_site_mixture          = p_node + offset;
			std::vector<double>::const_iterator    	p_site_mixture_left     = p_left + offset;
			std::vector<double>::const_iterator    	p_site_mixture_right    = p_right + offset;

            // get the pointers for this mixture category and this site
        	const double *       tp_a    = tp_begin;
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
                
            }
            
        }
        
    }
}


template<class charType, class treeType>
void RevBayesCore::GeneralCharEvoModel<charType, treeType>::computeTipLikelihood(const TopologyNode &node, size_t nodeIndex)
{
	std::vector<double>::iterator p_node = this->partialLikelihoods.begin() + this->activeLikelihood[nodeIndex]*this->activeLikelihoodOffset + nodeIndex*this->nodeOffset;
    
    std::string name = node.getName();
    const std::vector<bool> &gap_node = this->gapMatrix[name];
    const std::vector<unsigned long> &char_node = this->charMatrix[name];


    // compute the transition probabilities
    updateTransitionProbabilities( nodeIndex, node.getBranchLength() );

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
                        
                    }
                    
                }
                
            }
            
        }
        
    }
}


template<class charType, class treeType>
void RevBayesCore::GeneralCharEvoModel<charType, treeType>::resizeLikelihoodVectors( void ) {

    // we resize the partial likelihood vectors to the new dimensions
    size_t numSiteRates = this->numSiteRates;
    size_t numPatterns = this->numPatterns;
    size_t numChars = this->numChars;
    size_t numNodes = this->tau->getValue().getNumberOfNodes();

    this->partialLikelihoods.clear();
    this->partialLikelihoods.resize(2*numNodes*numSiteRates*numPatterns*numChars);

    this->transitionProbMatrices = std::vector<TransitionProbabilityMatrix>(this->numSiteRates, TransitionProbabilityMatrix(this->numChars) );

    // set the offsets for easier iteration through the likelihood vector
    this->activeLikelihoodOffset      =  numNodes*numSiteRates*numPatterns*numChars;
    this->nodeOffset                  =  numSiteRates*numPatterns*numChars;
    this->mixtureOffset               =  numPatterns*numChars;
}

template<class charType, class treeType>
void RevBayesCore::GeneralCharEvoModel<charType, treeType>::updateTransitionProbabilities(size_t nodeIdx, double brlen) {

    // first, get the rate matrix for this branch
    const RateMatrix *rm = this->getRateMatrix(nodeIdx);
    double branchTime = this->getClockRate(nodeIdx);
    const std::vector<double> &r = this->getSiteRates();

	for (size_t i = 0; i < r.size(); ++i)
	{
		rm->calculateTransitionProbabilities( branchTime * r[i] * brlen, this->transitionProbMatrices[i] );
	}

}


template<class charType, class treeType>
void RevBayesCore::GeneralCharEvoModel<charType, treeType>::swapParameter(const DagNode *oldP, const DagNode *newP) {
    
	AbstractBranchHeterogeneousModel<charType, treeType>::swapParameter(oldP,newP);
	AbstractSiteHomogeneousMixtureModel<charType, treeType>::swapParameter(oldP,newP);
    
}

template<class charType, class treeType>
void RevBayesCore::GeneralCharEvoModel<charType, treeType>::touchSpecialization( DagNode* affecter ) {
    
    AbstractBranchHeterogeneousModel<charType, treeType>::touchSpecialization(affecter);
    
}

template<class charType, class treeType>
void RevBayesCore::GeneralCharEvoModel<charType, treeType>::redrawValue( void ) {

    // delete the old value first
    delete this->value;

    // create a new character data object
    this->value = new DiscreteCharacterData<charType>();

    // create a vector of taxon data
    std::vector< DiscreteTaxonData<charType> > taxa = std::vector< DiscreteTaxonData< charType > >( this->tau->getValue().getNumberOfNodes(), DiscreteTaxonData<charType>() );

    // first, simulate the per site rates
    RandomNumberGenerator* rng = GLOBAL_RNG;
    std::vector<size_t> perSiteRates;
    for ( size_t i = 0; i < this->numSites; ++i )
    {
        // draw the state
        double u = rng->uniform01();
        size_t rateIndex = (int)(u*this->numSiteRates);
        perSiteRates.push_back( rateIndex );
    }

    // simulate the root sequence
    const std::vector< double > &stationaryFreqs = this->getRootFrequencies();
    DiscreteTaxonData< charType > &root = taxa[ this->tau->getValue().getRoot().getIndex() ];
    for ( size_t i = 0; i < this->numSites; ++i )
    {
        // create the character
        charType c;
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

        // add the character to the sequence
        root.addCharacter( c );
    }

    // recursively simulate the sequences
    simulate( this->tau->getValue().getRoot(), taxa, perSiteRates );

    // add the taxon data to the character data
    for (size_t i = 0; i < this->tau->getValue().getNumberOfTips(); ++i)
    {
        this->value->addTaxonData( taxa[i] );
    }

    // compress the data and initialize internal variables
    this->compress();

}

template<class charType, class treeType>
void RevBayesCore::GeneralCharEvoModel<charType, treeType>::simulate( const TopologyNode &node, std::vector< DiscreteTaxonData< charType > > &taxa, const std::vector<size_t> &perSiteRates) {

    // get the children of the node
    const std::vector<TopologyNode*>& children = node.getChildren();

    // get the sequence of this node
    size_t nodeIndex = node.getIndex();
    const DiscreteTaxonData< charType > &parent = taxa[ nodeIndex ];

    // simulate the sequence for each child
    RandomNumberGenerator* rng = GLOBAL_RNG;
    for (std::vector< TopologyNode* >::const_iterator it = children.begin(); it != children.end(); ++it)
    {
        const TopologyNode &child = *(*it);

        // update the transition probability matrix
        updateTransitionProbabilities( child.getIndex(), child.getBranchLength() );

        DiscreteTaxonData< charType > &taxon = taxa[ child.getIndex() ];
        for ( size_t i = 0; i < this->numSites; ++i )
        {
            // get the ancestral character for this site
            unsigned long parentState = parent.getCharacter( i ).getState();
            size_t p = 0;
            while ( parentState != 1 )
            {
                // shift to the next state
                parentState >>= 1;
                // increase the index
                ++p;
            }

            double *freqs = this->transitionProbMatrices[ perSiteRates[i] ][ p ];

            // create the character
            charType c;
            c.setToFirstState();
            // draw the state
            double u = rng->uniform01();
            while ( true )
            {
                u -= *freqs;

                if ( u > 0.0 )
                {
                    ++c;
                    ++freqs;
                }
                else
                {
                    break;
                }
            }

            // add the character to the sequence
            taxon.addCharacter( c );
        }

        if ( child.isTip() )
        {
            taxon.setTaxonName( child.getName() );
        }
        else
        {
            // recursively simulate the sequences
            simulate( child, taxa, perSiteRates );
        }

    }

}

template<class charType, class treeType>
const std::vector< RevBayesCore::DiscreteTaxonData<charType> >& RevBayesCore::GeneralCharEvoModel<charType, treeType>::getMapping(void) {

	redrawMapping();
	return mapping;
}

template<class charType, class treeType>
void RevBayesCore::GeneralCharEvoModel<charType, treeType>::redrawMapping(void)
{
    this->computeLnProbability();

	mapping = std::vector< DiscreteTaxonData< charType > >( this->tau->getValue().getNumberOfNodes(), DiscreteTaxonData<charType>() );

	// first, simulate the per site rates
	RandomNumberGenerator* rng = GLOBAL_RNG;

	std::vector<size_t> perSiteRates;
	for ( size_t i = 0; i < this->numSites; ++i )
	{
		// draw the state
		double u = rng->uniform01();
		size_t rateIndex = (int)(u*this->numSiteRates);
		perSiteRates.push_back( rateIndex );
	}

	// compute the ln probability by recursively calling the probability calculation for each node
	const TopologyNode &root = this->tau->getValue().getRoot();

	// we start with the root and then traverse down the tree
	size_t rootIndex = root.getIndex();

	// start by filling the likelihood vector for the two children of the root
	const TopologyNode &left = root.getChild(0);
	size_t leftIndex = left.getIndex();
	simulateMapping( left, leftIndex, perSiteRates );
	const TopologyNode &right = root.getChild(1);
	size_t rightIndex = right.getIndex();
	simulateMapping( right, rightIndex , perSiteRates);

	const std::vector<double> &f = this->getRootFrequencies();

	DiscreteTaxonData< charType > &rootData = mapping[ rootIndex ];
	for ( size_t i = 0; i < this->numSites; ++i )
	{
		size_t offset = perSiteRates[i]*this->mixtureOffset + this->patternMap[i]*this->siteOffset;

		std::vector<double>::const_iterator   p_left  = this->partialLikelihoods.begin() + this->activeLikelihood[leftIndex]*this->activeLikelihoodOffset + leftIndex*this->nodeOffset + offset;
		std::vector<double>::const_iterator   p_right = this->partialLikelihoods.begin() + this->activeLikelihood[rightIndex]*this->activeLikelihoodOffset + rightIndex*this->nodeOffset + offset;

		double sum = 0.0;
		for (size_t c1 = 0; c1 < this->numChars; ++c1)
		{
			sum += p_left[c1]*p_right[c1]*f[c1];
		}


		// create the character
		charType c;
		c.setToFirstState();
		// draw the state
		double u = rng->uniform01();

		double tmp = 0.0;
		while(tmp < u*sum)
		{
			unsigned long c1 = c.getState()-1;
			tmp += p_left[c1]*p_right[c1]*f[c1];
			c++;
		}
		c--;

		rootData.addCharacter(c);
	}

}

template<class charType, class treeType>
void RevBayesCore::GeneralCharEvoModel<charType, treeType>::simulateMapping(const TopologyNode &node, size_t nodeIndex, std::vector<size_t> perSiteRates ) {

    // first, simulate the per site rates
    RandomNumberGenerator* rng = GLOBAL_RNG;

    if(node.isTip()){
    	for ( size_t i = 0; i < this->numSites; ++i )
    		mapping[nodeIndex].addCharacter(this->value->getTaxonData(nodeIndex).getCharacter(i));
    	return;
    }

    const TopologyNode &left = node.getChild(0);
	size_t leftIndex = left.getIndex();
	simulateMapping( left, leftIndex, perSiteRates );
	const TopologyNode &right = node.getChild(1);
	size_t rightIndex = right.getIndex();
	simulateMapping( right, rightIndex , perSiteRates);


    // simulate the sequence at this node
    DiscreteTaxonData< charType > &nodeData = mapping[ nodeIndex ];
    for ( size_t i = 0; i < this->numSites; ++i )
    {
        size_t offset = perSiteRates[i]*this->mixtureOffset + this->patternMap[i]*this->siteOffset;

        std::vector<double>::const_iterator   p_left  = this->partialLikelihoods.begin() + this->activeLikelihood[leftIndex]*this->activeLikelihoodOffset + leftIndex*this->nodeOffset + offset;
	    std::vector<double>::const_iterator   p_right = this->partialLikelihoods.begin() + this->activeLikelihood[rightIndex]*this->activeLikelihoodOffset + rightIndex*this->nodeOffset + offset;

		double sum = 0.0;
		for (size_t c1 = 0; c1 < this->numChars; ++c1)
		{
			sum += p_left[c1]*p_right[c1];
		}

		// create the character
		charType c;
		c.setToFirstState();
		// draw the state
		double u = rng->uniform01();

		double tmp = 0.0;
		while(tmp < u*sum)
		{
			tmp += p_left[c.getState()-1]*p_right[c.getState()-1];
			c++;
		}
		c--;

		nodeData.addCharacter(c);
    }

}



#endif

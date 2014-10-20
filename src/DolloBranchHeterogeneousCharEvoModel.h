#ifndef DolloBranchHeterogeneousCharEvoModel_H
#define DolloBranchHeterogeneousCharEvoModel_H

#include "AbstractSiteHomogeneousMixtureCharEvoModel.h"
#include "DnaState.h"
#include "RateMatrix.h"
#include "RbVector.h"
#include "TopologyNode.h"
#include "TransitionProbabilityMatrix.h"
#include "TreeChangeEventListener.h"
#include "TypedDistribution.h"

namespace RevBayesCore {
    
    template<class charType, class treeType>
    class DolloBranchHeterogeneousCharEvoModel : public AbstractSiteHomogeneousMixtureCharEvoModel<charType, treeType> {
        
    public:
        DolloBranchHeterogeneousCharEvoModel(const TypedDagNode< treeType > *p, const TypedDagNode<Topology> *t, size_t nChars, bool c, size_t nSites, charType& astate);
        DolloBranchHeterogeneousCharEvoModel(const DolloBranchHeterogeneousCharEvoModel &n);                                                                                                //!< Copy constructor
        virtual                                            ~DolloBranchHeterogeneousCharEvoModel(void);                                                                   //!< Virtual destructor
        
        // public member functions
        DolloBranchHeterogeneousCharEvoModel*         		clone(void) const;                                                                          //!< Create an independent clone
        void                                                setClockRate(const TypedDagNode< double > *r);
		void                                                setClockRate(const TypedDagNode< std::vector< double > > *r);
		void                                                setRateMatrix(const TypedDagNode< RateMatrix > *rm);
		void                                                setRateMatrix(const TypedDagNode< RbVector< RateMatrix > > *rm);
		void                                                setRootFrequencies(const TypedDagNode< std::vector< double > > *f);
		void                                                setSiteRates(const TypedDagNode< std::vector< double > > *r);
		void                                                swapParameter(const DagNode *oldP, const DagNode *newP);

		//virtual void                                        redrawValue(void);
        
    protected:

	   void                                                computeRootLikelihood(size_t root, size_t l, size_t r);
	   void                                                computeInternalNodeLikelihood(const TopologyNode &n, size_t nIdx, size_t l, size_t r);
	   void                                                computeTipLikelihood(const TopologyNode &node, size_t nIdx);
	   const std::vector<double>&                          getRootFrequencies(void);
// (not needed)        void                                                keepSpecialization(DagNode* affecter);
// (not needed)        void                                                restoreSpecialization(DagNode *restorer);
	   void                                                touchSpecialization(DagNode *toucher);
	   void                                                updateTransitionProbabilities(size_t nodeIdx, double brlen);
        
        
    private:
        void                                                computeAncestral(void);
        bool                                                computeAncestralMap(const TopologyNode& node, size_t site);
        void                                                computeAncestralNodes(const TopologyNode& node, size_t site);

        // members
		const TypedDagNode< double >*                       homogeneousClockRate;
		const TypedDagNode< std::vector< double > >*        heterogeneousClockRates;
		const TypedDagNode< RateMatrix >*                   homogeneousRateMatrix;
		const TypedDagNode< RbVector< RateMatrix > >*       heterogeneousRateMatrices;
		const TypedDagNode< std::vector< double > >*        rootFrequencies;
		const TypedDagNode< std::vector< double > >*        siteRates;

		double												omega;

		std::vector<std::vector<std::vector<double> > >		ui;
		std::vector<double>									totalmass;


		// flags specifying which model variants we use
		bool                                                branchHeterogeneousClockRates;
		bool                                                branchHeterogeneousSubstitutionMatrices;
		bool                                                rateVariationAcrossSites;

        std::vector<std::vector<size_t> >                   ancestralNodes;
        std::map<size_t,bool>								ancestralMap;
        charType											absentstate;

        const TypedDagNode<Topology>*						topology;
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
RevBayesCore::DolloBranchHeterogeneousCharEvoModel<charType, treeType>::DolloBranchHeterogeneousCharEvoModel(const TypedDagNode< treeType > *p, const TypedDagNode<Topology> *t, size_t nChars, bool c, size_t nSites, charType& astate) : AbstractSiteHomogeneousMixtureCharEvoModel<charType, treeType>(  p, nChars, 1, c, nSites ) ,
	absentstate(astate), topology(t), omega(0.0)
{

	// initialize with default parameters
	homogeneousClockRate        = new ConstantNode<double>("clockRate", new double(1.0) );
	heterogeneousClockRates     = NULL;
	homogeneousRateMatrix       = new ConstantNode<RateMatrix>("rateMatrix", new RateMatrix_JC( nChars ) );
	heterogeneousRateMatrices   = NULL;
	rootFrequencies             = NULL;
	siteRates                   = NULL;


	// flags specifying which model variants we use
	branchHeterogeneousClockRates               = false;
	branchHeterogeneousSubstitutionMatrices     = false;
	rateVariationAcrossSites                    = false;

	// add the parameters to the parents list
	this->addParameter( homogeneousClockRate );
	this->addParameter( homogeneousRateMatrix );
	this->addParameter( topology );

	totalmass.resize(this->tau->getValue().getNumberOfNodes());
	ui.resize(this->tau->getValue().getNumberOfNodes());
	for(size_t i = 0; i < ui.size(); i ++){
		ui[i].resize(2);
		for(size_t m = 0; m < 2; m++){
			ui[i][m].resize(this->numSiteRates);
		}
	}

    //this->redrawValue();

}


template<class charType, class treeType>
RevBayesCore::DolloBranchHeterogeneousCharEvoModel<charType, treeType>::DolloBranchHeterogeneousCharEvoModel(const DolloBranchHeterogeneousCharEvoModel &d) : AbstractSiteHomogeneousMixtureCharEvoModel<charType, treeType>( d ),
	absentstate(d.absentstate), ancestralNodes(d.ancestralNodes), topology(d.topology)
{
	// parameters are automatically copied
	// initialize with default parameters
	homogeneousClockRate        = d.homogeneousClockRate;
	heterogeneousClockRates     = d.heterogeneousClockRates;
	homogeneousRateMatrix       = d.homogeneousRateMatrix;
	heterogeneousRateMatrices   = d.heterogeneousRateMatrices;
	rootFrequencies             = d.rootFrequencies;
	siteRates                   = d.siteRates;


	// flags specifying which model variants we use
	branchHeterogeneousClockRates               = d.branchHeterogeneousClockRates;
	branchHeterogeneousSubstitutionMatrices     = d.branchHeterogeneousSubstitutionMatrices;
	rateVariationAcrossSites                    = d.rateVariationAcrossSites;

	totalmass.clear();
	ui.clear();

	totalmass.resize(this->tau->getValue().getNumberOfNodes());
	ui.resize(this->tau->getValue().getNumberOfNodes());
	for(size_t i = 0; i < ui.size(); i ++){
		ui[i].clear();
		ui[i].resize(2);
		for(size_t m = 0; m < 2; m++){
			ui[i][m].clear();
			ui[i][m].resize(this->numSiteRates);
		}
	}
}


template<class charType, class treeType>
RevBayesCore::DolloBranchHeterogeneousCharEvoModel<charType, treeType>::~DolloBranchHeterogeneousCharEvoModel( void ) {
    // We don't delete the parameters, because they might be used somewhere else too. The model needs to do that!
    
}


template<class charType, class treeType>
RevBayesCore::DolloBranchHeterogeneousCharEvoModel<charType, treeType>* RevBayesCore::DolloBranchHeterogeneousCharEvoModel<charType, treeType>::clone( void ) const {
    
    return new DolloBranchHeterogeneousCharEvoModel<charType, treeType>( *this );
}



template<class charType, class treeType>
void RevBayesCore::DolloBranchHeterogeneousCharEvoModel<charType, treeType>::computeRootLikelihood( size_t root, size_t left, size_t right)
{
    // reset the likelihood
    this->lnProb = 0.0;

    // get the root frequencies
    const std::vector<double> &f                    = getRootFrequencies();
    std::vector<double>::const_iterator f_end       = f.end();
    std::vector<double>::const_iterator f_begin     = f.begin();

    // create a vector for the per mixture likelihoods
    // we need this vector to sum over the different mixture likelihoods
    std::vector<double> per_mixture_Likelihoods = std::vector<double>(this->numPatterns,0.0);

	for (size_t site = 0; site < this->numPatterns; site++){
		// iterate over all nodes in the ancestral set for this site pattern
		for (size_t anc = 0; anc < ancestralNodes[site].size(); anc++){
			size_t idx = ancestralNodes[site][anc];
			const TopologyNode &node = this->tau->getValue().getNode(idx);

			size_t nleft = node.getChild(0).getIndex();
			size_t nright = node.getChild(1).getIndex();

			// get the pointers to the partial likelihoods of the left and right subtree of this ancestral node
			const double* p_site_left   = this->partialLikelihoods + this->activeLikelihood[nleft]*this->activeLikelihoodOffset + nleft*this->nodeOffset + site*this->siteOffset;
			const double* p_site_right  = this->partialLikelihoods + this->activeLikelihood[nright]*this->activeLikelihoodOffset + nright*this->nodeOffset + site*this->siteOffset;

			// iterate over all mixture categories
			for (size_t mixture = 0; mixture < this->numSiteRates; ++mixture)
			{
					const double*   p_site_mixture_left = p_site_left;
				    const double*   p_site_mixture_right = p_site_right;

				    double prob_birth = 1.0;
					if(this->rateVariationAcrossSites == true){
						double r = this->siteRates->getValue()[mixture];
						prob_birth /= r;
						if(!node.isRoot())
							prob_birth *= (1.0 - exp(-r*node.getBranchLength()));
					}else{
						if(!node.isRoot())
							prob_birth = 1.0 - exp(-node.getBranchLength());
					}

					//if(prob_birth > 1.0)
					//	std::cerr << prob_birth << std::endl;

				    // temporary variable storing the likelihood
					double tmp = 0.0;
					// get the pointer to the stationary frequencies
					std::vector<double>::const_iterator f_j             = f_begin;
					// get the pointers to the likelihoods for this site and mixture category
					const double* p_site_left_j   = p_site_mixture_left;
					const double* p_site_right_j  = p_site_mixture_right;
					// iterate over all starting states
					for (; f_j != f_end; ++f_j)
					{
						// add the probability of starting from this state
						tmp += *p_site_left_j * *p_site_right_j * *f_j;
						// increment pointers
						++p_site_left_j; ++p_site_right_j;
					}

					// add the likelihood for this mixture category
					per_mixture_Likelihoods[site] += tmp*prob_birth;

					// increment the pointers to the next mixture category
					p_site_mixture_left+=this->mixtureOffset; p_site_mixture_right+=this->mixtureOffset;
			} //end-for over all mixtures (=rate categories)
		} // end-for over ancestral node set
    } // end-for over all sites (=patterns)

    // sum the log-likelihoods for all sites together
    std::vector< size_t >::const_iterator patterns = this->patternCounts.begin();
    for (size_t site = 0; site < this->numPatterns; ++site, ++patterns)
    {
        this->lnProb += log( per_mixture_Likelihoods[site] ) * *patterns;
    }

    //Nicholls and Gray 2003
    omega = 0;
    for (size_t node = 0; node < totalmass.size(); node++){
    	omega += totalmass[node];
    }

    double tl = this->tau->getValue().getNode(left).getBranchLength();
    double tr = this->tau->getValue().getNode(right).getBranchLength();

    totalmass[root] = 0;
    for (size_t mixture = 0; mixture < this->numSiteRates; ++mixture){
		double r = 1.0;
		if(this->rateVariationAcrossSites == true)
			r = this->siteRates->getValue()[mixture];

		ui[root][0][mixture] = (1 - exp(-r*tl)*(1 - ui[left][0][mixture] ) )*(1 - exp(-r*tr)*(1 - ui[right][0][mixture]) );
		ui[root][1][mixture] = (1 - exp(-r*tl)*(1 - ui[left][0][mixture] ) )*exp(-r*tr)*ui[right][1][mixture]
							  +(1 - exp(-r*tr)*(1 - ui[right][0][mixture]) )*exp(-r*tl)*ui[left][1][mixture];

		totalmass[root] += (1 - ui[root][0][mixture] - ui[root][1][mixture])/r;
	}

	omega += totalmass[root];

	//Alekseyenko et al. 2008
	//Integrated over lambda
    this->lnProb -= log(this->numSites) + this->numSites*log(omega);
}


template<class charType, class treeType>
void RevBayesCore::DolloBranchHeterogeneousCharEvoModel<charType, treeType>::computeInternalNodeLikelihood(const TopologyNode &node, size_t nodeIndex, size_t left, size_t right)
{
    
    // compute the transition probability matrix
    updateTransitionProbabilities( nodeIndex, node.getBranchLength() );

    // get the pointers to the partial likelihoods for this node and the two descendant subtrees
    const double*   p_left  = this->partialLikelihoods + this->activeLikelihood[left]*this->activeLikelihoodOffset + left*this->nodeOffset;
    const double*   p_right = this->partialLikelihoods + this->activeLikelihood[right]*this->activeLikelihoodOffset + right*this->nodeOffset;
    double*         p_node  = this->partialLikelihoods + this->activeLikelihood[nodeIndex]*this->activeLikelihoodOffset + nodeIndex*this->nodeOffset;

    // iterate over all mixture categories
    for (size_t mixture = 0; mixture < this->numSiteRates; ++mixture)
    {
        // the transition probability matrix for this mixture category
        const double*                       tp_begin    = this->transitionProbMatrices[mixture].theMatrix;

        // get the pointers to the likelihood for this mixture category
        size_t offset = mixture*this->mixtureOffset;
        double*          p_site_mixture          = p_node + offset;
        const double*    p_site_mixture_left     = p_left + offset;
        const double*    p_site_mixture_right    = p_right + offset;
        // compute the per site probabilities
        for (size_t site = 0; site < this->numPatterns ; ++site)
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

    // total mass computation for internal node
    // Nicholls and Gray 2003
    double tl = this->tau->getValue().getNode(left).getBranchLength();
    double tr = this->tau->getValue().getNode(right).getBranchLength();
    double t0 = this->tau->getValue().getNode(nodeIndex).getBranchLength();

    totalmass[nodeIndex] = 0;
    for (size_t mixture = 0; mixture < this->numSiteRates; ++mixture){
		double r = 1.0;
		if(this->rateVariationAcrossSites == true)
			r = this->siteRates->getValue()[mixture];

		ui[nodeIndex][0][mixture] = (1 - exp(-r*tl)*(1 - ui[left][0][mixture] ) )*(1 - exp(-r*tr)*(1 - ui[right][0][mixture]) );
		ui[nodeIndex][1][mixture] = (1 - exp(-r*tl)*(1 - ui[left][0][mixture] ) )*exp(-r*tr)*ui[right][1][mixture]
								   +(1 - exp(-r*tr)*(1 - ui[right][0][mixture]) )*exp(-r*tl)*ui[left][1][mixture];

		totalmass[nodeIndex] += (1 - ui[nodeIndex][0][mixture] - ui[nodeIndex][1][mixture])*(1-exp(-r*t0))/r;
    }

}


template<class charType, class treeType>
void RevBayesCore::DolloBranchHeterogeneousCharEvoModel<charType, treeType>::computeTipLikelihood(const TopologyNode &node, size_t nodeIndex)
{

    double* p_node = this->partialLikelihoods + this->activeLikelihood[nodeIndex]*this->activeLikelihoodOffset + nodeIndex*this->nodeOffset;

    std::string name = node.getName();
    const std::vector<bool> &gap_node = this->gapMatrix[name];
    const std::vector<unsigned long> &char_node = this->charMatrix[name];

    // compute the transition probabilities
    updateTransitionProbabilities( nodeIndex, node.getBranchLength() );
    
    double*   p_mixture      = p_node;
    
    // iterate over all mixture categories
    for (size_t mixture = 0; mixture < this->numSiteRates; ++mixture)
    {
        // the transition probability matrix for this mixture category
        const double*                       tp_begin    = this->transitionProbMatrices[mixture].theMatrix;

        // get the pointer to the likelihoods for this site and mixture category
        double*     p_site_mixture      = p_mixture;

        // iterate over all sites
        for (size_t site = 0; site != this->numPatterns; ++site)
        {
            
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
                        const double* d  = tp_begin+(this->numChars*c1);

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

                } // end-for over all possible initial character for the branch

            } // end-if a gap state

            // increment the pointers to next site
            p_site_mixture+=this->siteOffset;

        } // end-for over all sites/patterns in the sequence

        // increment the pointers to next mixture category
        p_mixture+=this->mixtureOffset;

    } // end-for over all mixture categories
    
    for(size_t mixture = 0; mixture < this->numSiteRates; mixture++){
		ui[nodeIndex][0][mixture] = 0;
		ui[nodeIndex][1][mixture] = 1;
    }
    totalmass[nodeIndex] = 0;
}


template<class charType, class treeType>
const std::vector<double>& RevBayesCore::DolloBranchHeterogeneousCharEvoModel<charType, treeType>::getRootFrequencies( void ) {

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
void RevBayesCore::DolloBranchHeterogeneousCharEvoModel<charType, treeType>::updateTransitionProbabilities(size_t nodeIdx, double brlen) {

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

    // second, get the clock rate for the branch
    double branchTime;
    if ( this->branchHeterogeneousClockRates == true )
    {
        branchTime = this->heterogeneousClockRates->getValue()[nodeIdx] * brlen;
    }
    else
    {
        branchTime = this->homogeneousClockRate->getValue() * brlen;
    }
    
    // and finally compute the per site rate transition probability matrix
    if ( this->rateVariationAcrossSites == true )
    {
        const std::vector<double> &r = this->siteRates->getValue();
        for (size_t i = 0; i < this->numSiteRates; ++i)
        {
            rm->calculateTransitionProbabilities( branchTime * r[i], this->transitionProbMatrices[i] );
        }
    }
    else
    {
        rm->calculateTransitionProbabilities( branchTime, this->transitionProbMatrices[0] );
    }

}



template<class charType, class treeType>
void RevBayesCore::DolloBranchHeterogeneousCharEvoModel<charType, treeType>::setClockRate(const TypedDagNode< double > *r) {
    
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
void RevBayesCore::DolloBranchHeterogeneousCharEvoModel<charType, treeType>::setClockRate(const TypedDagNode< std::vector< double > > *r) {
    
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
void RevBayesCore::DolloBranchHeterogeneousCharEvoModel<charType, treeType>::setRateMatrix(const TypedDagNode< RateMatrix > *rm) {

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
void RevBayesCore::DolloBranchHeterogeneousCharEvoModel<charType, treeType>::setRateMatrix(const TypedDagNode< RbVector< RateMatrix > > *rm) {

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
void RevBayesCore::DolloBranchHeterogeneousCharEvoModel<charType, treeType>::setRootFrequencies(const TypedDagNode< std::vector< double > > *f) {

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
void RevBayesCore::DolloBranchHeterogeneousCharEvoModel<charType, treeType>::setSiteRates(const TypedDagNode< std::vector< double > > *r) {
    
    // remove the old parameter first
    if ( siteRates != NULL )
    {
        this->removeParameter( siteRates );
        siteRates = NULL;
    }
    
    if ( r != NULL )
    {
        // set the value
        rateVariationAcrossSites = true;
        siteRates = r;
        this->numSiteRates = r->getValue().size();
        this->resizeLikelihoodVectors();
    }
    else
    {
        // set the value
        rateVariationAcrossSites = false;
        siteRates = NULL;
        this->numSiteRates = 1;
        this->resizeLikelihoodVectors();
        
    }
    
    // add the parameter
    this->addParameter( siteRates );

    // redraw the current value
    if ( this->dagNode != NULL && !this->dagNode->isClamped() )
    {
        this->redrawValue();
    }
}




template<class charType, class treeType>
void RevBayesCore::DolloBranchHeterogeneousCharEvoModel<charType, treeType>::swapParameter(const DagNode *oldP, const DagNode *newP) {
    
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
    else if (oldP == siteRates)
    {
        siteRates = static_cast<const TypedDagNode< std::vector< double > >* >( newP );
    }
    else if (oldP == topology)
    {
    	topology = static_cast<const TypedDagNode<Topology>* >( newP );
    }
    else
    {
        AbstractSiteHomogeneousMixtureCharEvoModel<charType, treeType>::swapParameter(oldP,newP);
    }
    
}

template<class charType, class treeType>
void RevBayesCore::DolloBranchHeterogeneousCharEvoModel<charType, treeType>::touchSpecialization( DagNode* affecter ) {
    
    // if the topology wasn't the culprit for the touch, then we just flag everything as dirty
    if ( affecter == heterogeneousClockRates )
    {
        const std::set<size_t> &indices = heterogeneousClockRates->getTouchedElementIndices();

        // maybe all of them have been touched or the flags haven't been set properly
        if ( indices.size() == 0 )
        {
            // just delegate the call
            AbstractSiteHomogeneousMixtureCharEvoModel<charType, treeType>::touchSpecialization( affecter );
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
            AbstractSiteHomogeneousMixtureCharEvoModel<charType, treeType>::touchSpecialization( affecter );
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
    else if ( affecter == topology )
    {
    	    computeAncestral();

	}
    else
	{

    	AbstractSiteHomogeneousMixtureCharEvoModel<charType, treeType>::touchSpecialization( affecter );
    }
    
}

template<class charType, class treeType>
void RevBayesCore::DolloBranchHeterogeneousCharEvoModel<charType, treeType>::computeAncestral( void ) {

	ancestralNodes.clear();
	ancestralNodes.resize(this->numPatterns);

	 const TopologyNode &root = this->tau->getValue().getRoot();

	for (size_t site = 0; site < this->numPatterns; ++site)
	{
		ancestralNodes[site].clear();

		computeAncestralMap(root,site);
		computeAncestralNodes(root,site);
	}
}

template<class charType, class treeType>
bool RevBayesCore::DolloBranchHeterogeneousCharEvoModel<charType, treeType>::computeAncestralMap(const TopologyNode& node, size_t site) {
	if(node.isTip()){
		AbstractTaxonData& taxon = this->getValue().getTaxonData( node.getName() );
		CharacterState &c = taxon.getCharacter(site);
		if(c != absentstate){
			ancestralMap[node.getIndex()] = true;
		}else{
			ancestralMap[node.getIndex()] = false;
		}
	}else{
		bool pleft = computeAncestralMap(node.getChild(0),site);
		bool pright = computeAncestralMap(node.getChild(1),site);
		ancestralMap[node.getIndex()] = pleft || pright;
	}

	return ancestralMap[node.getIndex()];
}

template<class charType, class treeType>
void RevBayesCore::DolloBranchHeterogeneousCharEvoModel<charType, treeType>::computeAncestralNodes(const TopologyNode& node, size_t site) {
	if(!node.isTip()){
		bool pleft = ancestralMap[node.getChild(0).getIndex()];
		bool pright = ancestralMap[node.getChild(1).getIndex()];
		if(pleft != pright){
			ancestralNodes[site].push_back(node.getIndex());
			if(pleft){
				computeAncestralNodes(node.getChild(0),site);
			}else{
				computeAncestralNodes(node.getChild(1),site);
			}
		}else if(pleft && pright){
			ancestralNodes[site].push_back(node.getIndex());
		}
	}
}

/*
template<class charType, class treeType>
void RevBayesCore::DolloBranchHeterogeneousCharEvoModel<charType, treeType>::redrawValue( void ) {

    // delete the old value first
    delete this->value;

    // create a new character data object
    this->value = new DiscreteCharacterData<charType>();

    // create a vector of taxon data
    std::vector< DiscreteTaxonData<charType> > taxa = std::vector< DiscreteTaxonData< charType > >( tau->getValue().getNumberOfNodes(), DiscreteTaxonData<charType>() );

    // first, simulate the per site rates
    RandomNumberGenerator* rng = GLOBAL_RNG;
    std::vector<size_t> perSiteRates;
    for ( size_t i = 0; i < numSites; ++i )
    {
        // draw the state
        double u = rng->uniform01();
        size_t rateIndex = (int)(u*numSiteRates);
        perSiteRates.push_back( rateIndex );
    }

    // simulate the root sequence
    const std::vector< double > &stationaryFreqs = getRootFrequencies();
    DiscreteTaxonData< charType > &root = taxa[ tau->getValue().getRoot().getIndex() ];
    for ( size_t i = 0; i < numSites; ++i )
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
    simulate( tau->getValue().getRoot(), taxa, perSiteRates );

    // add the taxon data to the character data
    for (size_t i = 0; i < tau->getValue().getNumberOfTips(); ++i)
    {
        this->value->addTaxonData( taxa[i] );
    }

    // compress the data and initialize internal variables
    this->compress();

}
*/

#endif

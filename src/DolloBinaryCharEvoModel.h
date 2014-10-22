#ifndef DolloBinaryCharEvoModel_H
#define DolloBinaryCharEvoModel_H

#include "BinaryCharEvoModel.h"
#include "DnaState.h"
#include "RateMatrix.h"
#include "RbVector.h"
#include "TopologyNode.h"
#include "TransitionProbabilityMatrix.h"
#include "TreeChangeEventListener.h"
#include "TypedDistribution.h"

namespace RevBayesCore {
    
    template<class treeType>
    class DolloBinaryCharEvoModel : public BinaryCharEvoModel<treeType> {
        
    public:
    	DolloBinaryCharEvoModel(const TypedDagNode< treeType > *p, const TypedDagNode<Topology> *t, bool c, size_t nSites, int type = 0);
        DolloBinaryCharEvoModel(const DolloBinaryCharEvoModel &n);                                                                                                //!< Copy constructor
        virtual                                            ~DolloBinaryCharEvoModel(void);                                                                   //!< Virtual destructor
        
        // public member functions
        DolloBinaryCharEvoModel*         					clone(void) const;
		void                                                swapParameter(const DagNode *oldP, const DagNode *newP);

		//virtual void                                        redrawValue(void);
        
    protected:

	   void                                                	computeRootLikelihood(size_t root, size_t l, size_t r);
	   void                                                	computeRootCorrection(size_t root, size_t l, size_t r);

	   void                                        			touchSpecialization(DagNode *toucher);
        
    private:
        void                                                computeAncestral(void);
        bool                                                computeAncestralMap(const TopologyNode& node, size_t site);
        void                                                computeAncestralNodes(const TopologyNode& node, size_t site);

        std::vector<std::vector<size_t> >                   ancestralNodes;
        std::map<size_t,bool>								ancestralMap;

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

template<class treeType>
RevBayesCore::DolloBinaryCharEvoModel<treeType>::DolloBinaryCharEvoModel(const TypedDagNode< treeType > *p, const TypedDagNode<Topology> *t, bool c, size_t nSites, int type) :
	BinaryCharEvoModel<treeType>(p, c , nSites, type),
	AbstractCharEvoModel<StandardState, treeType>(p, 2, c , nSites),
	topology(t)
{

	// initialize with default parameters
	this->addParameter( topology );
}


template<class treeType>
RevBayesCore::DolloBinaryCharEvoModel<treeType>::DolloBinaryCharEvoModel(const DolloBinaryCharEvoModel &d) :
	BinaryCharEvoModel<treeType>(d),
	AbstractCharEvoModel<StandardState, treeType>(d),
	ancestralNodes(d.ancestralNodes), topology(d.topology)
{
}


template<class treeType>
RevBayesCore::DolloBinaryCharEvoModel<treeType>::~DolloBinaryCharEvoModel( void ) {
    // We don't delete the parameters, because they might be used somewhere else too. The model needs to do that!
    
}


template<class treeType>
RevBayesCore::DolloBinaryCharEvoModel<treeType>* RevBayesCore::DolloBinaryCharEvoModel<treeType>::clone( void ) const {
    
    return new DolloBinaryCharEvoModel<treeType>( *this );
}

template<class treeType>
void RevBayesCore::DolloBinaryCharEvoModel<treeType>::computeRootLikelihood( size_t root, size_t left, size_t right)
{
    // reset the likelihood
    this->lnProb = 0.0;

    // get the root frequencies
    const std::vector<double> &f                    = this->getRootFrequencies();
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
			std::vector<double>::const_iterator p_site_left   = this->partialLikelihoods.begin() + this->activeLikelihood[nleft]*this->activeLikelihoodOffset + nleft*this->nodeOffset + site*this->siteOffset;
			std::vector<double>::const_iterator p_site_right  = this->partialLikelihoods.begin() + this->activeLikelihood[nright]*this->activeLikelihoodOffset + nright*this->nodeOffset + site*this->siteOffset;

			// iterate over all mixture categories
			for (size_t mixture = 0; mixture < this->numSiteRates; ++mixture)
			{
					size_t offset = mixture*this->mixtureOffset;

					std::vector<double>::const_iterator   p_site_left_j = p_site_left + offset;
					std::vector<double>::const_iterator   p_site_right_j = p_site_right + offset;

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
			}
		}
    }

    // sum the log-likelihoods for all sites together
    std::vector< size_t >::const_iterator patterns = this->patternCounts.begin();
    for (size_t site = 0; site < this->numPatterns; ++site, ++patterns)
    {
        this->lnProb += log( per_mixture_Likelihoods[site] ) * *patterns;
    }

    computeRootCorrection(root,left,right);

    this->lnProb -= this->lnCorrection;
}

template<class treeType>
void RevBayesCore::DolloBinaryCharEvoModel<treeType>::computeRootCorrection( size_t root, size_t left, size_t right)
{
	// reset the likelihood
	this->lnCorrection = 0.0;

	if(this->numCorrectionSites == 0)
		return;

	// get the root frequencies
	const std::vector<double> &f                    = this->getRootFrequencies();
	std::vector<double>::const_iterator f_end       = f.end();
	std::vector<double>::const_iterator f_begin     = f.begin();

	// create a vector for the per mixture likelihoods
	// we need this vector to sum over the different mixture likelihoods
	std::vector<double> per_mixture_Likelihoods = std::vector<double>(this->numCorrectionSites,0.0);

	// iterate over all sites
	for (size_t site = 0; site < this->numCorrectionSites; ++site){

		for (size_t anc = 0; anc < ancestralNodes[site].size(); anc++){
			size_t idx = ancestralNodes[site][anc];
			const TopologyNode &node = this->tau->getValue().getNode(idx);

			size_t nleft = node.getChild(0).getIndex();
			size_t nright = node.getChild(1).getIndex();

			// get the pointers to the partial likelihoods of the left and right subtree of this ancestral node
			std::vector<double>::const_iterator p_site_left   = this->partialLikelihoods.begin() + this->activeLikelihood[nleft]*this->activeLikelihoodOffset + nleft*this->nodeOffset + (site+this->numPatterns)*this->siteOffset;
			std::vector<double>::const_iterator p_site_right  = this->partialLikelihoods.begin() + this->activeLikelihood[nright]*this->activeLikelihoodOffset + nright*this->nodeOffset + (site+this->numPatterns)*this->siteOffset;

			for (size_t mixture = 0; mixture < this->numSiteRates; ++mixture)
			{

				size_t offset = mixture*this->mixtureOffset;

				std::vector<double>::const_iterator   p_site_left_j = p_site_left + offset;
				std::vector<double>::const_iterator   p_site_right_j = p_site_right + offset;

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

				double tmp = 0.0;

				// get the pointer to the stationary frequencies
				std::vector<double>::const_iterator f_j             = f_begin;

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
			}

		}

	} // end-for over all mixtures (=rate categories)

	// sum the log-likelihoods for all sites together
	for (size_t site = 0; site < this->numCorrectionSites; ++site)
	{
		this->lnCorrection += per_mixture_Likelihoods[site];
	}
	// normalize the log-probability
	this->lnCorrection /= this->numSiteRates;
	this->lnCorrection = this->numSites*log(1-this->lnCorrection);
}


template<class treeType>
void RevBayesCore::DolloBinaryCharEvoModel<treeType>::swapParameter(const DagNode *oldP, const DagNode *newP) {
    
    if (oldP == topology)
    {
    	topology = static_cast<const TypedDagNode<Topology>* >( newP );
    }
    else
    {
    	BinaryCharEvoModel<treeType>::swapParameter(oldP,newP);
    }
    
}

template<class treeType>
void RevBayesCore::DolloBinaryCharEvoModel<treeType>::touchSpecialization( DagNode* affecter ) {
    
    // if the topology wasn't the culprit for the touch, then we just flag everything as dirty
    if ( affecter == topology || affecter == this->dagNode )
    {
    	computeAncestral();

	}
    else
	{

    	BinaryCharEvoModel<treeType>::touchSpecialization(affecter);
    }
    
}

template<class treeType>
void RevBayesCore::DolloBinaryCharEvoModel<treeType>::computeAncestral( void ) {

	ancestralNodes.clear();
	ancestralNodes.resize(this->numPatterns+this->numCorrectionSites);

	const TopologyNode &root = this->tau->getValue().getRoot();

	for (size_t site = 0; site < this->numPatterns+this->numCorrectionSites; ++site)
	{
		ancestralNodes[site].clear();

		computeAncestralMap(root,site);
		computeAncestralNodes(root,site);
	}
}

template<class treeType>
bool RevBayesCore::DolloBinaryCharEvoModel<treeType>::computeAncestralMap(const TopologyNode& node, size_t site) {
	if(node.isTip()){
		if(site < this->numPatterns){
			AbstractTaxonData& taxon = this->getValue().getTaxonData( node.getName() );
			std::string c = taxon.getCharacter(site).getStringValue();
			if(c == "1"){
				ancestralMap[node.getIndex()] = true;
			}else{
				ancestralMap[node.getIndex()] = false;
			}
		}else{
			unsigned long c = this->correctionCharMatrix[node.getName()][site-this->numPatterns];
			if(c == 1){
				ancestralMap[node.getIndex()] = true;
			}else{
				ancestralMap[node.getIndex()] = false;
			}
		}
	}else{
		bool pleft = computeAncestralMap(node.getChild(0),site);
		bool pright = computeAncestralMap(node.getChild(1),site);
		ancestralMap[node.getIndex()] = pleft || pright;
	}

	return ancestralMap[node.getIndex()];
}

template<class treeType>
void RevBayesCore::DolloBinaryCharEvoModel<treeType>::computeAncestralNodes(const TopologyNode& node, size_t site) {
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
template<class treeType>
void RevBayesCore::DolloBinaryCharEvoModel<treeType>::redrawValue( void ) {

	this->dagNode->getLnProbability();
    // delete the old value first
    //delete this->value;

    // create a new character data object
    //this->value = new DiscreteCharacterData<charType>();

    // create a vector of taxon data
    //std::vector< DiscreteTaxonData<charType> > taxa = std::vector< DiscreteTaxonData< charType > >( tau->getValue().getNumberOfNodes(), DiscreteTaxonData<charType>() );

    // first, simulate a birth for each character
    // by sampling a node in proportion to its totalmass
    RandomNumberGenerator* rng = GLOBAL_RNG;
    std::vector<size_t> perSiteBirthNodes;
    //std::cerr << omega << std::endl;
    for ( size_t i = 0; i < this->numSites; ++i )
    {
    	double u = rng->uniform01();
    	double total = 0;
    	size_t i = 0;
    	while(total < u*omega){
    		total += totalmass[i++];
    	}
    	perSiteBirthNodes.push_back(i-1);
    }

    // next, simulate the per site rates
    std::vector<size_t> perSiteRates;
    for ( size_t i = 0; i < this->numSites; ++i )
    {
        // draw the state
        double u = rng->uniform01();
        size_t rateIndex = (int)(u*this->numSiteRates);
        perSiteRates.push_back( rateIndex );
    }

    // simulate the root sequence
    const std::vector< double > &stationaryFreqs = getRootFrequencies();
    DiscreteTaxonData< charType > &root = taxa[ tau->getValue().getRoot().getIndex() ];
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

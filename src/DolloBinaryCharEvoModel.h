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

		virtual void                                        redrawValue(void);
        
    protected:

	   void                                                	computeRootLikelihood(size_t root, size_t l, size_t r);
	   void                                                 computeRootCorrection(size_t root, size_t l, size_t r);
	   void                                                 computeInternalNodeCorrection(const TopologyNode &n, size_t nIdx, size_t l, size_t r);
	   void                                                 computeTipCorrection(const TopologyNode &node, size_t nIdx);

	   void                                        			touchSpecialization(DagNode *toucher);
	   void                                         		setCorrectionPatterns();
        
    private:
        void                                                computeAncestral(void);
        bool                                                computeAncestralMap(const TopologyNode& node, size_t site);
        void                                                computeAncestralNodes(const TopologyNode& node, size_t site);

        size_t 												simulate( const TopologyNode &node, std::vector<StandardState> &taxa, size_t birthNode, size_t rateIndex);


        std::vector<std::vector<size_t> >                   ancestralNodes;
        std::map<size_t,bool>								ancestralMap;

        const TypedDagNode<Topology>*						topology;
        double												omega;
        std::vector<double>									totalmass;
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

	this->numCorrectionSites = 2*(type > 0);

	if(this->numCorrectionSites > 0)
		totalmass.resize(this->tau->getValue().getNumberOfNodes());

	this->resizeLikelihoodVectors();
}


template<class treeType>
RevBayesCore::DolloBinaryCharEvoModel<treeType>::DolloBinaryCharEvoModel(const DolloBinaryCharEvoModel &d) :
	BinaryCharEvoModel<treeType>(d),
	AbstractCharEvoModel<StandardState, treeType>(d),
	ancestralNodes(d.ancestralNodes), topology(d.topology), totalmass(d.totalmass)
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

    this->lnProb -= log(this->numSites) + this->lnCorrection;
}

template<class treeType>
void RevBayesCore::DolloBinaryCharEvoModel<treeType>::setCorrectionPatterns()
{
	this->correctionCharMatrix.clear();
	this->correctionGapMatrix.clear();
}

template<class treeType>
void RevBayesCore::DolloBinaryCharEvoModel<treeType>::computeTipCorrection(const TopologyNode &node, size_t nodeIndex)
{
	if(this->numCorrectionSites == 0)
		return;

	std::vector<double>::iterator p_node = this->partialLikelihoods.begin() + this->activeLikelihood[nodeIndex]*this->activeLikelihoodOffset + nodeIndex*this->nodeOffset;

	// iterate over all mixture categories
	totalmass[nodeIndex] = 0;

	for (size_t mixture = 0; mixture < this->numSiteRates; ++mixture)
	{

		size_t offset = mixture*this->mixtureOffset + this->numPatterns*this->siteOffset;

		std::vector<double>::iterator     	u      = p_node + offset;

		//Probability of absence in all leaves descending from this node
		u[0] = 0.0;

		//Probability of presence in only one leaf descending from this node
		u[1] = 1.0;

		//Probability of presence in all leaves descending from this node
		u[2] = 1.0;

		//Probability of presence in all but one leaves descending from this node
		u[3] = 0.0;

	}
}

template<class treeType>
void RevBayesCore::DolloBinaryCharEvoModel<treeType>::computeInternalNodeCorrection(const TopologyNode &node, size_t nodeIndex, size_t left, size_t right)
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
	totalmass[nodeIndex] = 0;

	for (size_t mixture = 0; mixture < this->numSiteRates; ++mixture)
	{
		//get the 1->1 transition probabilities for each branch
		double tr = 	node.getBranchLength();
		this->updateTransitionProbabilities( nodeIndex, tr );
		double pr    		= this->transitionProbMatrices[mixture][1][1];

		double tij = left_node.getBranchLength();
		this->updateTransitionProbabilities( left, tij );
		double pij    		= this->transitionProbMatrices[mixture][1][1];

		double tik = right_node.getBranchLength();
		this->updateTransitionProbabilities( right, tik );
		double pik    		= this->transitionProbMatrices[mixture][1][1];

		//get the rate modifier for this mixture category
		double r = 1.0;
		if(this->rateVariationAcrossSites == true)
			r = this->siteRates->getValue()[mixture];

		size_t offset = mixture*this->mixtureOffset + this->numPatterns*this->siteOffset;

		std::vector<double>::iterator     			u       = p_node + offset;
		std::vector<double>::const_iterator     	u_j  	= p_left + offset;
		std::vector<double>::const_iterator     	u_k 	= p_right + offset;

		//Probability of absence in all leaves descending from this node
		u[0] = (1 - pij*(1 - u_j[0]))*(1 - pik*(1 - u_k[0]));

		//Probability of presence in only one leaf descending from this node
		u[1] = (1 - pij*(1 - u_j[0]))*pik*u_k[1] + (1 - pik*(1 - u_k[0]))*pij*u_j[1];

		//Probability of presence in all leaves descending from this node
		u[2] 	= pij*u_j[2]*pik*u_k[2];

		bool jl = left_node.isTip();
		bool kl = right_node.isTip();

		//Probability of presence in all but one leaves descending from this node
		u[3] = (1-pij*(1-u_j[3]) + !jl*(pij-1))*pik*u_k[2]
			 + (1-pik*(1-u_k[3]) + !kl*(pik-1))*pij*u_j[2];

		double prob = 1.0;
		if(this->type & NO_ABSENT_SITES)
			prob -= u[0];
		if(this->type & NO_SINGLETON_GAINS)
			prob -= u[1];
		if(this->type & NO_PRESENT_SITES)
			prob -= u[2];
		//If both of this node's children are leaves, then u[1] = u[3]
		if((this->type & NO_SINGLETON_LOSSES) && !(jl && kl))
			prob -= u[3];

		// correct rounding errors
		if(prob < 0)
			prob = 0;

		//Nicholls and Gray 2003
		totalmass[nodeIndex] += prob*(1-pr)/r;
	}

}

template<class treeType>
void RevBayesCore::DolloBinaryCharEvoModel<treeType>::computeRootCorrection( size_t root, size_t left, size_t right)
{
	// reset the likelihood
	omega = 0.0;
	if(this->numCorrectionSites == 0){
		//If all site-patterns can, in principle, be observed
		//then omega is just the tree length
		if(this->rateVariationAcrossSites == true){
			std::vector<double> r = this->siteRates->getValue();
			for(size_t i = 0; i < r.size(); i++){
				omega += this->tau->getValue().getTreeLength()/r[i];
			}
		}else{
			omega = this->tau->getValue().getTreeLength();
		}
	}else{
		// get the pointers to the partial likelihoods for this node and the two descendant subtrees
		std::vector<double>::const_iterator   p_left  = this->partialLikelihoods.begin() + this->activeLikelihood[left]*this->activeLikelihoodOffset + left*this->nodeOffset;
		std::vector<double>::const_iterator   p_right = this->partialLikelihoods.begin() + this->activeLikelihood[right]*this->activeLikelihoodOffset + right*this->nodeOffset;
		std::vector<double>::iterator         p_node  = this->partialLikelihoods.begin() + this->activeLikelihood[root]*this->activeLikelihoodOffset + root*this->nodeOffset;

		const TopologyNode &node = this->tau->getValue().getNode(root);
		const TopologyNode &left_node = this->tau->getValue().getNode(left);
		const TopologyNode &right_node = this->tau->getValue().getNode(right);

		// iterate over all mixture categories
		totalmass[root] = 0;

		for (size_t mixture = 0; mixture < this->numSiteRates; ++mixture)
		{
			//get the 1->1 transition probabilities for each branch
			this->updateTransitionProbabilities( left, this->tau->getValue().getNode(left).getBranchLength() );
			double pij    		= this->transitionProbMatrices[mixture][1][1];

			this->updateTransitionProbabilities( right, this->tau->getValue().getNode(right).getBranchLength() );
			double pik    		= this->transitionProbMatrices[mixture][1][1];

			//get the rate modifier for this mixture category
			double r = 1.0;
			if(this->rateVariationAcrossSites == true)
				r = this->siteRates->getValue()[mixture];

			size_t offset = mixture*this->mixtureOffset + this->numPatterns*this->siteOffset;

			std::vector<double>::iterator     			u       = p_node + offset;
			std::vector<double>::const_iterator     	u_j  	= p_left + offset;
			std::vector<double>::const_iterator     	u_k 	= p_right + offset;

			//Probability of absence in all leaves descending from this node
			u[0] = (1 - pij*(1 - u_j[0]))*(1 - pik*(1 - u_k[0]));

			//Probability of presence in only one leaf descending from this node
			u[1] = (1 - pij*(1 - u_j[0]))*pik*u_k[1] + (1 - pik*(1 - u_k[0]))*pij*u_j[1];

			//Probability of presence in all leaves descending from this node
			u[2] 	= pij*u_j[2]*pik*u_k[2];

			bool jl = left_node.isTip();
			bool kl = right_node.isTip();

			//Probability of presence in all but one leaves descending from this node
			u[3] = (1-pij*(1-u_j[3]) + !jl*(pij-1))*pik*u_k[2]
				 + (1-pik*(1-u_k[3]) + !kl*(pik-1))*pij*u_j[2];

			double prob = 1.0;
			if(this->type & NO_ABSENT_SITES)
				prob -= u[0];
			if(this->type & NO_SINGLETON_GAINS)
				prob -= u[1];
			if(this->type & NO_PRESENT_SITES)
				prob -= u[2];
			//If both of this node's children are leaves, then u[1] = u[3]
			if((this->type & NO_SINGLETON_LOSSES) && !(jl && kl))
				prob -= u[3];

			// correct rounding errors
			if(prob < 0)
				prob = 0;

			//Nicholls and Gray
			totalmass[root] += prob/r;
		}

		for(size_t i = 0; i < this->tau->getValue().getNumberOfNodes(); i++){
			omega += totalmass[i];
		}
	}

	this->lnCorrection = this->numSites*log(omega);
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
	ancestralNodes.resize(this->numPatterns);

	const TopologyNode &root = this->tau->getValue().getRoot();

	for (size_t site = 0; site < this->numPatterns; ++site)
	{
		ancestralNodes[site].clear();

		computeAncestralMap(root,site);
		computeAncestralNodes(root,site);
	}
}

template<class treeType>
bool RevBayesCore::DolloBinaryCharEvoModel<treeType>::computeAncestralMap(const TopologyNode& node, size_t site) {
	if(node.isTip()){
			AbstractTaxonData& taxon = this->getValue().getTaxonData( node.getName() );
			std::string c = taxon.getCharacter(site).getStringValue();
			if(c == "1"){
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


template<class treeType>
void RevBayesCore::DolloBinaryCharEvoModel<treeType>::redrawValue( void ) {

	this->computeLnProbability();

    // delete the old value first
    delete this->value;

    // create a new character data object
    this->value = new DiscreteCharacterData<StandardState>();

    size_t numTips = this->tau->getValue().getNumberOfTips();

    for (size_t t = 0; t < numTips; ++t)
	{
		DiscreteTaxonData<StandardState> data;
		data.setTaxonName( this->tau->getValue().getNode(t).getName() );
		this->value->addTaxonData(data);
	}

    // first, simulate a birth for each character
    // by sampling a node in proportion to its totalmass
    RandomNumberGenerator* rng = GLOBAL_RNG;

    size_t cap = 0;
    for ( size_t i = 0; i < this->numSites; ++i )
    {
    	double u = rng->uniform01();
    	double total = 0;
    	size_t birthNode = 0;
    	while(total < u*omega){
    		total += totalmass[birthNode++];
    	}
    	birthNode--;

        // draw the state
        u = rng->uniform01();
        size_t rateIndex = (int)(u*this->numSiteRates);

        std::vector<StandardState> taxa(this->tau->getValue().getNumberOfNodes(), StandardState());

		// recursively simulate the sequences
		size_t numLeaves = simulate( this->tau->getValue().getRoot(), taxa, birthNode, rateIndex );

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
		for (size_t t = 0; t < this->tau->getValue().getNumberOfTips(); ++t)
		{
			this->value->getTaxonData(t).addCharacter(taxa[t]);
		}
    }

    std::cerr << (double)cap/this->numSites << std::endl;

    // compress the data and initialize internal variables
    this->compress();

}

template<class treeType>
size_t RevBayesCore::DolloBinaryCharEvoModel<treeType>::simulate( const TopologyNode &node, std::vector<StandardState> &taxa, size_t birthNode, size_t rateIndex) {

	size_t numLeaves = 0;
    // get the children of the node
    const std::vector<TopologyNode*>& children = node.getChildren();

    // get the sequence of this node
    size_t nodeIndex = node.getIndex();
    StandardState &parentState = taxa[ nodeIndex ];

    if(node.getIndex() == birthNode)
    	parentState++;

    if(parentState.getState() == 2 && node.isTip())
    	numLeaves++;

    // simulate the sequence for each child
    RandomNumberGenerator* rng = GLOBAL_RNG;
    for (std::vector< TopologyNode* >::const_iterator it = children.begin(); it != children.end(); ++it)
    {
        const TopologyNode &child = *(*it);

        // update the transition probability matrix
        this->updateTransitionProbabilities( child.getIndex(), child.getBranchLength() );

        StandardState &childState = taxa[ child.getIndex() ];

		double *pij = this->transitionProbMatrices[ rateIndex ][ 1 ];

		if(parentState.getState() == 2){
			if(rng->uniform01() < pij[1]){
				childState++;
				numLeaves += simulate( child, taxa, birthNode, rateIndex );
			}
		}else{
			numLeaves += simulate( child, taxa, birthNode, rateIndex );
		}
	}

    return numLeaves;

}

#endif

#ifndef DolloBinaryPartitionModel_H
#define DolloBinaryPartitionModel_H

#include "BinaryPartitionModel.h"
#include "DnaState.h"
#include "RateMatrix.h"
#include "RbVector.h"
#include "DistributionPoisson.h"
#include "DistributionGamma.h"
#include "RbStatisticsHelper.h"
#include "TopologyNode.h"
#include "TransitionProbabilityMatrix.h"
#include "TreeChangeEventListener.h"
#include "TypedDistribution.h"

namespace RevBayesCore {
    
    template<class treeType>
    class DolloBinaryPartitionModel : public BinaryPartitionModel<treeType> {
        
    public:
    	DolloBinaryPartitionModel(const TypedDagNode< treeType > *p, const TypedDagNode<Topology> *t, bool c, size_t nSites, int type = 0);
        DolloBinaryPartitionModel(const DolloBinaryPartitionModel &n);                                                                                                //!< Copy constructor
        virtual                                            ~DolloBinaryPartitionModel(void);                                                                   //!< Virtual destructor
        
        // public member functions
        DolloBinaryPartitionModel*         					clone(void) const;
		void                                                swapParameter(const DagNode *oldP, const DagNode *newP);

		virtual void                                        redrawValue(void);
        
    protected:

	   void                                                	computeRootLikelihood(size_t root, size_t l, size_t r);
	   void                                                 computeRootCorrections(size_t root, size_t l, size_t r);
	   void                                                 computeInternalNodeCorrections(const TopologyNode &n, size_t nIdx, size_t l, size_t r);
	   void                                                 computeTipCorrections(const TopologyNode &node, size_t nIdx);

	   void                                        			touchSpecialization(DagNode *toucher);
	   void                                         		setCorrectionPatterns();
        
    private:
        void                                                computeAncestral(void);
        bool                                                computeAncestralMap(const TopologyNode& node, size_t site);
        void                                                computeAncestralNodes(const TopologyNode& node, size_t site);

        size_t 												simulate( const TopologyNode &node, std::vector<StandardState> &taxa, size_t rateIndex);


        std::vector<std::vector<size_t> >                   ancestralNodes;
        std::map<size_t,bool>								ancestralMap;

        bool												ancestralNeedsRecomputing;

        const TypedDagNode<Topology>*									topology;
        std::vector<std::vector<std::vector<std::vector<double> > >	>	totalmassPerPartition;
        std::vector<double>												omegaPerPartition;

        std::vector<double>									probAbsencePerPartition;
        std::map<size_t, size_t>							patterns2partitions;
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
RevBayesCore::DolloBinaryPartitionModel<treeType>::DolloBinaryPartitionModel(const TypedDagNode< treeType > *p, const TypedDagNode<Topology> *t, bool c, size_t nSites, int type) :
	BinaryPartitionModel<treeType>(p, c , nSites, type),
	AbstractCharEvoModel<StandardState, treeType>(p, 2, c , nSites),
	topology(t),
	ancestralNodes(this->numPatterns, std::vector<size_t>())
{

	// initialize with default parameters
	this->addParameter( topology );

	ancestralNeedsRecomputing = true;
}


template<class treeType>
RevBayesCore::DolloBinaryPartitionModel<treeType>::DolloBinaryPartitionModel(const DolloBinaryPartitionModel &d) :
	BinaryPartitionModel<treeType>(d),
	AbstractCharEvoModel<StandardState, treeType>(d),
	ancestralNodes(d.ancestralNodes),
	topology(d.topology),
	totalmassPerPartition(d.totalmassPerPartition),
	omegaPerPartition(d.omegaPerPartition),
	probAbsencePerPartition(d.probAbsencePerPartition),
	patterns2partitions(d.patterns2partitions)
{
	ancestralNeedsRecomputing = true;
}


template<class treeType>
RevBayesCore::DolloBinaryPartitionModel<treeType>::~DolloBinaryPartitionModel( void ) {
    // We don't delete the parameters, because they might be used somewhere else too. The model needs to do that!
    
}


template<class treeType>
RevBayesCore::DolloBinaryPartitionModel<treeType>* RevBayesCore::DolloBinaryPartitionModel<treeType>::clone( void ) const {
    
    return new DolloBinaryPartitionModel<treeType>( *this );
}

template<class treeType>
void RevBayesCore::DolloBinaryPartitionModel<treeType>::computeRootLikelihood( size_t root, size_t left, size_t right)
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

    if(ancestralNeedsRecomputing)
    	computeAncestral();

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
						prob_birth *= 1-exp(-r*node.getBranchLength());
				}


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

	computeRootCorrections(root,left,right);

    // sum the log-likelihoods for all sites together
    std::vector< size_t >::const_iterator patterns = this->patternCounts.begin();

	for (size_t site = 0; site < this->numPatterns; ++site, ++patterns)
	{
		// if a site likelihood is zero, it means it is a constant absence site
		// therefore, the likelihood is equal to probAbsence
		if(per_mixture_Likelihoods[site] == 0.0)
			per_mixture_Likelihoods[site] = probAbsencePerPartition[patterns2partitions[site]];

		this->lnProb += log(per_mixture_Likelihoods[site])* *patterns;
	}

	// normalize mixtures
	this->lnProb -= this->numSites*log( this->numSiteRates );

	//apply correction
	for(size_t partition = 0; partition < this->numCorrectionPartitions; partition++){
		this->lnProb -= this->numSitesPerPartition[partition]*this->perPartitionLnCorrections[partition];
	}

    //reset probAbsence
    probAbsencePerPartition = std::vector<double>(this->numCorrectionPartitions, 0.0);
}



template<class treeType>
void RevBayesCore::DolloBinaryPartitionModel<treeType>::setCorrectionPatterns()
{
	this->numCorrectionPartitions = 0;
	this->numSitesPerPartition.clear();
	this->numGapsPerPartition.clear();
	this->perMixtureCorrections.clear();
	this->perPartitionLnCorrections.clear();
	this->patternsPerPartition.clear();

	std::map<std::string, size_t> gapPatterns;
	std::map<std::string, size_t> gapPatterns2partitions;
	for(size_t site = 0; site < this->numPatterns; site++){
		std::stringstream gapstream;
		for(std::map<std::string, std::vector<bool> >::iterator it = this->gapMatrix.begin(); it != this->gapMatrix.end(); it++)
			gapstream << it->second[site];

		std::string gaps = gapstream.str();

		if(gapPatterns.find(gaps) == gapPatterns.end()){
			if(count(gaps.begin(), gaps.end(), '1') == 0){
				this->patternsPerPartition.insert(this->patternsPerPartition.begin(),site);
				patterns2partitions[site] = 0;
				for(std::map<size_t, size_t >::iterator pp = patterns2partitions.begin(); pp != patterns2partitions.end(); pp++)
					pp->second++;
			}else{
				patterns2partitions[site] = this->patternsPerPartition.size();
				this->patternsPerPartition.push_back(site);
			}
			gapPatterns2partitions[gaps] = patterns2partitions[site];
		}else{
			patterns2partitions[site] = gapPatterns2partitions[gaps];
		}

		gapPatterns[gaps] += this->patternCounts[site];
	}

	this->numCorrectionPartitions = this->patternsPerPartition.size();

	for(std::map<std::string, size_t>::iterator it = gapPatterns.begin(); it != gapPatterns.end(); it++){
		size_t gaps = count(it->first.begin(),it->first.end(), '1');

		if(gaps > 0){
			this->numSitesPerPartition.push_back(it->second);
			this->numGapsPerPartition.push_back(gaps);
		}else{
			this->numSitesPerPartition.insert(this->numSitesPerPartition.begin(), it->second);
			this->numGapsPerPartition.insert(this->numGapsPerPartition.begin(), 0);
		}
	}

    totalmassPerPartition = std::vector<std::vector<std::vector<std::vector<double> > > >(this->numCorrectionPartitions,
    						std::vector<std::vector<std::vector<double> > >(this->tau->getValue().getNumberOfNodes(),
    						std::vector<std::vector<double> >(this->numSiteRates,
    						std::vector<double>(2,0.0))));

	omegaPerPartition = std::vector<double>(this->numCorrectionPartitions,0.0);
	probAbsencePerPartition = std::vector<double>(this->numCorrectionPartitions,0.0);

	this->numCorrectionSitesPerPartition = std::vector<size_t>(this->numCorrectionPartitions,5);
	this->perPartitionLnCorrections = std::vector<double>(this->numCorrectionPartitions,0.0);

	this->correctionCharMatrix.clear();
	this->correctionGapMatrix.clear();
}

template<class treeType>
void RevBayesCore::DolloBinaryPartitionModel<treeType>::computeTipCorrections(const TopologyNode &node, size_t nodeIndex)
{
	if(this->type == NONE)
		return;

	std::vector<double>::iterator p_node = this->partialLikelihoods.begin() + this->activeLikelihood[nodeIndex]*this->activeLikelihoodOffset + nodeIndex*this->nodeOffset;

	// iterate over all mixture categories

	for (size_t partition = 0; partition < this->numCorrectionPartitions; ++partition)
	{
		for (size_t mixture = 0; mixture < this->numSiteRates; ++mixture)
		{

			size_t offset = mixture*this->mixtureOffset + this->numPatterns*this->siteOffset + 5*partition;

			std::vector<double>::iterator     	u      = p_node + offset;

			//Probability of absence in all observed leaves descending from this node
			u[0] = 0.0;

			//Probability of presence in only one observed leaf descending from this node
			u[1] = 1.0;

			//Probability of presence in all observed leaves descending from this node
			u[2] = 1.0;

			//Probability of presence in all but one observed leaves descending from this node
			u[3] = 0.0;

			//boolean: are all leaves descending from this node gaps?
			u[4] = this->gapMatrix[node.getName()][this->patternsPerPartition[partition]];

			totalmassPerPartition[partition][nodeIndex][mixture][0] = 0;
			totalmassPerPartition[partition][nodeIndex][mixture][1] = 0;

		}
	}
}

template<class treeType>
void RevBayesCore::DolloBinaryPartitionModel<treeType>::computeInternalNodeCorrections(const TopologyNode &node, size_t nodeIndex, size_t left, size_t right)
{
	if(this->type == NONE)
		return;

	// get the pointers to the partial likelihoods for this node and the two descendant subtrees
	std::vector<double>::const_iterator   p_left  = this->partialLikelihoods.begin() + this->activeLikelihood[left]*this->activeLikelihoodOffset + left*this->nodeOffset;
	std::vector<double>::const_iterator   p_right = this->partialLikelihoods.begin() + this->activeLikelihood[right]*this->activeLikelihoodOffset + right*this->nodeOffset;
	std::vector<double>::iterator         p_node  = this->partialLikelihoods.begin() + this->activeLikelihood[nodeIndex]*this->activeLikelihoodOffset + nodeIndex*this->nodeOffset;

	const TopologyNode &left_node = this->tau->getValue().getNode(left);
	const TopologyNode &right_node = this->tau->getValue().getNode(right);

	// iterate over all mixture categories

	for (size_t partition = 0; partition < this->numCorrectionPartitions; ++partition)
	{
		for (size_t mixture = 0; mixture < this->numSiteRates; ++mixture)
		{
			size_t offset = mixture*this->mixtureOffset + this->numPatterns*this->siteOffset + 5*partition;

			std::vector<double>::iterator     			u       = p_node + offset;
			std::vector<double>::const_iterator     	u_j  	= p_left + offset;
			std::vector<double>::const_iterator     	u_k 	= p_right + offset;

			//boolean: are all leaves descending from this node gaps?
			u[4] = u_j[4] && u_k[4];

			if(u[4]){

				u[0] = 0.0;
				u[1] = 0.0;
				u[2] = 0.0;
				u[3] = 0.0;

				totalmassPerPartition[partition][nodeIndex][mixture][0] = 0.0;
				totalmassPerPartition[partition][nodeIndex][mixture][1] = 0.0;

			}else{
				//get the 1->1 transition probabilities for each branch
				this->updateTransitionProbabilities( nodeIndex, node.getBranchLength() );
				double pr    		= this->transitionProbMatrices[mixture][1][1];

				this->updateTransitionProbabilities( left, left_node.getBranchLength() );
				double pij    		= this->transitionProbMatrices[mixture][1][1];

				this->updateTransitionProbabilities( right, right_node.getBranchLength() );
				double pik    		= this->transitionProbMatrices[mixture][1][1];

				//Probability of absence in all observed leaves descending from this node
				u[0] = (1 - !u_j[4]*pij*(1 - u_j[0]))*(1 - !u_k[4]*pik*(1 - u_k[0]));

				//Probability of presence in only one observed leaf descending from this node
				u[1] = (1 - !u_j[4]*pij*(1 - u_j[0]))*pik*u_k[1]*!u_k[4]
					 + (1 - !u_k[4]*pik*(1 - u_k[0]))*pij*u_j[1]*!u_j[4];

				//Probability of presence in all observed leaves descending from this node
				u[2] 	= (pij*u_j[2] + u_j[4]*(1 - pij*u_j[2]))*(pik*u_k[2] + u_k[4]*(1 - pik*u_k[2]));

				bool jl = left_node.isTip();
				bool kl = right_node.isTip();

				//Probability of presence in all but one observed leaves descending from this node
				u[3] = (1 - !u_j[4]*(pij*(1-u_j[3]) + !jl*(pij-1)))*pik*u_k[2]*!u_k[4]
					 + (1 - !u_k[4]*(pik*(1-u_k[3]) + !kl*(pik-1)))*pij*u_j[2]*!u_j[4];

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

				double r = 1.0;
				if(this->rateVariationAcrossSites == true)
					r = this->siteRates->getValue()[mixture];

				//Nicholls and Gray 2003
				totalmassPerPartition[partition][nodeIndex][mixture][0] = prob;
				totalmassPerPartition[partition][nodeIndex][mixture][1] = (1-pr)/r;

				probAbsencePerPartition[partition] += u[0]*(1-pr)/r;
			}
		}
	}
}

template<class treeType>
void RevBayesCore::DolloBinaryPartitionModel<treeType>::computeRootCorrections( size_t root, size_t left, size_t right)
{
	omegaPerPartition = std::vector<double>(this->numCorrectionPartitions, 0.0);

	for (size_t partition = 0; partition < this->numCorrectionPartitions; ++partition)
	{

		if(this->type == NONE){
			//If all site-patterns can, in principle, be observed
			//then omega is just the tree length

			omegaPerPartition[partition] = this->tau->getValue().getTreeLength()*this->numSiteRates;
		}else{
			// get the pointers to the partial likelihoods for this node and the two descendant subtrees
			std::vector<double>::const_iterator   p_left  = this->partialLikelihoods.begin() + this->activeLikelihood[left]*this->activeLikelihoodOffset + left*this->nodeOffset;
			std::vector<double>::const_iterator   p_right = this->partialLikelihoods.begin() + this->activeLikelihood[right]*this->activeLikelihoodOffset + right*this->nodeOffset;
			std::vector<double>::iterator         p_node  = this->partialLikelihoods.begin() + this->activeLikelihood[root]*this->activeLikelihoodOffset + root*this->nodeOffset;

			const TopologyNode &left_node = this->tau->getValue().getNode(left);
			const TopologyNode &right_node = this->tau->getValue().getNode(right);

			// iterate over all mixture categories
			for (size_t mixture = 0; mixture < this->numSiteRates; ++mixture)
			{
				size_t offset = mixture*this->mixtureOffset + this->numPatterns*this->siteOffset;

				std::vector<double>::iterator     			u       = p_node + offset;
				std::vector<double>::const_iterator     	u_j  	= p_left + offset;
				std::vector<double>::const_iterator     	u_k 	= p_right + offset;

				//boolean: are all leaves descending from this node gaps?
				u[4] = u_j[4] && u_k[4];

				if(u[4]){

					u[0] = 0.0;
					u[1] = 0.0;
					u[2] = 0.0;
					u[3] = 0.0;

					totalmassPerPartition[partition][root][mixture][0] = 0.0;
					totalmassPerPartition[partition][root][mixture][1] = 0.0;

				}else{
					//get the 1->1 transition probabilities for each branch
					this->updateTransitionProbabilities( left, left_node.getBranchLength() );
					double pij    		= this->transitionProbMatrices[mixture][1][1];

					this->updateTransitionProbabilities( right, right_node.getBranchLength() );
					double pik    		= this->transitionProbMatrices[mixture][1][1];

					//get the rate modifier for this mixture category
					double r = 1.0;
					if(this->rateVariationAcrossSites == true)
						r = this->siteRates->getValue()[mixture];

					//Probability of absence in all observed leaves descending from this node
					u[0] = (1 - !u_j[4]*pij*(1 - u_j[0]))*(1 - !u_k[4]*pik*(1 - u_k[0]));

					//Probability of presence in only one observed leaf descending from this node
					u[1] = (1 - !u_j[4]*pij*(1 - u_j[0]))*pik*u_k[1]*!u_k[4]
						 + (1 - !u_k[4]*pik*(1 - u_k[0]))*pij*u_j[1]*!u_j[4];

					//Probability of presence in all observed leaves descending from this node
					u[2] 	= (pij*u_j[2] + u_j[4]*(1 - pij*u_j[2]))*(pik*u_k[2] + u_k[4]*(1 - pik*u_k[2]));

					bool jl = left_node.isTip();
					bool kl = right_node.isTip();

					//Probability of presence in all but one observed leaves descending from this node
					u[3] = (1 - !u_j[4]*(pij*(1-u_j[3]) + !jl*(pij-1)))*pik*u_k[2]*!u_k[4]
						 + (1 - !u_k[4]*(pik*(1-u_k[3]) + !kl*(pik-1)))*pij*u_j[2]*!u_j[4];

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

					totalmassPerPartition[partition][root][mixture][0] = prob;
					totalmassPerPartition[partition][root][mixture][1] = 1.0/r;

					// just probability of absence site
					probAbsencePerPartition[partition] += u[0]/r;
				}

				for(size_t i = 0; i < this->tau->getValue().getNumberOfNodes(); i++)
					for(size_t mixture = 0; mixture < this->numSiteRates; mixture++)
						omegaPerPartition[partition] += totalmassPerPartition[partition][i][mixture][0]*totalmassPerPartition[partition][i][mixture][1];
			}
		}

		this->perPartitionLnCorrections[partition] = log(omegaPerPartition[partition]) - log(this->numSiteRates);
	}
}

template<class treeType>
void RevBayesCore::DolloBinaryPartitionModel<treeType>::swapParameter(const DagNode *oldP, const DagNode *newP) {
    
    if (oldP == topology)
    {
    	topology = static_cast<const TypedDagNode<Topology>* >( newP );
    }
    else
    {
    	BinaryPartitionModel<treeType>::swapParameter(oldP,newP);
    }
    
}

template<class treeType>
void RevBayesCore::DolloBinaryPartitionModel<treeType>::touchSpecialization( DagNode* affecter ) {
    
    // if the topology wasn't the culprit for the touch, then we just flag everything as dirty
    if ( affecter == topology || affecter == this->dagNode )
    	ancestralNeedsRecomputing = true;

	if(affecter == this->siteRates){
		for(size_t partition = 0; partition < this->numCorrectionPartitions; partition++)
			for(size_t node = 0; node < this->tau->getValue().getNumberOfNodes(); node++)
				totalmassPerPartition[partition][node] = std::vector<std::vector<double> >(this->numSiteRates,std::vector<double>(2,0.0));
	}

	BinaryPartitionModel<treeType>::touchSpecialization(affecter);
    
}

template<class treeType>
void RevBayesCore::DolloBinaryPartitionModel<treeType>::computeAncestral( void ) {

	ancestralNodes = std::vector<std::vector<size_t> >(this->numPatterns,std::vector<size_t>());

	const TopologyNode &root = this->tau->getValue().getRoot();

	for (size_t site = 0; site < this->numPatterns; ++site)
	{
		ancestralNodes[site].clear();
		computeAncestralMap(root,site);
		computeAncestralNodes(root,site);
	}

	ancestralNeedsRecomputing = false;
}

template<class treeType>
bool RevBayesCore::DolloBinaryPartitionModel<treeType>::computeAncestralMap(const TopologyNode& node, size_t site) {
	if(node.isTip()){
		unsigned long c = this->charMatrix[node.getName()][site];

		if(c == 1 + this->usingAmbiguousCharacters){
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
void RevBayesCore::DolloBinaryPartitionModel<treeType>::computeAncestralNodes(const TopologyNode& node, size_t site) {
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
void RevBayesCore::DolloBinaryPartitionModel<treeType>::redrawValue( void ) {

	//this->computeLnProbability();

    // delete the old value first
    delete this->value;

    // create a new character data object
    this->value = new DiscreteCharacterData<StandardState>();

    size_t numTips = this->tau->getValue().getNumberOfTips();
    size_t numNodes = this->tau->getValue().getNumberOfNodes();

    std::vector< DiscreteTaxonData<StandardState> > taxa = std::vector< DiscreteTaxonData<StandardState> >(numTips, DiscreteTaxonData<StandardState>() );

    RandomNumberGenerator* rng = GLOBAL_RNG;

    if(this->type != NONE){
		// first sample a birth rate (lambda)
		// from the marginal posterior lambda | N ~ Gamma(N, omega )
    	// given lambda ~ 1/lambda
		double lambda = RbStatistics::Gamma::rv(this->numSites, omegaPerPartition.front(), *rng);

		// then resample numSites from Poisson( lambda*omega )
		this->numSites = RbStatistics::Poisson::rv( lambda*omegaPerPartition.front(), *rng);

		//std::cerr << "lambda: " << lambda << "\tN: " << this->numSites << std::endl;
    }

    std::vector<size_t> birthNodes;
    std::vector<size_t> perSiteRates;
    for(size_t site = 0; site < this->numSites; site++){
		double u = rng->uniform01()*omegaPerPartition.front();
		double total = 0.0;

		// simulate a birth for each character
		// by sampling nodes in proportion to their share of totalmass
		size_t birthNode = 0;
		while(total < u){
			for(size_t mixture = 0; mixture < this->numSiteRates; mixture++)
				total += totalmassPerPartition.front()[birthNode][mixture][0]*totalmassPerPartition.front()[birthNode][mixture][1];
			if(total < u)
				birthNode++;
		}

		// then sample a rate category conditional on survival from this node
		// by sampling in proportion to the per mixture surivival probs

		total = 0.0;
		for(size_t mixture = 0; mixture < this->numSiteRates; mixture++)
			total += totalmassPerPartition.front()[birthNode][mixture][0];


		u = rng->uniform01()*total;
		size_t rateIndex = 0;

		double tmp = 0.0;
		while(tmp < u){
			tmp += totalmassPerPartition.front()[birthNode][rateIndex][0];
			if(tmp < u)
				rateIndex++;
		}

		birthNodes.push_back(birthNode);
		perSiteRates.push_back(rateIndex);
    }

    // then simulate the data conditional on these birth nodes and site rates
    for ( size_t i = 0; i < this->numSites; ++i )
    {
        // draw the state
        size_t rateIndex = perSiteRates[i];
        const TopologyNode &birthNode = this->tau->getValue().getNode(birthNodes[i]);

        std::vector<StandardState> siteData(numNodes, StandardState(this->absentState,this->absentState+this->presentState));

        siteData[birthNodes[i]].setState(this->presentState);

		// recursively simulate the sequences
		size_t numLeaves = simulate( birthNode, siteData, rateIndex );

		//if(numLeaves > 0)
		//std::cerr << numLeaves << std::endl;
		// reject any site-patterns that are forbidden by our data acquisition settings
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
size_t RevBayesCore::DolloBinaryPartitionModel<treeType>::simulate( const TopologyNode &node, std::vector<StandardState> &taxa, size_t rateIndex) {

	size_t numLeaves = 0;
    // get the children of the node
    const std::vector<TopologyNode*>& children = node.getChildren();

	// simulate the sequence for each child
	RandomNumberGenerator* rng = GLOBAL_RNG;
	for (std::vector< TopologyNode* >::const_iterator it = children.begin(); it != children.end(); ++it)
	{
		const TopologyNode &child = *(*it);

		// update the transition probability matrix
		this->updateTransitionProbabilities( child.getIndex(), child.getBranchLength() );

		StandardState &childState = taxa[ child.getIndex() ];

		double pij = this->transitionProbMatrices[ rateIndex ][1][1];

		if(rng->uniform01() < pij){
			childState.setState(this->presentState);

			if(child.isTip())
				numLeaves++;
			else
				numLeaves += simulate( child, taxa, rateIndex );
		}
	}

    return numLeaves;

}

#endif

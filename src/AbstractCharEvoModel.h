#ifndef AbstractCharEvoModel_H
#define AbstractCharEvoModel_H

#include "AbstractCharacterData.h"
#include "DiscreteTaxonData.h"
#include "DnaState.h"
#include "RateMatrix.h"
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
    class AbstractCharEvoModel : public TypedDistribution< AbstractCharacterData >, public TreeChangeEventListener {
        
    public:
        // Note, we need the size of the alignment in the constructor to correctly simulate an initial state
        AbstractCharEvoModel(const TypedDagNode<treeType> *t, size_t nChars, bool c, size_t nSites);
        AbstractCharEvoModel(const AbstractCharEvoModel &n);                                                                                          //!< Copy constructor
        virtual                                                            ~AbstractCharEvoModel(void);                                                              //!< Virtual destructor
        
        // public member functions
        // pure virtual
        virtual AbstractCharEvoModel*                 						clone(void) const = 0;
        
        // virtual (you need to overwrite this method if you have additional parameters)
        virtual void                                                        swapParameter(const DagNode *oldP, const DagNode *newP);                                //!< Implementation of swaping paramoms
        
        // non-virtual
        virtual double                                                      computeLnProbability(void);
        void                                                                fireTreeChangeEvent(const TopologyNode &n);                                             //!< The tree has changed and we want to know which part.
        void                                                                setValue(AbstractCharacterData *v);                                                   //!< Set the current value, e.g. attach an observation (clamp)
        virtual void                                                        redrawValue(void) = 0;
        void                                                                reInitialized(void);
        const treeType*														getTree() const;
        
    protected:
        // helper method for this and derived classes
        void                                                                recursivelyFlagNodeDirty(const TopologyNode& n);
        virtual void                                                        resizeLikelihoodVectors(void) = 0;
        
        // virtual methods that may be overwritten, but then the derived class should call this methods
        virtual void                                                        keepSpecialization(DagNode* affecter);
        virtual void                                                        restoreSpecialization(DagNode *restorer);
        virtual void                                                        touchSpecialization(DagNode *toucher);
        
        // pure virtual methods
        virtual void                                                        computeRootLikelihood(size_t root, size_t l, size_t r) = 0;
        virtual void                                                        computeInternalNodeLikelihood(const TopologyNode &n, size_t nIdx, size_t l, size_t r) = 0;
        virtual void                                                        computeTipLikelihood(const TopologyNode &node, size_t nIdx) = 0;
        virtual void                                                        updateTransitionProbabilities(size_t nodeIdx, double brlen) = 0;
        
        // members
        double                                                              lnProb;
        size_t                                                              numSites;
        const size_t                                                        numChars;
        const TypedDagNode<treeType>*                                       tau;
        std::vector<TransitionProbabilityMatrix>                            transitionProbMatrices;
        
        // the likelihoods
        std::vector<double>                                                 partialLikelihoods;
        std::vector<size_t>                                                 activeLikelihood;
        
        // the data
        std::map<std::string,std::vector<unsigned long> >                   charMatrix;
        std::map<std::string,std::vector<bool> >                            gapMatrix;
        std::map<size_t,size_t>												patternMap;
        std::vector<size_t>                                                 patternCounts;
        size_t                                                              numPatterns;
        bool                                                                compressed;
        
        // convenience variables available for derived classes too
        std::vector<bool>                                                   changedNodes;
        std::vector<bool>                                                   dirtyNodes;

        // offsets for nodes
        size_t                                                              activeLikelihoodOffset;
        size_t                                                              nodeOffset;
        size_t																siteOffset;
        
        // flags
        bool                                                                usingAmbiguousCharacters;
        bool                                                                treatUnknownAsGap;
        bool                                                                treatAmbiguousAsGaps;
        
    protected:
        // private methods
        void                                                                compress(void);
        void                                                                fillLikelihoodVector(const TopologyNode &n, size_t nIdx);

    
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
RevBayesCore::AbstractCharEvoModel<charType, treeType>::AbstractCharEvoModel(const TypedDagNode<treeType> *t, size_t nChars, bool c, size_t nSites) : TypedDistribution< AbstractCharacterData >(  new DiscreteCharacterData<charType>() ),
    lnProb(0.0),
	numSites( nSites ),
    numChars( nChars ),
    tau( t ), 
    charMatrix(), 
    gapMatrix(),
    patternCounts(),
    numPatterns( numSites ),
    compressed( c ),
	activeLikelihood( std::vector<size_t>(tau->getValue().getNumberOfNodes(), 0) ),
    changedNodes( std::vector<bool>(tau->getValue().getNumberOfNodes(),false) ),
    dirtyNodes( std::vector<bool>(tau->getValue().getNumberOfNodes(), true) ),
    usingAmbiguousCharacters( true ),
    treatUnknownAsGap( true ),
    treatAmbiguousAsGaps( true )
{
    // add the parameters to the parents list
    this->addParameter( tau );
    tau->getValue().getTreeChangeEventHandler().addListener( this );
    
    this->resizeLikelihoodVectors();
}


template<class charType, class treeType>
RevBayesCore::AbstractCharEvoModel<charType, treeType>::AbstractCharEvoModel(const AbstractCharEvoModel &n) : TypedDistribution< AbstractCharacterData >( n ),
    lnProb(n.lnProb),
	numSites( n.numSites ),
    numChars( n.numChars ),
    tau( n.tau ), 
    transitionProbMatrices( n.transitionProbMatrices ),
    activeLikelihood( n.activeLikelihood ),
    charMatrix( n.charMatrix ), 
    gapMatrix( n.gapMatrix ), 
    patternCounts( n.patternCounts ),
	patternMap(n.patternMap),
    numPatterns( n.numPatterns ),
    compressed( n.compressed ),
    changedNodes( n.changedNodes ),
    dirtyNodes( n.dirtyNodes ),
    usingAmbiguousCharacters( n.usingAmbiguousCharacters ),
    treatUnknownAsGap( n.treatUnknownAsGap ),
    treatAmbiguousAsGaps( n.treatAmbiguousAsGaps ),
	activeLikelihoodOffset(n.activeLikelihoodOffset),
	partialLikelihoods(n.partialLikelihoods),
	nodeOffset(n.nodeOffset),
	siteOffset(n.siteOffset)
{
    // parameters are automatically copied
    
    tau->getValue().getTreeChangeEventHandler().addListener( this );
}


template<class charType, class treeType>
RevBayesCore::AbstractCharEvoModel<charType, treeType>::~AbstractCharEvoModel( void ) {
    // We don't delete the paramoms, because they might be used somewhere else too. The model needs to do that!
    
    // remove myself from the tree listeners
    if ( tau != NULL ) 
    {
        // TODO: this needs to be implemented (Sebastian)
        tau->getValue().getTreeChangeEventHandler().removeListener( this );
    }
}


template<class charType, class treeType>
RevBayesCore::AbstractCharEvoModel<charType, treeType>* RevBayesCore::AbstractCharEvoModel<charType, treeType>::clone( void ) const
{
    
    return new AbstractCharEvoModel<charType, treeType>( *this );
}

template<class charType, class treeType>
const treeType* RevBayesCore::AbstractCharEvoModel<charType, treeType>::getTree( void ) const
{

    return &(tau->getValue());
}


template<class charType, class treeType>
void RevBayesCore::AbstractCharEvoModel<charType, treeType>::compress( void )
{
    
//    compressed = false;
    
    charMatrix.clear();
    gapMatrix.clear();
    patternCounts.clear();
    patternMap.clear();
    numPatterns = 0;
    
    size_t totalNumSites = value->getNumberOfCharacters();
    numSites = value->getNumberOfIncludedCharacters();

    // check whether there are ambiguous characters (besides gaps)
    bool ambiguousCharacters = false;
    // find the unique site patterns and compute their respective frequencies
    std::vector<TopologyNode*> nodes = tau->getValue().getNodes();
    for (size_t site = 0; site < totalNumSites; ++site)
    {
    	if(!value->isCharacterExcluded(site)){
			for (std::vector<TopologyNode*>::iterator it = nodes.begin(); it != nodes.end(); ++it)
			{
				if ( (*it)->isTip() )
				{
					// \todo modify this so that the distribution is actually defined on discrete character data
					AbstractTaxonData& taxon = value->getTaxonData( (*it)->getName() );
					DiscreteCharacterState &c = static_cast<DiscreteCharacterState &>( taxon.getCharacter(site) );

					// if we treat unknown characters as gaps and this is an unknown character then we change it
					// because we might then have a pattern more
					if ( treatAmbiguousAsGaps && c.isAmbiguous() )
					{
						c.setGapState( true );
					}
					else if ( treatUnknownAsGap && c.getNumberOfStates() == c.getNumberObservedStates() )
					{
						c.setGapState( true );
					}
					else if ( !c.isGapState() && c.isAmbiguous() )
					{
						ambiguousCharacters = true;
						break;
					}
				}
			}

			// break the loop if there was an ambiguous character
			if ( ambiguousCharacters )
			{
				break;
			}
    	}
    }
    // set the global variable if we use ambiguous characters
    usingAmbiguousCharacters = ambiguousCharacters;

    
    std::vector<bool> unique(numSites, true);
    // compress the character matrix if we're asked to
    if ( compressed ) 
    {
        // find the unique site patterns and compute their respective frequencies
        std::map<std::string,size_t> patterns;
        size_t thisSite = 0;
        for (size_t site = 0; site < totalNumSites; ++site)
        {
        	if(!value->isCharacterExcluded(site)){
				// create the site pattern
				std::string pattern = "";
				for (std::vector<TopologyNode*>::iterator it = nodes.begin(); it != nodes.end(); ++it)
				{
					if ( (*it)->isTip() )
					{
						AbstractTaxonData& taxon = value->getTaxonData( (*it)->getName() );
						CharacterState &c = taxon.getCharacter(site);
						pattern += c.getStringValue();
					}
				}
				// check if we have already seen this site pattern
				std::map<std::string, size_t>::const_iterator index = patterns.find( pattern );
				if ( index != patterns.end() )
				{
					// we have already seen this pattern
					// increase the frequency counter
					patternCounts[ index->second ]++;
					patternMap[thisSite] = index->second;

					// obviously this site isn't unique nor the first encounter
					unique[site] = false;
				}
				else
				{
					patternMap[thisSite] = numPatterns;

					// create a new pattern frequency counter for this pattern
					patternCounts.push_back(1);

					// insert this pattern with the corresponding index in the map
					patterns.insert( std::pair<std::string,size_t>(pattern,numPatterns) );

					// increase the pattern counter
					numPatterns++;

					// flag that this site is unique (or the first occurence of this pattern)
					unique[site] = true;
				}
				thisSite++;
        	}
        }
    } 
    else 
    {
        // we do not compress
        numPatterns = numSites;
        patternCounts = std::vector<size_t>(numSites,1);
        for (size_t site = 0; site < numSites; ++site)
        	patternMap[site] = site;
    }
    
    
    // allocate and fill the cells of the matrices
    for (std::vector<TopologyNode*>::iterator it = nodes.begin(); it != nodes.end(); ++it) 
    {
        if ( (*it)->isTip() ) 
        {
            std::string name = (*it)->getName();
            AbstractTaxonData& taxon = value->getTaxonData( name );

            // resize the column
            charMatrix[name].resize(numPatterns);
            gapMatrix[name].resize(numPatterns);
            size_t patternIndex = 0;
            for (size_t site = 0; site < totalNumSites; ++site)
            {
            	if(!value->isCharacterExcluded(site)){
					// only add this site if it is unique
					if ( unique[site] )
					{
						charType &c = static_cast<charType &>( taxon.getCharacter(site) );
						gapMatrix[name][patternIndex] = c.isGapState();

						if ( ambiguousCharacters )
						{
							// we use the actual state
							charMatrix[name][patternIndex] = c.getState();
						}
						else
						{
							// we use the index of the state
							size_t index = 0;
							unsigned long state = c.getState();
							state >>= 1;

							while ( state != 0 ) // there are still observed states left
							{

								// remove this state from the observed states
								state >>= 1;

								// increment the index
								++index;
							} // end-while over all observed states for this character

							charMatrix[name][patternIndex] = index;
						}

						// increase the pattern index
						patternIndex++;
					}
            	}
            }
        }
    }

    // finally we resize the partial likelihood vectors to the new pattern counts
    this->resizeLikelihoodVectors();
    
}


template<class charType, class treeType>
double RevBayesCore::AbstractCharEvoModel<charType, treeType>::computeLnProbability( void )
{
    // compute the ln probability by recursively calling the probability calculation for each node
    const TopologyNode &root = tau->getValue().getRoot();
    
    // we start with the root and then traverse down the tree
    size_t rootIndex = root.getIndex();

    // only necessary if the root is actually dirty
    if ( dirtyNodes[rootIndex] )
    {
                
        // mark as computed
        dirtyNodes[rootIndex] = false;


        // start by filling the likelihood vector for the two children of the root
        const TopologyNode &left = root.getChild(0);
        size_t leftIndex = left.getIndex();
        fillLikelihoodVector( left, leftIndex );
        const TopologyNode &right = root.getChild(1);
        size_t rightIndex = right.getIndex();
        fillLikelihoodVector( right, rightIndex );
        // compute the likelihood of the root
        computeRootLikelihood( rootIndex, leftIndex, rightIndex );
    }

    return this->lnProb;
}


template<class charType, class treeType>
void RevBayesCore::AbstractCharEvoModel<charType, treeType>::fillLikelihoodVector(const TopologyNode &node, size_t nodeIndex)
{    
    
    // check for recomputation
    if ( dirtyNodes[nodeIndex] ) 
    {
        // mark as computed
        dirtyNodes[nodeIndex] = false;
        
        if ( node.isTip() ) 
        {
            // this is a tip node
            // compute the likelihood for the tip and we are done
            computeTipLikelihood(node, nodeIndex);
        } 
        else 
        {
            // this is an internal node
            // start by filling the likelihood vector for the two children of this node
            const TopologyNode &left = node.getChild(0);
            size_t leftIndex = left.getIndex();
            fillLikelihoodVector( left, leftIndex );
            const TopologyNode &right = node.getChild(1);
            size_t rightIndex = right.getIndex();
            fillLikelihoodVector( right, rightIndex );
            
            // now compute the likelihoods of this internal node
            computeInternalNodeLikelihood(node,nodeIndex,leftIndex,rightIndex);
        }
    }
}



template<class charType, class treeType>
void RevBayesCore::AbstractCharEvoModel<charType, treeType>::fireTreeChangeEvent( const RevBayesCore::TopologyNode &n ) {
    
    // call a recursive flagging of all node above (closer to the root) and including this node
    recursivelyFlagNodeDirty( n );
    
}


template<class charType, class treeType>
void RevBayesCore::AbstractCharEvoModel<charType, treeType>::keepSpecialization( DagNode* affecter ) {
    
    // reset all flags
    for (std::vector<bool>::iterator it = this->dirtyNodes.begin(); it != this->dirtyNodes.end(); ++it) 
    {
        (*it) = false;
    }
    
    for (std::vector<bool>::iterator it = this->changedNodes.begin(); it != this->changedNodes.end(); ++it) 
    {
        (*it) = false;
    }
    
}



template<class charType, class treeType>
void RevBayesCore::AbstractCharEvoModel<charType, treeType>::recursivelyFlagNodeDirty( const RevBayesCore::TopologyNode &n ) {
    
    // we need to flag this node and all ancestral nodes for recomputation
    size_t index = n.getIndex();
    
    // if this node is already dirty, the also all the ancestral nodes must have been flagged as dirty
    if ( !dirtyNodes[index] ) 
    {
        // the root doesn't have an ancestor
        if ( !n.isRoot() ) 
        {
            recursivelyFlagNodeDirty( n.getParent() );
        }
        
        // set the flag
        dirtyNodes[index] = true;
        
        // if we previously haven't touched this node, then we need to change the active likelihood pointer
        if ( !changedNodes[index] ) 
        {
            activeLikelihood[index] ^= 1;
            changedNodes[index] = true;
        }
        
    }
    
}


template<class charType, class treeType>
void RevBayesCore::AbstractCharEvoModel<charType, treeType>::reInitialized( void ) {
    
    // we need to recompress because the tree may have changed
    compress();
}


template<class charType, class treeType>
void RevBayesCore::AbstractCharEvoModel<charType, treeType>::restoreSpecialization( DagNode* affecter ) {
    
    // reset the flags
    for (std::vector<bool>::iterator it = dirtyNodes.begin(); it != dirtyNodes.end(); ++it) 
    {
        (*it) = false;
    }
    
    // restore the active likelihoods vector
    for (size_t index = 0; index < changedNodes.size(); ++index) 
    {
        // we have to restore, that means if we have changed the active likelihood vector
        // then we need to revert this change
        if ( changedNodes[index] ) 
        {
            activeLikelihood[index] = (activeLikelihood[index] == 0 ? 1 : 0);
        }
        
        // set all flags to false
        changedNodes[index] = false;
    }
    
}


template<class charType, class treeType>
void RevBayesCore::AbstractCharEvoModel<charType, treeType>::setValue(AbstractCharacterData *v) {
    
    // delegate to the parent class
    TypedDistribution< AbstractCharacterData >::setValue(v);
    
    this->compress();
}


template<class charType, class treeType>
void RevBayesCore::AbstractCharEvoModel<charType, treeType>::swapParameter(const DagNode *oldP, const DagNode *newP) {
    
    // we only have the topology here as the parameter
    if (oldP == tau) 
    {
        tau->getValue().getTreeChangeEventHandler().removeListener( this );
        tau = static_cast<const TypedDagNode<treeType>* >( newP );
        tau->getValue().getTreeChangeEventHandler().addListener( this );
    }
    
}


template<class charType, class treeType>
void RevBayesCore::AbstractCharEvoModel<charType, treeType>::touchSpecialization( DagNode* affecter ) {
    
    // if the topology wasn't the culprit for the touch, then we just flag everything as dirty
    if ( affecter != tau )
    {
        for (std::vector<bool>::iterator it = dirtyNodes.begin(); it != dirtyNodes.end(); ++it) 
        {
            (*it) = true;
        }
        
        // flip the active likelihood pointers
        for (size_t index = 0; index < changedNodes.size(); ++index) 
        {
            if ( !changedNodes[index] ) 
            {
                activeLikelihood[index] = (activeLikelihood[index] == 0 ? 1 : 0);
                changedNodes[index] = true;
            }
        }
    }
    
}

template<class charType, class treeType>
void RevBayesCore::AbstractCharEvoModel<charType, treeType>::resizeLikelihoodVectors( void ) {

    // we resize the partial likelihood vectors to the new dimensions
    size_t numNodes = tau->getValue().getNumberOfNodes();

    partialLikelihoods.clear();
    partialLikelihoods.resize(2*numNodes*numPatterns*numChars);

    transitionProbMatrices = std::vector<TransitionProbabilityMatrix>(1, TransitionProbabilityMatrix(numChars) );

    // set the offsets for easier iteration through the likelihood vector
    activeLikelihoodOffset      =  numNodes*numPatterns*numChars;
    nodeOffset                  =  numPatterns*numChars;
    siteOffset                  =  numChars;
}

#endif

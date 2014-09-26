 ///
//  TreeHistory.h
//  rb_mlandis
//
//  Created by Michael Landis on 3/28/14.
//  Copyright (c) 2014 Michael Landis. All rights reserved.
//

#ifndef __rb_mlandis__AbstractTreeHistoryCtmc__
#define __rb_mlandis__AbstractTreeHistoryCtmc__

#include "AbstractCharacterData.h"
#include "BranchHistory.h"
#include "DiscreteTaxonData.h"
#include "DiscreteCharacterData.h"
#include "DiscreteCharacterState.h"
#include "DnaState.h"
#include "RandomNumberFactory.h"
#include "RandomNumberGenerator.h"
#include "RateMatrix.h"
#include "TopologyNode.h"
#include "TransitionProbabilityMatrix.h"
#include "Tree.h"
#include "TreeChangeEventListener.h"
#include "TypedDistribution.h"

#include <cmath>

namespace RevBayesCore {

    template<class charType, class treeType>
    class AbstractTreeHistoryCtmc : public TypedDistribution< AbstractCharacterData >, public TreeChangeEventListener {
        
    public:
        // Note, we need the size of the alignment in the constructor to correctly simulate an initial state
        AbstractTreeHistoryCtmc(const TypedDagNode<treeType> *t, size_t nChars, size_t nSites, bool useAmbigChar=false);
        AbstractTreeHistoryCtmc(const AbstractTreeHistoryCtmc &n);                                                                           //!< Copy constructor
        virtual                                                            ~AbstractTreeHistoryCtmc(void);                                   //!< Virtual destructor
        
        // public member functions
        // pure virtual
        virtual AbstractTreeHistoryCtmc*                                    clone(void) const = 0;                                           //!< Create an independent clone
        virtual void                                                        redrawValue(void) = 0;
        virtual void                                                        initializeValue(void) = 0;
        virtual double                                                      samplePathStart(const TopologyNode& node, const std::set<size_t>& indexSet) = 0;
        virtual double                                                      samplePathEnd(const TopologyNode& node, const std::set<size_t>& indexSet) = 0;
        virtual double                                                      samplePathHistory(const TopologyNode& node, const std::set<size_t>& indexSet) = 0;
        
        // virtual (you need to overwrite this method if you have additional parameters)
        virtual void                                                        swapParameter(const DagNode *oldP, const DagNode *newP);         //!< Implementation of swaping paramoms
        
        // non-virtual
        double                                                              computeLnProbability(void);
        void                                                                fireTreeChangeEvent(const TopologyNode &n);                      //!< The tree has changed and we want to know which part.
//        BranchHistory&                                                      getHistory(size_t idx);
//        const BranchHistory&                                                getHistory(size_t idx) const;
        BranchHistory&                                                      getHistory(const TopologyNode& nd);
        const BranchHistory&                                                getHistory(const TopologyNode& nd) const;

        std::vector<BranchHistory*>                                         getHistories(void);
        const std::vector<BranchHistory*>&                                  getHistories(void) const;
//        void                                                                setHistory(const BranchHistory& bh, size_t idx);
        void                                                                setHistory(const BranchHistory& bh, const TopologyNode& nd);
        void                                                                setHistories(const std::vector<BranchHistory*>& bh);
        void                                                                setValue(AbstractCharacterData *v);                              //!< Set the current value, e.g. attach an observation (clamp)
        
        //virtual const std::vector<double>&                                  getTipProbs(const TopologyNode& nd) = 0;
        //virtual const std::vector<std::vector<double> >&                    getTipProbs(void) = 0;
        
        virtual void                                                        simulate(void);

        
    protected:
        // helper method for this and derived classes
        void                                                                flagNodeDirty(const TopologyNode& n);
        void                                                                resizeLikelihoodVectors(void);
        
        // virtual methods that may be overwritten, but then the derived class should call this methods
        virtual void                                                        keepSpecialization(DagNode* affecter);
        virtual void                                                        restoreSpecialization(DagNode *restorer);
        virtual void                                                        touchSpecialization(DagNode *toucher);
        
        // pure virtual methods
        // pure virtual methods
		virtual void                                                        computeRootLikelihood(size_t root, size_t l, size_t r) = 0;
		virtual void                                                        computeInternalNodeLikelihood(const TopologyNode &n, size_t nIdx, size_t l, size_t r) = 0;
		virtual void                                                        computeTipLikelihood(const TopologyNode &node, size_t nIdx) = 0;
        virtual const std::vector<double>&                                  getRootFrequencies(void) = 0;
        
        // members
        double                                                              lnProb;
        const size_t                                                        numChars;
        size_t                                                              numSites;
        size_t                                                              numSiteRates;
        const TypedDagNode<treeType>*                                       tau;
        
        // the likelihoods
        std::vector<size_t>                                                 activeLikelihood;
        std::vector<std::vector<double> >                                   historyLikelihoods;
        
        // the data
        std::vector<std::vector<unsigned long> >                            charMatrix;
        std::vector<std::vector<bool> >                                     gapMatrix;
        std::vector<BranchHistory*>                                         histories;
        
        // convenience variables available for derived classes too
        std::vector<bool>                                                   changedNodes;
        std::vector<bool>                                                   dirtyNodes;
        
        // flags
        bool                                                                usingAmbiguousCharacters;
        bool                                                                treatUnknownAsGap;
        bool                                                                treatAmbiguousAsGaps;
        bool                                                                tipsInitialized;
        
    private:
        // private methods
        void                                                                fillLikelihoodVector(const TopologyNode &n);
        void                                                                initializeHistoriesVector(void);
        virtual void                                                        simulate(const TopologyNode& node, DiscreteCharacterData< charType > &taxa);
        

        
        
    };

}

template<class charType, class treeType>
RevBayesCore::AbstractTreeHistoryCtmc<charType, treeType>::AbstractTreeHistoryCtmc(const TypedDagNode<treeType> *t, size_t nChars, size_t nSites, bool useAmbigChar) : TypedDistribution< AbstractCharacterData >(  new DiscreteCharacterData<charType>() ),
numChars( nChars ),
numSites( nSites ),
numSiteRates( 1 ),
tau( t ),
activeLikelihood( std::vector<size_t>(tau->getValue().getNumberOfNodes(), 0) ),
historyLikelihoods(),
charMatrix(),
gapMatrix(),
histories(),
changedNodes( std::vector<bool>(tau->getValue().getNumberOfNodes(),false) ),
dirtyNodes( std::vector<bool>(tau->getValue().getNumberOfNodes(), true) ),
usingAmbiguousCharacters( useAmbigChar ),
treatUnknownAsGap( true ),
treatAmbiguousAsGaps( true ),
tipsInitialized( false )
{
    
    // add the paramoms to the parents list
    this->addParameter( tau );
    tau->getValue().getTreeChangeEventHandler().addListener( this );
    
    // initialize histories
    initializeHistoriesVector();
    
}


template<class charType, class treeType>
RevBayesCore::AbstractTreeHistoryCtmc<charType, treeType>::AbstractTreeHistoryCtmc(const AbstractTreeHistoryCtmc &n) : TypedDistribution< AbstractCharacterData >( n ),
numChars( n.numChars ),
numSites( n.numSites ),
numSiteRates( n.numSiteRates ),
tau( n.tau ),
activeLikelihood( n.activeLikelihood ),
historyLikelihoods( n.historyLikelihoods ),
charMatrix( n.charMatrix ),
gapMatrix( n.gapMatrix ),
histories( n.histories ),
changedNodes( n.changedNodes ),
dirtyNodes( n.dirtyNodes ),
usingAmbiguousCharacters( n.usingAmbiguousCharacters ),
treatUnknownAsGap( n.treatUnknownAsGap ),
treatAmbiguousAsGaps( n.treatAmbiguousAsGaps ),
tipsInitialized( n.tipsInitialized )
{
    // parameters are automatically copied
    tau->getValue().getTreeChangeEventHandler().addListener( this );
    
}


template<class charType, class treeType>
RevBayesCore::AbstractTreeHistoryCtmc<charType, treeType>::~AbstractTreeHistoryCtmc( void ) {
    // We don't delete the paramoms, because they might be used somewhere else too. The model needs to do that!
    
    // remove myself from the tree listeners
    if ( tau != NULL )
    {
        // TODO: this needs to be implemented (Sebastian)
        tau->getValue().getTreeChangeEventHandler().removeListener( this );
    }
}


template<class charType, class treeType>
RevBayesCore::AbstractTreeHistoryCtmc<charType, treeType>* RevBayesCore::AbstractTreeHistoryCtmc<charType, treeType>::clone( void ) const
{
    
    return new AbstractTreeHistoryCtmc<charType, treeType>( *this );
}

template<class charType, class treeType>
double RevBayesCore::AbstractTreeHistoryCtmc<charType, treeType>::computeLnProbability( void )
{
    this->lnProb = 0.0;

    
    const std::vector<TopologyNode*>& nodes = tau->getValue().getNodes();

//    std::cout << "recompute lnProb    ";
//    for (size_t i = 0; i < nodes.size(); i++)
//    {
//        if (changedNodes[i] || dirtyNodes[i])
//        {
//            std::cout << i << ":" << (changedNodes[i] ? "1" : "0") << (dirtyNodes[i] ? "1" : "0") << " ";
//        }
//    }
//    std::cout << "\n";
    
    for (size_t i = 0; i < nodes.size(); i++)
    {
        const TopologyNode& nd = *nodes[i];
        size_t nodeIndex = nd.getIndex();
//        dirtyNodes[nodeIndex] = true; // to be commented
//        activeLikelihood[nodeIndex] = 0; // to be commented
        fillLikelihoodVector(nd);
        this->lnProb += historyLikelihoods[ activeLikelihood[nodeIndex] ][nodeIndex];
    }
    //std::cout << this->lnProb << "\n";
//    computeRootLikelihood(tau->getValue().getRoot());
    return this->lnProb;
}


template<class charType, class treeType>
void RevBayesCore::AbstractTreeHistoryCtmc<charType, treeType>::fillLikelihoodVector(const TopologyNode &node)
{
    size_t nodeIndex = node.getIndex();
    if (!dirtyNodes[nodeIndex])
        return;

    // compute
    double lnL = computeInternalNodeLikelihood(node);
    historyLikelihoods[ activeLikelihood[nodeIndex] ][nodeIndex] = lnL;
    
    // mark as computed
    dirtyNodes[nodeIndex] = false;

}


template<class charType, class treeType>
void RevBayesCore::AbstractTreeHistoryCtmc<charType, treeType>::fireTreeChangeEvent( const RevBayesCore::TopologyNode &n ) {
    
    // call a recursive flagging of all node above (closer to the root) and including this node
    flagNodeDirty(n);
    
//    size_t idx = n.getIndex();
//    std::cout << "fireTreeChangeEvent() " << idx << "  " << (changedNodes[idx] ? "1" : "0") << (dirtyNodes[idx] ? "1" : "0") << "\n";
}

template<class charType, class treeType>
const RevBayesCore::BranchHistory&  RevBayesCore::AbstractTreeHistoryCtmc<charType, treeType>::getHistory(const TopologyNode& nd) const
{
    return histories[nd.getIndex()];
}

template<class charType, class treeType>
RevBayesCore::BranchHistory&  RevBayesCore::AbstractTreeHistoryCtmc<charType, treeType>::getHistory(const TopologyNode& nd)
{
    return *histories[nd.getIndex()];
}


//template<class charType, class treeType>
//const RevBayesCore::BranchHistory&  RevBayesCore::AbstractTreeHistoryCtmc<charType, treeType>::getHistory(size_t idx) const
//{
//    
//    return histories[idx];
//}
//
//template<class charType, class treeType>
//RevBayesCore::BranchHistory&  RevBayesCore::AbstractTreeHistoryCtmc<charType, treeType>::getHistory(size_t idx)
//{
//    return *histories[idx];
//}

template<class charType, class treeType>
const std::vector<RevBayesCore::BranchHistory*>& RevBayesCore::AbstractTreeHistoryCtmc<charType, treeType>::getHistories(void) const
{
    return histories;
}

template<class charType, class treeType>
std::vector<RevBayesCore::BranchHistory*> RevBayesCore::AbstractTreeHistoryCtmc<charType, treeType>::getHistories(void)
{
    return histories;
}



template<class charType, class treeType>
void RevBayesCore::AbstractTreeHistoryCtmc<charType, treeType>::initializeHistoriesVector( void ) {
    
    std::vector<TopologyNode*> nodes = tau->getValue().getNodes();
    histories.resize(nodes.size());
    for (size_t i = 0; i < nodes.size(); i++)
    {
        TopologyNode* nd = nodes[i];
        histories[nd->getIndex()] = new BranchHistory(numSites,numChars,nd->getIndex());
    }
    
    historyLikelihoods.resize(2);
    for (size_t i = 0; i < 2; i++)
        historyLikelihoods[i].resize(nodes.size(), 0.0);

}


template<class charType, class treeType>
void RevBayesCore::AbstractTreeHistoryCtmc<charType, treeType>::keepSpecialization( DagNode* affecter ) {
    
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
void RevBayesCore::AbstractTreeHistoryCtmc<charType, treeType>::flagNodeDirty( const RevBayesCore::TopologyNode &n ) {
    
    // we need to flag this node and all ancestral nodes for recomputation
    size_t index = n.getIndex();
    
    // if this node is already dirty, the also all the ancestral nodes must have been flagged as dirty
    if ( !dirtyNodes[index] )
    {
        // set the flag
        dirtyNodes[index] = true;
        
        // if we previously haven't touched this node, then we need to change the active likelihood pointer
        if ( !changedNodes[index] )
        {
            activeLikelihood[index] = (activeLikelihood[index] == 0 ? 1 : 0);
            //activeLikelihood[index] = 0;
            changedNodes[index] = true;
        }
        
    }
    
}

template<class charType, class treeType>
void RevBayesCore::AbstractTreeHistoryCtmc<charType, treeType>::resizeLikelihoodVectors( void ) {

    // we resize the partial likelihood vectors to the new dimensions
    delete [] partialLikelihoods;
    partialLikelihoods = new double[2*tau->getValue().getNumberOfNodes()*numSiteRates*numPatterns*numChars];

    transitionProbMatrices = std::vector<TransitionProbabilityMatrix>(numSiteRates, TransitionProbabilityMatrix(numChars) );

    // set the offsets for easier iteration through the likelihood vector
    activeLikelihoodOffset      =  tau->getValue().getNumberOfNodes()*numSiteRates*numPatterns*numChars;
    nodeOffset                  =  numSiteRates*numPatterns*numChars;
    mixtureOffset               =  numPatterns*numChars;
    siteOffset                  =  numChars;
}

template<class charType, class treeType>
void RevBayesCore::AbstractTreeHistoryCtmc<charType, treeType>::restoreSpecialization( DagNode* affecter ) {
    
    // reset the flags
    for (std::vector<bool>::iterator it = dirtyNodes.begin(); it != dirtyNodes.end(); ++it)
    {
        (*it) = false;
    }
    
    //std::cout << "affecter " << affecter->getName() << "\n";
    
    // restore the active likelihoods vector
    for (size_t index = 0; index < changedNodes.size(); ++index)
    {
        // we have to restore, that means if we have changed the active likelihood vector
        // then we need to revert this change
        if ( changedNodes[index] )
        {

            activeLikelihood[index] = (activeLikelihood[index] == 0 ? 1 : 0);
            //activeLikelihood[index] = 0;
//            if (affecter->getName() == "ctmc") std::cout << index << " " << activeLikelihood[index] << "\n";
        }
        
        // set all flags to false
        changedNodes[index] = false;
    }
    
    
    return;
}

template<class charType, class treeType>
void RevBayesCore::AbstractTreeHistoryCtmc<charType, treeType>::setHistory(const BranchHistory& bh, const TopologyNode& nd)
{
    histories[ nd.getIndex() ] = new BranchHistory(bh);
}


template<class charType, class treeType>
void RevBayesCore::AbstractTreeHistoryCtmc<charType, treeType>::setHistories(const std::vector<BranchHistory*>& bh)
{
    for (size_t i = 0; i < bh.size(); i++)
        histories[i] = bh[i];

}

template<class charType, class treeType>
void RevBayesCore::AbstractTreeHistoryCtmc<charType, treeType>::setValue(AbstractCharacterData *v) {
    
    // delegate to the parent class
    TypedDistribution< AbstractCharacterData >::setValue(v);
}


template<class charType, class treeType>
void RevBayesCore::AbstractTreeHistoryCtmc<charType, treeType>::simulate(void)
{
    ;
}

template<class charType, class treeType>
//void RevBayesCore::AbstractTreeHistoryCtmc<charType, treeType>::simulate( const TopologyNode &node, std::vector< DiscreteTaxonData< charType > > &taxa)
void RevBayesCore::AbstractTreeHistoryCtmc<charType, treeType>::simulate( const TopologyNode &node, DiscreteCharacterData< charType > &taxa)
{
//    
//    // get the children of the node
//    const std::vector<TopologyNode*>& children = node.getChildren();
//    
//    // get the sequence of this node
//    size_t nodeIndex = node.getIndex();
//    const DiscreteTaxonData< charType > &parent = taxa[ nodeIndex ];
//    
//    // simulate the sequence for each child
//    RandomNumberGenerator* rng = GLOBAL_RNG;
//    for (std::vector< TopologyNode* >::const_iterator it = children.begin(); it != children.end(); ++it)
//    {
//        const TopologyNode &child = *(*it);
//        
//        // update the transition probability matrix
//        updateTransitionProbabilities( child.getIndex(), child.getBranchLength() );
//        
//        DiscreteTaxonData< charType > &taxon = taxa[ child.getIndex() ];
//        for ( size_t i = 0; i < numSites; ++i )
//        {
//            // get the ancestral character for this site
//            unsigned long parentState = parent.getCharacter( i ).getState();
//            size_t p = 0;
//            while ( parentState != 1 )
//            {
//                // shift to the next state
//                parentState >>= 1;
//                // increase the index
//                ++p;
//            }
//            
//            double *freqs = transitionProbMatrices[ perSiteRates[i] ][ p ];
//            
//            // create the character
//            charType c;
//            c.setToFirstState();
//            // draw the state
//            double u = rng->uniform01();
//            while ( true )
//            {
//                u -= *freqs;
//                
//                if ( u > 0.0 )
//                {
//                    ++c;
//                    ++freqs;
//                }
//                else
//                {
//                    break;
//                }
//            }
//            
//            // add the character to the sequence
//            taxon.addCharacter( c );
//        }
//        
//        if ( child.isTip() )
//        {
//            taxon.setTaxonName( child.getName() );
//        }
//        else
//        {
//            // recursively simulate the sequences
//            simulate( child, taxa, perSiteRates );
//        }
//        
//    }
    
}


template<class charType, class treeType>
void RevBayesCore::AbstractTreeHistoryCtmc<charType, treeType>::swapParameter(const DagNode *oldP, const DagNode *newP) {
    
    // we only have the topology here as the parameter
    if (oldP == tau)
    {
        tau->getValue().getTreeChangeEventHandler().removeListener( this );
        tau = static_cast<const TypedDagNode<treeType>* >( newP );
        tau->getValue().getTreeChangeEventHandler().addListener( this );
    }
    
}


template<class charType, class treeType>
void RevBayesCore::AbstractTreeHistoryCtmc<charType, treeType>::touchSpecialization( DagNode* affecter ) {
    
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
                //activeLikelihood[index] = 0;
                changedNodes[index] = true;
            }
        }
    }
    else
    {
        
        true;
        ;
    }
}

#endif /* defined(__rb_mlandis__TreeHistory__) */

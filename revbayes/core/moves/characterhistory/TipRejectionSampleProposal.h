//
//  TipRejectionSampleProposal.h
//  rb_mlandis
//
//  Created by Michael Landis on 5/7/14.
//  Copyright (c) 2014 Michael Landis. All rights reserved.
//

#ifndef __rb_mlandis__TipRejectionSampleProposal__
#define __rb_mlandis__TipRejectionSampleProposal__

#include "BranchHistory.h"
#include "DeterministicNode.h"
#include "DiscreteCharacterData.h"
#include "DistributionBinomial.h"
#include "PathRejectionSampleProposal.h"
#include "Proposal.h"
#include "RandomNumberFactory.h"
#include "RandomNumberGenerator.h"
#include "RateMap.h"
#include "RbException.h"
#include "StochasticNode.h"
//#include "TransitionProbability.h"
#include "TypedDagNode.h"

#include <cmath>
#include <iostream>
#include <set>
#include <string>

namespace RevBayesCore {
    
    /**
     * The scaling operator.
     *
     * A scaling proposal draws a random uniform number u ~ unif(-0.5,0.5)
     * and scales the current vale by a scaling factor
     * sf = exp( lambda * u )
     * where lambda is the tuning parameter of the Proposal to influence the size of the proposals.
     *
     * @copyright Copyright 2009-
     * @author The RevBayes Development Core Team (Sebastian Hoehna)
     * @since 2009-09-08, version 1.0
     *
     */
    
    template<class charType, class treeType>
    class TipRejectionSampleProposal : public Proposal {
        
    public:
        TipRejectionSampleProposal( StochasticNode<AbstractCharacterData> *n, StochasticNode<treeType>* t, DeterministicNode<RateMap> *q, double l, TopologyNode* nd=NULL );                                                                //!<  constructor
        
        // Basic utility functions
        void                            assignNode(TopologyNode* nd);
        void                            assignSiteIndexSet(const std::set<size_t>& s);
        TipRejectionSampleProposal*     clone(void) const;                                                                  //!< Clone object
        void                            cleanProposal(void);
        double                          doProposal(void);                                                                   //!< Perform proposal
        const std::vector<DagNode*>&    getNodes(void) const;                                                               //!< Get the vector of DAG nodes this proposal is working on
        const std::string&              getProposalName(void) const;                                                        //!< Get the name of the proposal for summary printing
        void                            printParameterSummary(std::ostream &o) const;                                       //!< Print the parameter summary
        void                            prepareProposal(void);                                                              //!< Prepare the proposal
        double                          sampleTipCharacters(const std::set<size_t>& indexSet);    //!< Sample the characters at the node
        void                            swapNode(DagNode *oldN, DagNode *newN);                                             //!< Swap the DAG nodes on which the Proposal is working on
        void                            tune(double r);                                                                     //!< Tune the proposal to achieve a better acceptance/rejection ratio
        void                            undoProposal(void);                                                                 //!< Reject the proposal
        
    protected:
        
        // parameters
        StochasticNode<AbstractCharacterData>*  ctmc;
        StochasticNode<treeType>*               tau;
        DeterministicNode<RateMap>*             qmap;
        std::vector<DagNode*>                   nodes;
        
        // dimensions
        size_t                                  numNodes;
        size_t                                  numCharacters;
        size_t                                  numStates;
                
        // proposal
        std::vector<unsigned>                   storedNodeState;
        TopologyNode*                           node;
        std::set<size_t>                        siteIndexSet;
        double                                  storedLnProb;
        double                                  proposedLnProb;
        
        PathRejectionSampleProposal<charType,treeType>* nodeProposal;
        TransitionProbabilityMatrix nodeTpMatrix;
        
        double                                  lambda;
        int                                     monitorIndex;
        
        // flags
        bool                                    fixNodeIndex;
        bool                                    sampleNodeIndex;
        bool                                    sampleSiteIndexSet;
        
    };
    
}



/**
 * Constructor
 *
 * Here we simply allocate and initialize the Proposal object.
 */
template<class charType, class treeType>
RevBayesCore::TipRejectionSampleProposal<charType, treeType>::TipRejectionSampleProposal( StochasticNode<AbstractCharacterData> *n, StochasticNode<treeType> *t, DeterministicNode<RateMap>* q, double l, TopologyNode* nd) : Proposal(),
ctmc(n),
tau(t),
qmap(q),
numNodes(t->getValue().getNumberOfNodes()),
numCharacters(n->getValue().getNumberOfCharacters()),
numStates(static_cast<const DiscreteCharacterState&>(n->getValue().getCharacter(0,0)).getNumberOfStates()),
node(nd),
nodeTpMatrix(numStates),
lambda(l),
sampleNodeIndex(true),
sampleSiteIndexSet(true)

{
    monitorIndex = -120;
    
    nodes.push_back(ctmc);
    nodes.push_back(tau);
    nodes.push_back(qmap);
    
    nodeProposal = new PathRejectionSampleProposal<charType,treeType>(n,t,q,l,nd);
    
    fixNodeIndex = (node != NULL);
}


template<class charType, class treeType>
void RevBayesCore::TipRejectionSampleProposal<charType, treeType>::cleanProposal( void )
{
//    std::cout << "ACCEPT " << node->getIndex() << "\n";
    nodeProposal->cleanProposal();
}

/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'B'.
 *
 * \return A new copy of the proposal.
 */
template<class charType, class treeType>
RevBayesCore::TipRejectionSampleProposal<charType, treeType>* RevBayesCore::TipRejectionSampleProposal<charType, treeType>::clone( void ) const
{
    return new TipRejectionSampleProposal( *this );
}

template<class charType, class treeType>
void RevBayesCore::TipRejectionSampleProposal<charType, treeType>::assignNode(TopologyNode* nd)
{
    node = nd;
    sampleNodeIndex = false;
}

template<class charType, class treeType>
void RevBayesCore::TipRejectionSampleProposal<charType, treeType>::assignSiteIndexSet(const std::set<size_t>& s)
{
    siteIndexSet = s;
    sampleSiteIndexSet = false;
}


/**
 * Get Proposals' name of object
 *
 * \return The Proposals' name.
 */
template<class charType, class treeType>
const std::string& RevBayesCore::TipRejectionSampleProposal<charType, treeType>::getProposalName( void ) const
{
    static std::string name = "TipRejectionSampleProposal";
    
    return name;
}


/**
 * Get the vector of nodes on which this proposal is working on.
 *
 * \return  Const reference to a vector of nodes pointer on which the proposal operates.
 */
template<class charType, class treeType>
const std::vector<RevBayesCore::DagNode*>& RevBayesCore::TipRejectionSampleProposal<charType, treeType>::getNodes( void ) const
{
    
    return nodes;
}


/**
 * Perform the Proposal.
 *
 * A scaling Proposal draws a random uniform number u ~ unif(-0.5,0.5)
 * and scales the current vale by a scaling factor
 * sf = exp( lambda * u )
 * where lambda is the tuning parameter of the Proposal to influence the size of the proposals.
 *
 * \return The hastings ratio.
 */
template<class charType, class treeType>
double RevBayesCore::TipRejectionSampleProposal<charType, treeType>::doProposal( void )
{
    proposedLnProb = 0.0;

    double proposedLnProbRatio = 0.0;
    
    AbstractTreeHistoryCtmc<charType, treeType>* p = static_cast< AbstractTreeHistoryCtmc<charType, treeType>* >(&ctmc->getDistribution());
//    p->getHistory(*node).print();
    
    // update 1x pathEnd and 1x pathHistory values
    proposedLnProbRatio += sampleTipCharacters(siteIndexSet);
    proposedLnProbRatio += nodeProposal->doProposal();

//    p->getHistory(*node).print();
    return proposedLnProbRatio;
}


/**
 *
 */
template<class charType, class treeType>
void RevBayesCore::TipRejectionSampleProposal<charType, treeType>::prepareProposal( void )
{
    
    storedLnProb = 0.0;
    
    size_t numTips = tau->getValue().getNumberOfTips();
    if (sampleNodeIndex && !fixNodeIndex)
    {
        size_t idx = GLOBAL_RNG->uniform01() * numTips;
        node = &tau->getValue().getNode(idx);
    }
    
    if (sampleSiteIndexSet)
    {
        siteIndexSet.clear();
        siteIndexSet.insert(GLOBAL_RNG->uniform01() * numCharacters); // at least one is inserted
        for (size_t i = 0; i < numCharacters; i++)
        {
            if (GLOBAL_RNG->uniform01() < lambda)
            {
                siteIndexSet.insert(i);
            }
        }
//        std::cout << "sites ";
//        for (std::set<size_t>::iterator it = siteIndexSet.begin(); it != siteIndexSet.end(); it++)
//        {
//            std::cout << *it << " ";
//        }
//        std::cout << "\n";
    }
    
    AbstractTreeHistoryCtmc<charType, treeType>* p = static_cast< AbstractTreeHistoryCtmc<charType, treeType>* >(&ctmc->getDistribution());
    
    nodeProposal->assignNode(node);
    nodeProposal->assignSiteIndexSet(siteIndexSet);
    nodeProposal->prepareProposal();
    
    // store node state values
    storedNodeState.clear();
    
    storedNodeState.resize(numCharacters,0);
    const std::vector<CharacterEvent*>& nodeState = p->getHistory(*node).getChildCharacters();
    for (std::set<size_t>::iterator it = siteIndexSet.begin(); it != siteIndexSet.end(); it++)
    {
        unsigned s = nodeState[*it]->getState();
        storedNodeState[*it] = s;
    }

    sampleNodeIndex = true;
    sampleSiteIndexSet = true;
    
}


/**
 * Print the summary of the Proposal.
 *
 * The summary just contains the current value of the tuning parameter.
 * It is printed to the stream that it passed in.
 *
 * \param[in]     o     The stream to which we print the summary.
 */
template<class charType, class treeType>
void RevBayesCore::TipRejectionSampleProposal<charType, treeType>::printParameterSummary(std::ostream &o) const
{
    o << "lambda = " << lambda;
}


template<class charType, class treeType>
double RevBayesCore::TipRejectionSampleProposal<charType, treeType>::sampleTipCharacters(const std::set<size_t>& indexSet)
{
    double lnP = 0.0;

    
    if (node->isTip())
    {
        
        AbstractTreeHistoryCtmc<charType, treeType>* p = static_cast< AbstractTreeHistoryCtmc<charType, treeType>* >(&ctmc->getDistribution());
        
        qmap->getValue().calculateTransitionProbabilities(*node, nodeTpMatrix);
        
        std::vector<BranchHistory*> histories = p->getHistories();
        
        // for sampling probs
        const std::vector<CharacterEvent*>& nodeParentState = histories[node->getIndex()]->getParentCharacters();
        const std::vector<double>& tipProbs = p->getTipProbs(node->getIndex());

        // to update
        std::vector<CharacterEvent*> nodeChildState = histories[node->getIndex()]->getChildCharacters();
                
        for (std::set<size_t>::iterator it = siteIndexSet.begin(); it != siteIndexSet.end(); it++)
        {
            unsigned int ancS = nodeParentState[*it]->getState();
            
            double u = GLOBAL_RNG->uniform01();
            double g0 = nodeTpMatrix[ancS][0] * (1.0 - tipProbs[*it]);
            double g1 = nodeTpMatrix[ancS][1] * tipProbs[*it];
            
            unsigned int s = 0;
            if (u < g1 / (g0 + g1))
                s = 1;
            
            nodeChildState[*it]->setState(s);
        }
    }
    
    return lnP;
}

/**
 * Reject the Proposal.
 *
 * Since the Proposal stores the previous value and it is the only place
 * where complex undo operations are known/implement, we need to revert
 * the value of the ctmc/DAG-node to its original value.
 */
template<class charType, class treeType>
void RevBayesCore::TipRejectionSampleProposal<charType, treeType>::undoProposal( void )
{
    // swap current value and stored value
    AbstractTreeHistoryCtmc<charType, treeType>* p = static_cast< AbstractTreeHistoryCtmc<charType, treeType>* >(&ctmc->getDistribution());
    const std::vector<BranchHistory*>& histories = p->getHistories();
    
    // restore path state
//    std::cout << "REJECT " << node->getIndex() << "\n";
    nodeProposal->undoProposal();
    
    BranchHistory* bh = &p->getHistory(*node);
//    bh->print();

    // restore node state
    std::vector<CharacterEvent*> nodeChildState = histories[node->getIndex()]->getChildCharacters();
    
    for (std::set<size_t>::iterator it = siteIndexSet.begin(); it != siteIndexSet.end(); it++)
    {
        unsigned s = storedNodeState[*it];
        nodeChildState[*it]->setState(s);
    }
}


/**
 * Swap the current ctmc for a new one.
 *
 * \param[in]     oldN     The old ctmc that needs to be replaced.
 * \param[in]     newN     The new ctmc.
 */
template<class charType, class treeType>
void RevBayesCore::TipRejectionSampleProposal<charType, treeType>::swapNode(DagNode *oldN, DagNode *newN)
{
    
    if (oldN == ctmc)
    {
        ctmc = static_cast<StochasticNode<AbstractCharacterData>* >(newN) ;
    }
    else if (oldN == tau)
    {
        tau = static_cast<StochasticNode<treeType>* >(newN);
    }
    else if (oldN == qmap)
    {
        qmap = static_cast<DeterministicNode<RateMap>* >(newN);
    }
    
    nodeProposal->swapNode(oldN, newN);

}


/**
 * Tune the Proposal to accept the desired acceptance ratio.
 */
template<class charType, class treeType>
void RevBayesCore::TipRejectionSampleProposal<charType, treeType>::tune( double rate )
{
    ; // do nothing
}

#endif /* defined(__rb_mlandis__TipRejectionSampleProposal__) */

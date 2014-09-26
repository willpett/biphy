//
//  PathRejectionSampleProposal.h
//  rb_mlandis
//
//  Created by Michael Landis on 5/7/14.
//  Copyright (c) 2014 Michael Landis. All rights reserved.
//

#ifndef __rb_mlandis__PathRejectionSampleProposal__
#define __rb_mlandis__PathRejectionSampleProposal__

#include "BranchHistory.h"
#include "DeterministicNode.h"
#include "DiscreteCharacterData.h"
#include "DistributionBinomial.h"
#include "Proposal.h"
#include "RandomNumberFactory.h"
#include "RandomNumberGenerator.h"
#include "RateMap.h"
#include "RbException.h"
#include "StochasticNode.h"
#include "TopologyNode.h"
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
    class PathRejectionSampleProposal : public Proposal {
        
    public:
        PathRejectionSampleProposal( StochasticNode<AbstractCharacterData> *n, StochasticNode<treeType>* t, DeterministicNode<RateMap> *q, double l, TopologyNode* nd=NULL);   //!<  constructor
        PathRejectionSampleProposal( const PathRejectionSampleProposal& p );
        
        // Basic utility functions
        void                            assignNode(TopologyNode* nd);
        void                            assignSiteIndexSet(const std::set<size_t>& s);
        double                          computeLnProposal(const TopologyNode& nd, const BranchHistory& bh);
        void                            cleanProposal(void);
        PathRejectionSampleProposal*    clone(void) const;                                                                  //!< Clone object
        double                          doProposal(void);                                                                   //!< Perform proposal
        const std::vector<DagNode*>&    getNodes(void) const;                                                               //!< Get the vector of DAG nodes this proposal is working on
        const std::string&              getProposalName(void) const;                                                        //!< Get the name of the proposal for summary printing
        void                            printParameterSummary(std::ostream &o) const;                                       //!< Print the parameter summary
        void                            prepareProposal(void);                                                              //!< Prepare the proposal
        void                            swapNode(DagNode *oldN, DagNode *newN);                                             //!< Swap the DAG nodes on which the Proposal is working on
        void                            tune(double r);                                                                     //!< Tune the proposal to achieve a better acceptance/rejection ratio
        void                            undoProposal(void);                                                                 //!< Reject the proposal
        
    protected:
        
        void                            fillStateCounts(std::vector<CharacterEvent*> s, unsigned* counts);
        
        // parameters
        StochasticNode<AbstractCharacterData>*  ctmc;
        StochasticNode<treeType>*               tau;
        DeterministicNode<RateMap>*             qmap;
        
        //BranchHistory*                          storedValue;
        std::multiset<CharacterEvent*,CharacterEventCompare> storedHistory;
        std::multiset<CharacterEvent*,CharacterEventCompare> proposedHistory;
        
        TopologyNode*                           node;
        std::set<size_t>                        siteIndexSet;
        
        double                                  storedLnProb;
        double                                  proposedLnProb;
        
        size_t                                  numNodes;
        size_t                                  numCharacters;
        size_t                                  numStates;
        
        double                                  lambda;
        bool                                    fixNodeIndex;
        bool                                    sampleNodeIndex;
        bool                                    sampleSiteIndexSet;
        
        bool                                    printDebug;
        
        std::vector<DagNode*>                   nodes;
    };
    
}



/**
 * Constructor
 *
 * Here we simply allocate and initialize the Proposal object.
 */
template<class charType, class treeType>
RevBayesCore::PathRejectionSampleProposal<charType, treeType>::PathRejectionSampleProposal( StochasticNode<AbstractCharacterData> *n, StochasticNode<treeType> *t, DeterministicNode<RateMap>* q, double l, TopologyNode* nd) : Proposal(),
    ctmc(n),
    tau(t),
    qmap(q),
    node(nd),
    lambda(l),
    sampleNodeIndex(true),
    sampleSiteIndexSet(true)
{
    nodes.push_back(ctmc);
    nodes.push_back(tau);
    nodes.push_back(qmap);
    
    numNodes = t->getValue().getNumberOfNodes();
    numCharacters = n->getValue().getNumberOfCharacters();
    numStates = static_cast<const DiscreteCharacterState&>(n->getValue().getCharacter(0,0)).getNumberOfStates();
 
    printDebug = false;
    
    fixNodeIndex = (node != NULL);
}

template<class charType, class treeType>
RevBayesCore::PathRejectionSampleProposal<charType, treeType>::PathRejectionSampleProposal(const PathRejectionSampleProposal& m) : Proposal()
{
    if (this != &m)
    {
        ctmc = m.ctmc;
        tau = m.tau;
        qmap = m.qmap;
        node = m.node;
        lambda = m.lambda;
        sampleNodeIndex = m.sampleNodeIndex;
        sampleSiteIndexSet  = m.sampleSiteIndexSet;
        
        numNodes =  m.numNodes;
        numCharacters = m.numCharacters;
        numStates = m.numStates;
            
        printDebug = m.printDebug;
        fixNodeIndex = m.fixNodeIndex;
    }
}

template<class charType, class treeType>
void RevBayesCore::PathRejectionSampleProposal<charType, treeType>::cleanProposal( void )
{
    

    AbstractTreeHistoryCtmc<charType,treeType>& p = static_cast< AbstractTreeHistoryCtmc<charType, treeType>& >(ctmc->getDistribution());
    BranchHistory* bh = &p.getHistory(*node);

    if (printDebug) { std::cout << node->getIndex() << " ACCEPT\n"; bh->print(); std::cout << "------\n"; }
    
    // delete old events
    std::multiset<CharacterEvent*,CharacterEventCompare>::iterator it_h;
    for (it_h = storedHistory.begin(); it_h != storedHistory.end(); it_h++)
    {
//         std::cout << "clean deleting" <<  *it_h << "\n";
        delete *it_h;
    }
    storedHistory.clear();
    proposedHistory.clear();
}

/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'B'.
 *
 * \return A new copy of the proposal.
 */
template<class charType, class treeType>
RevBayesCore::PathRejectionSampleProposal<charType, treeType>* RevBayesCore::PathRejectionSampleProposal<charType, treeType>::clone( void ) const
{
    return new PathRejectionSampleProposal( *this );
}

template<class charType, class treeType>
void RevBayesCore::PathRejectionSampleProposal<charType, treeType>::assignNode(TopologyNode* nd)
{
    node = nd;
    sampleNodeIndex = false;
}

template<class charType, class treeType>
void RevBayesCore::PathRejectionSampleProposal<charType, treeType>::assignSiteIndexSet(const std::set<size_t>& s)
{
    siteIndexSet = s;
    sampleSiteIndexSet = false;
}


template<class charType, class treeType>
double RevBayesCore::PathRejectionSampleProposal<charType, treeType>::computeLnProposal(const TopologyNode& nd, const BranchHistory& bh)
{
    double lnP = 0.0;

    std::vector<CharacterEvent*> currState = bh.getParentCharacters();
    const std::multiset<CharacterEvent*,CharacterEventCompare>& history = bh.getHistory();
    std::multiset<CharacterEvent*,CharacterEventCompare>::iterator it_h;

    unsigned counts[numStates];
    for (size_t i = 0; i < numStates; i++)
        counts[i] = 0;
    fillStateCounts(currState, counts);

    const treeType& tree = tau->getValue();
    //const TopologyNode& node = tree.getNode(bh.getIndex());
    double branchLength = nd.getBranchLength();
    if (nd.isRoot())
    {
        return 0.0;
        branchLength = 1000.0;
    }
    
    
    double currAge = (!nd.isRoot() ? nd.getParent().getAge() : 10e200);
    
    // get sampling ratemap
    const RateMap& rm = qmap->getValue();
    
    // stepwise events
    double t = 0.0;
    double dt = 0.0;
    for (it_h = history.begin(); it_h != history.end(); it_h++)
    {
        // next event time
        double idx = (*it_h)->getIndex();
        dt = (*it_h)->getTime() - t;
    
        double tr = rm.getRate(nd, currState, *it_h, counts, currAge);
        double sr = rm.getSumOfRates(nd, currState, counts, currAge);
        
        // lnP for stepwise events for p(x->y)
        lnP += log(tr) - sr * dt * branchLength;
        
        // update counts
        counts[ currState[idx]->getState() ] -= 1;
        counts[ (*it_h)->getState() ] += 1;
        
        // update state
        currState[idx] = *it_h;
        t += dt;
        currAge -= dt * branchLength;
    }
    // lnL for final non-event
    double sr = rm.getSumOfRates(nd, currState, counts, currAge);
    lnP += -sr * (1.0 - t) * branchLength;
    
    if (nd.isRoot())
        return 0.0;
    
    return lnP;
}

template<class charType, class treeType>
void RevBayesCore::PathRejectionSampleProposal<charType, treeType>::fillStateCounts(std::vector<CharacterEvent*> s, unsigned int *counts)
{
    for (size_t i = 0; i < s.size(); i++)
        counts[ s[i]->getState() ] += 1;
}


/**
 * Get Proposals' name of object
 *
 * \return The Proposals' name.
 */
template<class charType, class treeType>
const std::string& RevBayesCore::PathRejectionSampleProposal<charType, treeType>::getProposalName( void ) const
{
    static std::string name = "PathRejectionSampleProposal";
    
    return name;
}


/**
 * Get the vector of nodes on which this proposal is working on.
 *
 * \return  Const reference to a vector of nodes pointer on which the proposal operates.
 */
template<class charType, class treeType>
const std::vector<RevBayesCore::DagNode*>& RevBayesCore::PathRejectionSampleProposal<charType, treeType>::getNodes( void ) const
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
double RevBayesCore::PathRejectionSampleProposal<charType, treeType>::doProposal( void )
{

    //if (nodeIndex == 5) { std::cout << "path old " << nodeIndex << "\n";    storedValue->print(); }
    
    AbstractTreeHistoryCtmc<charType,treeType>& p = static_cast< AbstractTreeHistoryCtmc<charType, treeType>& >(ctmc->getDistribution());
    proposedHistory.clear();
    
    // get model parameters
    const treeType& tree = this->tau->getValue();
    double branchLength = node->getBranchLength(); //tree.getBranchLength(node.getIndex());
    if (node->isRoot())
    {
        return 0.0;
        branchLength = 1000.0;//2*tree.getTreeLength();
    }
    
    const RateMap& rm = qmap->getValue();

    // clear characters
    BranchHistory* bh = &p.getHistory(*node);
    if (printDebug) { std::cout << "BEFORE PATH PROPOSAL\n"; bh->print(); }
    
    //proposedValue->clearEvents( siteIndexSet );
    
    // reject sample path history
    std::vector<CharacterEvent*> parentVector = bh->getParentCharacters();
    std::vector<CharacterEvent*> childVector =  bh->getChildCharacters();
    
//    std::cout << "path  " << node->getIndex() << " ";
//    for (size_t i = 0; i < numCharacters; i++)
//    {
//        std::cout << parentVector[i]->getState();
//    }
//    std::cout << "\n";
    //bh->print();
    
    for (std::set<size_t>::iterator it = siteIndexSet.begin(); it != siteIndexSet.end(); it++)
    {
        std::set<CharacterEvent*> tmpHistory;
        unsigned ct = 0;
        bool reject = false;
        
        unsigned int currState = parentVector[*it]->getState();
        unsigned int endState = childVector[*it]->getState();
        do
        {
//            std::cout << node->getIndex() << " " << *it << " " << ct++ << "\n";
            // delete previously rejected events
            tmpHistory.clear();
                      
            // proceed with rejection sampling
            currState = parentVector[*it]->getState();
            double t = 0.0;
            do
            {
                unsigned int nextState = (currState == 1 ? 0 : 1);
                double r = rm.getSiteRate(*node, currState, nextState);
                
                t += RbStatistics::Exponential::rv(r * branchLength, *GLOBAL_RNG);
                if (t < 1.0)
                {
                    currState = nextState;
                    CharacterEvent* evt = new CharacterEvent(*it, nextState, t);
                    tmpHistory.insert(evt);
                }
                else if (currState != endState)
                {
                    for (std::set<CharacterEvent*>::iterator it_h = tmpHistory.begin(); it_h != tmpHistory.end(); it_h++)
                        delete *it_h;
                }
                
//                std::cout << node->getIndex() << " " << *it << " " << t << " " << currState << " " << nextState << " " << endState << "\n";
            }
            while(t < 1.0);
//            std::cout << "tmpHistory.size() " << tmpHistory.size() << "\n";
        }
        while (currState != endState);
        
        for (std::set<CharacterEvent*>::iterator it = tmpHistory.begin(); it != tmpHistory.end(); it++)
        {
            proposedHistory.insert(*it);
//            history.insert(*it);
        }
    }
    
    // assign values back to model for likelihood
    bh->updateHistory(proposedHistory, siteIndexSet);
    //bh->print();

//    double r1=rm.getSiteRate(*node, 0, 1);
//    double r0=rm.getSiteRate(*node, 1,0);
//    double sum_r_ns = (r0+r1)/2 * numCharacters;
//    std::cout << node->getIndex() << " " << branchLength << " " << r0 << " " << r1 << " " << sum_r_ns << " " << branchLength*sum_r_ns << " " << bh->getNumEvents() << "\n";
    
    if (printDebug) { std::cout << "AFTER PATH PROPOSAL\n"; bh->print(); std::cout << "\n"; }
    
    // return hastings ratio
    proposedLnProb = computeLnProposal(*node, *bh);
    
    return storedLnProb - proposedLnProb;
}


/**
 *
 */
template<class charType, class treeType>
void RevBayesCore::PathRejectionSampleProposal<charType, treeType>::prepareProposal( void )
{
    storedHistory.clear();
    proposedHistory.clear();
    
    storedLnProb = 0.0;
    AbstractTreeHistoryCtmc<charType,treeType>& p = static_cast< AbstractTreeHistoryCtmc<charType, treeType>& >(ctmc->getDistribution());
    
    if (sampleNodeIndex && !fixNodeIndex)
    {
        size_t nodeIndex = GLOBAL_RNG->uniform01() * numNodes;
        node = &tau->getValue().getNode(nodeIndex);
    }
    if (node->isRoot()) return;

    // Enable to use dirty lnL flagging.
    p.fireTreeChangeEvent(*node);
    
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
    }

    sampleNodeIndex = true;
    sampleSiteIndexSet = true;

    BranchHistory* bh = &p.getHistory(*node);
    const std::multiset<CharacterEvent*,CharacterEventCompare>& history = bh->getHistory();

    // store history for events in siteIndexSet
    std::multiset<CharacterEvent*,CharacterEventCompare>::iterator it_h;
    for (it_h = history.begin(); it_h != history.end(); it_h++)
    {
        if (siteIndexSet.find( (*it_h)->getIndex() ) != siteIndexSet.end())
        {
//            if (node->getIndex() == 41)
//                std::cout << "yow " << (*it_h)->getTime() << " " << (*it_h)->getIndex() << " " << (*it_h)->getState() << "\n";

            storedHistory.insert(*it_h);
            //storedHistory.insert(new CharacterEvent(**it_h));
        }
    }
    
    storedLnProb = computeLnProposal(*node, *bh);

    if (printDebug) { std::cout << node->getIndex() << " PREPARE proposal\n"; bh->print(); }
    
    //    storedValue = new BranchHistory( p.getHistory(nodeIndex) );
    //    std::cout << "--------------------\n"; std::cout << "before proposal\n";    storedValue->print();    std::cout << "--------------------\n";
 //   storedLnProb = computeLnProposal(*storedValue);
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
void RevBayesCore::PathRejectionSampleProposal<charType, treeType>::printParameterSummary(std::ostream &o) const
{
    o << "lambda = " << lambda;
}


/**
 * Reject the Proposal.
 *
 * Since the Proposal stores the previous value and it is the only place
 * where complex undo operations are known/implement, we need to revert
 * the value of the ctmc/DAG-node to its original value.
 */
template<class charType, class treeType>
void RevBayesCore::PathRejectionSampleProposal<charType, treeType>::undoProposal( void )
{
    // swap current value and stored value
    AbstractTreeHistoryCtmc<charType, treeType>* p = static_cast< AbstractTreeHistoryCtmc<charType, treeType>* >(&ctmc->getDistribution());
    BranchHistory* bh = &p->getHistory(*node);
 
    if (printDebug) { std::cout << node->getIndex() << " BEFORE REJECT\n"; bh->print(); }
    
    bh->updateHistory(storedHistory, siteIndexSet);
    
    // delete new events
    std::multiset<CharacterEvent*,CharacterEventCompare>::iterator it_h;
    for (it_h = proposedHistory.begin(); it_h != proposedHistory.end(); it_h++)
    {
        delete *it_h;
    }
    storedHistory.clear();
    proposedHistory.clear();
    
    if (printDebug) { std::cout << node->getIndex() << " AFTER REJECT\n"; bh->print(); std::cout << "------\n"; }
    
}


/**
 * Swap the current ctmc for a new one.
 *
 * \param[in]     oldN     The old ctmc that needs to be replaced.
 * \param[in]     newN     The new ctmc.
 */
template<class charType, class treeType>
void RevBayesCore::PathRejectionSampleProposal<charType, treeType>::swapNode(DagNode *oldN, DagNode *newN)
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
    
}


/**
 * Tune the Proposal to accept the desired acceptance ratio.
 */
template<class charType, class treeType>
void RevBayesCore::PathRejectionSampleProposal<charType, treeType>::tune( double rate )
{
    ; // do nothing
}

#endif /* defined(__rb_mlandis__PathRejectionSampleProposal__) */

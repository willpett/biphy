//
//  TreeCharacterHistoryNhxMonitor.h
//  rb_mlandis
//
//  Created by Michael Landis on 5/6/14.
//  Copyright (c) 2014 Michael Landis. All rights reserved.
//

#ifndef rb_mlandis_TreeCharacterHistoryNhxMonitor_h
#define rb_mlandis_TreeCharacterHistoryNhxMonitor_h


#include "Monitor.h"
#include "BranchHistory.h"
#include "StochasticNode.h"
#include "TypedDagNode.h"
#include "TimeTree.h"

#include <fstream>
#include <iostream>
#include <string>
#include <vector>

namespace RevBayesCore {
    
    template<class charType, class treeType>
    class TreeCharacterHistoryNhxMonitor : public Monitor {
        
    public:
        // Constructors and Destructors
        TreeCharacterHistoryNhxMonitor(StochasticNode<AbstractCharacterData>* s, TypedDagNode<treeType> *t, std::vector<std::vector<double> > gc, int g, int mg, int burn, const std::string &fname, const std::string &del, bool pp=true, bool l=true, bool pr=true, bool ap=false, bool sm=true, bool sr=true);
        
        TreeCharacterHistoryNhxMonitor(const TreeCharacterHistoryNhxMonitor& f);
        
        // basic methods
        TreeCharacterHistoryNhxMonitor*   clone(void) const;                                                  //!< Clone the object
        
        // Monitor functions
        void                                    monitor(long gen);                                                  //!< Monitor at generation gen
        void                                    swapNode(DagNode *oldN, DagNode *newN);
        
        // FileMonitor functions
        void                                    closeStream(void);                                                  //!< Close stream after finish writing
        void                                    openStream(void);                                                   //!< Open the stream for writing
        void                                    printHeader(void);                                                  //!< Print header
        std::vector<unsigned int>               getChildCharacterCounts(int idx);
        std::vector<unsigned int>               getParentCharacterCounts(int idx);
        long                                    getNumSamples(void);
        
    private:
        std::string                             buildExtendedNewick();
        std::string                             buildExtendedNewick(TopologyNode* n);
        std::string                             buildCharacterHistoryString(TopologyNode* n, std::string brEnd="child");
        void                                    updateCharacterCounts(TopologyNode* n, std::string brEnd="child");
        std::string                             buildNhxString(void);
        
        // the stream to print
        std::fstream                            outStream;
        
        // parameters
        StochasticNode<AbstractCharacterData>*    variable;
        TypedDagNode<treeType>*                 tree;
        std::set<DagNode *>                     nodeVariables;
        std::vector<std::vector<double> >       geographicCoordinates;
        
        std::vector<std::vector<unsigned int> > parentCharacterCounts;
        std::vector<std::vector<unsigned int> > childCharacterCounts;
        
        size_t numHistories;
        size_t numCharacters;
        
        std::string                         filename;
        std::string                         separator;
        bool                                posterior;
        bool                                prior;
        bool                                likelihood;
        bool                                append;
        bool                                showMetadata;
        bool                                showRates;
        long                                numSamples;
        long                                maxGen;
        long                                burn;
        
    };
    
}

#include "AbstractTreeHistoryCtmc.h"
#include "Model.h"
#include "RbException.h"

#include <sstream>

/* Constructor */
template<class charType, class treeType>
RevBayesCore::TreeCharacterHistoryNhxMonitor<charType, treeType>::TreeCharacterHistoryNhxMonitor(StochasticNode<AbstractCharacterData>* s, TypedDagNode<treeType>* t, std::vector<std::vector<double> > gc, int g, int mg, int b, const std::string &fname, const std::string &del, bool pp, bool l, bool pr, bool ap, bool sm, bool sr) :
Monitor(g,t),
outStream(),
variable( s ),
tree( t ),
geographicCoordinates(gc),
filename( fname ),
separator( del ),
posterior( pp ),
prior( pr ),
likelihood( l ),
append(ap),
showMetadata(sm),
showRates(sr),
numSamples(0),
maxGen(mg),
burn(b) {
    
    std::cout << g << " " << mg << "\n";
    
    nodes.push_back(s);
    nodes.push_back(t);
    //for (size_t i = 0; i < branchHistories.size(); i++)
    //    nodes.push_back(branchHistories[i]);
    
//    numHistories = bh.size();
//    numCharacters = bh[0]->getValue().getNumCharacters();
    
    parentCharacterCounts.resize(numHistories);
    childCharacterCounts.resize(numHistories);
    for (size_t i = 0; i < numHistories; i++)
    {
        parentCharacterCounts[i].resize(numCharacters,0);
        childCharacterCounts[i].resize(numCharacters,0);
    }
}

template<class charType, class treeType>
RevBayesCore::TreeCharacterHistoryNhxMonitor<charType, treeType>::TreeCharacterHistoryNhxMonitor(const TreeCharacterHistoryNhxMonitor &m) :
Monitor( m ),
outStream( ),
variable( m.variable ),
tree( m.tree ),
nodeVariables( m.nodeVariables ),
geographicCoordinates(m.geographicCoordinates),
parentCharacterCounts(m.parentCharacterCounts),
childCharacterCounts(m.childCharacterCounts),
numHistories(m.numHistories),
numCharacters(m.numCharacters),
showMetadata(m.showMetadata),
showRates(m.showRates),
numSamples(m.numSamples),
maxGen(m.maxGen),
burn(m.burn) {
    
    filename    = m.filename;
    separator   = m.separator;
    prior       = m.prior;
    posterior   = m.posterior;
    likelihood  = m.likelihood;
    append      = m.append;
}


/* Clone the object */
template<class charType, class treeType>
RevBayesCore::TreeCharacterHistoryNhxMonitor<charType, treeType>* RevBayesCore::TreeCharacterHistoryNhxMonitor<charType, treeType>::clone(void) const {
    
    return new TreeCharacterHistoryNhxMonitor<charType,treeType>(*this);
}

template<class charType, class treeType>
void RevBayesCore::TreeCharacterHistoryNhxMonitor<charType, treeType>::closeStream() {
    outStream.close();
}

template<class charType, class treeType>
std::string RevBayesCore::TreeCharacterHistoryNhxMonitor<charType, treeType>::buildExtendedNewick( void ) {
    //tree->getValue().getRoot().setNewickNeedsRefreshing(true);
    numSamples++;
    std::string newick = buildExtendedNewick( &tree->getValue().getRoot() );
    return newick;
}

template<class charType, class treeType>
void RevBayesCore::TreeCharacterHistoryNhxMonitor<charType, treeType>::updateCharacterCounts(TopologyNode* n, std::string brEnd)
{
    
    AbstractTreeHistoryCtmc<charType, treeType>* p = static_cast< AbstractTreeHistoryCtmc<charType, treeType>* >(&variable->getDistribution());
    const BranchHistory& bh = p->getHistory(*n);

    std::vector<CharacterEvent*> characters;
    if (brEnd=="child")
    {
        characters = bh.getChildCharacters();
        for (size_t i = 0; i < characters.size(); i++)
            childCharacterCounts[n->getIndex()][i] += characters[i]->getState();
    }
    else if (brEnd=="parent")
    {
        characters = bh.getParentCharacters();
        for (size_t i = 0; i < characters.size(); i++)
            parentCharacterCounts[n->getIndex()][i] += characters[i]->getState();
    }
    
}

template<class charType, class treeType>
std::string RevBayesCore::TreeCharacterHistoryNhxMonitor<charType, treeType>::buildCharacterHistoryString(TopologyNode* n, std::string brEnd)
{
    size_t nd_idx = n->getIndex();
    AbstractTreeHistoryCtmc<charType, treeType>* p = static_cast< AbstractTreeHistoryCtmc<charType, treeType>* >(&variable->getDistribution());
    const BranchHistory& bh = p->getHistory(*n);

    std::vector<CharacterEvent*> characters;
    std::stringstream ss;
    
    if (brEnd=="child")
    {
        characters = bh.getChildCharacters();
        for (size_t i = 0; i < characters.size(); i++)
        {
            if (i != 0)
                ss << ",";
            ss << (double)childCharacterCounts[nd_idx][i] / numSamples;
        }
    }
    else if (brEnd=="parent")
    {
        characters = bh.getParentCharacters();
        for (size_t i = 0; i < characters.size(); i++)
        {
            if (i != 0)
                ss << ",";
            ss << (double)parentCharacterCounts[nd_idx][i] / numSamples;
        }
    }
    
    return ss.str();
}

/* Build newick string */
template<class charType, class treeType>
std::string RevBayesCore::TreeCharacterHistoryNhxMonitor<charType, treeType>::buildExtendedNewick( TopologyNode* n ) {
    // create the newick string
    std::stringstream o;
    
    // extended data is only found on admixture nodes
    std::string additionalInfo = "";
    
    // loop over admixture nodes per branch
    std::stringstream characterStream;
    
    double br = 1.0;
    
    if (showMetadata)
    {
        updateCharacterCounts(n,"child");
        
        characterStream << "[";
        
        // character history
        characterStream << "&ch={" << buildCharacterHistoryString(n,"child") << "}";
        
        // ... whatever else
        characterStream << "]";
        
        additionalInfo = characterStream.str();
    }
    
    // test whether this is a internal or external node
    if (n->isTip()) {
        // this is a tip so we just return the name of the node
        o << n->getIndex() << additionalInfo << ":" << n->getBranchLength();
    }
    
    else {
        o << "(";
        for (size_t i=0; i<(n->getNumberOfChildren()-1); i++) {
            o << buildExtendedNewick( &n->getChild(i) ) << ",";
        }
        o << buildExtendedNewick( &n->getChild(n->getNumberOfChildren()-1) ) << ")" << additionalInfo << ":" << n->getBranchLength() * br;
    }
    
    return o.str();
}


/** Monitor value at generation gen */
template<class charType, class treeType>
void RevBayesCore::TreeCharacterHistoryNhxMonitor<charType, treeType>::monitor(long gen) {
    
    // get the printing frequency
    int samplingFrequency = printgen;
    
    if (gen % samplingFrequency == 0 && gen != maxGen) { // && gen >= burn) {
        
        std::string en = buildExtendedNewick();
        std::cout << "\n" << en << "\n\n";
        
    }
    
    if (gen == maxGen)
    {
        std::string nhxStr = buildNhxString();
        outStream << nhxStr << std::endl;
    }
}

template<class charType, class treeType>
std::string RevBayesCore::TreeCharacterHistoryNhxMonitor<charType, treeType>::buildNhxString(void)
{
    std::stringstream nhxStrm;
    
    // begin nexus file
    nhxStrm << "#NEXUS" << "\n\n";
    
    // Nhx settings block
    nhxStrm << "Begin Nhx;\n";
    nhxStrm << "\tdrawtype pie\n";
    nhxStrm << "\tmodeltype biogeography\n";
    nhxStrm << "\tareatype discrete\n";
    nhxStrm << "\tmaptype clean\n";
    nhxStrm << "\tpieslicestyle even\n";
    nhxStrm << "\tpiefillstyle outwards\n";
    nhxStrm << "\ttimestart -" << tree->getValue().getRoot().getAge() << "\n";
    nhxStrm << "\tmarkerradius " << 200 << "\n";
    nhxStrm << "\tminareaval " << 0.1 << "\n";
    nhxStrm << "End;\n\n";
    
    // bayarea-fig block
    nhxStrm << "Begin bayarea-fig;\n";
    nhxStrm << "\tmapheight\t100\n";
    nhxStrm << "\tmapwidth\t150\n";
    nhxStrm << "\tcanvasheight\t2000\n";
    nhxStrm << "\tcanvaswidth\t1000\n";
    nhxStrm << "\tminareaval\t0.15\n";
    nhxStrm << "\tareacolors black\n";
    nhxStrm << "\tareatypes";
    for (unsigned i = 0; i < numCharacters; i++)
        nhxStrm << " 1";
    nhxStrm << "\n";
    nhxStrm << "\tareanames Default\n";
    nhxStrm << "End;\n\n";
    
    // taxa block
    nhxStrm << "Begin taxa;\n";
    nhxStrm << "\tDimensions ntax=" << tree->getValue().getNumberOfTips() << ";\n";
    nhxStrm << "\tTaxlabels\n";
    for (unsigned i = 0; i < tree->getValue().getNumberOfNodes(); i++)
    {
        TopologyNode* p = &tree->getValue().getNode(i);
        if (p->isTip())
        {
            nhxStrm << "\t\t" << p->getName() << "\n";
        }
    }
    nhxStrm << "\t\t;\n";
    nhxStrm << "End;\n\n";
    
    // geo block
    nhxStrm << "Begin geo;\n";
    nhxStrm << "\tDimensions ngeo=" << numCharacters << ";\n";
    nhxStrm << "\tCoords\n";
    
    for (unsigned i = 0; i < numCharacters; i++)
    {
        nhxStrm << "\t\t" << i << "\t" << geographicCoordinates[i][0] << "\t" << geographicCoordinates[i][1];
        if (i < (numCharacters - 1))
            nhxStrm << ",";
        nhxStrm << "\n";
    }
    nhxStrm << "\t\t;\n";
    nhxStrm << "End;\n\n";
    
    // tree block
    nhxStrm << "Begin trees;\n";
    nhxStrm << "\tTranslate\n";
    for (unsigned i = 0; i < tree->getValue().getNumberOfNodes(); i++)
    {
        TopologyNode* p = &tree->getValue().getNode(i);
        if (p->isTip())
        {
            nhxStrm << "\t\t" << p->getIndex() << "\t" << p->getName();
            if (i < (tree->getValue().getNumberOfNodes() - 1))
                nhxStrm << ",";
            nhxStrm << "\n";
        }
    }
    nhxStrm << "\t\t;\n";
    
    // write tree string
    std::string treeStr = "";
    treeStr = buildExtendedNewick(); //buildExtendedNewick(&tree->getValue().getRoot());
    std::cout << treeStr << "\n";
    std::cout << "nhxStr\n" << treeStr << "\n";
    nhxStrm << "tree TREE1 = " << treeStr << ";\n";
    nhxStrm << "End;\n";
    
    std::cout << "[";
    for (size_t i = 0; i < numCharacters; i++)
    {
        if (i != 0)
            std::cout << ",";
        std::cout << (double)childCharacterCounts[0][i] / numSamples;
    }
    std::cout << "]\n";
    
    return nhxStrm.str();
}


/** open the file stream for printing */
template<class charType, class treeType>
void RevBayesCore::TreeCharacterHistoryNhxMonitor<charType, treeType>::openStream(void) {
    
    // open the stream to the file
    if (append)
        outStream.open( filename.c_str(), std::fstream::out | std::fstream::app);
    else
        outStream.open( filename.c_str(), std::fstream::out);
}

/** Print header for monitored values */
template<class charType, class treeType>
void RevBayesCore::TreeCharacterHistoryNhxMonitor<charType, treeType>::printHeader() {
    
    // do nothing
    ;
}

template<class charType, class treeType>
std::vector<unsigned int> RevBayesCore::TreeCharacterHistoryNhxMonitor<charType, treeType>::getChildCharacterCounts(int idx)
{
    return childCharacterCounts[idx];
}

template<class charType, class treeType>
std::vector<unsigned int> RevBayesCore::TreeCharacterHistoryNhxMonitor<charType, treeType>::getParentCharacterCounts(int idx)
{
    return parentCharacterCounts[idx];
}

template<class charType, class treeType>
long RevBayesCore::TreeCharacterHistoryNhxMonitor<charType, treeType>::getNumSamples(void)
{
    return numSamples;
}

template<class charType, class treeType>
void RevBayesCore::TreeCharacterHistoryNhxMonitor<charType, treeType>::swapNode(DagNode *oldN, DagNode *newN) {
    
    bool found = false;
    if ( oldN == tree )
    {
        tree = static_cast< TypedDagNode<treeType> *>(newN);
        found = true;
    }
    else if ( oldN == variable )
    {
        variable = static_cast<StochasticNode<AbstractCharacterData>* >(newN);
        found = true;
    }
    
//    for (size_t i = 0; i < branchHistories.size(); i++)
//    {
//        if (oldN == branchHistories[i])
//        {
//            branchHistories[i] = static_cast<StochasticNode<BranchHistory>* >(newN);
//            found = true;
//        }
//    }
//    
//    /*
//     if (found == false)
//     {
//     // error catching
//     if ( nodeVariables.find(oldN) == nodeVariables.end() ) {
//     throw RbException("Cannot replace DAG node in this monitor because the monitor doesn't hold this DAG node.");
//     }
//     
//     nodeVariables.erase( oldN );
//     nodeVariables.insert( newN );
//     }
//     */
    
    // delegate to base class
    Monitor::swapNode(oldN, newN);
}


#endif

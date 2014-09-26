//
//  PhylowoodNhxMonitor.h
//  rb_mlandis
//
//  Created by Michael Landis on 10/16/13.
//  Copyright (c) 2013 Michael Landis. All rights reserved.
//

#ifndef __rb_mlandis__PhylowoodNhxMonitor__
#define __rb_mlandis__PhylowoodNhxMonitor__

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
    
    class PhylowoodNhxMonitor : public Monitor {
        
    public:
        // Constructors and Destructors
        PhylowoodNhxMonitor(TypedDagNode<TimeTree> *t, std::vector< StochasticNode< BranchHistory >* > bh, std::vector<std::vector<double> > gc, int g, int mg, int burn, const std::string &fname, const std::string &del, bool pp=true, bool l=true, bool pr=true, bool ap=false, bool sm=true, bool sr=true);
        
        // new PhylowoodNhxMonitor( tau, bh_vector_stochastic, 10, filepath + "rb.tree_chars.txt", "\t"));
        
        PhylowoodNhxMonitor(const PhylowoodNhxMonitor& f);
        
        // basic methods
        PhylowoodNhxMonitor*          clone(void) const;                                                  //!< Clone the object
        
        // Monitor functions
        void                                monitor(long gen);                                                  //!< Monitor at generation gen
        void                                swapNode(DagNode *oldN, DagNode *newN);
        
        // FileMonitor functions
        void                                closeStream(void);                                                  //!< Close stream after finish writing
        void                                openStream(void);                                                   //!< Open the stream for writing
        void                                printHeader(void);                                                  //!< Print header
        std::vector<unsigned int>           getChildCharacterCounts(int idx);
        std::vector<unsigned int>           getParentCharacterCounts(int idx);
        long                                getNumSamples(void);
        
    private:
        std::string                         buildExtendedNewick();
        std::string                         buildExtendedNewick(TopologyNode* n);
        std::string                         buildCharacterHistoryString(TopologyNode* n, std::string brEnd="child");
        void updateCharacterCounts(TopologyNode* n, std::string brEnd="child");
        std::string                         buildNhxString(void);
        
        // the stream to print
        std::fstream                        outStream;
        
        // parameters
        TypedDagNode<TimeTree>*             tree;
        std::vector<StochasticNode<BranchHistory>* > branchHistories;
        std::set<DagNode *>                 nodeVariables;
        std::vector<std::vector<double> > geographicCoordinates;
        
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


#endif /* defined(__rb_mlandis__PhylowoodNhxMonitor__) */

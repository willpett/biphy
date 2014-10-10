/**
 * @file
 * This file contains the declaration of the Gtr model test class. 
 * This test runs an MCMC analysis on a phylogenetic model using the GTR subsitution model.
 * The data is read as an alignment from file. The tree is estimated as well!
 *
 *
 * @brief Declaration of the Gtr+G model test class
 *
 * (c) Copyright 2009-
 * @date Last modified: $Date$
 * @author The RevBayes Development Core Team
 * @license GPL version 3
 * @since Version 1.0, 2012-07-05
 *
 * $Id$
 */


#ifndef TestBranchHeterogeneousBinaryModel_H
#define TestBranchHeterogeneousBinaryModel_H

#include <string>
#include <vector>

namespace RevBayesCore {
    
    class Tree;
    class AbstractCharacterData;
    
    class TestBranchHeterogeneousBinaryModel {
        
    public:
    	TestBranchHeterogeneousBinaryModel(const std::string &datafile,
    										const std::string &treefile,
    										const std::string &name,
    										const std::string &outgroupfile,
    										int branchprior,
    										bool ras,
    										int heterogeneous,
    										bool dollo,
    										int mixture,
    										bool rootprior,
    										double rootmin,
    										double rootmax,
    										int every,
    										int until,
    										int numchains,
    										int swapInterval,
    										double deltaTemp,
    										double sigmaTemp,
    										bool saveall,
    										bool nexus);
    	TestBranchHeterogeneousBinaryModel(const std::string &name, const std::string &cvfile, bool ppred);
        virtual                                ~TestBranchHeterogeneousBinaryModel(void);                                                            //!< Virtual destructor
        
        bool                                    run();
        void                                    open();
        void                                    save();
        
    private:
        
        // members
        std::string                             dataFile;
        std::string								name;
        std::string                             treeFile;
        std::string                             outgroupFile;
        std::string                             cvfile;

        int										branchprior;
        bool									ras;
        int										heterogeneous;
        bool									dollo;
        int 									mixture;
        bool									ppred;
        bool									rootprior;
        double                                  rootmin;
        double                                  rootmax;
        int                                     every;
        int                                     until;
        int                                     numChains;
        int                                     swapInterval;
        double                                  deltaTemp;
        double                                  sigmaTemp;
        bool									saveall;
        bool									nexus;
        bool									readstream;
        
    };
    
}

#endif

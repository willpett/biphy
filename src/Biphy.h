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


#ifndef Biphy_H
#define Biphy_H

#include <string>
#include <vector>
#include "BinaryPartitionModel.h"
#include "ParallelMcmcmc.h"
    
class Biphy {
        
    public:

    	enum ModelType {DOLLO, HOMOGENEOUS, HIERARCHICAL, MIXTURE, DPP};
		enum RootPrior {FREE, RIGID, TRUNCATED};
		enum BranchPrior {EXPONENTIAL, DIRICHLET};


    	Biphy(const std::string &datafile,
    										const std::string &treefile,
    										const std::string &name,
    										const std::string &outgroupfile,
    										ModelType modeltype,
											BranchPrior branchprior,
											RootPrior rootprior,
											int correction,
    										int dgam,
    										int mixture,
    										double rootmin,
    										double rootmax,
    										int every,
    										int until,
    										int numchains,
    										int swapInterval,
    										double delta,
    										double sigma,
    										bool saveall,
    										bool nexus
											);
    	Biphy(const std::string &name, const std::string &cvfile, bool ppred, bool dolloMapping = false);
    	Biphy(const std::string &name);
        virtual                                ~Biphy(void);                                                            //!< Virtual destructor
        
        void									init();
        void                                    open();
        void                                    save();
        void                                    run();
        
    private:
        
        // members
        std::string                             dataFile;
        std::string								name;
        std::string                             treeFile;
        std::string                             outgroupFile;
        std::string                             cvfile;

        ModelType 								modeltype;
		BranchPrior 							branchprior;
		RootPrior 								rootprior;
		int 									correction;

        int										dgam;
        int 									mixture;
        bool									ppred;
        double                                  rootmin;
        double                                  rootmax;
        int                                     every;
        int                                     until;
        int                                     numChains;
        int                                     swapInterval;
        double                                  delta;
        double                                  sigma;
        bool									saveall;
        bool									nexus;
        bool									readstream;
        bool									restart;
        bool									dolloMapping;
        
        RevBayesCore::ParallelMcmcmc* 			mcmc;

    };

#endif

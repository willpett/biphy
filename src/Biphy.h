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
#include <cstdlib>
#include <iostream>

#include <unistd.h>

#include "BinaryCharacterData.h"
#include "ConstantNode.h"
#include "ContinuousStochasticNode.h"
#include "DeterministicNode.h"
#include "ParallelMcmcmc.h"
#include "Tree.h"

struct ModelPrior {
    enum Type {DOLLO, MK, HOMOGENEOUS, HIERARCHICAL, MIXTURE};
};

struct RootPrior {
    enum Type {FREE, RIGID, TRUNCATED};
};

struct BranchPrior {
    enum Type {DEFAULT, FIXED, STRICT, EXPONENTIAL, DIRICHLET};
};

class Biphy {
        
    public:

    	Biphy(const std::string name,
                const std::string datafile,
				const std::string cvfile,
                const std::string treefile,
                const std::string outgroupfile,
                ModelPrior::Type modeltype,
                BranchPrior::Type branchprior,
                RootPrior::Type rootprior,
                int correction,
                int dgam,
				int dbeta,
				bool asymmbeta,
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
                bool nexus,
				bool percoding
			);
    	
    	/* stream reading mode */
    	Biphy(const std::string name, const std::string cvfile, int ppred, bool dolloMapping, bool site, bool ancestral, int burnin = 0, int every = 1);
    	
    	/* restart mode */
    	Biphy(const std::string name, size_t stepping = 0);
        
        void									init();
        void                                    run();

    protected:

        void                                    open();
        void                                    save();

        virtual void                            precheck();
        virtual void                            readInputFiles();

        virtual void                            printConfiguration();
        virtual void                            initModel();
        virtual void                            initMCMC();
        
        // members
        std::string                             dataFile;
        std::string								name;
        std::string                             treeFile;

        int										correction;
        int										dgam;
        int										dbeta;
        int                                     every;
        int                                     until;
        int                                     numChains;
        int                                     swapInterval;
        double                                  delta;
        double                                  sigma;
        bool									saveall;
        bool									readstream;
        bool									restart;

        Model* 								    model;
        ParallelMcmcmc* 						mcmc;

        std::vector<Move*> 					    moves;
		std::vector<Monitor*> 				    monitors;
		std::vector<DagNode*> 				    monitoredNodes;

        std::vector<Tree*> 		                trees;
        BinaryCharacterData* 	                data;

        ConstantNode<double>*					    one;
        ConstantNode<double>*					    zero;

        ContinuousStochasticNode*					shape;
        DeterministicNode<std::vector<double> >*	site_rates;
        
        // members
        std::string                             outgroupFile;
        std::string                             cvfile;

        ModelPrior::Type                        modeltype;
        BranchPrior::Type                       branchprior;
        RootPrior::Type                         rootprior;

        int                                     mixture;
        int                                    ppred;
        double                                  rootmin;
        double                                  rootmax;
        bool                                    nexus;
        bool									percoding;
        bool                                    dolloMapping;
        bool                                    perSiteLnProbs;
        bool                                    ancestral;
        bool									asymmbeta;
        size_t                                  steppingStones;

        Clade                                   outgroup;

        BinaryCharacterData*                    cvdata;

};

#endif

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


#ifndef MultiBiphy_H
#define MultiBiphy_H

#include <string>
#include <vector>
#include <cstdlib>
#include <iostream>
#include "ParallelMcmcmc.h"
    
class MultiBiphy {
        
    public:

    	enum ModelType {DOLLO, REVERSIBLE};
		enum RatePrior {HOMOGENEOUS, HIERARCHICAL};

    	MultiBiphy(const std::string datafile,
					const std::string treefile,
					const std::string name,
					ModelType modeltype,
					RatePrior A_prior,
					RatePrior B_prior,
					RatePrior A_species,
					RatePrior B_species,
					int correction,
					int dgam,
					int every,
					int until,
					int numchains,
					int swapInterval,
					double delta,
					double sigma,
					bool saveall,
					bool ancestral
					);
    	MultiBiphy(const std::string name, bool ancestral);
    	MultiBiphy(const std::string name);
        
        virtual void							init() = 0;
        void                                    open();
        void                                    save();
        void                            		run();
        
    protected:
        
        // members
        std::string                             dataFile;
        std::string								name;
        std::string                             treeFile;

        ModelType 								modeltype;
		RatePrior 								A_prior;
		RatePrior 								B_prior;
		RatePrior 								A_species;
		RatePrior 								B_species;

		int 									correction;

        int										dgam;
        int                                     every;
        int                                     until;
        int                                     numChains;
        int                                     swapInterval;
        double                                  delta;
        double                                  sigma;
        bool									saveall;
        bool									readstream;
        bool									restart;
        bool									ancestral;
        
        RevBayesCore::ParallelMcmcmc* 			mcmc;

};

#endif

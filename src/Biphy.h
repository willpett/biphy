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

#include "BinaryCharEvoModel.h"

#include "ParallelMcmcmc.h"
#include "AbstractCharacterData.h"
#include "BranchLengthTree.h"
#include "ConstantNode.h"
#include "ContinuousStochasticNode.h"
    
class Biphy {
        
    public:

    	Biphy(const std::string name,
			const std::string datafile,
			const std::string treefile,
			int correction,
			int dgam,
			int every,
			int until,
			int numchains,
			int swapInterval,
			double delta,
			double sigma,
			bool saveall
			);
    	Biphy(const std::string &name);
        
        void									init();
        void                                    run();

    protected:

        void                                    open();
        void                                    save();
        
        virtual void                            openParams(std::ifstream&);
        virtual void                            saveParams(std::ofstream&);

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
        bool									missing;
        int                                     every;
        int                                     until;
        int                                     numChains;
        int                                     swapInterval;
        double                                  delta;
        double                                  sigma;
        bool									saveall;
        bool									readstream;
        bool									restart;
        
        std::string symbols;

        RevBayesCore::Model* 								model;
        RevBayesCore::ParallelMcmcmc* 						mcmc;

        std::vector<RevBayesCore::Move*> 					moves;
		std::vector<RevBayesCore::Monitor*> 				monitors;
		std::vector<RevBayesCore::DagNode*> 				monitoredNodes;

        std::vector<RevBayesCore::BranchLengthTree*> 		trees;
        std::vector<RevBayesCore::AbstractCharacterData*> 	data;

        RevBayesCore::ConstantNode<double>*					one;
        RevBayesCore::ConstantNode<double>*					zero;

        RevBayesCore::ContinuousStochasticNode*					shape;
        RevBayesCore::DeterministicNode<std::vector<double> >*	site_rates;
        RevBayesCore::StochasticNode<double>*					xi;

    };

#endif

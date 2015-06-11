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

#include <iomanip>
#include "Biphy.h"
    
class MultiBiphy : public Biphy {
        
    public:

    	enum ModelType {DOLLO, REVERSIBLE};
		enum RatePrior {HOMOGENEOUS, HIERARCHICAL};

		MultiBiphy(const std::string name,
					const std::string datafile,
					const std::string treefile,
					ModelType modeltype,
					RatePrior A_prior,
					RatePrior B_prior,
					RatePrior A_species,
					RatePrior B_species,
					bool ancestral,
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
    	MultiBiphy(const std::string name, bool ancestral);
    	MultiBiphy(const std::string name);
        
    	virtual void		readInputFiles();
		virtual void		printConfiguration();
		virtual void		initModel();
		virtual void		initMCMC();

		virtual RevBayesCore::BinaryCharEvoModel<RevBayesCore::BranchLengthTree>*		initSubModel(size_t index) = 0;

		void                openParams(std::ifstream&);
		void                saveParams(std::ofstream&);
        
    protected:

		std::vector<RevBayesCore::StochasticNode<RevBayesCore::AbstractCharacterData>* > dataNodes;

        ModelType 								modeltype;
		RatePrior 								A_prior;
		RatePrior 								B_prior;
		RatePrior 								A_species;
		RatePrior 								B_species;

        bool									ancestral;

        std::vector<std::map<size_t, size_t> >	node2speciesMaps;
        std::map<size_t, size_t>				speciesIndex;
        std::vector<size_t>						rootSpecies;

		std::stringstream 						tmp_name;
		size_t numSpecies;
		size_t numTrees;
		size_t wt;
		size_t ws;

};

#endif

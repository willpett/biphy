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


#ifndef MultiBiphyLambdaMu_H
#define MultiBiphyLambdaMu_H

#include "MultiBiphy.h"
    
class MultiBiphyLambdaMu : public MultiBiphy {
        
    public:
		MultiBiphyLambdaMu(const std::string datafile,
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
		MultiBiphyLambdaMu(const std::string name, bool ancestral);
		MultiBiphyLambdaMu(const std::string name);

        void									init();

    };

#endif

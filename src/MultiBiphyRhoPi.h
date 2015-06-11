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


#ifndef MultiBiphyRhoPi_H
#define MultiBiphyRhoPi_H

#include "MultiBiphy.h"

namespace RevBayesCore {

	class MultiBiphyRhoPi : public MultiBiphy {

		public:

			MultiBiphyRhoPi(const std::string name,
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
			MultiBiphyRhoPi(const std::string name, bool ancestral);
			MultiBiphyRhoPi(const std::string name);

			void		printConfiguration();
			void		initModel( void );
			RevBayesCore::BinaryCharEvoModel<RevBayesCore::BranchLengthTree>*	initSubModel(size_t index);
			void		initMCMC( void );

		protected:
			RatePrior rho_prior;
			RatePrior rho_species;
			RatePrior pi_prior;
			RatePrior pi_species;

			// base frequencies prior
			DeterministicNode< RateMatrix >* q;
			std::vector<const TypedDagNode< RateMatrix >* > q_vector;

			TypedDagNode<std::vector<double> >*			rf;

			std::vector<StochasticNode<double >* > 		pi_vector;
			std::vector<DeterministicNode<double >* > 	pi0_vector;
			std::vector<StochasticNode<double >* > 		rho_vector;

			std::vector<StochasticNode<double >* > 		pi_x_vector;
			std::vector<DeterministicNode<double >* > 	pi0_x_vector;
			std::vector<StochasticNode<double >* >		rho_x_vector;

			std::stringstream tmp_name;

			TypedDagNode<double> *pi;
			DeterministicNode<double> *pi0;

			StochasticNode<double> *alpha;
			StochasticNode<double> *beta;

			StochasticNode<double> *alpha_x;
			StochasticNode<double> *beta_x;

			StochasticNode<double> *rho;
			StochasticNode<double> *rho_variance;

			StochasticNode<double> *rho_x_shape;

	};
}

#endif

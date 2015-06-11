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
namespace RevBayesCore {

	class MultiBiphyLambdaMu : public MultiBiphy {

		public:
			MultiBiphyLambdaMu(const std::string name,
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
			MultiBiphyLambdaMu(const std::string name, bool ancestral);
			MultiBiphyLambdaMu(const std::string name);

			void		printConfiguration();
			void		initModel( void );
			RevBayesCore::BinaryCharEvoModel<RevBayesCore::BranchLengthTree>*	initSubModel(size_t index);
			void		initMCMC( void );

		protected:
			RatePrior lambda_prior;
			RatePrior lambda_species;
			RatePrior mu_prior;
			RatePrior mu_species;

			StochasticNode<double> *					lambda;
			StochasticNode<double> *					mu;
			TypedDagNode<double>*					pi;

			std::vector<StochasticNode<double > *> 		lambda_vector;
			std::vector<StochasticNode<double > *> 		mu_vector;

			std::vector<StochasticNode<double > *> 		lambda_x_vector;
			std::vector<StochasticNode<double > *> 		mu_x_vector;

			std::vector<DeterministicNode<double > *> 	pi_vector;

			StochasticNode<double>*						lambda_variance;
			StochasticNode<double>*						mu_variance;

			StochasticNode<double>*						lambda_x_shape;
			StochasticNode<double>*						mu_x_shape;


			DeterministicNode< RateMatrix >* q = NULL;
			std::vector<const TypedDagNode< RateMatrix >* > q_vector;

			TypedDagNode<std::vector<double> > *rf;
			std::vector<TypedDagNode<std::vector<double> >* > rf_vector;
	};
}

#endif

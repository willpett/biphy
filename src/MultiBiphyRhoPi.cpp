#include "AbstractCharacterData.h"
#include "BinaryAddition.h"
#include "BinaryDivision.h"
//#include "BinaryDolloCompatibleMonitor.h"
#include "BinaryMultiplication.h"
#include "BinaryMissingCharEvoModel.h"
#include "BinarySubtraction.h"
#include "BetaDistribution.h"
#include "BetaSimplexMove.h"
#include "ConstantNode.h"
#include "ConstantFunction.h"
#include "DeterministicNode.h"
#include "DolloBinaryMissingCharEvoModel.h"
#include "ExponentialDistribution.h"
#include "FileMonitor.h"
#include "FreeBinaryRateMatrixFunction.h"
#include "FreeBinaryRateMatrixVectorFunction.h"
#include "GammaDistribution.h"
#include "ParallelMcmcmc.h"
#include "Mcmc.h"
#include "Model.h"
#include "Monitor.h"
#include "Move.h"
#include "NclReader.h"
#include "NHXTreeReader.h"
#include "NegativeBinomialDistribution.h"
#include "NewickTreeMonitor.h"
#include "NexusTreeMonitor.h"
#include "ParallelMcmcmc.h"
#include "PoissonDistribution.h"
#include "QuantileFunction.h"
#include "RbSettings.h"
#include "RbStatisticsHelper.h"
#include "RbVectorFunction.h"
#include "ScaleMove.h"
#include "ScreenMonitor.h"
#include "SimplexSingleElementScale.h"
#include "SlidingMove.h"
#include "VectorFunction.h"

#include <unistd.h>
#include "MultiBiphyRhoPi.h"

using namespace RevBayesCore;

MultiBiphyRhoPi::MultiBiphyRhoPi(const std::string n,
						const std::string df,
						const std::string t,
						ModelType mt,
						RatePrior Ap,
						RatePrior Bp,
						RatePrior As,
						RatePrior Bs,
						bool anc,
						int c,
						int d,
						int e,
						int u,
						int num,
						int swa,
						double de,
						double si,
						bool sav) :
		MultiBiphy(n,
		df,
		t,
		mt,
		Ap,
		Bp,
		As,
		Bs,
		anc,
		c,
		d,
		e,
		u,
		num,
		swa,
		de,
		si,
		sav)
{
	rho_prior 	= A_prior;
	rho_species = A_species;
	pi_prior 	= B_prior;
	pi_species 	= B_species;
}

MultiBiphyRhoPi::MultiBiphyRhoPi(const std::string n, bool anc) : MultiBiphy(n,anc)
{
}

MultiBiphyRhoPi::MultiBiphyRhoPi(const std::string n) : MultiBiphy(n)
{
}

void MultiBiphyRhoPi::printConfiguration( void ) {

	std::cout << std::endl;

    if(pi_prior == HOMOGENEOUS && pi_species == HOMOGENEOUS){
    	std::cout << "pi ~ Beta(1,1)\n";
    }else{
    	std::cout << "alpha ~ exp(1)\n";
		std::cout << "beta ~ exp(1)\n";
		std::cout << "pi_i ~ Beta(alpha,beta)\n";

		if(pi_prior == HIERARCHICAL && pi_species == HIERARCHICAL){
			std::cout << std::endl;

			std::cout << "species modulators:\n";
			std::cout << "alpha_x ~ exp(1)\n";
			std::cout << "beta_x ~ exp(1)\n";
			std::cout << "pi_x_i ~ Beta(alpha_x,beta_x)\n";
    	}
    }

    std::cout << "\n";

    if(rho_prior == HOMOGENEOUS && rho_species == HOMOGENEOUS){
		std::cout << "rho ~ exp(1)\n";
	}else{
		std::cout << "rho_mean ~ exp(1)\n";
		std::cout << "rho_variance ~ exp(1)\n";
		std::cout << "rho_i ~ gamma with mean = mu_mean and variance = mu_variance\n";

		if(rho_prior == HIERARCHICAL && rho_species == HIERARCHICAL){
			std::cout << std::endl;

			std::cout << "species multipliers:\n";
			std::cout << "rho_x_shape ~ exp(1)\n";
			std::cout << "rho_x_i ~ gamma(mu_x_shape, mu_x_shape)\n";
		}
	}

    MultiBiphy::printConfiguration();
}

void MultiBiphyRhoPi::initModel( void ) {

    //////////////////////
    // first the priors //
    //////////////////////

	Biphy::initModel();

	/* lambda prior */
	std::string n;
	if(pi_prior == HOMOGENEOUS){
		if(pi_species == HIERARCHICAL){
			alpha = new StochasticNode<double>("alpha", new ExponentialDistribution(one) );
			beta = new StochasticNode<double>("beta", new ExponentialDistribution(one) );

			DeterministicNode<double>* denom = new DeterministicNode<double>("denom", new BinaryAddition<double,double,double>(alpha,beta));
			pi = new DeterministicNode<double>("mean_pi", new BinaryDivision<double,double,double>(alpha,denom));

			pi0 = new DeterministicNode<double>("mean_pi0", new BinarySubtraction<double,double,double>(one,pi));

			for(std::map<size_t, size_t>::iterator it = speciesIndex.begin(); it != speciesIndex.end(); it++){
				tmp_name.str("");
				tmp_name << "pi(" << setfill('0') << setw(ws) << it->first << ")";
				pi_vector.push_back(new StochasticNode<double>( tmp_name.str(), new BetaDistribution(alpha,beta) ) );

				tmp_name.str("");
				tmp_name << "pi0(" << setfill('0') << setw(wt) << it->first << ")";
				pi0_vector.push_back(new DeterministicNode<double>(tmp_name.str(), new BinarySubtraction<double,double,double>(one,pi_vector.back())));
			}
		}else{
			pi = new StochasticNode<double>("pi", new BetaDistribution(one,one) );
			pi0 = new DeterministicNode<double>("pi0", new BinarySubtraction<double,double,double>(one,pi));
		}
	}else{
		alpha = new StochasticNode<double>("alpha", new ExponentialDistribution(one) );
		beta = new StochasticNode<double>("beta", new ExponentialDistribution(one) );

		DeterministicNode<double>* denom = new DeterministicNode<double>("mean_lambda_squared", new BinaryAddition<double,double,double>(alpha,beta));
		pi = new DeterministicNode<double>("mean_pi", new BinaryDivision<double,double,double>(alpha,denom));

		pi0 = new DeterministicNode<double>("mean_pi0", new BinarySubtraction<double,double,double>(one,pi));

		for (size_t i = 0 ; i < numTrees ; i++ ) {
			tmp_name.str("");
			tmp_name << "pi(" << setfill('0') << setw(wt) << i << ")";
			pi_vector.push_back(new StochasticNode<double>( tmp_name.str(), new BetaDistribution(alpha,beta) ) );

			tmp_name.str("");
			tmp_name << "pi0(" << setfill('0') << setw(wt) << i << ")";
			pi0_vector.push_back(new DeterministicNode<double>(tmp_name.str(), new BinarySubtraction<double,double,double>(one,pi_vector.back())));
		}
		if(pi_species == HIERARCHICAL){
			alpha_x = new StochasticNode<double>("alpha_x", new ExponentialDistribution(one) );
			beta_x = new StochasticNode<double>("beta_x", new ExponentialDistribution(one) );

			for(std::map<size_t, size_t>::iterator it = speciesIndex.begin(); it != speciesIndex.end(); it++){
				size_t j = it->first;
				tmp_name.str("");
				tmp_name << "pi_x(" << setfill('0') << setw(ws) << j << ")";
				pi_x_vector.push_back(new StochasticNode<double>( tmp_name.str(), new BetaDistribution(alpha_x,beta_x) ) );

				tmp_name.str("");
				tmp_name << "pi0_x(" << setfill('0') << setw(wt) << j << ")";
				pi0_x_vector.push_back(new DeterministicNode<double>(tmp_name.str(), new BinarySubtraction<double,double,double>(one,pi_x_vector.back())));
			}
		}
	}

	/* rho prior */
	if(rho_prior == HOMOGENEOUS){
		if(rho_species == HIERARCHICAL)
			n = "mean_rho";
		else
			n = "rho";
		rho = new StochasticNode<double>(n, new ExponentialDistribution(one) );
		if(rho_species == HIERARCHICAL){
			DeterministicNode<double>* mean_squared = new DeterministicNode<double>("mean_rho_squared", new BinaryMultiplication<double,double,double>(rho,rho));
			rho_variance = new StochasticNode<double>("rho_variance", new ExponentialDistribution(one) );

			// alpha = (μ/σ)^2, beta = μ/σ^2
			DeterministicNode<double>* alpha = new DeterministicNode<double>("rho_alpha", new BinaryDivision<double,double,double>(mean_squared,rho_variance));
			DeterministicNode<double>* beta = new DeterministicNode<double>("rho_beta", new BinaryDivision<double,double,double>(rho,rho_variance));

			for(std::map<size_t, size_t>::iterator it = speciesIndex.begin(); it != speciesIndex.end(); it++){
				tmp_name.str("");
				tmp_name << "rho(" << setfill('0') << setw(ws) << it->first << ")";
				rho_vector.push_back(new StochasticNode<double>( tmp_name.str(), new GammaDistribution(alpha,beta) ) );
			}
		}
	}else{
		rho = new StochasticNode<double>("mean_rho", new ExponentialDistribution(one) );

		DeterministicNode<double>* mean_squared = new DeterministicNode<double>("mean_rho_squared", new BinaryMultiplication<double,double,double>(rho,rho));
		rho_variance = new StochasticNode<double>("rho_variance", new ExponentialDistribution(one) );

		// alpha = (μ/σ)^2, beta = μ/σ^2
		DeterministicNode<double>* alpha = new DeterministicNode<double>("rho_alpha", new BinaryDivision<double,double,double>(mean_squared,rho_variance));
		DeterministicNode<double>* beta = new DeterministicNode<double>("rho_beta", new BinaryDivision<double,double,double>(rho,rho_variance));

		for (size_t i = 0 ; i < numTrees ; i++ ) {
			tmp_name.str("");
			tmp_name << "rho(" << setfill('0') << setw(wt) << i << ")";
			rho_vector.push_back(new StochasticNode<double>( tmp_name.str(), new GammaDistribution(alpha,beta) ) );
		}
		if(rho_species == HIERARCHICAL){
			rho_x_shape = new StochasticNode<double>("rho_x_shape", new ExponentialDistribution(one) );
			for(std::map<size_t, size_t>::iterator it = speciesIndex.begin(); it != speciesIndex.end(); it++){
				tmp_name.str("");
				tmp_name << "rho_x(" << setfill('0') << setw(ws) << it->first << ")";
				rho_x_vector.push_back(new StochasticNode<double>( tmp_name.str(), new GammaDistribution(rho_x_shape,rho_x_shape) ) );
			}
		}
	}

	std::vector<const TypedDagNode<double> *> rf_vec(2);
	if(modeltype == DOLLO){
		q = new DeterministicNode<RateMatrix>( "q", new FreeBinaryRateMatrixFunction(zero) );

		rf_vec[1] = zero;
		rf_vec[0] = one;
		rf = new DeterministicNode< std::vector< double > >( "rf", new VectorFunction< double >( rf_vec ) );
	}else{
		if(pi_prior == HOMOGENEOUS && pi_species == HOMOGENEOUS){
			q = new DeterministicNode<RateMatrix>( "q", new FreeBinaryRateMatrixFunction(pi) );

			rf_vec[1] = pi;
			rf_vec[0] = pi0;
			rf = new DeterministicNode< std::vector< double > >( "rf", new VectorFunction< double >( rf_vec ) );
		}else if(pi_prior == HOMOGENEOUS){
			for(size_t i = 0; i < numSpecies; i++){
				tmp_name.str("");
				tmp_name << "q(" << setfill('0') << setw(ws) << i << ")";
				q_vector.push_back(new DeterministicNode<RateMatrix>( tmp_name.str(), new FreeBinaryRateMatrixFunction(pi_vector[i]) ));
			}
		}
	}

	std::cout << "prior okay\n";

	MultiBiphy::initModel();
}

BinaryCharEvoModel<BranchLengthTree>* MultiBiphyRhoPi::initSubModel( size_t i ) {

		const TypedDagNode<double>* pi_i;
		const DeterministicNode<double>* pi0_i;
		if(pi_prior == HIERARCHICAL){
			pi_i = pi_vector[i];
			pi0_i = pi0_vector[i];
		}else{
			pi_i = pi;
			pi0_i = pi0;
		}

		const StochasticNode<double>* rho_i;
		if(rho_prior == HIERARCHICAL)
			rho_i = rho_vector[i];
		else
			rho_i = rho;

		size_t numNodes = trees[i]->getNumberOfNodes();


		tmp_name.str("");
		tmp_name << "psi(" << setfill('0') << setw(wt) << i << ")";
		ConstantNode<BranchLengthTree> *psi = new ConstantNode<BranchLengthTree>( tmp_name.str(), trees[i] );

		BinaryCharEvoModel<BranchLengthTree> *charModel;

		//std::cout << "built model " << i << std::endl;
		// check for species specific rates
		std::vector<TypedDagNode<double>* > scaled_pi;
		std::vector<TypedDagNode<double>* > scaled_pi0;

		if(pi_species == HIERARCHICAL){
			if(pi_prior == HIERARCHICAL){
				for (size_t j = 0 ; j < numSpecies; j++ ) {
					tmp_name.str("");
					tmp_name << "a(" << setfill('0') << setw(wt) << i << "," << setfill('0') << setw(ws) << j << ")";
					DeterministicNode<double>* a = new DeterministicNode<double>(tmp_name.str(), new BinaryMultiplication<double,double,double>(pi_i,pi_x_vector[j]));

					tmp_name.str("");
					tmp_name << "b(" << setfill('0') << setw(wt) << i << "," << setfill('0') << setw(ws) << j << ")";
					DeterministicNode<double>* b = new DeterministicNode<double>(tmp_name.str(), new BinaryMultiplication<double,double,double>(pi0_i,pi0_x_vector[j]));

					tmp_name.str("");
					tmp_name << "ab(" << setfill('0') << setw(wt) << i << "," << setfill('0') << setw(ws) << j << ")";
					DeterministicNode<double>* ab = new DeterministicNode<double>(tmp_name.str(), new BinaryAddition<double,double,double>(a,b));

					tmp_name.str("");
					tmp_name << "pi(" << setfill('0') << setw(wt) << i << "," << setfill('0') << setw(ws) << j << ")";
					scaled_pi.push_back(new DeterministicNode<double>(tmp_name.str(), new BinaryDivision<double,double,double>(a,ab)));

					tmp_name.str("");
					tmp_name << "pi0(" << setfill('0') << setw(wt) << i << "," << setfill('0') << setw(ws) << j << ")";
					scaled_pi0.push_back(new DeterministicNode<double>(tmp_name.str(), new BinarySubtraction<double,double,double>(one,scaled_pi.back())));
				}
			}else{
				for (size_t j = 0 ; j < numSpecies; j++ ) {
					scaled_pi.push_back(pi_vector[speciesIndex[node2speciesMaps[i][j]]]);
					scaled_pi0.push_back(pi0_vector[speciesIndex[node2speciesMaps[i][j]]]);
				}
			}
		}

		std::vector<const TypedDagNode<double> *> scaled_rho;

		if(rho_species == HIERARCHICAL){
			for (size_t j = 0 ; j < numSpecies; j++ ) {
				if(rho_prior == HIERARCHICAL){
					tmp_name.str("");
					tmp_name << "rho(" << setfill('0') << setw(wt) << i << "," << setfill('0') << setw(ws) << j << ")";
					scaled_rho.push_back(new DeterministicNode<double>(tmp_name.str(), new BinaryMultiplication<double,double,double>(rho_i,rho_x_vector[j])));
				}else{
					scaled_rho.push_back(rho_vector[j]);
				}
			}
		}

		// modeltype switch
		if(modeltype == REVERSIBLE){
			/* REVERSIBLE MODEL */

			charModel = new BinaryMissingCharEvoModel<BranchLengthTree>(psi, true, data[i]->getNumberOfIncludedCharacters(), correction);
			if(missing)
				((BinaryMissingCharEvoModel<BranchLengthTree>*)charModel)->setMissingRate(xi);

			//std::cout << "built model " << i << std::endl;
			// check for species specific rates
			if(pi_prior == HIERARCHICAL && pi_species == HIERARCHICAL){
				for (size_t j = 0 ; j < numSpecies; j++ ) {
					tmp_name.str("");
					tmp_name << "q(" << setfill('0') << setw(wt) << i << "," << setfill('0') << setw(ws) << j << ")";
					q_vector.push_back(new DeterministicNode<RateMatrix>( tmp_name.str(), new FreeBinaryRateMatrixFunction(scaled_pi[j]) ));
				}
			}else if(pi_prior == HIERARCHICAL && pi_species == HOMOGENEOUS){
				q = new DeterministicNode<RateMatrix>( "q", new FreeBinaryRateMatrixFunction(pi_i) );
			}

			if(pi_species == HOMOGENEOUS){
				charModel->setRateMatrix( q );

				if(pi_prior == HIERARCHICAL){
					std::vector<const TypedDagNode<double> *> rf_vec(2);
					rf_vec[1] = pi_i;
					rf_vec[0] = pi0_i;

					tmp_name.str("");
					tmp_name << "rf(" << setfill('0') << setw(wt) << i << ")";
					rf = new DeterministicNode< std::vector< double > >( tmp_name.str(), new VectorFunction< double >( rf_vec ) );
				}
			}else{
				std::vector<const TypedDagNode< RateMatrix >* > qs;

				for (size_t j = 0 ; j < numNodes; j++ )
					qs.push_back(q_vector[speciesIndex[node2speciesMaps[i][j]]]);

				tmp_name.str("");
				tmp_name << "q_vector(" << setfill('0') << setw(wt) << i << ")";
				DeterministicNode< RbVector< RateMatrix > >* qs_node = new DeterministicNode< RbVector< RateMatrix > >( tmp_name.str(), new RbVectorFunction<RateMatrix>(qs) );
				charModel->setRateMatrix( qs_node );

				std::vector<const TypedDagNode<double> *> rf_vec(2);
				if(pi_prior == HOMOGENEOUS){
					rf_vec[1] = pi_vector[speciesIndex[rootSpecies[i]]];
					rf_vec[0] = pi0_vector[speciesIndex[rootSpecies[i]]];
				}else{
					rf_vec[1] = scaled_pi[speciesIndex[rootSpecies[i]]];
					rf_vec[0] = scaled_pi0[speciesIndex[rootSpecies[i]]];
				}

				tmp_name.str("");
				tmp_name << "rf(" << setfill('0') << setw(wt) << i << ")";
				rf = new DeterministicNode< std::vector< double > >( tmp_name.str(), new VectorFunction< double >( rf_vec ) );
			}

			if(rho_species == HIERARCHICAL){
				std::vector<const TypedDagNode<double >* > rhos;

				tmp_name.str("");
				tmp_name << "rhos(" << setfill('0') << setw(wt) << i << ")";
				for (size_t j = 0 ; j < numNodes; j++ )
					rhos.push_back(scaled_rho[speciesIndex[node2speciesMaps[i][j]]]);

				TypedDagNode<std::vector<double> >* rho_clock = new DeterministicNode< std::vector< double > >( tmp_name.str(), new VectorFunction<double>( rhos ) );

				charModel->setClockRate(rho_clock);
			}else{
				charModel->setClockRate(rho_i);
			}
		}else{
			/* DOLLO MODEL */

			tmp_name.str("");
			tmp_name << "tau(" << setfill('0') << setw(wt) << i << ")";
			ConstantNode<Topology> *tau = new ConstantNode<Topology>( tmp_name.str(), (Topology*)&(psi->getValue().getTopology()) );

			charModel = new DolloBinaryMissingCharEvoModel<BranchLengthTree>(psi, tau, true, data[i]->getNumberOfIncludedCharacters(), correction);
			if(missing)
				((DolloBinaryMissingCharEvoModel<BranchLengthTree>*)charModel)->setMissingRate(xi);

			/* set clock rates */
			charModel->setRateMatrix( q );
		}

		charModel->setRootFrequencies( rf );

		if(dgam > 1)
			charModel->setSiteRates( site_rates );

		return charModel;
}

void MultiBiphyRhoPi::initMCMC( void ) {

	monitoredNodes.push_back(pi);

	if(pi_prior == HOMOGENEOUS && pi_species == HOMOGENEOUS){
		moves.push_back( new BetaSimplexMove((StochasticNode<double>*)pi, 1.0, true, 1.0 ) );
	}else{
		for (size_t i = 0 ; i < pi_vector.size() ; i++ ){
			moves.push_back( new BetaSimplexMove(pi_vector[i], 1.0, true, 1.0 ) );
			if(pi_prior == HOMOGENEOUS)
				monitoredNodes.push_back(pi_vector[i]);
		}

		if(pi_prior == HIERARCHICAL && pi_species == HIERARCHICAL){
			moves.push_back( new ScaleMove(alpha_x, 1.0, true, 1.0 ) );
			moves.push_back( new ScaleMove(beta_x, 1.0, true, 1.0 ) );
			for (size_t i = 0 ; i < numSpecies ; i++ ){
				moves.push_back( new BetaSimplexMove(pi_x_vector[i], 1.0, true, 1.0 ) );
				monitoredNodes.push_back(pi_x_vector[i]);
			}
		}

		moves.push_back( new ScaleMove(alpha, 1.0, true, 1.0 ) );
		moves.push_back( new ScaleMove(beta, 1.0, true, 1.0 ) );

		monitoredNodes.push_back(alpha);
		monitoredNodes.push_back(beta);

		if(pi_prior == HIERARCHICAL && pi_species == HIERARCHICAL){
			monitoredNodes.push_back(alpha_x);
			monitoredNodes.push_back(beta_x);
		}
	}

	moves.push_back( new ScaleMove(rho, 1.0, true, 1.0 ) );
	monitoredNodes.push_back(rho);

	if(!(rho_prior == HOMOGENEOUS && rho_species == HOMOGENEOUS)){
		moves.push_back( new ScaleMove(rho_variance, 1.0, true, 1.0 ) );
		monitoredNodes.push_back(rho_variance);

		for (size_t i = 0 ; i < rho_vector.size() ; i++ ){
			moves.push_back( new ScaleMove(rho_vector[i], 1.0, true, 1.0 ) );
			if(rho_prior == HOMOGENEOUS)
				monitoredNodes.push_back(rho_vector[i]);
		}

		if(rho_prior == HIERARCHICAL && rho_species == HIERARCHICAL){
			for (size_t i = 0 ; i < numSpecies ; i++ ){
				moves.push_back( new ScaleMove(rho_x_vector[i], 1.0, true, 1.0 ) );
				monitoredNodes.push_back(rho_x_vector[i]);
			}
			moves.push_back( new ScaleMove(rho_x_shape, 1.0, true, 1.0 ) );
			monitoredNodes.push_back(rho_x_shape);
		}
	}

	MultiBiphy::initMCMC();
}

#include "AbstractCharacterData.h"
#include "BinaryAddition.h"
#include "BinaryDivision.h"
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
#include "LnCorrectionFunction.h"
#include "ParallelMcmcmc.h"
#include "MappingMonitor.h"
#include "Mcmc.h"
#include "Model.h"
#include "Monitor.h"
#include "GibbsMove.h"
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
#include "SumFunction.h"
#include "VectorFunction.h"

#include "MultiBiphyLambdaMu.h"

using namespace RevBayesCore;

MultiBiphyLambdaMu::MultiBiphyLambdaMu(const std::string n,
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
	lambda_prior 	= A_prior;
	lambda_species 	= A_species;
	mu_prior 		= B_prior;
	mu_species 		= B_species;
}

MultiBiphyLambdaMu::MultiBiphyLambdaMu(const std::string n, bool anc) : MultiBiphy(n,anc)
{
}

MultiBiphyLambdaMu::MultiBiphyLambdaMu(const std::string n) : MultiBiphy(n)
{
}

void MultiBiphyLambdaMu::printConfiguration( void ) {

	std::cout << std::endl;

    if(lambda_prior == HOMOGENEOUS && lambda_species == HOMOGENEOUS){
    	std::cout << "lambda ~ exp(1)\n";
    }else{
    	std::cout << "mean_lambda ~ exp(1)\n";
		std::cout << "lambda_variance ~ exp(1)\n";
		std::cout << "lambda_i ~ gamma of mean = mean_lambda and var = lambda_variance\n\n";

		if(lambda_prior == HIERARCHICAL && lambda_species == HIERARCHICAL){
			std::cout << "species multipliers:\n";
			std::cout << "lambda_x_shape ~ exp(1)\n";
			std::cout << "lambda_x_i ~ gamme(lambda_x_shape,lambda_x_shape)\n\n";
    	}
    }

    std::cout << "\n";

    if(mu_prior == HOMOGENEOUS && mu_species == HOMOGENEOUS){
		std::cout << "mu ~ exp(1)\n";
	}else{
		std::cout << "mean_mu ~ exp(1)\n";
		std::cout << "mu_variance ~ exp(1)\n";
		std::cout << "mu_i ~ gamma of mean = mean_mu and var = mu_variance\n\n";

		if(mu_prior == HIERARCHICAL && mu_species == HIERARCHICAL){
			std::cout << "species multipliers:\n";
			std::cout << "mu_x_shape ~ exp(1)\n";
			std::cout << "mu_x_i ~ gamme(mu_x_shape,mu_x_shape)\n\n";
		}
	}

    MultiBiphy::printConfiguration();
}

void MultiBiphyLambdaMu::initModel( void ) {

    //////////////////////
    // first the priors //
    //////////////////////

	Biphy::initModel();

	/* lambda prior */
	if(lambda_prior == HOMOGENEOUS){
		if(lambda_species == HIERARCHICAL){
			lambda = new StochasticNode<double>("mean_lambda", new ExponentialDistribution(one) );
			lambda_variance = new StochasticNode<double>("lambda_variance", new ExponentialDistribution(one) );

			DeterministicNode<double>* mean_squared = new DeterministicNode<double>("mean_squared_lambda", new BinaryMultiplication<double,double,double>(lambda,lambda) );
			DeterministicNode<double>* alpha = new DeterministicNode<double>("lambda_alpha", new BinaryDivision<double,double,double>(mean_squared,lambda_variance) );
			DeterministicNode<double>* beta = new DeterministicNode<double>("lambda_beta", new BinaryDivision<double,double,double>(lambda_variance,lambda) );

			for(std::map<size_t, size_t>::iterator it = speciesIndex.begin(); it != speciesIndex.end(); it++){
				tmp_name.str("");
				tmp_name << "lambda(" << setfill('0') << setw(ws) << it->first << ")";
				lambda_vector.push_back(new StochasticNode<double>( tmp_name.str(), new GammaDistribution(alpha,beta) ) );
			}
		}else{
			lambda = new StochasticNode<double>("lambda", new ExponentialDistribution(one) );
		}
	}else{
		lambda = new StochasticNode<double>("mean_lambda", new ExponentialDistribution(one) );
		lambda_variance = new StochasticNode<double>("lambda_variance", new ExponentialDistribution(one) );

		DeterministicNode<double>* mean_squared = new DeterministicNode<double>("mean_squared_lambda", new BinaryMultiplication<double,double,double>(lambda,lambda) );
		DeterministicNode<double>* alpha = new DeterministicNode<double>("lambda_alpha", new BinaryDivision<double,double,double>(mean_squared,lambda_variance) );
		DeterministicNode<double>* beta = new DeterministicNode<double>("lambda_beta", new BinaryDivision<double,double,double>(lambda_variance,lambda) );

		for (size_t i = 0 ; i < numTrees ; i++ ) {
			tmp_name.str("");
			tmp_name << "lambda(" << setfill('0') << setw(ws) << i << ")";
			lambda_vector.push_back(new StochasticNode<double>( tmp_name.str(), new GammaDistribution(alpha,beta) ) );
		}
		if(lambda_species == HIERARCHICAL){
			lambda_x_shape = new StochasticNode<double>("lambda_x_shape", new ExponentialDistribution(one) );

			for(std::map<size_t, size_t>::iterator it = speciesIndex.begin(); it != speciesIndex.end(); it++){
				size_t j = it->first;
				tmp_name.str("");
				tmp_name << "lambda_x(" << setfill('0') << setw(ws) << j << ")";
				lambda_x_vector.push_back(new StochasticNode<double>( tmp_name.str(), new GammaDistribution(lambda_x_shape,lambda_x_shape) ) );
			}
		}
	}

	/* mu prior */
	if(mu_prior == HOMOGENEOUS){
		if(mu_species == HIERARCHICAL){
			mu = new StochasticNode<double>("mean_mu", new ExponentialDistribution(one) );
			mu_variance = new StochasticNode<double>("mu_variance", new ExponentialDistribution(one) );

			DeterministicNode<double>* mean_squared = new DeterministicNode<double>("mean_squared_mu", new BinaryMultiplication<double,double,double>(mu,mu) );
			DeterministicNode<double>* alpha = new DeterministicNode<double>("mu_alpha", new BinaryDivision<double,double,double>(mean_squared,mu_variance) );
			DeterministicNode<double>* beta = new DeterministicNode<double>("mu_beta", new BinaryDivision<double,double,double>(mu_variance,mu) );

			for(std::map<size_t, size_t>::iterator it = speciesIndex.begin(); it != speciesIndex.end(); it++){
				tmp_name.str("");
				tmp_name << "mu(" << setfill('0') << setw(ws) << it->first << ")";
				mu_vector.push_back(new StochasticNode<double>( tmp_name.str(), new GammaDistribution(alpha,beta) ) );
			}
		}else{
			mu = new StochasticNode<double>("mu", new ExponentialDistribution(one) );
		}
	}else{
		mu = new StochasticNode<double>("mean_mu", new ExponentialDistribution(one) );
		mu_variance = new StochasticNode<double>("mu_variance", new ExponentialDistribution(one) );

		DeterministicNode<double>* mean_squared = new DeterministicNode<double>("mean_squared_mu", new BinaryMultiplication<double,double,double>(mu,mu) );
		DeterministicNode<double>* alpha = new DeterministicNode<double>("mu_alpha", new BinaryDivision<double,double,double>(mean_squared,mu_variance) );
		DeterministicNode<double>* beta = new DeterministicNode<double>("mu_beta", new BinaryDivision<double,double,double>(mu_variance,mu) );

		for (size_t i = 0 ; i < numTrees ; i++ ) {
			tmp_name.str("");
			tmp_name << "mu(" << setfill('0') << setw(ws) << i << ")";
			mu_vector.push_back(new StochasticNode<double>( tmp_name.str(), new GammaDistribution(alpha,beta) ) );
		}
		if(mu_species == HIERARCHICAL){
			mu_x_shape = new StochasticNode<double>("mu_x_shape", new ExponentialDistribution(one) );

			for(std::map<size_t, size_t>::iterator it = speciesIndex.begin(); it != speciesIndex.end(); it++){
				size_t j = it->first;
				tmp_name.str("");
				tmp_name << "mu_x(" << setfill('0') << setw(ws) << j << ")";
				mu_x_vector.push_back(new StochasticNode<double>( tmp_name.str(), new GammaDistribution(mu_x_shape,mu_x_shape) ) );
			}
		}
	}

	// base frequencies prior

	if(modeltype == DOLLO){
		q = new DeterministicNode<RateMatrix>( "q", new FreeBinaryRateMatrixFunction(zero) );

		std::vector<const TypedDagNode<double> *> rf_vec(2);
		rf_vec[1] = zero;
		rf_vec[0] = one;
		rf = new DeterministicNode< std::vector< double > >( "rf", new VectorFunction< double >( rf_vec ) );
	}else{
		if(lambda_prior == HOMOGENEOUS && lambda_species == HOMOGENEOUS && mu_prior == HOMOGENEOUS && mu_species == HOMOGENEOUS){
			DeterministicNode<double>* denom = new DeterministicNode<double>( "denom", new BinaryAddition<double,double,double>(lambda,mu) );
			pi = new DeterministicNode<double>( "pi", new BinaryDivision<double,double,double>(lambda,denom) );
			q = new DeterministicNode<RateMatrix>( "q", new FreeBinaryRateMatrixFunction(pi) );

			DeterministicNode<double>* pi0 = new DeterministicNode<double>( "pi0", new BinarySubtraction<double,double,double>(one,pi) );

			std::vector<const TypedDagNode<double> *> rf_vec(2);
			rf_vec[1] = pi;
			rf_vec[0] = pi0;
			rf = new DeterministicNode< std::vector< double > >( "rf", new VectorFunction< double >( rf_vec ) );
		}else if(lambda_prior == HOMOGENEOUS && mu_prior == HOMOGENEOUS){
			for(size_t i = 0; i < numSpecies; i++){
				TypedDagNode<double>* tmp_lambda;
				if(lambda_species == HOMOGENEOUS)
					tmp_lambda = lambda;
				else
					tmp_lambda = lambda_vector[i];

				TypedDagNode<double>* tmp_mu;
				if(mu_species == HOMOGENEOUS)
					tmp_mu = mu;
				else
					tmp_mu = mu_vector[i];

				tmp_name.str("");
				tmp_name << "denom(" << setfill('0') << setw(ws) << i << ")";
				DeterministicNode<double>* denom = new DeterministicNode<double>( tmp_name.str(), new BinaryAddition<double,double,double>(tmp_lambda,tmp_mu) );

				tmp_name.str("");
				tmp_name << "pi(" << setfill('0') << setw(ws) << i << ")";
				pi_vector.push_back(new DeterministicNode<double>( tmp_name.str(), new BinaryDivision<double,double,double>(tmp_lambda,denom) ));

				tmp_name.str("");
				tmp_name << "pi0(" << setfill('0') << setw(ws) << i << ")";
				DeterministicNode<double>* pi0 = new DeterministicNode<double>( tmp_name.str(), new BinarySubtraction<double,double,double>(one,pi_vector.back()) );

				std::vector<const TypedDagNode<double> *> rf_vec(2);
				rf_vec[1] = pi_vector.back();
				rf_vec[0] = pi0;

				tmp_name.str("");
				tmp_name << "rf(" << setfill('0') << setw(ws) << i << ")";
				rf_vector.push_back(new DeterministicNode< std::vector< double > >( tmp_name.str(), new VectorFunction< double >( rf_vec ) ));

				tmp_name.str("");
				tmp_name << "q(" << setfill('0') << setw(ws) << i << ")";
				q_vector.push_back(new DeterministicNode<RateMatrix>( tmp_name.str(), new FreeBinaryRateMatrixFunction(pi_vector.back()) ));
			}
		}
	}

	std::cerr << "prior okay\n";

	MultiBiphy::initModel();
}

BinaryCharEvoModel<BranchLengthTree>* MultiBiphyLambdaMu::initSubModel( size_t i ) {

	const TypedDagNode<double>* lambda_i;
	if(lambda_prior == HIERARCHICAL){
		lambda_i = lambda_vector[i];
	}else{
		lambda_i = lambda;
	}

	const TypedDagNode<double>* mu_i;
	if(mu_prior == HIERARCHICAL){
		mu_i = mu_vector[i];
	}else{
		mu_i = mu;
	}

	size_t numNodes = trees[i]->getNumberOfNodes();


	tmp_name.str("");
	tmp_name << "psi(" << setfill('0') << setw(wt) << i << ")";
	ConstantNode<BranchLengthTree> *psi = new ConstantNode<BranchLengthTree>( tmp_name.str(), trees[i] );

	BinaryCharEvoModel<BranchLengthTree> *charModel;

	//std::cout << "built model " << i << std::endl;
	// check for species specific rates
	std::vector<TypedDagNode<double>* > scaled_lambda;
	std::vector<TypedDagNode<double>* > scaled_mu;

	if(!(lambda_prior == HOMOGENEOUS && mu_prior == HOMOGENEOUS) && !(lambda_species == HOMOGENEOUS && mu_species == HOMOGENEOUS)){
		for(size_t j = 0; j < numSpecies; j++){
			const TypedDagNode<double>* tmp_lambda;
			if(lambda_prior == HIERARCHICAL && lambda_species == HIERARCHICAL){
				tmp_name.str("");
				tmp_name << "lambda(" << setfill('0') << setw(wt) << i << "," << setfill('0') << setw(ws) << j << ")";
				tmp_lambda = new DeterministicNode<double>( tmp_name.str(), new BinaryMultiplication<double,double,double>(lambda_i,lambda_x_vector[j]) );
			}else if(lambda_species == HIERARCHICAL){
				tmp_lambda = lambda_vector[j];
			}else{
				tmp_lambda = lambda_i;
			}

			const TypedDagNode<double>* tmp_mu;
			if(mu_prior == HIERARCHICAL && mu_species == HIERARCHICAL){
				tmp_name.str("");
				tmp_name << "mu(" << setfill('0') << setw(wt) << i << "," << setfill('0') << setw(ws) << j << ")";
				tmp_mu = new DeterministicNode<double>( tmp_name.str(), new BinaryMultiplication<double,double,double>(mu_i,mu_x_vector[j]) );
			}else if(mu_species == HIERARCHICAL){
				tmp_mu = mu_vector[j];
			}else{
				tmp_mu = mu_i;
			}

			tmp_name.str("");
			tmp_name << "denom(" << setfill('0') << setw(wt) << i << "," << setfill('0') << setw(ws) << j << ")";
			DeterministicNode<double>* denom = new DeterministicNode<double>( tmp_name.str(), new BinaryAddition<double,double,double>(tmp_lambda,tmp_mu) );

			tmp_name.str("");
			tmp_name << "pi(" << setfill('0') << setw(wt) << i << "," << setfill('0') << setw(ws) << j << ")";
			pi_vector.push_back(new DeterministicNode<double>( tmp_name.str(), new BinaryDivision<double,double,double>(tmp_lambda,denom) ));

			tmp_name.str("");
			tmp_name << "pi0(" << setfill('0') << setw(wt) << i << "," << setfill('0') << setw(ws) << j << ")";
			DeterministicNode<double>* pi0 = new DeterministicNode<double>( tmp_name.str(), new BinarySubtraction<double,double,double>(one,pi_vector.back()) );

			std::vector<const TypedDagNode<double> *> rf_vec(2);
			rf_vec[1] = pi_vector.back();
			rf_vec[0] = pi0;

			tmp_name.str("");
			tmp_name << "rf(" << setfill('0') << setw(wt) << i << "," << setfill('0') << setw(ws) << j << ")";
			rf_vector.push_back(new DeterministicNode< std::vector< double > >( tmp_name.str(), new VectorFunction< double >( rf_vec ) ));

			tmp_name.str("");
			tmp_name << "q(" << setfill('0') << setw(wt) << i << "," << setfill('0') << setw(ws) << j << ")";
			q_vector.push_back(new DeterministicNode<RateMatrix>( tmp_name.str(), new FreeBinaryRateMatrixFunction(pi_vector.back()) ));
		}
	}else if(!(lambda_prior == HOMOGENEOUS && mu_prior == HOMOGENEOUS)){
		tmp_name.str("");
		tmp_name << "denom(" << setfill('0') << setw(wt) << i << ")";
		DeterministicNode<double>* denom = new DeterministicNode<double>( tmp_name.str(), new BinaryAddition<double,double,double>(lambda_i,mu_i) );

		tmp_name.str("");
		tmp_name << "pi(" << setfill('0') << setw(wt) << i << ")";
		pi_vector.push_back(new DeterministicNode<double>( tmp_name.str(), new BinaryDivision<double,double,double>(lambda_i,denom) ));

		tmp_name.str("");
		tmp_name << "pi0(" << setfill('0') << setw(wt) << i << ")";
		DeterministicNode<double>* pi0 = new DeterministicNode<double>( tmp_name.str(), new BinarySubtraction<double,double,double>(one,pi_vector.back()) );

		std::vector<const TypedDagNode<double> *> rf_vec(2);
		rf_vec[1] = pi_vector.back();
		rf_vec[0] = pi0;

		rf = new DeterministicNode< std::vector< double > >( tmp_name.str(), new VectorFunction< double >( rf_vec ) );

		tmp_name.str("");
		tmp_name << "q(" << setfill('0') << setw(wt) << i << ")";
		q = new DeterministicNode<RateMatrix>( tmp_name.str(), new FreeBinaryRateMatrixFunction(pi_vector.back()) );
	}

	// modeltype switch
	if(modeltype == REVERSIBLE){
		/* REVERSIBLE MODEL */

		charModel = new BinaryMissingCharEvoModel<BranchLengthTree>(psi, true, data[i]->getNumberOfIncludedCharacters(), correction);
		if(missing)
			((BinaryMissingCharEvoModel<BranchLengthTree>*)charModel)->setMissingRate(xi);

		if(!(lambda_species == HOMOGENEOUS && mu_species == HOMOGENEOUS)){
			std::vector<const TypedDagNode< RateMatrix >* > qs;

			for (size_t j = 0 ; j < numNodes; j++ )
				qs.push_back(q_vector[speciesIndex[node2speciesMaps[i][j]]]);

			tmp_name.str("");
			tmp_name << "q_vector(" << setfill('0') << setw(wt) << i << ")";
			DeterministicNode< RbVector< RateMatrix > >* qs_node = new DeterministicNode< RbVector< RateMatrix > >( tmp_name.str(), new RbVectorFunction<RateMatrix>(qs) );
			charModel->setRateMatrix( qs_node );

			rf = rf_vector[speciesIndex[rootSpecies[i]]];
		}else{
			charModel->setRateMatrix( q );
		}
	}

	charModel->setRootFrequencies( rf );

	if(dgam > 1)
		charModel->setSiteRates( site_rates );

	return charModel;
}

void MultiBiphyLambdaMu::initMCMC( void ) {

	moves.push_back( new ScaleMove(lambda, 1.0, true, 1.0 ) );
	moves.push_back( new ScaleMove(mu, 1.0, true, 1.0 ) );
	monitoredNodes.push_back(lambda);
	monitoredNodes.push_back(mu);

	/* lambda moves */
	if(lambda_prior == HOMOGENEOUS){
		if(lambda_species == HIERARCHICAL){
			moves.push_back( new ScaleMove(lambda_variance, 1.0, true, 1.0 ) );
			monitoredNodes.push_back(lambda_variance);

			for(size_t i = 0; i < lambda_vector.size(); i++){
				moves.push_back( new ScaleMove(lambda_vector[i], 1.0, true, 1.0 ) );
				monitoredNodes.push_back(lambda_vector[i]);
			}
		}
	}else{
		moves.push_back( new ScaleMove(lambda_variance, 1.0, true, 1.0 ) );
		monitoredNodes.push_back(lambda_variance);

		for(size_t i = 0; i < lambda_vector.size(); i++){
			moves.push_back( new ScaleMove(lambda_vector[i], 1.0, true, 1.0 ) );
		}
		if(lambda_species == HIERARCHICAL){
			moves.push_back( new ScaleMove(lambda_x_shape, 1.0, true, 1.0 ) );
			monitoredNodes.push_back(lambda_x_shape);

			for(size_t i = 0; i < lambda_x_vector.size(); i++){
				moves.push_back( new ScaleMove(lambda_x_vector[i], 1.0, true, 1.0 ) );
				monitoredNodes.push_back(lambda_x_vector[i]);
			}
		}
	}

	/* mu moves */
	if(mu_prior == HOMOGENEOUS){
		if(mu_species == HIERARCHICAL){
			moves.push_back( new ScaleMove(mu_variance, 1.0, true, 1.0 ) );
			monitoredNodes.push_back(mu_variance);

			for(size_t i = 0; i < mu_vector.size(); i++){
				moves.push_back( new ScaleMove(mu_vector[i], 1.0, true, 1.0 ) );
				monitoredNodes.push_back(mu_vector[i]);
			}
		}
	}else{
		moves.push_back( new ScaleMove(mu_variance, 1.0, true, 1.0 ) );
		monitoredNodes.push_back(mu_variance);

		for(size_t i = 0; i < mu_vector.size(); i++){
			moves.push_back( new ScaleMove(mu_vector[i], 1.0, true, 1.0 ) );
		}
		if(mu_species == HIERARCHICAL){
			moves.push_back( new ScaleMove(mu_x_shape, 1.0, true, 1.0 ) );
			monitoredNodes.push_back(mu_x_shape);

			for(size_t i = 0; i < mu_x_vector.size(); i++){
				moves.push_back( new ScaleMove(mu_x_vector[i], 1.0, true, 1.0 ) );
				monitoredNodes.push_back(mu_x_vector[i]);
			}
		}
	}

	MultiBiphy::initMCMC();
}

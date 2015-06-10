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

#include <unistd.h>
#include "MultiBiphyLambdaMu.h"

using namespace RevBayesCore;

MultiBiphyLambdaMu::MultiBiphyLambdaMu(const std::string df,
						const std::string n,
						const std::string t,
						ModelType mt,
						RatePrior Ap,
						RatePrior Bp,
						RatePrior As,
						RatePrior Bs,
						int c,
						int d,
						int e,
						int u,
						int num,
						int swa,
						double de,
						double si,
						bool sav,
						bool anc) :
		MultiBiphy(df,
		n,
		t,
		mt,
		Ap,
		Bp,
		As,
		Bs,
		c,
		d,
		e,
		u,
		num,
		swa,
		de,
		si,
		sav,
		anc)
{
	init();
	save();
}

MultiBiphyLambdaMu::MultiBiphyLambdaMu(const std::string n, bool anc) : MultiBiphy(n,anc)
{
	init();
}

MultiBiphyLambdaMu::MultiBiphyLambdaMu(const std::string n) : MultiBiphy(n)
{
	open();
	init();
}

void MultiBiphyLambdaMu::init( void ) {

	RevBayesCore::RandomNumberGenerator* rng = RevBayesCore::GLOBAL_RNG;
	rng->setSeed(std::vector<unsigned int>(2,time(NULL)*getpid()));
    
    /* First, we read in the data */
    // the matrix
	RbSettings::userSettings().setPrintNodeIndex(false);

	std::vector<AbstractCharacterData*> data;
	data = NclReader::getInstance().readMatrices(dataFile);
	if(data.empty()){
		std::cerr << "Error: failed to read datafile " << dataFile << std::endl;
		exit(1);
	}

	bool missing = false;
	size_t excluded = 0;
	size_t numSites = 0;
	for(size_t i = 0; i < data.size(); i++){
		if(data[i]->getDatatype() != "Standard"){
			std::cerr << "Error: incompatible datatype '" << data[i]->getDatatype() << "' for datafile " << dataFile << std::endl;
			exit(1);
		}
		for(size_t c = 0; c < data[i]->getNumberOfCharacters(); c++){
			bool present = false;
			for(size_t t = 0; t < data[i]->getNumberOfTaxa(); t++){
				StandardState state = ((DiscreteCharacterData<StandardState> *)data[i])->getCharacter(t,c);
				missing = missing || state.isGapState();
				present = present || state.getStringValue() == "1";
			}
			if(!present && correction == RevBayesCore::NO_ABSENT_SITES){
				data[i]->excludeCharacter(c);
				excluded++;
			}
		}
		numSites += data[i]->getNumberOfIncludedCharacters();
	}

	if(excluded)
		std::cout << "excluded " << excluded << " sites\n";

	std::cout << "included " << numSites << " sites\n";

	std::vector<BranchLengthTree*> trees;
	NHXTreeReader reader;
	trees = *(reader.readBranchLengthTrees(treeFile));
	if(trees.size() != data.size()){
		std::cerr << "Error: number of matrices (" << data.size() << ") does not match number of trees (" << trees.size() << ")\n";
		exit(1);
	}
	std::cout << "read " << trees.size() << " trees" << std::endl;

	std::vector<size_t> rootSpecies;

	std::map<size_t, size_t> counts;
	std::vector<std::map<size_t, size_t> > node2speciesMaps;
	for(size_t i = 0; i < trees.size(); i++){
		std::vector<TopologyNode*> nodes = trees[i]->getNodes();

		std::map<size_t, size_t> smap;
		for(size_t index = 0; index < nodes.size(); index++){
			size_t sIndex = size_t(nodes[index]->getBranchParameter("S"));

			smap[index] = sIndex;
			counts[sIndex]++;
			if(nodes[index]->isRoot()){
				rootSpecies.push_back(sIndex);
			}
		}

		node2speciesMaps.push_back(smap);
	}
	std::map<size_t, size_t> speciesIndex;

	size_t numSpecies = 0;
	for(std::map<size_t, size_t>::iterator it = counts.begin(); it != counts.end(); it++)
		speciesIndex[it->first] = numSpecies++;

	// speciesIndex[node2speciesMaps[gene][node]] = species parameter index

	std::cout << std::endl;
    if(modeltype == DOLLO){
		std::cout << "dollo model";
	}else{
		std::cout << "reversible model\n";
	}

    std::cout << "\n";

    RatePrior lambda_prior 		= A_prior;
    RatePrior lambda_species 	= A_species;
    RatePrior mu_prior 			= B_prior;
    RatePrior mu_species 		= B_species;

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

    if(dgam > 1){
		std::cout << "\n";

		std::cout << "rates across sites ~ discrete gamma with " << dgam << " rate categories\n";
		std::cout << "\tmean 1 and variance 1/alpha^2\n";
		std::cout << "alpha ~ exponential of mean 1\n";
	}

    if(correction != NONE){
    	std::cout << "\n";
    	std::cout << "correction for unobservable site-patterns:\n";
    	if(correction == NO_UNINFORMATIVE){
    		std::cout << "\tuninformative sites\n";
    	}else{
    		if((correction & NO_CONSTANT_SITES) == NO_CONSTANT_SITES)
				std::cout << "\tconstant sites\n";
			else if(correction & NO_ABSENT_SITES)
					std::cout << "\tconstant absence\n";
			else if(correction & NO_PRESENT_SITES)
					std::cout << "\tconstant presence\n";

			if((correction & NO_SINGLETONS) == NO_SINGLETONS)
				std::cout << "\tsingletons\n";
			else if(correction & NO_SINGLETON_GAINS)
				std::cout << "\tsingleton gains\n";
			else if(correction & NO_SINGLETON_LOSSES)
				std::cout << "\tsingleton losses\n";
    	}
    }

    if(missing)
    	std::cout << "proportion of missing data ~ Beta(1,1)\n";

    std::cout << "\n";

    //////////////////////
    // first the priors //
    //////////////////////
    
    // constant nodes
    ConstantNode<double> *one = new ConstantNode<double>("one", new double(1.0) );
    ConstantNode<double> *zero = new ConstantNode<double>("zero", new double(0.0) );

    StochasticNode<double> *lambda;
    StochasticNode<double> *mu;
    DeterministicNode<double> *	pi;

	std::vector<StochasticNode<double > *> lambda_vector;
	std::vector<StochasticNode<double > *> mu_vector;

	std::vector<StochasticNode<double > *> lambda_x_vector;
	std::vector<StochasticNode<double > *> mu_x_vector;

	std::vector<DeterministicNode<double > *> pi_vector;

	//number of branches
	size_t numTrees = trees.size();
	size_t wt = numTrees > 0 ? (int) log10 ((double) numTrees) + 1 : 1;
	size_t ws = numSpecies > 0 ? (int) log10 ((double) numSpecies) + 1 : 1;

	std::stringstream tmp_name;

	StochasticNode<double> *lambda_variance;
	StochasticNode<double> *mu_variance;

	StochasticNode<double> *lambda_x_shape;
	StochasticNode<double> *mu_x_shape;

	std::vector<const TypedDagNode<int>* > M_vector;

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

	// RAS prior
	ContinuousStochasticNode *shape;
	std::vector<const TypedDagNode<double>* > gamma_rates = std::vector<const TypedDagNode<double>* >();
	DeterministicNode<std::vector<double> > *site_rates;
	DeterministicNode<std::vector<double> > *site_rates_norm;
	if(dgam > 1){
		shape = new ContinuousStochasticNode("shape", new ExponentialDistribution(one) );
		for(size_t cat = 0; cat < dgam; cat++){
			std::stringstream name;
			std::stringstream value_name;
			name << "q";
			name << cat+1;
			value_name << name.str() << "_value";
			gamma_rates.push_back( new DeterministicNode<double>(value_name.str(), new QuantileFunction(new ConstantNode<double>(name.str(), new double((cat+1.0/2.0)/dgam) ), new GammaDistribution(shape, shape) ) ));
		}
		site_rates = new DeterministicNode<std::vector<double> >( "site_rates", new VectorFunction<double>(gamma_rates) );
		//site_rates_norm = new DeterministicNode<std::vector<double> >( "site_rates_norm", new NormalizeVectorFunction(site_rates) );
	}

	//missing data rate prior
	StochasticNode<double> *xi;
	if(missing)
		xi = new StochasticNode<double>( "xi", new BetaDistribution( one,one ) );

	// base frequencies prior
	DeterministicNode< RateMatrix >* q = NULL;
	std::vector<const TypedDagNode< RateMatrix >* > q_vector;

	TypedDagNode<std::vector<double> > *rf;
	std::vector<TypedDagNode<std::vector<double> >* > rf_vector;

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

	std::cout << "prior okay\n";

	// Dollo:		q = fixed, phi = constant
	// RevH:		q = fixed, phi = fixed
	// RevHSpecies:	q_vector, phi_vector = species, phi = root phi

	std::vector<StochasticNode<AbstractCharacterData>* > dataNodes;

	std::cout << "Initializing " << trees.size() << " likelihood functions\n";
	for(size_t i = 0; i < trees.size(); i++){

		try{
		//std::cerr << trees[i]->getTipNames() << std::endl;

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

		TypedDagNode< std::vector< double > >* xis;
		if(missing){
			std::vector<const TypedDagNode<double > *> xi_vector;
			for(size_t k = 0; k < trees[i]->getNumberOfTips(); k++){
				tmp_name.str("");
				tmp_name << "xi(" << setfill('0') << setw(wt) << i << "," << setw(ws) << k << ")";
				xi_vector.push_back(new DeterministicNode<double>( tmp_name.str(), new ConstantFunction< double >( xi ) ));
			}

			tmp_name.str("");
			tmp_name << "xis(" << setfill('0') << setw(wt) << i << ")";
			xis = new DeterministicNode< std::vector< double > >( tmp_name.str(), new VectorFunction<double>( xi_vector ) );
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

			if(missing)
				charModel = new BinaryMissingCharEvoModel<BranchLengthTree>(psi, xis, true, data[i]->getNumberOfIncludedCharacters(), correction);
			else
				charModel = new BinaryCharEvoModel<BranchLengthTree>(psi, true, data[i]->getNumberOfIncludedCharacters(), correction);

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

		tmp_name.str("");
		tmp_name << "S(" << setfill('0') << setw(wt) << i << ")";
		StochasticNode< AbstractCharacterData >* charactermodel = new StochasticNode< AbstractCharacterData >(tmp_name.str(), charModel );
		charactermodel->clamp( data[i] );

		dataNodes.push_back(charactermodel);

		if(correction != NONE){
			tmp_name.str("");
			tmp_name << "N(" << setfill('0') << setw(wt) << i << ")";
			ConstantNode<int>* N = new ConstantNode<int>(tmp_name.str(), new int(data[i]->getNumberOfIncludedCharacters() + 1));

			tmp_name.str("");
			tmp_name << "p(" << setfill('0') << setw(wt) << i << ")";
			DeterministicNode<double>* p = new DeterministicNode<double>(tmp_name.str(), new LnCorrectionFunction<StandardState, BranchLengthTree>(dataNodes.back()) );

			tmp_name.str("");
			tmp_name << "M(" << setfill('0') << setw(wt) << i << ")";
			M_vector.push_back(new StochasticNode<int>(tmp_name.str(), new NegativeBinomialDistribution(N, p) ));
		}

		//std::cerr << charactermodel->getLnProbability() << std::endl;
		}catch(RbException& e){
			std::cerr << "Error building model for " << data[i]->getFileName() << std::endl;
			std::cerr << e.getMessage() << std::endl;
		}
	}
	std::cout << "done\n";
    
    // add the moves
    std::vector<Move*> moves;
	std::vector<Monitor*> monitors;
	std::vector<DagNode*> monitoredNodes;

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

	if(dgam > 1){
		moves.push_back( new ScaleMove(shape, 1.0, true, 1.0) );
		monitoredNodes.push_back( shape );
		//monitoredNodes.push_back(site_rates_norm);
	}

	if(correction != NONE){
		for(size_t i = 0; i < M_vector.size(); i++){
			moves.push_back( new GibbsMove<int>((StochasticNode<int>*)M_vector[i], 1.0 ) );
			//monitoredNodes.push_back((StochasticNode<int>*)M_vector[i] );
		}

		DeterministicNode<std::vector<int> >* M_vector_node = new DeterministicNode< std::vector<int> >( "M_vector", new VectorFunction<int>( M_vector ) );
		DeterministicNode<int>* sum = new DeterministicNode<int>("M", new SumFunction<int>(M_vector_node));

		monitoredNodes.push_back( sum );
	}

	if(missing){
		moves.push_back( new BetaSimplexMove(xi, 1.0, true, 2.0 ) );
		monitoredNodes.push_back(xi);
	}

	Model myModel = Model(one);
	std::cout << "model okay\n";

	std::vector<DagNode*> nodes = myModel.getDagNodes();
	//for(size_t i = 0; i < nodes.size(); i++)
	//	std::cerr << nodes[i]->getName() << std::endl;

	bool useParallelMcmcmc = (numChains > 1);

	if(ancestral){
		for(size_t i = 0; i < dataNodes.size(); i++){
			std::string treelist = name+"."+data[i]->getFileName()+".treelist";
			if(!restart)
				remove(treelist.c_str());
			monitors.push_back( new MappingMonitor<StandardState,BranchLengthTree>( dataNodes[i], every, name+"."+data[i]->getFileName()+".treelist", useParallelMcmcmc || restart) );
		}
	}

	monitors.push_back( new FileMonitor( monitoredNodes, every, name+".trace", "\t", false, true, false, useParallelMcmcmc || restart, useParallelMcmcmc, false ) );

	//std::cout << "init mcmc\n";
	double startingHeat = 1.0;
	mcmc = new ParallelMcmcmc(myModel, moves, monitors, name+".stream", "random", every, numChains, numChains, swapInterval, delta, sigma, startingHeat, saveall);
	std::cout << "mcmc initalized\n";
}

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
#include "MultiBiphyRhoPi.h"

using namespace RevBayesCore;

MultiBiphyRhoPi::MultiBiphyRhoPi(const std::string df,
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

MultiBiphyRhoPi::MultiBiphyRhoPi(const std::string n, bool anc) : MultiBiphy(n,anc)
{
	init();
}

MultiBiphyRhoPi::MultiBiphyRhoPi(const std::string n) : MultiBiphy(n)
{
	open();
	init();
}

void MultiBiphyRhoPi::init( void ) {

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

    RatePrior rho_prior 	= A_prior;
	RatePrior rho_species 	= A_species;
	RatePrior pi_prior 		= B_prior;
	RatePrior pi_species 	= B_species;

    if(pi_prior == HOMOGENEOUS && pi_species == HOMOGENEOUS){
    	std::cout << "pi ~ Beta(1,1)\n";
    }else{
    	std::cout << "alpha ~ exp(1)\n";
		std::cout << "beta ~ exp(1)\n";
		std::cout << "pi_i ~ Beta(alpha,beta)\n\n";

		if(pi_prior == HIERARCHICAL && pi_species == HIERARCHICAL){
			std::cout << "species modulators:\n";
			std::cout << "alpha_x ~ exp(1)\n";
			std::cout << "beta_x ~ exp(1)\n";
			std::cout << "pi_x_i ~ Beta(alpha_x,beta_x)\n\n";
    	}
    }

    std::cout << "\n";

    if(rho_prior == HOMOGENEOUS && rho_species == HOMOGENEOUS){
		std::cout << "rho ~ exp(1)\n";
	}else{
		std::cout << "rho_mean ~ exp(1)\n";
		std::cout << "rho_variance ~ exp(1)\n";
		std::cout << "rho_i ~ gamma with mean = mu_mean and variance = mu_variance\n\n";

		if(rho_prior == HIERARCHICAL && rho_species == HIERARCHICAL){
			std::cout << "species multipliers:\n";
			std::cout << "rho_x_shape ~ exp(1)\n";
			std::cout << "rho_x_i ~ gamma(mu_x_shape, mu_x_shape)\n";
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

	std::vector<StochasticNode<double > *> pi_vector;
	std::vector<DeterministicNode<double > *> pi0_vector;
	std::vector<StochasticNode<double > *> rho_vector;

	std::vector<StochasticNode<double > *> pi_x_vector;
	std::vector<DeterministicNode<double > *> pi0_x_vector;
	std::vector<StochasticNode<double > *> rho_x_vector;

	//number of branches
	size_t numTrees = trees.size();
	size_t wt = numTrees > 0 ? (int) log10 ((double) numTrees) + 1 : 1;
	size_t ws = numSpecies > 0 ? (int) log10 ((double) numSpecies) + 1 : 1;

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

	std::vector<const TypedDagNode<int>* > M_vector;

	/* lambda prior */
	std::string n;
	if(pi_prior == HOMOGENEOUS){
		if(pi_species == HIERARCHICAL){
			alpha = new StochasticNode<double>("alpha", new ExponentialDistribution(one) );
			beta = new StochasticNode<double>("beta", new ExponentialDistribution(one) );

			DeterministicNode<double>* denom = new DeterministicNode<double>("mean_lambda_squared", new BinaryAddition<double,double,double>(alpha,beta));
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

	std::vector<const TypedDagNode<double> *> rf_vec(2);
	TypedDagNode<std::vector<double> > *rf;

	if(modeltype == DOLLO){
		q = new DeterministicNode<RateMatrix>( "q", new FreeBinaryRateMatrixFunction(zero) );

		rf_vec[1] = zero;
		rf_vec[0] = one;
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

	// Dollo:		q = fixed, phi = constant
	// RevH:		q = fixed, phi = fixed
	// RevHSpecies:	q_vector, phi_vector = species, phi = root phi

	std::vector<StochasticNode<AbstractCharacterData>* > dataNodes;

	std::cout << "Initializing " << trees.size() << " likelihood functions\n";
	for(size_t i = 0; i < trees.size(); i++){

		try{
		//std::cerr << trees[i]->getTipNames() << std::endl;

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
				tmp_name.str("");
				tmp_name << "rho(" << setfill('0') << setw(wt) << i << "," << setfill('0') << setw(ws) << j << ")";

				if(rho_prior == HIERARCHICAL)
					scaled_rho.push_back(new DeterministicNode<double>(tmp_name.str(), new BinaryMultiplication<double,double,double>(rho_i,rho_x_vector[j])));
				else
					scaled_rho.push_back(rho_vector[speciesIndex[node2speciesMaps[i][j]]]);
			}
		}

		// modeltype switch
		if(modeltype == REVERSIBLE){
			/* REVERSIBLE MODEL */

			if(missing)
				charModel = new BinaryMissingCharEvoModel<BranchLengthTree>(psi, xis, true, data[i]->getNumberOfIncludedCharacters(), correction);
			else
				charModel = new BinaryCharEvoModel<BranchLengthTree>(psi, true, data[i]->getNumberOfIncludedCharacters(), correction);

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

			std::vector<const TypedDagNode<double> *> scaled_rho;
			if(rho_species == HIERARCHICAL){
				tmp_name.str("");
				tmp_name << "rhos(" << setfill('0') << setw(wt) << i << ")";
				TypedDagNode<std::vector<double> >* rho_clock = new DeterministicNode< std::vector< double > >( tmp_name.str(), new VectorFunction<double>( scaled_rho ) );

				charModel->setClockRate(rho_clock);
			}else{
				charModel->setClockRate(rho_i);
			}
		}else{
			/* DOLLO MODEL */

			tmp_name.str("");
			tmp_name << "tau(" << setfill('0') << setw(wt) << i << ")";
			ConstantNode<Topology> *tau = new ConstantNode<Topology>( tmp_name.str(), (Topology*)&(psi->getValue().getTopology()) );

			if(missing)
				charModel = new DolloBinaryMissingCharEvoModel<BranchLengthTree>(psi, tau, xis, true, data[i]->getNumberOfIncludedCharacters(), correction);
			else
				charModel = new DolloBinaryCharEvoModel<BranchLengthTree>(psi, tau, true, data[i]->getNumberOfIncludedCharacters(), correction);

			/* set clock rates */

			charModel->setRateMatrix( q );
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
			DeterministicNode<double>* p = new DeterministicNode<double>(tmp_name.str(), new LnCorrectionFunction<StandardState, BranchLengthTree>(charactermodel));

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

	if(dgam > 1){
		moves.push_back( new ScaleMove(shape, 1.0, true, 1.0) );
		monitoredNodes.push_back( shape );
		//monitoredNodes.push_back(site_rates_norm);
	}

	if(correction != NONE){
		for(size_t i = 0; i < M_vector.size(); i++){
			moves.push_back( new GibbsMove<int>((StochasticNode<int>*)M_vector[i], 1.0 ) );
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
	//	std::cerr << nodes[i]->getLnProbability() << std::endl;

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

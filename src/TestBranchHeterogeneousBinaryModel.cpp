#include "BinaryDivision.h"
#include "BinaryMultiplication.h"
#include "BinarySubtraction.h"
#include "BetaDistribution.h"
#include "BetaSimplexMove.h"
#include "Clade.h"
#include "CladeReader.h"
#include "ConstantNode.h"
#include "CrossValidationScoreMonitor.h"
#include "DeterministicNode.h"
#include "DirichletDistribution.h"
#include "DirichletProcessPriorDistribution.h"
#include "DolloBranchHeterogeneousCharEvoModel.h"
#include "DPPAllocateAuxGibbsMove.h"
#include "DPPGibbsConcentrationMove.h"
#include "DppNumTablesStatistic.h"
#include "DPPBetaSimplexMove.h"
#include "ExponentialDistribution.h"
#include "ExtendedNewickTreeMonitor.h"
#include "FileMonitor.h"
#include "FreeBinaryRateMatrixFunction.h"
#include "FreeBinaryRateMatrixVectorFunction.h"
#include "GammaDistribution.h"
#include "GeneralBranchHeterogeneousCharEvoModel.h"
#include "ParallelMcmcmc.h"
#include "MixtureDistribution.h"
#include "MixtureAllocationMove.h"
#include "Mcmc.h"
#include "MeanFunction.h"
#include "Model.h"
#include "ModelStreamMonitor.h"
#include "ModelStreamReader.h"
#include "Monitor.h"
#include "Move.h"
#include "MultiMove.h"
#include "NclReader.h"
#include "NewickTreeMonitor.h"
#include "NexusTreeMonitor.h"
#include "PosteriorPredictiveStateFrequencyMonitor.h"
#include "QuantileFunction.h"
#include "RbSettings.h"
#include "RbStatisticsHelper.h"
#include "RbVectorFunction.h"
#include "ScaleMove.h"
#include "ScreenMonitor.h"
#include "SimplexSingleElementScale.h"
#include "SlidingMove.h"
#include "SubtreePruneRegraft.h"
#include "TestBranchHeterogeneousBinaryModel.h"
#include "TimeTree.h"
#include "TreeAssemblyFunction.h"
#include "TreeLengthStatistic.h"
#include "TreeScale.h"
#include "TruncatedDistributionUnnormalized.h"
#include "UniformDistribution.h"
#include "UniformConstrainedTimeTreeDistribution.h"
#include "UniformRootedTopologyDistribution.h"
#include "UniformTimeTreeDistribution.h"
#include "VectorFunction.h"
#include "VectorScaleFunction.h"


using namespace RevBayesCore;

TestBranchHeterogeneousBinaryModel::TestBranchHeterogeneousBinaryModel(const std::string &datafile,
									const std::string &name,
									const std::string &treefile,
									const std::string &outgroupfile,
									int branchprior,
									bool ras,
									int heterogeneous,
									bool dollo,
									int mixture,
									bool rootprior,
									double rootmin,
									double rootmax,
									int every,
									int until,
									int numchains,
									int swapInterval,
									double deltaTemp,
									double sigmaTemp,
									bool saveall,
									bool nexus) :
		dataFile( datafile ),
		name( name ),
		treeFile( treefile ),
		outgroupFile( outgroupfile ),
		cvfile("None"),
		branchprior(branchprior),
		ras( ras),
		heterogeneous( heterogeneous ),
		dollo(dollo),
		mixture( mixture),
		ppred(false),
		rootprior( rootprior),
		rootmin( rootmin),
		rootmax( rootmax),
		every( every ),
		until( until ),
		numChains( numchains),
		swapInterval( swapInterval ),
		deltaTemp( deltaTemp ),
		sigmaTemp( sigmaTemp ),
		saveall( saveall),
		nexus(nexus)
{
    save();
    readstream = false;
}

TestBranchHeterogeneousBinaryModel::TestBranchHeterogeneousBinaryModel(const std::string &name, const std::string &cvfile, bool ppred) :
		name( name ),
		cvfile(cvfile),
		ppred(ppred)
{
    open();
    readstream = true;
    every = 1;
}

TestBranchHeterogeneousBinaryModel::~TestBranchHeterogeneousBinaryModel() {
    // nothing to do
}

void TestBranchHeterogeneousBinaryModel::open( void ) {
	std::ifstream is((name + ".param").c_str());
	if (!is)        {
		std::cerr << "error : cannot find file : " << name << ".param\n";
		exit(1);
	}

	is >> dataFile;
	is >> name;
	is >> treeFile;
	is >> outgroupFile;
	is >> heterogeneous;
	is >> dollo;
	is >> mixture;
	is >> ras;
	is >> rootprior;
	is >> rootmin;
	is >> rootmax;
	is >> every;
	is >> until;
	is >> numChains;
	is >> swapInterval;
	is >> deltaTemp;
	is >> sigmaTemp;
	is >> saveall;
	is >> nexus;

	is.close();
}

void TestBranchHeterogeneousBinaryModel::save( void ) {
	std::ofstream os((name + ".param").c_str());

	os << dataFile << "\n";
	os << name << "\n";
	os << treeFile << "\n";
	os << outgroupFile << "\n";
	os << heterogeneous << "\n";
	os << dollo << "\n";
	os << mixture << "\n";
	os << ras << "\n";
	os << rootprior << "\n";
	os << rootmin << "\n";
	os << rootmax << "\n";
	os << every << "\n";
	os << until << "\n";
	os << numChains << "\n";
	os << swapInterval << "\n";
	os << deltaTemp << "\n";
	os << sigmaTemp << "\n";
	os << saveall << "\n";
	os << nexus << "\n";

	os.close();
}



bool TestBranchHeterogeneousBinaryModel::run( void ) {
    
    /* First, we read in the data */
    // the matrix
	RbSettings::userSettings().setPrintNodeIndex(false);

    std::vector<AbstractCharacterData*> data = NclReader::getInstance().readMatrices(dataFile);
    std::cout << "Read " << data.size() << " matrices." << std::endl;

    std::vector<AbstractCharacterData*> cvdata;
    if(cvfile != "None"){
		cvdata = NclReader::getInstance().readMatrices(cvfile);
		std::cout << "Cross-validating " << cvdata.size() << " matrices." << std::endl;
    }

    std::vector<TimeTree*> trees;
    Topology fixedTree;
    bool treeFixed = false;
    if(treeFile != "None"){
    	trees = NclReader::getInstance().readTimeTrees(treeFile);
    	treeFixed = trees.size() > 0;
    	if(trees.size() > 0){
    		fixedTree = trees[0]->getTopology();
    		std::cout << "Using fixed tree:\n";
    		std::cout << fixedTree.getNewickRepresentation() << std::endl;
    		std::cout << std::endl;
    	}
    }
   // std::cout << trees[0]->getNewickRepresentation() << std::endl;
    
    //construct constraints
    Clade outgroup(std::vector<std::string>(),0);
    if(outgroupFile != "None"){
    	CladeReader reader(outgroupFile);
    	outgroup = reader.getClades()[0];
    	std::cout << "Outgroup:\n";
    	std::cout << outgroup << std::endl;
    	std::cout << std::endl;
    }
    
    if(heterogeneous){
    	std::cout << "branch-wise constant time-heterogeneous model:\n";
    	switch(heterogeneous){
			case 2:
				std::cout << "branch frequencies = pi(i) ~ DPP(G,cp)\n";
				std::cout << "cp ~ Gamma(a,b)\n";
				std::cout << "a,b ~ iid exponential of mean 1";
				std::cout << "G = Beta(alpha,beta)\n";
				break;
			case 3:
				std::cout << "branch frequencies = pi(i) ~ mixture of " << mixture << " Beta(alpha,beta) distributions\n";
				break;
			default:
				std::cout << "branch frequencies = pi(i) ~ iid Beta(alpha,beta)\n";
				break;
    	}
	}else if(dollo){
		std::cout << "dollo model\n";
	}else{
		std::cout << "time-homogeneous model:\n";
	}
    if(heterogeneous){
    	if(rootprior){
			std::cout << "root frequency  = phi ~ Beta(alpha,beta) truncated on (" << rootmin << "," << rootmax << ")\n";
		}else{
			switch(heterogeneous){
				case 2:
					std::cout << "root frequency = phi ~ DPP(G,cp)\n";
					break;
				case 3:
					std::cout << "root frequency = phi ~ mixture of " << mixture << " Beta(alpha,beta) distributions\n";
					break;
				default:
					std::cout << "root frequency = phi ~ Beta(alpha,beta)\n";
					break;
			}
			std::cout << "alpha, beta ~ iid exponential of mean 1\n";
		}
	}else if(!dollo){
		if(rootprior){
			std::cout << "root frequency  = phi ~ Uniform(" << rootmin << "," << rootmax << ")\n";
		}else{
			std::cout << "root frequency  = phi ~ Uniform(0,1)\n";
		}
	}

    std::cout << "\n";

    if(branchprior == 0){
		std::cout << "branch lengths ~ exponential of mean mu\n";
		std::cout << "mu ~ exponential of mean 0.1\n";
    }else if(branchprior == 1){
    	std::cout << "branch lengths ~ dirichlet with concentration parameter = 1\n";
    	std::cout << "tree length ~ exponential of mean 1\n";
    }

    if(ras){
		std::cout << "\n";

		std::cout << "ras ~ gamma of mean 1 and variance 1/lambda^2\n";
		std::cout << "lambda ~ exponential of mean 1\n";
	}

    std::cout << "\n";

    //////////////////////
    // first the priors //
    //////////////////////
    
    // constant nodes
    ConstantNode<int> *idx = new ConstantNode<int>("idx", new int(1) );
    ConstantNode<double> *one = new ConstantNode<double>("one", new double(1.0) );
    ConstantNode<double> *ten = new ConstantNode<double>("ten", new double(10.0) );
    ConstantNode<double> *zero = new ConstantNode<double>("zero", new double(0.0) );

    //number of branches
	size_t numBranches = 2*data[0]->getNumberOfTaxa() - 2;

	//mixture vectors
	std::vector<const TypedDagNode<double > *> pi_cats;
	TypedDagNode< std::vector< double > >* pi_mix;
	std::vector<const TypedDagNode<double> *> pi_stat;

	//dpp priors
	ContinuousStochasticNode *dpA;
	ContinuousStochasticNode *dpB;
	StochasticNode<double> *cp;
	DeterministicNode<int> *numCats;

	TypedDagNode< std::vector< double > >* pi_vector;

	std::vector<const TypedDagNode< RateMatrix >* > qs;
	DeterministicNode< RateMatrix >* q = NULL;
	DeterministicNode< RbVector< RateMatrix > >* qs_node = NULL;

	// declaring a vector of clock rates
	std::vector<const TypedDagNode<double> *> branchRates;
	std::vector< ContinuousStochasticNode *> branchRates_nonConst;

    // base frequencies hyperprior
    StochasticNode<double> *alpha;
    StochasticNode<double> *beta;

    // base frequencies prior
    TypedDagNode<double> *phi;
    if(dollo){
    	phi = zero;
    }else{
		if(rootprior){
			if(heterogeneous){
				alpha = new StochasticNode<double>( "alpha", new ExponentialDistribution(one) );
				beta = new StochasticNode<double>( "beta", new ExponentialDistribution(one) );
				phi = new StochasticNode<double>( "phi", new TruncatedDistributionUnnormalized( new BetaDistribution( alpha, beta ), new ConstantNode<double>("rootmin", new double(rootmin) ), new ConstantNode<double>("rootmax", new double(rootmax) ) ) );
			}else{
				phi = new StochasticNode<double>( "phi", new UniformDistribution( new ConstantNode<double>("rootmin", new double(rootmin) ), new ConstantNode<double>("rootmax", new double(rootmax) ) ) );
			}
		}else{
			if(heterogeneous){
				alpha = new StochasticNode<double>( "alpha", new ExponentialDistribution(one) );
				beta = new StochasticNode<double>( "beta", new ExponentialDistribution(one) );
				phi = new StochasticNode<double>( "phi", new BetaDistribution( alpha,beta ) );
			}else{
				phi = new StochasticNode<double>( "phi", new BetaDistribution( one,one ) );
			}
		}
    }

    // branch frequency prior
	if(heterogeneous){
		if(heterogeneous == 1){
			for (unsigned int i = 0 ; i < numBranches ; i++ ) {
				std::ostringstream pi_name;
				pi_name << "pi(" << i << ")";
				pi_stat.push_back(new StochasticNode<double>( pi_name.str(), new BetaDistribution(alpha,beta) ) );

				std::ostringstream q_name;
				q_name << "q(" << i << ")";
				qs.push_back(new DeterministicNode<RateMatrix>( q_name.str(), new FreeBinaryRateMatrixFunction(pi_stat[i]) ));
			}
			pi_vector = new DeterministicNode< std::vector< double > >( "pi", new VectorFunction<double>( pi_stat ) );
			qs_node = new DeterministicNode< RbVector< RateMatrix > >( "q_vector", new RbVectorFunction<RateMatrix>(qs) );
		}else if(heterogeneous == 2){
			// Setting up the hyper prior on the concentration parameter
			// This hyperprior is fully conditional on the DPP using a gamma distribution
			//    the parameters of the gamma distribution are set so that the expectation of the hyperprior = meanCP
			//    where meanCP = dpA / dpB
			dpA  = new ContinuousStochasticNode("dp_a", new ExponentialDistribution(one));
			dpB  = new ContinuousStochasticNode("dp_b", new ExponentialDistribution(one));
			cp = new StochasticNode<double>("dpp.cp", new GammaDistribution(dpA, dpB) );

			// G_0 is an Beta distribution
			TypedDistribution<double> *g = new BetaDistribution(alpha,beta);

			// branchRates ~ DPP(g, cp, numBranches)
			pi_vector = new StochasticNode<std::vector<double> >("pi", new DirichletProcessPriorDistribution<double>(g, cp, (int)numBranches + 1) );
			// a deterministic node for calculating the number of rate categories (required for the Gibbs move on cp)
			numCats = new DeterministicNode<int>("dpp.k", new DppNumTablesStatistic<double>((StochasticNode<std::vector<double> >*)pi_vector) );

			delete phi;
			phi = new DeterministicNode<double>("phi", new VectorIndexOperator<double>((StochasticNode<std::vector<double> >*)pi_vector,idx));

			qs_node = new DeterministicNode< RbVector< RateMatrix > >( "q_vector", new FreeBinaryRateMatrixVectorFunction(pi_vector,true) );
		}else if(heterogeneous == 3){
			ConstantNode<std::vector<double> > *probs = new ConstantNode<std::vector<double> >("probs", new std::vector<double>(mixture,1.0) );
			//probs = new DeterministicNode<std::vector<double> >( "probs", new DirichletDistribution(conc));
			for(size_t i = 0; i < mixture; i++){
				std::ostringstream pi_name;
				pi_name << "pi(" << i << ")";
				pi_cats.push_back(new StochasticNode<double >( pi_name.str(), new BetaDistribution(alpha,beta) ) );
			}
			pi_mix = new DeterministicNode< std::vector< double > >( "pi_mix", new VectorFunction< double >( pi_cats ) );
			for (unsigned int i = 0 ; i < numBranches ; i++ ) {
				std::ostringstream q_name;
				q_name << "q(" << i << ")";
				pi_stat.push_back(new StochasticNode<double>( q_name.str()+"_pi", new MixtureDistribution<double>(pi_mix,probs) ));
				qs.push_back(new DeterministicNode<RateMatrix>( q_name.str(), new FreeBinaryRateMatrixFunction(pi_stat[i]) ));
			}

			delete phi;
			phi = new StochasticNode<double>( "phi", new MixtureDistribution<double>(pi_mix,probs) );

			pi_vector = new DeterministicNode< std::vector< double > >( "pi", new VectorFunction<double>( pi_stat ) );
			qs_node = new DeterministicNode< RbVector< RateMatrix > >( "q_vector", new RbVectorFunction<RateMatrix>(qs) );
		}

	}else{
		q = new DeterministicNode<RateMatrix>( "q", new FreeBinaryRateMatrixFunction(phi) );
	}

	// branch length hyperprior
	ContinuousStochasticNode *mu;
	DeterministicNode<double> *rec_mu;
	TypedDagNode< std::vector< double > >* br_vector;

	StochasticNode< std::vector< double > >* br_times;
	DeterministicNode<double>* tree_length_alpha;
	TypedDagNode<double>* tree_length;

	// branch length prior
	if(branchprior == 0){
		mu = new ContinuousStochasticNode("mu", new ExponentialDistribution(ten) );
		rec_mu = new DeterministicNode<double>( "rec_mu", new BinaryDivision<double,double,double>(one,mu));
		int w = numBranches > 0 ? (int) log10 ((double) numBranches) + 1 : 1;
		for (unsigned int i = 0 ; i < numBranches ; i++ ) {
			std::ostringstream br_name;
			br_name << "br(" << setfill('0') << setw(w) << i << ")";

			ContinuousStochasticNode* tmp_branch_rate = new ContinuousStochasticNode( br_name.str(), new ExponentialDistribution(rec_mu));
			branchRates.push_back( tmp_branch_rate );
			branchRates_nonConst.push_back( tmp_branch_rate );
		}
		br_vector = new DeterministicNode< std::vector< double > >( "br_vector", new VectorFunction< double >( branchRates ) );
	}else if(branchprior == 1){
		ConstantNode<std::vector<double> > *conc = new ConstantNode<std::vector<double> >("conc", new std::vector<double>(numBranches,1.0) );
		tree_length = new StochasticNode<double>("length", new GammaDistribution(one,one));
		br_times = new StochasticNode< std::vector< double > >( "br_times", new DirichletDistribution( conc));
		br_vector = new DeterministicNode< std::vector< double > >( "br_vector", new VectorScaleFunction(br_times,tree_length) );
	}

	// RAS prior
    ContinuousStochasticNode *lambda;
    std::vector<const TypedDagNode<double>* > gamma_rates = std::vector<const TypedDagNode<double>* >();
    DeterministicNode<std::vector<double> > *site_rates;
    DeterministicNode<std::vector<double> > *site_rates_norm;
    if(ras){
		lambda = new ContinuousStochasticNode("lambda", new ExponentialDistribution(one) );
		gamma_rates.push_back( new DeterministicNode<double>("q4_value", new QuantileFunction(new ConstantNode<double>("q1", new double(0.125) ), new GammaDistribution(lambda, lambda) ) ));
		gamma_rates.push_back( new DeterministicNode<double>("q4_value", new QuantileFunction(new ConstantNode<double>("q2", new double(0.375) ), new GammaDistribution(lambda, lambda) ) ));
		gamma_rates.push_back( new DeterministicNode<double>("q4_value", new QuantileFunction(new ConstantNode<double>("q3", new double(0.625) ), new GammaDistribution(lambda, lambda) ) ));
		gamma_rates.push_back( new DeterministicNode<double>("q4_value", new QuantileFunction(new ConstantNode<double>("q4", new double(0.875) ), new GammaDistribution(lambda, lambda) ) ));
		site_rates = new DeterministicNode<std::vector<double> >( "site_rates", new VectorFunction<double>(gamma_rates) );
		site_rates_norm = new DeterministicNode<std::vector<double> >( "site_rates_norm", new NormalizeVectorFunction(site_rates) );
    }
            
    // tree prior
    std::vector<std::string> names = data[0]->getTaxonNames();
    //StochasticNode<TimeTree> *tau = new StochasticNode<TimeTree>( "tau", new UniformConstrainedTimeTreeDistribution(one,names,constraints,outgroup) );
    StochasticNode<Topology> *tau = new StochasticNode<Topology>( "tau", new UniformRootedTopologyDistribution(names, std::vector<Clade>(), outgroup) );
    if(treeFixed){
    	tau->setValue(fixedTree, true);
    	tau->setIgnoreRedraw(true);
    }
    
    DeterministicNode<BranchLengthTree> *psi = new DeterministicNode<BranchLengthTree>( "psi", new TreeAssemblyFunction(tau, br_vector) );
    std::cout << "tree okay\n";

    AbstractSiteHomogeneousMixtureCharEvoModel<StandardState, BranchLengthTree> *charModel;
    DolloBranchHeterogeneousCharEvoModel<StandardState, BranchLengthTree> *DcharModel;
    GeneralBranchHeterogeneousCharEvoModel<StandardState, BranchLengthTree> *GcharModel;
    if(dollo){
    	StandardState absorbingstate("01");
		absorbingstate.setState("0");
		DcharModel = new DolloBranchHeterogeneousCharEvoModel<StandardState, BranchLengthTree>(psi, 2, true, data[0]->getNumberOfCharacters(), absorbingstate);

    }else{
    	GcharModel = new GeneralBranchHeterogeneousCharEvoModel<StandardState, BranchLengthTree>(psi, 2, true, data[0]->getNumberOfCharacters());
    }


    std::vector<const TypedDagNode<double> *> rf_vec(2);
	TypedDagNode<std::vector<double> > *rf;

    if(dollo){
    	rf_vec[1] = one;
    	rf_vec[0] = zero;
    	rf = new DeterministicNode< std::vector< double > >( "rf", new VectorFunction< double >( rf_vec ) );
		if(heterogeneous){
			DcharModel->setRateMatrix( qs_node );
			DcharModel->setRootFrequencies( rf );
		}else{
			DcharModel->setRateMatrix( q );
			DcharModel->setRootFrequencies( rf );
		}
		DcharModel->setClockRate( one );
		if(ras)
			DcharModel->setSiteRates( site_rates_norm );

		charModel = DcharModel;
    }else{
    	rf_vec[1] = phi;
    	rf_vec[0] = new DeterministicNode<double>( "pi0", new BinarySubtraction<double,double,double>(one,phi));
    	rf = new DeterministicNode< std::vector< double > >( "rf", new VectorFunction< double >( rf_vec ) );
    	if(heterogeneous){
			GcharModel->setRateMatrix( qs_node );
			GcharModel->setRootFrequencies( rf );
		}else{
			GcharModel->setRateMatrix( q );
			GcharModel->setRootFrequencies( rf );
		}
		GcharModel->setClockRate( one );
		if(ras)
			GcharModel->setSiteRates( site_rates_norm );

		charModel = GcharModel;
    }

    StochasticNode< AbstractCharacterData > *charactermodel = new StochasticNode< AbstractCharacterData >("S", charModel );
    charactermodel->clamp( data[0] );
    std::cout << "data ok\n";

    
    /* add the moves */
    std::vector<Move*> moves;
	std::vector<Monitor*> monitors;
	std::vector<DagNode*> monitoredNodes;

	if(readstream){
		if(ppred)
			monitors.push_back( new PosteriorPredictiveStateFrequencyMonitor( charactermodel, every, name+".ppred") );
		if(cvdata.size() > 0)
			monitors.push_back( new CrossValidationScoreMonitor( charactermodel, cvdata[0], every, name+".cv") );

		//monitoredNodes.push_back(pi_vector);
		//monitors.push_back( new FileMonitor( monitoredNodes, every, name+".out", "\t", false, true, false, false, false, false ) );

		Model myModel = Model(charactermodel);
		std::cout << "model okay\n";

		ModelStreamReader stream(myModel,monitors,every,name+".chain");
		std::cout << "reading from stream\n";
		stream.run(until);
	}else{
		//moves.push_back( new NearestNeighborInterchange( tau, 2.0 ) );
		if(!treeFixed)
			moves.push_back( new SubtreePruneRegraft( tau, 5.0 ) );

		for (size_t i = 0 ; i < numBranches ; i ++ ) {
			if(heterogeneous == 1){
				moves.push_back( new BetaSimplexMove((StochasticNode<double>*)pi_stat[i], 1.0, true, 2.0 ) );
			}else if(heterogeneous == 3){
				moves.push_back( new MixtureAllocationMove<double>((StochasticNode<double>*)pi_stat[i], 2.0 ) );
			}
			if(branchprior == 0)
				moves.push_back( new ScaleMove(branchRates_nonConst[i], 2.0, true, 1.0 ) );
		}

		if(branchprior == 0){
			tree_length = new DeterministicNode<double >("length", new TreeLengthStatistic<BranchLengthTree>(psi) );
		}else if(branchprior == 1){
			moves.push_back(new ScaleMove((StochasticNode<double>*)tree_length, 1.0, true, 1.0));
			moves.push_back(new SimplexSingleElementScale(br_times, 2.0, true, (int)numBranches/3));
		}
		monitoredNodes.push_back( tree_length );

		if(heterogeneous == 3){
			moves.push_back( new MixtureAllocationMove<double>((StochasticNode<double>*)phi, 2.0 ) );
		}else if(heterogeneous != 2 && !dollo){
			moves.push_back( new BetaSimplexMove((StochasticNode<double>*)phi, 1.0, true, 5.0 ) );
		}
		if(!dollo)
			monitoredNodes.push_back( phi );

		DeterministicNode<double>* meanpi;
		if(heterogeneous){
			meanpi = new DeterministicNode<double >( "mean_pi", new MeanFunction(pi_vector) );
			monitoredNodes.push_back( meanpi );
		}

		if(branchprior==0){
			moves.push_back( new ScaleMove(mu, 1.0, true, 1.0) );
			monitoredNodes.push_back( mu );
		}

		if(ras){
			moves.push_back( new ScaleMove(lambda, 1.0, true, 1.0) );
			monitoredNodes.push_back( lambda );
			//monitoredNodes.push_back(site_rates_norm);
		}

		if(heterogeneous){
			moves.push_back( new ScaleMove(alpha, 1.0, true, 1.0) );
			monitoredNodes.push_back( alpha );

			moves.push_back( new ScaleMove(beta, 1.0, true, 1.0) );
			monitoredNodes.push_back( beta );
		}

		if(heterogeneous == 2){
			monitoredNodes.push_back(numCats);
			monitoredNodes.push_back(cp);

			moves.push_back( new ScaleMove(dpA, 1.0, true, 1.0) );
			monitoredNodes.push_back(dpA);

			moves.push_back( new ScaleMove(dpB, 1.0, true, 1.0) );
			monitoredNodes.push_back(dpB);

			moves.push_back( new DPPBetaSimplexMove( (StochasticNode<std::vector<double> >*)pi_vector, 1.0 , 1.0 ) );
			moves.push_back( new DPPAllocateAuxGibbsMove<double>( (StochasticNode<std::vector<double> >*)pi_vector, 4, 2.0 ) );
			moves.push_back( new DPPGibbsConcentrationMove<double>( cp, numCats, dpA, dpB, (int)numBranches + 1, 2.0 ) );
		}else if(heterogeneous == 3){
			std::vector<Move*> mixmoves;
			for (size_t i = 0 ; i < mixture; i ++ ) {
				//mixmoves.push_back( new ScaleMove((StochasticNode<double>*)pi_cats[i], 1.0, true, 2.0 ) );
				moves.push_back( new BetaSimplexMove((StochasticNode<double>*)pi_cats[i], 1.0, true, 2.0 ) );
				monitoredNodes.push_back((StochasticNode<double>*)pi_cats[i]);
			}
			//moves.push_back(new MultiMove(mixmoves,one,2.0,true));
		}

		//monitoredNodes.push_back(pi_vector);

		bool useParallelMcmcmc = (numChains > 1);
		monitors.push_back( new FileMonitor( monitoredNodes, every, name+".trace", "\t", false, true, false, useParallelMcmcmc, useParallelMcmcmc, false ) );
		if(!treeFixed){
			monitors.push_back( new NewickTreeMonitor( psi, every, name+".treelist", useParallelMcmcmc) );
			if(nexus){
				std::set<TypedDagNode<std::vector<double> > *> piset;
				if(heterogeneous)
					piset.insert(pi_vector);
				monitors.push_back( new NexusTreeMonitor( psi, piset, every, name+".treelist.nex", useParallelMcmcmc) );
				//monitors.push_back( new ExtendedNewickTreeMonitor( tau, piset, every, name+".treelist", "\t", false, false, false, useParallelMcmcmc) );
			}
		}

		Model myModel = Model(charactermodel);
		std::cout << "model okay\n";

		if(saveall)
			monitors.push_back( new ModelStreamMonitor( myModel, every, name+".chain", (numChains > 1)) );

		std::cerr << charactermodel->getLnProbability() << std::endl;

		int numProcesses = numChains;
		double startingHeat = 1.0;
		ParallelMcmcmc myPmc3(myModel, moves, monitors, "random", numChains, numProcesses, swapInterval, deltaTemp, sigmaTemp, startingHeat);
		std::cout << "running\n";
		myPmc3.run(until);
	}

    return true;
}

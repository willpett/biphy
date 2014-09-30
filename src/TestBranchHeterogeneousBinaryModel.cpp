#include "BinaryDivision.h"
#include "BinaryMultiplication.h"
#include "BinarySubtraction.h"
#include "BetaDistribution.h"
#include "Clade.h"
#include "CladeReader.h"
#include "ConstantNode.h"
#include "CrossValidationScoreMonitor.h"
#include "DeterministicNode.h"
#include "DirichletDistribution.h"
#include "DirichletProcessPriorDistribution.h"
#include "DPPAllocateAuxGibbsMove.h"
#include "DPPGibbsConcentrationMove.h"
#include "DppNumTablesStatistic.h"
#include "DPPScaleCatValsMove.h"
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
#include "Model.h"
#include "ModelMonitor.h"
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
#include "MeanFunction.h"

using namespace RevBayesCore;

TestBranchHeterogeneousBinaryModel::TestBranchHeterogeneousBinaryModel(const std::string &datafile,
									const std::string &name,
									const std::string &treefile,
									const std::string &outgroupfile,
									const std::string &cvfile,
									bool heterogeneous,
									int mixture,
									bool dpp,
									bool ppred,
									bool rootprior,
									double rootmin,
									double rootmax,
									int every,
									int until,
									int numchains,
									int swapInterval,
									double deltaTemp,
									double sigmaTemp) :
		dataFile( datafile ),
		name( name ),
		treeFile( treefile ),
		outgroupFile( outgroupfile ),
		cvfile( cvfile ),
		heterogeneous( heterogeneous ),
		mixture( mixture),
		dpp( dpp),
		ppred( ppred ),
		rootprior( rootprior),
		rootmin( rootmin),
		rootmax( rootmax),
		every( every ),
		until( until ),
		numChains( numchains),
		swapInterval( swapInterval ),
		deltaTemp( deltaTemp ),
		sigmaTemp( sigmaTemp )
{
    
}

TestBranchHeterogeneousBinaryModel::~TestBranchHeterogeneousBinaryModel() {
    // nothing to do
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
    	if(dpp){
    		std::cout << "branch frequencies = pi(i) ~ DPP(cp)\n";
    	}else if(mixture > 1){
    		std::cout << "branch frequencies = pi(i) ~ mixture of " << mixture << " Beta(alpha,beta) distributions\n";
    	}else{
    		std::cout << "branch frequencies = pi(i) ~ iid Beta(alpha,beta)\n";
    	}
	}else{
		std::cout << "time-homogeneous model:\n";
	}
    if(rootprior){
		std::cout << "root frequency  = phi ~ Beta(alpha,beta) truncated on (" << rootmin << "," << rootmax << ")\n";
	}else if(dpp){
		std::cout << "root frequency = phi ~ DPP(cp)\n";
	}else if(mixture > 1){
		std::cout << "root frequency = phi ~ mixture of " << mixture << " Beta(alpha,beta) distributions\n";
	}else{
		std::cout << "root frequency = phi ~ Beta(alpha,beta)\n";
	}
    std::cout << "alpha, beta ~ iid exponential of mean 1\n";

    std::cout << "\n";

    std::cout << "branch lengths ~ exponential of mean mu\n";
    std::cout << "mu ~ exponential of mean 0.1\n";

    std::cout << "\n";

    std::cout << "ras ~ gamma of mean 1 and variance 1/lambda^2\n";
    std::cout << "lambda ~ exponential of mean 1\n";

    std::cout << "\n";

    //////////////////////
    // first the priors //
    //////////////////////
    
    // constant nodes
    ConstantNode<int> *idx = new ConstantNode<int>("idx", new int(1) );
    ConstantNode<double> *one = new ConstantNode<double>("one", new double(1.0) );
    ConstantNode<double> *ten = new ConstantNode<double>("ten", new double(10.0) );

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
    StochasticNode<double> *alpha = new StochasticNode<double>( "alpha", new ExponentialDistribution(one) );
    StochasticNode<double> *beta = new StochasticNode<double>( "beta", new ExponentialDistribution(one) );

    // base frequencies prior
    TypedDagNode<double> *phi;
    if(rootprior){
    	phi = new StochasticNode<double>( "phi", new TruncatedDistributionUnnormalized( new BetaDistribution( alpha, beta ), new ConstantNode<double>("rootmin", new double(rootmin) ), new ConstantNode<double>("rootmax", new double(rootmax) ) ) );
    }else{
    	phi = new StochasticNode<double>( "phi", new BetaDistribution( alpha,beta ) );
    }

    // branch frequency prior
	if(heterogeneous){
		if(dpp){
			// The prior mean number of rate categories induces an expected concentration parameter of ~meanCP
			double priorMean = 12.0;
			double meanCP = RbStatistics::Helper::dppConcParamFromNumTables(priorMean, (double)numBranches + 1);

			// Setting up the hyper prior on the concentratino parameter
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
		}else{
			if(mixture > 1){
				ConstantNode<std::vector<double> > *probs = new ConstantNode<std::vector<double> >("conc", new std::vector<double>(mixture,1.0) );
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
			}else{
				for (unsigned int i = 0 ; i < numBranches ; i++ ) {
					std::ostringstream pi_name;
					pi_name << "pi(" << i << ")";
					pi_stat.push_back(new StochasticNode<double>( pi_name.str(), new BetaDistribution(alpha,beta) ) );

					std::ostringstream q_name;
					q_name << "q(" << i << ")";
					qs.push_back(new DeterministicNode<RateMatrix>( q_name.str(), new FreeBinaryRateMatrixFunction(pi_stat[i]) ));
				}
			}
			pi_vector = new DeterministicNode< std::vector< double > >( "pi", new VectorFunction<double>( pi_stat ) );
			qs_node = new DeterministicNode< RbVector< RateMatrix > >( "q_vector", new RbVectorFunction<RateMatrix>(qs) );
		}
	}else{
		if(mixture > 1){
			std::cerr << "error: branch mixture frequencies only allowed with time-heterogeneous model\n";
			exit(1);
		}
		q = new DeterministicNode<RateMatrix>( "q", new FreeBinaryRateMatrixFunction(phi) );
	}

	// branch length hyperprior
	ContinuousStochasticNode *mu = new ContinuousStochasticNode("mu", new ExponentialDistribution(ten) );
	DeterministicNode<double> *rec_mu = new DeterministicNode<double>( "rec_mu", new BinaryDivision<double,double,double>(one,mu));

	// branch length prior
	for (unsigned int i = 0 ; i < numBranches ; i++ ) {
		std::ostringstream br_name;
		br_name << "br(" << i << ")";
		ContinuousStochasticNode* tmp_branch_rate = new ContinuousStochasticNode( br_name.str(), new ExponentialDistribution(rec_mu));
		branchRates.push_back( tmp_branch_rate );
		branchRates_nonConst.push_back( tmp_branch_rate );
	}
	DeterministicNode< std::vector< double > >* br_vector = new DeterministicNode< std::vector< double > >( "br_vector", new VectorFunction< double >( branchRates ) );

	// RAS prior
    ContinuousStochasticNode *lambda = new ContinuousStochasticNode("lambda", new ExponentialDistribution(one) );
    ConstantNode<double> *q1 = new ConstantNode<double>("q1", new double(0.125) );
    DeterministicNode<double> *q1_value = new DeterministicNode<double>("q1_value", new QuantileFunction(q1, new GammaDistribution(lambda, lambda) ) );
    ConstantNode<double> *q2 = new ConstantNode<double>("q2", new double(0.375) );
    DeterministicNode<double> *q2_value = new DeterministicNode<double>("q2_value", new QuantileFunction(q2, new GammaDistribution(lambda, lambda) ) );
    ConstantNode<double> *q3 = new ConstantNode<double>("q3", new double(0.625) );
    DeterministicNode<double> *q3_value = new DeterministicNode<double>("q3_value", new QuantileFunction(q3, new GammaDistribution(lambda, lambda) ) );
    ConstantNode<double> *q4 = new ConstantNode<double>("q4", new double(0.875) );
    DeterministicNode<double> *q4_value = new DeterministicNode<double>("q4_value", new QuantileFunction(q4, new GammaDistribution(lambda, lambda) ) );
    std::vector<const TypedDagNode<double>* > gamma_rates = std::vector<const TypedDagNode<double>* >();
    gamma_rates.push_back(q1_value);
    gamma_rates.push_back(q2_value);
    gamma_rates.push_back(q3_value);
    gamma_rates.push_back(q4_value);
    
    DeterministicNode<std::vector<double> > *site_rates = new DeterministicNode<std::vector<double> >( "site_rates", new VectorFunction<double>(gamma_rates) );
    DeterministicNode<std::vector<double> > *site_rates_norm = new DeterministicNode<std::vector<double> >( "site_rates_norm", new NormalizeVectorFunction(site_rates) );
            
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
    GeneralBranchHeterogeneousCharEvoModel<StandardState, BranchLengthTree> *charModel = new GeneralBranchHeterogeneousCharEvoModel<StandardState, BranchLengthTree>(psi, 2, true, data[0]->getNumberOfCharacters());
    //GeneralBranchHeterogeneousCharEvoModel<StandardState, TimeTree> *charModel = new GeneralBranchHeterogeneousCharEvoModel<StandardState, TimeTree>(tau, 2, true, data[0]->getNumberOfCharacters() );

    std::vector<const TypedDagNode<double> *> rf_vec(2);
	rf_vec[1] = phi;
	rf_vec[0] = new DeterministicNode<double>( "pi0", new BinarySubtraction<double,double,double>(one,phi));
	TypedDagNode<std::vector<double> > *rf = new DeterministicNode< std::vector< double > >( "rf", new VectorFunction< double >( rf_vec ) );

    if(heterogeneous){
    	charModel->setRateMatrix( qs_node );
    	charModel->setRootFrequencies( rf );
    }else{
    	charModel->setRateMatrix( q );
    	charModel->setRootFrequencies( rf );
    }
    charModel->setClockRate( one );
    charModel->setSiteRates( site_rates_norm );

    StochasticNode< AbstractCharacterData > *charactermodel = new StochasticNode< AbstractCharacterData >("S", charModel );
    charactermodel->clamp( data[0] );
    std::cout << "data ok\n";
    
    
    /* add the moves */
    std::vector<Move*> moves;
	std::vector<Monitor*> monitors;
	std::vector<DagNode*> monitoredNodes;

    //moves.push_back( new NearestNeighborInterchange( tau, 2.0 ) );
    if(!treeFixed)
    	moves.push_back( new SubtreePruneRegraft( tau, 5.0 ) );

    DeterministicNode<double>* treeLength = new DeterministicNode<double >( "length", new TreeLengthStatistic<BranchLengthTree>(psi) );
    monitoredNodes.push_back( treeLength );

    DeterministicNode<double>* meanpi;
	if(heterogeneous){
		meanpi = new DeterministicNode<double >( "mean_pi", new MeanFunction(pi_vector) );
		monitoredNodes.push_back( meanpi );
	}

    moves.push_back( new ScaleMove(mu, 1.0, true, 1.0) );
    monitoredNodes.push_back( mu );

    moves.push_back( new ScaleMove(lambda, 1.0, true, 1.0) );
    monitoredNodes.push_back( lambda );

    if(mixture > 1){
    	moves.push_back( new MixtureAllocationMove<double>((StochasticNode<double>*)phi, 2.0 ) );
    }else if(!dpp){
    	moves.push_back( new ScaleMove((StochasticNode<double>*)phi, 1.0, true, 5.0 ) );
    }
    monitoredNodes.push_back( phi );

    moves.push_back( new ScaleMove(alpha, 1.0, true, 1.0) );
    monitoredNodes.push_back( alpha );

    moves.push_back( new ScaleMove(beta, 1.0, true, 1.0) );
	monitoredNodes.push_back( beta );

	for (size_t i = 0 ; i < numBranches ; i ++ ) {
		if(heterogeneous && !mixture && !dpp){
			moves.push_back( new ScaleMove((StochasticNode<double>*)pi_stat[i], 1.0, true, 2.0 ) );
		}else if(mixture > 1){
			moves.push_back( new MixtureAllocationMove<double>((StochasticNode<double>*)pi_stat[i], 2.0 ) );
		}
		moves.push_back( new ScaleMove(branchRates_nonConst[i], 2.0, true, 1.0 ) );
	}

	if(heterogeneous){
		if(dpp){
			monitoredNodes.push_back(numCats);
			monitoredNodes.push_back(cp);

			moves.push_back( new ScaleMove(dpA, 1.0, true, 1.0) );
			monitoredNodes.push_back(dpA);

			moves.push_back( new ScaleMove(dpB, 1.0, true, 1.0) );
			monitoredNodes.push_back(dpB);

			moves.push_back( new DPPScaleCatValsMove( (StochasticNode<std::vector<double> >*)pi_vector, log(2.0), 2.0 ) );
			moves.push_back( new DPPAllocateAuxGibbsMove<double>( (StochasticNode<std::vector<double> >*)pi_vector, 4, 2.0 ) );
			moves.push_back( new DPPGibbsConcentrationMove<double>( cp, numCats, dpA, dpB, (int)numBranches + 1, 2.0 ) );
		}else if(mixture > 1){
			std::vector<Move*> mixmoves;
			for (size_t i = 0 ; i < mixture; i ++ ) {
				//mixmoves.push_back( new ScaleMove((StochasticNode<double>*)pi_cats[i], 1.0, true, 2.0 ) );
				moves.push_back( new ScaleMove((StochasticNode<double>*)pi_cats[i], 1.0, true, 2.0 ) );
				monitoredNodes.push_back((StochasticNode<double>*)pi_cats[i]);
			}
			//moves.push_back(new MultiMove(mixmoves,one,2.0,true));
		}
	}


    bool useParallelMcmcmc = (numChains > 1);
	monitors.push_back( new FileMonitor( monitoredNodes, every, name+".trace", "\t", false, true, false, useParallelMcmcmc, useParallelMcmcmc, false ) );
	if(!treeFixed){
		if(heterogeneous){
			std::set<TypedDagNode<std::vector<double> > *> piset;
			piset.insert(pi_vector);
			monitors.push_back( new NexusTreeMonitor( psi, piset, every, name+".treelist", useParallelMcmcmc) );
			//monitors.push_back( new ExtendedNewickTreeMonitor( tau, piset, every, name+".treelist", "\t", false, false, false, useParallelMcmcmc) );
		}else{
			monitors.push_back( new NewickTreeMonitor( psi, every, name+".treelist", useParallelMcmcmc) );
		}
	}
	if(ppred){
		monitors.push_back( new PosteriorPredictiveStateFrequencyMonitor( charactermodel, every, name+".ppred", useParallelMcmcmc) );
	}
	if(cvdata.size() > 0){
		monitors.push_back( new CrossValidationScoreMonitor( charactermodel, cvdata[0], every, name+".cv", useParallelMcmcmc) );
	}
	Model myModel = Model(charactermodel);
	std::cout << "model okay\n";

	int numProcesses = numChains;
	double startingHeat = 1.0;
	ParallelMcmcmc myPmc3(myModel, moves, monitors, "random", numChains, numProcesses, swapInterval, deltaTemp, sigmaTemp, startingHeat);
    std::cout << "running\n";
	myPmc3.run(until);
    
	for (std::vector<Move*>::iterator it = moves.begin(); it != moves.end(); ++it) {
		const Move *theMove = *it;
		delete theMove;
	}
	for (std::vector<Monitor*>::iterator it = monitors.begin(); it != monitors.end(); ++it) {
		const Monitor *theMonitor = *it;
		delete theMonitor;
	}
    
    return true;
}

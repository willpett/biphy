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
#include "Mcmc.h"
#include "Model.h"
#include "ModelMonitor.h"
#include "Monitor.h"
#include "Move.h"
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
#include "UniformDistribution.h"
#include "UniformConstrainedTimeTreeDistribution.h"
#include "UniformRootedTopologyDistribution.h"
#include "UniformTimeTreeDistribution.h"
#include "VectorFunction.h"

using namespace RevBayesCore;

TestBranchHeterogeneousBinaryModel::TestBranchHeterogeneousBinaryModel(const std::string &datafile,
									const std::string &name,
									const std::string &treefile,
									const std::string &outgroupfile,
									const std::string &cvfile,
									bool heterogeneous,
									bool ppred,
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
		ppred( ppred ),
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
    
    /* set up the model graph */
    
    //////////////////////
    // first the priors //
    //////////////////////
    
    // constant nodes
    ConstantNode<int> *idx = new ConstantNode<int>("idx", new int(2) );
    ConstantNode<double> *one = new ConstantNode<double>("one", new double(1.0) );
    ConstantNode<double> *ten = new ConstantNode<double>("ten", new double(10.0) );

    // base frequencies hyperprior
    StochasticNode<double> *bf_alpha = new StochasticNode<double>( "bf_alpha", new ExponentialDistribution(one) );
    StochasticNode<double> *bf_beta = new StochasticNode<double>( "bf_beta", new ExponentialDistribution(one) );

    std::vector< const TypedDagNode<double>* > bf_rootvec(2);
    bf_rootvec[0] = bf_alpha;
    bf_rootvec[1] = bf_beta;
    DeterministicNode<std::vector<double> > *bf = new DeterministicNode<std::vector<double> >( "bf", new VectorFunction< double >( bf_rootvec ) );

    // base frequencies prior
    TypedDagNode<std::vector<double> > *rf;
    if(!meanRootFreq)
    	rf = new StochasticNode<std::vector<double> >( "rf", new DirichletDistribution( bf ) );

	//Declaring a vector of matrices
	size_t numBranches = 2*data[0]->getNumberOfTaxa() - 2;
	std::vector< const TypedDagNode< RateMatrix >* > qs;
	std::vector<StochasticNode < std::vector<double> >* > pis;
	std::vector<const TypedDagNode<double> *> pi_stat;

    // declaring a vector of clock rates
	std::vector<const TypedDagNode<double> *> branchRates;
	std::vector< ContinuousStochasticNode *> branchRates_nonConst;
/*
	// The prior mean number of rate categories induces an expected concentration parameter of ~meanCP
	double priorMean = 12.0;
	double meanCP = RbStatistics::Helper::dppConcParamFromNumTables(priorMean, (double)numBranches);

	// Setting up the hyper prior on the concentratino parameter
	// This hyperprior is fully conditional on the DPP using a gamma distribution
	//    the parameters of the gamma distribution are set so that the expectation of the hyperprior = meanCP
	//    where meanCP = dpA / dpB
	ContinuousStochasticNode *dpA  = new ContinuousStochasticNode("dp_a", new ExponentialDistribution(one));
	ContinuousStochasticNode *dpB  = new ContinuousStochasticNode("dp_b", new ExponentialDistribution(one));
	StochasticNode<double> *cp = new StochasticNode<double>("DPP.cp", new GammaDistribution(dpA, dpB) );

	// G_0 is an Beta distribution
	TypedDistribution<double> *g = new BetaDistribution(one,one);

	// branchRates ~ DPP(g, cp, numBranches)
	StochasticNode<std::vector<double> > *pi_vector = new StochasticNode<std::vector<double> >("pi_vector", new DirichletProcessPriorDistribution<double>(g, cp, (int)numBranches) );

	// a deterministic node for calculating the number of rate categories (required for the Gibbs move on cp)
	DeterministicNode<int> *numCats = new DeterministicNode<int>("DPPNumCats", new DppNumTablesStatistic<double>(pi_vector) );
*/
	DeterministicNode< RateMatrix >* q = NULL;
	DeterministicNode< RbVector< RateMatrix > >* qs_node = NULL;
	if(!heterogeneous){
		q = new DeterministicNode<RateMatrix>( "q", new FreeBinaryRateMatrixFunction(rf) );
	}else{
		//qs_node = new DeterministicNode< RbVector< RateMatrix > >( "q", new FreeBinaryRateMatrixVectorFunction(pi_vector) );
	}

	ContinuousStochasticNode *rec_lambda = new ContinuousStochasticNode("rec_lambda", new ExponentialDistribution(ten) );
	DeterministicNode<double> *lambda = new DeterministicNode<double>( "lambda", new BinaryDivision<double,double,double>(one,rec_lambda));
    	for (unsigned int i = 0 ; i < numBranches ; i++ ) {
    		if(heterogeneous){
			std::ostringstream pi_name;
			pi_name << "pi(" << i << ")";
			pis.push_back(new StochasticNode<std::vector<double> >( pi_name.str(), new DirichletDistribution(bf) ) );
			pi_stat.push_back(new DeterministicNode<double>( pi_name.str()+"_1", new VectorIndexOperator<double>(pis.back(),idx)));
			std::ostringstream q_name;
			q_name << "q(" << i << ")";
			qs.push_back(new DeterministicNode<RateMatrix>( q_name.str(), new FreeBinaryRateMatrixFunction(pis[i]) ));
    		}

		std::ostringstream br_name;
		br_name << "br(" << i << ")";
		ContinuousStochasticNode* tmp_branch_rate = new ContinuousStochasticNode( br_name.str(), new ExponentialDistribution(lambda));
		branchRates.push_back( tmp_branch_rate );
		branchRates_nonConst.push_back( tmp_branch_rate );
	}
    DeterministicNode< std::vector< double > >* pi_vector = new DeterministicNode< std::vector< double > >( "pi", new VectorFunction< double >( pi_stat ) );
	DeterministicNode< std::vector< double > >* br_vector = new DeterministicNode< std::vector< double > >( "br_vector", new VectorFunction< double >( branchRates ) );

	if(heterogeneous){
		qs_node = new DeterministicNode< RbVector< RateMatrix > >( "q_vector", new RbVectorFunction<RateMatrix>(qs) );
	}

    ContinuousStochasticNode *alpha = new ContinuousStochasticNode("alpha", new ExponentialDistribution(one) );
    ConstantNode<double> *q1 = new ConstantNode<double>("q1", new double(0.125) );
    DeterministicNode<double> *q1_value = new DeterministicNode<double>("q1_value", new QuantileFunction(q1, new GammaDistribution(alpha, alpha) ) );
    ConstantNode<double> *q2 = new ConstantNode<double>("q2", new double(0.375) );
    DeterministicNode<double> *q2_value = new DeterministicNode<double>("q2_value", new QuantileFunction(q2, new GammaDistribution(alpha, alpha) ) );
    ConstantNode<double> *q3 = new ConstantNode<double>("q3", new double(0.625) );
    DeterministicNode<double> *q3_value = new DeterministicNode<double>("q3_value", new QuantileFunction(q3, new GammaDistribution(alpha, alpha) ) );
    ConstantNode<double> *q4 = new ConstantNode<double>("q4", new double(0.875) );
    DeterministicNode<double> *q4_value = new DeterministicNode<double>("q4_value", new QuantileFunction(q4, new GammaDistribution(alpha, alpha) ) );
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
    GeneralBranchHeterogeneousCharEvoModel<StandardState, BranchLengthTree> *charModel = new GeneralBranchHeterogeneousCharEvoModel<StandardState, BranchLengthTree>(psi, 2, true, data[0]->getNumberOfCharacters(), meanRootFreq);
    //GeneralBranchHeterogeneousCharEvoModel<StandardState, TimeTree> *charModel = new GeneralBranchHeterogeneousCharEvoModel<StandardState, TimeTree>(tau, 2, true, data[0]->getNumberOfCharacters() );

    if(heterogeneous){
    	charModel->setRateMatrix( qs_node );
    	charModel->setRootFrequencies( rf );
    }else{
    	charModel->setRateMatrix( q );
    	charModel->setRootFrequencies( rf );
    }
    //ContinuousStochasticNode* rho = new ContinuousStochasticNode( "rho", new ExponentialDistribution(one));
    charModel->setClockRate( one );
    charModel->setSiteRates( site_rates_norm );

    StochasticNode< AbstractCharacterData > *charactermodel = new StochasticNode< AbstractCharacterData >("S", charModel );
    charactermodel->clamp( data[0] );
    std::cout << "data ok\n";
    
    
    /* add the moves */
    std::vector<Move*> moves;
    //moves.push_back( new NearestNeighborInterchange( tau, 2.0 ) );
    if(!treeFixed)
    	moves.push_back( new SubtreePruneRegraft( tau, 5.0 ) );

    moves.push_back( new SimplexSingleElementScale((StochasticNode<std::vector<double> >*)rf, 1.0, true, 1.0 ) );

    moves.push_back( new ScaleMove(alpha, 1.0, true, 1.0) );
    moves.push_back( new ScaleMove(bf_alpha, 1.0, true, 1.0) );
    moves.push_back( new ScaleMove(bf_beta, 1.0, true, 1.0) );
    moves.push_back( new ScaleMove(rec_lambda, 1.0, true, 1.0) );
    //moves.push_back( new ScaleMove(dpA, 1.0, true, 1.0) );
    //moves.push_back( new ScaleMove(dpB, 1.0, true, 1.0) );
    //moves.push_back( new DPPScaleCatValsMove( pi_vector, log(2.0), 2.0 ) );
    //moves.push_back( new DPPAllocateAuxGibbsMove<double>( pi_vector, 4, 2.0 ) );
    //moves.push_back( new DPPGibbsConcentrationMove<double>( cp, numCats, dpA, dpB, (int)numBranches, 2.0 ) );
 
    for (unsigned int i = 0 ; i < numBranches ; i ++ ) {
    	if(heterogeneous){
    		moves.push_back( new SimplexSingleElementScale( pis[i], 1.0, true, 2.0 ) );
    	}
    	moves.push_back( new ScaleMove(branchRates_nonConst[i], 2.0, true, 1.0 ) );
    }
    std::cout << "moves okay\n";
    bool useParallelMcmcmc = (numChains > 1);

    /* add the monitors */
	std::vector<Monitor*> monitors;
	std::vector<DagNode*> monitoredNodes;
	DeterministicNode<double>* treeLength = new DeterministicNode<double >( "length", new TreeLengthStatistic<BranchLengthTree>(psi) );
	monitoredNodes.push_back( treeLength );
	monitoredNodes.push_back( alpha );
	monitoredNodes.push_back( lambda );
	monitoredNodes.push_back( bf_alpha );
	monitoredNodes.push_back( bf_beta );
	//monitoredNodes.push_back( numCats );
	//monitoredNodes.push_back( cp );
	//monitoredNodes.push_back( dpA );
	//monitoredNodes.push_back( dpB );
	monitoredNodes.push_back( rf );

	monitors.push_back( new FileMonitor( monitoredNodes, every, name+".trace", "\t", false, true, false, useParallelMcmcmc, useParallelMcmcmc, false ) );
	if(!treeFixed){
		/*
		if(heterogeneous){
			std::set<TypedDagNode<std::vector<double> > *> piset;
			piset.insert(pi_vector);
			monitors.push_back( new NexusTreeMonitor( psi, piset, every, name+".treelist", useParallelMcmcmc) );
			//monitors.push_back( new ExtendedNewickTreeMonitor( tau, piset, every, name+".treelist", "\t", false, false, false, useParallelMcmcmc) );
		}else{
		*/
			monitors.push_back( new NewickTreeMonitor( psi, every, name+".treelist", useParallelMcmcmc) );
		//}
	}
	if(ppred){
		monitors.push_back( new PosteriorPredictiveStateFrequencyMonitor( charactermodel, every, name+".ppred", useParallelMcmcmc) );
	}
	if(cvdata.size() > 0){
		monitors.push_back( new CrossValidationScoreMonitor( charactermodel, cvdata[0], every, name+".cv", useParallelMcmcmc) );
	}
	std::cout << "monitors okay\n";
    	Model myModel = Model(charactermodel);
    	std::cout << "model okay\n";

	int numProcesses = numChains;
	double startingHeat = 1.0;
    	ParallelMcmcmc myPmc3(myModel, moves, monitors, "random", numChains, numProcesses, swapInterval, deltaTemp, sigmaTemp, startingHeat);
    	std::cout << "running\n";
	myPmc3.run(until);
	myPmc3.printOperatorSummary();
    
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

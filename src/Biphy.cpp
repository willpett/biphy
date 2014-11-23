#include "AbstractCharacterData.h"
#include "Biphy.h"
#include "BinaryCharEvoModel.h"
#include "BinaryDivision.h"
#include "BinaryDolloCompatibleMonitor.h"
#include "BinaryMultiplication.h"
#include "BinarySubtraction.h"
#include "BetaDistribution.h"
#include "BetaSimplexMove.h"
#include "Clade.h"
#include "CladeReader.h"
#include "ConstantNode.h"
#include "ConstantFunction.h"
#include "CrossValidationScoreMonitor.h"
#include "DeterministicNode.h"
#include "DirichletDistribution.h"
#include "DirichletProcessPriorDistribution.h"
#include "DolloBinaryCharEvoModel.h"
#include "DPPAllocateAuxGibbsMove.h"
#include "DPPGibbsConcentrationMove.h"
#include "DppNumTablesStatistic.h"
#include "DPPBetaSimplexMove.h"
#include "ExponentialDistribution.h"
#include "FastaWriter.h"
#include "FileMonitor.h"
#include "FreeBinaryRateMatrixFunction.h"
#include "FreeBinaryRateMatrixVectorFunction.h"
#include "GammaDistribution.h"
#include "ParallelMcmcmc.h"
#include "MixtureDistribution.h"
#include "MixtureAllocationMove.h"
#include "Mcmc.h"
#include "MeanFunction.h"
#include "Model.h"
#include "Monitor.h"
#include "Move.h"
#include "MultiMove.h"
#include "NclReader.h"
#include "NewickTreeMonitor.h"
#include "NexusTreeMonitor.h"
#include "ParallelMcmcmc.h"
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
#include "TimeTree.h"
#include "TreeAssemblyFunction.h"
#include "TreeLengthStatistic.h"
#include "TruncatedDistributionUnnormalized.h"
#include "UniformDistribution.h"
#include "UniformRootedTopologyDistribution.h"
#include "VectorFunction.h"
#include "VectorScaleFunction.h"


using namespace RevBayesCore;

Biphy::Biphy(const std::string &datafile,
									const std::string &name,
									const std::string &treefile,
									const std::string &outgroupfile,
									ModelType modeltype,
									BranchPrior branchprior,
									RootPrior rootprior,
									CorrectionType correction,
									int dgam,
									int mixture,
									double rootmin,
									double rootmax,
									int every,
									int until,
									int numchains,
									int swapInterval,
									double delta,
									double sigma,
									bool saveall,
									bool nexus) :
		dataFile( datafile ),
		name( name ),
		treeFile( treefile ),
		outgroupFile( outgroupfile ),
		cvfile("None"),
		modeltype(modeltype),
		branchprior(branchprior),
		rootprior(rootprior),
		correction(correction),
		dgam( dgam),
		mixture( mixture),
		ppred(false),
		rootmin( rootmin),
		rootmax( rootmax),
		every( every ),
		until( until ),
		numChains( numchains),
		swapInterval( swapInterval ),
		delta( delta ),
		sigma( sigma ),
		saveall( saveall),
		nexus(nexus),
		readstream(false),
		restart(false)
{
	init();
    save();
}

Biphy::Biphy(const std::string &name, const std::string &cvfile, bool ppred, bool dolloMapping) :
		name( name ),
		cvfile(cvfile),
		ppred(ppred),
		readstream(true),
		restart(false),
		dolloMapping(dolloMapping)
{
	every = 1;

	open();
	init();
}

Biphy::Biphy(const std::string &name) :
		name( name ), cvfile("None"), ppred(false), readstream(false), restart(true)
{
	open();
	init();
}

Biphy::~Biphy() {
    // nothing to do
}

void Biphy::open( void ) {
	std::ifstream is((name + ".param").c_str());
	if (!is)        {
		std::cerr << "error : cannot find file : " << name << ".param\n";
		exit(1);
	}
	int i;

	is >> dataFile;
	is >> name;
	is >> treeFile;
	is >> outgroupFile;
	is >> i;
	modeltype = static_cast<ModelType>(i);
	is >> i;
	branchprior = static_cast<BranchPrior>(i);
	is >> i;
	rootprior = static_cast<RootPrior>(i);
	is >> i;
	correction = static_cast<CorrectionType>(i);
	is >> mixture;
	is >> dgam;
	is >> rootmin;
	is >> rootmax;
	is >> every;
	is >> until;
	is >> numChains;
	is >> swapInterval;
	is >> delta;
	is >> sigma;
	is >> saveall;
	is >> nexus;

	is.close();
}

void Biphy::save( void ) {
	std::ofstream os((name + ".param").c_str());

	os << dataFile << "\n";
	os << name << "\n";
	os << treeFile << "\n";
	os << outgroupFile << "\n";
	os << modeltype << "\n";
	os << branchprior << "\n";
	os << rootprior << "\n";
	os << correction << "\n";
	os << mixture << "\n";
	os << dgam << "\n";
	os << rootmin << "\n";
	os << rootmax << "\n";
	os << every << "\n";
	os << until << "\n";
	os << numChains << "\n";
	os << swapInterval << "\n";
	os << delta << "\n";
	os << sigma << "\n";
	os << saveall << "\n";
	os << nexus << "\n";

	os.close();
}



void Biphy::init( void ) {
    
    /* First, we read in the data */
    // the matrix
	RbSettings::userSettings().setPrintNodeIndex(false);

	if(dolloMapping && modeltype == DOLLO)
		throw(RbException("Error: simulation of Dollo incompatible mappings under Dollo model doesn't make sense"));

	//const string format("phylip|standard|sequential");
    std::vector<AbstractCharacterData*> data = NclReader::getInstance().readMatrices(dataFile);

    // data checks
    if(data.empty()){
		//std::cerr << "Error: failed to read datafile" << std::endl;
		exit(1);
	}

    if(data[0]->getDatatype() != "Standard"){
    	std::cerr << "Error: incompatible datatype '" << data[0]->getDatatype() << "'" << std::endl;
    	exit(1);
    }

    std::string symbols = ((DiscreteCharacterData<StandardState> *)data[0])->getCharacter(0,0).getStateLabels();
	std::stringstream states;
	states << "character states:\t";
	for(size_t c = 0; c < symbols.size(); c++){
		if(c == 0)
			states << symbols[c] << " = absent";
		if(c == 1)
			states << symbols[c] << " = present" ;
		if(c != symbols.size() - 1)
			states << ", ";
	}
	std::cout << states.str() << std::endl;

    std::vector<AbstractCharacterData*> cvdata;
	if(cvfile != "None"){
		cvdata = NclReader::getInstance().readMatrices(cvfile);
		std::cout << "computing cross-validation scores for " << cvfile << std::endl;

		// data checks
		if(cvdata[0]->getDatatype() != "Standard"){
			std::cerr << "Error: incompatible datatype '" << data[0]->getDatatype() << "'" << std::endl;
			exit(1);
		}

		std::string cvsymbols = ((DiscreteCharacterData<StandardState> *)cvdata[0])->getCharacter(0,0).getStateLabels();
		if(cvsymbols != symbols){
			std::cerr << "Error: character states in test alignment don't match the original dataset" << std::endl;
			exit(1);
		}
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
    
    if(modeltype > HOMOGENEOUS){
    	std::cout << "branch-wise stationary frequencies:\n";
    	switch(modeltype){
			case DPP:
				std::cout << "pi(i) ~ DPP(G,cp)\n";
				std::cout << "cp ~ Gamma(1,1)\n";
				//std::cout << "a = b = 1~ exponential of mean 1\n";
				std::cout << "G = Beta(beta1,beta2)\n";
				break;
			case MIXTURE:
				std::cout << "pi(i) ~ Beta(beta1,beta2) mixture with " << mixture << " components\n";
				break;
			default:
				std::cout << "pi(i) ~ Beta(beta1,beta2)\n";
				break;
    	}
	}else if(modeltype == DOLLO){
		std::cout << "dollo model";
	}else{
		std::cout << "time-homogeneous model\n";
	}

    std::cout << "\n";

    if(modeltype != DOLLO){
		if(rootprior == RIGID)
			std::cout << "rigid ";
		std::cout << "root frequency:\n";
		if(modeltype > HOMOGENEOUS){
			switch(modeltype){
				case DPP:
					std::cout << "phi ~ DPP(G,cp)\n";
					break;
				case MIXTURE:
					std::cout << "phi ~ Beta(beta1,beta2) mixture with " << mixture << " components\n";
					break;
				default:
					if(rootprior == TRUNCATED)
						std::cout << "phi ~ Beta(beta1,beta2) truncated on (" << rootmin << "," << rootmax << ")\n";
					else
						std::cout << "phi ~ Beta(beta1,beta2)\n";
					break;
			}
			std::cout << "\nbeta1, beta2 ~ exponential of mean 1\n";
		}else if(modeltype != DOLLO){
			if(rootprior == TRUNCATED){
				std::cout << "phi ~ Uniform(" << rootmin << "," << rootmax << ")\n";
			}else{
				std::cout << "phi ~ Beta(1,1)\n";
			}
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

    if(dgam > 1){
		std::cout << "\n";

		std::cout << "rates ~ discrete gamma with " << dgam << " rate categories\n";
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

    std::cout << "\n";

    //////////////////////
    // first the priors //
    //////////////////////
    
    // constant nodes
    ConstantNode<double> *one = new ConstantNode<double>("one", new double(1.0) );
    ConstantNode<double> *ten = new ConstantNode<double>("ten", new double(10.0) );
    ConstantNode<double> *zero = new ConstantNode<double>("zero", new double(0.0) );

    //number of branches
	size_t numBranches = 2*data[0]->getNumberOfTaxa() - 2;
	size_t w = numBranches > 0 ? (int) log10 ((double) numBranches) + 1 : 1;

	if(modeltype == MIXTURE)
		w = mixture > 0 ? (int) log10 ((double) mixture) + 1 : 1;


    // tree prior
    std::vector<std::string> names = data[0]->getTaxonNames();
    //StochasticNode<TimeTree> *tau = new StochasticNode<TimeTree>( "tau", new UniformConstrainedTimeTreeDistribution(one,names,constraints,outgroup) );
    StochasticNode<Topology> *tau = new StochasticNode<Topology>( "tau", new UniformRootedTopologyDistribution(names, std::vector<Clade>(), outgroup) );
    if(treeFixed){
    	tau->setValue(fixedTree, true);
    	tau->setIgnoreRedraw(true);
    }

    // base frequencies hyperprior
	StochasticNode<double> *beta1;
	StochasticNode<double> *beta2;

	// base frequencies prior
	TypedDagNode<double> *phi;
	if(modeltype == DOLLO){
		phi = new DeterministicNode<double>("phi", new ConstantFunction<double>(zero));
	}else{
		if(rootprior == TRUNCATED){
			if(modeltype > HOMOGENEOUS){
				beta1 = new StochasticNode<double>( "beta1", new ExponentialDistribution(one) );
				beta2 = new StochasticNode<double>( "beta2", new ExponentialDistribution(one) );
				phi = new StochasticNode<double>( "phi", new TruncatedDistributionUnnormalized( new BetaDistribution( beta1, beta2 ), new ConstantNode<double>("rootmin", new double(rootmin) ), new ConstantNode<double>("rootmax", new double(rootmax) ) ) );
			}else{
				phi = new StochasticNode<double>( "phi", new UniformDistribution( new ConstantNode<double>("rootmin", new double(rootmin) ), new ConstantNode<double>("rootmax", new double(rootmax) ) ) );
			}
		}else{
			if(modeltype > HOMOGENEOUS){
				beta1 = new StochasticNode<double>( "beta1", new ExponentialDistribution(one) );
				beta2 = new StochasticNode<double>( "beta2", new ExponentialDistribution(one) );
				phi = new StochasticNode<double>( "phi", new BetaDistribution( beta1,beta2 ) );
			}else{
				phi = new StochasticNode<double>( "phi", new BetaDistribution( one,one ) );
			}
		}
	}

	//mixture vectors
	std::vector<const TypedDagNode<double > *> pi_cats;
	TypedDagNode< std::vector< double > >* pi_mix;
	std::vector<const TypedDagNode<double> *> pi_stat;

	//dpp priors
	TypedDagNode<double> *dpA;
	TypedDagNode<double> *dpB;
	StochasticNode<double> *cp;
	DeterministicNode<int> *numCats;

	TypedDagNode< std::vector< double > >* pi_vector;
	TypedDagNode< std::vector< double > >* dpp_vector;

	std::vector<const TypedDagNode< RateMatrix >* > qs;
	DeterministicNode< RateMatrix >* q = NULL;
	DeterministicNode< RbVector< RateMatrix > >* qs_node = NULL;

    const TopologyNode &root = tau->getValue().getRoot();
    size_t left = root.getChild(0).getIndex();
    size_t right = root.getChild(1).getIndex();

    // branch frequency prior
	if(modeltype == HIERARCHICAL){
		for (size_t i = 0 ; i < numBranches ; i++ ) {
			std::ostringstream pi_name;
			pi_name << "pi(" << setfill('0') << setw(w) << i << ")";
			std::ostringstream q_name;
			q_name << "q(" << setfill('0') << setw(w) << i << ")";
			if((i == left || i == right) && rootprior == RIGID)
				pi_stat.push_back(new DeterministicNode<double>(pi_name.str(), new ConstantFunction<double>(phi)));
			else
				pi_stat.push_back(new StochasticNode<double>( pi_name.str(), new BetaDistribution(beta1,beta2) ) );
			qs.push_back(new DeterministicNode<RateMatrix>( q_name.str(), new FreeBinaryRateMatrixFunction(pi_stat[i]) ));
		}
		pi_vector = new DeterministicNode< std::vector< double > >( "pi", new VectorFunction<double>( pi_stat ) );
		qs_node = new DeterministicNode< RbVector< RateMatrix > >( "q_vector", new RbVectorFunction<RateMatrix>(qs) );
	}else if(modeltype == DPP){
		// Setting up the hyper prior on the concentration parameter
		// This hyperprior is fully conditional on the DPP using a gamma distribution
		//    the parameters of the gamma distribution are set so that the expectation of the hyperprior = meanCP
		//    where meanCP = dpA / dpB
		double meandpp = 1;
		dpA  = new ConstantNode<double>("dp_a", new double(meandpp) );//new ContinuousStochasticNode("dp_a", new ExponentialDistribution(one));
		dpB  = new ConstantNode<double>("dp_b", new double(1.0) );//new ContinuousStochasticNode("dp_b", new ExponentialDistribution(one));
		cp = new StochasticNode<double>("dpp.cp", new GammaDistribution(dpA, dpB) );

		// G_0 is an Beta distribution
		TypedDistribution<double> *g = new BetaDistribution(beta1,beta2);

		// branchRates ~ DPP(g, cp, numBranches)
		dpp_vector = new StochasticNode<std::vector<double> >("dpp_vector", new DirichletProcessPriorDistribution<double>(g, cp, numBranches + 1 - 2*(rootprior == RIGID)) );
		// a deterministic node for calculating the number of rate categories (required for the Gibbs move on cp)
		numCats = new DeterministicNode<int>("dpp.k", new DppNumTablesStatistic<double>((StochasticNode<std::vector<double> >*)dpp_vector) );

		delete phi;
		ConstantNode<int> *idx = new ConstantNode<int>("idx", new int(numBranches + 1 - 2*(rootprior == RIGID)) );
		phi = new DeterministicNode<double>("phi", new VectorIndexOperator<double>((StochasticNode<std::vector<double> >*)dpp_vector,idx));

		if(rootprior == RIGID){
			size_t comp = 1;
			std::vector<const TypedDagNode<double>* > tmp;
			for(size_t i = 0; i < numBranches; i++){
				std::ostringstream tmp_name;
				tmp_name << comp;

				if(i == left || i == right){
					tmp.push_back(new DeterministicNode<double>("tmp"+tmp_name.str(),new ConstantFunction<double>(phi)));
				}else{
					ConstantNode<int> *tmp_idx = new ConstantNode<int>("idx"+tmp_name.str(), new int(comp) );
					tmp.push_back(new DeterministicNode<double>("tmp"+tmp_name.str(), new VectorIndexOperator<double>((StochasticNode<std::vector<double> >*)dpp_vector,tmp_idx)));
					comp++;
				}
			}
			tmp.push_back(new DeterministicNode<double>("tmp1",new ConstantFunction<double>(phi)));
			pi_vector = new DeterministicNode< std::vector< double > >( "pi", new VectorFunction<double>( tmp ) );
		}else{
			pi_vector = new DeterministicNode<std::vector< double > >( "pi", new ConstantFunction<std::vector< double > >( dpp_vector ) );
		}

		qs_node = new DeterministicNode< RbVector< RateMatrix > >( "q_vector", new FreeBinaryRateMatrixVectorFunction(pi_vector,true) );
	}else if(modeltype == MIXTURE){
		ConstantNode<std::vector<double> > *probs = new ConstantNode<std::vector<double> >("probs", new std::vector<double>(mixture,1.0/mixture) );
		for(size_t i = 0; i < mixture; i++){
			std::ostringstream pi_name;
			pi_name << "pi(" << setfill('0') << setw(w) << i << ")";
			pi_cats.push_back(new StochasticNode<double >( pi_name.str(), new BetaDistribution(beta1,beta2) ) );
		}
		pi_mix = new DeterministicNode< std::vector< double > >( "pi_mix", new VectorFunction< double >( pi_cats ) );

		delete phi;
		phi = new StochasticNode<double>( "phi", new MixtureDistribution<double>(pi_mix,probs) );

		for (size_t i = 0 ; i < numBranches ; i++ ) {
			std::ostringstream q_name;
			q_name << "q(" << setfill('0') << setw(w) << i << ")";
			if((i == left || i == right) && rootprior == RIGID){
				pi_stat.push_back(new DeterministicNode<double>( q_name.str()+"_pi", new ConstantFunction< double >( phi ) ));
			}else{
				pi_stat.push_back(new StochasticNode<double>( q_name.str()+"_pi", new MixtureDistribution<double>(pi_mix,probs) ));
			}
			qs.push_back(new DeterministicNode<RateMatrix>( q_name.str(), new FreeBinaryRateMatrixFunction(pi_stat[i]) ));
		}

		pi_vector = new DeterministicNode< std::vector< double > >( "pi", new VectorFunction<double>( pi_stat ) );
		qs_node = new DeterministicNode< RbVector< RateMatrix > >( "q_vector", new RbVectorFunction<RateMatrix>(qs) );
	}else{
		q = new DeterministicNode<RateMatrix>( "q", new FreeBinaryRateMatrixFunction(phi) );
	}

	// branch length hyperprior

	ContinuousStochasticNode *mu;
	DeterministicNode<double> *rec_mu;
	TypedDagNode< std::vector<double> >* br_vector;

	StochasticNode< std::vector<double> >* br_times;
	TypedDagNode<double>* tree_length;

	// declaring a vector of clock rates
	std::vector<const TypedDagNode<double> *> branchRates;
	std::vector< ContinuousStochasticNode *> branchRates_nonConst;

	// branch length prior
	if(branchprior == EXPONENTIAL){
		mu = new ContinuousStochasticNode("mu", new ExponentialDistribution(ten) );
		rec_mu = new DeterministicNode<double>( "rec_mu", new BinaryDivision<double,double,double>(one,mu));
		for (size_t i = 0 ; i < numBranches ; i++ ) {
			std::ostringstream br_name;
			br_name << "br(" << setfill('0') << setw(w) << i << ")";

			ContinuousStochasticNode* tmp_branch_rate = new ContinuousStochasticNode( br_name.str(), new ExponentialDistribution(rec_mu));
			branchRates.push_back( tmp_branch_rate );
			branchRates_nonConst.push_back( tmp_branch_rate );
		}
		br_vector = new DeterministicNode< std::vector< double > >( "br_vector", new VectorFunction< double >( branchRates ) );
	}else if(branchprior == DIRICHLET){
		ConstantNode<std::vector<double> > *conc = new ConstantNode<std::vector<double> >("conc", new std::vector<double>(numBranches,1.0) );
		tree_length = new StochasticNode<double>("length", new GammaDistribution(one,one));
		br_times = new StochasticNode< std::vector< double > >( "br_times", new DirichletDistribution( conc));
		br_vector = new DeterministicNode< std::vector< double > >( "br_vector", new VectorScaleFunction(br_times,tree_length) );
	}

	// RAS prior
    ContinuousStochasticNode *alpha;
    std::vector<const TypedDagNode<double>* > gamma_rates = std::vector<const TypedDagNode<double>* >();
    DeterministicNode<std::vector<double> > *site_rates;
    DeterministicNode<std::vector<double> > *site_rates_norm;
    if(dgam > 1){
		alpha = new ContinuousStochasticNode("alpha", new ExponentialDistribution(one) );
		for(size_t cat = 0; cat < dgam; cat++){
			std::stringstream name;
			std::stringstream value_name;
			name << "q";
			name << cat+1;
			value_name << name << "_value";
			gamma_rates.push_back( new DeterministicNode<double>(value_name.str(), new QuantileFunction(new ConstantNode<double>(name.str(), new double((cat+1.0/2.0)/dgam) ), new GammaDistribution(alpha, alpha) ) ));
		}
		site_rates = new DeterministicNode<std::vector<double> >( "site_rates", new VectorFunction<double>(gamma_rates) );
		site_rates_norm = new DeterministicNode<std::vector<double> >( "site_rates_norm", new NormalizeVectorFunction(site_rates) );
    }
    
    DeterministicNode<BranchLengthTree> *psi = new DeterministicNode<BranchLengthTree>( "psi", new TreeAssemblyFunction(tau, br_vector) );
    std::cout << "tree okay\n";

    BinaryCharEvoModel<BranchLengthTree> *charModel;
    StochasticNode<double> *birthrate;
    if(modeltype == DOLLO){
		charModel = new DolloBinaryCharEvoModel<BranchLengthTree>(psi, tau, true, data[0]->getNumberOfCharacters(), correction);
	}else{
		charModel = new BinaryCharEvoModel<BranchLengthTree>(psi, true, data[0]->getNumberOfCharacters(), correction);
    }

    std::vector<const TypedDagNode<double> *> rf_vec(2);
	TypedDagNode<std::vector<double> > *rf;

    if(modeltype == DOLLO){
    	rf_vec[1] = one;
    	rf_vec[0] = zero;
    	rf = new DeterministicNode< std::vector< double > >( "rf", new VectorFunction< double >( rf_vec ) );
    }else{
    	rf_vec[1] = phi;
    	rf_vec[0] = new DeterministicNode<double>( "pi0", new BinarySubtraction<double,double,double>(one,phi));
    	rf = new DeterministicNode< std::vector< double > >( "rf", new VectorFunction< double >( rf_vec ) );
    }
    if(modeltype > HOMOGENEOUS){
    	charModel->setRateMatrix( qs_node );
    	charModel->setRootFrequencies( rf );
	}else{
		charModel->setRateMatrix( q );
		charModel->setRootFrequencies( rf );
	}

	if(dgam > 1)
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
		moves.push_back( new SubtreePruneRegraft( tau, 5.0, rootprior == RIGID) );

	for (size_t i = 0 ; i < numBranches ; i ++ ) {
		if(!(i == left || i == right) || rootprior == FREE)
			if(modeltype == HIERARCHICAL){
				moves.push_back( new BetaSimplexMove((StochasticNode<double>*)pi_stat[i], 1.0, true, 2.0 ) );
			}else if(modeltype == MIXTURE){
				moves.push_back( new MixtureAllocationMove<double>((StochasticNode<double>*)pi_stat[i], 2.0 ) );
			}
		if(branchprior == EXPONENTIAL)
			moves.push_back( new ScaleMove(branchRates_nonConst[i], 1.0, true, 1.0 ) );
	}

	if(branchprior == EXPONENTIAL){
		tree_length = new DeterministicNode<double >("length", new TreeLengthStatistic<BranchLengthTree>(psi) );
	}else if(branchprior == DIRICHLET){
		moves.push_back(new ScaleMove((StochasticNode<double>*)tree_length, 1.0, true, 1.0));
		moves.push_back(new SimplexSingleElementScale(br_times, 2.0, true, (int)numBranches/3));
	}
	monitoredNodes.push_back( tree_length );

	if(modeltype > DOLLO && modeltype < DPP){
		if(modeltype == MIXTURE){
			moves.push_back( new MixtureAllocationMove<double>((StochasticNode<double>*)phi, 2.0 ) );
		}else{
			moves.push_back( new BetaSimplexMove((StochasticNode<double>*)phi, 1.0, true, 5.0 ) );
		}
	}
	if(modeltype != DOLLO)
		monitoredNodes.push_back( phi );

	DeterministicNode<double>* meanpi;
	if(modeltype > HOMOGENEOUS){
		meanpi = new DeterministicNode<double >( "mean_pi", new MeanFunction(pi_vector) );
		monitoredNodes.push_back( meanpi );
	}

	if(branchprior == EXPONENTIAL){
		moves.push_back( new ScaleMove(mu, 1.0, true, 1.0) );
		monitoredNodes.push_back( mu );
	}

	if(dgam > 1){
		moves.push_back( new ScaleMove(alpha, 1.0, true, 1.0) );
		monitoredNodes.push_back( alpha );
		//monitoredNodes.push_back(site_rates_norm);
	}

	if(modeltype > HOMOGENEOUS){
		moves.push_back( new ScaleMove(beta1, 1.0, true, 1.0) );
		monitoredNodes.push_back( beta1 );

		moves.push_back( new ScaleMove(beta2, 1.0, true, 1.0) );
		monitoredNodes.push_back( beta2 );

		if(modeltype == DPP){
			monitoredNodes.push_back(numCats);
			monitoredNodes.push_back(cp);

			//moves.push_back( new ScaleMove(dpA, 1.0, true, 1.0) );
			//monitoredNodes.push_back(dpA);

			//moves.push_back( new ScaleMove(dpB, 1.0, true, 1.0) );
			//monitoredNodes.push_back(dpB);

			moves.push_back( new DPPBetaSimplexMove( (StochasticNode<std::vector<double> >*)dpp_vector, 1.0 , 1.0 ) );
			moves.push_back( new DPPAllocateAuxGibbsMove<double>( (StochasticNode<std::vector<double> >*)dpp_vector, 4, 2.0 ) );
			moves.push_back( new DPPGibbsConcentrationMove<double>( cp, numCats, dpA, dpB, numBranches + 1, 2.0 ) );
		}else if(modeltype == MIXTURE){
			std::vector<Move*> mixmoves;
			for (size_t i = 0 ; i < mixture; i ++ ) {
				//mixmoves.push_back( new ScaleMove((StochasticNode<double>*)pi_cats[i], 1.0, true, 2.0 ) );
				moves.push_back( new BetaSimplexMove((StochasticNode<double>*)pi_cats[i], 1.0, true, 2.0 ) );
				monitoredNodes.push_back((StochasticNode<double>*)pi_cats[i]);
			}
			//moves.push_back(new MultiMove(mixmoves,one,2.0,true));
		}
	}

	Model myModel = Model(charactermodel);
	std::cout << "model okay\n";

	bool useParallelMcmcmc = (numChains > 1);

	if(readstream){
		if(ppred)
			monitors.push_back( new PosteriorPredictiveStateFrequencyMonitor( charactermodel, every, name+".ppred",useParallelMcmcmc) );
		if(cvdata.size() > 0)
			monitors.push_back( new CrossValidationScoreMonitor( charactermodel, cvdata[0], every, name+".cv",useParallelMcmcmc) );
		if(dolloMapping)
			monitors.push_back( new BinaryDolloCompatibleMonitor<BranchLengthTree>( charactermodel, every, name+".dollo.fa") );
	}else{
		monitors.push_back( new FileMonitor( monitoredNodes, every, name+".trace", "\t", false, true, false, useParallelMcmcmc || restart, useParallelMcmcmc, false ) );
		if(!treeFixed){
			monitors.push_back( new NewickTreeMonitor( psi, every, name+".treelist", useParallelMcmcmc || restart) );
			if(nexus){
				std::set<TypedDagNode<std::vector<double> > *> piset;
				if(modeltype > HOMOGENEOUS)
					piset.insert(pi_vector);
				monitors.push_back( new NexusTreeMonitor( psi, piset, every, name+".treelist.nex", useParallelMcmcmc || restart) );
				//monitors.push_back( new ExtendedNewickTreeMonitor( tau, piset, every, name+".treelist", "\t", false, false, false, useParallelMcmcmc) );
			}
		}
	}

	double startingHeat = 1.0;
	mcmc = new ParallelMcmcmc(myModel, moves, monitors, name+".stream", "random", every, numChains, numChains, swapInterval, delta, sigma, startingHeat, saveall);
}

void Biphy::run( void ) {
	if(readstream){
		if(!saveall){
			std::cerr << "Error: "+name+".stream does not contain complete output. run again with -s option\n\n";
			exit(1);
		}
		std::cout << "reading from stream\n";
		mcmc->readStream(until);
	}else{
		if(restart)
			mcmc->restore();

		std::cout << "running\n";
		mcmc->run(until);
	}

	std::cout << "done\n";
}

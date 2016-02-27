#include "Biphy.h"

#include "BetaDistribution.h"
#include "BetaSimplexMove.h"
#include "BinaryCharacterDataReader.h"
#include "BinaryDivision.h"
#include "BinaryDolloSubstitutionModel.h"
#include "BinarySubstitutionModel.h"
#include "ConstantFunction.h"
#include "CladeReader.h"
#include "CrossValidationScoreMonitor.h"
#include "DirichletDistribution.h"
#include "DiscreteBetaValuesFunction.h"
#include "EntropyFunction.h"
#include "ExponentialDistribution.h"
#include "FileMonitor.h"
#include "GammaDistribution.h"
#include "LogitFunction.h"
#include "MappingMonitor.h"
#include "MeanFunction.h"
#include "MixtureDistribution.h"
#include "MixtureAllocationMove.h"
#include "NexusTreeMonitor.h"
#include "NewickTreeMonitor.h"
#include "NewickTreeReader.h"
#include "NormalDistribution.h"
#include "NormalizeVectorFunction.h"
#include "PerSiteLnProbMonitor.h"
#include "PosteriorPredictiveStateFrequencyMonitor.h"
#include "QuantileFunction.h"
#include "RandomNumberFactory.h"
#include "RandomNumberGenerator.h"
#include "ScaleMove.h"
#include "Settings.h"
#include "SimplexMove.h"
#include "SimplexSingleElementScale.h"
#include "SlidingMove.h"
#include "SubtreePruneRegraft.h"
#include "SumFunction.h"
#include "TruncatedDistributionUnnormalized.h"
#include "UniformDistribution.h"
#include "UniformTopologyDistribution.h"
#include "VectorFunction.h"
#include "VectorScaleFunction.h"

#include <iomanip>

Biphy::Biphy(const std::string n,
        const std::string df,
		const std::string cv,
        const std::string t,
        const std::string out,
        ModelPrior::Type mt,
        BranchPrior::Type br,
        RootPrior::Type rt,
        int c,
        int d,
		int b,
		bool sb,
        int mix,
        double rm,
        double rx,
        int e,
        int u,
        int num,
        int swa,
        double de,
        double si,
        bool sav,
        bool nex) :
    dataFile( df ),
	cvfile(cv),
    name( n ),
    treeFile( t ),
    correction(c),
    dgam( d),
	dbeta(b),
	asymmbeta(sb),
    every( e ),
    until( u ),
    numChains( num),
    swapInterval( swa ),
    delta( de ),
    sigma( si ),
    saveall( sav),
    readstream(false),
    restart(false),
    mixture( mix),
    ppred(false),
    perSiteLnProbs(false),
    outgroupFile( out ),
    modeltype(mt),
    branchprior(br),
    rootprior(rt),
    nexus(nex),
    outgroup(std::vector<std::string>())
{
    save();
}


Biphy::Biphy(const std::string n, const std::string cv, bool pp, bool dm, bool site, bool anc) :
    name(n),
    cvfile(cv),
    ppred(pp),
    dolloMapping(dm),
    perSiteLnProbs(site),
    outgroup(std::vector<std::string>()),
    ancestral(anc)
{

    readstream = true;
    restart = false;

    open();

    every = 1;
}


Biphy::Biphy(const std::string n) :
        name( n ), readstream(false), restart(true), cvfile("None"), ppred(false), perSiteLnProbs(false), outgroup(std::vector<std::string>())
{
    open();
}


void Biphy::open( void ) {
    std::ifstream is((name + ".param").c_str());
    if (!is)        {
        std::cerr << "error : cannot find file : " << name << ".param\n";
        exit(1);
    }

    is >> dataFile;
    is >> name;
    is >> treeFile;
    is >> correction;
    is >> dgam;
    is >> dbeta;
    is >> asymmbeta;
    is >> every;
    is >> until;
    is >> numChains;
    is >> swapInterval;
    is >> delta;
    is >> sigma;
    is >> saveall;
    int i;
    is >> outgroupFile;
    is >> i;
    modeltype = static_cast<ModelPrior::Type>(i);
    is >> i;
    branchprior = static_cast<BranchPrior::Type>(i);
    is >> i;
    rootprior = static_cast<RootPrior::Type>(i);
    is >> mixture;
    is >> rootmin;
    is >> rootmax;
    is >> nexus;

    is.close();
}


void Biphy::save( void ) {
    remove((name+".param").c_str());
    std::ofstream os((name + ".param").c_str(), std::fstream::app);

    os << dataFile << "\n";
    os << name << "\n";
    os << treeFile << "\n";
    os << correction << "\n";
    os << dgam << "\n";
    os << dbeta << "\n";
    os << asymmbeta << "\n";
    os << every << "\n";
    os << until << "\n";
    os << numChains << "\n";
    os << swapInterval << "\n";
    os << delta << "\n";
    os << sigma << "\n";
    os << saveall << "\n";
    os << outgroupFile << "\n";
    os << modeltype << "\n";
    os << branchprior << "\n";
    os << rootprior << "\n";
    os << mixture << "\n";
    os << rootmin << "\n";
    os << rootmax << "\n";
    os << nexus << "\n";

    os.close();
}


void Biphy::init( void ) {

    precheck();
    readInputFiles();
    printConfiguration();
    initModel();
    initMCMC();
}



void Biphy::precheck( void ) {
    
    RandomNumberGenerator* rng = GLOBAL_RNG;
    rng->setSeed(std::vector<unsigned int>(2,time(NULL)*getpid()));

    /* First, we read in the data */
    // the matrix
    Settings::userSettings().setPrintNodeIndex(false);

}

void Biphy::readInputFiles( void ) {

    //const string format("phylip|standard|sequential");
    data = BinaryCharacterDataReader::getInstance().readMatrix(dataFile);
    
    // data checks
    if(data == NULL || data->getNumberOfTaxa() == 0){
        std::cerr << "Error: failed to read datafile" << std::endl;
        exit(1);
    }

    NewickTreeReader reader;
    if(treeFile != "None"){
        trees = *(reader.readTrees(treeFile));
        /*if(trees.size() != data.size()){
            std::cerr << "Error: number of matrices (" << data.size() << ") does not match number of trees (" << trees.size() << ")\n";
            exit(1);
        }*/
    }

    if(trees.size() > 0)
        std::cout << "Read " << trees.size() << " trees\n";
    
    if(cvfile != "None"){
        cvdata = BinaryCharacterDataReader::getInstance().readMatrix(cvfile);
        std::cout << "computing cross-validation scores for " << cvfile << std::endl;
    
        // data checks
        if(cvdata == NULL || cvdata->getNumberOfTaxa() == 0){
            std::cerr << "Error: failed to read cv datafile" << std::endl;
            exit(1);
        }
    }
    
    //construct constraints
    if(outgroupFile != "None"){
        CladeReader reader(outgroupFile);
        outgroup = reader.getClades()[0];
        std::cout << "Outgroup:\n";
        std::cout << outgroup << std::endl;
        std::cout << std::endl;
    }
}

void Biphy::printConfiguration( void ) {

    if(modeltype > ModelPrior::HOMOGENEOUS){
        std::cout << "branch-wise stationary frequencies:\n";
        switch(modeltype){
            case ModelPrior::MIXTURE:
                std::cout << "pi(i) ~ Beta(beta1,beta2) mixture with " << mixture << " components\n";
                break;
            default:
                std::cout << "pi(i) ~ Beta(beta1,beta2)\n";
                break;
        }
    }
    else if(modeltype == ModelPrior::DOLLO)
    {
        std::cout << "dollo model";
    }
    else if(modeltype == ModelPrior::HOMOGENEOUS)
    {
        std::cout << "homogeneous asymmetric reversible model\n";
    }
    else if(modeltype == ModelPrior::MK)
	{
		std::cout << "Mk model\n";
	}

    std::cout << "\n";

    if(modeltype != ModelPrior::DOLLO){
        if(rootprior == RootPrior::RIGID)
            std::cout << "rigid ";
        if(modeltype > ModelPrior::HOMOGENEOUS)
        {
        	std::cout << "root frequency:\n";
            switch(modeltype){
                case ModelPrior::MIXTURE:
                    std::cout << "phi ~ Beta(beta1,beta2) mixture with " << mixture << " components\n";
                    break;
                default:
                    if(rootprior == RootPrior::TRUNCATED)
                        std::cout << "phi ~ Beta(beta1,beta2) truncated on (" << rootmin << "," << rootmax << ")\n";
                    else
                        std::cout << "phi ~ Beta(beta1,beta2)\n";
                    break;
            }
            std::cout << "\nbeta1, beta2 ~ exponential of mean 1\n";
        }
        else if(rootprior == RootPrior::TRUNCATED)
        {
        	std::cout << "root frequency:\n";
        	std::cout << "phi ~ Uniform(" << rootmin << "," << rootmax << ")\n";
		}
		else if(modeltype == ModelPrior::HOMOGENEOUS)
		{
			std::cout << "root frequency:\n";
			std::cout << "phi ~ Beta(1,1)\n";
		}
    }

    std::cout << "\n";

    if(branchprior == BranchPrior::EXPONENTIAL)
    {
        std::cout << "branch lengths ~ iid exponential of mean mu\n";
        std::cout << "mu ~ exponential of mean 0.1\n";
    }
    else if(branchprior == BranchPrior::DIRICHLET)
    {
        std::cout << "branch lengths ~ dirichlet with concentration parameter = 1\n";
        std::cout << "tree length ~ exponential of mean 1\n";
    }
    else if(branchprior == BranchPrior::FIXED)
    {
    	std::cout << "fixed branch lengths\n";
    }
    else if(branchprior == BranchPrior::STRICT)
	{
		std::cout << "strict clock rate ~ exponential of mean 1\n";
	}


    if(dgam > 1){
    	std::cout << "\n";
        std::cout << "rates across sites ~ discrete gamma with " << dgam << " rate categories\n";
        std::cout << "\tmean 1 and variance 1/alpha^2\n";
        std::cout << "alpha ~ exponential of mean 1\n";
    }

    if(dbeta > 1){

        std::cout << "\n";

    	if(!asymmbeta)
    	{
			std::cout << "frequencies across sites ~ discrete Beta(beta,beta) with " << dgam << " rate categories\n";
			std::cout << "beta ~ exponential of mean 1\n";
    	}
    	else
    	{
    		std::cout << "frequencies across sites ~ discrete Beta(beta1,beta2) with " << dgam << " rate categories\n";
    		std::cout << "beta1 ~ exponential of mean 1\n";
    		std::cout << "beta2 ~ exponential of mean 1\n";
    	}
	}

    if(correction != AscertainmentBias::ALL){
        std::cout << "\n";
        std::cout << "correction for unobservable site-patterns:\n";
        if(correction == AscertainmentBias::INFORMATIVE){
            std::cout << "\tuninformative sites\n";
        }else{
            if((correction & AscertainmentBias::VARIABLE) == AscertainmentBias::VARIABLE)
                    std::cout << "\tconstant sites\n";
            else if(correction & AscertainmentBias::NOABSENCESITES)
                    std::cout << "\tconstant absence\n";
            else if(correction & AscertainmentBias::NOPRESENCESITES)
                    std::cout << "\tconstant presence\n";

            if((correction & AscertainmentBias::NOSINGLETONS) == AscertainmentBias::NOSINGLETONS)
                std::cout << "\tsingletons\n";
            else if(correction & AscertainmentBias::NOSINGLETONPRESENCE)
                std::cout << "\tsingleton presence\n";
            else if(correction & AscertainmentBias::NOSINGLETONABSENCE)
                std::cout << "\tsingleton absence\n";
        }
    }

    std::cout << "\n";
    
#if defined SSE_ENABLED
    std::cout << "using SSE likelihood calculator\n\n";
#elif defined AVX_ENABLED
    std::cout << "using AVX likelihood calculator\n\n";
#else
    std::cout << "using generic likelihood calculator\n\n";
#endif
}


void Biphy::initModel( void ) {

    // constant nodes
    one = new ConstantNode<double>("one", new double(1.0) );
    zero = new ConstantNode<double>("zero", new double(0.0) );

    // RAS prior
    std::vector<const TypedDagNode<double>* > gamma_rates = std::vector<const TypedDagNode<double>* >();
    //DeterministicNode<std::vector<double> > *site_rates_norm;
    if(dgam > 1){
        shape = new ContinuousStochasticNode("alpha", new ExponentialDistribution(one) );
        for(size_t cat = 0; cat < dgam; cat++){
            std::stringstream name;
            std::stringstream value_name;
            name << "q";
            name << cat+1;
            value_name << name.str() << "_value";
            gamma_rates.push_back( new DeterministicNode<double>(value_name.str(), new QuantileFunction(new ConstantNode<double>(name.str(), new double((cat+1.0/2.0)/dgam) ), new GammaDistribution(shape, shape) ) ));
        }
        DeterministicNode<std::vector<double> >* site_rates_unnorm = new DeterministicNode<std::vector<double> >( "site_rates_unnorm", new VectorFunction<double>(gamma_rates) );
        site_rates = new DeterministicNode<std::vector<double> >( "site_rates", new NormalizeVectorFunction(site_rates_unnorm) );
        //site_rates = new DeterministicNode<std::vector<double> >( "site_rates", new VectorFunction<double>(gamma_rates) );
    }

    // constant nodes
    ConstantNode<double> *ten = new ConstantNode<double>("ten", new double(10.0) );

    //number of branches
    bool rooted = modeltype != ModelPrior::HOMOGENEOUS;
    
    size_t numBranches = 2*data->getNumberOfTaxa() - 2 - !rooted;
    size_t numNodes = numBranches + 1;
    
    size_t w = numNodes > 0 ? (int) log10 ((double) numNodes) + 1 : 1;

    if(modeltype == ModelPrior::MIXTURE)
      w = mixture > 0 ? (int) log10 ((double) mixture) + 1 : 1;

    // tree prior
    std::vector<std::string> names = data->getTaxonNames();
    
    //StochasticNode<TimeTree> *tau = new StochasticNode<TimeTree>( "tau", new UniformConstrainedTimeTreeDistribution(one,names,constraints,outgroup) );
    StochasticNode<Tree> *tau = new StochasticNode<Tree>( "tau", new UniformTopologyDistribution(names, outgroup, rooted) );
    if(trees.size() > 0){
        tau->setValue(trees[0], true);
        tau->setIgnoreRedraw(true);
    }

    // base frequencies hyperprior
    StochasticNode<double> *beta1;
    StochasticNode<double> *beta2;

    // base frequencies prior
    ContinuousStochasticNode* logit;
    TypedDagNode<double> *phi;
    if(modeltype != ModelPrior::DOLLO)
    {
        if(rootprior == RootPrior::TRUNCATED)
        {
            if(modeltype > ModelPrior::HOMOGENEOUS)
            {
                beta1 = new StochasticNode<double>( "beta1", new ExponentialDistribution(one) );
                beta2 = new StochasticNode<double>( "beta2", new ExponentialDistribution(one) );
                phi = new StochasticNode<double>( "phi", new TruncatedDistributionUnnormalized( new BetaDistribution( beta1, beta2 ), new ConstantNode<double>("rootmin", new double(rootmin) ), new ConstantNode<double>("rootmax", new double(rootmax) ) ) );
            }
            else
            {
                phi = new StochasticNode<double>( "phi", new UniformDistribution( new ConstantNode<double>("rootmin", new double(rootmin) ), new ConstantNode<double>("rootmax", new double(rootmax) ) ) );
            }
        }
        else
        {
            if(modeltype > ModelPrior::HOMOGENEOUS)
            {
                beta1 = new StochasticNode<double>( "beta1", new ExponentialDistribution(one) );
                beta2 = new StochasticNode<double>( "beta2", new ExponentialDistribution(one) );
                if(modeltype == ModelPrior::HIERARCHICAL && rootprior == RootPrior::RIGID)
                    phi = new StochasticNode<double>( "phi", new BetaDistribution( beta1,beta2 ) );
            }
            else if(modeltype == ModelPrior::HOMOGENEOUS)
            {
            	/*ConstantNode<double>* var = new ConstantNode<double>("var", new double(1.66378) );
            	logit = new ContinuousStochasticNode("logit", new NormalDistribution(zero, var));
            	phi = new DeterministicNode<double>( "phi", new LogitFunction( logit, true ) );*/

            	phi = new StochasticNode<double>( "phi", new BetaDistribution( one,one ) );
            }
        }
    }

    //mixture vectors
    std::vector<const TypedDagNode<double > *> pi_cats;
    TypedDagNode< std::vector< double > >* pi_mix;
    std::vector<const TypedDagNode<double> *> pi_stat;

    TypedDagNode< std::vector< double > >* pi_vector;
    TypedDagNode< std::vector< double > >* mix_probs;

    const TopologyNode &root = tau->getValue().getRoot();
    size_t rootIndex = root.getIndex();
    size_t left = root.getChild(0).getIndex();
    size_t right = root.getChild(1).getIndex();

    // branch frequency prior
    if(modeltype == ModelPrior::MK && dbeta > 1)
	{
    	if(asymmbeta)
    	{
    		beta1 = new ContinuousStochasticNode("beta1", new ExponentialDistribution(one) );
    		beta2 = new ContinuousStochasticNode("beta2", new ExponentialDistribution(one) );
    	}
    	else
    	{
    		beta1 = new ContinuousStochasticNode("beta", new ExponentialDistribution(one) );
    		beta2 = beta1;
    	}

    	std::vector<double> breaks;
		for(size_t cat = 0; cat < dbeta - 1; cat++){
			breaks.push_back((cat + 1.0)/dbeta);
		}
		pi_vector = new DeterministicNode<std::vector<double> >( "pi", new DiscreteBetaValuesFunction(beta1, beta2, breaks) );
	}
    else if(modeltype == ModelPrior::HIERARCHICAL)
    {
        for (size_t i = 0 ; i < numNodes ; i++ )
        {
            std::ostringstream pi_name;
            pi_name << "pi(" << std::setfill('0') << std::setw(w) << i << ")";
            if( ((i == left || i == right) && rootprior == RootPrior::RIGID) || (i == rootIndex && (rootprior == RootPrior::RIGID || rootprior == RootPrior::TRUNCATED)) )
                pi_stat.push_back(new DeterministicNode<double>(pi_name.str(), new ConstantFunction<double>(phi)));
            else
                pi_stat.push_back(new StochasticNode<double>( pi_name.str(), new BetaDistribution(beta1,beta2) ) );
        }
        pi_vector = new DeterministicNode< std::vector< double > >( "pi", new VectorFunction<double>( pi_stat ) );
    }
    else if(modeltype == ModelPrior::MIXTURE)
	{
		//mix_probs = new ConstantNode<std::vector<double> >("probs", new std::vector<double>(mixture, 1.0/mixture) );
		ConstantNode<std::vector<double> > *conc_mix = new ConstantNode<std::vector<double> >("conc_mix", new std::vector<double>(mixture,1.0) );
		mix_probs = new StochasticNode< std::vector< double > >( "probs", new DirichletDistribution( conc_mix));
		//std::vector<double> breaks;
		//for(size_t cat = 0; cat < mixture - 1; cat++){
		for(size_t cat = 0; cat < mixture; cat++){
			std::ostringstream name;
			name << "pi(" << std::setfill('0') << std::setw(w) << cat << ")";
			pi_cats.push_back(new StochasticNode<double >( name.str(), new BetaDistribution(beta1,beta2) ) );

			//breaks.push_back((cat+1.0)/mixture);
		}
		pi_mix = new DeterministicNode<std::vector<double> >( "pi", new VectorFunction<double>(pi_cats) );
		//pi_mix = new DeterministicNode< std::vector< double > >( "pi", new DiscreteBetaValuesFunction( beta1, beta2, breaks ) );

		if(rootprior == RootPrior::RIGID)
			phi = new StochasticNode<double>( "phi", new MixtureDistribution<double>(pi_mix, mix_probs) );

		for (size_t i = 0 ; i < numNodes ; i++ ) {
			std::ostringstream q_name;
			q_name << "q(" << std::setfill('0') << std::setw(w) << i << ")";
			if( (i == left || i == right || i == rootIndex) && rootprior == RootPrior::RIGID )
				pi_stat.push_back(new DeterministicNode<double>( q_name.str()+"_pi", new ConstantFunction< double >( phi ) ));
			else
				pi_stat.push_back(new StochasticNode<double>( q_name.str()+"_pi", new MixtureDistribution<double>(pi_mix, mix_probs) ));
		}

		pi_vector = new DeterministicNode< std::vector< double > >( "pi_stat", new VectorFunction<double>( pi_stat ) );
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
    if(branchprior == BranchPrior::EXPONENTIAL || branchprior == BranchPrior::STRICT)
    {
        if(branchprior == BranchPrior::STRICT)
        	mu = new ContinuousStochasticNode("mu", new ExponentialDistribution(one) );
        else
        	mu = new ContinuousStochasticNode("mu", new ExponentialDistribution(ten) );

        if(branchprior == BranchPrior::EXPONENTIAL)
        {
        	rec_mu = new DeterministicNode<double>( "rec_mu", new BinaryDivision<double,double,double>(one,mu));
            for (size_t i = 0 ; i < numBranches ; i++ )
            {
                std::ostringstream br_name;
                br_name << "br(" << std::setfill('0') << std::setw(w) << i << ")";
    
                ContinuousStochasticNode* tmp_branch_rate = new ContinuousStochasticNode( br_name.str(), new ExponentialDistribution(rec_mu));
                branchRates.push_back( tmp_branch_rate );
                branchRates_nonConst.push_back( tmp_branch_rate );
            }

            br_vector = new DeterministicNode< std::vector< double > >( "br_vector", new VectorFunction< double >( branchRates ) );
        }
    }
    else if(branchprior == BranchPrior::DIRICHLET)
    {
        ConstantNode<std::vector<double> > *conc = new ConstantNode<std::vector<double> >("conc", new std::vector<double>(numBranches,1.0) );
        tree_length = new StochasticNode<double>("length", new GammaDistribution(one,one));
        br_times = new StochasticNode< std::vector< double > >( "br_times", new DirichletDistribution( conc));
        br_vector = new DeterministicNode< std::vector< double > >( "br_vector", new VectorScaleFunction<double>(br_times,tree_length) );
    }

    // Substitution model
    BinarySubstitutionModel* charModel;
    BinarySubstitutionModel* charModel2;
    
    if(modeltype == ModelPrior::DOLLO)
        charModel = new BinaryDolloSubstitutionModel(tau, AscertainmentBias::Coding(correction));
    else
    {
        charModel = new BinarySubstitutionModel(tau, AscertainmentBias::Coding(correction));
    }
    
    if(cvdata != NULL)
    	charModel->setVerbose(false);

    if(branchprior == BranchPrior::EXPONENTIAL || branchprior == BranchPrior::DIRICHLET)
    {
        charModel->setBranchLengths(br_vector);
    }
    else if(branchprior == BranchPrior::STRICT)
        charModel->setClockRate(mu);
    
    if(modeltype > ModelPrior::HOMOGENEOUS)
        charModel->setStationaryFrequency( pi_vector );
    else if(modeltype == ModelPrior::HOMOGENEOUS)
    {
        charModel->setStationaryFrequency( phi );
    }
    else if(modeltype == ModelPrior::MK && dbeta > 1)
	{
		charModel->setSiteFrequencies( pi_vector );
	}

    if(dgam > 1)
    {
    	charModel->setSiteRates( site_rates );
    }

    StochasticNode< BinaryCharacterData >* charactermodel = new StochasticNode< BinaryCharacterData >("S", charModel );
    charactermodel->clamp( data );

    if(trees.empty())
        moves.push_back( new SubtreePruneRegraft(tau, 0.3848 , rootprior == RootPrior::RIGID) );

    if(branchprior == BranchPrior::EXPONENTIAL || branchprior == BranchPrior::DIRICHLET)
	{
		if(branchprior == BranchPrior::EXPONENTIAL)
		{
			moves.push_back( new ScaleMove(mu, 0.3, true, 0.02) );
			monitoredNodes.push_back( mu );

			tree_length = new DeterministicNode<double >("length", new SumFunction<double>(br_vector) );
			monitoredNodes.push_back( tree_length );

			std::vector<Move*> br_moves;

			for (size_t i = 0 ; i < numBranches ; i ++ )
				moves.push_back( new ScaleMove(branchRates_nonConst[i], 2.0*log(1.6), true, 0.3846/float(numBranches) ) );
		}
		else if(branchprior == BranchPrior::DIRICHLET)
		{
			moves.push_back(new ScaleMove((StochasticNode<double>*)tree_length, 1.0, true, 0.1*(0.3846+0.1346+0.0577)));
			moves.push_back(new SimplexMove(br_times, 100.0, numBranches, 0.0, true, 0.9*(0.3846+0.1346+0.0577)/2));
			moves.push_back(new SimplexMove(br_times, 100.0, 2, 0.0, true, 0.9*(0.3846+0.1346+0.0577)/2));
			monitoredNodes.push_back( tree_length );
		}

	}
	else if(branchprior == BranchPrior::STRICT)
	{
		moves.push_back( new ScaleMove(mu, 0.3, true, 0.02) );
		monitoredNodes.push_back( mu );
	}

    if(modeltype > ModelPrior::HOMOGENEOUS)
	{
		for (size_t i = 0 ; i < numNodes ; i ++ )
		{
			if( (rootprior != RootPrior::TRUNCATED && rootprior != RootPrior::RIGID) ||
				(rootprior == RootPrior::RIGID && !(i == left || i == right || i == rootIndex)) ||
				(rootprior == RootPrior::TRUNCATED && i != rootIndex)
				)
			{
				if(modeltype == ModelPrior::HIERARCHICAL)
					moves.push_back( new BetaSimplexMove((StochasticNode<double>*)pi_stat[i], 100.0, true, 0.5*0.0096/float(numBranches) ) );
				else if(modeltype == ModelPrior::MIXTURE)
					moves.push_back( new MixtureAllocationMove<double>((StochasticNode<double>*)pi_stat[i], 0.02 ) );
			}
		}
	}

    if(modeltype == ModelPrior::HOMOGENEOUS || (modeltype > ModelPrior::HOMOGENEOUS && (rootprior == RootPrior::TRUNCATED || rootprior == RootPrior::RIGID)) )
	{
		if(modeltype == ModelPrior::MIXTURE)
			moves.push_back( new MixtureAllocationMove<double>((StochasticNode<double>*)phi, 0.02 ) );
		else
			moves.push_back( new BetaSimplexMove((StochasticNode<double>*)phi, 1000.0, true, 0.02 ) );

		if(modeltype == ModelPrior::HOMOGENEOUS || rootprior == RootPrior::RIGID || rootprior == RootPrior::TRUNCATED)
			monitoredNodes.push_back( phi );
	}

    DeterministicNode<double>* meanpi;
    if(modeltype == ModelPrior::HIERARCHICAL)
    {
        meanpi = new DeterministicNode<double >( "mean_pi", new MeanFunction(pi_vector) );
        monitoredNodes.push_back( meanpi );
    }

    if(modeltype > ModelPrior::HOMOGENEOUS || (modeltype == ModelPrior::MK && dbeta > 1))
    {
        moves.push_back( new ScaleMove(beta1, 0.3, true, 0.0096/4) );
        monitoredNodes.push_back( beta1 );

        if(modeltype > ModelPrior::HOMOGENEOUS || asymmbeta)
        {
        	moves.push_back( new ScaleMove(beta2, 0.3, true, 0.0096/4) );
        	monitoredNodes.push_back( beta2 );
        }

        if(modeltype == ModelPrior::MIXTURE){
			for (size_t i = 0 ; i < mixture; i ++ ) {
				moves.push_back( new BetaSimplexMove((StochasticNode<double>*)pi_cats[i], 1000.0, true, 0.02 ) );
			}
			moves.push_back( new SimplexMove((StochasticNode<std::vector<double> >*)mix_probs, 10.0, mixture, 0.0, true, 0.02 ) );

			//monitoredNodes.push_back(pi_mix);
			//monitoredNodes.push_back(mix_probs);
			monitoredNodes.push_back(new DeterministicNode<double >( "pi_ent", new EntropyFunction(pi_mix) ));
			monitoredNodes.push_back(new DeterministicNode<double >( "mix_ent", new EntropyFunction(mix_probs) ));
		}
    }
    
    if(dgam > 1)
    {
        moves.push_back( new ScaleMove(shape, 0.2, true, 0.02) );
        monitoredNodes.push_back( shape );
    }

    bool useParallelMcmcmc = (numChains > 1);

    if(readstream)
    {
        //tau->getValue().getTreeChangeEventHandler().clear();
        
        if(ppred)
            monitors.push_back( new PosteriorPredictiveStateFrequencyMonitor( charactermodel, every, name+".ppred", useParallelMcmcmc) );
        if(cvdata != NULL)
            monitors.push_back( new CrossValidationScoreMonitor( charactermodel, cvdata, every, name+".cv", useParallelMcmcmc) );
        if(perSiteLnProbs)
            monitors.push_back( new PerSiteLnProbMonitor( charactermodel, every, name+".lnprobs", useParallelMcmcmc) );
        if(ancestral)
            monitors.push_back( new MappingMonitor( charactermodel, every, name+".mapping", useParallelMcmcmc) );
        //if(dolloMapping)
        //  monitors.push_back( new BinaryDolloCompatibleMonitor<BranchLengthTree>( charactermodel, every, name+".dollo.fa") );
    }
    else
    {
        if(trees.empty() || branchprior == BranchPrior::EXPONENTIAL || branchprior == BranchPrior::DIRICHLET)
        {
            monitors.push_back( new NewickTreeMonitor( tau, every, name+".treelist", useParallelMcmcmc || restart) );
            if(nexus)
            {
                std::set<TypedDagNode<std::vector<double> > *> piset;
                if(modeltype > ModelPrior::HOMOGENEOUS)
                    piset.insert(pi_vector);
                monitors.push_back( new NexusTreeMonitor( tau, piset, every, name+".treelist.nex", useParallelMcmcmc || restart) );
                //monitors.push_back( new ExtendedNewickTreeMonitor( tau, piset, every, name+".treelist", "\t", false, false, false, useParallelMcmcmc) );
            }
        }
    }
    
}


void Biphy::initMCMC( void ) {
    model = new Model(one);
    std::cout << "model okay\n";

    bool useParallelMcmcmc = (numChains > 1);

    if(!readstream)
        monitors.push_back( new FileMonitor( monitoredNodes, every, name+".trace", "\t", false, true, false, useParallelMcmcmc || restart, useParallelMcmcmc, false ) );

    double startingHeat = 1.0;
    mcmc = new ParallelMcmcmc(model, moves, monitors, name+".stream", "single", every, numChains, numChains, swapInterval, delta, sigma, startingHeat, saveall);
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

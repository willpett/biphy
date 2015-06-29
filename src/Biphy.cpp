#include "Biphy.h"

#include "BetaDistribution.h"
#include "BetaSimplexMove.h"
#include "ConstantFunction.h"
#include "DeterministicNode.h"
#include "ExponentialDistribution.h"
#include "FileMonitor.h"
#include "GammaDistribution.h"
#include "NclReader.h"
#include "NewickTreeReader.h"
#include "NormalizeVectorFunction.h"
#include "RandomNumberFactory.h"
#include "RbSettings.h"
#include "RbStatisticsHelper.h"
#include "QuantileFunction.h"
#include "ScaleMove.h"
#include "VectorFunction.h"

#include <unistd.h>

using namespace RevBayesCore;

Biphy::Biphy(const std::string n,
				const std::string d,
				const std::string t,
				int c,
				int dg,
				int e,
				int u,
				int num,
				int swa,
				double de,
				double si,
				bool sav) :
		dataFile( d ),
		name( n ),
		treeFile( t ),
		correction(c),
		dgam( dg),
		every( e ),
		until( u ),
		numChains( num),
		swapInterval( swa ),
		delta( de ),
		sigma( si ),
		saveall( sav),
		readstream(false),
		restart(false)
{
    save();
}

Biphy::Biphy(const std::string name) :
		name( name ), readstream(false), restart(true)
{
	open();
}

void Biphy::open( void ) {
	std::ifstream is((name + ".param").c_str());
	if (!is)        {
		std::cerr << "error : cannot find file : " << name << ".param\n";
		exit(1);
	}

	openParams(is);

	is.close();
}

void Biphy::openParams( std::ifstream& is ) {

	is >> dataFile;
	is >> name;
	is >> treeFile;
	is >> correction;
	is >> missing;
	is >> dgam;
	is >> every;
	is >> until;
	is >> numChains;
	is >> swapInterval;
	is >> delta;
	is >> sigma;
	is >> saveall;
}

void Biphy::save( void ) {
	remove((name+".param").c_str());
	std::ofstream os((name + ".param").c_str(), std::fstream::app);

	saveParams(os);

	os.close();
}

void Biphy::saveParams( std::ofstream& os ) {
	os << dataFile << "\n";
	os << name << "\n";
	os << treeFile << "\n";
	os << correction << "\n";
	os << missing << "\n";
	os << dgam << "\n";
	os << every << "\n";
	os << until << "\n";
	os << numChains << "\n";
	os << swapInterval << "\n";
	os << delta << "\n";
	os << sigma << "\n";
	os << saveall << "\n";
}

void Biphy::init( void ) {

	precheck();
	readInputFiles();
	printConfiguration();
	initModel();
	initMCMC();
}



void Biphy::precheck( void ) {
    
	RevBayesCore::RandomNumberGenerator* rng = RevBayesCore::GLOBAL_RNG;
	rng->setSeed(std::vector<unsigned int>(2,time(NULL)*getpid()));

    /* First, we read in the data */
    // the matrix
	RbSettings::userSettings().setPrintNodeIndex(false);

}

void Biphy::readInputFiles( void ) {

	//const string format("phylip|standard|sequential");
	data = NclReader::getInstance().readMatrices(dataFile);

    // data checks
    if(data.empty()){
		//std::cerr << "Error: failed to read datafile" << std::endl;
		exit(1);
	}

    if(data[0]->getDatatype() != "Standard"){
    	std::cerr << "Error: incompatible datatype '" << data[0]->getDatatype() << "'" << std::endl;
    	exit(1);
    }

	size_t numTaxa = data.front()->getNumberOfTaxa();
	size_t excluded = 0;
	size_t numSites = 0;
    for(size_t i = 0; i < data.size(); i++){
		if(data[i]->getDatatype() != "Standard"){
			std::cerr << "Error: incompatible datatype '" << data[i]->getDatatype() << "' for datafile " << dataFile << std::endl;
			exit(1);
		}
		std::map<size_t,size_t> cts;
		//std::cout << "Read " << data[i]->getNumberOfCharacters() << " characters from datafile " << dataFile << std::endl;
		for(size_t c = 0; c < data[i]->getNumberOfCharacters(); c++){
			size_t present = 0;
			for(size_t t = 0; t < data[i]->getNumberOfTaxa(); t++){
				StandardState state = ((DiscreteCharacterData<StandardState> *)data[i])->getCharacter(t,c);
				missing = missing || state.isGapState();
				present += size_t(state.getStringValue() == "1");
			}
			cts[present]++;
			if( ((present == 0) 			&& (correction & RevBayesCore::NO_ABSENT_SITES)) ||
				((present == numTaxa) 		&& (correction & RevBayesCore::NO_PRESENT_SITES)) ||
				((present == 1) 			&& (correction & RevBayesCore::NO_SINGLETON_GAINS)) ||
				((present == numTaxa - 1)	&& (correction & RevBayesCore::NO_SINGLETON_LOSSES))
			){
				data[i]->excludeCharacter(c);
				excluded++;
			}
		}
		//for(std::map<size_t,size_t>::iterator it = cts.begin(); it != cts.end(); it++)
		//	std::cerr << it->first << "\t" << it->second << std::endl;
		numSites += data[i]->getNumberOfIncludedCharacters();
	}

    if(excluded > 0)
    	std::cout << "excluded " << excluded << " characters" << std::endl;

    symbols = ((DiscreteCharacterData<StandardState> *)data[0])->getCharacter(0,0).getStateLabels();
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

	NewickTreeReader reader;
    if(treeFile != "None"){
    	trees = *(reader.readBranchLengthTrees(treeFile));
    	if(trees.size() != data.size()){
			std::cerr << "Error: number of matrices (" << data.size() << ") does not match number of trees (" << trees.size() << ")\n";
			exit(1);
		}
    }

    if(trees.size() > 0)
    	std::cout << "Read " << trees.size() << " trees\n";
}

void Biphy::printConfiguration( void ) {

	std::cout << "\n";

	if(dgam > 1){
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
		std::cout << "\nproportion of missing data ~ Beta(1,1)\n";

	std::cout << "\n";
}

void Biphy::initModel( void ) {

	// constant nodes
	one = new ConstantNode<double>("one", new double(1.0) );
	zero = new ConstantNode<double>("zero", new double(0.0) );

	// RAS prior
	std::vector<const TypedDagNode<double>* > gamma_rates = std::vector<const TypedDagNode<double>* >();
	//DeterministicNode<std::vector<double> > *site_rates_norm;
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
		//DeterministicNode<std::vector<double> >* site_rates_unnorm = new DeterministicNode<std::vector<double> >( "site_rates_unnorm", new VectorFunction<double>(gamma_rates) );
		//site_rates = new DeterministicNode<std::vector<double> >( "site_rates", new NormalizeVectorFunction(site_rates_unnorm) );
		site_rates = new DeterministicNode<std::vector<double> >( "site_rates", new VectorFunction<double>(gamma_rates) );
	}

	if(missing)
		xi = new StochasticNode<double>( "xi", new BetaDistribution( one,one ) );

}

void Biphy::initMCMC( void ) {
	model = new Model(one);
	std::cout << "model okay\n";

	if(dgam > 1){
		moves.push_back( new ScaleMove(shape, 1.0, true, 1.0) );
		monitoredNodes.push_back( shape );
		//monitoredNodes.push_back(site_rates_norm);
	}

	if(missing){
		moves.push_back( new BetaSimplexMove(xi, 1.0, true, 2.0 ) );
		monitoredNodes.push_back(xi);
	}

	bool useParallelMcmcmc = (numChains > 1);

	if(!readstream)
		monitors.push_back( new FileMonitor( monitoredNodes, every, name+".trace", "\t", false, true, false, useParallelMcmcmc || restart, useParallelMcmcmc, false ) );

	double startingHeat = 1.0;
	mcmc = new ParallelMcmcmc(*model, moves, monitors, name+".stream", "random", every, numChains, numChains, swapInterval, delta, sigma, startingHeat, saveall);
	std::cout << "mcmc okay\n";
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

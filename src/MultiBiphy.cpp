#include <cstdlib>
#include <iostream>
#include "MultiBiphy.h"
    
MultiBiphy::MultiBiphy(const std::string df,
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
		dataFile( df ),
		name( n ),
		treeFile( t ),
		modeltype(mt),
		A_prior(Ap),
		B_prior(Bp),
		A_species(As),
		B_species(Bs),
		correction(c),
		dgam( d),
		every( e ),
		until( u ),
		numChains( num),
		swapInterval( swa ),
		delta( de ),
		sigma( si ),
		saveall( sav),
		ancestral( anc),
		readstream(false),
		restart(false)
{
}

MultiBiphy::MultiBiphy(const std::string n, bool anc) :
		name( n ),
		readstream(false),
		restart(true),
		ancestral(anc)
{
	open();

	if(anc){
		ancestral = anc;
		every = 1;
		readstream = true;
		restart = false;
	}
}

MultiBiphy::MultiBiphy(const std::string n) :
		name( n ), readstream(false), restart(true)
{

}

void MultiBiphy::open( void ) {
	std::ifstream is((name + ".param").c_str());
	if (!is)        {
		std::cerr << "error : cannot find file : " << name << ".param\n";
		exit(1);
	}
	int i;

	is >> dataFile;
	is >> name;
	is >> treeFile;
	is >> i;
	modeltype = static_cast<ModelType>(i);
	is >> i;
	A_prior = static_cast<RatePrior>(i);
	is >> i;
	A_species = static_cast<RatePrior>(i);
	is >> i;
	B_prior = static_cast<RatePrior>(i);
	is >> i;
	B_species = static_cast<RatePrior>(i);
	is >> correction;
	is >> dgam;
	is >> every;
	is >> until;
	is >> numChains;
	is >> swapInterval;
	is >> delta;
	is >> sigma;
	is >> saveall;
	is >> ancestral;

	is.close();
}

void MultiBiphy::save( void ) {
	std::ofstream os((name + ".param").c_str());

	os << dataFile << "\n";
	os << name << "\n";
	os << treeFile << "\n";
	os << modeltype << "\n";
	os << A_prior << "\n";
	os << A_species << "\n";
	os << B_prior << "\n";
	os << B_species << "\n";
	os << correction << "\n";
	os << dgam << "\n";
	os << every << "\n";
	os << until << "\n";
	os << numChains << "\n";
	os << swapInterval << "\n";
	os << delta << "\n";
	os << sigma << "\n";
	os << saveall << "\n";
	os << ancestral << "\n";

	os.close();
}

void MultiBiphy::run( void ) {
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

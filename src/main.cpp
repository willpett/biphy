#include "TestBranchHeterogeneousBinaryModel.h"
#include "RbException.h"
#include <iostream>
#include <stdlib.h>
#include <fstream>

bool fexists(const std::string& filename) {
  std::ifstream ifile(filename.c_str());
  return ifile;
}

void our_terminate (void);

namespace {
    static const bool SET_TERMINATE = std::set_terminate(our_terminate);
}

void our_terminate (void) { // try 1
    try { throw; }
    catch (RbException& e) {
    	std::cout << "RbException:\t" << e.getMessage() << '\n';
    }
    exit(1);
}


int main (int argc, const char * argv[])
{
		std::string datafile = "None";
		std::string treefile = "None";
		std::string cvfile = "None";
		std::string outgroupfile = "None";
		std::string name = "";

		int numChains = 1;
		int swapInterval = 1;
		double deltaTemp = 0.1;
		double sigmaTemp = 1;
		int branchprior = 0;
		int heterogeneous = 0;
		bool dollo = false;
		int mixture = 0;
		bool overwrite = false;
		bool ppred = false;
		bool rigidroot = false;
		bool rootprior = false;
		bool ras = false;
		double rootmin = 0.0;
		double rootmax = 1.0;
		bool saveall = false;
		bool nexus = false;

		int every = 1;
		int until = -1;
		try	{

			if (argc <= 2)	{
				throw(0);
			}

			int i = 1;
			while (i < argc)	{
				std::string s = argv[i];

				if (s == "-d")	{
					i++;
					datafile = argv[i];
				}else if (s == "-f")	{
					overwrite = true;
				}
				else if ((s == "-t") || (s == "-T"))	{
					i++;
					treefile = argv[i];
				}
				else if (s == "-cv")	{
					i++;
					cvfile = argv[i];
				}
				else if (s == "-o")	{
					i++;
					outgroupfile = argv[i];
				}else if (s == "-s"){
					saveall = true;
				}else if (s == "-ras"){
					ras = true;
				}else if (s == "-rr"){
					rigidroot = true;
				}
				else if (s == "-nh")	{
					heterogeneous = 1;
				}
				else if (s == "-e")	{
					nexus = true;
				}else if (s == "-dollo"){
					dollo = true;
				}else if (s == "-h")	{
					heterogeneous = 0;
				}else if (s == "-dpp")	{
					heterogeneous = 2;
				}else if (s == "-dir")	{
					branchprior = 1;
				}else if (s == "-n")	{
					i++;
					numChains = atoi(argv[i]);
				}else if (s == "-m")	{
					i++;
					mixture = atoi(argv[i]);
					if(mixture > 1){
						heterogeneous = 3;
					}else{
						mixture = 0;
					}
				}
				else if (s == "-delta")	{
					i++;
					deltaTemp = atof(argv[i]);
				}
				else if (s == "-sigma")	{
					i++;
					sigmaTemp = atof(argv[i]);
				}
				else if (s == "-si")	{
					i++;
					swapInterval = atoi(argv[i]);
				}else if (s == "-x")	{
					i++;
					if (i == argc) throw(0);
					every = atoi(argv[i]);
					i++;
					if (i == argc) throw(0);
					until = atoi(argv[i]);
				}else if (s == "-rp")	{
					i++;
					if (i == argc) throw(0);
					rootmin = atof(argv[i]);
					i++;
					if (i == argc) throw(0);
					rootmax = atof(argv[i]);
					rootprior = true;
				}else if (s == "-ppred")	{
					ppred = true;
				}else{
					if (i != (argc -1))	{
						throw(0);
					}
					name = argv[i];
				}
				i++;
			}

			if(name == "")
				throw(0);
			if(mixture > 1 && !heterogeneous)
				throw(0);
			if(rootprior && (heterogeneous > 1 || dollo))
				throw(0);
			if(heterogeneous && dollo)
				throw(0);
			if(rigidroot && !heterogeneous)
				throw(0);
		}
		catch(...)	{
			std::cerr << "biphy version 1.0\n";
			std::cerr << '\n';
			std::cerr << "usage: biphy -d <alignment> [-x <every> <until>] <chainname>\n\n";
			std::cerr << "Model options:\n";
			std::cerr << "\t-dollo\t\tdollo model\n";
			std::cerr << "\t-h\t\ttime-homogeneous binary substitution model (default)\n";
			std::cerr << "\t-nh\t\ttime-heterogeneous hierarchical beta model\n";
			std::cerr << "\t-m <int>\ttime-heterogeneous mixture model with <int> components\n";
			std::cerr << "\t-dpp\t\tdirichlet process prior on branch frequencies\n";
			std::cerr << "\t-ras\t\tdiscrete gamma rates across sites model\n";
			std::cerr << "\t-dir\t\tcompound dirichlet branch length prior\n\n";
			std::cerr << "Optional constraints:\n";
			std::cerr << "\t-t <file>\tfixed tree filename\n";
			std::cerr << "\t-o <file>\toutgroup clade file\n";
			std::cerr << "\t-rp <min> <max>\ttruncate root frequency prior (-nh or -h only)\n";
			std::cerr << "\t-rr\t\trigid root frequency (heterogeneous models only)\n\n";
			std::cerr << "MCMCMC options:\n";
			std::cerr << "\t-n <int>\tnumber of chains (default = 1)\n";
			std::cerr << "\t-delta <float>\t(default = 0.1)\n";
			std::cerr << "\t-sigma <float>\t(default = 1)\n\n";
			std::cerr << "Other options:\n";
			std::cerr << "\t-s\t\tsave entire output (default: disabled)\n";
			std::cerr << "\t-e\t\tsave nexus treefile output (default: disabled)\n";
			std::cerr << "\t-ppred\t\tposterior predictive simulation of tip frequencies\n";
			std::cerr << "\t-cv <file>\tcross-validation test alignment\n";
			exit(1);
		}

	RevBayesCore::TestBranchHeterogeneousBinaryModel *chain;

	if(datafile != "None"){
		if(datafile.at(0) != '.' && datafile.at(0) != '/'){
			datafile = "./"+datafile;
		}
		if(treefile != "None" && treefile.at(0) != '.' && treefile.at(0) != '/'){
			treefile = "./"+treefile;
		}
		if(outgroupfile != "None" && outgroupfile.at(0) != '.' && outgroupfile.at(0) != '/'){
			outgroupfile = "./"+outgroupfile;
		}
		if(cvfile != "None" && cvfile.at(0) != '.' && cvfile.at(0) != '/'){
			cvfile = "./"+cvfile;
		}

		if(fexists(name+".param") && !overwrite){
			std::cerr << "chain '" << name << "' exists. use overwrite option -f\n";
			exit(1);
		}else if(overwrite){
			if(saveall)
				remove((name+".chain").c_str());
			remove((name+".param").c_str());
			remove((name+".trace").c_str());
			remove((name+".treelist").c_str());
			if(nexus)
				remove((name+".treelist.nex").c_str());
			if(ppred)
				remove((name+".ppred").c_str());
			if(cvfile != "None")
				remove((name+".cv").c_str());
		}

		chain = new RevBayesCore::TestBranchHeterogeneousBinaryModel(datafile,name,treefile,outgroupfile,branchprior,ras,heterogeneous,dollo,mixture,rigidroot,rootprior,rootmin,rootmax,every,until,numChains,swapInterval,deltaTemp,sigmaTemp,saveall,nexus);
	}else{
		if(!fexists(name+".chain")){
			std::cerr << "chain '" << name << "' does not exist\n";
			exit(1);
		}
		chain = new RevBayesCore::TestBranchHeterogeneousBinaryModel(name,cvfile,ppred);
	}

	try
	{
    	chain->run();
	}
	catch (RbException& e)
	{
		std::cout << "Error:\t" << e.getMessage() << '\n';
		exit(1);
	}
    
    return 0;
}

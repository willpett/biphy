#include "TestBranchHeterogeneousBinaryModel.h"
#include <iostream>
#include <stdlib.h>
#include <fstream>

bool fexists(const std::string& filename) {
  std::ifstream ifile(filename.c_str());
  return ifile;
}

int main (int argc, const char * argv[])
{
	std::cerr << "biphy version 1.0\n";
	std::cerr << '\n';

	std::string datafile = "";
	std::string treefile = "None";
	std::string cvfile = "None";
	std::string outgroupfile = "None";
	std::string name = "";

		int numChains = 1;
		int swapInterval = 1;
		double deltaTemp = 0.1;
		double sigmaTemp = 1;
		bool heterogeneous = false;
		bool overwrite = false;
		bool ppred = false;
		bool rootprior = false;
		double rootmin = 0.0;
		double rootmax = 1.0;

		int every = 1;
		int until = -1;
		try	{

			if (argc == 1)	{
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
				}
				else if (s == "-nh")	{
					heterogeneous = true;
				}
				else if (s == "-h")	{
					heterogeneous = false;
				}else if (s == "-n")	{
					i++;
					numChains = atoi(argv[i]);
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
		}
		catch(...)	{
			std::cerr << "biphy -d <alignment> [-x <every> <until>] <chainname>\n\n";
			std::cerr << "Model options:\n";
			std::cerr << "\t-h\t\ttime-homogeneous binary substitution model (default)\n";
			std::cerr << "\t-nh\t\tnon-homogeneous binary substitution model\n";
			std::cerr << "Optional constraints:\n";
			std::cerr << "\t-t <file>\tfixed tree filename\n";
			std::cerr << "\t-o <file>\toutgroup clade file\n";
			std::cerr << "\t-rp <min> <max>\ttruncate root frequency prior\n\n";
			std::cerr << "MCMCMC options:\n";
			std::cerr << "\t-n <int>\tnumber of chains (default = 1)\n";
			std::cerr << "\t-delta <float>\t(default = 0.1)\n";
			std::cerr << "\t-sigma <float>\t(default = 1)\n\n";
			std::cerr << "Other options:\n";
			std::cerr << "\t-ppred\t\tposterior predictive simulation of tip frequencies\n";
			std::cerr << "\t-cv <file>\tcross-validation test alignment\n";
			exit(1);
		}

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

	if(fexists(name+".trace") && !overwrite){
		std::cerr << "chain '" << name << "' exists. use overwrite option -f\n";
		exit(1);
	}else if(overwrite){
		remove((name+".trace").c_str());
		remove((name+".treelist").c_str());
		if(ppred)
			remove((name+".ppred").c_str());
		if(cvfile != "None")
			remove((name+".cv").c_str());
	}
    RevBayesCore::TestBranchHeterogeneousBinaryModel t = RevBayesCore::TestBranchHeterogeneousBinaryModel(datafile,name,treefile,outgroupfile,cvfile,heterogeneous,ppred,rootprior,rootmin,rootmax,every,until,numChains,swapInterval,deltaTemp,sigmaTemp);
    t.run();
    
    return 0;
}

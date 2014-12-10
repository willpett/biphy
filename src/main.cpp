#include "RbException.h"
#include <iostream>
#include <stdlib.h>
#include <fstream>
#include "util.h"
#include "Biphy.h"

using namespace std;

int main (int argc, const char * argv[])
{
	string datafile = "None";
	string treefile = "None";
	string cvfile = "None";
	string outgroupfile = "None";
	string name = "";

	Biphy::ModelType modeltype = Biphy::HOMOGENEOUS;
	int mixture = 0;

	Biphy::BranchPrior branchprior = Biphy::EXPONENTIAL;

	int dgam = 4;

	Biphy::RootPrior rootprior = Biphy::FREE;
	double rootmin = 0.0;
	double rootmax = 1.0;

	int correction = RevBayesCore::NONE;

	int every = 1;
	int until = -1;
	bool overwrite = false;

	int numChains = 1;
	int swapInterval = 1;
	double delta = 0.05;
	double sigma = 1;

	bool saveall = false;
	bool nexus = false;
	bool ppred = false;
	bool dolloMapping = false;

	Biphy *chain = NULL;

	try	{
		if (argc == 2)	{
			name = argv[1];
			if (name == "-help" || name == "--help" || name == "-h")	{
				throw(0);
			}
			if(!fexists(name+".stream")){
				cerr << "run '" << name << "' does not exist\n";
				exit(1);
			}

			chain = new Biphy(name);
		}else if (argc == 1){
			throw(0);
		}else{

			int i = 1;
			while (i < argc)	{
				string s = argv[i];

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
				}else if (s == "-dgam"){
					i++;
					if (i == argc)	{
						cerr << "error in command: -dgam <int>\n\n";
						exit(1);
					}
					s = argv[i];
					if (! IsInt(s))	{
						cerr << "error in command: -dgam <int>\n\n";
						exit(1);
					}
					dgam = atoi(argv[i]);
				}else if (s == "-rigid"){
					rootprior = Biphy::RIGID;
				}
				else if (s == "-nh")	{
					modeltype = Biphy::HIERARCHICAL;
				}
				else if (s == "-e")	{
					nexus = true;
				}else if (s == "-dollo"){
					modeltype = Biphy::DOLLO;
				}else if (s == "-map"){
					dolloMapping = true;
				}else if (s == "-h")	{
					modeltype = Biphy::HOMOGENEOUS;
				}else if (s == "-dpp")	{
					modeltype = Biphy::DPP;
				}else if (s == "-ldir")	{
					branchprior = Biphy::DIRICHLET;
				}else if (s == "-lexp")	{
					branchprior = Biphy::EXPONENTIAL;
				}else if (s == "-n")	{
					i++;
					if (i == argc)	{
						cerr << "error in command: -n <int>\n\n";
						exit(1);
					}
					s = argv[i];
					if (! IsInt(s))	{
						cerr << "error in command: -n <int>\n\n";
						exit(1);
					}
					numChains = atoi(argv[i]);
				}else if (s == "-u")	{
					i++;
					if (i == argc)	{
						cerr << "error in command: -u <int>\n\n";
						exit(1);
					}
					s = argv[i];
					if (! IsInt(s))	{
						cerr << "error in command: -u <int>\n\n";
						exit(1);
					}
					correction = static_cast<RevBayesCore::CorrectionType>(atoi(argv[i]));
					if(correction < 0 || correction > 15)	{
						cerr << "error in command: -u <int>\n\n";
						exit(1);
					}
				}else if (s == "-m")	{
					i++;
					if (i == argc)	{
						cerr << "error in command: -m <int>\n\n";
						exit(1);
					}
					s = argv[i];
					if (! IsInt(s))	{
						cerr << "error in command: -m <int>\n\n";
						exit(1);
					}
					mixture = atoi(argv[i]);
					if(mixture > 1){
						modeltype = Biphy::MIXTURE;
					}else{
						modeltype = Biphy::HOMOGENEOUS;
						mixture = 0;
					}
				}
				else if (s == "-delta")	{
					i++;
					if (i == argc)	{
						cerr << "error in command: -delta <float>\n\n";
						exit(1);
					}
					s = argv[i];
					if (! IsFloat(s))	{
						cerr << "error in command: -delta <float>\n\n";
						exit(1);
					}
					delta = atof(argv[i]);
				}
				else if (s == "-sigma")	{
					i++;
					if (i == argc)	{
						cerr << "error in command: -sigma <float>\n\n";
						exit(1);
					}
					s = argv[i];
					if (! IsFloat(s))	{
						cerr << "error in command: -sigma <float>\n\n";
						exit(1);
					}
					sigma = atof(argv[i]);
				}
				else if (s == "-si")	{
					i++;
					if (i == argc)	{
						cerr << "error in command: -si <int>\n\n";
						exit(1);
					}
					s = argv[i];
					if (!IsInt(s))	{
						cerr << "error in command: -si <int>\n\n";
						exit(1);
					}
					swapInterval = atof(argv[i]);
				}else if (s == "-x")	{
					i++;
					if (i == argc)	{
						cerr << "error in command: -x <every> <until>\n\n";
						exit(1);
					}
					s = argv[i];
					if (! IsInt(s))	{
						cerr << "error in command: -x <every> <until>\n\n";
						exit(1);
					}
					every = atoi(argv[i]);
					i++;
					if (i != argc)	{
						s = argv[i];
						if (IsInt(s))	{
							until = atoi(argv[i]);
						}
						else	{
							i--;
						}
					}
				}else if (s == "-rp")	{
					i++;
					if (i == argc)	{
						cerr << "error in command: -rp <min> <max>\n\n";
						exit(1);
					}
					s = argv[i];
					if (! IsFloat(s))	{
						cerr << "error in command: -rp <min> <max>\n\n";
						exit(1);
					}
					rootmin = atof(argv[i]);
					i++;
					if (i == argc)	{
						cerr << "error in command: -rp <min> <max>\n\n";
						exit(1);
					}
					s = argv[i];
					if (! IsFloat(s))	{
						cerr << "error in command: -rp <min> <max>\n\n";
						exit(1);
					}
					rootmax = atof(argv[i]);
					if(rootmin < 0 || rootmax > 1.0 || rootmin > rootmax){
						cerr << "error in command: -rp <min> <max>\n\n";
						exit(1);
					}
					rootprior = Biphy::TRUNCATED;
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

			if(name == ""){
				cerr << "error: could not determine run name\n\n";
				exit(1);
			}

			if(mixture > 1 && modeltype != Biphy::MIXTURE)
				throw(0);

			if(rootprior == Biphy::TRUNCATED && (modeltype > Biphy::HIERARCHICAL || modeltype == Biphy::DOLLO)){
				cerr << "error: truncated root frequency only applies to -h or -nh\n\n";
				exit(1);
			}
			if(rootprior == Biphy::RIGID && (modeltype < Biphy::HIERARCHICAL))
				rootprior = Biphy::FREE;

			if(modeltype == Biphy::DOLLO)
				correction = (RevBayesCore::CorrectionType)(correction | RevBayesCore::NO_ABSENT_SITES);
		}
	}
	catch(RbException &e){
		throw(e);
	}catch(...){

		cerr << "biphy " << VERSION << "\n";
		cerr << '\n';
		cerr << "usage: biphy -d <data file> [-x <every> [<until>] ] <run name>\n";

		cerr << "\nBranch frequency prior:\n";
		cerr << "\t-dollo\t\tstochastic dollo model (enables -u 1)\n";
		cerr << "\t-h\t\ttime-homogeneous beta model (default)\n";
		cerr << "\t-nh\t\thierarchical beta model\n";
		cerr << "\t-m <int>\tmixture with <int> beta components\n";
		cerr << "\t-dpp\t\tdirichlet process prior\n";

		cerr << "\nBranch length prior:\n";
		cerr << "\t-lexp\t\thierarchical exponential prior (default)\n";
		cerr << "\t-ldir\t\tcompound dirichlet prior\n";

		cerr << "\nRates across sites prior:\n";
		cerr << "\t-dgam <int>\tdiscrete gamma model with <int> categories (default: 4)\n";
		cerr << "\t\t\t0 or 1 = constant rates model\n";

		cerr << "\nCorrections for unobservable site patterns:\n";
		cerr << "\t-u <int>\twhere <int> is one of:\n";
		cerr << "\t\t0 = no site patterns have been omitted (default)\n";
		cerr << "\t\t1 = constant absence sites have been removed\n";
		cerr << "\t\t2 = constant presence sites have been removed\n";
		cerr << "\t\t4 = singleton gains have been removed\n";
		cerr << "\t\t8 = singleton losses have been removed\n";

		cerr << "\tcombinations are achieved by adding the above values\n";
		cerr << "\te.g.  3 = constant sites have been removed\n";
		cerr << "\t     15 = uninformative sites have been removed\n";

		cerr << "\nOptional constraints:\n";
		cerr << "\t-t <file>\tfixed tree filename\n";
		cerr << "\t-o <file>\toutgroup clade file\n";
		cerr << "\t-rigid\trigid root frequency for heterogeneous models\n";
		cerr << "\t-rp <min> <max>\ttruncate root frequency prior on (min,max) (-nh or -h only)\n";

		cerr << "\nMCMCMC options:\n";
		cerr << "\t-n <int>\tnumber of chains (default = 1)\n";
		cerr << "\t-si <int>\tchain swap interval (default = 1)\n";
		cerr << "\t-delta <float>\t(default = 0.1)\n";
		cerr << "\t-sigma <float>\t(default = 1)\n";

		cerr << "\nOutput options:\n";
		cerr << "\t-s\t\tsave entire output (default: disabled)\n";
		cerr << "\t-e\t\tsave nexus treefile (default: disabled)\n";

		cerr << "\nModel-checking options:\n";
		cerr << "\t-ppred\t\tposterior predictive simulation of tip frequencies\n";
		cerr << "\t-cv <file>\tcross-validation test alignment\n";
		exit(1);
	}

	if(chain == NULL)
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

			if(fexists(name+".param") && !overwrite){
				cerr << "run '" << name << "' exists. use overwrite option -f\n";
				exit(1);
			}else if(overwrite){
				remove((name+".stream").c_str());
				remove((name+".param").c_str());
				remove((name+".trace").c_str());
				remove((name+".treelist").c_str());
				if(nexus)
					remove((name+".treelist.nex").c_str());
			}

			chain = new Biphy(datafile,name,treefile,outgroupfile,modeltype,branchprior,rootprior,correction,dgam,mixture,rootmin,rootmax,every,until,numChains,swapInterval,delta,sigma,saveall,nexus);
		}else{
			if(!fexists(name+".stream")){
				cerr << "run '" << name << "' does not exist\n";
				exit(1);
			}
			if(cvfile != "None" && cvfile.at(0) != '.' && cvfile.at(0) != '/'){
				cvfile = "./"+cvfile;
			}
			if(ppred)
				remove((name+".ppred").c_str());
			if(cvfile != "None")
				remove((name+".cv").c_str());
			if(dolloMapping)
				remove((name+".dollo.fa").c_str());
			chain = new Biphy(name,cvfile,ppred,dolloMapping);
		}

	try
	{
    	chain->run();
	}
	catch (RbException& e)
	{
		cerr << "Error:\t" << e.getMessage() << '\n';
		exit(1);
	}
    
    return 0;
}

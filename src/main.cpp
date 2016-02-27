#include <iostream>
#include <stdlib.h>
#include <fstream>

#include "Exception.h"
#include "Biphy.h"
#include "BinarySubstitutionModel.h"
#include "StringUtilities.h"

using namespace std;

int main (int argc, const char * argv[])
{
	string datafile = "None";
	string treefile = "None";
	string cvfile = "None";
	string outgroupfile = "None";
	string name = "";

	ModelPrior::Type modeltype = ModelPrior::HOMOGENEOUS;
	int mixture = 0;

	BranchPrior::Type branchprior = BranchPrior::DEFAULT;

	int dgam = 4;
	int dbeta = 0;

	RootPrior::Type rootprior = RootPrior::FREE;
	double rootmin = 0.0;
	double rootmax = 1.0;

	int correction = AscertainmentBias::ALL;

	int every = 1000;
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
	bool persite = false;
	bool ancestral = false;
	bool asymmbeta = false;

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
				}
				else if (s == "-f")	{
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
				else if (s == "-s"){
					saveall = true;
				}
				else if (s == "-a"){
                    ancestral = true;
                }
				else if (s == "-dgam"){
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
				}
				else if (s == "-dbeta"){
					i++;
					if (i == argc)	{
						cerr << "error in command: -dbeta <int>\n\n";
						exit(1);
					}
					s = argv[i];
					if (! IsInt(s))	{
						cerr << "error in command: -dbeta <int>\n\n";
						exit(1);
					}
					dbeta = atoi(argv[i]);
				}
				else if (s == "-asymmbeta"){
					asymmbeta = true;
				}
                else if (s == "-rigid"){
					rootprior = RootPrior::RIGID;
				}
				else if (s == "-e")	{
					nexus = true;
				}
				else if (s == "-site")	{
					persite = true;
				}
				else if (s == "-map"){
					dolloMapping = true;
				}
				else if (s == "-nh")	{
					modeltype = ModelPrior::HIERARCHICAL;
				}
				else if (s == "-dollo"){
					modeltype = ModelPrior::DOLLO;
				}
				else if (s == "-h")	{
					modeltype = ModelPrior::HOMOGENEOUS;
				}
				else if (s == "-mk")	{
					modeltype = ModelPrior::MK;
				}
				else if (s == "-ldir")	{
					branchprior = BranchPrior::DIRICHLET;
				}
				else if (s == "-lexp")	{
					branchprior = BranchPrior::EXPONENTIAL;
				}
				else if (s == "-lfixed")	{
					branchprior = BranchPrior::FIXED;
				}
				else if (s == "-lstrict")	{
					branchprior = BranchPrior::STRICT;
				}
				else if (s == "-n")	{
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
				}
				else if (s == "-u")	{
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
					correction = AscertainmentBias::Coding(atoi(argv[i]));
					if(correction < 0 || correction > 15)	{
						cerr << "error in command: -u <int>\n\n";
						exit(1);
					}
				}else if (s == "-mix")	{
					i++;
					if (i == argc)	{
						cerr << "error in command: -mix <int>\n\n";
						exit(1);
					}
					s = argv[i];
					if (! IsInt(s))	{
						cerr << "error in command: -mix <int>\n\n";
						exit(1);
					}
					mixture = atoi(argv[i]);
					if(mixture > 1){
						modeltype = ModelPrior::MIXTURE;
					}else{
						modeltype = ModelPrior::HOMOGENEOUS;
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
					if (!IsInt(s))	{
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
					if (!IsFloat(s))	{
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
					if (!IsFloat(s))	{
						cerr << "error in command: -rp <min> <max>\n\n";
						exit(1);
					}
					rootmax = atof(argv[i]);
					if(rootmin < 0 || rootmax > 1.0 || rootmin > rootmax){
						cerr << "error in command: -rp <min> <max>\n\n";
						exit(1);
					}
					rootprior = RootPrior::TRUNCATED;
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

			if(dbeta > 1 && modeltype != ModelPrior::MK)
			{
				cerr << "error: -dbeta only used with Mk model\n\n";
				exit(1);
			}

			if(rootprior == RootPrior::RIGID && (modeltype < ModelPrior::HIERARCHICAL))
				rootprior = RootPrior::FREE;

			if(rootprior == RootPrior::TRUNCATED && (modeltype == ModelPrior::MIXTURE))
			{
				cerr << "error: truncated root prior incompatible with branch mixture model\n\n";
				exit(1);
			}

			if((branchprior == BranchPrior::FIXED || branchprior == BranchPrior::STRICT) && treefile == "None")
			{
				cerr << "error: ";
				if(branchprior == BranchPrior::FIXED)
					cerr << "fixed branch length";
				else if(branchprior == BranchPrior::STRICT)
					cerr << "strict clock";
				cerr << " prior is used only with fixed input tree file\n\n";
				exit(1);
			}
			else if(branchprior == BranchPrior::DEFAULT)
			{
				if(treefile == "None")
					branchprior = BranchPrior::EXPONENTIAL;
				else
					branchprior = BranchPrior::FIXED;
			}
		}
	}
	catch(Exception &e){
		throw(e);
	}catch(...){

		cerr << "biphy " << VERSION << "\n";
		cerr << '\n';
		cerr << "usage: biphy -d <data file> [-x <every> [<until>] ] <run name>\n";

		cerr << "\nBranch frequency prior:\n";
		cerr << "\t-dollo\t\tirreversible stochastic dollo model\n";
		cerr << "\t-mk\t\tsymmetric reversible Mk model\n";
		cerr << "\t-h\t\tasymmetric reversible model (default)\n";
		cerr << "\t-nh\t\tbranch-heterogeneous asymmetric reversible model\n";
		//cerr << "\t-mix <int>\tbeta mixture with <int> components\n";

		cerr << "\nBranch length prior:\n";
		cerr << "\t-lexp\t\thierarchical exponential prior (default)\n";

		cerr << "\nClock models (require fixed input branch lengths via -t option):\n";
		cerr << "\t-lfixed\t\tfix the branch lengths at their input values (default)\n";
		cerr << "\t-lstrict\tstrict clock. same as -lfixed, but with a tree-length multiplier\n\n";

		cerr << "\nAcross sites mixtures:\n";
		cerr << "\t-dgam <int>\tdiscrete gamma rates mixture with <int> categories (default: 4)\n";
		cerr << "\t-dbeta <int>\tbeta frequency mixture with <int> categories (Mk model only)\n";
		cerr << "\t\t0 or 1 = across sites mixture disabled\n";
		cerr << "\t-asymmbeta\tasymmetric beta mixture\n\n";

		cerr << "\nCorrections for unobservable site patterns:\n";
		cerr << "\t-u <int>\twhere <int> is one of:\n";
		cerr << "\t\t0 = all site patterns are observable (default)\n";
		cerr << "\t\t1 = no constant absence sites\n";
		cerr << "\t\t2 = no constant presence sites\n";
		cerr << "\t\t4 = no singleton presence sites\n";
		cerr << "\t\t8 = no singleton absence sites\n";

		cerr << "\tcombinations are achieved by adding the above values\n";
		cerr << "\te.g.  3 = no constant sites\n";
		cerr << "\t     15 = no uninformative sites\n";

		cerr << "\nOptional constraints:\n";
		cerr << "\t-t <file>\tfixed tree filename\n";
		cerr << "\t-o <file>\toutgroup clade file\n";
		cerr << "\t-rigid\trigid root frequency for heterogeneous models\n";
		cerr << "\t-rp <min> <max>\ttruncate root frequency on (min,max) (asymmetric reversible models only)\n";

		cerr << "\nMCMCMC options:\n";
		cerr << "\t-n <int>\tnumber of chains (default = 1)\n";
		cerr << "\t-si <int>\tchain swap interval (default = 1)\n";
		cerr << "\t-delta <float>\t(default = 0.1)\n";
		cerr << "\t-sigma <float>\t(default = 1)\n";

		cerr << "\nOutput options:\n";
		cerr << "\t-s\t\tsave entire output\n";
		cerr << "\t-e\t\tsave nexus treefile\n";

		cerr << "\nStream-reading options:\n";
		cerr << "\t-ppred\t\tposterior predictive simulation of tip frequencies\n";
		cerr << "\t-cv <file>\tcross-validation test alignment\n";
		cerr << "\t-a\t\tsimulate ancestral states\n";

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

			chain = new Biphy(name,datafile,cvfile,treefile,outgroupfile,modeltype,branchprior,rootprior,correction,dgam,dbeta,asymmbeta,mixture,rootmin,rootmax,every,until,numChains,swapInterval,delta,sigma,saveall,nexus);
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
			if(persite)
			    remove((name+".lnprobs").c_str());
			if(dolloMapping)
				remove((name+".dollo.fa").c_str());
			if(ancestral)
			    remove((name+".mapping").c_str());
			chain = new Biphy(name,cvfile,ppred,dolloMapping,persite,ancestral);
		}

	try
	{
		chain->init();
		chain->run();
	}
	catch (Exception& e)
	{
		cerr << "Error:\t" << e.getMessage() << '\n';
		exit(1);
	}
    
    exit(0);
}

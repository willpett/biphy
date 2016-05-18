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

	BranchPrior::Type branchprior = BranchPrior::EXPONENTIAL;

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
	int swapInterval = 10;
	double delta = 0.05;
	double sigma = 0.3;

	bool saveall = true;
	bool nexus = false;
	int ppred = -1;
	bool dolloMapping = false;
	bool persite = false;
	bool ancestral = false;
	bool asymmbeta = false;

	size_t num_stones = 0;

	Biphy *chain = NULL;

	try	{
		if (argc == 2)	{
		    name = argv[1];
			if (name == "-help" || name == "--help" || name == "-h")	{
				throw 0;
			}
			if(!fexists(name+".stream")){
				throw Exception("run '"+name+"' does not exist");
			}

			chain = new Biphy(name,num_stones);
		}else if (argc == 1){
			throw 0;
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
					saveall = false;
				}
				else if (s == "-a"){
                    ancestral = true;
                }
				else if (s == "-dgam"){
					i++;
					if (i == argc)	{
					    throw Exception("-dgam <int>");
					}
					s = argv[i];
					if (! IsInt(s))	{
					    throw Exception("-dgam <int>");
					}
					dgam = atoi(argv[i]);
				}
				else if (s == "-ss"){
                    i++;
                    if (i == argc)  {
                        throw Exception("-ss <int>");
                    }
                    s = argv[i];
                    if (! IsInt(s)) {
                        throw Exception("-ss <int>");
                    }
                    num_stones = atoi(argv[i]);
                }
				else if (s == "-dbeta"){
					i++;
					if (i == argc)	{
					    throw Exception("-dbeta <int>");
					}
					s = argv[i];
					if (! IsInt(s))	{
					    throw Exception("-dbeta <int>");
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
					    throw Exception("-n <int>");
					}
					s = argv[i];
					if (! IsInt(s))	{
					    throw Exception("-n <int>");
					}
					numChains = atoi(argv[i]);
				}
				else if (s == "-u")	{
					i++;
					if (i == argc)	{
					    throw Exception("-u <int>");
					}
					s = argv[i];
					if (! IsInt(s))	{
					    throw Exception("-u <int>");
					}
					correction = AscertainmentBias::Coding(atoi(argv[i]));
					if(correction < 0 || (correction > 15 && correction != 17))	{
						throw Exception("-u <int>");
					}
				}else if (s == "-mix")	{
					i++;
					if (i == argc)	{
					    throw Exception("-mix <int>");
					}
					s = argv[i];
					if (! IsInt(s))	{
					    throw Exception("-mix <int>");
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
					    throw Exception("-delta <float>");
					}
					s = argv[i];
					if (! IsFloat(s))	{
					    throw Exception("-delta <float>");
					}
					delta = atof(argv[i]);
				}
				else if (s == "-sigma")	{
					i++;
					if (i == argc)	{
					    throw Exception("-sigma <float>");
					}
					s = argv[i];
					if (! IsFloat(s))	{
					    throw Exception("-sigma <float>");
					}
					sigma = atof(argv[i]);
				}
				else if (s == "-si")	{
					i++;
					if (i == argc)	{
					    throw Exception("-si <int>");
					}
					s = argv[i];
					if (!IsInt(s))	{
					    throw Exception("-si <int>");
					}
					swapInterval = atoi(argv[i]);
				}else if (s == "-x")	{
					i++;
					if (i == argc)	{
					    throw Exception("-x <every> <until>");
					}
					s = argv[i];
					if (!IsInt(s))	{
					    throw Exception("-x <every> <until>");
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
					    throw Exception("-rp <min> <max>");
					}
					s = argv[i];
					if (!IsFloat(s))	{
					    throw Exception("-rp <min> <max>");
					}
					rootmin = atof(argv[i]);
					i++;
					if (i == argc)	{
					    throw Exception("-rp <min> <max>");
					}
					s = argv[i];
					if (!IsFloat(s))	{
					    throw Exception("-rp <min> <max>");
					}
					rootmax = atof(argv[i]);
					if(rootmin < 0 || rootmax > 1.0 || rootmin > rootmax){
					    throw Exception("-rp <min> <max>");
					}
					rootprior = RootPrior::TRUNCATED;
				}else if (s == "-ppred")	{
					i++;
					if (i == argc)	{
					    throw Exception("-ppred <int>");
					}
					s = argv[i];
					if (!IsInt(s))	{
					    throw Exception("-ppred <int>");
					}
					ppred = atoi(argv[i]);
				}else{
					if (i != (argc -1))	{
						throw 0;
					}
					name = argv[i];
				}
				i++;
			}

			if(name == ""){
			    throw Exception("could not determine run name");
			}

			if(dbeta > 1 && modeltype != ModelPrior::MK)
			{
			    throw Exception("-dbeta only used with Mk model");
			}

			if(rootprior == RootPrior::RIGID && (modeltype < ModelPrior::HIERARCHICAL))
				rootprior = RootPrior::FREE;

			if(rootprior == RootPrior::TRUNCATED && (modeltype == ModelPrior::MIXTURE))
			{
			    throw Exception("truncated root prior incompatible with branch mixture model");
			}

			if((branchprior == BranchPrior::FIXED || branchprior == BranchPrior::STRICT) && treefile == "None")
			{
			    stringstream ss;
				ss << "";
				if(branchprior == BranchPrior::FIXED)
					ss << "fixed branch length";
				else if(branchprior == BranchPrior::STRICT)
					ss << "strict clock";
				ss << " prior is used only with input tree file";
				throw Exception(ss.str());
			}
		}
	}
	catch(Exception &e)
	{
	    cerr << "Error: " << e.getMessage() << std::endl;
	    return 1;
	}
	catch(...)
	{
	    stringstream ss;
		ss << "biphy " << VERSION << "\n";
		ss << '\n';
		ss << "usage: biphy -d <data file> [-x <every> <until> ] <run name>\n";

		ss << "\nBranch frequency prior:\n";
		ss << "\t-dollo\t\tirreversible stochastic dollo model\n";
		ss << "\t-mk\t\tsymmetric reversible Mk model\n";
		ss << "\t-h\t\tasymmetric reversible model (default)\n";
		ss << "\t-nh\t\tbranch-heterogeneous asymmetric reversible model\n";
		//ss << "\t-mix <int>\tbeta mixture with <int> components\n";

		ss << "\nBranch length prior:\n";
		ss << "\t-lexp\t\thierarchical exponential prior (default)\n";
		ss << "\t-lfixed\t\tfix the branch lengths at input values via the -t option\n";
		ss << "\t-lstrict\tsame as -lfixed with a uniform branch multiplier (clock_rate)\n\n";

		ss << "\nAcross sites mixtures:\n";
		ss << "\t-dgam <int>\tdiscrete gamma rates mixture with <int> categories (default: 4)\n";
		ss << "\t-dbeta <int>\tbeta frequency mixture with <int> categories (Mk model only)\n";
		ss << "\t\t0 or 1 = across sites mixture disabled\n";
		ss << "\t-asymmbeta\tasymmetric beta mixture\n\n";

		ss << "\nCorrections for unobservable site patterns:\n";
		ss << "\t-u <int>\twhere <int> is one of:\n";
		ss << "\t\t0 = all site patterns are observable (default)\n";
		ss << "\t\t1 = no constant absence sites\n";
		ss << "\t\t2 = no constant presence sites\n";
		ss << "\t\t4 = no singleton presence sites\n";
		ss << "\t\t8 = no singleton absence sites\n";

		ss << "\tcombinations are achieved by adding the above values\n";
		ss << "\te.g.  3 = no constant sites\n";
		ss << "\t     15 = no uninformative sites\n";

		ss << "\nOptional constraints:\n";
		ss << "\t-t <file>\tconstraint tree filename\n";
		ss << "\t-o <file>\toutgroup clade file\n";
		ss << "\t-rigid\t\trigid root frequency for heterogeneous models\n";
		ss << "\t-rp <min> <max>\ttruncate root frequency on (min,max) (asymmetric reversible models only)\n";

		ss << "\nMCMCMC options:\n";
		ss << "\t-n <int>\tnumber of chains (default = 1)\n";
		ss << "\t-si <int>\tchain swap interval (default = 10)\n";
		ss << "\t-delta <float>\t(default = 0.05)\n";

		ss << "\nOutput options:\n";
		ss << "\t-s\t\tdo not save entire output\n";
		ss << "\t-e\t\tsave nexus treefile\n";

		ss << "\nStream-reading options:\n";
		ss << "\t-ppred <int>\tposterior predictive simulation (0 or 1)\n";
		ss << "\t-cv <file>\tcross-validation test alignment\n";
		ss << "\t-a\t\tsimulate ancestral states\n";
		ss << "\t-ss <int>\tinitialize steppingstone sampler with <int> chains\n";

		cerr << ss.str();

		return 1;
	}

	try
	{
        if(chain == NULL)
            if(datafile == "None" && num_stones > 0)
            {
                if(!fexists(name+".param"))
                {
                    throw Exception("run '"+name+"' does not exist");
                }

                chain = new Biphy(name,num_stones);
            }
            else if(datafile != "None")
            {
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
                    throw Exception("run '"+name+"' exists. use overwrite option -f");
                }else if(overwrite){
                    remove((name+".stream").c_str());
                    remove((name+".param").c_str());
                    remove((name+".trace").c_str());
                    remove((name+".treelist").c_str());
                    if(nexus)
                        remove((name+".treelist.nex").c_str());
                }

                chain = new Biphy(name,datafile,cvfile,treefile,outgroupfile,modeltype,branchprior,rootprior,correction,dgam,dbeta,asymmbeta,mixture,rootmin,rootmax,every,until,numChains,swapInterval,delta,sigma,saveall,nexus);
            }
            else
            {
                if(!fexists(name+".stream")){
                    throw Exception("run '"+name+"' does not exist");
                }
                if(cvfile != "None" && cvfile.at(0) != '.' && cvfile.at(0) != '/'){
                    cvfile = "./"+cvfile;
                }
                if(ppred)
                {
                    remove((name+".ppred").c_str());
                }
                if(cvfile != "None")
                    remove((name+".cv").c_str());
                if(persite)
                    remove((name+".lnprobs").c_str());
                if(dolloMapping)
                    remove((name+".dollo.fa").c_str());
                if(ancestral)
                    remove((name+".mapping").c_str());
                chain = new Biphy(name,cvfile,ppred,dolloMapping,persite,ancestral,every,until);
            }

        chain->init();
        chain->run();
	}
	catch(Exception& e)
	{
	    cerr << "Error: " << e.getMessage() << std::endl;
	}
}

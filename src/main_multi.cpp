#include "RbException.h"
#include <iostream>
#include <stdlib.h>
#include <fstream>

#include "BinaryCharEvoModel.h"
#include "MultiBiphyRhoPi.h"
#include "MultiBiphyLambdaMu.h"
#include "util.h"

using namespace std;

int main (int argc, const char * argv[])
{
	string datafile = "None";
	string treefile = "None";
	string name = "";

	enum Parameterization { LAMBDA_MU, RHO_PI };

	MultiBiphy::ModelType modeltype = MultiBiphy::REVERSIBLE;
	MultiBiphy::RatePrior A_prior 	= MultiBiphy::HOMOGENEOUS;
	MultiBiphy::RatePrior B_prior 	= MultiBiphy::HOMOGENEOUS;
	MultiBiphy::RatePrior A_species = MultiBiphy::HOMOGENEOUS;
	MultiBiphy::RatePrior B_species = MultiBiphy::HOMOGENEOUS;

	Parameterization parameterization = LAMBDA_MU;
	int AB = 0;

	int dgam = 4;

	int correction = RevBayesCore::NONE;

	int every = 1;
	int until = -1;
	bool overwrite = false;

	int numChains = 1;
	int swapInterval = 1;
	double delta = 0.05;
	double sigma = 1;

	bool saveall = false;
	bool ancestral = false;

	MultiBiphy *chain = NULL;

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

			chain = new MultiBiphyRhoPi(name);
			chain = new MultiBiphyLambdaMu(name);
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
				}else if (s == "-a")	{
					ancestral = true;
				}
				else if ((s == "-t"))	{
					i++;
					treefile = argv[i];
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
				}else if (s == "-rev")	{
					modeltype = MultiBiphyRhoPi::REVERSIBLE;
				//}else if (s == "-dollo"){
				//	modeltype = MultiBiphyRhoPi::DOLLO;
				}else if (s == "-rho_pi")	{
					parameterization = RHO_PI;
					i++;
					if (i == argc)	{
						cerr << "error in command: -rho_pi <int>\n\n";
						exit(1);
					}
					s = argv[i];
					if (! IsInt(s))	{
						cerr << "error in command: -rho_pi <int>\n\n";
						exit(1);
					}
					AB = atoi(argv[i]);
				}else if (s == "-lambda_mu")	{
					parameterization = LAMBDA_MU;
					i++;
					if (i == argc)	{
						cerr << "error in command: -lambda_mu <int>\n\n";
						exit(1);
					}
					s = argv[i];
					if (! IsInt(s))	{
						cerr << "error in command: -lambda_mu <int>\n\n";
						exit(1);
					}
					AB = atoi(argv[i]);
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

			if(modeltype == MultiBiphyRhoPi::DOLLO)
				correction = (RevBayesCore::CorrectionType)(correction | RevBayesCore::NO_ABSENT_SITES);

			if(AB & 0x01)
				A_prior = MultiBiphy::HIERARCHICAL;
			if(AB & 0x02)
				B_prior = MultiBiphy::HIERARCHICAL;
			if(AB & 0x04)
				A_species = MultiBiphy::HIERARCHICAL;
			if(AB & 0x08)
				B_species = MultiBiphy::HIERARCHICAL;
		}
	}
	catch(RbException &e){
		throw(e);
	}catch(...){

		cerr << "multibiphy " << VERSION << "\n";
		cerr << '\n';
		cerr << "usage: multibiphy -t <tree file> -d <data file> [-x <every> [<until>] ] <run name>\n";

		cerr << "\nModel type:\n";
		//cerr << "\t-dollo\t\tstochastic dollo model (enables -u 1)\n";
		cerr << "\t-rev\t\treversible model (default)\n";

		cerr << "\nPrior parameterization:\n";
		cerr << "\t-rho_pi <int>\n";
		cerr << "\t-lambda_mu <int>\n";
		cerr << "\t\twhere -A_B <int> is one of:\n";
		cerr << "\t\t0 = homogeneous A and B (default)\n";
		cerr << "\t\t1 = heterogeneous A across trees\n";
		cerr << "\t\t2 = heterogeneous B across trees\n";
		cerr << "\t\t4 = heterogeneous A across species\n";
		cerr << "\t\t8 = heterogeneous B across species\n";

		cerr << "\tcombinations are achieved by adding the above values\n";
		cerr << "\te.g.  3 = heterogeneous A,B across trees\n";
		cerr << "\t     12 = heterogeneous A,B across species\n";

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

		cerr << "\t-a\tsimulate ancestral states";

		cerr << "\nMCMCMC options:\n";
		cerr << "\t-n <int>\tnumber of chains (default = 1)\n";
		cerr << "\t-si <int>\tchain swap interval (default = 1)\n";
		cerr << "\t-delta <float>\t(default = 0.1)\n";
		cerr << "\t-sigma <float>\t(default = 1)\n";

		cerr << "\nOutput options:\n";
		cerr << "\t-s\t\tsave entire output (default: disabled)\n";
		exit(1);
	}

	if(chain == NULL)
		if(datafile != "None" && treefile != "None"){
			if(datafile.at(0) != '.' && datafile.at(0) != '/'){
				datafile = "./"+datafile;
			}
			if(treefile.at(0) != '.' && treefile.at(0) != '/'){
				treefile = "./"+treefile;
			}

			if(fexists(name+".param") && !overwrite){
				cerr << "run '" << name << "' exists. use overwrite option -f\n";
				exit(1);
			}else if(overwrite){
				remove((name+".stream").c_str());
				remove((name+".param").c_str());
				remove((name+".trace").c_str());
			}


			if(parameterization == LAMBDA_MU){
				chain = new MultiBiphyLambdaMu(datafile,name,treefile,modeltype,A_prior,B_prior,A_species,B_species,correction,dgam,every,until,numChains,swapInterval,delta,sigma,saveall,ancestral);
				std::cerr << "done" << std::endl;
			}else{
				chain = new MultiBiphyRhoPi(datafile,name,treefile,modeltype,A_prior,B_prior,A_species,B_species,correction,dgam,every,until,numChains,swapInterval,delta,sigma,saveall,ancestral);
			}
		}else{
			if(!fexists(name+".stream")){
				cerr << "run '" << name << "' does not exist\n";
				exit(1);
			}
			if(parameterization == LAMBDA_MU)
				chain = new MultiBiphyLambdaMu(name,ancestral);
			else
				chain = new MultiBiphyRhoPi(name,ancestral);
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

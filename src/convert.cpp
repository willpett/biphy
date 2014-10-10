#include <iostream>
#include <stdlib.h>
#include <fstream>

#include "RbException.h"
#include "BranchLengthTree.h"
#include "TimeTree.h"
#include "NclReader.h"
#include "Topology.h"
#include "RbUtil.h"

using namespace RevBayesCore;

int main (int argc, const char * argv[])
{
		string datafile = "None";
		string treefile = "None";

		try	{

			if (argc <= 2)	{
				throw(0);
			}

			int i = 1;
			while (i < argc)	{
				string s = argv[i];

				if (s == "-d")	{
					i++;
					datafile = argv[i];
				}
				else if ((s == "-t") || (s == "-T"))	{
					i++;
					treefile = argv[i];
				}
				i++;
			}

			if(datafile == "None" || treefile == "None")
				throw(0);
		}
		catch(...)	{
			cerr << "usage: convert -d <tracefile> -t <treefile>\n";
			exit(1);
		}

	ifstream treelist( treefile.c_str());
	ifstream trace( datafile.c_str());

	//std::vector<BranchLengthTree*>* trees = NclReader::getInstance().readBranchLengthTrees(treefile,"nexus");

	size_t ncol = 10;
	vector<string> values(ncol);
	for(size_t i = 0; i < ncol ; i++)
		trace >> values[i];

	BranchLengthTree * tree = new BranchLengthTree();
	try{
		while(treelist >> *tree){
		//for(size_t t = 0; t < trees->size() ; t++){
			//tree = (*trees)[t];
			size_t numBranches = 2*tree->getNumberOfTips() - 2;

			for(size_t i = 0; i < ncol ; i++)
				trace >> values[i];

			cout << values[8] << endl;
			cout << values[9] << endl;
			for(size_t i = 0; i < numBranches ; i++)
				cout << tree->getBranchLength(i) << endl;

			cout << values[6] << endl;
			cout << values[5] << endl;
			cout << values[7] << endl;

			for(size_t i = 0; i < numBranches ; i++)
				cout << tree->getNode(i).getBranchParameter("pi_vector") << endl;

			tree->clearBranchParameters();
			cout << tree->getTopology() << endl;
		}
	}catch(RbException &e){
		exit(0);
		//cerr << "error:\t" << e.getMessage() << endl;
	}
    
    return 0;
}

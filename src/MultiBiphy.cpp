#include "MultiBiphy.h"
#include "NHXTreeReader.h"
#include "TopologyNode.h"
#include "DeterministicNode.h"
#include "StandardState.h"
#include "BranchLengthTree.h"
#include "LnCorrectionFunction.h"
#include "NegativeBinomialDistribution.h"
#include "GibbsMove.h"
#include "VectorFunction.h"
#include "SumFunction.h"
#include "MappingMonitor.h"

#include <algorithm>

using namespace RevBayesCore;

MultiBiphy::MultiBiphy(const std::string n,
						const std::string df,
						const std::string t,
						ModelType mt,
						RatePrior Ap,
						RatePrior Bp,
						RatePrior As,
						RatePrior Bs,
						bool anc,
						int c,
						int d,
						int e,
						int u,
						int num,
						int swa,
						double de,
						double si,
						bool sav) : Biphy(n,df,t,c,d,e,u,num,swa,de,si,sav),
	modeltype(mt), A_prior(Ap), B_prior(Bp), A_species(As), B_species(Bs), ancestral( anc)
{
	save();
}

MultiBiphy::MultiBiphy(const std::string n, bool anc) : Biphy(n),
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

MultiBiphy::MultiBiphy(const std::string n) : Biphy(n)
{
	open();
}

void MultiBiphy::printConfiguration( void ) {

	std::cout << std::endl;

    if(modeltype == DOLLO){
		std::cout << "dollo model";
	}else{
		std::cout << "reversible model\n";
	}

    Biphy::printConfiguration();
}

void MultiBiphy::readInputFiles( void ) {

	Biphy::readInputFiles();

	if(trees.size() > 0)
		for(size_t i = 0; i < trees.size(); i++)
			delete trees[i];

	NHXTreeReader reader;
	trees = *(reader.readBranchLengthTrees(treeFile));

	node2speciesMaps.clear();
	rootSpecies.clear();
	speciesIndex.clear();
	std::map<size_t, size_t> counts;

	size_t found = 0;
	std::vector<BranchLengthTree*> newtrees = trees;
	for(size_t i = 0; i < trees.size(); i++){
		std::vector<std::string> tree_names = trees[i]->getTipNames();
		std::sort(tree_names.begin(), tree_names.end());
		for(size_t j = 0; j < data.size(); j++){
			std::vector<std::string> data_names = data[j]->getTaxonNames();
			std::sort(data_names.begin(), data_names.end());
			if(std::equal(data_names.begin(), data_names.end(), tree_names.begin())){
				found++;
				newtrees[j] = trees[i];
				break;
			}
		}
		std::vector<TopologyNode*> nodes = trees[i]->getNodes();

		std::map<size_t, size_t> smap;
		for(size_t index = 0; index < nodes.size(); index++){
			size_t sIndex = size_t(nodes[index]->getBranchParameter("S"));

			smap[index] = sIndex;
			counts[sIndex]++;
			if(nodes[index]->isRoot()){
				rootSpecies.push_back(sIndex);
			}
		}

		node2speciesMaps.push_back(smap);
	}

	if(found != trees.size()){
		std::cerr << found << " Error: taxon labels in matrices do not match those in trees\n";
		exit(1);
	}
	trees = newtrees;

	numSpecies = 0;
	for(std::map<size_t, size_t>::iterator it = counts.begin(); it != counts.end(); it++)
		speciesIndex[it->first] = numSpecies++;

	numTrees = trees.size();
	wt = numTrees > 0 ? (int) log10 ((double) numTrees) + 1 : 1;
	ws = numSpecies > 0 ? (int) log10 ((double) numSpecies) + 1 : 1;
}

void MultiBiphy::initModel( void ) {
	size_t wt = trees.size() > 0 ? (int) log10 ((double) trees.size()) + 1 : 1;

	std::cout << "Initializing " << trees.size() << " likelihood functions\n";
	for(size_t i = 0; i < trees.size(); i++){
		try{
			BinaryCharEvoModel<BranchLengthTree> *charModel = initSubModel(i);

			tmp_name.str("");
			tmp_name << "S(" << std::setfill('0') << std::setw(wt) << i << ")";
			StochasticNode< AbstractCharacterData >* charactermodel = new StochasticNode< AbstractCharacterData >(tmp_name.str(), charModel );
			charactermodel->clamp( data[i] );

			dataNodes.push_back(charactermodel);

		//std::cerr << charactermodel->getLnProbability() << std::endl;
		}catch(RbException& e){
			std::cerr << "Error building model for " << data[i]->getFileName() << std::endl;
			std::cerr << e.getMessage() << std::endl;
		}
	}
}

void MultiBiphy::initMCMC( void ) {
	std::vector<const TypedDagNode<int>* > M_vector;

	if(correction != NONE && modeltype == REVERSIBLE){
		for(size_t i = 0; i < dataNodes.size(); i++){
			tmp_name.str("");
			tmp_name << "N(" << std::setfill('0') << std::setw(wt) << i << ")";
			ConstantNode<int>* N = new ConstantNode<int>(tmp_name.str(), new int(data[i]->getNumberOfIncludedCharacters() + 1));

			tmp_name.str("");
			tmp_name << "p(" << std::setfill('0') << std::setw(wt) << i << ")";
			DeterministicNode<double>* p = new DeterministicNode<double>(tmp_name.str(), new LnCorrectionFunction<StandardState, BranchLengthTree>(dataNodes.back()) );

			tmp_name.str("");
			tmp_name << "M(" << std::setfill('0') << std::setw(wt) << i << ")";
			M_vector.push_back(new StochasticNode<int>(tmp_name.str(), new NegativeBinomialDistribution(N, p) ));

			moves.push_back( new GibbsMove<int>((StochasticNode<int>*)M_vector.back(), 1.0 ) );
			//monitoredNodes.push_back((StochasticNode<int>*)M_vector[i] );
		}

		DeterministicNode<std::vector<int> >* M_vector_node = new DeterministicNode< std::vector<int> >( "M_vector", new VectorFunction<int>( M_vector ) );
		DeterministicNode<int>* sum = new DeterministicNode<int>("M", new SumFunction<int>(M_vector_node));

		monitoredNodes.push_back( sum );
	}

	bool useParallelMcmcmc = (numChains > 1);

	if(ancestral){
		for(size_t i = 0; i < dataNodes.size(); i++){
			std::string treelist = name+"."+data[i]->getFileName()+".treelist";
			if(!restart)
				remove(treelist.c_str());
			monitors.push_back( new MappingMonitor<StandardState,BranchLengthTree>( dataNodes[i], every, name+"."+data[i]->getFileName()+".treelist", useParallelMcmcmc || restart) );
		}
	}

	Biphy::initMCMC();
}

void MultiBiphy::openParams( std::ifstream& is ) {
	int i;
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
	is >> ancestral;
}

void MultiBiphy::saveParams( std::ofstream& os ) {

	os << modeltype << "\n";
	os << A_prior << "\n";
	os << A_species << "\n";
	os << B_prior << "\n";
	os << B_species << "\n";
	os << ancestral << "\n";
}

#include "ParallelMcmcmc.h"
#include "RandomNumberFactory.h"
#include "RandomNumberGenerator.h"
#include <iostream>
#include <vector>
#include <cmath>
#include <sstream>
#include <limits>
#include <stdlib.h>
#include "../../utils/Exception.h"

#if defined (USE_LIB_OPENMP)
    #include <omp.h>
#endif
 
#define DEBUG_PMC3 0

ParallelMcmcmc::ParallelMcmcmc(Model* m, const std::vector<Move*> &moves, const std::vector<Monitor*> &mons, std::string fn, std::string sT, int ev, int nc, int np, int si, double dt, double st, double sh, bool saveall) : Cloneable( ),
			filename(fn),
			numChains(nc),
			numProcesses(np),
			scheduleType(sT),
			every(ev),
			currentGeneration(0),
			swapInterval(si),
			delta(dt),
			sigma(st),
			startingHeat(sh),
			saveall(saveall)
{
    activeIndex = 0;
    
    for (size_t i = 0; i < numChains; i++)
    {
        // get chain heat
        double b = computeBeta(delta,sigma,i) * startingHeat;
        
        // create chains
        bool a = (i == 0 ? true : false);
        Mcmc* oneChain = new Mcmc(m, moves, mons, scheduleType, a, b, i);
    	oneChain->setChainIndex(i);
        oneChain->startMonitors();
        
        // add chain to team
        chains.push_back(oneChain);
        chainIdxByHeat.push_back(i);
    }

    numNodes = chains[0]->getNumNodes();
    //std::cout << "\n";
    
    
    // assign chains to processors
#if defined (USE_LIB_OPENMP)
    omp_set_num_threads((unsigned)numProcesses);
#endif
    
    if (numChains < numProcesses)
        numProcesses = numChains;
    
    chainsPerProcess.resize(numProcesses);
    for (size_t i = 0, j = 0; i < numChains; i++, j++)
    {
        if (j >= numProcesses)
            j = 0;
        chainsPerProcess[j].push_back(i);
    }

    for (size_t i = 0; i < numChains; i++)
    {
        
        ;//std::cout << i << ": " << chains[i]->getChainHeat() << " ";
        ;//std::cout << chains[i]->getLnPosterior() << "\n";
    }
}

ParallelMcmcmc::ParallelMcmcmc(const ParallelMcmcmc &m)
{
    sigma = m.sigma;
    delta = m.delta;
    startingHeat = m.startingHeat;
    numChains = m.numChains;
    numProcesses = m.numProcesses;
    numNodes = m.numNodes;
    filename = m.filename;
    every = m.every;
    saveall = m.saveall;

    chainIdxByHeat = m.chainIdxByHeat;
    swapInterval = m.swapInterval;
    activeIndex = m.activeIndex;
    scheduleType = m.scheduleType;

    chainsPerProcess.clear();
    for (size_t i = 0; i < m.chainsPerProcess.size(); i++)
    {
        //chainsPerProcess[i].push_back(m.chainsPerProcess[i]);
        chainsPerProcess.push_back(m.chainsPerProcess[i]);
    }
    
    chains.clear();
    for (size_t i = 0; i < m.chains.size(); i++)
        chains.push_back(new Mcmc( *(m.chains[i])) );
    
    currentGeneration = m.currentGeneration;
}

ParallelMcmcmc::~ParallelMcmcmc(void)
{
    for (size_t i = 0; i < chains.size(); i++)
    {
        
        delete chains[i];
    }
    chains.clear();
}

void ParallelMcmcmc::initialize(void)
{
    
}

double ParallelMcmcmc::computeBeta(double d, double s, size_t idx)
{
    // MJL: May want other distributions of beta in the future
    return pow(1 + d, -pow(idx,s));
}

void ParallelMcmcmc::burnin(int g, int ti)
{
    
}

ParallelMcmcmc* ParallelMcmcmc::clone(void) const
{
    return new ParallelMcmcmc(*this);
}

void ParallelMcmcmc::printOperatorSummary(void) const
{
    for (size_t i = 0; i < numChains; i++)
    {
        //std::cout << "\nChain " << i << ":\n";
        chains[i]->printOperatorSummary();
    }
}

unsigned int ParallelMcmcmc::getCurrentGeneration(void)
{
    return currentGeneration;
}

void ParallelMcmcmc::run(int generations)
{
	if(stream.is_open())
		stream.close();

	stream.open( filename.c_str(), std::fstream::out | std::fstream::app);

    // print file header
    if (currentGeneration == 0){
        chains[0]->monitor(0);
		stream << currentGeneration << "\n";
		toStream(stream);
    }

    // run chain
    for (int i = currentGeneration+1; i <= generations || generations == -1; i += swapInterval)
    {
        // start parallel job per block of swapInterval cycles
        size_t np = numProcesses; // in fact, used by the macro below
        int pid = 0;
        
        #pragma omp parallel default(shared) private(np, pid)
        {
            #ifdef USE_LIB_OPENMP
                pid = omp_get_thread_num();
            #endif
            
            // Create process per chain
            for (size_t j = 0; j < chainsPerProcess[pid].size(); j++)
            {
                // get chain index from job vector
                size_t chainIdx = chainsPerProcess[pid][j];
                
                // Advance cycles in blocks of size swapInterval
                for (size_t k = 0; k < swapInterval && (i+k) <= generations; k++)
                {
                    // advance chain j by a single cycle
                    chains[chainIdx]->nextCycle(true);
/*
                    // monitor chain activeIndex only
                    //if (chainIdx == activeIndex)
                    if (chains[chainIdx]->isChainActive() )
                    {

            			//std::cerr << "monitoring " << chainIdx << std::endl;
                        //chains[activeIndex]->
                    	//std::cout << "monitoring " << chainIdx << std::endl;
                        chains[chainIdx]->monitor(i+k);
                    }
                    //std::cout << chainIdx << "    lnPosterior  " << chains[chainIdx]->getLnPosterior() << " chainHeat  " << chains[chainIdx]->getChainHeat() << "\n";
                    //chains[chainIdx]->monitor(i+k);
                     *
                     */
                }
            }
            
            // synchronize chains
            #pragma omp barrier

        } // processor job end

        currentGeneration += swapInterval;

        swapChains();

		chains[chainIdxByHeat[0]]->monitor(currentGeneration);

		if(currentGeneration % every*swapInterval == 0){
			//std::cerr << "-----\n";
			std::stringstream output;
			output << currentGeneration << "\n";
			toStream(output);

			if(!saveall){
				stream.close();
				stream.open( filename.c_str(), std::fstream::trunc | std::fstream::out);
			}
			stream << output.str();
		}
    }
}

// MJL: allow swapChains to take a swap function -- e.g. pairwise swap for 1..n-1
void ParallelMcmcmc::swapChains(void)
{
    size_t numChains = chains.size();
    double lnProposalRatio = 0.0;
    
    // exit if there is only one chain
    if (numChains < 2)
        return;
    
    size_t numAccepted = 0;
    
    //for (size_t i = 1; i < numChains; i++)
    for (size_t i = numChains-1; i > 0; i--)
    {
        
        // swap adjacent chains
        size_t j = chainIdxByHeat[i-1];
        size_t k = chainIdxByHeat[i];
        
        // compute exchange ratio
        double bj = chains[j]->getChainHeat();
        double bk = chains[k]->getChainHeat();
        double lnPj = chains[j]->getLnPosterior();
        double lnPk = chains[k]->getLnPosterior();
        double lnR = bj * (lnPk - lnPj) + bk * (lnPj - lnPk) + lnProposalRatio;
        
        // determine whether we accept or reject the chain swap
        bool accept = false;
        double u = GLOBAL_RNG->uniform01();
        if (lnR >= 0)
            accept = true;
        else if (lnR < -100)
            accept = false;
        else if (u < exp(lnR))
            accept = true;
        else
            accept = false;
        
        if (accept == true)
            numAccepted++;
        
        // test override
        //accept = true;
#if DEBUG_PMC3
        std::cout << "\nbj " << bj << "; bk " << bk << "; lnPj " << lnPj << "; lnPk " << lnPk << "\n";
        std::cout << "bj*(lnPk-lnPj) " << bj*(lnPk-lnPj) << "; bk*(lnPj-lnPk) " << bk*(lnPj-lnPk) << "\n";
        std::cout << "swapChains()\t" << j << " <--> " << k << "  " << lnR << "\n";
        std::cout << u << "  " << exp(lnR) << "  " << (accept ? "accept\n" : "reject\n");
#endif
        
        // on accept, swap beta values and active chains
        if (accept)
        {
            // swap chain heat values
            chains[j]->setChainHeat(bk);
            chains[k]->setChainHeat(bj);
            
            //size_t tmpIdx = j;
            chainIdxByHeat[i-1] = k;
            chainIdxByHeat[i] = j;

            // swap active chain
            if (activeIndex == j)
            {
                activeIndex = (int)k;
                chains[j]->setChainActive(false);
                chains[k]->setChainActive(true);
            }
        }
        //std::cout << "activeIndex " << activeIndex << "\n";
    }
    
#if DEBUG_PMC3
   
    int nc = (numChains < 10 || true ? numChains : 10);
    for (int j = 0; j < nc; j++)
    {
        int i = chainIdxByHeat[j];
        std::cout << i << " " << chains[i]->getChainHeat() << " * " << chains[i]->getLnPosterior() << " = " << chains[i]->getChainHeat() * chains[i]->getLnPosterior() << "\n";
        //chains[i]->getModelLnProbability() << " " << (chains[i]->isChainActive() ? "*" : "") << (i == activeIndex ? "#" : "") << "\n";
    }
    std::cout << "freq accepted: " << (double)numAccepted/(numChains-1) << "\n";
    
    std::cout << "\n";
# endif
    
}

bool ParallelMcmcmc::restore(void){
	if(stream.is_open())
		stream.close();

	stream.open(filename.c_str(), std::fstream::in);

	size_t pos = stream.tellg();
	size_t lastsample = pos;

	while(!stream.eof()){
		stream >> currentGeneration;
		fromStream(stream,false);

		lastsample = pos;
		pos = stream.tellg();
	}

	stream.clear();
	stream.seekg(lastsample);

	stream >> currentGeneration;
	fromStream(stream);

	//exit(1);

	return true;
}

void ParallelMcmcmc::fromStream(std::istream& is, bool keep, bool keepCold){
	int index = 0;
	for ( size_t i = 0; i < numChains; i++){
		if(numChains > 1){
			is >> index;
			if(!is)
				throw(Exception("premature end of stream"));
		}

		chains[index]->fromStream(is, keep || (keepCold && i == 0));

		chains[index]->setChainActive(i == 0);
		chains[index]->setChainHeat(computeBeta(delta,sigma,i));

		chainIdxByHeat[i] = index;
		if(i == 0)
			activeIndex = index;

		//if(keep)
		//	std::cerr << index << "\t" << chains[index]->getChainHeat() << " * " << chains[index]->getModelLnProbability() << "\t= " << chains[index]->getChainHeat()*chains[index]->getModelLnProbability() << std::endl;
	}
}

void ParallelMcmcmc::toStream(std::ostream& os){
	std::stringstream sample;
	sample.precision(std::numeric_limits<double>::digits10+2);
	for ( size_t i = 0; i < numChains; i++){
		if(numChains > 1)
			sample << chainIdxByHeat[i] << "\n";
		chains[chainIdxByHeat[i]]->toStream(sample);
	}
	os << sample.str();
	os.flush();
}

void ParallelMcmcmc::readStream(size_t burnin)
{
	if(stream.is_open())
		stream.close();

	stream.open(filename.c_str(), std::fstream::in);

	// print file header

    for (int k=0; k==k; k++)
    {
		stream >> currentGeneration;
		if(!stream)
			throw(Exception("premature end of stream"));

		fromStream(stream, false, true);
		//std::cerr << ".";
		if(k == burnin)
			chains[chainIdxByHeat[0]]->monitor(0);
		else if(k > burnin)
			chains[chainIdxByHeat[0]]->monitor(currentGeneration);

		if(stream.eof())
			break;
    }
}

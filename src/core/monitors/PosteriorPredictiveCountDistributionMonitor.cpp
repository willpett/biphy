/**
 * @file
 * This file contains the implementation of PosteriorPredictiveCountDistributionMonitor, used to save
 * information to file about the tree and additional variables at nodes, e.g. branch rates, 
 * in extended Newick format.
 *
 * @brief Implementation of PosteriorPredictiveCountDistributionMonitor
 *
 * (c) Copyright 2009- under GPL version 3
 * @date Last modified: $Date: 2012-11-26 05:34:51 -0800 (Mon, 26 Nov 2012) $
 * @author The RevBayes Development Core Team
 * @license GPL version 3
 * @version 1.0
 * @since 2012-06-21, version 1.0
 *
 * $Id: FileMonitor.cpp 1867 2012-11-26 13:34:51Z hoehna $
 */


#include "PosteriorPredictiveCountDistributionMonitor.h"
#include "DagNode.h"
#include "Model.h"
#include "Monitor.h"
#include "Exception.h"
#include "BinarySubstitutionModel.h"

#include <algorithm>

/* Constructor */
PosteriorPredictiveCountDistributionMonitor::PosteriorPredictiveCountDistributionMonitor(StochasticNode<BinaryCharacterData> *t, int g, const std::string &fname, bool ap) :
Monitor(g,t), outStream(), data( t ), filename( fname ), append(ap)
{

}


PosteriorPredictiveCountDistributionMonitor::PosteriorPredictiveCountDistributionMonitor(const PosteriorPredictiveCountDistributionMonitor &m) :
		Monitor( m ), outStream(), data( m.data )
{
    
    filename    = m.filename;
    append      = m.append;
    
    //if (m.outStream.is_open())
    //	openStream();
}


/* Clone the object */
PosteriorPredictiveCountDistributionMonitor* PosteriorPredictiveCountDistributionMonitor::clone(void) const {
    
    return new PosteriorPredictiveCountDistributionMonitor(*this);
}


void PosteriorPredictiveCountDistributionMonitor::closeStream() {
    outStream.close();
}


/** Monitor value at generation gen */
void PosteriorPredictiveCountDistributionMonitor::monitor(long gen) {
    
    size_t numTaxa = data->getValue().getNumberOfTaxa();

	// print observed values
	if(gen == 0)
	{
		std::vector<int> frequencies = data->getValue().computeCountDistribution();
		for(size_t i=0;i<=numTaxa;i++){
			std::cout << frequencies[i] << "\t";
		}
		std::cout << std::endl;
	}

    // get the printing frequency
    int samplingFrequency = printgen;

    if (gen % samplingFrequency == 0) {
        //BinaryCharacterData* oldData = data->getValue().clone();
    	data->redraw();
    	BinarySubstitutionModel* model = dynamic_cast<BinarySubstitutionModel*>(&(data->getDistribution()));
    	std::vector<int> frequencies = model->getCountDistribution();
    	//std::vector<std::string> newNames = newData->getTaxonNames();
    	for(size_t i=0;i<=numTaxa;i++){
			//size_t idx = std::find(newNames.begin(),newNames.end(),taxonNames[i]) - newNames.begin();
			outStream << frequencies[i] << "\t";
    	}
    	outStream << std::endl;
    	//data->clamp(oldData);
    }
}


/** open the file stream for printing */
void PosteriorPredictiveCountDistributionMonitor::openStream(void) {
    
    // open the stream to the file
    if (append)
        outStream.open( filename.c_str(), std::fstream::out | std::fstream::app);
    else
        outStream.open( filename.c_str(), std::fstream::out);
}

/** Print header for monitored values */
void PosteriorPredictiveCountDistributionMonitor::printHeader() {

	size_t numTaxa = data->getValue().getNumberOfTaxa();

	for(size_t i=0;i<=numTaxa;i++){
		outStream << i;
		if(i != numTaxa)
			outStream << "\t";
	}
	outStream << std::endl;
}


void PosteriorPredictiveCountDistributionMonitor::swapNode(DagNode *oldN, DagNode *newN) {
	if ( oldN == data ) {
		data = static_cast< StochasticNode<BinaryCharacterData> *>( newN );
	}
    // delegate to base class
    Monitor::swapNode(oldN, newN);
}



/**
 * @file
 * This file contains the implementation of PosteriorPredictiveStateFrequencyMonitor, used to save
 * information to file about the tree and additional variables at nodes, e.g. branch rates, 
 * in extended Newick format.
 *
 * @brief Implementation of PosteriorPredictiveStateFrequencyMonitor
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


#include "PosteriorPredictiveStateFrequencyMonitor.h"
#include "DagNode.h"
#include "Model.h"
#include "Monitor.h"
#include "Exception.h"

#include <algorithm>

/* Constructor */
PosteriorPredictiveStateFrequencyMonitor::PosteriorPredictiveStateFrequencyMonitor(StochasticNode<BinaryCharacterData> *t, int g, const std::string &fname, bool ap) :
Monitor(g,t), outStream(), data( t ), filename( fname ), append(ap)
{
    for(size_t i=0;i<data->getValue().getNumberOfTaxa();i++){
    	taxonNames.push_back(data->getValue().getTaxonNameWithIndex(i));
    }
}


PosteriorPredictiveStateFrequencyMonitor::PosteriorPredictiveStateFrequencyMonitor(const PosteriorPredictiveStateFrequencyMonitor &m) :
		Monitor( m ), outStream(), data( m.data ), taxonNames( m.taxonNames)
{
    
    filename    = m.filename;
    append      = m.append;
    
    //if (m.outStream.is_open())
    //  openStream();
}


/* Clone the object */
PosteriorPredictiveStateFrequencyMonitor* PosteriorPredictiveStateFrequencyMonitor::clone(void) const {
    
    return new PosteriorPredictiveStateFrequencyMonitor(*this);
}


void PosteriorPredictiveStateFrequencyMonitor::closeStream() {
    outStream.close();
}


/** Monitor value at generation gen */
void PosteriorPredictiveStateFrequencyMonitor::monitor(long gen) {
    
    // get the printing frequency
    int samplingFrequency = printgen;
    
    if (gen % samplingFrequency == 0) {
        BinaryCharacterData* oldData = data->getValue().clone();
    	data->redraw();
    	BinaryCharacterData* newData = &data->getValue();
    	std::vector<double> frequencies = newData->computeStateFrequencies();
    	//std::vector<std::string> newNames = newData->getTaxonNames();
    	for(size_t i=0;i<taxonNames.size();i++){
			//size_t idx = std::find(newNames.begin(),newNames.end(),taxonNames[i]) - newNames.begin();
			outStream << frequencies[i] << "\t";
    	}
    	outStream << std::endl;
    	data->clamp(oldData);
    }
}


/** open the file stream for printing */
void PosteriorPredictiveStateFrequencyMonitor::openStream(void) {
    
    // open the stream to the file
    if (append)
        outStream.open( filename.c_str(), std::fstream::out | std::fstream::app);
    else
        outStream.open( filename.c_str(), std::fstream::out);
}

/** Print header for monitored values */
void PosteriorPredictiveStateFrequencyMonitor::printHeader() {
	for(size_t i=0;i<taxonNames.size();i++){
		outStream << taxonNames[i] << "\t";
	}
	outStream << std::endl;
}


void PosteriorPredictiveStateFrequencyMonitor::swapNode(DagNode *oldN, DagNode *newN) {
	if ( oldN == data ) {
		data = static_cast< StochasticNode<BinaryCharacterData> *>( newN );
	}
    // delegate to base class
    Monitor::swapNode(oldN, newN);
}



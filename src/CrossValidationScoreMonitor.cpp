/**
 * @file
 * This file contains the implementation of CrossValidationScoreMonitor, used to save
 * information to file about the tree and additional variables at nodes, e.g. branch rates, 
 * in extended Newick format.
 *
 * @brief Implementation of CrossValidationScoreMonitor
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


#include "CrossValidationScoreMonitor.h"
#include "AbstractDiscreteCharacterData.h"
#include "DagNode.h"
#include "Model.h"
#include "Monitor.h"
#include "RbException.h"

using namespace RevBayesCore;

/* Constructor */
CrossValidationScoreMonitor::CrossValidationScoreMonitor(StochasticNode<AbstractCharacterData> *t, AbstractCharacterData* test, int g, const std::string &fname, bool ap) : Monitor(g,t), outStream(), data( t ), test( test ), filename( fname ), append(ap) {
    if(dynamic_cast<AbstractDiscreteCharacterData*>(t->getValue().clone()) == 0){
    	throw RbException("CrossValidationScoreMonitor requires DiscreteCharacterData");
    }
}


CrossValidationScoreMonitor::CrossValidationScoreMonitor(const CrossValidationScoreMonitor &m) : Monitor( m ), outStream(), data( m.data ), test(m.test){
    
    filename    = m.filename;
    append      = m.append;
    
    if (m.outStream.is_open())
        openStream();
}


/* Clone the object */
CrossValidationScoreMonitor* CrossValidationScoreMonitor::clone(void) const {
    
    return new CrossValidationScoreMonitor(*this);
}


void CrossValidationScoreMonitor::closeStream() {
    outStream.close();
}


/** Monitor value at generation gen */
void CrossValidationScoreMonitor::monitor(long gen) {
    
    // get the printing frequency
    int samplingFrequency = printgen;
    
    if (gen % samplingFrequency == 0) {
    	AbstractCharacterData* learned = data->getValue().clone();
    	data->clamp(test->clone());
    	double cvScore = data->getLnProbability();
    	outStream << gen << "\t" << cvScore << std::endl;
    	data->clamp(learned);
    }
}


/** open the file stream for printing */
void CrossValidationScoreMonitor::openStream(void) {
    
    // open the stream to the file
    if (append)
        outStream.open( filename.c_str(), std::fstream::out | std::fstream::app);
    else
        outStream.open( filename.c_str(), std::fstream::out);
}

/** Print header for monitored values */
void CrossValidationScoreMonitor::printHeader() {
	outStream << "gen\tCV" << std::endl;
}


void CrossValidationScoreMonitor::swapNode(DagNode *oldN, DagNode *newN) {
	if ( oldN == data ) {
		data = static_cast< StochasticNode<AbstractCharacterData> *>( newN );
	}
    // delegate to base class
    Monitor::swapNode(oldN, newN);
}



/**
 * @file
 * This file contains the implementation of ModelStreamMonitor, used to save
 * information to file about the tree and additional variables at nodes, e.g. branch rates,
 * in extended Newick format.
 *
 * @brief Implementation of ModelStreamMonitor
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


#include "ModelStreamMonitor.h"
#include "AbstractDiscreteCharacterData.h"
#include "DagNode.h"
#include "Model.h"
#include "Monitor.h"
#include "RbException.h"

#include <algorithm>

using namespace RevBayesCore;

/* Constructor */
ModelStreamMonitor::ModelStreamMonitor(const Model& m, int g, const std::string &fname, bool ap) : Monitor(g,m.getDagNodes()), nodeNames(), outStream(), filename( fname ), append(ap) {

	getOrderedStochasticNodes();
}


ModelStreamMonitor::ModelStreamMonitor(const ModelStreamMonitor &m) : Monitor( m ), nodeNames(m.nodeNames), outStream() {

    filename    = m.filename;
    append      = m.append;

    if (m.outStream.is_open())
        openStream();

    getOrderedStochasticNodes();
}


/* Clone the object */
ModelStreamMonitor* ModelStreamMonitor::clone(void) const {

    return new ModelStreamMonitor(*this);
}


void ModelStreamMonitor::closeStream() {
    outStream.close();
}


/** Monitor value at generation gen */
void ModelStreamMonitor::monitor(long gen) {

    // get the printing frequency
    int samplingFrequency = printgen;

    class valuetype;

    if (gen % samplingFrequency == 0) {
    	std::set<std::string>::iterator it;
    	for ( it=nodeNames.begin() ; it != nodeNames.end(); it++ ){
    		if(!nodeMap[*it]->isClamped())
    			outStream << nodeMap[*it] << "\n";
    	}
    	outStream.flush();
    }
}


/** open the file stream for printing */
void ModelStreamMonitor::openStream(void) {

    // open the stream to the file
    if (append)
        outStream.open( filename.c_str(), std::fstream::out | std::fstream::app);
    else
        outStream.open( filename.c_str(), std::fstream::out);
}

/** Creates a vector of stochastic nodes, starting from the source nodes to the sink nodes */
void ModelStreamMonitor::getOrderedStochasticNodes(void) {


	std::vector<DagNode *>::const_iterator it;
	for ( it=nodes.begin() ; it != nodes.end(); it++ ){
		if((*it)->isStochastic()){
			nodeNames.insert((*it)->getName());
			nodeMap[(*it)->getName()] = (*it);
		}
	}
}

void ModelStreamMonitor::swapNode(DagNode *oldN, DagNode *newN)
{
    // error catching
    std::set<std::string>::iterator it = find(nodeNames.begin(), nodeNames.end(), oldN->getName());

    if (it != nodeNames.end())
		nodeMap[*it] = newN;

    Monitor::swapNode(oldN,newN);
}



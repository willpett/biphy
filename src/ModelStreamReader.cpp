/**
 * @file
 * This file contains the implementation of ModelStreamReader, used to save
 * information to file about the tree and additional variables at nodes, e.g. branch rates,
 * in extended Newick format.
 *
 * @brief Implementation of ModelStreamReader
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


#include "ModelStreamReader.h"
#include "AbstractDiscreteCharacterData.h"
#include "DagNode.h"
#include "Model.h"
#include "Monitor.h"
#include "RbException.h"

#include <algorithm>

using namespace RevBayesCore;

/* Constructor */
ModelStreamReader::ModelStreamReader(const Model& m, const std::vector<Monitor*> &mons, int g, const std::string &fname) :
		model(m),
		generation(0),
		monitors(),
		nodeNames(),
		inStream(),
		filename( fname ),
		printgen( g)
{
	replaceDag(mons);
	initializeStream();
	initializeMonitors();
}


ModelStreamReader::ModelStreamReader(const ModelStreamReader &m) :
		model(m.model),
		generation(m.generation),
		monitors(),
		nodeNames(),
		inStream(),
		filename(m.filename),
		printgen( m.printgen)
{

    if (m.inStream.is_open())
        openStream();

    const std::vector<Monitor*>& mons = m.monitors;
    // create an independent copy of the model, monitors and moves
	replaceDag(mons);

	initializeStream();
	initializeMonitors();
}

/**
 * Destructor. Frees the DAG nodes (the model), monitor and the move schedule.
 */
ModelStreamReader::~ModelStreamReader(void) {

    for (std::vector<Monitor*>::iterator it = monitors.begin(); it != monitors.end(); ++it)
    {
        Monitor *theMonitor = (*it);
        delete theMonitor;
    }
}


/* Clone the object */
ModelStreamReader* ModelStreamReader::clone(void) const {

    return new ModelStreamReader(*this);
}


void ModelStreamReader::closeStream() {
    inStream.close();
}

void ModelStreamReader::run(size_t kIterations) {
	if ( generation == 0 )
    {
        // Monitor
        startMonitors();
    }

	for (int k=1; k<=kIterations; k++)
	{
		if(!nextSample())
			break;
		// Monitor
		monitor(generation);

	}
}

/** Monitor value at generation gen */
bool ModelStreamReader::nextSample() {

    // get the printing frequency
    int samplingFrequency = printgen;

    std::set<std::string>::iterator it;
	for ( it=nodeNames.begin() ; it != nodeNames.end(); it++ ){
		if(!nodeMap[*it]->isClamped()){
			if(inStream >> nodeMap[*it])
				nodeMap[*it]->keep();
			else
				return false;
		}
	}

    generation++;
    return true;
}

void ModelStreamReader::initializeMonitors(void)
{
    for (size_t i=0; i<monitors.size(); i++)
    {
        monitors[i]->setModel( &model );
    }
}

void ModelStreamReader::initializeStream(void)
{
	getOrderedStochasticNodes();

	openStream();
}

/** Creates a vector of stochastic nodes, starting from the source nodes to the sink nodes */
void ModelStreamReader::getOrderedStochasticNodes(void) {

	std::vector<DagNode *>::const_iterator it;
	for ( it=model.getDagNodes().begin() ; it != model.getDagNodes().end(); it++ ){
		if((*it)->isStochastic()){
			nodeNames.insert((*it)->getName());
			nodeMap[(*it)->getName()] = (*it);
		}
	}
}

void ModelStreamReader::monitor(unsigned long g)
{

    // Monitor
    for (size_t i = 0; i < monitors.size(); i++)
    {
    	//std::cout << monitors[i] << std::endl;
        monitors[i]->monitor(g);
    }

}

void ModelStreamReader::replaceDag(const std::vector<Monitor *> &mons)
{

    // we need to replace the DAG nodes of the monitors and moves
    const std::vector<DagNode*>& modelNodes = model.getDagNodes();

    for (std::vector<Monitor*>::const_iterator it = mons.begin(); it != mons.end(); ++it)
    {
        Monitor *theMonitor = (*it)->clone();
        std::vector<DagNode*> nodes = theMonitor->getDagNodes();
        for (std::vector<DagNode*>::const_iterator j = nodes.begin(); j != nodes.end(); ++j)
        {

            // error checking
            if ( (*j)->getName() == "" )
                throw RbException( "Unable to connect monitor to DAG copy because variable name was lost");

            DagNode* theNewNode = NULL;
            for (std::vector<DagNode*>::const_iterator k = modelNodes.begin(); k != modelNodes.end(); ++k)
            {
                if ( (*k)->getName() == (*j)->getName() )
                {
                    theNewNode = *k;
                    break;
                }
            }
            // error checking
            if ( theNewNode == NULL )
            {
                throw RbException("Cannot find node with name '" + (*j)->getName() + "' in the model but received a monitor working on it.");
            }

            // now swap the node
            theMonitor->swapNode( *j, theNewNode );
        }
        monitors.push_back( theMonitor );
    }
}


/** open the file stream for printing */
void ModelStreamReader::openStream(void) {

    // open the stream to the file
    inStream.open( filename.c_str(), std::fstream::in);
}

void ModelStreamReader::startMonitors( void ) {

    /* Open the output file and print headers */
    for (size_t i=0; i<monitors.size(); i++)
    {

        // open filestream for each monitor
        monitors[i]->openStream();
        monitors[i]->printHeader();
    }
}



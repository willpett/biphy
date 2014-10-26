/**
 * @file
 * This file contains the declaration of a BinaryDolloCompatibleMonitor, used to save information
 * to a file about the tree and additional variable on nodes (or branches) using
 * the extended Newick formay.
 *
 * @brief Declaration of BinaryDolloCompatibleMonitor
 *
 * (c) Copyright 2009- under GPL version 3
 * @date Last modified: $Date: 2012-11-26 05:34:51 -0800 (Mon, 26 Nov 2012) $
 * @author The RevBayes Development Core Team
 * @license GPL version 3
 * @version 1.0
 * @since 2012-12-11, version 1.0
 *
 * $Id: BinaryDolloCompatibleMonitor.h 1867 2012-11-26 13:34:51Z hoehna $
 */

#ifndef BinaryDolloCompatibleMonitor_H
#define BinaryDolloCompatibleMonitor_H

#include "Monitor.h"
#include "DolloBinaryCharEvoModel.h"
#include "AbstractCharacterData.h"
#include "StochasticNode.h"
#include "FastaWriter.h"

#include <fstream>
#include <iostream>
#include <string>
#include <vector>

namespace RevBayesCore {
    
	template<class treeType>
    class BinaryDolloCompatibleMonitor : public Monitor {
        
    public:
        // Constructors and Destructors
        BinaryDolloCompatibleMonitor(StochasticNode<AbstractCharacterData> *t, int g, const std::string &fname);
        //!< Constructor with set of DAG node
        BinaryDolloCompatibleMonitor(const BinaryDolloCompatibleMonitor& f);
        
        // basic methods
        BinaryDolloCompatibleMonitor*       clone(void) const;                                                  //!< Clone the object
        
        // Monitor functions
        void                                monitor(long gen);                                                  //!< Monitor at generation gen
        void                                swapNode(DagNode *oldN, DagNode *newN);

        // FileMonitor functions
        void                                closeStream(void);                                                  //!< Close stream after finish writing
        void                                openStream(void);                                                   //!< Open the stream for writing
        void                                printHeader(void);                                                  //!< Print header
        
    private:        
        
        // parameters
        StochasticNode<AbstractCharacterData>*      data;
        std::string                         		filename;
        
    };

template<class treeType>
BinaryDolloCompatibleMonitor<treeType>::BinaryDolloCompatibleMonitor(StochasticNode<AbstractCharacterData> *t, int g, const std::string &fname) : Monitor(g,t), data( t ), filename( fname ) {
    if(dynamic_cast<BinaryCharEvoModel<treeType>*>(t->getDistribution().clone()) == 0){
    	throw RbException("BinaryDolloCompatibleMonitor requires BinaryCharEvoModel");
    }

    if(dynamic_cast<DolloBinaryCharEvoModel<treeType>*>(t->getDistribution().clone())){
		throw RbException("BinaryDolloCompatibleMonitor requires BinaryCharEvoModel");
	}
}

template<class treeType>
BinaryDolloCompatibleMonitor<treeType>::BinaryDolloCompatibleMonitor(const BinaryDolloCompatibleMonitor &m) : Monitor( m ), data( m.data ){

    filename    = m.filename;
}


/* Clone the object */
template<class treeType>
BinaryDolloCompatibleMonitor<treeType>* BinaryDolloCompatibleMonitor<treeType>::clone(void) const {

    return new BinaryDolloCompatibleMonitor<treeType>(*this);
}

template<class treeType>
void BinaryDolloCompatibleMonitor<treeType>::closeStream() {
}


/** Monitor value at generation gen */
template<class treeType>
void BinaryDolloCompatibleMonitor<treeType>::monitor(long gen) {

    // get the printing frequency
    int samplingFrequency = printgen;

    std::stringstream str;
    str << gen;
    if (gen % samplingFrequency == 0) {
    	DiscreteTaxonData<StandardState>* seq = dynamic_cast<BinaryCharEvoModel<treeType>&>(data->getDistribution()).getDolloCompatible();
    	seq->setTaxonName(str.str());
    	DiscreteCharacterData<StandardState> data;
    	data.addTaxonData(*seq);

    	FastaWriter writer;
    	writer.writeData(filename,data,true);
    }
}


/** open the file stream for printing */
template<class treeType>
void BinaryDolloCompatibleMonitor<treeType>::openStream(void) {

}

/** Print header for monitored values */
template<class treeType>
void BinaryDolloCompatibleMonitor<treeType>::printHeader() {
}


template<class treeType>
void BinaryDolloCompatibleMonitor<treeType>::swapNode(DagNode *oldN, DagNode *newN) {
	if ( oldN == data ) {
		data = static_cast< StochasticNode<AbstractCharacterData> *>( newN );
	}
    // delegate to base class
    Monitor::swapNode(oldN, newN);
}

}

#endif


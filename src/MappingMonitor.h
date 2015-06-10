/**
 * @file
 * This file contains the declaration of a MappingMonitor, used to save information
 * to a file about the tree and additional variable on nodes (or branches) using
 * the extended Newick formay.
 *
 * @brief Declaration of MappingMonitor
 *
 * (c) Copyright 2009- under GPL version 3
 * @date Last modified: $Date: 2012-11-26 05:34:51 -0800 (Mon, 26 Nov 2012) $
 * @author The RevBayes Development Core Team
 * @license GPL version 3
 * @version 1.0
 * @since 2012-12-11, version 1.0
 *
 * $Id: MappingMonitor.h 1867 2012-11-26 13:34:51Z hoehna $
 */

#ifndef MappingMonitor_H
#define MappingMonitor_H

#include "Monitor.h"
#include "AbstractCharacterData.h"
#include "StochasticNode.h"

#include "GeneralCharEvoModel.h"
#include "AbstractDiscreteCharacterData.h"
#include "RbException.h"

#include <fstream>
#include <iostream>
#include <string>
#include <vector>

namespace RevBayesCore {

	template<class charType, class treeType>
    class MappingMonitor : public Monitor {
        
    public:
        // Constructors and Destructors
        MappingMonitor(StochasticNode<AbstractCharacterData> *t, int g, const std::string &fname, bool ap=false);
        //!< Constructor with set of DAG node
        MappingMonitor(const MappingMonitor& f);
        
        // basic methods
        MappingMonitor*          clone(void) const;                                                  //!< Clone the object
        
        // Monitor functions
        void                                monitor(long gen);                                                  //!< Monitor at generation gen
        void                                swapNode(DagNode *oldN, DagNode *newN);

        // FileMonitor functions
        void                                closeStream(void);                                                  //!< Close stream after finish writing
        void                                openStream(void);                                                   //!< Open the stream for writing
        void                                printHeader(void);                                                  //!< Print header
        
    private:        
        // the stream to print
        std::fstream                        outStream;
        
        // parameters
        StochasticNode<AbstractCharacterData>*      data;
        std::string                         filename;
        bool                                append;
        
    };

}

/* Constructor */
template<class charType, class treeType>
RevBayesCore::MappingMonitor<charType, treeType>::MappingMonitor(StochasticNode<AbstractCharacterData> *t, int g, const std::string &fname, bool ap) : Monitor(g,t), outStream(), data( t ), filename( fname ), append(ap) {
    if(dynamic_cast<AbstractDiscreteCharacterData*>(&(t->getValue())) == 0){
    	throw RbException("MappingMonitor requires DiscreteCharacterData");
    }
    if(dynamic_cast<GeneralCharEvoModel<charType, treeType>*>(&(t->getDistribution())) == 0){
		throw RbException("MappingMonitor requires GeneralCharEvolModel");
	}
}

template<class charType, class treeType>
RevBayesCore::MappingMonitor<charType, treeType>::MappingMonitor(const MappingMonitor &m) : Monitor( m ), outStream(), data( m.data ) {

    filename    = m.filename;
    append      = m.append;

    if (m.outStream.is_open())
        openStream();
}


/* Clone the object */
template<class charType, class treeType>
RevBayesCore::MappingMonitor<charType, treeType>* RevBayesCore::MappingMonitor<charType, treeType>::clone(void) const {

    return new MappingMonitor(*this);
}

template<class charType, class treeType>
void RevBayesCore::MappingMonitor<charType, treeType>::closeStream() {
    outStream.close();
}


/** Monitor value at generation gen */
template<class charType, class treeType>
void RevBayesCore::MappingMonitor<charType, treeType>::monitor(long gen) {

    // get the printing frequency
    int samplingFrequency = printgen;

    if (gen % samplingFrequency == 0) {
    	GeneralCharEvoModel<charType,treeType> model = dynamic_cast<GeneralCharEvoModel<charType, treeType>& >(data->getDistribution());
    	const std::vector< DiscreteTaxonData<charType> >& mapping = model.getMapping();
    	treeType* tree = model.getTree()->clone();

    	std::vector<TopologyNode*> nodes = tree->getNodes();
    	for(size_t i = 0; i < nodes.size(); i++){
    		std::string name = nodes[i]->getName();
    		name += "[&sequence=";
    		const RevBayesCore::DiscreteTaxonData<charType>& d = mapping[i];
    		for(size_t c = 0; c < d.size(); c++){
    			name += d[c].getStringValue();
    		}
    		name += "]";
    		nodes[i]->setName(name);
    	}

    	tree->clearBranchParameters();
    	outStream << tree->getNewickRepresentation() << std::endl;
    	delete tree;
    }
}


/** open the file stream for printing */
template<class charType, class treeType>
void RevBayesCore::MappingMonitor<charType, treeType>::openStream(void) {

    // open the stream to the file
    if (append)
        outStream.open( filename.c_str(), std::fstream::out | std::fstream::app);
    else
        outStream.open( filename.c_str(), std::fstream::out);
}

/** Print header for monitored values */
template<class charType, class treeType>
void RevBayesCore::MappingMonitor<charType, treeType>::printHeader() {
	//for(size_t i=0;i<taxonNames.size();i++){
	//	outStream << taxonNames[i] << "\t";
	//}
	//outStream << std::endl;
}

template<class charType, class treeType>
void RevBayesCore::MappingMonitor<charType, treeType>::swapNode(DagNode *oldN, DagNode *newN) {
	if ( oldN == data ) {
		data = static_cast< StochasticNode<AbstractCharacterData> *>( newN );
	}
    // delegate to base class
    Monitor::swapNode(oldN, newN);
}

#endif


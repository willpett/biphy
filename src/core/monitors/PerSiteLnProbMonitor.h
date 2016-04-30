/**
 * @file
 * This file contains the declaration of a PerSiteLnProbMonitor, used to save information
 * to a file about the tree and additional variable on nodes (or branches) using
 * the extended Newick formay.
 *
 * @brief Declaration of PerSiteLnProbMonitor
 *
 * (c) Copyright 2009- under GPL version 3
 * @date Last modified: $Date: 2012-11-26 05:34:51 -0800 (Mon, 26 Nov 2012) $
 * @author The RevBayes Development Core Team
 * @license GPL version 3
 * @version 1.0
 * @since 2012-12-11, version 1.0
 *
 * $Id: PerSiteLnProbMonitor.h 1867 2012-11-26 13:34:51Z hoehna $
 */

#ifndef PerSiteLnProbMonitor_H
#define PerSiteLnProbMonitor_H

#include "Monitor.h"
#include "BinaryCharacterData.h"
#include "StochasticNode.h"
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include "Exception.h"


class PerSiteLnProbMonitor : public Monitor {
    
public:
    // Constructors and Destructors
    PerSiteLnProbMonitor(StochasticNode<BinaryCharacterData > *t, int g, const std::string &fname, bool ap=false);
    //!< Constructor with set of DAG node
    PerSiteLnProbMonitor(const PerSiteLnProbMonitor& f);
    
    // basic methods
    PerSiteLnProbMonitor*          clone(void) const;                                                  //!< Clone the object
    
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
    StochasticNode<BinaryCharacterData >*    data;
    std::string                                        filename;
    bool                                               append;
    
};


PerSiteLnProbMonitor::PerSiteLnProbMonitor(StochasticNode<BinaryCharacterData > *t, int g, const std::string &fname, bool ap) : 
    Monitor(g,t), outStream(), data( t ), filename( fname ), append(ap) 
{
    if(dynamic_cast<BinarySubstitutionModel *>(&(t->getDistribution())) == 0){
        throw Exception("PerSiteLnProbMonitor requires BinarySubstitutionModel");
    }
    
    dynamic_cast<BinarySubstitutionModel& >(data->getDistribution()).setVerbose(false);
}



PerSiteLnProbMonitor::PerSiteLnProbMonitor(const PerSiteLnProbMonitor &m) : 
        Monitor( m ), outStream(), data( m.data )
{

    filename    = m.filename;
    append      = m.append;

    if (m.outStream.is_open())
        openStream();
}


/* Clone the object */


PerSiteLnProbMonitor* PerSiteLnProbMonitor::clone(void) const {

    return new PerSiteLnProbMonitor(*this);
}



void PerSiteLnProbMonitor::closeStream() {
    outStream.close();
}


/** Monitor value at generation gen */


void PerSiteLnProbMonitor::monitor(long gen) {

    // get the printing frequency
    int samplingFrequency = printgen;

    if (gen % samplingFrequency == 0) {
        BinarySubstitutionModel model1 = dynamic_cast<BinarySubstitutionModel& >(data->getDistribution());
        
        RealVector probs1 = model1.getPerSiteLnProbs();

        size_t numSites = model1.getNumPatterns();
       
        for(size_t i = 0; i < numSites; i++){
            outStream << probs1[i];
            if(i != numSites - 1)
                outStream << "\t";
        }
        outStream << std::endl;
    }
}


/** open the file stream for printing */


void PerSiteLnProbMonitor::openStream(void) {

    // open the stream to the file
    if (append)
        outStream.open( filename.c_str(), std::fstream::out | std::fstream::app);
    else
        outStream.open( filename.c_str(), std::fstream::out);
}

/** Print header for monitored values */


void PerSiteLnProbMonitor::printHeader() {
    BinarySubstitutionModel model = dynamic_cast<BinarySubstitutionModel& >(data->getDistribution());

    size_t numSites = model.getNumPatterns();
    for(size_t i=0; i<numSites; i++)
    {
        outStream << i;
        if(i != numSites - 1)
            outStream << "\t";
    }
    outStream << std::endl;
}


void PerSiteLnProbMonitor::swapNode(DagNode *oldN, DagNode *newN) {
    if ( oldN == data ) {
        data = static_cast< StochasticNode<BinaryCharacterData > *>( newN );
    }
    // delegate to base class
    Monitor::swapNode(oldN, newN);
}

#endif


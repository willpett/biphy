/**
 * @file
 * This file contains the implementation of FileMonitor, used to save
 * information to file about the monitoring of a variable DAG node.
 *
 * @brief Implementation of FileMonitor
 *
 * (c) Copyright 2009- under GPL version 3
 * @date Last modified: $Date$
 * @author The RevBayes Development Core Team
 * @license GPL version 3
 * @version 1.0
 * @since 2012-06-21, version 1.0
 *
 * $Id$
 */


#include "FileMonitor.h"
#include "DagNode.h"
#include "Mcmc.h"
#include "Model.h"
#include "Monitor.h"

using namespace RevBayesCore;

/* Constructor */
FileMonitor::FileMonitor(DagNode *n, int g, const std::string &fname, const std::string &del, bool pp, bool l, bool pr, bool ap, bool ci, bool ch) : Monitor(g,n), outStream(), filename( fname ), separator( del ), posterior( pp ), prior( pr ), likelihood( l ), append(ap), chainIdx(ci), chainHeat(ch) {
    
}


/* Constructor */
FileMonitor::FileMonitor(const std::set<DagNode *> &n, int g, const std::string &fname, const std::string &del, bool pp, bool l, bool pr, bool ap, bool ci, bool ch) : Monitor(g,n), outStream(), filename( fname ), separator( del ), posterior( pp ), prior( pr ), likelihood( l ), append(ap), chainIdx(ci), chainHeat(ch) {
    
}

FileMonitor::FileMonitor(const std::vector<DagNode *> &n, int g, const std::string &fname, const std::string &del, bool pp, bool l, bool pr, bool ap, bool ci, bool ch) : Monitor(g,n), outStream(), filename( fname ), separator( del ), posterior( pp ), prior( pr ), likelihood( l ), append(ap), chainIdx(ci), chainHeat(ch) {
    
}


FileMonitor::FileMonitor(const FileMonitor &f) : Monitor( f ), outStream() {
    
    filename    = f.filename;
    separator   = f.separator;
    prior       = f.prior;
    posterior   = f.posterior;
    likelihood  = f.likelihood;
    append      = f.append;
    chainIdx    = f.chainIdx;
    chainHeat   = f.chainHeat;
    
    if (f.outStream.is_open())
        openStream();
}


/* Clone the object */
FileMonitor* FileMonitor::clone(void) const {

    return new FileMonitor(*this);
}


void FileMonitor::closeStream() {
    outStream.close();
}


/** Monitor value at generation gen */
void FileMonitor::monitor(long gen) {

    // get the printing frequency
    long samplingFrequency = printgen;
    
    if (gen % samplingFrequency == 0) {
        // print the iteration number first
        outStream << gen;
        
        if ( posterior ) {
            // add a separator before every new element
            outStream << separator;
            
            const std::vector<DagNode*> &n = model->getDagNodes();
            double pp = 0.0;
            for (std::vector<DagNode*>::const_iterator it = n.begin(); it != n.end(); ++it)
            {
                pp += (*it)->getLnProbability();
            }
            outStream << pp;
        }
        
        if ( likelihood ) {
            // add a separator before every new element
            outStream << separator;
            
            const std::vector<DagNode*> &n = model->getDagNodes();
            double pp = 0.0;
            for (std::vector<DagNode*>::const_iterator it = n.begin(); it != n.end(); ++it) {
                if ( (*it)->isClamped() ) {
                    pp += (*it)->getLnProbability();
                }
            }
            outStream << pp;
        }
        
        if ( prior ) {
            // add a separator before every new element
            outStream << separator;
            
            const std::vector<DagNode*> &n = model->getDagNodes();
            double pp = 0.0;
            for (std::vector<DagNode*>::const_iterator it = n.begin(); it != n.end(); ++it) {
                if ( !(*it)->isClamped() ) {
                    pp += (*it)->getLnProbability();
                }
            }
            outStream << pp;
        }
        
        if ( chainIdx ) {
            // add a separator before every new element
            outStream << separator;
            outStream << mcmc->getChainIndex();
            
        }
        
        if (chainHeat)
        {
            outStream << separator;
            outStream << mcmc->getChainHeat();
        }
        
        for (std::vector<DagNode*>::iterator i = nodes.begin(); i != nodes.end(); ++i) {
            // add a separator before every new element
            outStream << separator;
            
            // get the node
            DagNode *node = *i;
            
            // print the value
            node->printValue(outStream, separator);
        }

        outStream << std::endl;
    
    }
}


/** open the file stream for printing */
void FileMonitor::openStream(void) {
    
    // open the stream to the file
    if (append)
        outStream.open( filename.c_str(), std::fstream::out | std::fstream::app);
    else
        outStream.open( filename.c_str(), std::fstream::out);    
}

/** Print header for monitored values */
void FileMonitor::printHeader() {
  
    // print one column for the iteration number
    outStream << "gen";
    
    if ( posterior ) {
        // add a separator before every new element
        outStream << separator;
        outStream << "Posterior";
    }
    
    if ( likelihood ) {
        // add a separator before every new element
        outStream << separator;
        outStream << "lnL";
    }
    
    if ( prior ) {
        // add a separator before every new element
        outStream << separator;
        outStream << "pr";
    }
    
    if ( chainIdx ) {
        outStream << separator;
        outStream << "chIdx";
    }
    
    if (chainHeat)
    {
        outStream << separator;
        outStream << "chBeta";
    }
    
    for (std::vector<DagNode *>::const_iterator it=nodes.begin(); it!=nodes.end(); it++) {
        // add a separator before every new element
        outStream << separator;
        
         const DagNode* theNode = *it;
        
        // print the header
        if (theNode->getName() != "")
            theNode->printName(outStream,separator);
        else
            outStream << "Unnamed";
    }
    
    outStream << std::endl;
}



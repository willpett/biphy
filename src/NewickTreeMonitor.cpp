/**
 * @file
 * This file contains the implementation of NewickTreeMonitor, used to save
 * information to file about the tree and additional variables at nodes, e.g. branch rates, 
 * in extended Newick format.
 *
 * @brief Implementation of NewickTreeMonitor
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


#include "NewickTreeMonitor.h"
#include "DagNode.h"
#include "Model.h"
#include "Monitor.h"

#include "Exception.h"

/* Constructor */
NewickTreeMonitor::NewickTreeMonitor(TypedDagNode<Tree> *t, int g, const std::string &fname, bool ap) : Monitor(g,t), outStream(), tree( t ), filename( fname ), append(ap) {
    
}


/* Constructor */
NewickTreeMonitor::NewickTreeMonitor(TypedDagNode<Tree> *t, const std::set<TypedDagNode< std::vector<double> > *> &n, int g, const std::string &fname, bool ap) : Monitor(g,t), outStream(), tree( t ), nodeVariables( n ), filename( fname ), append(ap) {
//    this->nodes.insert( tree );
    
    for (std::set<TypedDagNode< std::vector<double> > *>::iterator it = nodeVariables.begin(); it != nodeVariables.end(); ++it) {
        this->nodes.push_back( *it );
    }
}


NewickTreeMonitor::NewickTreeMonitor(const NewickTreeMonitor &m) : Monitor( m ), outStream(), tree( m.tree ), nodeVariables( m.nodeVariables ) {
    
    filename    = m.filename;
    append      = m.append;
    
    //if (m.outStream.is_open())
    //    openStream();
}


/* Clone the object */
NewickTreeMonitor* NewickTreeMonitor::clone(void) const {
    
    return new NewickTreeMonitor(*this);
}


void NewickTreeMonitor::closeStream() {
    outStream.close();
}


/** Monitor value at generation gen */
void NewickTreeMonitor::monitor(long gen) {
    
    // get the printing frequency
    int samplingFrequency = printgen;
    
    if (gen % samplingFrequency == 0) {
        tree->getValue().clearBranchParameters();
        //for (std::set<TypedDagNode< std::vector<double> > *>::iterator it = nodeVariables.begin(); it != nodeVariables.end(); ++it) {
        //    tree->getValue().addBranchParameter((*it)->getName(), (*it)->getValue(), false);
        //}
            
        outStream << tree->getValue().getNewickRepresentation();
        outStream << ";" << std::endl;
        
    }
}


/** open the file stream for printing */
void NewickTreeMonitor::openStream(void) {
    
    // open the stream to the file
    if (append)
        outStream.open( filename.c_str(), std::fstream::out | std::fstream::app);
    else
        outStream.open( filename.c_str(), std::fstream::out);
}

/** Print header for monitored values */
void NewickTreeMonitor::printHeader() {
}


void NewickTreeMonitor::swapNode(DagNode *oldN, DagNode *newN) {
    
    TypedDagNode<std::vector<double> >* nodeVar = dynamic_cast< TypedDagNode<std::vector<double> > *>(oldN);
    if ( oldN == tree ) {
        tree = static_cast< TypedDagNode< Tree > *>( newN );
    } else if ( nodeVar != NULL ) {
        // error catching
        if ( nodeVariables.find(nodeVar) == nodeVariables.end() ) {
            throw Exception("Cannot replace DAG node with name\"" + oldN->getName() + "\" in this extended newick monitor because the monitor doesn't hold this DAG node.");
        }
        
        nodeVariables.erase( nodeVar );
        nodeVariables.insert( static_cast< TypedDagNode<std::vector<double> > *>(newN) );
    }
    
    // delegate to base class
    Monitor::swapNode(oldN, newN);
}



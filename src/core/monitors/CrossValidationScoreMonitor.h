/**
 * @file
 * This file contains the declaration of a CrossValidationScoreMonitor, used to save information
 * to a file about the tree and additional variable on nodes (or branches) using
 * the extended Newick formay.
 *
 * @brief Declaration of CrossValidationScoreMonitor
 *
 * (c) Copyright 2009- under GPL version 3
 * @date Last modified: $Date: 2012-11-26 05:34:51 -0800 (Mon, 26 Nov 2012) $
 * @author The RevBayes Development Core Team
 * @license GPL version 3
 * @version 1.0
 * @since 2012-12-11, version 1.0
 *
 * $Id: CrossValidationScoreMonitor.h 1867 2012-11-26 13:34:51Z hoehna $
 */

#ifndef CrossValidationScoreMonitor_H
#define CrossValidationScoreMonitor_H

#include "Monitor.h"
#include "BinaryCharacterData.h"
#include "StochasticNode.h"

#include <fstream>
#include <iostream>
#include <string>
#include <vector>

class CrossValidationScoreMonitor : public Monitor {
    
public:
    // Constructors and Destructors
    CrossValidationScoreMonitor(StochasticNode<BinaryCharacterData > *t, BinaryCharacterData* test, int g, const std::string &fname, bool ap=false);
    //!< Constructor with set of DAG node
    CrossValidationScoreMonitor(const CrossValidationScoreMonitor& f);
    
    // basic methods
    CrossValidationScoreMonitor*          clone(void) const;                                                  //!< Clone the object
    
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
    StochasticNode<BinaryCharacterData >*      data;
    BinaryCharacterData*				test;
    std::string                         filename;
    bool                                append;
    
};

#endif


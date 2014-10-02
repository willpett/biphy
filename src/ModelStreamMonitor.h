/**
 * @file
 * This file contains the declaration of a ModelStreamMonitor, used to save information
 * to a file about the tree and additional variable on nodes (or branches) using
 * the extended Newick formay.
 *
 * @brief Declaration of ModelStreamMonitor
 *
 * (c) Copyright 2009- under GPL version 3
 * @date Last modified: $Date: 2012-11-26 05:34:51 -0800 (Mon, 26 Nov 2012) $
 * @author The RevBayes Development Core Team
 * @license GPL version 3
 * @version 1.0
 * @since 2012-12-11, version 1.0
 *
 * $Id: ModelStreamMonitor.h 1867 2012-11-26 13:34:51Z hoehna $
 */

#ifndef ModelStreamMonitor_H
#define ModelStreamMonitor_H

#include "Monitor.h"
#include "AbstractCharacterData.h"
#include "StochasticNode.h"

#include <fstream>
#include <iostream>
#include <string>
#include <vector>

namespace RevBayesCore {

    class ModelStreamMonitor : public Monitor {

    public:
        // Constructors and Destructors
        ModelStreamMonitor(const Model& m, int g, const std::string &fname, bool ap=false);
        ModelStreamMonitor(const ModelStreamMonitor& f);

        // basic methods
        ModelStreamMonitor*          		clone(void) const;                                                  //!< Clone the object
        virtual void                      	swapNode(DagNode *oldN, DagNode *newN);

        // Monitor functions
        void                                monitor(long gen);                                                  //!< Monitor at generation gen

        // FileMonitor functions
        void                                closeStream(void);                                                  //!< Close stream after finish writing
        void                                openStream(void);                                                   //!< Open the stream for writing

    private:
        void                                getOrderedStochasticNodes(void);

        std::set<std::string>				nodeNames;
        std::map<std::string,DagNode*>		nodeMap;
        // the stream to print
        std::fstream                        outStream;

        // parameters
        std::string                         filename;
        bool                                append;

    };

}

#endif


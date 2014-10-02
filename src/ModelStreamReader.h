/**
 * @file
 * This file contains the declaration of a ModelStreamReader, used to save information
 * to a file about the tree and additional variable on nodes (or branches) using
 * the extended Newick formay.
 *
 * @brief Declaration of ModelStreamReader
 *
 * (c) Copyright 2009- under GPL version 3
 * @date Last modified: $Date: 2012-11-26 05:34:51 -0800 (Mon, 26 Nov 2012) $
 * @author The RevBayes Development Core Team
 * @license GPL version 3
 * @version 1.0
 * @since 2012-12-11, version 1.0
 *
 * $Id: ModelStreamReader.h 1867 2012-11-26 13:34:51Z hoehna $
 */

#ifndef ModelStreamReader_H
#define ModelStreamReader_H

#include "Model.h"
#include "Monitor.h"
#include "AbstractCharacterData.h"
#include "StochasticNode.h"

#include <fstream>
#include <iostream>
#include <string>
#include <vector>

namespace RevBayesCore {

    class ModelStreamReader : Cloneable {

    public:
        // Constructors and Destructors
        ModelStreamReader(const Model& m, const std::vector<Monitor*> &mons, int g, const std::string &fname);
        ModelStreamReader(const ModelStreamReader& f);
        virtual                             ~ModelStreamReader(void);

        // basic methods
        ModelStreamReader*          		clone(void) const;                                                  //!< Clone the object
        void                                run(size_t g);
        void                                monitor(unsigned long g);

        // Monitor functions
        bool                                nextSample();                                                  //!< Monitor at generation gen

        // FileMonitor functions
        void                                closeStream(void);                                                  //!< Close stream after finish writing
        void                                openStream(void);                                                   //!< Open the stream for writing

    protected:
        void                                getOrderedStochasticNodes(void);
        void                                initializeStream(void);                                                                  //!< Initialize objects for mcmc sampling
		void                                initializeMonitors(void);                                                               //!< Assign model and mcmc ptrs to monitors
		void                                replaceDag(const std::vector<Monitor*> &mons);
		void                                startMonitors(void);

    private:

        Model                               model;
        size_t								generation;

        std::set<std::string>				nodeNames;
		std::map<std::string,DagNode*>		nodeMap;
        std::fstream                        inStream;

        // parameters
        std::string                         filename;
        int									printgen;
        std::vector<Monitor*>               monitors;

    };

}

#endif


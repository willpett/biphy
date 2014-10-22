/**
 * @file
 * This file contains the declaration of the Gtr model test class. 
 * This test runs an MCMC analysis on a phylogenetic model using the GTR subsitution model.
 * The data is read as an alignment from file. The tree is estimated as well!
 *
 *
 * @brief Declaration of the Gtr+G model test class
 *
 * (c) Copyright 2009-
 * @date Last modified: $Date$
 * @author The RevBayes Development Core Team
 * @license GPL version 3
 * @since Version 1.0, 2012-07-05
 *
 * $Id$
 */


#ifndef Biphy_H
#define Biphy_H

#include <string>
#include <vector>
    
    class Biphy {
        
    public:
    	Biphy(const std::string &datafile,
    										const std::string &treefile,
    										const std::string &name,
    										const std::string &outgroupfile,
    										int branchprior,
    										bool ras,
    										int heterogeneous,
    										bool dollo,
    										int mixture,
    										bool rigidroot,
    										bool rootprior,
    										double rootmin,
    										double rootmax,
    										int every,
    										int until,
    										int numchains,
    										int swapInterval,
    										double deltaTemp,
    										double sigmaTemp,
    										bool saveall,
    										bool nexus,
											int correctionType);
    	Biphy(const std::string &name, const std::string &cvfile, bool ppred);
    	Biphy(const std::string &name);
        virtual                                ~Biphy(void);                                                            //!< Virtual destructor
        
        bool                                    run();
        void                                    open();
        void                                    save();
        
    private:
        
        // members
        std::string                             dataFile;
        std::string								name;
        std::string                             treeFile;
        std::string                             outgroupFile;
        std::string                             cvfile;

        int										branchprior;
        bool									ras;
        int										heterogeneous;
        bool									dollo;
        int 									mixture;
        bool									ppred;
        bool									rigidroot;
        bool									rootprior;
        double                                  rootmin;
        double                                  rootmax;
        int                                     every;
        int                                     until;
        int                                     numChains;
        int                                     swapInterval;
        double                                  deltaTemp;
        double                                  sigmaTemp;
        bool									saveall;
        bool									nexus;
        int										correctionType;
        bool									readstream;
        bool									restart;
        
    };

#endif

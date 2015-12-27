/**
 * @file
 * This file contains the declaration of BinaryCharacterDataReader, which is the wrapper class for the NCL
 * for reading character matrices, sequences, and trees.
 
 * @brief Declaration of BinaryCharacterDataReader
 *
 * (c) Copyright 2009- under GPL version 3
 * @date Last modified: $Date$
 * @author The RevBayes Development Core Team
 * @license GPL version 3
 * @version 1.0
 * @since 2009-08-16, version 1.0
 * @extends DAGNode
 *
 * $Id$
 */

#ifndef BinaryCharacterDataReader_H
#define BinaryCharacterDataReader_H

#include "BinaryCharacterData.h"
#include "ContinuousBinaryCharacterData.h"
#include "TopologyNode.h"

#include <map>
#include <set>
#include <string>
#include <vector>
#include "FileManager.h"

class Tree;

class BinaryCharacterDataReader{
    
public:
    
    void                                        addWarning(std::string s) { warningsSummary.insert(s); }                        //!< Add a warning to the warnings vector
    void                                        clearWarnings(void) { warningsSummary.clear(); }                                //!< Clear all of the warnings from the warnings vector
    size_t                                      getNumWarnings(void) { return warningsSummary.size(); }                         //!< Return the number of warnings
    std::set<std::string>&                      getWarnings(void) { return warningsSummary; }                                   //!< Get a reference to the warnings vector
    static BinaryCharacterDataReader&           getInstance(void);
    
    // file type methods
    bool                                        isFastaFile(std::string& fn);
    bool                                        isPastaFile(std::string& fn);
    bool                                        isNexusFile(const std::string& fn);                                             //!< Checks if the file is in NEXUS format
    bool                                        isPhylipFile(std::string& fn, bool& isInterleaved);         //!< Checks if the file is in Phylip format
    
    // alignment stuff
    BinaryCharacterData*                        readMatrix(std::string fn);
    BinaryCharacterData*                        ReadNexus(std::string filespec);
    BinaryCharacterData*                        ReadPhylipSequential(std::string filespec);
    BinaryCharacterData*                        ReadFasta(std::string filespec);
    BinaryCharacterData*                        ReadPhylip(std::string filespec, int repeattaxa);
    ContinuousBinaryCharacterData*              ReadPasta(std::string filespec);
private:
    BinaryCharacterDataReader(void) { }                                                             //!< Default constructor
    BinaryCharacterDataReader(const BinaryCharacterDataReader& r) { }                                               //!< Copy constructor
    virtual                                    ~BinaryCharacterDataReader(void) { }                                                             //!< Destructor

    bool                                        fileExists(const char *fn) const;                                               //!< Returns whether a file exists
    std::string                                 findFileNameFromPath(const std::string& fp) const;                              //!< Returns the file name from a file path
    void                                        formatError(FileManager& fm, std::string& errorStr);                        //!< Attempt to determine the type of data
    
    std::set<std::string>                       warningsSummary;                                                                //!< A vector that contains the warnings that acumulate
};

#endif

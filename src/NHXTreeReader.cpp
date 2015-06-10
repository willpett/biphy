#include "BranchLengthTree.h"
#include "RbException.h"

#include <fstream>
#include "NHXTreeReader.h"
#include "NHXConverter.h"

using namespace RevBayesCore;


/**
 * Default constructor.
 * 
 * The default constructor does nothing except allocating the object.
 */
NHXTreeReader::NHXTreeReader()
{
    
}



/**
 *
 */
std::vector<BranchLengthTree*>* NHXTreeReader::readBranchLengthTrees(std::string const &fn)
{
    /* Open file */
    std::ifstream inFile( fn.c_str() );
    
    if ( !inFile )
        throw RbException( "Could not open file \"" + fn + "\"" );
    
    /* Initialize */
    std::vector<BranchLengthTree*>* trees = new std::vector<BranchLengthTree*>();
    std::string commandLine;
    
    /* line-processing loop */
    while ( inFile.good() ) 
    {
        
        // Read a line
        std::string line;
        getline( inFile, line );
        
        // skip empty lines
        if (line.length() == 0) 
        {
            continue;
        }
                
        
        NHXConverter c;
        BranchLengthTree *blTree = c.convertFromNHX( line );

        trees->push_back( blTree );
    }
    
    return trees;
}

std::vector<Topology*>* NHXTreeReader::readTopologies(std::string const &fn)
{
    /* Open file */
    std::ifstream inFile( fn.c_str() );

    if ( !inFile )
        throw RbException( "Could not open file \"" + fn + "\"" );

    /* Initialize */
    std::vector<Topology*>* trees = new std::vector<Topology*>();
    std::string commandLine;

    /* line-processing loop */
    while ( inFile.good() )
    {

        // Read a line
        std::string line;
        getline( inFile, line );

        // skip empty lines
        if (line.length() == 0)
        {
            continue;
        }


        NHXConverter c;
        Topology *blTree = c.convertTopologyFromNHX( line );

        trees->push_back( blTree );
    }

    return trees;
}

#include "Tree.h"
#include "NewickConverter.h"
#include "NewickTreeReader.h"
#include <fstream>
#include "Exception.h"

/**
 * Default constructor.
 * 
 * The default constructor does nothing except allocating the object.
 */
NewickTreeReader::NewickTreeReader() 
{
    
}



/**
 *
 */
std::vector<Tree*>* NewickTreeReader::readTrees(std::string const &fn) 
{
    /* Open file */
    std::ifstream inFile( fn.c_str() );
    
    if ( !inFile )
        throw Exception( "Could not open file \"" + fn + "\"" );
    
    /* Initialize */
    std::vector<Tree*>* trees = new std::vector<Tree*>();
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
                
        
        NewickConverter c;
        Tree *blTree = c.convertFromNewick( line );

        trees->push_back( blTree );
    }
    
    return trees;
}

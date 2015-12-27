#include "BinaryCharacterData.h"
#include "FastaWriter.h"

/**
 * Default constructor.
 * 
 * The default constructor does nothing except allocating the object.
 */
FastaWriter::FastaWriter( void ) 
{
    
}


/**
 * This method simply writes a character data object into a file in Fasta format.
 *
 * \param[in]   fileName    The name of the file into which the objects is to be written.
 * \param[in]   data        The character data object which is written out.
 */
void FastaWriter::writeData(std::string const &fileName, const BinaryCharacterData &data, bool append)
{
    
    // the filestream object
    std::fstream outStream;
    
    // open the stream to the file
    if(append)
    	outStream.open( fileName.c_str(), std::fstream::out | std::fstream::app);
    else
    	outStream.open( fileName.c_str(), std::fstream::out );
    
    const std::vector<std::string> &taxonNames = data.getTaxonNames();
    for (std::vector<std::string>::const_iterator it = taxonNames.begin();  it != taxonNames.end(); ++it) 
    {
        outStream << ">" << *it << std::endl;
        const BinaryTaxonData* taxon = data.getTaxonData( *it );
        size_t nChars = taxon->getNumberOfCharacters();
        for (size_t i = 0; i < nChars; ++i) 
        {
            outStream << taxon->getCharacter( i );
        }
        outStream << std::endl;
    }
    
    // close the stream
    outStream.close();
}

#include "AbstractTaxonData.h"
#include "CharacterState.h"
#include "FastaWriter.h"

using namespace RevBayesCore;


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
void FastaWriter::writeData(std::string const &fileName, const AbstractCharacterData &data) 
{
    
    // the filestream object
    std::fstream outStream;
    
    // open the stream to the file
    outStream.open( fileName.c_str(), std::fstream::out );
    
    const std::vector<std::string> &taxonNames = data.getTaxonNames();
    for (std::vector<std::string>::const_iterator it = taxonNames.begin();  it != taxonNames.end(); ++it) 
    {
        outStream << ">" << *it << std::endl;
        const AbstractTaxonData &taxon = data.getTaxonData( *it );
        size_t nChars = taxon.getNumberOfCharacters();
        for (size_t i = 0; i < nChars; ++i) 
        {
            const CharacterState &c = taxon.getCharacter( i );  
            outStream << c.getStringValue();
        }
        outStream << std::endl;
    }
    
    // close the stream
    outStream.close();
}

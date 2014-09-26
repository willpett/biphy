#include "CladeReader.h"


using namespace RevBayesCore;


/**
 * Constructor. Here we read in immidiately the file and the we parse through each line 
 * and extract the taxon information.
 *
 * \param[in]     fn       The name of the file where the data is stored.
 * \param[in]     delim    The delimiter between the columns.
 */
CladeReader::CladeReader(const std::string &fn, char delim) : DelimitedDataReader( fn, delim )
{
    
    for (size_t i = 0; i < chars.size(); ++i) 
    {
        const std::vector<std::string>& line = chars[i];
        Clade c = Clade(line,0);
        
        clades.push_back( c );
    }
    
}


/**
 * Get the taxon information read from the file.
 *
 * \return The vector of taxa.
 */
const std::vector<Clade>& CladeReader::getClades( void ) const
{
    
    return clades;
}

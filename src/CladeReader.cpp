#include "CladeReader.h"


CladeReader::CladeReader(const std::string &fn, char delim) : DelimitedDataReader( fn, delim )
{
    
    for (size_t i = 0; i < chars.size(); ++i) 
    {
        const std::vector<std::string>& line = chars[i];
        Clade c = Clade(line);
        
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

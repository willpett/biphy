#ifndef CladeReader_H
#define CladeReader_H

#include "DelimitedDataReader.h"
#include "Clade.h"

namespace RevBayesCore {
    
    
    
    /**
     * Reader for taxon names mapped to sampling dates.
     *
     * This reader is a simple file reader of a delimited file, e.g., by tab-stops.
     * In the first column should be the taxon name and in the second column the date of sampling.
     *
     *
     * @copyright Copyright 2009-
     * @author The RevBayes Development Core Team (Sebastian Hoehna)
     * @since 2014-03-18, version 1.0
     *
     */
    class CladeReader : public DelimitedDataReader {
        
    public:
        
    	CladeReader(const std::string &fn, char d='\t');                    //!< Constructor
        
        const std::vector<Clade>&   getClades(void) const;                        //!< Get the taxa.

        
    protected:
        
        std::vector<Clade>          clades;
    };
}

#endif

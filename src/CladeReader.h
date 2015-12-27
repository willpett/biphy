#ifndef CladeReader_H
#define CladeReader_H

#include "DelimitedDataReader.h"
#include "Clade.h"

class CladeReader : public DelimitedDataReader {
    
public:
    
    CladeReader(const std::string &fn, char d='\t');                    //!< Constructor
    
    const std::vector<Clade>&   getClades(void) const;                        //!< Get the taxa.

    
protected:
    
    std::vector<Clade>          clades;
};

#endif

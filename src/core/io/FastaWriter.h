#ifndef FastaWriter_H
#define FastaWriter_H

#include "BinaryCharacterData.h"

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>

class FastaWriter {
    
public:
    FastaWriter();
    
    void                    writeData(const std::string& fn, const BinaryCharacterData &d, bool append = false);
    
    
};


#endif

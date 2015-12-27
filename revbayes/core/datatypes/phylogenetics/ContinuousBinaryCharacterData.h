#ifndef ContinuousBinaryCharacterData_H
#define ContinuousBinaryCharacterData_H

#include "BinaryCharacterData.h"

#include <map>
#include <set>
#include <string>
#include <vector>

class ContinuousBinaryCharacterData : public BinaryCharacterData {
    
public:
    
    std::string                         getDatatype(void) const { return "Continuous"; }
    
    ContinuousBinaryCharacterData*      clone(void) const;
    
};

#endif


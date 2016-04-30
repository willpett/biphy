#include "ContinuousBinaryCharacterData.h"

ContinuousBinaryCharacterData* ContinuousBinaryCharacterData::clone( void ) const 
{
    
    return new ContinuousBinaryCharacterData(*this);
}


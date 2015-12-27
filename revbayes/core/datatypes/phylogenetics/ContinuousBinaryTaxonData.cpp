#include "ContinuousBinaryTaxonData.h"


ContinuousBinaryTaxonData::ContinuousBinaryTaxonData(void) : BinaryTaxonData(),
        sequence() 
{
    
}


/**
 * Constructor with taxon name.
 * Does nothing except instanciating the object.
 */

ContinuousBinaryTaxonData::ContinuousBinaryTaxonData(const std::string &tname) : BinaryTaxonData(tname),
        sequence() 
{
    
}

ContinuousBinaryTaxonData* ContinuousBinaryTaxonData::clone( void ) const
{

    return new ContinuousBinaryTaxonData(*this);
}



/**
 * Subscript const operator for convenience access.
 *
 * \param[in]    i    The position of the character.
 *
 * \return            A non-const reference to the character
 */

RealNumber ContinuousBinaryTaxonData::operator[](size_t i) const
{
    
    if (i >= sequence.size()){
        throw Exception("Index out of bounds");
    }
    
    return sequence[i];
}


/**
 * Push back a new character.
 * 
 * \param[in]    newChar    The new character.
 */

void ContinuousBinaryTaxonData::addCharacter( RealNumber newChar ) 
{
    
    sequence.push_back( newChar );
    gaps.push_back( false );
    
}


/**
 * Get-operator for convenience access.
 *
 * \param[in]    i    The position of the character.
 *
 * \return            A non-const reference to the character
 */

RealNumber ContinuousBinaryTaxonData::getCharacter(size_t index) const 
{
    
    if (index >= sequence.size()){
        throw Exception("Index out of bounds");
    }
    
    return sequence[index];
}

size_t ContinuousBinaryTaxonData::getNumberOfCharacters(void) const 
{
    
    return sequence.size();
}


void ContinuousBinaryTaxonData::setCharacter(size_t index, RealNumber val) 
{
    
    if (index >= sequence.size()){
        throw Exception("Index out of bounds");
    }
    
    sequence[index] = val;
}


RealNumber ContinuousBinaryTaxonData::getElement(size_t i) const 
{
    
    if (i >= sequence.size()){
        throw Exception("Index out of bounds");
    }
    
    return sequence[i];
}

RealNumber ContinuousBinaryTaxonData::computeStateFrequency() const 
{
    
    RealNumber p = 0;
    size_t len = 0;
    
    for(size_t i = 0; i < gaps.size(); i++)
    {
        if(!gaps[i])
        {
            p += sequence[i];
            len++;
        }
    }
    
    return p/len;
}

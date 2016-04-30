#include "DiscreteBinaryTaxonData.h"

DiscreteBinaryTaxonData::DiscreteBinaryTaxonData(void) : 
    BinaryTaxonData(), 
    sequence()
{
    
}

DiscreteBinaryTaxonData::~DiscreteBinaryTaxonData()
{

}



DiscreteBinaryTaxonData* DiscreteBinaryTaxonData::clone( void ) const
{

    return new DiscreteBinaryTaxonData(*this);
}


/**
 * Constructor with taxon name.
 * Does nothing except instanciating the object.
 */

DiscreteBinaryTaxonData::DiscreteBinaryTaxonData(const std::string &tname) : 
    BinaryTaxonData(tname), 
    sequence()
{
    
}


/**
 * Subscript const operator for convenience access.
 *
 * \param[in]    i    The position of the character.
 *
 * \return            A non-const reference to the character
 */

RealNumber DiscreteBinaryTaxonData::operator[](size_t i) const
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

void DiscreteBinaryTaxonData::addCharacter( RealNumber newChar ) 
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

RealNumber DiscreteBinaryTaxonData::getCharacter(size_t index) const 
{
    
    if (index >= sequence.size()){
        throw Exception("Index out of bounds");
    }
    
    return sequence[index];
}

void DiscreteBinaryTaxonData::setCharacter(size_t index, RealNumber val) 
{
    
    if (index >= sequence.size()){
        throw Exception("Index out of bounds");
    }
    
    sequence[index] = val;
}


RealNumber DiscreteBinaryTaxonData::getElement(size_t i) const 
{
    
    if (i >= sequence.size()){
        throw Exception("Index out of bounds");
    }
    
    return sequence[i];
}


/**
 * Get the number of character stored in this object
 *
 * \return            The number of characters
 */

size_t DiscreteBinaryTaxonData::getNumberOfCharacters(void) const 
{
    
    return sequence.size();
}

RealNumber DiscreteBinaryTaxonData::computeStateFrequency() const 
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

    //std::cerr << p << "\t" << len << "\t" << p/RealNumber(len) << std::endl;

    return p/RealNumber(len);
}

#include "BinaryTaxonData.h"

BinaryTaxonData::BinaryTaxonData(void) : 
    taxonName("")
{
    
}

BinaryTaxonData::~BinaryTaxonData()
{

}


/**
 * Constructor with taxon name.
 * Does nothing except instanciating the object.
 */

BinaryTaxonData::BinaryTaxonData(const std::string &tname) : 
    taxonName(tname)
{
    
}

bool BinaryTaxonData::getGap(size_t index) const 
{
    
    if (index >= gaps.size()){
        throw Exception("Index out of bounds");
    }
    
    return gaps[index];
}

void BinaryTaxonData::setGap(size_t index, bool val) 
{
    
    if (index >= gaps.size()){
        throw Exception("Index out of bounds");
    }
    
    gaps[index] = val;
}


/**
 * Get the name of the taxon.
 *
 * \return            The taxon's name.
 */

const std::string& BinaryTaxonData::getTaxonName(void) const 
{
    
    return taxonName;
}


/**
 * Set the name of the taxon.
 *
 * \param[in]    tn    The new name of the taxon.
 */

void BinaryTaxonData::setTaxonName(std::string tn) 
{
    
    taxonName = tn;
}


/**
 * Get the size of the taxon which is the same as the number of characters.
 *
 * \return            The number of characters.
 */

size_t BinaryTaxonData::size(void) const 
{
    
    return gaps.size();
}

size_t BinaryTaxonData::getNonGapLength() const 
{
    
    size_t len = 0;
    
    for(size_t i = 0; i < gaps.size(); i++)
    {
        len += !gaps[i];
    }
    
    return len;
}

#include "BinaryCharacterData.h"


BinaryCharacterData::BinaryCharacterData()
{

}

BinaryCharacterData::~BinaryCharacterData()
{
	clear();
}

BinaryCharacterData::BinaryCharacterData(const BinaryCharacterData& m)
{
	sequenceNames     = m.sequenceNames;
	deletedTaxa       = m.deletedTaxa;
	deletedCharacters = m.deletedCharacters;
	fileName	      = m.fileName;

	for(size_t i=0; i < m.taxonData.size(); i++)
	{
		taxonData.push_back(m.taxonData[i]->clone());
	}
}


/** 
 * Index (const) operator to access a TaxonData object at position i.
 *
 * \param[in]    i    The position of the TaxonData object.
 *
 * \return            The TaxonData object at position i.
 */

BinaryTaxonData* BinaryCharacterData::operator[]( const size_t i ) 
{
    
    return getTaxonData( i );
}



/**
 * Less than operator.
 */

bool BinaryCharacterData::operator<(const BinaryCharacterData &x) const 
{
    
    return sequenceNames.size() < x.sequenceNames.size();
}


/** 
 * Add a sequence (TaxonData) to the character data object.
 *
 * \param[in]    obsd    The TaxonData object that should be added.
 */

void BinaryCharacterData::addTaxonData(BinaryTaxonData* obs) 
{
    
    // add the sequence name to the list
    sequenceNames.push_back( obs->getTaxonName() );
    
    // add the sequence also as a member so that we can access it by name
    taxonData.push_back(obs);
    
}



/** 
 * Clear the object, that is, remove all TaxonData elements.
 */

void BinaryCharacterData::clear( void ) 
{
    
    sequenceNames.clear();
    for(size_t i = 0; i < taxonData.size(); i++)
        delete taxonData[i];
    
    taxonData.clear();
    
}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'B'.
 *
 * \return A new copy of the object. 
 */

BinaryCharacterData* BinaryCharacterData::clone( void ) const 
{
    
    return new BinaryCharacterData(*this);
}



/**
 * Compute the state frequencies per site.
 *
 * \return       A matrix of site frequencies where each column is a site a each row the frequency of a character.
 */

std::vector<double> BinaryCharacterData::computeStateFrequencies( void ) const 
{
    
    bool tmp;
    size_t numSequences = sequenceNames.size();
    std::vector<double> m;
    for (size_t i = 0; i < numSequences; ++i) 
    {
        if(isTaxonExcluded(i))
            continue;
        
        const BinaryTaxonData* seq = getTaxonData(i);
        double pi = seq->computeStateFrequency();
        
        m.push_back(pi);
    }
    
    return m;
}

std::vector<int> BinaryCharacterData::computeCountDistribution( void ) const
{

    bool tmp;
    size_t numSequences = sequenceNames.size();
    std::vector<int> m = std::vector<int>(numSequences + 1, 0);
    for(size_t c = 0; c < this->getNumberOfCharacters(); c++)
    {
		if(isCharacterExcluded(c))
			continue;

		int count = 0;
		for (size_t i = 0; i < numSequences; ++i)
		{
			if(isTaxonExcluded(i))
				continue;

			const BinaryTaxonData* seq = getTaxonData(i);
			count += ((*seq)[c] > 0.0);
		}

		m[count]++;
    }

    return m;
}


/** 
 * Exclude a character.
 * We don't actually delete the character but mark it for exclusion.
 *
 * \param[in]    i    The position of the character that will be excluded.
 */

void BinaryCharacterData::excludeCharacter(size_t i) 
{
    
    if (i >= getTaxonData( 0 )->size() ) 
    {
        std::stringstream o;
        o << "Only " << getNumberOfCharacters() << " characters in matrix";
        throw Exception( o.str() );
    }
    
    
    deletedCharacters.insert( i );
    
}


/** 
 * Exclude a taxon.
 * We don't actually delete the taxon but instead mark it for exclusion.
 *
 * \param[in]    i    The index of the taxon that will be excluded.
 */

void BinaryCharacterData::excludeTaxon(size_t i) 
{

    if (i >= taxonData.size()) 
    {
        std::stringstream o;
        o << "Only " << taxonData.size() << " taxa in matrix";
        throw Exception( o.str() );
    }
    
    deletedTaxa.insert( i );
}


/** 
 * Exclude a taxon.
 * We don't actually delete the taxon but instead mark it for exclusion.
 *
 * \param[in]    s    The name of the taxon that will be excluded.
 */

void BinaryCharacterData::excludeTaxon(std::string& s) 
{
    
    for (size_t i = 0; i < getNumberOfTaxa(); i++) 
    {
        if (s == sequenceNames[i] ) 
        {
            deletedTaxa.insert( i );
            break;
        }
    }
    
}

bool BinaryCharacterData::getGap( size_t tn, size_t cn ) const 
{
    
    if ( cn >= getNumberOfCharacters() ){
        throw Exception( "Character index out of range" );
    }
    
    return getTaxonData( tn )->getGap(cn);
}


/**
 * Get the data type of the character stored in this object.
 *
 * \return      The data type (e.g. DNA, RNA or Standard).
 */

std::string BinaryCharacterData::getDatatype(void) const 
{
    return "Binary";
}


/**
 * Get the file name from whcih the character data object was read in.
 *
 * \return    The original file name.
 */

const std::string& BinaryCharacterData::getFileName(void) const 
{
    
    return fileName;
}


/** 
 * Get the number of characters in taxon data object. 
 * This i regardless of whether the character are included or excluded.
 * For simplicity we assume that all taxon data objects contain the same number
 * of character and thus we simply return the number from the first taxon data object.
 *
 * \return    The total number of characters
 */

size_t BinaryCharacterData::getNumberOfCharacters(void) const 
{
    
    if (getNumberOfTaxa() > 0) 
    {
        return getTaxonData(0)->getNumberOfCharacters();
    }
    
    return 0;
}


/** 
 * Get the number of characters in the i-th taxon data object. 
 * This i regardless of whether the character are included or excluded.
 *
 * \param[in]    i     The index of the taxon data object.
 *
 * \return             The total number of characters
 */

size_t BinaryCharacterData::getNumberOfCharacters(size_t idx) const {
    
    if (getNumberOfTaxa() > 0) 
    {
        return getTaxonData(idx)->getNumberOfCharacters();
    }
    
    return 0;
}


/** 
 * Get the number of characters in taxon data object. 
 * This i regardless of whether the character are included or excluded.
 * For simplicity we assume that all taxon data objects contain the same number
 * of character and thus we simply return the number from the first taxon data object.
 *
 * \return    The total number of characters
 */

size_t BinaryCharacterData::getNumberOfIncludedCharacters(void) const {
    
    if (getNumberOfTaxa() > 0) 
    {
        return getTaxonData(0)->getNumberOfCharacters() - deletedCharacters.size();
    }
    return 0;
}


/** 
 * Get the number of included characters in the i-th taxon data object.
 *
 * \param[in]    i     The index of the taxon data object.
 *
 * \return             The total number of characters
 */

size_t BinaryCharacterData::getNumberOfIncludedCharacters(size_t idx) const {
    
    if (getNumberOfTaxa() > 0) 
    {
        return getTaxonData(idx)->getNumberOfCharacters() - deletedCharacters.size();
    }
    
    return 0;
}


/**
 * Get the number of taxa currently stored in this object.
 *
 * \return       The number of taxa.
 */

size_t BinaryCharacterData::getNumberOfTaxa(void) const 
{
    
    return sequenceNames.size();
}


/** 
 * Get the taxon data object with index tn.
 *
 * \return     A const reference to the taxon data object at position tn.
 */

const BinaryTaxonData* BinaryCharacterData::getTaxonData( size_t tn ) const 
{
    
    if ( tn >= getNumberOfTaxa() ){
        throw Exception( "Taxon index out of range" );
    }
    
    return taxonData[tn];
    
}


/** 
 * Get the taxon data object at position tn.
 *
 * \return     A non-const reference to the taxon data object at position tn.
 */

BinaryTaxonData* BinaryCharacterData::getTaxonData( size_t tn ) 
{
    
    if ( tn >= getNumberOfTaxa() ){
        throw Exception( "Taxon index out of range" );
    }
    
    return taxonData[tn];
    
}


/** 
 * Get the taxon data object with name tn.
 *
 * \return     A non-const reference to the taxon data object with name tn.
 */

const BinaryTaxonData* BinaryCharacterData::getTaxonData( const std::string &tn ) const 
{
    
    if ( tn == "" )
        throw Exception("Ambiguous taxon name.");
    
    for(size_t i = 0; i < sequenceNames.size(); i++)
    {
        if(sequenceNames[i] == tn)
            return taxonData[i];
    }
    
    throw Exception("Cannot find taxon '" + tn + "' in the CharacterData matrix.");
}


/** 
 * Get the taxon data object with name tn.
 *
 * \return     A const reference to the taxon data object with name tn.
 */

BinaryTaxonData* BinaryCharacterData::getTaxonData( const std::string &tn ) 
{
    if ( tn == "" )
        throw Exception("Ambiguous taxon name.");
    
    for(size_t i = 0; i < sequenceNames.size(); i++)
    {
        if(sequenceNames[i] == tn)
            return taxonData[i];
    }
    
    throw Exception("Cannot find taxon '" + tn + "' in the CharacterData matrix.");
    
}


/**
 * Get the names of all taxa.
 *
 * \return     A vector of all taxon names.
 */

const std::vector<std::string>& BinaryCharacterData::getTaxonNames( void ) const 
{
    
    return sequenceNames;
}



/** 
 * Get the taxon name with index idx.
 *
 * \param[in]    idx    The position of the taxon.
 *
 * \return              The name of the taxon.
 */

const std::string& BinaryCharacterData::getTaxonNameWithIndex( size_t idx ) const 
{
    
    return sequenceNames[idx];
}



/** 
 * Get the index of the taxon with name s.
 *
 * \param[in]    s    The name of the taxon.
 *
 * \return            The index of the taxon.
 */

size_t BinaryCharacterData::indexOfTaxonWithName( std::string& s ) const 
{
    
    // search through all names
    for (size_t i=0; i<sequenceNames.size(); i++) 
    {
        if (s == sequenceNames[i] ) 
        {
            return i;
        }
    }
    
    return -1;
}


/** 
 * Is this character pattern constant at site idx?
 * 
 * \param[in]   idx    The site at which we want to know if it is constant?
 */

bool BinaryCharacterData::isCharacterConstant(size_t idx) const 
{
    
    RealNumber f;
    bool first = true;
    for ( size_t i=0; i<getNumberOfTaxa(); i++ ) 
    {
        if ( isTaxonExcluded(i) == false ) 
        {
            if ( first )
            {
                f = taxonData[i]->getCharacter(idx );
                first = false;
            }
            else 
            {
                RealNumber s = taxonData[i]->getCharacter(idx );
                if ( f != s )
                    return false;
            }
        }
    }
    
    return true;
}


/** 
 * Is the character excluded?
 *
 * \param[in]    i   The position of the character.
 */

bool BinaryCharacterData::isCharacterExcluded(size_t i) const 
{
    
    std::set<size_t>::const_iterator it = deletedCharacters.find( i );
    if ( it != deletedCharacters.end() )
        return true;
    
    return false;
}


/** 
 * Does the character have missing or ambiguous characters?
 *
 * \param[in]    idx    The position of the character in question.
 */

bool BinaryCharacterData::isCharacterMissingOrAmbiguous(size_t idx) const 
{
    
    for ( size_t i=0; i<getNumberOfTaxa(); i++ )
    {
        if ( isTaxonExcluded(i) == false )
        {
            if ( getGap( i, idx ) )
                return true;
        }
    }
    
    return false;
}



/**
 * Is the taxon excluded.
 *
 * \param[in]    idx    The position of the taxon in question.
 */

bool BinaryCharacterData::isTaxonExcluded(size_t i) const 
{
    
    std::set<size_t>::const_iterator it = deletedTaxa.find( i );
    if ( it != deletedTaxa.end() )
        return true;
    
    return false;
}


/** 
 * Is the taxon excluded?
 *
 * \param[in]    s    The name of the taxon in question.
 */

bool BinaryCharacterData::isTaxonExcluded(std::string& s) const 
{
    
    size_t i = indexOfTaxonWithName(s);
    std::set<size_t>::const_iterator it = deletedTaxa.find( i );
    if ( it != deletedTaxa.end() )
        return true;
    
    return false;
}


/** 
 * Calculates and returns the number of constant characters.
 */

size_t BinaryCharacterData::numConstantPatterns( void ) const 
{
    
    size_t nc = 0;
    for (size_t i=0; i<getNumberOfCharacters(); i++)
    {
        if ( isCharacterExcluded(i) == false && isCharacterConstant(i) == true )
            nc++;
    }
    
    return nc;
}


/** 
 * Returns the number of characters with missing or ambiguous data
 */

size_t BinaryCharacterData::numMissAmbig(void) const 
{
    
    size_t nma = 0;
    for (size_t i=0; i<getNumberOfCharacters(); i++)
    {
        if ( isCharacterExcluded(i) == false && isCharacterMissingOrAmbiguous(i) == true )
            nma++;
    }
    
    return nma;
}


/** 
 * Restore a character. We simply do not mark the character as excluded anymore.
 *
 * \param[in]    i    The position of the character to restore.
 */

void BinaryCharacterData::restoreCharacter(size_t i) 
{
    
    if (i >= getNumberOfCharacters() ){
        throw Exception( "Character index out of range" );
    }
    deletedCharacters.erase( i );
    
}


/** 
 * Restore a taxon. We simply do not mark the taxon as excluded anymore
 *
 *
 * \param[in]    i    The position of the taxon in question.
 */

void BinaryCharacterData::restoreTaxon(size_t i) 
{
    
    if ( i >= getNumberOfTaxa() )
        return;
    
    deletedTaxa.erase( i );
    
}


/** 
 * Restore a taxon. We simply do not mark the taxon as excluded anymore.
 *
 * \param[in]    s    The name of the taxon in question.
 */

void BinaryCharacterData::restoreTaxon(std::string& s) 
{
    
    size_t i = indexOfTaxonWithName( s );
    
    deletedTaxa.erase( i );
    
}


/**
 * Set the original file name for this character data object.
 *
 * \param[in]    fn    The new file name.
 */

void BinaryCharacterData::setFileName(const std::string& fn) 
{
    
    fileName = fn;
    
}


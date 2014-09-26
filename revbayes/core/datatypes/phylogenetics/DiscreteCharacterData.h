#ifndef DiscreteCharacterData_H
#define DiscreteCharacterData_H

#include "AbstractDiscreteCharacterData.h"
#include "DiscreteCharacterState.h"
#include "MatrixReal.h"
#include "DiscreteTaxonData.h"

#include <map>
#include <set>
#include <string>
#include <vector>

namespace RevBayesCore {
    
    /**
     * Template class for all discrete character data objects.
     *
     * The discrete character data class is template class for all discrete character data objects.
     * Here the data is stored via a set of DiscreteTaxonData objects.
     * Several convenience functions are provided to access, store and delete the data.
     *
     * @copyright Copyright 2009-
     * @author The RevBayes Development Core Team (Sebastian Hoehna)
     * @since 2013-04-15, version 1.0
     */
    template<class charType>
    class DiscreteCharacterData : public AbstractDiscreteCharacterData {
        
    public:
        DiscreteCharacterData();                                                                                        //!< Default constructor
        
        // Overloaded operators
        const DiscreteTaxonData<charType>&  operator[](size_t i) const;                                                 //!< Subscript operator (const)
        bool                                operator<(const DiscreteCharacterData& x) const;                            //!< Less than operator
        
        // implemented methods of the Cloneable interface
        DiscreteCharacterData<charType>*    clone(void) const;
        
        // Container functions
        void                                clear(void);
        
        // CharacterData functions
        void                                addTaxonData(const AbstractTaxonData &obs);                                 //!< Add taxon data
        void                                addTaxonData(const AbstractDiscreteTaxonData &obs);                         //!< Add discrete taxon data
        void                                addTaxonData(const DiscreteTaxonData<charType> &obs);                               //!< Add taxon data
        MatrixReal                          computeStateFrequencies(void) const;
        void                                excludeCharacter(size_t i);                                                 //!< Exclude character
        void                                excludeTaxon(size_t i);                                                     //!< Exclude taxon
        void                                excludeTaxon(std::string& s);                                               //!< Exclude taxon
        const charType&                     getCharacter(size_t tn, size_t cn) const;                                   //!< Return a reference to a character element in the character matrix
        std::string                         getDatatype(void) const;
        const std::string&                  getFileName(void) const;                                                    //!< Returns the name of the file the data came from
        size_t                              getNumberOfCharacters(void) const;                                          //!< Number of characters
        size_t                              getNumberOfCharacters(size_t idx) const;                                    //!< Number of characters for a specific taxon
        size_t                              getNumberOfIncludedCharacters(void) const;                                          //!< Number of characters
        size_t                              getNumberOfIncludedCharacters(size_t idx) const;                                    //!< Number of characters for a specific taxon
        size_t                              getNumberOfStates(void) const;                                              //!< Get the number of states for the characters in this matrix
        size_t                              getNumberOfTaxa(void) const;                                                //!< Number of taxa
        DiscreteTaxonData<charType>&        getTaxonData(size_t tn);                                                    //!< Return a reference to a sequence in the character matrix
        const DiscreteTaxonData<charType>&  getTaxonData(size_t tn) const;                                              //!< Return a reference to a sequence in the character matrix
        DiscreteTaxonData<charType>&        getTaxonData(const std::string &tn);                                        //!< Return a reference to a sequence in the character matrix
        const DiscreteTaxonData<charType>&  getTaxonData(const std::string &tn) const;                                  //!< Return a reference to a sequence in the character matrix
        const std::vector<std::string>&     getTaxonNames(void) const;                                                  //!< Get the names of the taxa
        const std::string&                  getTaxonNameWithIndex(size_t idx) const;                                    //!< Returns the idx-th taxon name
        bool                                isCharacterExcluded(size_t i) const;                                        //!< Is the character excluded
        bool                                isHomologyEstablished(void) const;                                          //!< Returns whether the homology of the characters has been established
        bool                                isTaxonExcluded(size_t i) const;                                            //!< Is the taxon excluded
        bool                                isTaxonExcluded(std::string& s) const;                                      //!< Is the taxon excluded
        void                                restoreCharacter(size_t i);                                                 //!< Restore character
        void                                restoreTaxon(size_t i);                                                     //!< Restore taxon
        void                                restoreTaxon(std::string& s);                                               //!< Restore taxon
        void                                setFileName(const std::string &fn);                                         //!< Set the file name
        void                                setHomologyEstablished(bool tf);                                            //!< Set whether the homology of the characters has been established
        
    protected:
        // Utility functions
        size_t                              indexOfTaxonWithName(std::string& s) const;                                 //!< Get the index of the taxon
        bool                                isCharacterConstant(size_t idx) const;                                      //!< Is the idx-th character a constant pattern?
        bool                                isCharacterMissingOrAmbiguous(size_t idx) const;                            //!< Does the character have missing or ambiguous data?
        size_t                              numConstantPatterns(void) const;                                            //!< The number of constant patterns
        size_t                              numMissAmbig(void) const;                                                   //!< The number of patterns with missing or ambiguous characters
        
        // Member variables
        std::set<size_t>                    deletedTaxa;                                                                //!< Set of deleted taxa
        std::set<size_t>                    deletedCharacters;                                                          //!< Set of deleted characters
        std::string                         fileName;                                                                   //!< The path/filename from where this matrix originated
        std::vector<std::string>            sequenceNames;                                                              //!< names of the sequences
        size_t                              sequenceLength;                                                             //!< The length of each sequence
        bool                                homologyEstablished;                                                        //!< Whether the homology of the characters has been established
        
        std::map<std::string, DiscreteTaxonData<charType> >    taxonMap;
        
    };
    
    
}


#include "DiscreteCharacterState.h"
#include "DiscreteTaxonData.h"
#include "RbException.h"

#include <cmath>
#include <fstream>
#include <sstream>
#include <string>


/**
 * Default constructor,
 * Does nothing except instanciating the object.
 */
template<class charType>
RevBayesCore::DiscreteCharacterData<charType>::DiscreteCharacterData() 
{

}



/** 
 * Index (const) operator to access a TaxonData object at position i.
 *
 * \param[in]    i    The position of the TaxonData object.
 *
 * \return            The TaxonData object at position i.
 */
template<class charType>
const RevBayesCore::DiscreteTaxonData<charType>& RevBayesCore::DiscreteCharacterData<charType>::operator[]( const size_t i ) const 
{
    
    return getTaxonData( i );
}



/**
 * Less than operator.
 */
template<class charType>
bool RevBayesCore::DiscreteCharacterData<charType>::operator<(const DiscreteCharacterData<charType> &x) const 
{
    
    return sequenceNames.size() < x.sequenceNames.size();
}


/** 
 * Add a sequence (TaxonData) to the character data object.
 *
 * \param[in]    obsd    The TaxonData object that should be added.
 */
template<class charType>
void RevBayesCore::DiscreteCharacterData<charType>::addTaxonData(const AbstractTaxonData &obsd) 
{
    
#ifdef ASSERTIONS_ALL
    if ( dynamic_cast<const DiscreteTaxonData<charType>* >( &obsd ) == NULL ) 
    {
        throw RbException("Inserting wrong character type into CharacterData!!!");
    }
#endif
    
    // delegate the call to the specialized method
    addTaxonData( static_cast<const DiscreteTaxonData<charType>& >( obsd ) );
    
}


/** 
 * Add a sequence (TaxonData) to the character data object.
 *
 * \param[in]    obsd    The TaxonData object that should be added.
 */
template<class charType>
void RevBayesCore::DiscreteCharacterData<charType>::addTaxonData(const AbstractDiscreteTaxonData &obsd) 
{
    
#ifdef ASSERTIONS_ALL
    if ( dynamic_cast<const DiscreteTaxonData<charType>* >( &obsd ) == NULL ) 
    {
        throw RbException("Inserting wrong character type into CharacterData!!!");
    }
#endif
    
    // delegate the call to the specialized method
    addTaxonData( static_cast<const DiscreteTaxonData<charType>& >( obsd ) );
    
}


/** 
 * Add a sequence (TaxonData) to the character data object.
 *
 * \param[in]    obsd    The TaxonData object that should be added.
 */
template<class charType>
void RevBayesCore::DiscreteCharacterData<charType>::addTaxonData(const DiscreteTaxonData<charType> &obs) 
{
    
    // add the sequence name to the list
    sequenceNames.push_back( obs.getTaxonName() );
    
    // add the sequence also as a member so that we can access it by name
    taxonMap.insert( std::pair<std::string, DiscreteTaxonData<charType> >( obs.getTaxonName(), obs ) );
    
}



/** 
 * Clear the object, that is, remove all TaxonData elements.
 */
template<class charType>
void RevBayesCore::DiscreteCharacterData<charType>::clear( void ) 
{
    
    sequenceNames.clear();
    taxonMap.clear();
    
}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'B'.
 *
 * \return A new copy of the object. 
 */
template<class charType>
RevBayesCore::DiscreteCharacterData<charType>* RevBayesCore::DiscreteCharacterData<charType>::clone( void ) const 
{
    
    return new DiscreteCharacterData<charType>(*this);
}



/**
 * Compute the state frequencies per site.
 *
 * \return       A matrix of site frequencies where each column is a site a each row the frequency of a character.
 */
template<class charType>
RevBayesCore::MatrixReal RevBayesCore::DiscreteCharacterData<charType>::computeStateFrequencies( void ) const 
{
    
    charType tmp;
    std::string labels = tmp.getStateLabels();
    size_t alphabetSize = tmp.getNumberOfStates();
    size_t numSequences = this->sequenceNames.size();
    MatrixReal m(numSequences,alphabetSize);
    for (size_t i = 0; i < numSequences; ++i) 
    {
        const DiscreteTaxonData<charType>& seq = this->getTaxonData(i);
        size_t l = seq.size();
        double nonGapSeqLength = 0.0;
        std::vector<double> stateCounts(alphabetSize);
        for (size_t j = 0; j < l; ++j) 
        {
            const charType& c = seq[j];
            
            if ( !c.isGapState() ) 
            {
                nonGapSeqLength++;
                
                unsigned long state = c.getState();
                double numObservedStates = c.getNumberObservedStates();
                
                size_t index = 0;
                do 
                {
                    
                    if ( (state & 1) == 1 ) {
                        // add a uniform probability of having observed each of the ambiguous characters
                        stateCounts[index] += 1.0 / numObservedStates;
                    }
                    state >>= 1;
                    ++index;
                } while ( state != 0 );
            }
        } // finished loop over sequence
        
        // set the observed state frequencies for this sequence into the matrix
        std::vector<double> &observedFreqs = m[i];
        for (size_t j = 0; j < alphabetSize; ++j) {
            observedFreqs[j] = stateCounts[j] / nonGapSeqLength;
        }
    }
    
    return m;
}


/** 
 * Exclude a character.
 * We don't actually delete the character but mark it for exclusion.
 *
 * \param[in]    i    The position of the character that will be excluded.
 */
template<class charType>
void RevBayesCore::DiscreteCharacterData<charType>::excludeCharacter(size_t i) 
{
    
    if (i >= getTaxonData( 0 ).size() ) 
    {
        std::stringstream o;
        o << "Only " << getNumberOfCharacters() << " characters in matrix";
        throw RbException( o.str() );
    }
    
    
    deletedCharacters.insert( i );
    
}


/** 
 * Exclude a taxon.
 * We don't actually delete the taxon but instead mark it for exclusion.
 *
 * \param[in]    i    The index of the taxon that will be excluded.
 */
template<class charType>
void RevBayesCore::DiscreteCharacterData<charType>::excludeTaxon(size_t i) 
{
    
    if (i >= taxonMap.size()) 
    {
        std::stringstream o;
        o << "Only " << taxonMap.size() << " taxa in matrix";
        throw RbException( o.str() );
    }
    
    deletedTaxa.insert( i );
}


/** 
 * Exclude a taxon.
 * We don't actually delete the taxon but instead mark it for exclusion.
 *
 * \param[in]    s    The name of the taxon that will be excluded.
 */
template<class charType>
void RevBayesCore::DiscreteCharacterData<charType>::excludeTaxon(std::string& s) 
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



/** 
 * Get the cn-th character of the tn-th taxon.
 *
 * \param[in]    tn     The index/position of the taxon.
 * \param[in]    cn     The position of the character.
 *
 * \return              The cn-th character of the tn-th taxon. 
 */
template<class charType>
const charType& RevBayesCore::DiscreteCharacterData<charType>::getCharacter( size_t tn, size_t cn ) const 
{
    
    if ( cn >= getNumberOfCharacters() )
        throw RbException( "Character index out of range" );
    
    return getTaxonData( tn )[cn];
}


/**
 * Get the data type of the character stored in this object.
 *
 * \return      The data type (e.g. DNA, RNA or Standard).
 */
template<class charType>
std::string RevBayesCore::DiscreteCharacterData<charType>::getDatatype(void) const 
{
    
    std::string dt = "";
    if ( sequenceNames.size() > 0 ) 
    {
        const DiscreteTaxonData<charType> &t = getTaxonData( sequenceNames[0] );
        if ( t.size() > 0 ) 
        {
            dt = t[0].getDatatype();
        }
        
    }
    
    return dt;
}


/**
 * Get the file name from whcih the character data object was read in.
 *
 * \return    The original file name.
 */
template<class charType>
const std::string& RevBayesCore::DiscreteCharacterData<charType>::getFileName(void) const 
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
template<class charType>
size_t RevBayesCore::DiscreteCharacterData<charType>::getNumberOfCharacters(void) const 
{
    
    if (getNumberOfTaxa() > 0) 
    {
        return getTaxonData(0).getNumberOfCharacters();
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
template<class charType>
size_t RevBayesCore::DiscreteCharacterData<charType>::getNumberOfCharacters(size_t idx) const {
    
    if (getNumberOfTaxa() > 0) 
    {
        return getTaxonData(idx).getNumberOfCharacters();
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
template<class charType>
size_t RevBayesCore::DiscreteCharacterData<charType>::getNumberOfIncludedCharacters(void) const {
    
    if (getNumberOfTaxa() > 0) 
    {
        return getTaxonData(0).getNumberOfCharacters() - deletedCharacters.size();
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
template<class charType>
size_t RevBayesCore::DiscreteCharacterData<charType>::getNumberOfIncludedCharacters(size_t idx) const {
    
    if (getNumberOfTaxa() > 0) 
    {
        return getTaxonData(idx).getNumberOfCharacters() - deletedCharacters.size();
    }
    
    return 0;
}


/** 
 * Get the number of states for the characters in this object. 
 * We assume that all of the characters in the matrix are of the same
 * type and have the same number of potential states. 
 *
 * \return      The number of states for the characters.
 */
template<class charType>
size_t RevBayesCore::DiscreteCharacterData<charType>::getNumberOfStates(void) const 
{
    
    // Get the first character in the matrix
    if ( getNumberOfTaxa() == 0 )
        return 0;
    
    const DiscreteTaxonData<charType>& sequence = getTaxonData( 0 );
    if ( sequence.getNumberOfCharacters() == 0 )
        return 0;
    
    return sequence[0].getNumberOfStates();
}


/**
 * Get the number of taxa currently stored in this object.
 *
 * \return       The number of taxa.
 */
template<class charType>
size_t RevBayesCore::DiscreteCharacterData<charType>::getNumberOfTaxa(void) const 
{
    
    return sequenceNames.size();
}


/** 
 * Get the taxon data object with index tn.
 *
 * \return     A const reference to the taxon data object at position tn.
 */
template<class charType>
const RevBayesCore::DiscreteTaxonData<charType>& RevBayesCore::DiscreteCharacterData<charType>::getTaxonData( size_t tn ) const 
{
    
    if ( tn >= getNumberOfTaxa() )
        throw RbException( "Taxon index out of range" );
    
    const std::string& name = sequenceNames[tn];
    const typename std::map<std::string, DiscreteTaxonData<charType> >::const_iterator& i = taxonMap.find( name ); 
    
    if (i != taxonMap.end() ) 
    {
        return i->second;
    }
    else 
    {
        throw RbException("Cannot find taxon '" + name + "' in the CharacterData matrix.");
    }
    
}


/** 
 * Get the taxon data object at position tn.
 *
 * \return     A non-const reference to the taxon data object at position tn.
 */
template<class charType>
RevBayesCore::DiscreteTaxonData<charType>& RevBayesCore::DiscreteCharacterData<charType>::getTaxonData( size_t tn ) 
{
    
    if ( tn >= getNumberOfTaxa() )
        throw RbException( "Taxon index out of range" );
    
    const std::string& name = sequenceNames[tn];
    const typename std::map<std::string, DiscreteTaxonData<charType> >::iterator& i = taxonMap.find( name ); 
    
    if (i != taxonMap.end() ) 
    {
        return i->second;
    }
    else 
    {
        throw RbException("Cannot find taxon '" + name + "' in the CharacterData matrix.");
    }
    
}


/** 
 * Get the taxon data object with name tn.
 *
 * \return     A non-const reference to the taxon data object with name tn.
 */
template<class charType>
const RevBayesCore::DiscreteTaxonData<charType>& RevBayesCore::DiscreteCharacterData<charType>::getTaxonData( const std::string &tn ) const 
{
    
    if ( tn == "" ) 
    {
        throw RbException("Ambiguous taxon name.");
    }
    
    const typename std::map<std::string, DiscreteTaxonData<charType> >::const_iterator& i = taxonMap.find(tn); 
    
    if (i != taxonMap.end() ) 
    {
        return i->second;
    }
    else 
    {
        throw RbException("Cannot find taxon '" + tn + "' in the CharacterData matrix.");
    }
    
}


/** 
 * Get the taxon data object with name tn.
 *
 * \return     A const reference to the taxon data object with name tn.
 */
template<class charType>
RevBayesCore::DiscreteTaxonData<charType>& RevBayesCore::DiscreteCharacterData<charType>::getTaxonData( const std::string &tn ) 
{
    
    
    if ( tn == "" ) 
    {
        throw RbException("Ambiguous taxon name.");
    }
    
    const typename std::map<std::string, DiscreteTaxonData<charType> >::iterator& i = taxonMap.find(tn); 
    
    if (i != taxonMap.end() ) 
    {
        return i->second;
    }
    else 
    {
        
        throw RbException("Cannot find taxon '" + tn + "' in the CharacterData matrix.");
    }
    
}


/**
 * Get the names of all taxa.
 *
 * \return     A vector of all taxon names.
 */
template<class charType>
const std::vector<std::string>& RevBayesCore::DiscreteCharacterData<charType>::getTaxonNames( void ) const 
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
template<class charType>
const std::string& RevBayesCore::DiscreteCharacterData<charType>::getTaxonNameWithIndex( size_t idx ) const 
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
template<class charType>
size_t RevBayesCore::DiscreteCharacterData<charType>::indexOfTaxonWithName( std::string& s ) const 
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
template<class charType>
bool RevBayesCore::DiscreteCharacterData<charType>::isCharacterConstant(size_t idx) const 
{
    
    const CharacterState* f = NULL;
    for ( size_t i=0; i<getNumberOfTaxa(); i++ ) 
    {
        if ( isTaxonExcluded(i) == false ) 
        {
            if ( f == NULL )
                f = &getCharacter( i, idx );
            else 
            {
                const CharacterState* s = &getCharacter( i , idx );
                if ( (*f) != (*s) )
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
template<class charType>
bool RevBayesCore::DiscreteCharacterData<charType>::isCharacterExcluded(size_t i) const 
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
template<class charType>
bool RevBayesCore::DiscreteCharacterData<charType>::isCharacterMissingOrAmbiguous(size_t idx) const 
{
    
    for ( size_t i=0; i<getNumberOfTaxa(); i++ )
    {
        if ( isTaxonExcluded(i) == false )
        {
            const CharacterState& c = getCharacter( i, idx );
            if ( c.isGapState() == true || c.isAmbiguous() == true )
                return true;
        }
    }
    
    return false;
}


/**
 * Is the homology established, i.e., is the character data object aligned?
 *
 * \return     True/False whether the homology was established.
 */
template<class charType>
bool RevBayesCore::DiscreteCharacterData<charType>::isHomologyEstablished(void) const 
{
    
    return homologyEstablished;
}


/**
 * Is the taxon excluded.
 *
 * \param[in]    idx    The position of the taxon in question.
 */
template<class charType>
bool RevBayesCore::DiscreteCharacterData<charType>::isTaxonExcluded(size_t i) const 
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
template<class charType>
bool RevBayesCore::DiscreteCharacterData<charType>::isTaxonExcluded(std::string& s) const 
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
template<class charType>
size_t RevBayesCore::DiscreteCharacterData<charType>::numConstantPatterns( void ) const 
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
template<class charType>
size_t RevBayesCore::DiscreteCharacterData<charType>::numMissAmbig(void) const 
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
template<class charType>
void RevBayesCore::DiscreteCharacterData<charType>::restoreCharacter(size_t i) 
{
    
    if (i >= getNumberOfCharacters() )
        throw RbException( "Character index out of range" );
    
    deletedCharacters.erase( i );
    
}


/** 
 * Restore a taxon. We simply do not mark the taxon as excluded anymore
 *
 *
 * \param[in]    i    The position of the taxon in question.
 */
template<class charType>
void RevBayesCore::DiscreteCharacterData<charType>::restoreTaxon(size_t i) 
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
template<class charType>
void RevBayesCore::DiscreteCharacterData<charType>::restoreTaxon(std::string& s) 
{
    
    size_t i = indexOfTaxonWithName( s );
    
    deletedTaxa.erase( i );
    
}


/**
 * Set the original file name for this character data object.
 *
 * \param[in]    fn    The new file name.
 */
template<class charType>
void RevBayesCore::DiscreteCharacterData<charType>::setFileName(const std::string& fn) 
{
    
    fileName = fn;
    
}


/**
 * Set whether the homology has been established.
 *
 * \param[in]    tf    Whether the homology has been established.
 */
template<class charType>
void RevBayesCore::DiscreteCharacterData<charType>::setHomologyEstablished(bool tf) 
{
    
    homologyEstablished = tf;
    
}



#endif


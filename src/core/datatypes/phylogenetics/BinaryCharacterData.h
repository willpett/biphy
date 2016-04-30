#ifndef BinaryCharacterData_H
#define BinaryCharacterData_H

#include "BinaryTaxonData.h"

#include <map>
#include <set>
#include <string>
#include <vector>

#include "numerics.h"

class BinaryCharacterData {
    
public:
    BinaryCharacterData();
    BinaryCharacterData(const BinaryCharacterData&);
    virtual ~BinaryCharacterData();
    
    // Overloaded operators
    BinaryTaxonData*                    operator[](size_t i);                                                 //!< Subscript operator (const)
    bool                                operator<(const BinaryCharacterData& x) const;                            //!< Less than operator
    
    // implemented methods of the Cloneable interface
    BinaryCharacterData*                clone(void) const;
    
    // Container functions
    void                                clear(void);
    
    // CharacterData functions
    void                                addTaxonData(BinaryTaxonData* obs);                                 //!< Add taxon data
    std::vector<double>                 computeStateFrequencies(void) const;
    std::vector<int>                 	computeCountDistribution(void) const;
    void                                excludeCharacter(size_t i);                                                 //!< Exclude character
    void                                excludeTaxon(size_t i);                                                     //!< Exclude taxon
    void                                excludeTaxon(std::string& s);
    bool                                getGap(size_t tn, size_t cn) const;
    virtual std::string                 getDatatype(void) const;
    const std::string&                  getFileName(void) const; 
    size_t                              getNumberOfCharacters(void) const;                                          //!< Number of characters
    size_t                              getNumberOfCharacters(size_t idx) const;                                    //!< Number of characters for a specific taxon
    size_t                              getNumberOfIncludedCharacters(void) const;                                          //!< Number of characters
    size_t                              getNumberOfIncludedCharacters(size_t idx) const;                                    //!< Number of characters for a specific taxon
    size_t                              getNumberOfTaxa(void) const;                                                //!< Number of taxa
    BinaryTaxonData*                    getTaxonData(size_t tn);                                                    //!< Return a reference to a sequence in the character matrix
    const BinaryTaxonData*              getTaxonData(size_t tn) const;                                              //!< Return a reference to a sequence in the character matrix
    BinaryTaxonData*                    getTaxonData(const std::string &tn);                                        //!< Return a reference to a sequence in the character matrix
    const BinaryTaxonData*              getTaxonData(const std::string &tn) const;                                  //!< Return a reference to a sequence in the character matrix
    const std::vector<std::string>&     getTaxonNames(void) const;                                                  //!< Get the names of the taxa
    const std::string&                  getTaxonNameWithIndex(size_t idx) const;                                    //!< Returns the idx-th taxon name
    bool                                isCharacterExcluded(size_t i) const;                                         //!< Returns whether the homology of the characters has been established
    bool                                isTaxonExcluded(size_t i) const;                                            //!< Is the taxon excluded
    bool                                isTaxonExcluded(std::string& s) const;                                      //!< Is the taxon excluded
    void                                restoreCharacter(size_t i);                                                 //!< Restore character
    void                                restoreTaxon(size_t i);                                                     //!< Restore taxon
    void                                restoreTaxon(std::string& s);                                               //!< Restore taxon
    void                                setFileName(const std::string &fn);                                          //!< Set whether the homology of the characters has been established
    
protected:
    // Utility functions
    size_t                              indexOfTaxonWithName(std::string& s) const;                                 //!< Get the index of the taxon
    bool                                isCharacterConstant(size_t idx) const;                                      //!< Is the idx-th character a constant pattern?
    bool                                isCharacterMissingOrAmbiguous(size_t idx) const;                            //!< Does the character have missing or ambiguous data?
    size_t                              numConstantPatterns(void) const;                                            //!< The number of constant patterns
    size_t                              numMissAmbig(void) const; 
            
    
    std::set<size_t>                    deletedTaxa;                                                                //!< Set of deleted taxa
    std::set<size_t>                    deletedCharacters;                                                          //!< Set of deleted characters
    std::string                         fileName;                                                                   //!< The path/filename from where this matrix originated
    std::vector<std::string>            sequenceNames;
    
    std::vector<BinaryTaxonData* >      taxonData;
    
};

#endif


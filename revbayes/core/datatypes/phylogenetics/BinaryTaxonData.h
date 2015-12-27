#ifndef BinaryTaxonData_H
#define BinaryTaxonData_H

#include <string>
#include <vector>

#include "numerics.h"
#include "Exception.h"
#include "Options.h"

class BinaryTaxonData {

public:
    BinaryTaxonData(void);                                                                                    //!< Set type spec of container from type of elements
    BinaryTaxonData(const std::string &tname);
    virtual ~BinaryTaxonData();
    
    virtual BinaryTaxonData*                clone(void) const = 0;

    virtual RealNumber                      operator[](size_t i) const = 0;                                         //!< Const index op
                   
    // TaxonData functions
    virtual void                            addCharacter(RealNumber newChar ) = 0;                             //!< Push back a new character
    virtual RealNumber                      getCharacter(size_t index) const = 0;
    bool                                    getGap(size_t index) const;
    void                                    setGap(size_t index, bool gap);
    virtual void                            setCharacter(size_t index, RealNumber val) = 0;
    virtual RealNumber                      getElement(size_t i) const = 0;                                         //!< Const index op
    virtual size_t                          getNumberOfCharacters(void) const = 0;                                  //!< How many characters
    const std::string&                      getTaxonName(void) const;                                           //!< Return the name of the character vector
    void                                    setTaxonName(std::string tn);                                       //!< Set the taxon name
    size_t                                  size(void) const;
    
    size_t                                  getNonGapLength() const;
    virtual RealNumber                      computeStateFrequency() const = 0;
    
protected:
    std::string                             taxonName;                                                          //!< Name of the taxon for this vector of characters               
    std::vector<bool>                       gaps;

};


#endif

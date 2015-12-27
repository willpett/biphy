#ifndef DiscreteBinaryTaxonData_H
#define DiscreteBinaryTaxonData_H

#include <string>
#include <vector>

#include "numerics.h"
#include "Exception.h"
#include "Options.h"
#include "BinaryTaxonData.h"

class DiscreteBinaryTaxonData : public BinaryTaxonData {

public:
    DiscreteBinaryTaxonData(void);                                                                                    //!< Set type spec of container from type of elements
    DiscreteBinaryTaxonData(const std::string &tname);
    virtual ~DiscreteBinaryTaxonData();
    
    DiscreteBinaryTaxonData*        clone(void) const;

    virtual RealNumber              operator[](size_t i) const;                                         //!< Const index op
                   
    // TaxonData functions
    void                            addCharacter(RealNumber newChar );                             //!< Push back a new character
    RealNumber                      getCharacter(size_t index) const;
    void                            setCharacter(size_t index, RealNumber val);
    RealNumber                      getElement(size_t i) const;                                         //!< Const index op
    size_t                          getNumberOfCharacters(void) const;                                  //!< How many characters
    
    RealNumber                      computeStateFrequency() const;
    
protected:
    std::vector<bool>                       sequence;

};


#endif

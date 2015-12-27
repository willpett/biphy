#ifndef ContinuousBinaryTaxonData_H
#define ContinuousBinaryTaxonData_H

#include <string>
#include <vector>

#include "numerics.h"
#include "Exception.h"
#include "Options.h"

#include "BinaryTaxonData.h"

class ContinuousBinaryTaxonData : public BinaryTaxonData {

public:
    ContinuousBinaryTaxonData(void);                                                                                    //!< Set type spec of container from type of elements
    ContinuousBinaryTaxonData(const std::string &tname);                                                                //!< Set type spec of container from type of elements

    ContinuousBinaryTaxonData*        		clone(void) const;

    RealNumber                              operator[](size_t i) const;                                         //!< Const index op
                   
    // TaxonData functions
    void                                    addCharacter(RealNumber newChar );                             //!< Push back a new character
    RealNumber                              getCharacter(size_t index) const;
    void                                    setCharacter(size_t index, RealNumber val);
    RealNumber                              getElement(size_t i) const;
    size_t                                  getNumberOfCharacters(void) const;
    
    RealNumber                              computeStateFrequency() const;
    
private:              
    std::vector<RealNumber>                 sequence;

};

#endif

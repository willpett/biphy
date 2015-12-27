#ifndef Clade_H
#define Clade_H

#include <vector>
#include <string>

class Clade  {
    
    public:
        Clade(void);                                                //! Default constructor: empty clade of age 0.0
        Clade(const std::string &t);                        //!< Default constructor with optional index
        Clade(const std::vector<std::string> &n);               //!< Default constructor with optional index
    
    virtual                                    ~Clade() {}
    
    std::vector<std::string>::const_iterator          begin(void) const;
    std::vector<std::string>::iterator                begin(void);
    std::vector<std::string>::const_iterator          end(void) const;
    std::vector<std::string>::iterator                end(void);
    // overloaded operators
    bool                                        operator==(const Clade &t) const;
    bool                                        operator!=(const Clade &t) const;
    bool                                        operator<(const Clade &t) const;
    bool                                        operator<=(const Clade &t) const;
    bool                                        operator>(const Clade &t) const;
    bool                                        operator>=(const Clade &t) const;

    
    // Basic utility functions
    Clade*                                      clone(void) const;                                          //!< Clone object
    
    // public methods
    void                                        addTaxon(const std::string &t);                                        //!< Get the age of this clade.
    std::vector<std::string>&                   getTaxa(void);                                              //!< Get the taxon names.
    const std::vector<std::string>&             getTaxa(void) const;                                 //!< Get a single taxon name.
    const std::string&                          getTaxonName(size_t i) const;                         //!< Set a single taxon's age.
    size_t                                      size(void) const;                                           //!< Get the number of taxa.
    std::string                                 toString(void) const;                                       //!< Convert this value into a string.
    
    // public TopologyNode functions
    
private: 
    
    // members
    std::vector<std::string>                    taxa;
    
};

// Global functions using the class
std::ostream&                       operator<<(std::ostream& o, const Clade& x);

#endif

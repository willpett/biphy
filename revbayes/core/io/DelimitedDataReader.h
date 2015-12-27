#ifndef DelimitedDataReader_H
#define DelimitedDataReader_H

#include <vector>
#include <string>

class DelimitedDataReader {
    
public:
    
    
protected:
    // protected constructor to disallow construction
    DelimitedDataReader(const std::string &fn, char d='\t');

    
    // protected methods only callable for derived classes
    void                                                readData(void);
    const std::vector<std::vector<std::string> >&       getChars(void);

    
    // protected member only accessible for derived classes
    std::string                                         filename;
    char                                                delimiter;
    std::vector<std::vector<std::string> >              chars;
    
};
#endif 

#ifndef NewickTreeReader_H
#define NewickTreeReader_H

#include <string>
#include <vector>

class NewickTreeReader {
    
public:
    NewickTreeReader();                                                                                     //!< Default constructor.
    
    std::vector<Tree*>*     readTrees(const std::string &fn);                   //!< Read a set of trees with branch lengths in newick format from a file.
    
};

#endif

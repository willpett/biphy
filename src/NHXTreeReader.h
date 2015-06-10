#ifndef NHXTreeReader_H
#define NHXTreeReader_H

#include <string>
#include <vector>

namespace RevBayesCore {
    
    class BranchLengthTree;
    
    /**
     * NHX tree reader.
     *
     * The NHX tree reader provides convenience methods for reading trees in NHX format.
     *
     *
     * @copyright Copyright 2009-
     * @author The RevBayes Development Core Team (Sebastian Hoehna)
     * @since 2014-01-29, version 1.0
     *
     */
    class NHXTreeReader {
        
    public:
        NHXTreeReader();                                                                                     //!< Default constructor.
        
        std::vector<BranchLengthTree*>*     readBranchLengthTrees(const std::string &fn);                   //!< Read a set of trees with branch lengths in NHX format from a file.
        std::vector<Topology*>*         	readTopologies(const std::string &fn);

    };
    
}

#endif

/**
 * @file
 * This file contains the declaration of our tree interface, a light-weight class 
 * that holds the pointer to the root node of a tree and provides some convinience functions.
 *
 * The derived classes differ mainly in which type of topology nodes they store.
 *
 * @todo: If the derived classes only differ in the topology nodes, then we could use templates for the derived classes instead!
 *
 * @brief Declaration of the Tree
 *
 * (c) Copyright 2009- under GPL version 3
 * @date Last modified: $Date: 2012-07-05 16:47:08 +0200 (Thu, 05 Jul 2012) $
 * @author The RevBayes Development Core Team
 * @license GPL version 3
 * @version 1.0
 * @since 2011-04-13, version 1.0
 *
 * $Id: Tree.h 1651 2012-07-05 14:47:08Z hoehna $
 */

#ifndef Tree_H
#define Tree_H

#include "Cloneable.h"
#include "TreeChangeEventHandler.h"

#include <vector>
#include <string>

namespace RevBayesCore {
    
    class Topology;
    class TopologyNode;

    class Tree : public Cloneable {
    
    public:       
        virtual                                            ~Tree(void);                                                                                                        //!< Destructor
        
        // overloaded operators
        bool                                                operator==(const Tree &t) const;
        bool                                                operator!=(const Tree &t) const;
        bool                                                operator<(const Tree &t) const;
        
        // virtual basic utility functions
        virtual Tree*                                       clone(void) const = 0;                                                      //!< Clone object
    
        // virtual Tree functions
        virtual double                                      getAge(size_t idx) const = 0;
        virtual double                                      getBranchLength(size_t idx) const = 0;
        virtual double                                      getTime(size_t idx) const = 0;

        // public Tree methods
        void                                                addBranchParameter(const std::string &n, const std::vector<double> &p, bool io);
        void                                                addNodeParameter(const std::string &n, const std::vector<double> &p, bool io);
        void                                                clearBranchParameters(void);
        std::vector<std::string>                            getTipNames() const;
        TopologyNode&                                       getNode(size_t idx);                                                    //!< Get the node at index
        const TopologyNode&                                 getNode(size_t idx) const;                                              //!< Get the node at index
        std::vector<TopologyNode*>                          getNodes(void) const;                                                   //!< Get a pointer to the nodes in the Tree
//        TopologyNode&                                       getNodeByIndex(size_t idx);
//        const TopologyNode&                                 getNodeByIndex(size_t idx) const;
//        std::vector<TopologyNode*>                          getNodesByIndex(void) const;
        size_t                                              getNumberOfInteriorNodes(void) const;                                   //!< Get the number of nodes in the Tree
        size_t                                              getNumberOfNodes(void) const;                                           //!< Get the number of nodes in the Tree
        size_t                                              getNumberOfTips(void) const;                                            //!< Get the number of tip nodes in the Tree
        const TopologyNode&                                 getInteriorNode(size_t indx) const;                                     //!< Get a pointer to interior node i
        const std::string&                                  getNewickRepresentation() const;                                        //!< Get the newick representation of this Tree
        TopologyNode&                                       getRoot(void);                                                          //!< Get a pointer to the root node of the Tree
        const TopologyNode&                                 getRoot(void) const;                                                    //!< Get a pointer to the root node of the Tree
        TopologyNode&                                       getTipNode(size_t indx);                                                //!< Get a pointer to tip node i
        const TopologyNode&                                 getTipNode(size_t indx) const;                                          //!< Get a pointer to tip node i
        double                                              getTmrca(const TopologyNode &n);
        TreeChangeEventHandler&                             getTreeChangeEventHandler(void) const;                                      //!< Get the change-event handler for this tree
        const Topology&                                     getTopology(void) const;                                                //!< Get the topology of this tree
        bool                                                hasSameTopology( const Tree &t ) const;                                 //!< Has this tree the same topology?
        bool                                                hasSameTopology( const Topology &t ) const;                             //!< Has this tree the same topology?
        bool                                                isBinary(void) const;                                                   //!< Is the Tree rooted
        bool                                                isRooted(void) const;                                                   //!< Is the Tree rooted
//        void                                                setRooted(bool tf);                                                     //!< Set the rootedness of the Tree
//        void                                                setRoot(TopologyNode* r);                                               //!< Set the root and bootstrap the Tree from it
        void                                                setTopology(const Topology *t, bool owns);
        
    protected:
        // preventing instantiation of the base class
        Tree(void);                                                                                                         //!< Default constructor
        Tree(const Tree& t);                                                                                                //!< Copy constructor
        Tree&                                               operator=(const Tree& t);

        // pure virtual methods you have to overwrite
        virtual void                                        resizeElementVectors(size_t n) = 0;

        // protected members
        const Topology*                                     topology;
        mutable TreeChangeEventHandler                      changeEventHandler;
        bool                                                ownsTopology;
    };

    // Global functions using the class
    std::ostream&                       operator<<(std::ostream& o, const Tree& x);                                         //!< Overloaded output operator

}

#endif


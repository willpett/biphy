#ifndef Model_H
#define Model_H

#include <map>
#include <set>
#include <vector>

#include "Cloneable.h"

class DagNode;

/**
 * \brief The Model class holds its independent copy of the model graph (DAG nodes) contained in the model
 * and provides methods for convenient access.
 *
 * A Model object holds the model graph which is an independent copy
 * of the DAG nodes contained in the model. The model graph may or may not be connected.
 * The model graph is obtained by providing a set of DAG nodes to the model and the model
 * object extracts all DAG nodes that are connected to these provided nodes.
 * An independent copy of these nodes is then obtained. 
 * The model object simply provides access to these DAG nodes.
 *
 * \copyright (c) Copyright 2009-2013 (GPL version 3)
 * \author The RevBayes Development Core Team (Sebastian Hoehna)
 * \since Version 1.0, 2012-06-21
 *
 */
class Model : public Cloneable {
    
    public:
                                                                Model(DagNode* source);                                   //!< Constructor from a single DAG node from which the model graph is extracted.
                                                                Model(std::set<DagNode*> sources);                  //!< Constructor from a set of DAG nodes from each of which the model graph is extracted.
                                                                Model(const Model &m);                                          //!< Copy constructor.
    virtual                                                    ~Model(void);                                                    //!< Destructor.

    
    // overloaded operators
    Model&                                                      operator=(const Model& x);                                      //!< Assignment operator.
    
    
    // convenience methods
    virtual Model*                                              clone(void) const;                                              //!< Clone function. This is similar to the copy constructor but useful in inheritance.

    // getters and setters
    std::vector<DagNode*>&                                      getDagNodes(void);                                              //!< Non-constant getter function of the set of DAG nodes contained in the model graph.
    const std::vector<DagNode*>&                                getDagNodes(void) const;                                        //!< Constant getter function of the map between the pointer of the original DAG nodes to the pointers of the copied DAG nodes.  

private:
    
    // private methods
    void                                                        addSourceNode(DagNode *sourceNode);                       //!< Add a source node, extract the model graph and create and indepedent copy of it.
    void                                                        cloneSourceNode(DagNode *sourceNode);
    // members
    std::vector<DagNode*>                                       nodes;
    std::set<DagNode*>                                          sources;                                                        //!< Set of source nodes for the model graph.
};

#endif

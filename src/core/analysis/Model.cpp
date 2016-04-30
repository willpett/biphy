#include "Model.h"

#include "../utils/Exception.h"
#include "DagNode.h"


/**
 * Constructor from a single DAG node.
 * The model graph is extracted by obtaining all DAG nodes connected to the provide source node.
 * The entire model graph is copied and a map between the pointers to the original DAG nodes and
 * the copied DAG nodes is created for convinience access.
 *
 * \param[in]    source    The DAG node from which the model graph is extracted.
 */
Model::Model(DagNode *source) 
{
    
    // add this node to the source nodes and build model graph
    addSourceNode( source );
    
}


/**
 * Constructor from a set of DAG nodes.
 * The model graph is extracted by obtaining all DAG nodes that are connected to either of the provide source nodes.
 * The entire model graph is copied and a map between the pointers to the original DAG nodes and
 * the copied DAG nodes is created for convinience access.
 *
 * \param[in]    sources    The set of DAG nodes from which the model graph is extracted.
 */
Model::Model(std::set<DagNode*> s) : 
    sources() 
{
    
    // iterate over all sources
    for (std::set<DagNode*>::const_iterator it = s.begin(); it != s.end(); ++it) 
    {
        // add this node and build model graph
        addSourceNode( *it );
    }
    
}


/**
 * Copy constructor. We instantiate the model from the previously stored source nodes. 
 *
 * \param[in]    m    The model object to copy.
 */
Model::Model(const Model &m) : sources() 
{
    
    // iterate over all sources
    for (std::set<DagNode*>::const_iterator it = m.sources.begin(); it != m.sources.end(); ++it) 
    {
        // add this node and build model graph
        cloneSourceNode( *it );
    }
    
}

/**
 * Destructor.
 * We have created new copied of the DAG nodes so we need to delete these here again.
 */
Model::~Model( void ) 
{
    
    // delete each DAG node from the copied model graph.
    for (std::vector<DagNode*>::iterator it = nodes.begin(); it != nodes.end(); ++it) 
    {
        DagNode *theNode = *it;
        if ( theNode->decrementReferenceCount() == 0 )
        {
            delete theNode;
        }
    }
    
    while ( !sources.empty() ) 
    {
        std::set<DagNode*>::iterator theNode = sources.begin();
        sources.erase( theNode );
        
        if ( (*theNode)->decrementReferenceCount() == 0)
        {
            delete *theNode;
        }
    }
    
}


/**
 * Assignment operator.
 * We have stored previously all the source nodes so we can simply construct the model graph by extracting
 * the DAG again from these source nodes. Note that we are acutally making a copy of the original graph again
 * which may cause that we get a different model graph (e.g. if the graph has change in the meanwhile).
 * The reason behind this approach is that we keep a valid map from the pointers of the original DAG nodes to
 * the pointers to the copied DAG nodes.
 *
 * \todo Check that this copy constructor is good or if we should use a different mechanism for the map between the
 * original DAG nodes and the copied DAG nodes.
 *
 * \param[in]    m    The model object to copy.
 */
Model& Model::operator=(const Model &x) 
{
    
    if ( this != &x )
    {
        // first remove all DAG nodes
        for (std::vector<DagNode*>::iterator it = nodes.begin(); it != nodes.end(); ++it) 
        {
            DagNode *theNode = *it;
            if ( theNode->decrementReferenceCount() == 0 )
            {
                delete theNode;
            }
        }
        
        // empty the source nodes
        while ( !sources.empty() ) 
        {
            std::set<DagNode*>::iterator theNode = sources.begin();
            sources.erase( theNode );
            
            if ( (*theNode)->decrementReferenceCount() == 0)
            {
                delete *theNode;
            }
        }
        
        // iterate over all sources
        for (std::set<DagNode*>::const_iterator it = x.sources.begin(); it != x.sources.end(); ++it) 
        {
            // add this node and build model graph
            cloneSourceNode( *it );
        }
    }
    
    return *this;
}

void Model::addSourceNode( DagNode *sourceNode) 
{
    
    // check that the source node is a valid pointer
    if (sourceNode == NULL)
        throw Exception("Cannot instantiate a model with a NULL DAG node.");
    
    // add the source node to our set of sources
    sourceNode->incrementReferenceCount();
    sources.insert( sourceNode );
    
    // we don't really know which nodes are new in our nodes map.
    // therefore we empty the nodes map and fill it again.
    for (std::vector<DagNode*>::iterator it = nodes.begin(); it != nodes.end(); ++it) 
        (*it)->decrementReferenceCount();
    
    nodes.clear();
    std::set<const DagNode*> nodeset;
    sourceNode->collectDAG(nodeset);
    
    std::set<const DagNode*>::iterator i = nodeset.begin();
    while ( i != nodeset.end() ) 
    {
        // get the copied node
        DagNode* theNewNode = const_cast<DagNode*>(*i);
        
        // increment the iterator;
        ++i;
        
        // increment the reference count to the new node
        theNewNode->incrementReferenceCount();
            
        // insert in direct access vector
        nodes.push_back( theNewNode );
    }
    
    
}

void Model::cloneSourceNode( DagNode *sourceNode) 
{
    // check that the source node is a valid pointer
    if (sourceNode == NULL)
        throw Exception("Cannot instantiate a model with a NULL DAG node.");
    
    std::map<const DagNode*, DagNode* > nodesMap;
    // copy the entire graph connected to the source node
    // only if the node is not contained already in the nodesMap will it be copied.
    sourceNode->cloneDAG(nodesMap);
    
    // add the source node to our set of sources
    DagNode *theNewSource = nodesMap[sourceNode];
    theNewSource->incrementReferenceCount();
    sources.insert( theNewSource );
        
    /* insert new nodes into direct access vector */
    std::map<const DagNode*, DagNode* >::iterator i = nodesMap.begin();
    
    // we don't really know which nodes are new in our nodes map.
    // therefore we empty the nodes map and fill it again.
    for (std::vector<DagNode*>::iterator it = nodes.begin(); it != nodes.end(); ++it) 
    {
        
        DagNode *theNode = *it;
        theNode->decrementReferenceCount();
        
    }
    nodes.clear();
    while ( i != nodesMap.end() ) 
    {
        // get the copied node
        DagNode* theNewNode = (*i).second;
        
        // increment the iterator;
        ++i;
        
        // increment the reference count to the new node
        theNewNode->incrementReferenceCount();
            
        // insert in direct access vector
        nodes.push_back( theNewNode );
    }
}

/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'B'.
 *
 * \return A new copy of the model. 
 */
Model* Model::clone( void ) const 
{
    
    return new Model( *this );
}


/**
 * Constant getter function for the vector of DAG nodes of the model.
 *
 * \return Vector of DAG nodes constituting to this model.
 */
const std::vector<DagNode *>& Model::getDagNodes( void ) const 
{
    
    return nodes;
}


/**
 * Non-constant getter function for the vector of DAG nodes of the model.
 *
 * \return Vector of DAG nodes constituting to this model.
 */
std::vector<DagNode *>& Model::getDagNodes( void ) 
{
    
    return nodes;
}

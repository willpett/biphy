#include "TreeAssemblyFunction.h"
#include "RbException.h"

using namespace RevBayesCore;

TreeAssemblyFunction::TreeAssemblyFunction(const TypedDagNode<Topology> *t, const TypedDagNode<std::vector<double> > *b) : TypedFunction<BranchLengthTree>( new BranchLengthTree() ), tau( t ), brlen( b ) {
    // add the lambda parameter as a parent
    addParameter( tau );
    addParameter( brlen );
    
    value->setTopology( &(tau->getValue()), false );
    
    update();
}


TreeAssemblyFunction::TreeAssemblyFunction(const TreeAssemblyFunction &n) : TypedFunction<BranchLengthTree>( n ), tau( n.tau ), brlen( n.brlen ) {
    // no need to add parameters, happens automatically
}


TreeAssemblyFunction::~TreeAssemblyFunction( void ) {
    // We don't delete the parameters, because they might be used somewhere else too. The model needs to do that!
}



TreeAssemblyFunction* TreeAssemblyFunction::clone( void ) const {
    return new TreeAssemblyFunction( *this );
}


void TreeAssemblyFunction::keep(DagNode *affecter) {
    
//    touchedNodeIndices.clear();
}


void TreeAssemblyFunction::restore(DagNode *restorer) {
    
//    touchedNodeIndices.clear();
}


void TreeAssemblyFunction::touch(DagNode *toucher) {
    
    if ( toucher == brlen ) {
        const std::set<size_t> &touchedIndices = toucher->getTouchedElementIndices();
        touchedNodeIndices.insert(touchedIndices.begin(), touchedIndices.end());
    }
}


void TreeAssemblyFunction::update( void ) {
    
    if ( touchedNodeIndices.size() > 0 ) {
        const std::vector<double> &v = brlen->getValue();
        for (std::set<size_t>::iterator it = touchedNodeIndices.begin(); it != touchedNodeIndices.end(); ++it) {
            value->setBranchLength(*it, v[*it]);
        }
        touchedNodeIndices.clear();
    } else {
        const std::vector<double> &v = brlen->getValue();
        for (size_t i = 0; i < v.size(); ++i) {
            value->setBranchLength(i, v[i]);
        }
    }
}



void TreeAssemblyFunction::swapParameterInternal(const DagNode *oldP, const DagNode *newP) {
    if (oldP == tau) {
        tau = static_cast<const TypedDagNode<Topology>* >( newP );
        value->setTopology( &(tau->getValue()), false );
    }
    else if (oldP == brlen) {
        brlen = static_cast<const TypedDagNode<std::vector<double> >* >( newP );
    }
}



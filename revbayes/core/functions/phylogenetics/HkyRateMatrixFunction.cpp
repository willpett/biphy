#include "HkyRateMatrixFunction.h"
#include "RbException.h"

using namespace RevBayesCore;

HkyRateMatrixFunction::HkyRateMatrixFunction(const TypedDagNode<double> *k, const TypedDagNode<std::vector<double> > *bf) : TypedFunction<RateMatrix>( new RateMatrix_HKY() ), baseFrequencies( bf ), kappa( k ) {
    // add the lambda parameter as a parent
    addParameter( baseFrequencies );
    addParameter( kappa );
    
    update();
}


HkyRateMatrixFunction::HkyRateMatrixFunction(const HkyRateMatrixFunction &n) : TypedFunction<RateMatrix>( n ), baseFrequencies( n.baseFrequencies ), kappa( n.kappa ) {
    // no need to add parameters, happens automatically
}


HkyRateMatrixFunction::~HkyRateMatrixFunction( void ) {
    // We don't delete the parameters, because they might be used somewhere else too. The model needs to do that!
}



HkyRateMatrixFunction* HkyRateMatrixFunction::clone( void ) const {
    return new HkyRateMatrixFunction( *this );
}


void HkyRateMatrixFunction::update( void ) {
    // get the information from the arguments for reading the file
    double k = kappa->getValue();
    const std::vector<double>& f = baseFrequencies->getValue();
    
    
    // set the base frequencies
    static_cast< RateMatrix_HKY* >(value)->setStationaryFrequencies( f );
    static_cast< RateMatrix_HKY* >(value)->setKappa( k );
    
    value->updateMatrix();
}



void HkyRateMatrixFunction::swapParameterInternal(const DagNode *oldP, const DagNode *newP) {
    if (oldP == baseFrequencies) {
        baseFrequencies = static_cast<const TypedDagNode<std::vector<double> >* >( newP );
    }
    else if (oldP == kappa) {
        kappa = static_cast<const TypedDagNode< double >* >( newP );
    }
}



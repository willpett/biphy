#include "CoalaFunction.h"
#include "RbException.h"
#include "RbMathVector.h"

using namespace RevBayesCore;

CoalaFunction::CoalaFunction(const TypedDagNode<std::vector<double> > *coords, const MatrixReal &ca, const std::vector<double> &cw) : TypedFunction< std::vector<double> >( new std::vector<double>(coords->getValue().size()) ), coordinates( coords ), coa( ca ), colWeights( cw ) {
    // add the coordinates parameter as a parent
    addParameter( coordinates );
    
    update();
}


CoalaFunction::CoalaFunction(const CoalaFunction &n) : TypedFunction< std::vector<double> >( n ), coordinates( n.coordinates ), coa( n.coa ), colWeights( n.colWeights ) {
    // no need to add parameters, happens automatically
}


CoalaFunction::~CoalaFunction( void ) {
    // We don't delete the parameters, because they might be used somewhere else too. The model needs to do that!
}



CoalaFunction* CoalaFunction::clone( void ) const {
    return new CoalaFunction( *this );
}


void CoalaFunction::update( void ) {
    
    // get the information from the arguments for reading the file
    const std::vector<double>& c = coordinates->getValue();
    
    // Now, frequencies are computed from the vector of coordinates and the transpose of the principal axes matrix (P_):
    std::vector<double> tmpFreqs = coa * c;
    for (unsigned int i = 0; i < tmpFreqs.size(); i++)
    {
        tmpFreqs[i] = (tmpFreqs[i] + 1) * colWeights[i];
    }
    *value = tmpFreqs;
//    std::vector<double> &freq = *value;
//    freq = tmpFreqs;
//    
//    // Frequencies are not allowed to be lower than 10^-3 or higher than 0.5:
//    bool norm = false;
//    for (unsigned int i = 0; i < 20; i++)
//    {
//        if (freq[i] < 0.001)
//        {
//            norm = true;
//            freq[i] = 0.001;
//        }
//        if (freq[i] > 0.2)
//        {
//            norm = true;
//            freq[i] = 0.2;
//        }
//    }
//    if (norm == true)
//    {
//        RbMath::normalize( freq, 1.0 );
//    }
}



void CoalaFunction::swapParameterInternal(const DagNode *oldP, const DagNode *newP) {
    
    if ( oldP == coordinates ) 
    {
        coordinates = static_cast<const TypedDagNode<std::vector<double> >* >( newP );
    }
}



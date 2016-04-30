#include "MappingMonitor.h"
#include "BinarySubstitutionModel.h"

#include <sstream>

MappingMonitor::MappingMonitor(StochasticNode<BinaryCharacterData> *t, int g, const std::string &fname, bool ap) : Monitor(g,t), outStream(), data( t ), filename( fname ), append(ap) {
    if(dynamic_cast<BinarySubstitutionModel*>(&(t->getDistribution())) == 0){
        throw Exception("MappingMonitor requires BinarySubstitutionModel");
    }
}


MappingMonitor::MappingMonitor(const MappingMonitor &m) : Monitor( m ), outStream(), data( m.data ) {

    filename    = m.filename;
    append      = m.append;

    //if (m.outStream.is_open())
    //    openStream();
}


/* Clone the object */

MappingMonitor* MappingMonitor::clone(void) const {

    return new MappingMonitor(*this);
}


void MappingMonitor::closeStream() {
    outStream.close();
}


/** Monitor value at generation gen */

void MappingMonitor::monitor(long gen) {

    // get the printing frequency
    int samplingFrequency = printgen;

    if (gen % samplingFrequency == 0) {
        BinarySubstitutionModel* model = dynamic_cast<BinarySubstitutionModel* >(&data->getDistribution());
        const std::vector< DiscreteBinaryTaxonData >& mapping = model->getMapping();

        std::vector<TopologyNode*> nodes = model->getTree()->getNodes();

        for(size_t i = 0; i < nodes.size(); i++)
        {
            const DiscreteBinaryTaxonData& d = mapping[i];
            for(size_t c = 0; c < d.size(); c++)
            	outStream << d[c];

            if(i != nodes.size() - 1)
                outStream << "\t";
        }

        outStream << std::endl;

        //tree->clearBranchParameters();
        //outStream << tree->getNewickRepresentation();
        //outStream << std::endl;
        //delete tree;
    }
}


/** open the file stream for printing */

void MappingMonitor::openStream(void) {

    // open the stream to the file
    if (append)
        outStream.open( filename.c_str(), std::fstream::out | std::fstream::app);
    else
        outStream.open( filename.c_str(), std::fstream::out);
}

/** Print header for monitored values */

void MappingMonitor::printHeader() {
    BinarySubstitutionModel* model = dynamic_cast<BinarySubstitutionModel* >(&data->getDistribution());
    const std::vector< DiscreteBinaryTaxonData >& mapping = model->getMapping();

    size_t nodes = model->getTree()->getNumberOfNodes();
    for(size_t i=0; i<nodes; i++)
    {
        outStream << i;
        if(i != nodes - 1)
            outStream << "\t";
    }
    outStream << std::endl;
}


void MappingMonitor::swapNode(DagNode *oldN, DagNode *newN) {
    if ( oldN == data ) {
        data = static_cast< StochasticNode<BinaryCharacterData> *>( newN );
    }
    // delegate to base class
    Monitor::swapNode(oldN, newN);
}

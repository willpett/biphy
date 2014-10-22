#ifndef BinaryCharEvoModel_H
#define BinaryCharEvoModel_H

#include "AbstractSiteCorrectionModel.h"
#include "StandardState.h"
#include "RateMatrix.h"
#include "RbVector.h"
#include "TopologyNode.h"
#include "TransitionProbabilityMatrix.h"
#include "TreeChangeEventListener.h"
#include "TypedDistribution.h"

namespace RevBayesCore {
    
	enum CorrectionType { NONE 					= 0x00,
						  NO_ABSENT_SITES		= 0x01,
						  NO_PRESENT_SITES		= 0x02,
						  NO_CONSTANT_SITES		= 0x03,
						  NO_SINGLETON_GAINS 	= 0x04,
						  NO_SINGLETON_LOSSES	= 0x08,
						  NO_UNINFORMATIVE		= 0x0F
						};

    template<class treeType>
    class BinaryCharEvoModel : public AbstractSiteCorrectionModel<StandardState, treeType> {

    public:
    	BinaryCharEvoModel(const TypedDagNode< treeType > *p, bool c, size_t nSites, int type = 0);
        BinaryCharEvoModel(const BinaryCharEvoModel &n);                                                                                                //!< Copy constructor
        virtual                                            ~BinaryCharEvoModel(void);                                                                   //!< Virtual destructor
        
        // public member functions
        BinaryCharEvoModel*         						clone(void) const;
        
        void                                         		setCorrectionPatterns(AbstractCharacterData* data);

    protected:
        
        
    private:

	   int													type;

    };
    
}


#include "ConstantNode.h"
#include "DiscreteCharacterState.h"
#include "RandomNumberFactory.h"
#include "TopologyNode.h"
#include "TransitionProbabilityMatrix.h"

#include <cmath>
#include <cstring>

template<class treeType>
RevBayesCore::BinaryCharEvoModel<treeType>::BinaryCharEvoModel(const TypedDagNode< treeType > *t, bool c, size_t nSites, int type) :
	AbstractCharEvoModel<StandardState, treeType>(t, 2, c, nSites),
	AbstractSiteCorrectionModel<StandardState, treeType>(t, 2, c , nSites), type(type)
{
	setCorrectionPatterns(NULL);
}


template<class treeType>
RevBayesCore::BinaryCharEvoModel<treeType>::BinaryCharEvoModel(const BinaryCharEvoModel &d) :
	AbstractCharEvoModel<StandardState, treeType>( d ),
	AbstractSiteCorrectionModel<StandardState, treeType>(d), type(d.type)
{
}


template<class treeType>
RevBayesCore::BinaryCharEvoModel<treeType>::~BinaryCharEvoModel( void ) {
    // We don't delete the parameters, because they might be used somewhere else too. The model needs to do that!
    
}


template<class treeType>
RevBayesCore::BinaryCharEvoModel<treeType>* RevBayesCore::BinaryCharEvoModel<treeType>::clone( void ) const {
    
    return new BinaryCharEvoModel<treeType>( *this );
}

template<class treeType>
void RevBayesCore::BinaryCharEvoModel<treeType>::setCorrectionPatterns(AbstractCharacterData* data){

	size_t numTips = this->tau->getValue().getNumberOfTips();

	this->numCorrectionSites = (bool)(type & NO_ABSENT_SITES)
							 + (bool)(type & NO_PRESENT_SITES)
							 + numTips*(bool)(type & NO_SINGLETON_GAINS)
							 + numTips*(bool)(type & NO_SINGLETON_LOSSES);

	this->resizeLikelihoodVectors();

	std::vector<TopologyNode*> nodes = this->tau->getValue().getNodes();
	// allocate and fill the cells of the matrices

	size_t nonSingletons = (bool)(type & NO_ABSENT_SITES) + (bool)(type & NO_PRESENT_SITES);
	size_t i = 0;
	for (std::vector<TopologyNode*>::iterator it = nodes.begin(); it != nodes.end(); ++it)
	{
		if ( (*it)->isTip() )
		{
			std::string name = (*it)->getName();

			// resize the column
			this->correctionCharMatrix[name].resize(this->numCorrectionSites);
			this->correctionGapMatrix[name].resize(this->numCorrectionSites);

			for (size_t site = 0; site < this->numCorrectionSites; ++site)
			{
				this->correctionGapMatrix[name][site] = false;
				if(site < nonSingletons){
					this->correctionCharMatrix[name][site] = site + !(type & NO_ABSENT_SITES);
				}else{
					if(site - nonSingletons < numTips)
						if(site - nonSingletons == i){
							this->correctionCharMatrix[name][site] = !(type & NO_SINGLETON_LOSSES);
						}else{
							this->correctionCharMatrix[name][site] = (bool)(type & NO_SINGLETON_LOSSES);
						}
					else
						if(site - nonSingletons == numTips + i){
							this->correctionCharMatrix[name][site] = 1;
						}else{
							this->correctionCharMatrix[name][site] = 0;
						}
				}
			}
			i++;
		}
	}
}


#endif

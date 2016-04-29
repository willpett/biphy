#ifndef BinaryDolloSubstitutionModel_H
#define BinaryDolloSubstitutionModel_H

#include "BinarySubstitutionModel.h"


class BinaryDolloSubstitutionModel : public BinarySubstitutionModel {
    
public:
    // Note, we need the size of the alignment in the constructor to correctly simulate an initial state
    BinaryDolloSubstitutionModel(const TypedDagNode<Tree> *t, AscertainmentBias::Coding c = AscertainmentBias::ALL);
    BinaryDolloSubstitutionModel(const BinaryDolloSubstitutionModel &n);
    virtual                                                            ~BinaryDolloSubstitutionModel(void);
    
    // public member functions
    // pure virtual
    virtual BinaryDolloSubstitutionModel*                 				clone(void) const;
    // non-virtual
    virtual void                                                        redrawValue(void);
    const std::vector< DiscreteBinaryTaxonData >&                       getMapping(void);
    
    void                                                                setStationaryFrequency(const TypedDagNode< std::vector< double > > *r);
    void                                                                setStationaryFrequency(const TypedDagNode< double > *f);
    RealNumber                                                          getStationaryFrequency(size_t n);
    
protected:
    // helper method for this and derived classes
    RealVector                                                          integrationFactors;
    size_t                                                              activeIntegrationOffset;
    std::vector<std::vector<size_t> >                                   maskNodeObservationCounts;
    
    // virtual methods that may be overwritten, but then the derived class should call this methods
    void                                                                resizeLikelihoodVectors(void);
    virtual void                                                        getIncludedSiteIndices();
    
    // pure virtual methods
    virtual void                                                        computeNodeLikelihood(const TopologyNode &n, size_t nIdx, size_t l, size_t r);
    virtual void                                                        computeNodeLikelihood(const TopologyNode &n, size_t nIdx, size_t l, size_t r, size_t m);
           
    virtual void                                                        computeNodeCorrection(const TopologyNode &n, size_t nIdx, size_t l, size_t r);
    virtual void                                                        computeNodeCorrection(const TopologyNode &n, size_t nIdx, size_t l, size_t r, size_t m);
    virtual void                                                        computeNodeCorrection(const TopologyNode &node, size_t nIdx);
    
    virtual RealNumber                                                  sumRootLikelihood( void );
    virtual RealNumber                                                  sumUncorrectedRootLikelihood( void );
    
    RealNumber                                                          getAncestralNodeWeights(const TopologyNode &node, size_t pattern, std::vector<RealNumber>* weights = NULL, std::vector<size_t>* nodes = NULL);
        
    virtual void                                                        updateTransitionProbabilities();
    
    void                                                                simulateDolloMapping( const TopologyNode &node, size_t pattern, std::vector<bool> &data, size_t rateIndex, std::pair<size_t, size_t>& charCounts);
 
protected:
    // private methods
    void                                                                simulate( const TopologyNode &node, std::vector<RealNumber> &data, size_t rateIndex, std::pair<size_t, size_t>& charCounts);
    
    void                                                                scale(size_t i, size_t l, size_t r);
    void                                                                scale(size_t i, size_t l, size_t r, size_t m);
    void                                                                scaleCorrection(size_t i, size_t l, size_t r);
    void                                                                scaleCorrection(size_t i, size_t l, size_t r, size_t m);


};

#endif

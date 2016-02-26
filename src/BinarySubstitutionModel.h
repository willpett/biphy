#ifndef BinarySubstitutionModel_H
#define BinarySubstitutionModel_H

#include "BinaryCharacterData.h"
#include "DiscreteBinaryTaxonData.h"
#include "Tree.h"
#include "TreeChangeEventListener.h"
#include "TypedDistribution.h"
#include <memory.h>
#include <cmath>

#include "MathFunctions.h"

struct AscertainmentBias {
    
  enum Coding { ALL                 = 0x00,
                NOABSENCESITES      = 0x01,
                NOPRESENCESITES     = 0x02,
                VARIABLE            = 0x03,
                NOSINGLETONPRESENCE = 0x04,
                NOSINGLETONABSENCE  = 0x08,
                NOSINGLETONS        = 0x0C,
                INFORMATIVE         = 0x0F
              };
};


class BinarySubstitutionModel : public TypedDistribution< BinaryCharacterData >, public TreeChangeEventListener {
    
public:
    // Note, we need the size of the alignment in the constructor to correctly simulate an initial state
    BinarySubstitutionModel(const TypedDagNode<Tree> *t, AscertainmentBias::Coding c = AscertainmentBias::ALL);
    BinarySubstitutionModel(const BinarySubstitutionModel &n);
    virtual                                                            ~BinarySubstitutionModel(void);
    
    // public member functions
    // pure virtual
    virtual BinarySubstitutionModel*                 				    clone(void) const;
    
    // virtual (you need to overwrite this method if you have additional parameters)
    virtual void                                                        swapParameter(const DagNode *oldP, const DagNode *newP);
    
    // non-virtual
    virtual double                                                      computeLnProbability(void);
    void                                                                fireTreeChangeEvent(const TopologyNode &n);
    void                                                                setValue(BinaryCharacterData *v);
    virtual void                                                        redrawValue(void);
    virtual const std::vector< DiscreteBinaryTaxonData >&               getMapping(void);
    size_t                                                              getNumSites(void) const;
    size_t                                                              getNumPatterns(void) const;
    const Tree*												            getTree() const;
    
    RealNumber                                                          getLnCorrection() const;
    RealVector                                                          getPerSiteLnProbs() const;
    
    void                                                                setClockRate(const TypedDagNode< double > *r);
    void                                                                setClockRate(const TypedDagNode< std::vector< double > > *r);
    void                                                                setBranchLengths(const TypedDagNode< std::vector< double > > *r);
    virtual void                                                        setStationaryFrequency(const TypedDagNode< std::vector< double > > *r);
    virtual void                                                        setStationaryFrequency(const TypedDagNode< double > *f);
    virtual void                                                        setSiteRates(const TypedDagNode< std::vector< double > > *r);
    virtual void                                                        setSiteFrequencies(const TypedDagNode< std::vector< double > > *r);
    
    void                                                                setVerbose(bool v);
    void                                                                setUseScaling(bool s);
    double                                                              getClockRate(size_t n);
    double                                                              getBranchLength(size_t n);
    virtual RealNumber                                                  getStationaryFrequency(size_t n, size_t mixture);
    
protected:
    // helper method for this and derived classes
    void                                                                recursivelyFlagNodeDirty(const TopologyNode& n);
    virtual void                                                        resizeLikelihoodVectors(void);
    
    // virtual methods that may be overwritten, but then the derived class should call this methods
    virtual void                                                        keepSpecialization(DagNode* affecter);
    virtual void                                                        restoreSpecialization(DagNode *restorer);
    virtual void                                                        touchSpecialization(DagNode *toucher);
    
    // pure virtual methods
    virtual void                                                        computeNodeLikelihood(const TopologyNode &n, size_t nIdx, size_t l, size_t r);
    virtual void                                                        computeNodeLikelihood(const TopologyNode &n, size_t nIdx, size_t l, size_t r, size_t m);
    virtual void                                                        computeNodeLikelihood(const TopologyNode &node, size_t nIdx);
    
    virtual void                                                        computeNodeCorrection(const TopologyNode &n, size_t nIdx, size_t l, size_t r);
    virtual void                                                        computeNodeCorrection(const TopologyNode &n, size_t nIdx, size_t l, size_t r, size_t m);
    virtual void                                                        computeNodeCorrection(const TopologyNode &node, size_t nIdx);
    
    virtual RealNumber                                                  sumRootLikelihood( void );
    virtual RealNumber                                                  sumUncorrectedRootLikelihood( void );
    
    virtual void                                                        updateTransitionProbabilities();

    virtual void                                                        getIncludedSiteIndices();
    virtual bool                                                        isSitePatternCompatible( std::pair<size_t,size_t>, size_t gap );
    
    const TypedDagNode<Tree>*                                           tau;
    
    AscertainmentBias::Coding                                           coding;
    // members
    RealNumber                                                          lnProb;
    size_t                                                              numSites;
    size_t                                                              numNodes;
    size_t                                                              numTaxa;
    size_t                                                              numPatterns;
#ifdef SIMD_ENABLED
    size_t                                                              numAllocatedPatterns;
    size_t                                                              numSIMDBlocks;
#endif
    bool                                                                verbose;
    size_t                                                              scalingDensity;
    
    // RAS model
    size_t                                                              mixtureOffset;
    size_t                                                              rateOffset;
    size_t                                                              numSiteRates;
    size_t                                                              numSiteFrequencies;
    
    // Correction model
    
    size_t                                                              N;
    RealNumber                                                          perSiteCorrection;

    size_t                                                              numCorrectionMasks;
    size_t                                                              activeCorrectionOffset;
    size_t                                                              correctionNodeOffset;
    size_t                                                              correctionRateOffset;
    size_t                                                              correctionMixtureOffset;
    
    // full model
    size_t                                                              activeLikelihoodOffset;
    size_t                                                              nodeOffset;
    size_t                                                              siteOffset;
    
    size_t                                                              tNodeOffset;
    size_t                                                              tRateOffset;
    size_t                                                              tMixtureOffset;
    size_t                                                              tActiveOffset;
    
    bool                                                                rateVariationAcrossSites;
    bool                                                                frequencyVariationAcrossSites;
    bool                                                                branchHeterogeneousClockRates;
    bool                                                                branchHeterogeneousFrequencies;
    
    // convenience variables available for derived classes too
    std::set<size_t>                                                    changedNodes;
    std::vector<bool>                                                   dirtyNodes;
    std::set<size_t>                                                    touchedNodes;
    
    const TypedDagNode< std::vector< double > >*                        siteRates;
    
    // Branch model
    const TypedDagNode< double >*                                       homogeneousClockRate;
    const TypedDagNode< std::vector< double > >*                        heterogeneousClockRates;
    const TypedDagNode< std::vector< double > >*                        branchLengths;
    const TypedDagNode< std::vector< double > >*                        heterogeneousFrequencies;
    const TypedDagNode< double >*                                       homogeneousFrequency;
    
    RealVector                                                          correctionLikelihoods;
    
    std::vector<std::vector<bool> >                                     correctionMaskMatrix;
    std::vector<size_t>                                                 correctionMaskCounts;
    std::vector<size_t>                                                 maskObservationCounts;
    RealVector                                                          perMaskCorrections;

    RealVector                                                          perMixtureCorrections;
    
    RealVector                                                          perNodeCorrectionLogScalingFactors;
    size_t                                                              activeCorrectionScalingOffset;
    RealVector                                                          perNodeSiteLogScalingFactors;
    size_t                                                              activeScalingOffset;
    
    bool                                                                useScaling;
    
    // the likelihoods
    RealVector                                                          per_site_Likelihoods;
    RealVector                                                          partialLikelihoods;
    std::vector<bool>                                                   activeLikelihood;
    std::vector<bool>                                                   activeProbability;
    
    // the data
    std::vector<size_t>                                                 siteIndices;
    std::vector<size_t>                                                 patternCounts;
    std::vector<size_t>                                                 site2pattern;
    std::vector<size_t>                                                 pattern2site;
    std::vector<size_t>                                                 site2mask;
    RealVector                                                          transitionProbabilities;
    
    bool                                                                continuous;
    
    
    
protected:
    // private methods
    virtual void                                                        compress(void);
    void                                                                fillLikelihoodVector(const TopologyNode &n, size_t nIdx);
    virtual void                                                        simulate( const TopologyNode &node, std::vector<RealNumber> &data, size_t rateIndex, size_t freqIndex, std::pair<size_t, size_t>& charCounts);
    
    virtual void                                                        scale(size_t i, size_t l, size_t r);
    virtual void                                                        scale(size_t i, size_t l, size_t r, size_t m);
    virtual void                                                        scaleCorrection(size_t i, size_t l, size_t r);
    virtual void                                                        scaleCorrection(size_t i, size_t l, size_t r, size_t m);
    
    std::vector< DiscreteBinaryTaxonData >                              mapping;
    
private:
    virtual void                                                        simulateMapping(const TopologyNode& parent, const TopologyNode& node, std::vector<size_t>& perSiteRates, std::vector<size_t>& perSiteFreqs);


};

#endif

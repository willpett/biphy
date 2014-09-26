#ifndef NucleotideBranchHeterogeneousCharEvoModel_H
#define NucleotideBranchHeterogeneousCharEvoModel_H

#include "AbstractSiteHomogeneousMixtureCharEvoModel.h"
#include "DnaState.h"
#include "RateMatrix.h"
#include "RbVector.h"
#include "TopologyNode.h"
#include "TransitionProbabilityMatrix.h"
#include "TreeChangeEventListener.h"
#include "TypedDistribution.h"

namespace RevBayesCore {
    
    template<class charType, class treeType>
    class NucleotideBranchHeterogeneousCharEvoModel : public AbstractSiteHomogeneousMixtureCharEvoModel<charType, treeType> {
        
    public:
        NucleotideBranchHeterogeneousCharEvoModel(const TypedDagNode< treeType > *t, bool c, size_t nSites);
        NucleotideBranchHeterogeneousCharEvoModel(const NucleotideBranchHeterogeneousCharEvoModel &n);                                                                                                //!< Copy constructor
        virtual                                            ~NucleotideBranchHeterogeneousCharEvoModel(void);                                                                   //!< Virtual destructor
        
        // public member functions
        NucleotideBranchHeterogeneousCharEvoModel*          clone(void) const;                                                                          //!< Create an independent clone
        void                                                setClockRate(const TypedDagNode< double > *r);
        void                                                setClockRate(const TypedDagNode< std::vector< double > > *r);
        void                                                setRateMatrix(const TypedDagNode< RateMatrix > *rm);
        void                                                setRateMatrix(const TypedDagNode< RbVector< RateMatrix > > *rm);
        void                                                setRootFrequencies(const TypedDagNode< std::vector< double > > *f);
        void                                                setSiteRates(const TypedDagNode< std::vector< double > > *r);
        void                                                swapParameter(const DagNode *oldP, const DagNode *newP);                                    //!< Implementation of swaping parameters
        
    protected:
        
        void                                                computeRootLikelihood(size_t root, size_t l, size_t r);
        void                                                computeInternalNodeLikelihood(const TopologyNode &n, size_t nIdx, size_t l, size_t r);
        void                                                computeTipLikelihood(const TopologyNode &node, size_t nIdx);
        const std::vector<double>&                          getRootFrequencies(void);
        // (not needed)        void                                                keepSpecialization(DagNode* affecter);
        // (not needed)        void                                                restoreSpecialization(DagNode *restorer);
        void                                                touchSpecialization(DagNode *toucher);
        void                                                updateTransitionProbabilities(size_t nodeIdx, double brlen);
        
        
    private:        
        // members
        const TypedDagNode< double >*                       homogeneousClockRate;
        const TypedDagNode< std::vector< double > >*        heterogeneousClockRates;
        const TypedDagNode< RateMatrix >*                   homogeneousRateMatrix;
        const TypedDagNode< RbVector< RateMatrix > >*       heterogeneousRateMatrices;
        const TypedDagNode< std::vector< double > >*        rootFrequencies;
        const TypedDagNode< std::vector< double > >*        siteRates;
        
        
        // flags specifying which model variants we use
        bool                                                branchHeterogeneousClockRates;
        bool                                                branchHeterogeneousSubstitutionMatrices;
        bool                                                rateVariationAcrossSites;
    };
    
}


#include "ConstantNode.h"
#include "DiscreteCharacterState.h"
#include "RateMatrix_JC.h"
#include "RandomNumberFactory.h"
#include "TopologyNode.h"
#include "TransitionProbabilityMatrix.h"

#include <cmath>
#include <cstring>
#if defined (SSE_ENABLED)
#include <xmmintrin.h>
#include <emmintrin.h>
#include <pmmintrin.h>
#elif defined ( AVX_ENABLED )
#include <xmmintrin.h>
#include <emmintrin.h>
#include <pmmintrin.h>
#include <immintrin.h>
#endif

template<class charType, class treeType>
RevBayesCore::NucleotideBranchHeterogeneousCharEvoModel<charType, treeType>::NucleotideBranchHeterogeneousCharEvoModel(const TypedDagNode<treeType> *t, bool c, size_t nSites) : AbstractSiteHomogeneousMixtureCharEvoModel<charType, treeType>(  t, 4, 1, c, nSites ) {
    
    // initialize with default parameters
    homogeneousClockRate        = new ConstantNode<double>("clockRate", new double(1.0) );
    heterogeneousClockRates     = NULL;
    homogeneousRateMatrix       = new ConstantNode<RateMatrix>("rateMatrix", new RateMatrix_JC( 4 ) );
    heterogeneousRateMatrices   = NULL;
    rootFrequencies             = NULL;
    siteRates                   = NULL;
    
    
    // flags specifying which model variants we use
    branchHeterogeneousClockRates               = false;
    branchHeterogeneousSubstitutionMatrices     = false;
    rateVariationAcrossSites                    = false;
    
    // add the parameters to the parents list
    this->addParameter( homogeneousClockRate );
    this->addParameter( homogeneousRateMatrix );
    
    this->redrawValue();
    
}


template<class charType, class treeType>
RevBayesCore::NucleotideBranchHeterogeneousCharEvoModel<charType, treeType>::NucleotideBranchHeterogeneousCharEvoModel(const NucleotideBranchHeterogeneousCharEvoModel &d) : AbstractSiteHomogeneousMixtureCharEvoModel<charType, treeType>( d ) {
    // parameters are automatically copied
    // initialize with default parameters
    homogeneousClockRate        = d.homogeneousClockRate;
    heterogeneousClockRates     = d.heterogeneousClockRates;
    homogeneousRateMatrix       = d.homogeneousRateMatrix;
    heterogeneousRateMatrices   = d.heterogeneousRateMatrices;
    rootFrequencies             = d.rootFrequencies;
    siteRates                   = d.siteRates;
    
    
    // flags specifying which model variants we use
    branchHeterogeneousClockRates               = d.branchHeterogeneousClockRates;
    branchHeterogeneousSubstitutionMatrices     = d.branchHeterogeneousSubstitutionMatrices;
    rateVariationAcrossSites                    = d.rateVariationAcrossSites;
    
}


template<class charType, class treeType>
RevBayesCore::NucleotideBranchHeterogeneousCharEvoModel<charType, treeType>::~NucleotideBranchHeterogeneousCharEvoModel( void ) {
    // We don't delete the parameters, because they might be used somewhere else too. The model needs to do that!
    
}


template<class charType, class treeType>
RevBayesCore::NucleotideBranchHeterogeneousCharEvoModel<charType, treeType>* RevBayesCore::NucleotideBranchHeterogeneousCharEvoModel<charType, treeType>::clone( void ) const {
    
    return new NucleotideBranchHeterogeneousCharEvoModel<charType, treeType>( *this );
}



template<class charType, class treeType>
void RevBayesCore::NucleotideBranchHeterogeneousCharEvoModel<charType, treeType>::computeRootLikelihood( size_t root, size_t left, size_t right) 
{
    
    // reset the likelihood
    this->lnProb = 0.0;
    
    // get the root frequencies
    const std::vector<double> &f                    = getRootFrequencies();
    
//    
//#   if defined ( SSE_ENABLED )
//    
//    __m128d f01 = _mm_set_pd(f[0],f[1]);
//    __m128d f23 = _mm_set_pd(f[2],f[3]);
//    double *tmp = new double [2];
//    
//#   elif defined ( AVX_ENABLED )
//    
//#   else
    
    double  f0  = f[0];
    double  f1  = f[1];
    double  f2  = f[2];
    double  f3  = f[3];
    
//#   endif
    
    // get the pointers to the partial likelihoods of the left and right subtree
    const double* p_left   = this->partialLikelihoods + this->activeLikelihood[left]*this->activeLikelihoodOffset + left*this->nodeOffset;
    const double* p_right  = this->partialLikelihoods + this->activeLikelihood[right]*this->activeLikelihoodOffset + right*this->nodeOffset;
    
    // create a vector for the per mixture likelihoods
    // we need this vector to sum over the different mixture likelihoods
    std::vector<double> per_mixture_Likelihoods = std::vector<double>(this->numPatterns,0.0);
    
    // get pointers the likelihood for both subtrees
    const double*   p_mixture_left     = p_left;
    const double*   p_mixture_right    = p_right;
    // iterate over all mixture categories
    for (size_t mixture = 0; mixture < this->numSiteRates; ++mixture) 
    {
        
        // get pointers to the likelihood for this mixture category
        const double*   p_site_mixture_left     = p_mixture_left;
        const double*   p_site_mixture_right    = p_mixture_right;
        // iterate over all sites
        for (size_t site = 0; site < this->numPatterns; ++site)
        {
            
//#           if defined (SSE_ENABLED)
//            
//            __m128d a01 = _mm_load_pd(p_site_mixture_left);
//            __m128d a23 = _mm_load_pd(p_site_mixture_left+2);
//            
//            __m128d b01 = _mm_load_pd(p_site_mixture_right);
//            __m128d b23 = _mm_load_pd(p_site_mixture_right+2);
//            
//            __m128d p01 = _mm_mul_pd(a01,b01);
//            __m128d p23 = _mm_mul_pd(a23,b23);
//            
//            __m128d pf01 = _mm_mul_pd(p01,f01);
//            __m128d pf23 = _mm_mul_pd(p23,f23);
//            
//            
//            __m128d sum_0123 = _mm_hadd_pd(pf01,pf23);
//            _mm_store_pd(tmp,sum_0123);
//            per_mixture_Likelihoods[site] += tmp[0] + tmp[1];
//            
//#           elif defined( AVX_ENABLED )
//            
//            
//#           else
            
            double tmp = p_site_mixture_left[0] * p_site_mixture_right[0] * f0;
            tmp += p_site_mixture_left[1] * p_site_mixture_right[1] * f1;
            tmp += p_site_mixture_left[2] * p_site_mixture_right[2] * f2;
            tmp += p_site_mixture_left[3] * p_site_mixture_right[3] * f3;
            
            // add the likelihood for this mixture category
            per_mixture_Likelihoods[site] += tmp;
            
//#           endif
            
            // increment the pointers to the next site
            p_site_mixture_left+=this->siteOffset; p_site_mixture_right+=this->siteOffset;
            
        } // end-for over all sites (=patterns)
        
        // increment the pointers to the next mixture category
        p_mixture_left+=this->mixtureOffset; p_mixture_right+=this->mixtureOffset;
        
    } // end-for over all mixtures (=rate categories)
    
    // sum the log-likelihoods for all sites together
    std::vector< size_t >::const_iterator patterns = this->patternCounts.begin();
    for (size_t site = 0; site < this->numPatterns; ++site, ++patterns)
    {
        this->lnProb += log( per_mixture_Likelihoods[site] ) * *patterns;
    }
    // normalize the log-probability
    this->lnProb -= log( this->numSiteRates ) * this->numSites;
    
//#   if defined (SSE_ENABLED)
//    delete [] tmp;
//#   endif
    
}


template<class charType, class treeType>
void RevBayesCore::NucleotideBranchHeterogeneousCharEvoModel<charType, treeType>::computeInternalNodeLikelihood(const TopologyNode &node, size_t nodeIndex, size_t left, size_t right) 
{   
    
    // compute the transition probability matrix
    updateTransitionProbabilities( nodeIndex, node.getBranchLength() );
    
#   if defined ( SSE_ENABLED )
    
    double* p_left   = this->partialLikelihoods + this->activeLikelihood[left]*this->activeLikelihoodOffset + left*this->nodeOffset;
    double* p_right  = this->partialLikelihoods + this->activeLikelihood[right]*this->activeLikelihoodOffset + right*this->nodeOffset;
    double* p_node   = this->partialLikelihoods + this->activeLikelihood[nodeIndex]*this->activeLikelihoodOffset + nodeIndex*this->nodeOffset;
    //    __m128d* p_left   = (__m128d *) this->partialLikelihoods + this->activeLikelihood[left]*this->activeLikelihoodOffset + left*this->nodeOffset;
    //    __m128d* p_right  = (__m128d *) this->partialLikelihoods + this->activeLikelihood[right]*this->activeLikelihoodOffset + right*this->nodeOffset;
    //    __m128d* p_node   = (__m128d *) this->partialLikelihoods + this->activeLikelihood[nodeIndex]*this->activeLikelihoodOffset + nodeIndex*this->nodeOffset;

#   elif defined ( AVX_ENABLED )

    double* p_left   = this->partialLikelihoods + this->activeLikelihood[left]*this->activeLikelihoodOffset + left*this->nodeOffset;
    double* p_right  = this->partialLikelihoods + this->activeLikelihood[right]*this->activeLikelihoodOffset + right*this->nodeOffset;
    double* p_node   = this->partialLikelihoods + this->activeLikelihood[nodeIndex]*this->activeLikelihoodOffset + nodeIndex*this->nodeOffset;

#   else

    // get the pointers to the partial likelihoods for this node and the two descendant subtrees
    const double*   p_left  = this->partialLikelihoods + this->activeLikelihood[left]*this->activeLikelihoodOffset + left*this->nodeOffset;
    const double*   p_right = this->partialLikelihoods + this->activeLikelihood[right]*this->activeLikelihoodOffset + right*this->nodeOffset;
    double*         p_node  = this->partialLikelihoods + this->activeLikelihood[nodeIndex]*this->activeLikelihoodOffset + nodeIndex*this->nodeOffset;

#   endif
    
    // iterate over all mixture categories
    for (size_t mixture = 0; mixture < this->numSiteRates; ++mixture)
    {
        // the transition probability matrix for this mixture category
        const double* tp_begin = this->transitionProbMatrices[mixture].theMatrix;
        
        // get the pointers to the likelihood for this mixture category
        size_t offset = mixture*this->mixtureOffset;
        
#       if defined ( SSE_ENABLED )
        
        double*          p_site_mixture          = p_node + offset;
        const double*    p_site_mixture_left     = p_left + offset;
        const double*    p_site_mixture_right    = p_right + offset;
        
        __m128d tp_a_ac = _mm_load_pd(tp_begin);
        __m128d tp_a_gt = _mm_load_pd(tp_begin+2);
        __m128d tp_c_ac = _mm_load_pd(tp_begin+4);
        __m128d tp_c_gt = _mm_load_pd(tp_begin+6);
        __m128d tp_g_ac = _mm_load_pd(tp_begin+8);
        __m128d tp_g_gt = _mm_load_pd(tp_begin+10);
        __m128d tp_t_ac = _mm_load_pd(tp_begin+12);
        __m128d tp_t_gt = _mm_load_pd(tp_begin+14);
        
#       elif defined ( AVX_ENABLED )
        
        double*          p_site_mixture          = p_node + offset;
        const double*    p_site_mixture_left     = p_left + offset;
        const double*    p_site_mixture_right    = p_right + offset;
        
        __m256d tp_a = _mm256_load_pd(tp_begin);
        __m256d tp_c = _mm256_load_pd(tp_begin+4);
        __m256d tp_g = _mm256_load_pd(tp_begin+8);
        __m256d tp_t = _mm256_load_pd(tp_begin+12);
        
#       else

        double*          p_site_mixture          = p_node + offset;
        const double*    p_site_mixture_left     = p_left + offset;
        const double*    p_site_mixture_right    = p_right + offset;

#       endif

        // compute the per site probabilities
        for (size_t site = 0; site < this->numPatterns ; ++site)
        {
            
#           if defined ( SSE_ENABLED )
            
            __m128d a01 = _mm_load_pd(p_site_mixture_left);
            __m128d a23 = _mm_load_pd(p_site_mixture_left+2);
            
            __m128d b01 = _mm_load_pd(p_site_mixture_right);
            __m128d b23 = _mm_load_pd(p_site_mixture_right+2);
            
            __m128d p01 = _mm_mul_pd(a01,b01);
            __m128d p23 = _mm_mul_pd(a23,b23);
            
            __m128d a_ac = _mm_mul_pd(p01, tp_a_ac   );
            __m128d a_gt = _mm_mul_pd(p23, tp_a_gt );
            __m128d a_acgt = _mm_hadd_pd(a_ac,a_gt);
            
            __m128d c_ac = _mm_mul_pd(p01, tp_c_ac );
            __m128d c_gt = _mm_mul_pd(p23, tp_c_gt );
            __m128d c_acgt = _mm_hadd_pd(c_ac,c_gt);
            
            //            *p_site_mixture = _mm_hadd_pd(a_acgt,c_acgt);
            __m128d ac = _mm_hadd_pd(a_acgt,c_acgt);
            _mm_store_pd(p_site_mixture,ac);
            
            
            __m128d g_ac = _mm_mul_pd(p01, tp_g_ac  );
            __m128d g_gt = _mm_mul_pd(p23, tp_g_gt );
            __m128d g_acgt = _mm_hadd_pd(g_ac,g_gt);
            
            __m128d t_ac = _mm_mul_pd(p01, tp_t_ac );
            __m128d t_gt = _mm_mul_pd(p23, tp_t_gt );
            __m128d t_acgt = _mm_hadd_pd(t_ac,t_gt);
            
            //            p_site_mixture[2] = _mm_hadd_pd(g_acgt,t_acgt);
            __m128d gt = _mm_hadd_pd(g_acgt,t_acgt);
            _mm_store_pd(p_site_mixture+2,gt);
 
#           elif defined ( AVX_ENABLED )
 
            __m256d a = _mm256_load_pd(p_site_mixture_left);
            __m256d b = _mm256_load_pd(p_site_mixture_right);
            __m256d p = _mm_mul_pd(a,b);
            
            __m256d a_acgt = _mm256_mul_pd(p, tp_a );
            __m256d c_acgt = _mm256_mul_pd(p, tp_c );
            __m256d g_acgt = _mm256_mul_pd(p, tp_g );
            __m256d t_acgt = _mm256_mul_pd(p, tp_t );
            
            __m256d ac   = _mm256_hadd_pd(a_acgt,c_acgt);
            __m256d gt   = _mm256_hadd_pd(g_acgt,t_acgt)
            
            __m256d acgt = _mm256_hadd_pd(ac,gt);

            _mm256_store_pd(p_site_mixture,acgt);

#           else

            double p0 = p_site_mixture_left[0] * p_site_mixture_right[0];
            double p1 = p_site_mixture_left[1] * p_site_mixture_right[1];
            double p2 = p_site_mixture_left[2] * p_site_mixture_right[2];
            double p3 = p_site_mixture_left[3] * p_site_mixture_right[3];
            
            double sum = p0 * tp_begin[0];
            sum += p1 * tp_begin[1];
            sum += p2 * tp_begin[2];
            sum += p3 * tp_begin[3];
            
            p_site_mixture[0] = sum;
            
            sum = p0 * tp_begin[4];
            sum += p1 * tp_begin[5];
            sum += p2 * tp_begin[6];
            sum += p3 * tp_begin[7];
            
            p_site_mixture[1] = sum;
            
            sum = p0 * tp_begin[8];
            sum += p1 * tp_begin[9];
            sum += p2 * tp_begin[10];
            sum += p3 * tp_begin[11];
            
            p_site_mixture[2] = sum;
            
            sum = p0 * tp_begin[12];
            sum += p1 * tp_begin[13];
            sum += p2 * tp_begin[14];
            sum += p3 * tp_begin[15];
            
            p_site_mixture[3] = sum;

#           endif
            
            // increment the pointers to the next site
            p_site_mixture_left+=this->siteOffset; p_site_mixture_right+=this->siteOffset; p_site_mixture+=this->siteOffset;

                        
        } // end-for over all sites (=patterns)
        
    } // end-for over all mixtures (=rate-categories)
    
}


template<class charType, class treeType>
void RevBayesCore::NucleotideBranchHeterogeneousCharEvoModel<charType, treeType>::computeTipLikelihood(const TopologyNode &node, size_t nodeIndex) 
{    
    
    double* p_node = this->partialLikelihoods + this->activeLikelihood[nodeIndex]*this->activeLikelihoodOffset + nodeIndex*this->nodeOffset;
    
    const std::vector<bool> &gap_node = this->gapMatrix[nodeIndex];
    const std::vector<unsigned long> &char_node = this->charMatrix[nodeIndex];
    
    // compute the transition probabilities
    updateTransitionProbabilities( nodeIndex, node.getBranchLength() );
    
    double*   p_mixture      = p_node;
    
    // iterate over all mixture categories
    for (size_t mixture = 0; mixture < this->numSiteRates; ++mixture)
    {
        // the transition probability matrix for this mixture category
        const double*                       tp_begin    = this->transitionProbMatrices[mixture].theMatrix;
        
        // get the pointer to the likelihoods for this site and mixture category
        double*     p_site_mixture      = p_mixture;
        
        // iterate over all sites
        for (size_t site = 0; site != this->numPatterns; ++site) 
        {
            
            // is this site a gap?
            if ( gap_node[site] ) 
            {
                // since this is a gap we need to assume that the actual state could have been any state
                p_site_mixture[0] = 1.0;
                p_site_mixture[1] = 1.0;
                p_site_mixture[2] = 1.0;
                p_site_mixture[3] = 1.0;
                
            } 
            else // we have observed a character
            {
                
                // get the original character
                unsigned long org_val = char_node[site];
                                    
                if ( this->usingAmbiguousCharacters )
                {                
                    
                    double p0 = 0.0;
                    double p1 = 0.0;
                    double p2 = 0.0;
                    double p3 = 0.0;
                    
                    if ( (org_val & 1) == 1 )
                    {
                        p0 = tp_begin[0];
                        p1 = tp_begin[4];
                        p2 = tp_begin[8];
                        p3 = tp_begin[12];
                    }
                    
                    if ( (org_val & 2) == 2 )
                    {
                        p0 += tp_begin[1];
                        p1 += tp_begin[5];
                        p2 += tp_begin[9];
                        p3 += tp_begin[13];
                    }
                    
                    if ( (org_val & 4) == 4 )
                    {
                        p0 += tp_begin[2];
                        p1 += tp_begin[6];
                        p2 += tp_begin[10];
                        p3 += tp_begin[14];
                    }
                    
                    if ( (org_val & 8) == 8 )
                    {
                        p0 += tp_begin[3];
                        p1 += tp_begin[7];
                        p2 += tp_begin[11];
                        p3 += tp_begin[15];
                    }
                    
                    p_site_mixture[0] = p0;
                    p_site_mixture[1] = p1;
                    p_site_mixture[2] = p2;
                    p_site_mixture[3] = p3;
                    
                } 
                else // no ambiguous characters in use
                {
                        
                    // store the likelihood
                    p_site_mixture[0] = tp_begin[org_val];
                    p_site_mixture[1] = tp_begin[4+org_val];
                    p_site_mixture[2] = tp_begin[8+org_val];
                    p_site_mixture[3] = tp_begin[12+org_val];
                        
                }
                
            } // end-if a gap state
            
            // increment the pointers to next site
            p_site_mixture+=this->siteOffset; 
            
        } // end-for over all sites/patterns in the sequence
        
        // increment the pointers to next mixture category
        p_mixture+=this->mixtureOffset;
        
    } // end-for over all mixture categories
    
}



template<class charType, class treeType>
const std::vector<double>& RevBayesCore::NucleotideBranchHeterogeneousCharEvoModel<charType, treeType>::getRootFrequencies( void ) {
    
    if ( branchHeterogeneousSubstitutionMatrices || rootFrequencies != NULL ) 
    {
        return rootFrequencies->getValue();
    } 
    else 
    {
        return homogeneousRateMatrix->getValue().getStationaryFrequencies();
    }
    
}


template<class charType, class treeType>
void RevBayesCore::NucleotideBranchHeterogeneousCharEvoModel<charType, treeType>::updateTransitionProbabilities(size_t nodeIdx, double brlen) {
    
    // first, get the rate matrix for this branch
    const RateMatrix *rm;
    if ( this->branchHeterogeneousSubstitutionMatrices == true ) 
    {
        rm = &this->heterogeneousRateMatrices->getValue()[nodeIdx];
    } 
    else 
    {
        rm = &this->homogeneousRateMatrix->getValue();
    }
    
    // second, get the clock rate for the branch
    double branchTime;
    if ( this->branchHeterogeneousClockRates == true ) 
    {
        branchTime = this->heterogeneousClockRates->getValue()[nodeIdx] * brlen;
    } 
    else 
    {
        branchTime = this->homogeneousClockRate->getValue() * brlen;
    }
    
    // and finally compute the per site rate transition probability matrix
    if ( this->rateVariationAcrossSites == true ) 
    {
        const std::vector<double> &r = this->siteRates->getValue();
        for (size_t i = 0; i < this->numSiteRates; ++i)
        {
            rm->calculateTransitionProbabilities( branchTime * r[i], this->transitionProbMatrices[i] );
        }
    } 
    else 
    {
        rm->calculateTransitionProbabilities( branchTime, this->transitionProbMatrices[0] );
    }
    
}



template<class charType, class treeType>
void RevBayesCore::NucleotideBranchHeterogeneousCharEvoModel<charType, treeType>::setClockRate(const TypedDagNode< double > *r) {
    
    // remove the old parameter first
    if ( homogeneousClockRate != NULL ) 
    {
        this->removeParameter( homogeneousClockRate );
        homogeneousClockRate = NULL;
    }
    else // heterogeneousClockRate != NULL
    {
        this->removeParameter( heterogeneousClockRates );
        heterogeneousClockRates = NULL;
    }
    
    // set the value
    branchHeterogeneousClockRates = false;
    homogeneousClockRate = r;
    
    // add the parameter
    this->addParameter( homogeneousClockRate );
    
    // redraw the current value
    if ( this->dagNode == NULL || !this->dagNode->isClamped() ) 
    {
        this->redrawValue();
    }
    
}



template<class charType, class treeType>
void RevBayesCore::NucleotideBranchHeterogeneousCharEvoModel<charType, treeType>::setClockRate(const TypedDagNode< std::vector< double > > *r) {
    
    // remove the old parameter first
    if ( homogeneousClockRate != NULL ) 
    {
        this->removeParameter( homogeneousClockRate );
        homogeneousClockRate = NULL;
    }
    else // heterogeneousClockRate != NULL
    {
        this->removeParameter( heterogeneousClockRates );
        heterogeneousClockRates = NULL;
    }
    
    // set the value
    branchHeterogeneousClockRates = true;
    heterogeneousClockRates = r;
    
    // add the parameter
    this->addParameter( heterogeneousClockRates );
    
    // redraw the current value
    if ( this->dagNode == NULL || !this->dagNode->isClamped() ) 
    {
        this->redrawValue();
    }
    
}


template<class charType, class treeType>
void RevBayesCore::NucleotideBranchHeterogeneousCharEvoModel<charType, treeType>::setRateMatrix(const TypedDagNode< RateMatrix > *rm) {
    
    // remove the old parameter first
    if ( homogeneousRateMatrix != NULL ) 
    {
        this->removeParameter( homogeneousRateMatrix );
        homogeneousRateMatrix = NULL;
    }
    else // heterogeneousRateMatrix != NULL
    {
        this->removeParameter( heterogeneousRateMatrices );
        heterogeneousRateMatrices = NULL;
    }
    
    // set the value
    branchHeterogeneousSubstitutionMatrices = false;
    homogeneousRateMatrix = rm;
    
    // add the parameter
    this->addParameter( homogeneousRateMatrix );
    
    // redraw the current value
    if ( this->dagNode == NULL || !this->dagNode->isClamped() ) 
    {
        this->redrawValue();
    }
    
}


template<class charType, class treeType>
void RevBayesCore::NucleotideBranchHeterogeneousCharEvoModel<charType, treeType>::setRateMatrix(const TypedDagNode< RbVector< RateMatrix > > *rm) {
    
    // remove the old parameter first
    if ( homogeneousRateMatrix != NULL ) 
    {
        this->removeParameter( homogeneousRateMatrix );
        homogeneousRateMatrix = NULL;
    }
    else // heterogeneousRateMatrix != NULL
    {
        this->removeParameter( heterogeneousRateMatrices );
        heterogeneousRateMatrices = NULL;
    }
    
    // set the value
    branchHeterogeneousSubstitutionMatrices = true;
    heterogeneousRateMatrices = rm;
    
    // add the parameter
    this->addParameter( heterogeneousRateMatrices );
    
    // redraw the current value
    if ( this->dagNode == NULL || !this->dagNode->isClamped() ) 
    {
        this->redrawValue();
    }
    
}


template<class charType, class treeType>
void RevBayesCore::NucleotideBranchHeterogeneousCharEvoModel<charType, treeType>::setRootFrequencies(const TypedDagNode< std::vector< double > > *f) {
    
    // remove the old parameter first
    if ( rootFrequencies != NULL ) 
    {
        this->removeParameter( rootFrequencies );
        rootFrequencies = NULL;
    }
    
    if ( f != NULL )
    {
        // set the value
        //        branchHeterogeneousSubstitutionMatrices = true;
        rootFrequencies = f;
        
        // add the parameter
        this->addParameter( rootFrequencies );
    }
    else
    {
        branchHeterogeneousSubstitutionMatrices = false;
    }
    
    // redraw the current value
    if ( this->dagNode == NULL || !this->dagNode->isClamped() ) 
    {
        this->redrawValue();
    }
}


template<class charType, class treeType>
void RevBayesCore::NucleotideBranchHeterogeneousCharEvoModel<charType, treeType>::setSiteRates(const TypedDagNode< std::vector< double > > *r) {
    
    // remove the old parameter first
    if ( siteRates != NULL ) 
    {
        this->removeParameter( siteRates );
        siteRates = NULL;
    }
    
    if ( r != NULL ) 
    {
        // set the value
        rateVariationAcrossSites = true;
        siteRates = r;
        this->numSiteRates = r->getValue().size();
        this->resizeLikelihoodVectors();
    }
    else
    {
        // set the value
        rateVariationAcrossSites = false;
        siteRates = NULL;
        this->numSiteRates = 1;
        this->resizeLikelihoodVectors();
        
    }
    
    // add the parameter
    this->addParameter( siteRates );
    
    // redraw the current value
    if ( this->dagNode == NULL || !this->dagNode->isClamped() ) 
    {
        this->redrawValue();
    }
}




template<class charType, class treeType>
void RevBayesCore::NucleotideBranchHeterogeneousCharEvoModel<charType, treeType>::swapParameter(const DagNode *oldP, const DagNode *newP) {
    
    if (oldP == homogeneousClockRate) 
    {
        homogeneousClockRate = static_cast<const TypedDagNode< double >* >( newP );
    }
    else if (oldP == heterogeneousClockRates) 
    {
        heterogeneousClockRates = static_cast<const TypedDagNode< std::vector< double > >* >( newP );
    }
    else if (oldP == homogeneousRateMatrix) 
    {
        homogeneousRateMatrix = static_cast<const TypedDagNode< RateMatrix >* >( newP );
    }
    else if (oldP == heterogeneousRateMatrices) 
    {
        heterogeneousRateMatrices = static_cast<const TypedDagNode< RbVector< RateMatrix > >* >( newP );
    }
    else if (oldP == rootFrequencies) 
    {
        rootFrequencies = static_cast<const TypedDagNode< std::vector< double > >* >( newP );
    }
    else if (oldP == siteRates) 
    {
        siteRates = static_cast<const TypedDagNode< std::vector< double > >* >( newP );
    }
    else 
    {
        AbstractSiteHomogeneousMixtureCharEvoModel<charType, treeType>::swapParameter(oldP,newP);
    }
    
}

template<class charType, class treeType>
void RevBayesCore::NucleotideBranchHeterogeneousCharEvoModel<charType, treeType>::touchSpecialization( DagNode* affecter ) {
    
    // if the topology wasn't the culprit for the touch, then we just flag everything as dirty
    if ( affecter == heterogeneousClockRates ) 
    {
        const std::set<size_t> &indices = heterogeneousClockRates->getTouchedElementIndices();
        
        // maybe all of them have been touched or the flags haven't been set properly
        if ( indices.size() == 0 ) 
        {
            // just delegate the call
            AbstractSiteHomogeneousMixtureCharEvoModel<charType, treeType>::touchSpecialization( affecter );
        } 
        else 
        {
            const std::vector<TopologyNode *> &nodes = this->tau->getValue().getNodes();
            // flag recomputation only for the nodes
            for (std::set<size_t>::iterator it = indices.begin(); it != indices.end(); ++it) 
            {
                this->recursivelyFlagNodeDirty( *nodes[*it] );
            } 
        }
    }
    else if ( affecter == heterogeneousRateMatrices )
    {
        
        const std::set<size_t> &indices = heterogeneousRateMatrices->getTouchedElementIndices();
        
        // maybe all of them have been touched or the flags haven't been set properly
        if ( indices.size() == 0 ) 
        {
            // just delegate the call
            AbstractSiteHomogeneousMixtureCharEvoModel<charType, treeType>::touchSpecialization( affecter );
        } 
        else 
        {
            const std::vector<TopologyNode *> &nodes = this->tau->getValue().getNodes();
            // flag recomputation only for the nodes
            for (std::set<size_t>::iterator it = indices.begin(); it != indices.end(); ++it) 
            {
                this->recursivelyFlagNodeDirty( *nodes[*it] );
            } 
        }
    }
    else if ( affecter == rootFrequencies )
    {
        
        const TopologyNode &root = this->tau->getValue().getRoot();
        this->recursivelyFlagNodeDirty( root );
    }
    else
    {
        AbstractSiteHomogeneousMixtureCharEvoModel<charType, treeType>::touchSpecialization( affecter );
    }
    
}


#endif

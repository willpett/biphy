#include "BinaryDolloSubstitutionModel.h"

#include "BinaryTaxonData.h"
#include "ContinuousBinaryCharacterData.h"
#include "ContinuousBinaryTaxonData.h"
#include "DistributionGamma.h"
#include "DistributionPoisson.h"
#include "RandomNumberGenerator.h"
#include "RandomNumberFactory.h"
#include "TopologyNode.h"

#include "numerics.h"

BinaryDolloSubstitutionModel::BinaryDolloSubstitutionModel(const TypedDagNode<Tree> *t, AscertainmentBias::Coding c) : BinarySubstitutionModel(t,c)
{
    this->numCorrectionMasks = 1;
}


BinaryDolloSubstitutionModel::BinaryDolloSubstitutionModel(const BinaryDolloSubstitutionModel &n) : BinarySubstitutionModel(n),
    integrationFactors(n.integrationFactors),
    activeIntegrationOffset(n.activeIntegrationOffset),
    maskNodeObservationCounts(n.maskNodeObservationCounts)
{   
}



BinaryDolloSubstitutionModel::~BinaryDolloSubstitutionModel( void ) {
    
}

BinaryDolloSubstitutionModel* BinaryDolloSubstitutionModel::clone( void ) const
{
    
    return new BinaryDolloSubstitutionModel( *this );
}

void BinaryDolloSubstitutionModel::getIncludedSiteIndices( void )
{
    correctionMaskCounts.clear();
    maskObservationCounts.clear();
    maskNodeObservationCounts.clear();
    correctionMaskMatrix.clear();
    site2mask.clear();
    std::fill(countDistribution.begin(), countDistribution.end(), 0);

    // find the unique site patterns and compute their respective frequencies
    std::vector<TopologyNode*> nodes = tau->getValue().getNodes();
    size_t tips = tau->getValue().getNumberOfTips();

    // create a vector with the correct site indices
    // some of the sites may have been excluded
    siteIndices.clear();
    size_t siteIndex = 0;
    size_t incompatible = 0;

    std::map<std::string, size_t> maskIndices;

    if(coding != AscertainmentBias::ALL)
    {
        // if we are using a correction, add a correction mask with 0 gaps.
        // it is required when simulating data, but not necessarily used
        // in computing the likelihood (e.g. if all sites have at least one gap)

        std::string gapless(tips, ' ');
        std::vector<bool> mask(tips, false);

        maskIndices[gapless] = 0;

        correctionMaskCounts.push_back(0);
        maskNodeObservationCounts.push_back(std::vector<size_t>(nodes.size(), 1));
        maskObservationCounts.push_back(tips);
        correctionMaskMatrix.push_back(mask);
    }

    for (size_t i = 0; i < this->value->getNumberOfCharacters(); ++i)
    {
        std::pair<size_t, size_t> charCounts;
        size_t numGap = 0;

        std::string mask = "";
        std::vector<bool> maskData;
        std::vector<size_t> observationCounts(nodes.size(), 1);
        for (std::vector<TopologyNode*>::iterator it = nodes.begin(); it != nodes.end(); ++it)
        {
            if ( (*it)->isTip() )
            {
                const BinaryTaxonData* taxon = this->value->getTaxonData( (*it)->getName() );
                RealNumber c = taxon->getCharacter(siteIndex);

                bool gap = taxon->getGap(siteIndex);

                if(gap)
                {
                    if(coding != AscertainmentBias::ALL)
                        mask += "-";
                    numGap++;
                }
                else
                {
                    if(coding != AscertainmentBias::ALL)
                        mask += " ";

                    if(c == 1.0)
                        charCounts.second++;
                    else if(c == 0.0)
                        charCounts.first++;
                }

                if(coding != AscertainmentBias::ALL)
                    maskData.push_back(gap);

                observationCounts[(*it)->getIndex()] = !gap;
            }
        }

        if( !isSitePatternCompatible(charCounts, numGap) )
        {
            incompatible++;
        }
        else
        {
            countDistribution[charCounts.second]++;

            siteIndices.push_back(siteIndex);

            if(coding != AscertainmentBias::ALL)
            {
                // increase the count for this mask
                std::map<std::string, size_t>::iterator it = maskIndices.find(mask);
                if(it != maskIndices.end())
                {
                    correctionMaskCounts[it->second]++;

                    site2mask.push_back(maskIndices[mask]);
                }
                else
                {
                    maskIndices[mask] = correctionMaskCounts.size();

                    site2mask.push_back(maskIndices[mask]);

                    correctionMaskCounts.push_back(1);
                    maskObservationCounts.push_back(tips - numGap);
                    maskNodeObservationCounts.push_back(observationCounts);
                    correctionMaskMatrix.push_back(maskData);
                }
            }
        }

        siteIndex++;
    }

    // resize our datset to account for the newly excluded characters
    numSites = siteIndices.size();
    N = numSites;

    // warn if we have to exclude some incompatible characters
    if(incompatible > 0 && verbose)
        std::cerr << "Warning: There are " << incompatible << " characters incompatible with the specified coding bias. These characters will be excluded." << std::endl;

    // readjust the number of correction sites to account for masked sites
    if(coding != AscertainmentBias::ALL)
        numCorrectionMasks = correctionMaskCounts.size();
}

void BinaryDolloSubstitutionModel::resizeLikelihoodVectors( void ) {

    // we resize the partial likelihood vectors to the new dimensions
    perMixtureCorrections = RealVector(numSiteRates*numNodes, 0.0);
    perMaskCorrections    = RealVector(numCorrectionMasks, 0.0);
            
    if(coding != AscertainmentBias::ALL)
    {
        correctionRateOffset 	= numCorrectionMasks*4;
        correctionMixtureOffset = correctionRateOffset;
        correctionNodeOffset    = numSiteRates*correctionMixtureOffset;
        activeCorrectionOffset  = numNodes*correctionNodeOffset;

        correctionLikelihoods = RealVector(2*activeCorrectionOffset, 0.0);
        
        if(useScaling)
        {
            activeCorrectionScalingOffset = numNodes*numCorrectionMasks;
            perNodeCorrectionLogScalingFactors = RealVector(2*activeCorrectionScalingOffset, 0.0);
        }
        else
        {
            perNodeCorrectionLogScalingFactors.clear();
        }
    }
    
    if(numPatterns == 0)
        return;
    
    //resize the likelihood vectors
#ifdef SIMD_ENABLED
    if(!useScaling)
        siteOffset              = 3*REALS_PER_SIMD_REGISTER;
    else
        siteOffset              = 2*REALS_PER_SIMD_REGISTER;
    
    numSIMDBlocks               = size_t((numPatterns - 1)/REALS_PER_SIMD_REGISTER) + 1;
    mixtureOffset               = numSIMDBlocks*siteOffset;
    
    numAllocatedPatterns        = numSIMDBlocks*REALS_PER_SIMD_REGISTER;
#else
    siteOffset                  = 2;
    mixtureOffset               = numPatterns*siteOffset;
    
    size_t numAllocatedPatterns = numPatterns;
#endif
    rateOffset 					= mixtureOffset;
    nodeOffset                  = numSiteRates*mixtureOffset;
    activeLikelihoodOffset      = numNodes*nodeOffset;
    
    partialLikelihoods   = RealVector(2 * activeLikelihoodOffset, 0.0);
    per_site_Likelihoods = RealVector(numAllocatedPatterns, 0.0);
    
    if(useScaling)
    {
        activeScalingOffset = numNodes*numAllocatedPatterns;
        perNodeSiteLogScalingFactors = RealVector(2 * activeScalingOffset, 0.0);
    }
    else
    {
        perNodeSiteLogScalingFactors.clear();
    }
    
    // reset the transitionProbability vector
    tMixtureOffset = 2;
    tRateOffset = tMixtureOffset;
    tNodeOffset = numSiteRates*tMixtureOffset;
    tActiveOffset = (numNodes - 1)*tNodeOffset;
    
    transitionProbabilities = RealVector(2 * tActiveOffset, 0.0);
    
    integrationFactors = RealVector(2 * numNodes * numSiteRates, 0.0);
    activeIntegrationOffset = numNodes*numSiteRates;
    
    // Reinitialize the tip data
    for(size_t nodeIndex = 0; nodeIndex < tau->getValue().getNumberOfTips(); nodeIndex++)
    {
        const TopologyNode& node = tau->getValue().getNode(nodeIndex);
        // this is a tip node
        // compute the likelihood for the tip and we are done
        BinarySubstitutionModel::computeNodeLikelihood(node, nodeIndex);
        
        if(coding != AscertainmentBias::ALL)
            computeNodeCorrection(node, nodeIndex);
    }
}


void BinaryDolloSubstitutionModel::computeNodeLikelihood(const TopologyNode &node, size_t nodeIndex, size_t left, size_t right)
{   
    // get the pointers to the partial likelihoods for this node and the two descendant subtrees
    RealVector::const_iterator  p_left  = partialLikelihoods.begin() + activeLikelihood[left]*activeLikelihoodOffset      + left*nodeOffset;
    RealVector::const_iterator  p_right = partialLikelihoods.begin() + activeLikelihood[right]*activeLikelihoodOffset     + right*nodeOffset;
    RealVector::iterator        p_node  = partialLikelihoods.begin() + activeLikelihood[nodeIndex]*activeLikelihoodOffset + nodeIndex*nodeOffset;
    
    RealVector::iterator        pi_left   = transitionProbabilities.begin() + activeProbability[left]*tActiveOffset  + left*tNodeOffset;
    RealVector::iterator        pi_right  = transitionProbabilities.begin() + activeProbability[right]*tActiveOffset + right*tNodeOffset;

#ifdef SIMD_ENABLED   
    RealVector::iterator        int_left  = integrationFactors.begin() + activeProbability[left]*activeIntegrationOffset  + left*numSiteRates;  
    RealVector::iterator        int_right = integrationFactors.begin() + activeProbability[right]*activeIntegrationOffset + right*numSiteRates;                    
            
    SIMDRegister          m1, m2, m3, m4, mIntL, mIntR, zero;

        zero = SIMD_SET1 (0.0);
        
#endif
    
    // iterate over all mixture categories
    for (size_t mixture = 0; mixture < numSiteRates; ++mixture)
    {   
        size_t tOffset = mixture*tMixtureOffset;
        
        RealVector::iterator    t_left   = pi_left  + tOffset;
        RealVector::iterator    t_right  = pi_right + tOffset;
        
#ifdef SIMD_ENABLED   
        
        if(!useScaling)
        {
            mIntL = SIMD_LOAD1 (&int_left[mixture]);
            mIntR = SIMD_LOAD1 (&int_right[mixture]);
        }
        
        for (size_t site = 0; site < numSIMDBlocks ; ++site)
        {
            size_t offset = mixture*mixtureOffset + site*siteOffset;
            
            SIMDRegister *          p_site_mixture    = (SIMDRegister *)&*(p_node  + offset);
            SIMDRegister *     p_site_mixture_left    = (SIMDRegister *)&*(p_left  + offset);
            SIMDRegister *    p_site_mixture_right    = (SIMDRegister *)&*(p_right + offset);

            /* compute the ancestral map */
            p_site_mixture[0] = SIMD_MUL (p_site_mixture_left[0], p_site_mixture_right[0]);
            
            /* compute the survival probabilities */
            m1 = SIMD_LOAD1 (&t_left[0]);
            m1 = SIMD_MUL (m1, p_site_mixture_left[0]);
            
            m2 = SIMD_LOAD1 (&t_left[1]);
            m2 = SIMD_MUL (m2, p_site_mixture_left[1]);
            
            m3 = SIMD_ADD (m1, m2);
            
            m1 = SIMD_LOAD1 (&t_right[0]);
            m1 = SIMD_MUL (m1, p_site_mixture_right[0]);
            
            m2 = SIMD_LOAD1 (&t_right[1]);
            m2 = SIMD_MUL (m2, p_site_mixture_right[1]);
            
            m1 = SIMD_ADD (m1, m2);
            p_site_mixture[1] = SIMD_MUL (m1, m3);
            
            if(!useScaling)
            {
                // sum up the the integrated likelihood terms for ancestral nodes
                
                //get the term from the left node
                m1 = SIMD_MUL (mIntL, p_site_mixture_left[1]);
                m1 = SIMD_ADD (m1, p_site_mixture_left[2]);
                
                //get the term from the right node
                m2 = SIMD_MUL (mIntR, p_site_mixture_right[1]);
                m2 = SIMD_ADD (m2, p_site_mixture_right[2]);
                
                // if only one child has descendants, it can potentially be an ancestral node
                //
                // if p_site_mixture_left[0] == 0 and p_site_mixture_right[0] > 0.0  , then p_site_mixture[2] = m1
                // if p_site_mixture_left[0] > 0  and p_site_mixture_right[0] == 0.0 , then p_site_mixture[2] = m2
                // if p_site_mixture_left[0] > 0  and p_site_mixture_right[0] > 0.0  , then p_site_mixture[2] = m1 + m2
                // if p_site_mixture_left[0] == 0 and p_site_mixture_right[0] == 0.0 , then p_site_mixture[2] = 0.0
                m3 = SIMD_CMPNEQ (p_site_mixture_left[0], zero);
                m4 = SIMD_CMPNEQ (p_site_mixture_right[0], zero);
                
                m3 = SIMD_AND(m3, m2);
                m4 = SIMD_AND(m4, m1);
                
                p_site_mixture[2] = SIMD_ADD (m3, m4);
                                
            }
#else   
        // compute the per site probabilities
        for (size_t site = 0; site < numPatterns ; ++site)
        {
            size_t offset = mixture*mixtureOffset + site*siteOffset;
            
            RealVector::iterator          p_site_mixture          = p_node + offset;
            RealVector::const_iterator    p_site_mixture_left     = p_left + offset;
            RealVector::const_iterator    p_site_mixture_right    = p_right + offset;
            
            p_site_mixture[0] = ( p_site_mixture_left[0]  * p_site_mixture_right[0] );
            
            p_site_mixture[1] = ( p_site_mixture_left[0]  * t_left[0]  + p_site_mixture_left[1]  * t_left[1] )
                              * ( p_site_mixture_right[0] * t_right[0] + p_site_mixture_right[1] * t_right[1]);
            
            /*//if only one child has descendants, it can potentially be an ancestral node
            if(p_site_mixture_left[0] != p_site_mixture_right[0])
            {
                // only one of these terms will be non-zero
                p_site_mixture[2] = p_site_mixture_left[0]  * (p_site_mixture_right[2] + p_site_mixture_right[1] * rightIntegrationFactor)
                                  + p_site_mixture_right[0] * (p_site_mixture_left[2]  + p_site_mixture_left[1]  * leftIntegrationFactor);
            }
            else
            {
                p_site_mixture[2] = 0.0;
            }*/
#endif            
        } // end-for over all sites (=patterns)
        
    } // end-for over all mixtures (=rate-categories)
    
}



void BinaryDolloSubstitutionModel::computeNodeLikelihood(const TopologyNode &node, size_t nodeIndex, size_t left, size_t right, size_t middle)
{
    // get the pointers to the partial likelihoods for this node and the two descendant subtrees
    RealVector::const_iterator   p_left   = partialLikelihoods.begin() + activeLikelihood[left]*activeLikelihoodOffset + left*nodeOffset;
    RealVector::const_iterator   p_right  = partialLikelihoods.begin() + activeLikelihood[right]*activeLikelihoodOffset + right*nodeOffset;
    RealVector::const_iterator   p_middle = partialLikelihoods.begin() + activeLikelihood[middle]*activeLikelihoodOffset + middle*nodeOffset;
    RealVector::iterator         p_node   = partialLikelihoods.begin() + activeLikelihood[nodeIndex]*activeLikelihoodOffset + nodeIndex*nodeOffset;
    
    RealVector::iterator    pi_left   = transitionProbabilities.begin() + activeProbability[left]*tActiveOffset + left*tNodeOffset;
    RealVector::iterator    pi_right  = transitionProbabilities.begin() + activeProbability[right]*tActiveOffset + right*tNodeOffset;
    RealVector::iterator    pi_middle = transitionProbabilities.begin() + activeProbability[middle]*tActiveOffset + middle*tNodeOffset;
    
#ifdef SIMD_ENABLED
    RealVector::iterator        int_left  = integrationFactors.begin() + activeProbability[left]*activeIntegrationOffset   + left*numSiteRates;  
    RealVector::iterator        int_right = integrationFactors.begin() + activeProbability[right]*activeIntegrationOffset  + right*numSiteRates; 
    RealVector::iterator       int_middle = integrationFactors.begin() + activeProbability[middle]*activeIntegrationOffset + middle*numSiteRates; 
        
    SIMDRegister          m1, m2, m3, m4, m5, mIntL, mIntR, mIntM, zero;

    zero = SIMD_SET1 (0.0);
        
#endif    
    
    // iterate over all mixture categories
    for (size_t mixture = 0; mixture < numSiteRates; ++mixture)
    {   
        size_t tOffset = mixture*tMixtureOffset;
                
        RealVector::iterator    t_left   = pi_left   + tOffset;
        RealVector::iterator    t_right  = pi_right  + tOffset;
        RealVector::iterator    t_middle = pi_middle + tOffset;
        
#ifdef SIMD_ENABLED
        
        if(!useScaling)
        {
            mIntL = SIMD_LOAD1 (&int_left[mixture]);
            mIntR = SIMD_LOAD1 (&int_right[mixture]);
            mIntM = SIMD_LOAD1 (&int_middle[mixture]);
        }
        
        for (size_t site = 0; site < numSIMDBlocks ; ++site)
        {
            size_t offset = mixture*mixtureOffset + site*siteOffset;
            
            SIMDRegister *          p_site_mixture    = (SIMDRegister *)&*(p_node  + offset);
            SIMDRegister *     p_site_mixture_left    = (SIMDRegister *)&*(p_left  + offset);
            SIMDRegister *    p_site_mixture_right    = (SIMDRegister *)&*(p_right + offset);
            SIMDRegister *   p_site_mixture_middle    = (SIMDRegister *)&*(p_middle + offset);

            /*compute the ancestral map */
            p_site_mixture[0] = SIMD_MUL (p_site_mixture_left[0], p_site_mixture_right[0]);
            p_site_mixture[0] = SIMD_MUL (p_site_mixture[0], p_site_mixture_middle[0]);
            
            /*compute the survival probabilities */
            m1 = SIMD_LOAD1 (&t_left[0]);
            m1 = SIMD_MUL (m1, p_site_mixture_left[0]);
            
            m2 = SIMD_LOAD1 (&t_left[1]);
            m2 = SIMD_MUL (m2, p_site_mixture_left[1]);
            
            m3 = SIMD_ADD (m1, m2);
            
            
            m1 = SIMD_LOAD1 (&t_right[0]);
            m1 = SIMD_MUL (m1, p_site_mixture_right[0]);
            
            m2 = SIMD_LOAD1 (&t_right[1]);
            m2 = SIMD_MUL (m2, p_site_mixture_right[1]);
            
            m1 = SIMD_ADD (m1, m2);
            m3 = SIMD_MUL (m1, m3);
            
                        
            m1 = SIMD_LOAD1 (&t_middle[0]);
            m1 = SIMD_MUL (m1, p_site_mixture_middle[0]);
            
            m2 = SIMD_LOAD1 (&t_middle[1]);
            m2 = SIMD_MUL (m2, p_site_mixture_middle[1]);
            
            m1 = SIMD_ADD (m1, m2);
            p_site_mixture[1] = SIMD_MUL (m1, m3);
            
            if(!useScaling)
            {
                /*compute the integrated likelihood */
                m1 = SIMD_MUL (mIntL, p_site_mixture_left[1]);
                m1 = SIMD_ADD (m1, p_site_mixture_left[2]);
                
                m2 = SIMD_CMPNEQ (p_site_mixture_right[0], zero);
                m3 = SIMD_CMPNEQ (p_site_mixture_middle[0], zero);
                
                m4 = SIMD_AND(m2, m3);
                m1 = SIMD_AND(m4, m1);
                
                m2 = SIMD_MUL (mIntR, p_site_mixture_right[1]);
                m2 = SIMD_ADD (m2, p_site_mixture_right[2]);
                
                m4 = SIMD_CMPNEQ (p_site_mixture_left[0], zero);
                                
                m5 = SIMD_AND(m3, m4);
                m2 = SIMD_AND(m5, m2);
                
                p_site_mixture[2] = SIMD_ADD(m1, m2);
                
                m1 = SIMD_MUL (mIntM, p_site_mixture_middle[1]);
                m1 = SIMD_ADD (m1, p_site_mixture_middle[2]);
                
                m5 = SIMD_AND(m2, m4);
                m5 = SIMD_AND(m5, m1);
                
                p_site_mixture[2] = SIMD_ADD(p_site_mixture[2], m5);
            }
#else        
                                
        // compute the per site probabilities
        for (size_t site = 0; site < numPatterns ; ++site)
        {
            size_t offset = mixture*mixtureOffset + site*siteOffset;
            
            RealVector::iterator          p_site_mixture          = p_node   + offset;
            RealVector::const_iterator    p_site_mixture_left     = p_left   + offset;
            RealVector::const_iterator    p_site_mixture_right    = p_right  + offset;
            RealVector::const_iterator    p_site_mixture_middle   = p_middle + offset;
            
            p_site_mixture[0] = ( p_site_mixture_left[0] * p_site_mixture_right[0] * p_site_mixture_middle[0]);
            
            p_site_mixture[1] = ( p_site_mixture_left[0]   * t_left[0]   + p_site_mixture_left[1]   * t_left[1]  )
                              * ( p_site_mixture_right[0]  * t_right[0]  + p_site_mixture_right[1]  * t_right[1] )
                              * ( p_site_mixture_middle[0] * t_middle[0] + p_site_mixture_middle[1] * t_middle[1]);
            
            /*  if only one child has descendants, then it can potentially be an ancestral node
            if(p_site_mixture_left[0] + p_site_mixture_right[0] + p_site_mixture_middle[0] == 2.0)
            {
                // only one of these terms will be non-zero
                p_site_mixture[2] = (p_site_mixture_right[0] * p_site_mixture_middle[0]) * (p_site_mixture_left[2]   + p_site_mixture_left[1]   * leftIntegrationFactor)
                                  + (p_site_mixture_left[0]  * p_site_mixture_middle[0]) * (p_site_mixture_right[2]  + p_site_mixture_right[1]  * rightIntegrationFactor)
                                  + (p_site_mixture_left[0]  * p_site_mixture_right[0])  * (p_site_mixture_middle[2] + p_site_mixture_middle[1] * middleIntegrationFactor);
                                  
            }
            else
            {
                p_site_mixture[2] = 0.0;
            }*/
#endif            
        } // end-for over all sites (=patterns)
        
    } // end-for over all mixtures (=rate-categories)
    
}
    




void BinaryDolloSubstitutionModel::computeNodeCorrection(const TopologyNode &node, size_t nodeIndex)
{
    for(size_t active = 0; active < 2; active++)
    {
        RealVector::iterator p_node = correctionLikelihoods.begin() + active*activeCorrectionOffset + nodeIndex*correctionNodeOffset;
    
        // iterate over all mixture categories
        for (size_t mixture = 0; mixture < numSiteRates; ++mixture)
        {
            for(size_t mask = 0; mask < numCorrectionMasks; mask++)
            {
                size_t offset = mixture*correctionMixtureOffset + mask*4;
    
                RealVector::iterator         u      = p_node + offset;
                
                bool gap = correctionMaskMatrix[mask][nodeIndex];
                
                RealVector::iterator         uC = u;
                RealVector::iterator         uI = uC + 2;
                                    
                for(size_t c = 0; c < 2; c++)
                {
                                
                    // Probability of constant state c this tip
                    // when the state at this tip is 1
                    uC[c] = (c == 1) && !gap;
                    
                    // Probability of invert singleton state c this tip
                    // when the state at this tip is 1
                    uI[c] = (c != 1) && !gap;
                }
            }
    
        }
    }
}

void BinaryDolloSubstitutionModel::computeNodeCorrection(const TopologyNode &node, size_t nodeIndex, size_t left, size_t right, size_t middle)
{
    for(size_t mask = 0; mask < numCorrectionMasks; mask++)
    {
        maskNodeObservationCounts[mask][nodeIndex] = maskNodeObservationCounts[mask][left] + maskNodeObservationCounts[mask][right] + maskNodeObservationCounts[mask][middle];
    }

    // get the pointers to the partial likelihoods for this node and the two descendant subtrees
    RealVector::const_iterator   p_left   = correctionLikelihoods.begin() + activeLikelihood[left]*activeCorrectionOffset + left*correctionNodeOffset;
    RealVector::const_iterator   p_right  = correctionLikelihoods.begin() + activeLikelihood[right]*activeCorrectionOffset + right*correctionNodeOffset;
    RealVector::const_iterator   p_middle = correctionLikelihoods.begin() + activeLikelihood[middle]*activeCorrectionOffset + middle*correctionNodeOffset;
    RealVector::iterator         p_node   = correctionLikelihoods.begin() + activeLikelihood[nodeIndex]*activeCorrectionOffset + nodeIndex*correctionNodeOffset;

    RealVector::iterator    pi_left   = transitionProbabilities.begin() + activeProbability[left]*tActiveOffset + left*tNodeOffset;
    RealVector::iterator    pi_right  = transitionProbabilities.begin() + activeProbability[right]*tActiveOffset + right*tNodeOffset;
    RealVector::iterator    pi_middle = transitionProbabilities.begin() + activeProbability[middle]*tActiveOffset + middle*tNodeOffset;
    
    // iterate over all mixture categories
    for (size_t mixture = 0; mixture < numSiteRates; ++mixture)
    {
        size_t tOffset = mixture*tMixtureOffset;
        
        RealVector::iterator    t_left   = pi_left   + tOffset;
        RealVector::iterator    t_right  = pi_right  + tOffset;
        RealVector::iterator    t_middle = pi_middle  + tOffset;
                
        for(size_t mask = 0; mask < numCorrectionMasks; mask++){

            size_t offset = mixture*correctionMixtureOffset + mask*4;
                        
            RealVector::iterator         C_i = p_node  + offset;
            RealVector::iterator         I_i = C_i + 2;
            
            RealVector::const_iterator   C_j = p_left  + offset;
            RealVector::const_iterator   I_j = C_j + 2;
            
            RealVector::const_iterator   C_k = p_right  + offset;
            RealVector::const_iterator   I_k = C_k + 2;
            
            RealVector::const_iterator   C_l = p_middle + offset;
            RealVector::const_iterator   I_l = C_l + 2;
            
            RealNumber P11j =  t_left[1];
            RealNumber P11k = t_right[1];
            RealNumber P11l = t_middle[1];
            
            RealNumber P10j =  t_left[0];
            RealNumber P10k = t_right[0];
            RealNumber P10l = t_middle[0];
                
            C_i[0] = (P10j + P11j*C_j[0]) * (P10k + P11k*C_k[0]) * (P10l + P11l*C_l[0]);
            C_i[1] = P11j*C_j[1] * P11k*C_k[1] * P11l*C_l[1];
                    
            I_i[0] = P11j*I_j[0] * (P10k + P11k*C_k[0]) * (P10l + P11l*C_l[0])
                   + (P10j + P11j*C_j[0]) * P11k*I_k[0] * (P10l + P11l*C_l[0])
                   + (P10j + P11j*C_j[0]) * (P10k + P11k*C_k[0]) * P11l*I_l[0];
            
            bool jTip = left < numTaxa; 
            bool kTip = right < numTaxa;
            bool lTip = middle < numTaxa;
            
            I_i[1] = (P11j*I_j[1] + P10j*jTip) * P11k*C_k[1] * P11l*C_l[1] 
                   + P11j*C_j[1] * (P11k*I_j[1] + P10k*kTip) * P11l*C_l[1] 
                   + P11j*C_j[1] * P11k*C_k[1] * (P11l*I_l[1] + P10l*lTip);
        }
    }
}



void BinaryDolloSubstitutionModel::computeNodeCorrection(const TopologyNode &node, size_t nodeIndex, size_t left, size_t right)
{
    for(size_t mask = 0; mask < numCorrectionMasks; mask++)
    {
        maskNodeObservationCounts[mask][nodeIndex] = maskNodeObservationCounts[mask][left] + maskNodeObservationCounts[mask][right];
    }

    // get the pointers to the partial likelihoods for this node and the two descendant subtrees
    RealVector::const_iterator   p_left  = correctionLikelihoods.begin() + activeLikelihood[left]*activeCorrectionOffset + left*correctionNodeOffset;
    RealVector::const_iterator   p_right = correctionLikelihoods.begin() + activeLikelihood[right]*activeCorrectionOffset + right*correctionNodeOffset;
    RealVector::iterator         p_node  = correctionLikelihoods.begin() + activeLikelihood[nodeIndex]*activeCorrectionOffset + nodeIndex*correctionNodeOffset;

    RealVector::iterator    pi_left   = transitionProbabilities.begin() + activeProbability[left]*tActiveOffset + left*tNodeOffset;
    RealVector::iterator    pi_right  = transitionProbabilities.begin() + activeProbability[right]*tActiveOffset + right*tNodeOffset;

    // iterate over all mixture categories
    for (size_t mixture = 0; mixture < numSiteRates; ++mixture)
    {
        size_t tOffset = mixture*tMixtureOffset;
        
        RealVector::iterator    t_left   = pi_left   + tOffset;
        RealVector::iterator    t_right  = pi_right  + tOffset;

        RealNumber P11j =  t_left[1];
        RealNumber P11k = t_right[1];
        
        RealNumber P10j =  t_left[0];
        RealNumber P10k = t_right[0];
        
        for(size_t mask = 0; mask < numCorrectionMasks; mask++){

            size_t offset = mixture*correctionMixtureOffset + mask*4;
            
            RealVector::iterator         C_i = p_node  + offset;
            RealVector::iterator         I_i = C_i + 2;
            
            RealVector::const_iterator   C_j = p_left  + offset;
            RealVector::const_iterator   I_j = C_j + 2;
            
            RealVector::const_iterator   C_k = p_right  + offset;
            RealVector::const_iterator   I_k = C_k + 2;
                
            C_i[0] = (P10j + P11j*C_j[0]) * (P10k + P11k*C_k[0]);
            C_i[1] = P11j*C_j[1] * P11k*C_k[1];
                    
            I_i[0] = P11j*I_j[0] * (P10k + P11k*C_k[0])
                   + (P10j + P11j*C_j[0]) * P11k*I_k[0];
            
            bool jTip = left < numTaxa; 
            bool kTip = right < numTaxa;
            
            I_i[1] = (P11j*I_j[1] + P10j*jTip) * P11k*C_k[1] 
                   + P11j*C_j[1] * (P11k*I_k[1] + P10k*kTip);
        }
    }
}


RealNumber BinaryDolloSubstitutionModel::sumRootLikelihood( void )
{
    RealNumber sumPartialProbs = sumUncorrectedRootLikelihood();
    
    std::fill(perMaskCorrections.begin(), perMaskCorrections.end(), 0.0);
    std::fill(perMixtureCorrections.begin(), perMixtureCorrections.end(), 0.0);

    double sampling = 1.0 - getSamplingRate();
            
    if(coding == AscertainmentBias::ALL)
    {
        for(size_t nodeIndex = 0; nodeIndex < numNodes; nodeIndex++)
        {
            RealVector::iterator integrationFactor = integrationFactors.begin() + activeProbability[nodeIndex]*activeIntegrationOffset + nodeIndex*numSiteRates;
            
            for (size_t mixture = 0; mixture < numSiteRates; ++mixture)
            {
                perMixtureCorrections[nodeIndex*numSiteRates + mixture] = integrationFactor[mixture];
                
                perMaskCorrections[0] += integrationFactor[mixture];
            }
        }
        
        perMaskCorrections[0] = log(perMaskCorrections[0]) - log(RealNumber(numSiteRates));
        
        sumPartialProbs -= perMaskCorrections[0]*numSites;
    }
    else
    {
        // iterate over each correction mask
        for(size_t mask = 0; mask < numCorrectionMasks; mask++)
        {   
            for(size_t nodeIndex = 0; nodeIndex < numNodes; nodeIndex++)
            {
                // get the root node
                const TopologyNode &node = this->tau->getValue().getNode(nodeIndex);
                
                RealVector::const_iterator p_node = correctionLikelihoods.begin() + activeLikelihood[nodeIndex] * activeCorrectionOffset  + nodeIndex*correctionNodeOffset;
                        
                RealNumber logScalingFactor = useScaling ? perNodeCorrectionLogScalingFactors[activeLikelihood[nodeIndex]*activeCorrectionScalingOffset + nodeIndex*numCorrectionMasks + mask] : 0.0;
        	RealNumber one = useScaling ? exp(-logScalingFactor) : 1.0;
         
                RealVector::iterator integrationFactor = integrationFactors.begin() + activeProbability[nodeIndex]*activeIntegrationOffset + nodeIndex*numSiteRates;
                            
                // iterate over all mixture categories
                for (size_t mixture = 0; mixture < numSiteRates; ++mixture)
                {    
                    size_t offset = mixture*correctionMixtureOffset + mask*4;
        
                    // constant site pattern likelihoods
                    RealVector::const_iterator         uC_i = p_node   + offset;
                    // invert singleton likelihoods
                    RealVector::const_iterator         uI_i = uC_i + 2;
                    
                    RealNumber prob = 0.0;
                        
                    if(coding & AscertainmentBias::NOABSENCESITES)
                        prob += uC_i[0];
                    
                    if(coding & AscertainmentBias::NOSINGLETONPRESENCE)
                        prob += uI_i[0];
                    else
                    	prob += uI_i[0]*sampling;
                    
                    // if there is only one observed tip, then don't double-count singleton gains
                    if((coding & AscertainmentBias::NOPRESENCESITES) && maskNodeObservationCounts[mask][nodeIndex] == maskObservationCounts[mask] && maskObservationCounts[mask] > 1)
                        prob += uC_i[1];
                
                    // if there are only two observed tips, then don't double-count singleton gains
                    // if there is only one observed tip, then don't double-count absence sites
                    if((coding & AscertainmentBias::NOSINGLETONABSENCE) && maskObservationCounts[mask] > 2)
                    {
                        if(maskNodeObservationCounts[mask][nodeIndex] == maskObservationCounts[mask])
                            prob += uI_i[1];
                        else if(maskNodeObservationCounts[mask][nodeIndex] == maskObservationCounts[mask] - 1)
                            prob += uC_i[1];
                    }
                    
                    // impose a per-mixture boundary
                    if(prob <= 0.0 || prob > one)
                    {
                        prob = Constants::Real::nan;
                    }

                    RealNumber mixprob = integrationFactor[mixture] * (1.0 - prob / one);
                    
                    if(mask == 0)
                        perMixtureCorrections[nodeIndex*numSiteRates + mixture] = mixprob;
                    
                    perMaskCorrections[mask] += mixprob;
                }
            }
            
            // normalize the log-probability
            perMaskCorrections[mask] = log(perMaskCorrections[mask]/numSiteRates);
            
            // apply the correction for this correction mask
            sumPartialProbs -= perMaskCorrections[mask]*correctionMaskCounts[mask];
        }
    }
    
    // TODO: do we do this if there is no ascertainment bias (is N a random variable)?
    sumPartialProbs -= log(numSites);

    return sumPartialProbs;
}


RealNumber BinaryDolloSubstitutionModel::sumUncorrectedRootLikelihood( void )
{
    //reset the per site likelihood
    std::fill(per_site_Likelihoods.begin(), per_site_Likelihoods.end(), 0.0);
        
    const TopologyNode& root = tau->getValue().getRoot();
    size_t rootIndex = root.getIndex();
    
    double logSampling = log(getSamplingRate());

    double sumPartialProbs = 0.0;
    
#ifdef SIMD_ENABLED
    RealVector::iterator p_node            = partialLikelihoods.begin() + activeLikelihood[rootIndex]  * activeLikelihoodOffset  + rootIndex*nodeOffset;
    
    if(!useScaling)
    {
        RealVector::iterator integrationFactor = integrationFactors.begin() + activeProbability[rootIndex] * activeIntegrationOffset + rootIndex*numSiteRates;

        SIMDRegister * mTotals = (SIMDRegister *)&per_site_Likelihoods[0];
        SIMDRegister m1,m2,mInt;//,mask,one,zero;
    
        for (size_t mixture = 0; mixture < numSiteRates; mixture++)
        {
            mInt = SIMD_LOAD1 (&integrationFactor[mixture]); // log-scaled if useScaling
        
            for (size_t site = 0; site < numSIMDBlocks; site++)
            { 
                // the root node is always an ancestral node
                SIMDRegister * p_site_mixture = (SIMDRegister *)&*(p_node + mixture*mixtureOffset + site*siteOffset); 

                m1 = SIMD_MUL(mInt, p_site_mixture[1]);
                m1 = SIMD_ADD(m1, p_site_mixture[2]);
                
                mTotals[site] = SIMD_ADD(m1, mTotals[site]);
            }
        }

        for (size_t pattern = 0; pattern < numPatterns; pattern++)
        {
            // if we are using scaling, then per_site_Likelihoods will already be log-scaled at this point
            per_site_Likelihoods[pattern] = log( per_site_Likelihoods[pattern] / numSiteRates );   
            
            per_site_Likelihoods[pattern] += (pattern2numPresent[pattern] == 1) ? logSampling : 0.0;
	    
            sumPartialProbs += per_site_Likelihoods[pattern]*patternCounts[pattern];
        }
    }
    else
    {
        for (size_t pattern = 0; pattern < numPatterns; ++pattern)
        {
            // if we are using scaling, then we have to log-sum-exp the integrated node probs
            std::vector<RealNumber> integratedNodeProbs;
                        
            RealNumber max = getAncestralNodeWeights(root, pattern, &integratedNodeProbs);
            
            per_site_Likelihoods[pattern] = Math::log_sum_exp(integratedNodeProbs, max) - log(RealNumber(numSiteRates)); 
            
            per_site_Likelihoods[pattern] += (pattern2numPresent[pattern] == 1) ? logSampling : 0.0;
	    
            sumPartialProbs += per_site_Likelihoods[pattern]*patternCounts[pattern];
        }
    }
#else
    for (size_t pattern = 0; pattern < numPatterns; ++pattern)
    {
        // if we are using scaling, then we have to log-sum-exp the integrated node probs
        if ( useScaling )
        {
            std::vector<RealNumber> logIntegratedNodeProbs;
            
            RealNumber max = getAncestralNodeWeights(root, pattern, &logIntegratedNodeProbs);
            
            per_site_Likelihoods[pattern] = Math::log_sum_exp(logIntegratedNodeProbs, max); 
        }
        else
        {
            //otherwise, we can just add them up
            RealNumber prob = getAncestralNodeWeights(root, pattern);
            
            per_site_Likelihoods[pattern] = log(prob);
        }
        
        //normalize the log prob
        per_site_Likelihoods[pattern] -= log(RealNumber(numSiteRates));

        per_site_Likelihoods[pattern] += (pattern2numPresent[pattern] == 1) ? logSampling : 0.0;
            
        sumPartialProbs += per_site_Likelihoods[pattern]*patternCounts[pattern];
    }
#endif
    
    return sumPartialProbs;
}


RealNumber BinaryDolloSubstitutionModel::getAncestralNodeWeights(const TopologyNode &node, size_t pattern, std::vector<RealNumber>* weights, std::vector<size_t>* nodes)
{   
    //reset the per site likelihood
#ifdef SIMD_ENABLED
    size_t patternOffset = (pattern % REALS_PER_SIMD_REGISTER) + size_t(pattern/REALS_PER_SIMD_REGISTER)*siteOffset;
    size_t idx = REALS_PER_SIMD_REGISTER;
#else
    size_t patternOffset = pattern*siteOffset;
    size_t idx = 1;
#endif
    
    RealNumber patternProb = 0.0;
    
    if(useScaling)
        patternProb = -std::numeric_limits<RealNumber>::infinity();
        
    size_t nodeIndex = node.getIndex();
    
    RealVector::const_iterator p_node  = partialLikelihoods.begin() + activeLikelihood[nodeIndex] * activeLikelihoodOffset  + nodeIndex*nodeOffset + patternOffset;
    
#ifdef SIMD_ENABLED
    RealNumber logScalingFactor = useScaling ? perNodeSiteLogScalingFactors[activeLikelihood[nodeIndex]*activeScalingOffset + nodeIndex*numAllocatedPatterns + pattern] : 0.0;
#else
    RealNumber logScalingFactor = useScaling ? perNodeSiteLogScalingFactors[activeLikelihood[nodeIndex]*activeScalingOffset + nodeIndex*numPatterns + pattern] : 0.0;
#endif
    
    RealVector::iterator integrationFactor = integrationFactors.begin() + activeProbability[nodeIndex]*activeIntegrationOffset + nodeIndex*numSiteRates;
            
    //otherwise, it is an ancestral node so we add the integrated likelihood
    for (size_t mixture = 0; mixture < numSiteRates; ++mixture)
    {            
        size_t offset = mixture*mixtureOffset;
        
        // get the pointers to the partial likelihoods of the left and right subtree
        RealVector::const_iterator p_mixture        = p_node  + offset;
        
        if ( useScaling )
        {
            RealNumber prob = log(p_mixture[idx]*integrationFactor[mixture]) + logScalingFactor;
            
            patternProb = std::max(prob, patternProb);
            
            if(weights != NULL)
                weights->push_back( prob );
        }
        else
        {
            patternProb += p_mixture[idx]*integrationFactor[mixture];
        }
    }
    
    if(nodes != NULL)
        nodes->push_back(nodeIndex);
    
    if(weights != NULL && !useScaling)
        weights->push_back(patternProb);
    
    if(node.isTip())
        return patternProb;
    
    const TopologyNode& leftChild = node.getChild(0);
    size_t left = leftChild.getIndex();
    
    const TopologyNode& rightChild = node.getChild(1);
    size_t right = rightChild.getIndex();

    RealVector::const_iterator p_left  = partialLikelihoods.begin() + activeLikelihood[left] * activeLikelihoodOffset  + left*nodeOffset + patternOffset;
    RealVector::const_iterator p_right  = partialLikelihoods.begin() + activeLikelihood[right] * activeLikelihoodOffset  + right*nodeOffset + patternOffset;
    
    //if some descendents of both children are present, then this is the last ancestral node
    if(p_left[0] != p_right[0])
    {
        if(p_right[0])
        {
            if(useScaling)
                patternProb = std::max(patternProb, getAncestralNodeWeights(leftChild, pattern, weights, nodes));
            else
                patternProb += getAncestralNodeWeights(leftChild, pattern, weights, nodes);
        }
        else if(p_left[0])
        {
            if(useScaling)
                patternProb = std::max(patternProb, getAncestralNodeWeights(rightChild, pattern, weights, nodes));
            else
                patternProb += getAncestralNodeWeights(rightChild, pattern, weights, nodes);
        }
    }
    
    return patternProb;
}

void BinaryDolloSubstitutionModel::scale( size_t nodeIndex, size_t left, size_t right )
{
    RealVector::iterator p_node = partialLikelihoods.begin() + activeLikelihood[nodeIndex]*activeLikelihoodOffset + nodeIndex*nodeOffset;
    
#ifdef SIMD_ENABLED
          SIMDRegister* p_scaler        = (SIMDRegister*)&*(perNodeSiteLogScalingFactors.begin() + activeLikelihood[nodeIndex]*activeScalingOffset + nodeIndex*numAllocatedPatterns);
    const SIMDRegister* p_scaler_left   = (SIMDRegister*)&*(perNodeSiteLogScalingFactors.begin() + activeLikelihood[left]*activeScalingOffset      + left*numAllocatedPatterns);
    const SIMDRegister* p_scaler_right  = (SIMDRegister*)&*(perNodeSiteLogScalingFactors.begin() + activeLikelihood[right]*activeScalingOffset     + right*numAllocatedPatterns);
#else
               RealVector::iterator p_scaler   = perNodeSiteLogScalingFactors.begin() + activeLikelihood[nodeIndex]*activeScalingOffset + nodeIndex*numPatterns;
    RealVector::const_iterator p_scaler_left   = perNodeSiteLogScalingFactors.begin() + activeLikelihood[left]*activeScalingOffset      + left*numPatterns;
    RealVector::const_iterator p_scaler_right  = perNodeSiteLogScalingFactors.begin() + activeLikelihood[right]*activeScalingOffset     + right*numPatterns;
#endif
    
    if ( useScaling == true && nodeIndex % scalingDensity == 0 )
    {
#ifdef SIMD_ENABLED
        // iterate over all mixture categories
        for (size_t site = 0; site < numSIMDBlocks ; ++site)
        {
            // the max probability
            SIMDRegister max;
            
            // compute the per site probabilities
            for (size_t mixture = 0; mixture < numSiteRates; ++mixture)
            {
                // get the pointers to the likelihood for this mixture category
                size_t offset = mixture*mixtureOffset + site*siteOffset;

                SIMDRegister *          p_site_mixture = (SIMDRegister*)&*(p_node + offset);
                if(mixture == 0)
                    max = p_site_mixture[1];
                else
                    max = SIMD_MAX(max, p_site_mixture[1]);
            }
            
            p_scaler[site] = SIMD_ADD(p_scaler_left[site], p_scaler_right[site]);
            p_scaler[site] = SIMD_ADD(p_scaler[site], SIMD_LOG(max));

            // compute the per site probabilities
            for (size_t mixture = 0; mixture < numSiteRates; ++mixture)
            {
                // get the pointers to the likelihood for this mixture category
                size_t offset = mixture*mixtureOffset + site*siteOffset;

                SIMDRegister *          p_site_mixture = (SIMDRegister*)&*(p_node + offset);
                p_site_mixture[1] = SIMD_DIV(p_site_mixture[1], max);

            }
#else
        // iterate over all mixture categories
        for (size_t site = 0; site < numPatterns ; ++site)
        {   
            // the max probability
            RealNumber max = 0.0;

            // compute the per site probabilities
            for (size_t mixture = 0; mixture < numSiteRates; ++mixture)
            {
                // get the pointers to the likelihood for this mixture category
                size_t offset = mixture*mixtureOffset + site*siteOffset;

                RealVector::iterator          p_site_mixture          = p_node + offset;

                if ( p_site_mixture[1] > max )
                {
                    max = p_site_mixture[1];
                }

            }

            p_scaler[site] = p_scaler_left[site] + p_scaler_right[site] + log(max);


            // compute the per site probabilities
            for (size_t mixture = 0; mixture < numSiteRates; ++mixture)
            {
                // get the pointers to the likelihood for this mixture category
                size_t offset = mixture*mixtureOffset + site*siteOffset;

                RealVector::iterator          p_site_mixture          = p_node + offset;

                p_site_mixture[1] /= max;

            }
#endif

        }
    }
    else if ( useScaling == true )
    {
        // iterate over all mixture categories
#ifdef SIMD_ENABLED
        for (size_t site = 0; site < numSIMDBlocks ; ++site)
        {
            p_scaler[site] = SIMD_ADD(p_scaler_left[site], p_scaler_right[site]);
        }
#else
        for (size_t site = 0; site < numPatterns ; ++site)
        {
            p_scaler[site] = p_scaler_left[site] + p_scaler_right[site];
        }
#endif

    }
}


void BinaryDolloSubstitutionModel::scale( size_t nodeIndex, size_t left, size_t right, size_t middle )
{
    RealVector::iterator p_node = partialLikelihoods.begin() + activeLikelihood[nodeIndex]*activeLikelihoodOffset + nodeIndex*nodeOffset;
    
#ifdef SIMD_ENABLED
    size_t scaleOffset = numSIMDBlocks*REALS_PER_SIMD_REGISTER;
    
          SIMDRegister* p_scaler        = (SIMDRegister*)&*(perNodeSiteLogScalingFactors.begin() + activeLikelihood[nodeIndex]*activeScalingOffset + nodeIndex*numAllocatedPatterns);
    const SIMDRegister* p_scaler_left   = (SIMDRegister*)&*(perNodeSiteLogScalingFactors.begin() + activeLikelihood[left]*activeScalingOffset      + left*numAllocatedPatterns);
    const SIMDRegister* p_scaler_right  = (SIMDRegister*)&*(perNodeSiteLogScalingFactors.begin() + activeLikelihood[right]*activeScalingOffset     + right*numAllocatedPatterns);
    const SIMDRegister* p_scaler_middle = (SIMDRegister*)&*(perNodeSiteLogScalingFactors.begin() + activeLikelihood[middle]*activeScalingOffset    + middle*numAllocatedPatterns);
#else
               RealVector::iterator p_scaler   = perNodeSiteLogScalingFactors.begin() + activeLikelihood[nodeIndex]*activeScalingOffset + nodeIndex*numPatterns;
    RealVector::const_iterator p_scaler_left   = perNodeSiteLogScalingFactors.begin() + activeLikelihood[left]*activeScalingOffset      + left*numPatterns;
    RealVector::const_iterator p_scaler_right  = perNodeSiteLogScalingFactors.begin() + activeLikelihood[right]*activeScalingOffset     + right*numPatterns;
    RealVector::const_iterator p_scaler_middle = perNodeSiteLogScalingFactors.begin() + activeLikelihood[middle]*activeScalingOffset    + middle*numPatterns;
#endif
    
    if ( useScaling == true && nodeIndex % scalingDensity == 0 )
    {
#ifdef SIMD_ENABLED
        // iterate over all mixture categories
        for (size_t site = 0; site < numSIMDBlocks ; ++site)
        {
            // the max probability
            SIMDRegister max;
            
            // compute the per site probabilities
            for (size_t mixture = 0; mixture < numSiteRates; ++mixture)
            {
                // get the pointers to the likelihood for this mixture category
                size_t offset = mixture*mixtureOffset + site*siteOffset;

                SIMDRegister *          p_site_mixture = (SIMDRegister*)&*(p_node + offset);
                
                if(mixture == 0)
                    max = p_site_mixture[1];
                else
                    max = SIMD_MAX(max, p_site_mixture[1]);
            }
            p_scaler[site] = SIMD_ADD(p_scaler_left[site], p_scaler_right[site]);
            p_scaler[site] = SIMD_ADD(p_scaler[site], p_scaler_middle[site]);
            p_scaler[site] = SIMD_ADD(p_scaler[site], SIMD_LOG(max));

            // compute the per site probabilities
            for (size_t mixture = 0; mixture < numSiteRates; ++mixture)
            {
                // get the pointers to the likelihood for this mixture category
                size_t offset = mixture*mixtureOffset + site*siteOffset;

                SIMDRegister *          p_site_mixture = (SIMDRegister*)&*(p_node + offset);
                p_site_mixture[1] = SIMD_DIV(p_site_mixture[1], max);

            }
#else
        // iterate over all mixture categories
        for (size_t site = 0; site < numPatterns ; ++site)
        {   
            // the max probability
            RealNumber max = 0.0;

            // compute the per site probabilities
            for (size_t mixture = 0; mixture < numSiteRates; ++mixture)
            {
                // get the pointers to the likelihood for this mixture category
                size_t offset = mixture*mixtureOffset + site*siteOffset;

                RealVector::iterator          p_site_mixture          = p_node + offset;

                if ( p_site_mixture[1] > max )
                {
                    max = p_site_mixture[1];
                }

            }

            p_scaler[site] = p_scaler_left[site] + p_scaler_right[site] + p_scaler_middle[site] + log(max);


            // compute the per site probabilities
            for (size_t mixture = 0; mixture < numSiteRates; ++mixture)
            {
                // get the pointers to the likelihood for this mixture category
                size_t offset = mixture*mixtureOffset + site*siteOffset;

                RealVector::iterator          p_site_mixture          = p_node + offset;

                p_site_mixture[1] /= max;

            }
#endif

        }
    }
    else if ( useScaling == true )
    {
        // iterate over all mixture categories
#ifdef SIMD_ENABLED
        for (size_t site = 0; site < numSIMDBlocks ; ++site)
        {
            p_scaler[site] = SIMD_ADD(p_scaler_left[site], p_scaler_right[site]);
            p_scaler[site] = SIMD_ADD(p_scaler[site], p_scaler_middle[site]);
        }
#else
        for (size_t site = 0; site < numPatterns ; ++site)
        {
            p_scaler[site] = p_scaler_left[site] + p_scaler_right[site] + p_scaler_middle[site];
        }
#endif
    }
}


void BinaryDolloSubstitutionModel::scaleCorrection( size_t nodeIndex, size_t left, size_t right )
{

    RealVector::iterator p_node   = correctionLikelihoods.begin() + activeLikelihood[nodeIndex]*activeCorrectionOffset + nodeIndex*correctionNodeOffset;

               RealVector::iterator p_scaler  = perNodeCorrectionLogScalingFactors.begin() + activeLikelihood[nodeIndex]*activeCorrectionScalingOffset + nodeIndex*numCorrectionMasks;
    RealVector::const_iterator p_scaler_left  = perNodeCorrectionLogScalingFactors.begin() + activeLikelihood[left]*activeCorrectionScalingOffset      + left*numCorrectionMasks;
    RealVector::const_iterator p_scaler_right = perNodeCorrectionLogScalingFactors.begin() + activeLikelihood[right]*activeCorrectionScalingOffset     + right*numCorrectionMasks;
    
    if ( useScaling == true && nodeIndex % scalingDensity == 0 )
    {   
#ifdef SIMD_ENABLED
        for (size_t mask = 0; mask < numCorrectionMasks; ++mask)
        {
            
            SIMDRegister max;
            
            for (size_t mixture = 0; mixture < numSiteRates; ++mixture)
            {
                for (size_t c = 0; c < 4/REALS_PER_SIMD_REGISTER; ++c)
                {
                    size_t offset = mixture*correctionMixtureOffset + mask*4 + c*REALS_PER_SIMD_REGISTER;

                    SIMDRegister *   u_i  = (SIMDRegister *)&*(p_node  + offset);
                    
                    if(mixture == 0)
                        max = *u_i;
                    else
                        max = SIMD_MAX(max, *u_i);
                }
            }
            
            RealNumber maximum = 0.0;
            
            RealNumber* tmp = (RealNumber*)&max;
            
            for(size_t i = 0; i < REALS_PER_SIMD_REGISTER; i++)
                maximum = std::max(maximum, tmp[i]);

            p_scaler[mask] = p_scaler_left[mask] + p_scaler_right[mask] + log(maximum);
            
            max = SIMD_LOAD1(&maximum);
            
            // compute the per site probabilities
            for (size_t mixture = 0; mixture < numSiteRates; ++mixture)
            {
                for (size_t c = 0; c < 4/REALS_PER_SIMD_REGISTER; ++c)
                {
                    size_t offset = mixture*correctionMixtureOffset + mask*4 + c*REALS_PER_SIMD_REGISTER;

                    SIMDRegister *   u_i  = (SIMDRegister *)&*(p_node  + offset);
                    
                    *u_i = SIMD_DIV(*u_i, max);
                }
            }
        }
#else
        for (size_t mask = 0; mask < numCorrectionMasks; ++mask)
        {
            
            RealNumber max = 0.0;
            
            for (size_t mixture = 0; mixture < numSiteRates; ++mixture)
            {
                size_t offset = mixture*correctionMixtureOffset + mask*4;
                    
                RealVector::const_iterator   u_i  = p_node  + offset;
                
                for(size_t c = 0; c < 2; c++)
                {
                    RealVector::const_iterator   uC_i = u_i;
                    RealVector::const_iterator   uI_i = uC_i + 2;
    
                    max = std::max(max, std::max(uC_i[c], uI_i[c]));
                }
            }

            p_scaler[mask] = p_scaler_left[mask] + p_scaler_right[mask] + log(max);
    
            
            // compute the per site probabilities
            for (size_t mixture = 0; mixture < numSiteRates; ++mixture)
            {
                size_t offset = mixture*correctionMixtureOffset + mask*4;
                    
                RealVector::iterator   u_i  = p_node  + offset;
                
                for(size_t c = 0; c < 2; c++)
                {
                                            
                    RealVector::iterator   uC_i = u_i;
                    RealVector::iterator   uI_i = uC_i + 2;
    
                    uC_i[c] /= max;
                    uI_i[c] /= max;
                }
            }
        }
#endif
    }
    else if ( useScaling == true )
    {
        // iterate over all character states
        for(size_t mask = 0; mask < numCorrectionMasks; mask++)
        {
            p_scaler[mask] = p_scaler_left[mask] + p_scaler_right[mask];           
        }
    }
}


void BinaryDolloSubstitutionModel::scaleCorrection( size_t nodeIndex, size_t left, size_t right, size_t middle )
{
    RealVector::iterator p_node   = correctionLikelihoods.begin() + activeLikelihood[nodeIndex]*activeCorrectionOffset + nodeIndex*correctionNodeOffset;

               RealVector::iterator p_scaler  = perNodeCorrectionLogScalingFactors.begin() + activeLikelihood[nodeIndex]*activeCorrectionScalingOffset + nodeIndex*numCorrectionMasks;
    RealVector::const_iterator p_scaler_left  = perNodeCorrectionLogScalingFactors.begin() + activeLikelihood[left]*activeCorrectionScalingOffset      + left*numCorrectionMasks;
    RealVector::const_iterator p_scaler_right = perNodeCorrectionLogScalingFactors.begin() + activeLikelihood[right]*activeCorrectionScalingOffset     + right*numCorrectionMasks;
   RealVector::const_iterator p_scaler_middle = perNodeCorrectionLogScalingFactors.begin() + activeLikelihood[middle]*activeCorrectionScalingOffset    + middle*numCorrectionMasks;
    
    if ( useScaling == true && nodeIndex % scalingDensity == 0 )
    {   
#ifdef SIMD_ENABLED
        for (size_t mask = 0; mask < numCorrectionMasks; ++mask)
        {
            
            SIMDRegister max;
            
            for (size_t mixture = 0; mixture < numSiteRates; ++mixture)
            {
                for (size_t c = 0; c < 4/REALS_PER_SIMD_REGISTER; ++mixture)
                {
                    size_t offset = mixture*correctionMixtureOffset + mask*4 + c*REALS_PER_SIMD_REGISTER;

                    SIMDRegister *   u_i  = (SIMDRegister *)&*(p_node  + offset);

                    if(mixture == 0)
                        max = *u_i;
                    else
                        max = SIMD_MAX(max, *u_i);
                }
            }
            
            RealNumber maximum = 0.0;
            
            RealNumber* tmp = (RealNumber*)&max;
            
            for(size_t i = 0; i < REALS_PER_SIMD_REGISTER; i++)
                maximum = std::max(maximum, tmp[i]);

            p_scaler[mask] = p_scaler_left[mask] + p_scaler_right[mask] + p_scaler_middle[mask] + log(maximum);
            

            max = SIMD_LOAD1(&maximum);
            
            // compute the per site probabilities
            for (size_t mixture = 0; mixture < numSiteRates; ++mixture)
            {
                for (size_t c = 0; c < 4/REALS_PER_SIMD_REGISTER; ++c)
                {
                    size_t offset = mixture*correctionMixtureOffset + mask*4 + c*REALS_PER_SIMD_REGISTER;

                    SIMDRegister *   u_i  = (SIMDRegister *)&*(p_node  + offset);
                    
                    *u_i = SIMD_DIV(*u_i, max);
                }
            }
        }
#else
        for (size_t mask = 0; mask < numCorrectionMasks; ++mask)
        {
            
            RealNumber max = 0.0;
            
            for (size_t mixture = 0; mixture < numSiteRates; ++mixture)
            {
                size_t offset = mixture*correctionMixtureOffset + mask*4;
                    
                RealVector::const_iterator   u_i  = p_node  + offset;
                
                for(size_t c = 0; c < 2; c++)
                {
                    RealVector::const_iterator   uC_i = u_i;
                    RealVector::const_iterator   uI_i = uC_i + 2;
    
                    max = std::max(max, std::max(uC_i[c], uI_i[c]));
                }
            }

            p_scaler[mask] = p_scaler_left[mask] + p_scaler_right[mask] + p_scaler_middle[mask] + log(max);
    
            
            // compute the per site probabilities
            for (size_t mixture = 0; mixture < numSiteRates; ++mixture)
            {
                size_t offset = mixture*correctionMixtureOffset + mask*4;
                    
                RealVector::iterator   u_i  = p_node  + offset;
                
                for(size_t c = 0; c < 2; c++)
                {
                    RealVector::iterator   uC_i = u_i;
                    RealVector::iterator   uI_i = uC_i + 2;
    
                    uC_i[c] /= max;
                    uI_i[c] /= max;
                }
            }
        }
#endif
    }
    else if ( useScaling == true )
    {
        // iterate over all character states
        for(size_t mask = 0; mask < numCorrectionMasks; mask++)
        {
            p_scaler[mask] = p_scaler_left[mask] + p_scaler_right[mask] + p_scaler_middle[mask];           
        }
    }
}


void BinaryDolloSubstitutionModel::updateTransitionProbabilities() {

    size_t count = 0;
    for(std::set<size_t>::iterator it = touchedNodes.begin(); it != touchedNodes.end(); it++)
    {
        size_t nodeIndex = *it;
        
        RealNumber rate = getClockRate(nodeIndex);
        
        RealVector::iterator integrationFactor = integrationFactors.begin() + activeProbability[nodeIndex]*activeIntegrationOffset + nodeIndex*numSiteRates;
        
        // only compute the integration factor for the root
        if(nodeIndex == numNodes - 1)
        {
            for (size_t mixture = 0; mixture < numSiteRates; ++mixture)
            {
                // compute the integration factor for this node
                RealNumber r = 1.0;
                if(rateVariationAcrossSites)
                    r = siteRates->getValue()[mixture];
                
                integrationFactor[mixture] = 1.0/(rate * r);
            }
        }
        else
        {
        
            RealNumber brlen = getBranchLength(nodeIndex);
        
            RealVector::iterator p_node = transitionProbabilities.begin() + activeProbability[nodeIndex]*tActiveOffset + nodeIndex*tNodeOffset;
                                   
            for (size_t mixture = 0; mixture < numSiteRates; ++mixture)
            {            
                RealNumber r = 1.0;
                if(rateVariationAcrossSites)
                    r = siteRates->getValue()[mixture];
                
                RealNumber expPart = exp( - rate * brlen * r );
                
                expPart = (expPart <= 0.0 || expPart >= 1.0) ? Constants::Real::nan : expPart; 
                
		RealVector::iterator p_node_mixture = p_node + mixture*tMixtureOffset;
        
                p_node_mixture[0] = 1.0 - expPart;
                p_node_mixture[1] = expPart;
                
                // compute the integration factor for this node
                integrationFactor[mixture] = p_node_mixture[0]/(rate * r);
            }
        }
        
        count++;
    }

}


RealNumber BinaryDolloSubstitutionModel::getStationaryFrequency( size_t nodeIndex ) {
    
    throw(Exception("Stationary frequency does not apply to Dollo model"));
}


void BinaryDolloSubstitutionModel::setStationaryFrequency(const TypedDagNode< double > *r)
{

    throw(Exception("Stationary frequency does not apply to Dollo model"));

}


void BinaryDolloSubstitutionModel::setStationaryFrequency(const TypedDagNode< std::vector< double > > *f) {

    throw(Exception("Stationary frequency does not apply to Dollo model"));

}


void BinaryDolloSubstitutionModel::redrawValue( void ) {

    if(numSites == 0)
        return; 
    
    this->dagNode->touch();
    updateTransitionProbabilities();

    // delete the old value first
    delete this->value;

    // create a new character data object
    if(continuous)
        this->value = new ContinuousBinaryCharacterData();
    else
        this->value = new BinaryCharacterData();

    RandomNumberGenerator* rng = GLOBAL_RNG;
    
    std::vector< BinaryTaxonData* > taxa;

    for(size_t i = 0; i < numTaxa; i++)
	{
		if(continuous)
			taxa.push_back( new ContinuousBinaryTaxonData("") );
		else
			taxa.push_back( new DiscreteBinaryTaxonData("") );
	}

    std::vector<size_t> rateIndices;
    std::vector<size_t> birthNodes;

    if(coding != AscertainmentBias::ALL)
    {
        /*
        // first sample the raw character birth rate from the marginal posterior 
        // lambda | N ~ Gamma(N, exp(lnCorrection) )
        double lambda = Statistics::Gamma::rv(N, exp(perMaskCorrections[0]), *rng);
    
        // then sample the observed number of characters (numSites) from the likelihood:
        // numSites | lambda ~ Poisson(lambda * exp(lnCorrection) )
        numSites = Statistics::Poisson::rv( lambda * exp(perMaskCorrections[0]), *rng);
        */

        for ( size_t i = 0; i < numSites; i++ )
        {
			double u = rng->uniform01()*perMaskCorrections[0];
			double total = 0.0;

			// simulate a birth for this character
			// by sampling nodes in proportion to the per node survival probs
			size_t birthNode = 0;
			size_t rateIndex = 0;
			while(total < u)
			{
				total += perMixtureCorrections[birthNode*numSiteRates + rateIndex];

				if(total < u)
				{
					rateIndex++;
					if(rateIndex == numSiteRates)
					{
						rateIndex = 0;
						birthNode++;
					}
				}
			}

			rateIndices.push_back(rateIndex);
			birthNodes.push_back(birthNode);
        }
    }

    double sampling = getSamplingRate();
    std::fill(countDistribution.begin(), countDistribution.end(), 0);

    // then sample site-patterns using rejection sampling,
    // rejecting those that match the unobservable ones.
    for ( size_t i = 0; i < numSites; i++ )
    {
        std::vector<RealNumber> siteData(numNodes, 0.0);

        const TopologyNode& root = tau->getValue().getNode(birthNodes[i]);
        
        // recursively simulate the sequences
        std::pair<size_t, size_t> charCounts;
        simulate( root, siteData, rateIndices[i], charCounts);

        charCounts.first = numTaxa - charCounts.second;

        if( !isSitePatternCompatible(charCounts, 0) )
        {
            i--;
            continue;
        }
        else if(samplingRate != NULL && charCounts.second == 1 && rng->uniform01() >= sampling)
		{
			i--;
			continue;
		}

        countDistribution[charCounts.second]++;

        // add the taxon data to the character data
        for (size_t t = 0; t < numTaxa; ++t)
        {
            taxa[t]->addCharacter(siteData[t]);
        }
    }

    // add the taxon data to the character data
    for (size_t i = 0; i < tau->getValue().getNumberOfTips(); ++i)
    {
        taxa[i]->setTaxonName(tau->getValue().getNode(i).getName());
        this->value->addTaxonData( taxa[i] );
    }
}


void BinaryDolloSubstitutionModel::simulate( const TopologyNode &node, std::vector<RealNumber> &data, size_t rateIndex, std::pair<size_t, size_t>& charCounts) {

    // get the children of the node
    const std::vector<TopologyNode*>& children = node.getChildren();

    // get the sequence of this node
    size_t nodeIndex = node.getIndex();

    // simulate the sequence for each child
    RandomNumberGenerator* rng = GLOBAL_RNG;
    for (std::vector< TopologyNode* >::const_iterator it = children.begin(); it != children.end(); ++it)
    {
        const TopologyNode &child = *(*it);
        size_t childIndex = child.getIndex();

        RealVector::iterator pi = transitionProbabilities.begin() + activeProbability[childIndex]*tActiveOffset + childIndex*tNodeOffset + rateIndex*tMixtureOffset;
        
        RealNumber& childState = data[ childIndex ];
        
        childState = 1.0;
        
        if(child.isTip())
        {
            //random survival at tips
            if(continuous)
            {
                childState = pi[1];
            }
            else
            {
                double u = rng->uniform01();
                if(u < pi[0])
                    childState = 0.0;   
            }
            
            // count if survived at tip
            if(childState > 0.0)
                charCounts.second++;
        }
        else
        {
            // random survival
            double u = rng->uniform01();
            if(u < pi[1])
                simulate( child, data, rateIndex, charCounts);
            else
                childState = 0.0;
        }
    }

}

const std::vector< DiscreteBinaryTaxonData >& BinaryDolloSubstitutionModel::getMapping(void) {

    computeLnProbability();

    mapping = std::vector< DiscreteBinaryTaxonData >( tau->getValue().getNumberOfNodes(), DiscreteBinaryTaxonData() );

    // first, simulate the per site rates
    RandomNumberGenerator* rng = GLOBAL_RNG;
    
    const TopologyNode& root = tau->getValue().getRoot();

    for ( size_t site = 0; site < numSites; ++site )
    {
        std::vector<RealNumber> nodeWeights;
        std::vector<size_t> ancestralNodes;
        
        RealNumber total = 0.0;
        
        std::vector<RealNumber>::iterator weight;
        std::vector<size_t>::iterator birthNode;
                
        // simulate a birth for this character
        // by sampling ancestral nodes in proportion to the per node survival probs       
        if(useScaling)
        {
            RealNumber max = getAncestralNodeWeights(root, site2pattern[site], &nodeWeights, &ancestralNodes);
            
            weight    = nodeWeights.begin();
            birthNode = ancestralNodes.begin();
                    
            total = Math::log_sum_exp(nodeWeights, max);
            
            
            double u = rng->uniform01();
            
            double tmp = 0.0;
            
            for(size_t mixture =0; mixture < numSiteRates; mixture++)
            {
                tmp += exp(*weight - total);
                weight++;
            }
            
            while(tmp <= u)
            {
                birthNode++;
                
                for(size_t mixture =0; mixture < numSiteRates; mixture++)
                {
                    tmp += exp(*weight - total);
                    weight++;
                }
            }
                    
        }
        else
        {
            total = getAncestralNodeWeights(root, site2pattern[site], &nodeWeights, &ancestralNodes);
            
            weight = nodeWeights.begin();
            birthNode = ancestralNodes.begin();
                    
            double u = rng->uniform01()*total;
            
            total = *weight;
            
            while(total <= u)
            {
                birthNode++;
                weight++;
                
                total += *weight;
            }
        }
        
        size_t nodeIndex = *birthNode;
        
#ifdef SIMD_ENABLED
        size_t patternOffset = (site2pattern[site] % REALS_PER_SIMD_REGISTER) + size_t(site2pattern[site]/REALS_PER_SIMD_REGISTER)*siteOffset;
        size_t idx = REALS_PER_SIMD_REGISTER;
#else
        size_t patternOffset = site2pattern[site]*siteOffset;
        size_t idx = 1;
#endif
        
        // get the pointers to the partial likelihoods for this node
        RealVector::iterator p_node            = partialLikelihoods.begin() + activeLikelihood[nodeIndex]*activeLikelihoodOffset  + nodeIndex*nodeOffset + patternOffset;
        RealVector::iterator integrationFactor = integrationFactors.begin() + activeProbability[nodeIndex]*activeIntegrationOffset + nodeIndex*numSiteRates;    
            
        
        // then sample a rate category conditional on survival from this node
        // by sampling in proportion to the per mixture integrated surivival probs
        total = 0.0;
        for(size_t mixture = 0; mixture < numSiteRates; mixture++)
        {
            RealVector::iterator p_node_mixture = p_node + mixture*mixtureOffset;
            total += p_node_mixture[1]*integrationFactor[mixture];
        }

        double u = rng->uniform01()*total;
        size_t rateIndex = 0;

        double tmp = 0.0;
        while(tmp < u){
            RealVector::iterator p_node_mixture = p_node + rateIndex*mixtureOffset;
            
            tmp += p_node_mixture[1]*integrationFactor[rateIndex];
            if(tmp < u)
                rateIndex++;
        }
        
        // create empty site data absent everywhere but this node for now
        std::vector<bool> siteData(numNodes, false);

        const TopologyNode& birthRoot = tau->getValue().getNode(nodeIndex);
        
        siteData[ birthRoot.getIndex() ] = true;
        
        const std::vector<TopologyNode*>& children = birthRoot.getChildren();
                
        std::pair<size_t, size_t> charCounts;
        for (std::vector< TopologyNode* >::const_iterator it = children.begin(); it != children.end(); ++it)
            simulateDolloMapping( *(*it), site2pattern[site], siteData, rateIndex, charCounts);
        
        // add the taxon data to the character data
        for (size_t nodeIndex = 0; nodeIndex < numNodes; ++nodeIndex)
        {
            mapping[nodeIndex].addCharacter(siteData[nodeIndex]);
        }
    }
        
    return mapping;
}



void BinaryDolloSubstitutionModel::simulateDolloMapping( const TopologyNode &node, size_t pattern, std::vector<bool> &data, size_t rateIndex, std::pair<size_t, size_t>& charCounts)
{
    RandomNumberGenerator* rng = GLOBAL_RNG;
    
    size_t nodeIndex = node.getIndex();
    
#ifdef SIMD_ENABLED
    size_t patternOffset = (pattern % REALS_PER_SIMD_REGISTER) + size_t(pattern/REALS_PER_SIMD_REGISTER)*siteOffset;
    size_t idx = REALS_PER_SIMD_REGISTER;
#else
    size_t patternOffset = pattern*siteOffset;
    size_t idx = 1;
#endif
    // get the pointers to the partial likelihoods for this node
    RealVector::iterator p_node  = partialLikelihoods.begin() + activeLikelihood[nodeIndex]*activeLikelihoodOffset  + nodeIndex*nodeOffset + rateIndex*mixtureOffset + patternOffset;
        
    RealVector::iterator pi = transitionProbabilities.begin() + activeProbability[nodeIndex]*tActiveOffset + nodeIndex*tNodeOffset + nodeIndex*tMixtureOffset;
    
    data[ nodeIndex ] = true;
    
    RealNumber p1 = p_node[idx]*pi[1];
    RealNumber p0 = p_node[0]*pi[0];
            
    if(p0 > 0.0)
    {
        double u = rng->uniform01()*(p0 + p1);
        if(u <= p0)
            data[ nodeIndex ] = false;
    }
    
    if(node.isTip())
    {
        if(data[ nodeIndex ])
            charCounts.second++;
        else
            charCounts.first++;
    }
    else
    {
        const std::vector<TopologyNode*>& children = node.getChildren();
        
        for (std::vector< TopologyNode* >::const_iterator it = children.begin(); it != children.end(); ++it)
            simulateDolloMapping( *(*it), pattern, data, rateIndex, charCounts);
    }
}

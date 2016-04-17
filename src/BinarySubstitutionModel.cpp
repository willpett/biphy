#include "BinarySubstitutionModel.h"

#include "BinaryTaxonData.h"
#include "ContinuousBinaryCharacterData.h"
#include "ContinuousBinaryTaxonData.h"
#include "DistributionBinomial.h"
#include "DistributionNegativeBinomial.h"
#include "RandomNumberGenerator.h"
#include "RandomNumberFactory.h"
#include "TopologyNode.h"

#include "numerics.h"

BinarySubstitutionModel::BinarySubstitutionModel(const TypedDagNode<Tree> *t, AscertainmentBias::Coding c) : TypedDistribution< BinaryCharacterData >(  new BinaryCharacterData() ),
    coding(c),
    lnProb( 0.0 ),
    numNodes( t->getValue().getNumberOfNodes() ),
    numTaxa( t->getValue().getNumberOfTips() ),
    numSites( 0 ),
    numSiteRates( 1 ),
	numSiteFrequencies( 1 ),
    tau( t ),
    activeLikelihood( std::vector<bool>(numNodes, false) ),
    activeProbability( std::vector<bool>(numNodes, false) ),
    useScaling( false ),
    numPatterns( 0 ),
    dirtyNodes( std::vector<bool>(numNodes, true) ),
    numCorrectionMasks(0),
    verbose(true),
    continuous(false),
    scalingDensity(10)
{
    for(size_t i = 0; i < numNodes - 1; i++)
    {
        touchedNodes.insert(i);
    }
    
    // initialize with default parameters
    homogeneousFrequency      = NULL;
    homogeneousClockRate      = NULL;
    heterogeneousClockRates   = NULL;
    heterogeneousFrequencies  = NULL;
    branchLengths             = NULL;
    siteRates                 = NULL;
    samplingRate			  = NULL;
        
    // add the parameters to the parents list
    this->addParameter( tau );
    tau->getValue().getTreeChangeEventHandler().addListener( this );
    
    // flags specifying which model variants we use
    branchHeterogeneousClockRates   = false;
    branchHeterogeneousFrequencies  = false;
    rateVariationAcrossSites        = false;
    frequencyVariationAcrossSites   = false;

    countDistribution = std::vector<int>(numTaxa + 1, 0);
}


BinarySubstitutionModel::BinarySubstitutionModel(const BinarySubstitutionModel &n) : TypedDistribution< BinaryCharacterData >( n ),
    coding(n.coding),
    
    siteIndices(n.siteIndices),
    patternCounts( n.patternCounts ),
    pattern2site(n.pattern2site),
	pattern2numPresent(n.pattern2numPresent),
    site2pattern(n.site2pattern),
    site2mask(n.site2mask),
#ifdef SIMD_ENABLED
    numSIMDBlocks(n.numSIMDBlocks),
    numAllocatedPatterns(n.numAllocatedPatterns),
#endif
    numPatterns( n.numPatterns ),
    numSiteRates(n.numSiteRates),
	numSiteFrequencies(n.numSiteFrequencies),
    numNodes(n.numNodes),
    numSites( n.numSites ),
    numTaxa(n.numTaxa),
    
    nodeOffset(n.nodeOffset),
	rateOffset(n.rateOffset),
    mixtureOffset(n.mixtureOffset),
    siteOffset(n.siteOffset),
    activeLikelihoodOffset(n.activeLikelihoodOffset),
    
    lnProb(n.lnProb),
    activeLikelihood( n.activeLikelihood ),
    partialLikelihoods(n.partialLikelihoods),
    per_site_Likelihoods(n.per_site_Likelihoods),

    useScaling(n.useScaling),
    
    // parameters
    tau( n.tau ),
	samplingRate( n.samplingRate ),
    siteRates(n.siteRates),
    homogeneousClockRate(n.homogeneousClockRate),
    heterogeneousClockRates(n.heterogeneousClockRates),
    branchHeterogeneousClockRates(n.branchHeterogeneousClockRates),
    homogeneousFrequency(n.homogeneousFrequency),
    heterogeneousFrequencies(n.heterogeneousFrequencies),
    branchLengths(n.branchLengths),
    
    branchHeterogeneousFrequencies(n.branchHeterogeneousFrequencies),
    rateVariationAcrossSites(n.rateVariationAcrossSites),
	frequencyVariationAcrossSites(n.frequencyVariationAcrossSites),

    activeProbability( n.activeProbability ),
    transitionProbabilities( n.transitionProbabilities ),
    tActiveOffset(n.tActiveOffset),
    tNodeOffset(n.tNodeOffset),
    tRateOffset(n.tRateOffset),
	tMixtureOffset(n.tMixtureOffset),
    
    N(n.N),
    numCorrectionMasks(n.numCorrectionMasks),
	correctionMaskOffset(n.correctionMaskOffset),
    activeCorrectionOffset(n.activeCorrectionOffset),
    correctionNodeOffset(n.correctionNodeOffset),
	correctionRateOffset(n.correctionRateOffset),
    correctionMixtureOffset(n.correctionMixtureOffset),
    correctionMaskMatrix(n.correctionMaskMatrix),
    correctionMaskCounts(n.correctionMaskCounts),
    maskObservationCounts(n.maskObservationCounts),
    
    perSiteCorrection(n.perSiteCorrection),
	perCodingProbs(n.perCodingProbs),
    perMaskCorrections(n.perMaskCorrections),
    perMixtureCorrections(n.perMixtureCorrections),
    correctionLikelihoods(n.correctionLikelihoods),
    
    perNodeCorrectionLogScalingFactors(n.perNodeCorrectionLogScalingFactors),
    activeCorrectionScalingOffset(n.activeCorrectionScalingOffset),
    perNodeSiteLogScalingFactors(n.perNodeSiteLogScalingFactors),
    activeScalingOffset(n.activeScalingOffset),

    dirtyNodes( n.dirtyNodes ),
    touchedNodes(n.touchedNodes),
    changedNodes( n.changedNodes ),
    
    verbose(n.verbose),
    continuous(n.continuous),
    scalingDensity(n.scalingDensity),
	countDistribution(n.countDistribution)
{   
    tau->getValue().getTreeChangeEventHandler().addListener( this );
}


void BinarySubstitutionModel::resizeLikelihoodVectors( void ) {

    // we resize the partial likelihood vectors to the new dimensions
    if(coding != AscertainmentBias::ALL)
    {
    	correctionMaskOffset		= 8;
    	correctionRateOffset 		= numCorrectionMasks*correctionMaskOffset;
        correctionMixtureOffset 	= numSiteRates*correctionRateOffset;
        correctionNodeOffset    	= numSiteFrequencies*correctionMixtureOffset;
        activeCorrectionOffset  	= numNodes*correctionNodeOffset;

        correctionLikelihoods = RealVector(2*activeCorrectionOffset, 0.0);
        perMixtureCorrections = RealVector(numSiteRates*numSiteFrequencies*numCorrectionMasks, 0.0);
        perMaskCorrections    = RealVector(numCorrectionMasks, 0.0);
        
        if(useScaling)
        {
            activeCorrectionScalingOffset = numCorrectionMasks*numNodes;
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
    numSIMDBlocks               = size_t((numPatterns - 1)/REALS_PER_SIMD_REGISTER) + 1;
    siteOffset                  = 2*REALS_PER_SIMD_REGISTER;
    rateOffset               	= numSIMDBlocks*siteOffset;
    mixtureOffset				= numSiteRates*rateOffset;
    
    numAllocatedPatterns        = numSIMDBlocks*REALS_PER_SIMD_REGISTER;
#else
    siteOffset                  = 2;
    rateOffset               	= numPatterns*siteOffset;
    mixtureOffset               = numSiteRates*rateOffset;
    
    size_t numAllocatedPatterns = numPatterns;
#endif
    nodeOffset                  = numSiteFrequencies*mixtureOffset;
    activeLikelihoodOffset      = numNodes*nodeOffset;
    
    partialLikelihoods   = RealVector(2 * activeLikelihoodOffset, 0.0);
    per_site_Likelihoods = RealVector(numAllocatedPatterns, 0.0);
    
    if(useScaling)
    {
        activeScalingOffset           = numAllocatedPatterns*numNodes;
        perNodeSiteLogScalingFactors = RealVector(2 * activeScalingOffset, 0.0);
    }
    else
    {
        perNodeSiteLogScalingFactors.clear();
    }
    
    // reset the transitionProbability vector
    tRateOffset    	= 4;
    tMixtureOffset 	= numSiteRates*tRateOffset;
    tNodeOffset 	= numSiteFrequencies*tMixtureOffset;
    tActiveOffset 	= (numNodes - 1)*tNodeOffset;
    
    transitionProbabilities = RealVector(2 * tActiveOffset, 0.0);
    
    // Reinitialize the tip data
    for(size_t nodeIndex = 0; nodeIndex < tau->getValue().getNumberOfTips(); nodeIndex++)
    {
        const TopologyNode& node = tau->getValue().getNode(nodeIndex);
        // this is a tip node
        // compute the likelihood for the tip and we are done
        computeNodeLikelihood(node, nodeIndex);
        
        if(coding != AscertainmentBias::ALL)
            computeNodeCorrection(node, nodeIndex);
    }
}



BinarySubstitutionModel::~BinarySubstitutionModel( void ) {
    // We don't delete the params, because they might be used somewhere else too. The model needs to do that!
    
    // remove myself from the tree listeners
    if ( tau != NULL ) 
    {
        // TODO: this needs to be implemented (Sebastian)
        tau->getValue().getTreeChangeEventHandler().removeListener( this );
    }
}




BinarySubstitutionModel* BinarySubstitutionModel::clone( void ) const
{
    
    return new BinarySubstitutionModel( *this );
}



const Tree* BinarySubstitutionModel::getTree( void ) const
{

    return &(tau->getValue());
}


void BinarySubstitutionModel::getIncludedSiteIndices( void )
{
    correctionMaskCounts.clear();
    maskObservationCounts.clear();
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
        maskObservationCounts.push_back(tips);
        correctionMaskMatrix.push_back(mask);
    }

    for (size_t i = 0; i < this->value->getNumberOfCharacters(); ++i)
    {
        std::pair<size_t, size_t> charCounts;
        size_t numGap = 0;

        std::string mask = "";
        std::vector<bool> maskData;
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

bool BinarySubstitutionModel::isSitePatternCompatible( std::pair<size_t,size_t> charCounts, size_t numGap)
{
    size_t zero = charCounts.first;
    size_t one  = charCounts.second;

    bool compatible = true;
    
    if(zero == (numTaxa - numGap) && (coding & AscertainmentBias::NOABSENCESITES) )
    {
        compatible = false;
    }
    else if(one == (numTaxa - numGap)  && (coding & AscertainmentBias::NOPRESENCESITES) )
    {
        compatible = false;
    }
    else if(zero == 1 && one == (numTaxa - numGap - 1) && (coding & AscertainmentBias::NOSINGLETONABSENCE) )
    {
        compatible = false;
    }
    else if(zero == (numTaxa - numGap - 1) && one == 1 && (coding & AscertainmentBias::NOSINGLETONPRESENCE) )
    {
        compatible = false;
    }
    else if(continuous && zero == (numTaxa - numGap - 1) && (coding & AscertainmentBias::NOSINGLETONPRESENCE) && (coding & AscertainmentBias::NOABSENCESITES) )
    {
        compatible = false;
    }
    else if(continuous && one == (numTaxa - numGap - 1) && (coding & AscertainmentBias::NOSINGLETONABSENCE) && (coding & AscertainmentBias::NOPRESENCESITES) )
    {
        compatible = false;
    }
    
    return compatible;
}


void BinarySubstitutionModel::compress( void )
{

    patternCounts.clear();
    pattern2site.clear();
    site2pattern.clear();
    numPatterns = 0;

    // resize the matrices
    numTaxa = tau->getValue().getNumberOfTips();
    
    // create a vector with the correct site indices
    // some of the sites may have been excluded
    getIncludedSiteIndices();
    
    if(numSites == 0)
        return;
    
    size_t siteIndex = siteIndices.back() + 1;
    
    // find the unique site patterns and compute their respective frequencies
    std::vector<TopologyNode*> nodes = tau->getValue().getNodes();

    std::vector<bool> unique(numSites, true);

    // compress the character matrix if we're asked to
    // find the unique site patterns and compute their respective frequencies
    std::map<std::string,size_t> patterns;
    for (size_t site = 0; site < numSites; ++site)
    {
        // create the site pattern
        std::string pattern = "";
        size_t numPresent = 0;
        for (std::vector<TopologyNode*>::iterator it = nodes.begin(); it != nodes.end(); ++it)
        {
            if ( (*it)->isTip() )
            {
                const BinaryTaxonData* taxon = this->value->getTaxonData( (*it)->getName() );
                
                std::stringstream ss;
                RealNumber c = taxon->getCharacter(siteIndices[site]);
                ss << c;
                pattern += ss.str();

                numPresent += (c > 0.0);
            }
        }
        // check if we have already seen this site pattern
        std::map<std::string, size_t>::const_iterator index = patterns.find( pattern );
        if ( index != patterns.end() )
        {
            // we have already seen this pattern
            // increase the frequency counter
            patternCounts[ index->second ]++;

            // obviously this site isn't unique nor the first encounter
            unique[site] = false;
            
            site2pattern.push_back(index->second);
        }
        else
        {
            // create a new pattern frequency counter for this pattern
            patternCounts.push_back(1);
            pattern2numPresent.push_back(numPresent);

            // insert this pattern with the corresponding index in the map
            patterns.insert( std::pair<std::string,size_t>(pattern,numPatterns) );

            // remember which site this pattern comes from
            pattern2site.push_back(site);
            site2pattern.push_back(numPatterns);
                        
            // increase the pattern counter
            numPatterns++;

            // flag that this site is unique (or the first occurence of this pattern)
            unique[site] = true;
        }
    }
    
    // finally we resize the partial likelihood vectors to the new pattern counts
    resizeLikelihoodVectors();
}



double BinarySubstitutionModel::computeLnProbability( void )
{

    updateTransitionProbabilities();
    
    // compute the ln probability by recursively calling the probability calculation for each node
    const TopologyNode &root = tau->getValue().getRoot();

    // we start with the root and then traverse down the tree
    size_t rootIndex = root.getIndex();

    // only necessary if the root is actually dirty
    if ( dirtyNodes[rootIndex] )
    {

        // start by filling the likelihood vector for the children of the root
        if ( root.getNumberOfChildren() == 2 ) // rooted trees have two children for the root
        {
            const TopologyNode &left = root.getChild(0);
            size_t leftIndex = left.getIndex();
            fillLikelihoodVector( left, leftIndex );
            const TopologyNode &right = root.getChild(1);
            size_t rightIndex = right.getIndex();
            fillLikelihoodVector( right, rightIndex );

            computeNodeLikelihood( root, rootIndex, leftIndex, rightIndex );
            scale(rootIndex, leftIndex, rightIndex);
            
            if(coding != AscertainmentBias::ALL)
            {
                computeNodeCorrection( root, rootIndex, leftIndex, rightIndex );
                scaleCorrection(rootIndex, leftIndex, rightIndex);
            }

        }
        else if ( root.getNumberOfChildren() == 3 ) // unrooted trees have three children for the root
        {
            const TopologyNode &left = root.getChild(0);
            size_t leftIndex = left.getIndex();
            fillLikelihoodVector( left, leftIndex );
            const TopologyNode &right = root.getChild(1);
            size_t rightIndex = right.getIndex();
            fillLikelihoodVector( right, rightIndex );
            const TopologyNode &middle = root.getChild(2);
            size_t middleIndex = middle.getIndex();
            fillLikelihoodVector( middle, middleIndex );

            computeNodeLikelihood( root, rootIndex, leftIndex, rightIndex, middleIndex );
            scale(rootIndex, leftIndex, rightIndex, middleIndex);
            
            if(coding != AscertainmentBias::ALL)
            {
                computeNodeCorrection( root, rootIndex, leftIndex, rightIndex, middleIndex );
                scaleCorrection(rootIndex, leftIndex, rightIndex, middleIndex);
            }

        }
        else
        {
            std::cerr << root.getNumberOfChildren() << std::endl;
        	throw Exception("The root node has an unexpected number of children. Only 2 (for rooted trees) or 3 (for unrooted trees) are allowed.");
        }


        // sum the partials up
        lnProb = sumRootLikelihood();

    }
    //std::cerr << lnProb << std::endl; 
    return lnProb;
}




void BinarySubstitutionModel::fillLikelihoodVector(const TopologyNode &node, size_t nodeIndex)
{    
    if ( !node.isTip() && dirtyNodes[nodeIndex]) 
    {   
        // mark as computed
        dirtyNodes[nodeIndex] = false;
        
        // this is an internal node
        // start by filling the likelihood vector for the two children of this node
        const TopologyNode &left = node.getChild(0);
        size_t leftIndex = left.getIndex();
        fillLikelihoodVector( left, leftIndex );
        const TopologyNode &right = node.getChild(1);
        size_t rightIndex = right.getIndex();
        fillLikelihoodVector( right, rightIndex );
        
        // now compute the likelihoods of this internal node
        computeNodeLikelihood(node,nodeIndex,leftIndex,rightIndex);
        scale(nodeIndex,leftIndex,rightIndex);
        
        if(coding != AscertainmentBias::ALL)
        {
            computeNodeCorrection(node,nodeIndex,leftIndex,rightIndex);
            scaleCorrection(nodeIndex, leftIndex, rightIndex);
        }
    }
}



void BinarySubstitutionModel::computeNodeLikelihood(const TopologyNode &node, size_t nodeIndex, size_t left, size_t right)
{   
    // get the pointers to the partial likelihoods for this node and the two descendant subtrees
    RealVector::const_iterator  p_left  = partialLikelihoods.begin() + activeLikelihood[left]*activeLikelihoodOffset      + left*nodeOffset;
    RealVector::const_iterator  p_right = partialLikelihoods.begin() + activeLikelihood[right]*activeLikelihoodOffset     + right*nodeOffset;
    RealVector::iterator        p_node  = partialLikelihoods.begin() + activeLikelihood[nodeIndex]*activeLikelihoodOffset + nodeIndex*nodeOffset;
    
    RealVector::iterator        pi_left   = transitionProbabilities.begin() + activeProbability[left]*tActiveOffset  + left*tNodeOffset;
    RealVector::iterator        pi_right  = transitionProbabilities.begin() + activeProbability[right]*tActiveOffset + right*tNodeOffset;

#ifdef SIMD_ENABLED
    SIMDRegister          m1, m2, m3, mrf0, mrf1;
#endif

    for (size_t freq = 0; freq < numSiteFrequencies; ++freq)
    {
    	// get the root frequencies
		RealNumber rf1 = getStationaryFrequency(nodeIndex, freq);
		RealNumber rf0 = 1.0 - rf1;

#ifdef SIMD_ENABLED
		if(node.isRoot())
		{
#ifdef AVX_ENABLED
			mrf0 = _mm256_broadcast_ss (&rf0);
			mrf1 = _mm256_broadcast_ss (&rf1);
#else
			mrf0 = _mm_load1_pX (&rf0);
			mrf1 = _mm_load1_pX (&rf1);
#endif
		}
#endif
    
		// iterate over all mixture categories
		for (size_t rate = 0; rate < numSiteRates; ++rate)
		{
			size_t tOffset = rate*tRateOffset + freq*tMixtureOffset;

			RealVector::iterator    t_left   = pi_left  + tOffset;
			RealVector::iterator    t_right  = pi_right + tOffset;

#ifdef SIMD_ENABLED
			for (size_t site = 0; site < numSIMDBlocks ; ++site)
			{
				size_t offset = rate*rateOffset + freq*mixtureOffset + site*siteOffset;

				SIMDRegister *          p_site_mixture    = (SIMDRegister *)&*(p_node  + offset);
				SIMDRegister *     p_site_mixture_left    = (SIMDRegister *)&*(p_left  + offset);
				SIMDRegister *    p_site_mixture_right    = (SIMDRegister *)&*(p_right + offset);
#ifdef AVX_ENABLED
				m1 = _mm256_broadcast_ss (&t_left[0]);
				m1 = _mm256_mul_ps (m1, p_site_mixture_left[0]);

				m2 = _mm256_broadcast_ss (&t_left[1]);
				m2 = _mm256_mul_ps (m2, p_site_mixture_left[1]);

				m3 = _mm256_add_ps (m1, m2);

				m1 = _mm256_broadcast_ss (&t_right[0]);
				m1 = _mm256_mul_ps (m1, p_site_mixture_right[0]);

				m2 = _mm256_broadcast_ss (&t_right[1]);
				m2 = _mm256_mul_ps (m2, p_site_mixture_right[1]);

				m1 = _mm256_add_ps (m1, m2);
				p_site_mixture[0] = _mm256_mul_ps (m1, m3);

				if(node.isRoot())
					p_site_mixture[0] = _mm256_mul_ps (mrf0, p_site_mixture[0]);

				m1 = _mm256_broadcast_ss (&t_left[2]);
				m1 = _mm256_mul_ps (m1, p_site_mixture_left[0]);

				m2 = _mm256_broadcast_ss (&t_left[3]);
				m2 = _mm256_mul_ps (m2, p_site_mixture_left[1]);

				m3 = _mm256_add_ps (m1, m2);


				m1 = _mm256_broadcast_ss (&t_right[2]);
				m1 = _mm256_mul_ps (m1, p_site_mixture_right[0]);

				m2 = _mm256_broadcast_ss (&t_right[3]);
				m2 = _mm256_mul_ps (m2, p_site_mixture_right[1]);

				m1 = _mm256_add_ps (m1, m2);
				p_site_mixture[1] = _mm256_mul_ps (m1, m3);

				if(node.isRoot())
					p_site_mixture[1] = _mm256_mul_ps (mrf1, p_site_mixture[1]);
#else
				m1 = _mm_load1_pX (&t_left[0]);
				m1 = _mm_mul_pX (m1, p_site_mixture_left[0]);

				m2 = _mm_load1_pX (&t_left[1]);
				m2 = _mm_mul_pX (m2, p_site_mixture_left[1]);

				m3 = _mm_add_pX (m1, m2);

				m1 = _mm_load1_pX (&t_right[0]);
				m1 = _mm_mul_pX (m1, p_site_mixture_right[0]);

				m2 = _mm_load1_pX (&t_right[1]);
				m2 = _mm_mul_pX (m2, p_site_mixture_right[1]);

				m1 = _mm_add_pX (m1, m2);
				p_site_mixture[0] = _mm_mul_pX (m1, m3);

				if(node.isRoot())
					p_site_mixture[0] = _mm_mul_pX (mrf0, p_site_mixture[0]);

				m1 = _mm_load1_pX (&t_left[2]);
				m1 = _mm_mul_pX (m1, p_site_mixture_left[0]);

				m2 = _mm_load1_pX (&t_left[3]);
				m2 = _mm_mul_pX (m2, p_site_mixture_left[1]);

				m3 = _mm_add_pX (m1, m2);


				m1 = _mm_load1_pX (&t_right[2]);
				m1 = _mm_mul_pX (m1, p_site_mixture_right[0]);

				m2 = _mm_load1_pX (&t_right[3]);
				m2 = _mm_mul_pX (m2, p_site_mixture_right[1]);

				m1 = _mm_add_pX (m1, m2);
				p_site_mixture[1] = _mm_mul_pX (m1, m3);

				if(node.isRoot())
					p_site_mixture[1] = _mm_mul_pX (mrf1, p_site_mixture[1]);
#endif
#else   
			// compute the per site probabilities
			for (size_t site = 0; site < numPatterns ; ++site)
			{
				size_t offset = rate*rateOffset + freq*mixtureOffset + site*siteOffset;

				RealVector::iterator          p_site_mixture          = p_node + offset;
				RealVector::const_iterator    p_site_mixture_left     = p_left + offset;
				RealVector::const_iterator    p_site_mixture_right    = p_right + offset;

				p_site_mixture[0] = ( p_site_mixture_left[0]  * t_left[0]  + p_site_mixture_left[1]  * t_left[1] )
								  * ( p_site_mixture_right[0] * t_right[0] + p_site_mixture_right[1] * t_right[1]);

				p_site_mixture[1] = ( p_site_mixture_left[0]  * t_left[2]  + p_site_mixture_left[1]  * t_left[3] )
								  * ( p_site_mixture_right[0] * t_right[2] + p_site_mixture_right[1] * t_right[3]);

				if(node.isRoot())
				{
					p_site_mixture[0] *= rf0;
					p_site_mixture[1] *= rf1;
				}
#endif            
			} // end-for over all sites (=patterns)

		} // end-for over all mixtures (=rate-categories)
	}
}



void BinarySubstitutionModel::computeNodeLikelihood(const TopologyNode &node, size_t nodeIndex, size_t left, size_t right, size_t middle)
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
    SIMDRegister          m1, m2, m3, mrf0, mrf1;
#endif

    for (size_t freq = 0; freq < numSiteFrequencies; ++freq)
    {
    	RealNumber rf1 = getStationaryFrequency(nodeIndex, freq);
		RealNumber rf0 = 1.0 - rf1;

#ifdef SIMD_ENABLED
		if(node.isRoot())
		{
#ifdef AVX_ENABLED
			mrf0 = _mm256_broadcast_ss (&rf0);
			mrf1 = _mm256_broadcast_ss (&rf1);
#else
			mrf0 = _mm_load1_pX (&rf0);
			mrf1 = _mm_load1_pX (&rf1);
#endif
		}
#endif    
    
		// iterate over all mixture categories
		for (size_t rate = 0; rate < numSiteRates; ++rate)
		{
			size_t tOffset = rate*tRateOffset + freq*tMixtureOffset;

			RealVector::iterator    t_left   = pi_left   + tOffset;
			RealVector::iterator    t_right  = pi_right  + tOffset;
			RealVector::iterator    t_middle = pi_middle + tOffset;

	#ifdef SIMD_ENABLED

			for (size_t site = 0; site < numSIMDBlocks ; ++site)
			{
				size_t offset = rate*rateOffset + freq*mixtureOffset + site*siteOffset;

				SIMDRegister *          p_site_mixture    = (SIMDRegister *)&*(p_node  + offset);
				SIMDRegister *     p_site_mixture_left    = (SIMDRegister *)&*(p_left  + offset);
				SIMDRegister *    p_site_mixture_right    = (SIMDRegister *)&*(p_right + offset);
				SIMDRegister *   p_site_mixture_middle    = (SIMDRegister *)&*(p_middle + offset);

	#ifdef AVX_ENABLED
				m1 = _mm256_broadcast_ss (&t_left[0]);
				m1 = _mm256_mul_ps (m1, p_site_mixture_left[0]);

				m2 = _mm256_broadcast_ss (&t_left[1]);
				m2 = _mm256_mul_ps (m2, p_site_mixture_left[1]);

				m3 = _mm256_add_ps (m1, m2);


				m1 = _mm256_broadcast_ss (&t_right[0]);
				m1 = _mm256_mul_ps (m1, p_site_mixture_right[0]);

				m2 = _mm256_broadcast_ss (&t_right[1]);
				m2 = _mm256_mul_ps (m2, p_site_mixture_right[1]);

				m1 = _mm256_add_ps (m1, m2);
				m3 = _mm256_mul_ps (m1, m3);


				m1 = _mm256_broadcast_ss (&t_middle[0]);
				m1 = _mm256_mul_ps (m1, p_site_mixture_middle[0]);

				m2 = _mm256_broadcast_ss (&t_middle[1]);
				m2 = _mm256_mul_ps (m2, p_site_mixture_middle[1]);

				m1 = _mm256_add_ps (m1, m2);
				p_site_mixture[0] = _mm256_mul_ps (m1, m3);

				if(node.isRoot())
					p_site_mixture[0] = _mm256_mul_ps (mrf0, p_site_mixture[0]);

				m1 = _mm256_broadcast_ss (&t_left[2]);
				m1 = _mm256_mul_ps (m1, p_site_mixture_left[0]);

				m2 = _mm256_broadcast_ss (&t_left[3]);
				m2 = _mm256_mul_ps (m2, p_site_mixture_left[1]);

				m3 = _mm256_add_ps (m1, m2);


				m1 = _mm256_broadcast_ss (&t_right[2]);
				m1 = _mm256_mul_ps (m1, p_site_mixture_right[0]);

				m2 = _mm256_broadcast_ss (&t_right[3]);
				m2 = _mm256_mul_ps (m2, p_site_mixture_right[1]);

				m1 = _mm256_add_ps (m1, m2);
				m3 = _mm256_mul_ps (m1, m3);


				m1 = _mm256_broadcast_ss (&t_middle[2]);
				m1 = _mm256_mul_ps (m1, p_site_mixture_middle[0]);

				m2 = _mm256_broadcast_ss (&t_middle[3]);
				m2 = _mm256_mul_ps (m2, p_site_mixture_middle[1]);

				m1 = _mm256_add_ps (m1, m2);
				p_site_mixture[1] = _mm256_mul_ps (m1, m3);

				if(node.isRoot())
					p_site_mixture[1] = _mm256_mul_ps (mrf1, p_site_mixture[1]);
	#else
				m1 = _mm_load1_pX (&t_left[0]);
				m1 = _mm_mul_pX (m1, p_site_mixture_left[0]);

				m2 = _mm_load1_pX (&t_left[1]);
				m2 = _mm_mul_pX (m2, p_site_mixture_left[1]);

				m3 = _mm_add_pX (m1, m2);


				m1 = _mm_load1_pX (&t_right[0]);
				m1 = _mm_mul_pX (m1, p_site_mixture_right[0]);

				m2 = _mm_load1_pX (&t_right[1]);
				m2 = _mm_mul_pX (m2, p_site_mixture_right[1]);

				m1 = _mm_add_pX (m1, m2);
				m3 = _mm_mul_pX (m1, m3);


				m1 = _mm_load1_pX (&t_middle[0]);
				m1 = _mm_mul_pX (m1, p_site_mixture_middle[0]);

				m2 = _mm_load1_pX (&t_middle[1]);
				m2 = _mm_mul_pX (m2, p_site_mixture_middle[1]);

				m1 = _mm_add_pX (m1, m2);
				p_site_mixture[0] = _mm_mul_pX (m1, m3);

				if(node.isRoot())
					p_site_mixture[0] = _mm_mul_pX (mrf0, p_site_mixture[0]);

				m1 = _mm_load1_pX (&t_left[2]);
				m1 = _mm_mul_pX (m1, p_site_mixture_left[0]);

				m2 = _mm_load1_pX (&t_left[3]);
				m2 = _mm_mul_pX (m2, p_site_mixture_left[1]);

				m3 = _mm_add_pX (m1, m2);


				m1 = _mm_load1_pX (&t_right[2]);
				m1 = _mm_mul_pX (m1, p_site_mixture_right[0]);

				m2 = _mm_load1_pX (&t_right[3]);
				m2 = _mm_mul_pX (m2, p_site_mixture_right[1]);

				m1 = _mm_add_pX (m1, m2);
				m3 = _mm_mul_pX (m1, m3);


				m1 = _mm_load1_pX (&t_middle[2]);
				m1 = _mm_mul_pX (m1, p_site_mixture_middle[0]);

				m2 = _mm_load1_pX (&t_middle[3]);
				m2 = _mm_mul_pX (m2, p_site_mixture_middle[1]);

				m1 = _mm_add_pX (m1, m2);
				p_site_mixture[1] = _mm_mul_pX (m1, m3);

				if(node.isRoot())
					p_site_mixture[1] = _mm_mul_pX (mrf1, p_site_mixture[1]);
	#endif
	#else
			// compute the per site probabilities
			for (size_t site = 0; site < numPatterns ; ++site)
			{
				size_t offset = rate*rateOffset + freq*mixtureOffset + site*siteOffset;

				RealVector::iterator          p_site_mixture          = p_node   + offset;
				RealVector::const_iterator    p_site_mixture_left     = p_left   + offset;
				RealVector::const_iterator    p_site_mixture_right    = p_right  + offset;
				RealVector::const_iterator    p_site_mixture_middle   = p_middle + offset;

				p_site_mixture[0] = ( p_site_mixture_left[0]   * t_left[0]    + p_site_mixture_left[1]   * t_left[1]  )
								  * ( p_site_mixture_right[0]  * t_right[0]   + p_site_mixture_right[1]  * t_right[1] )
								  * ( p_site_mixture_middle[0] * pi_middle[0] + p_site_mixture_middle[1] * pi_middle[1]);

				p_site_mixture[1] = ( p_site_mixture_left[0]   * t_left[2]   + p_site_mixture_left[1]   * t_left[3]  )
								  * ( p_site_mixture_right[0]  * t_right[2]  + p_site_mixture_right[1]  * t_right[3] )
								  * ( p_site_mixture_middle[0] * t_middle[2] + p_site_mixture_middle[1] * t_middle[3]);

				if(node.isRoot())
				{
					p_site_mixture[0] *= rf0;
					p_site_mixture[1] *= rf1;
				}
	#endif
			} // end-for over all sites (=patterns)

		} // end-for over all mixtures (=rate-categories)
	}
}


void BinarySubstitutionModel::computeNodeLikelihood(const TopologyNode &node, size_t nodeIndex) 
{    
    const BinaryTaxonData* td = this->value->getTaxonData(nodeIndex);
    
    for(size_t active = 0; active < 2; active++)
    {
        RealVector::iterator p_node = partialLikelihoods.begin() + active*activeLikelihoodOffset + nodeIndex*nodeOffset;
        
        // iterate over all mixture categories
        for (size_t freq = 0; freq < numSiteFrequencies; ++freq)
        {

			for (size_t rate = 0; rate < numSiteRates; ++rate)
			{
#ifdef SIMD_ENABLED
				// iterate over all sites

				for (size_t site = 0; site < numSIMDBlocks; site++)
				{
					RealVector::iterator     p_site_mixture      = p_node + rate*rateOffset + freq*mixtureOffset + site*siteOffset;

					for (size_t ss = 0; ss < REALS_PER_SIMD_REGISTER; ss++)
					{
						size_t pattern = site*REALS_PER_SIMD_REGISTER + ss;
						size_t index   = ss;

						if(pattern >= numPatterns)
							break;

						// is this site a gap?
						if ( td->getGap(siteIndices[pattern2site[pattern]]) )
						{
							// since this is a gap we need to assume that the actual state could have been any state
							p_site_mixture[index] = 1.0;
							p_site_mixture[index + REALS_PER_SIMD_REGISTER] = 1.0;
						}
						else // we have observed a character
						{
							// get the original character
							RealNumber c = td->getCharacter(siteIndices[pattern2site[pattern]]);

							p_site_mixture[index] = 1.0 - c;
							p_site_mixture[index + REALS_PER_SIMD_REGISTER] = c;

						}
					}
#else            
				// iterate over all sites
				for (size_t pattern = 0; pattern < numPatterns; pattern++)
				{
					RealVector::iterator     p_site_mixture      = p_node + rate*rateOffset + freq*mixtureOffset + pattern*siteOffset;

					// is this site a gap?
					if ( td->getGap(siteIndices[pattern2site[pattern]]) )
					{
						// since this is a gap we need to assume that the actual state could have been any state
						p_site_mixture[0] = 1.0;
						p_site_mixture[1] = 1.0;
					}
					else // we have observed a character
					{
						// get the original character
						RealNumber c = td->getCharacter(siteIndices[pattern2site[pattern]]);

						p_site_mixture[0] = 1.0 - c;
						p_site_mixture[1] = c;

					}
#endif                
				}

			} // end-for over all mixture categories
		}
    }
}



void BinarySubstitutionModel::computeNodeCorrection(const TopologyNode &node, size_t nodeIndex)
{
    for(size_t active = 0; active < 2; active++)
    {
        RealVector::iterator p_node = correctionLikelihoods.begin() + active*activeCorrectionOffset + nodeIndex*correctionNodeOffset;
    
        // iterate over all mixture categories
        for (size_t freq = 0; freq < numSiteFrequencies; ++freq)
        {

			for (size_t rate = 0; rate < numSiteRates; ++rate)
			{
				for(size_t mask = 0; mask < numCorrectionMasks; mask++)
				{
					size_t offset = rate*correctionRateOffset + freq*correctionMixtureOffset + mask*8;

					RealVector::iterator         u      = p_node + offset;

					bool gap = correctionMaskMatrix[mask][nodeIndex];

					for(size_t ci = 0; ci < 2; ci++)
					{
						RealVector::iterator         uC = u  + ci*2;
						RealVector::iterator         uI = uC + 4;

						for(size_t c = 0; c < 2; c++)
						{

							// Probability of constant state c this tip
							// when the state at this tip is ci
							uC[c] = (c == ci) && !gap;

							// Probability of invert singleton state c this tip
							// when the state at this tip is ci
							uI[c] = (c != ci) && !gap;
						}
					}
				}

			}
        }
    }
}



void BinarySubstitutionModel::computeNodeCorrection(const TopologyNode &node, size_t nodeIndex, size_t left, size_t right, size_t middle)
{
    // get the pointers to the partial likelihoods for this node and the two descendant subtrees
    RealVector::const_iterator   p_left   = correctionLikelihoods.begin() + activeLikelihood[left]*activeCorrectionOffset + left*correctionNodeOffset;
    RealVector::const_iterator   p_right  = correctionLikelihoods.begin() + activeLikelihood[right]*activeCorrectionOffset + right*correctionNodeOffset;
    RealVector::const_iterator   p_middle = correctionLikelihoods.begin() + activeLikelihood[middle]*activeCorrectionOffset + middle*correctionNodeOffset;
    RealVector::iterator         p_node   = correctionLikelihoods.begin() + activeLikelihood[nodeIndex]*activeCorrectionOffset + nodeIndex*correctionNodeOffset;

    RealVector::iterator    pi_left   = transitionProbabilities.begin() + activeProbability[left]*tActiveOffset + left*tNodeOffset;
    RealVector::iterator    pi_right  = transitionProbabilities.begin() + activeProbability[right]*tActiveOffset + right*tNodeOffset;
    RealVector::iterator    pi_middle = transitionProbabilities.begin() + activeProbability[middle]*tActiveOffset + middle*tNodeOffset;
    
    // iterate over all mixture categories
    for (size_t freq = 0; freq < numSiteFrequencies; ++freq)
    {
#ifdef SIMD_CORRECTION_ENABLED
		RealNumber rf1 = getStationaryFrequency(root);
		RealNumber rf0 = 1.0 - rf1;

		__m128 mrf = _mm_setr_ps(rf0, rf0, rf1, rf1);
#else
		std::vector<RealNumber> f(2, 1.0);
		f[1] = getStationaryFrequency(nodeIndex, freq);
		f[0] = 1.0 - f[1];
#endif
		for (size_t rate = 0; rate < numSiteRates; ++rate)
		{
			size_t tOffset = rate*tRateOffset + freq*tMixtureOffset;

#ifdef SIMD_CORRECTION_ENABLED
			__m128*    t_left   = (__m128*)&*(pi_left   + tOffset);
			__m128*    t_right  = (__m128*)&*(pi_right  + tOffset);
			__m128*    t_middle = (__m128*)&*(pi_middle  + tOffset);
#else
			RealVector::iterator    t_left   = pi_left   + tOffset;
			RealVector::iterator    t_right  = pi_right  + tOffset;
			RealVector::iterator    t_middle = pi_middle + tOffset;
#endif

			for(size_t mask = 0; mask < numCorrectionMasks; mask++){

				size_t offset = rate*correctionRateOffset + freq*correctionMixtureOffset + mask*8;

#ifdef SIMD_CORRECTION_ENABLED
				__m128* uC_i = (__m128*)&*(p_node  + offset);
				__m128* uC_j = (__m128*)&*(p_left  + offset);
				__m128* uC_k = (__m128*)&*(p_right + offset);
				__m128* uC_l = (__m128*)&*(p_middle + offset);

				__m128* uI_i = (__m128*)&*(p_node  + offset + 4);
				__m128* uI_j = (__m128*)&*(p_left  + offset + 4);
				__m128* uI_k = (__m128*)&*(p_right + offset + 4);
				__m128* uI_l = (__m128*)&*(p_middle + offset + 4);

				*uC_i = _mm_setzero_ps();
				*uI_i = _mm_setzero_ps();

				__m128 m1a,m1b,m1c,m2a,m2b,m2c,m3a,m3b,m3c,m4a,m4b,m4c,mA1,mA2,mA2,mA;

				for(size_t i = 0; i < 4; i++)
				{

					__m128* rC_j = i == 0 ? uI_j : uC_j;
					__m128* rC_k = i == 1 ? uI_k : uC_k;
					__m128* rC_l = i == 2 ? uI_l : uC_l;
					__m128* rC_i = i == 3 ? uC_i : uI_i;

					// compute for constant sites
					m1a = _mm_shuffle_ps(*t_left,   *t_right,  _MM_SHUFFLE(3, 2, 1, 0) );
					m1b = _mm_shuffle_ps(*t_right,  *t_middle, _MM_SHUFFLE(3, 2, 1, 0) );
					m1c = _mm_shuffle_ps(*t_middle, *t_left,   _MM_SHUFFLE(3, 2, 1, 0) );

					m2a = _mm_shuffle_ps(*rC_j, *rC_k, _MM_SHUFFLE(3, 0, 3, 0) );
					m2b = _mm_shuffle_ps(*rC_k, *rC_l, _MM_SHUFFLE(3, 0, 3, 0) );
					m2c = _mm_shuffle_ps(*rC_l, *rC_j, _MM_SHUFFLE(3, 0, 3, 0) );

					m3a = _mm_shuffle_ps(*t_left,   *t_right,  _MM_SHUFFLE(2, 3, 0, 1) );
					m3b = _mm_shuffle_ps(*t_right,  *t_middle, _MM_SHUFFLE(2, 3, 0, 1) );
					m3c = _mm_shuffle_ps(*t_middle, *t_left,   _MM_SHUFFLE(2, 3, 0, 1) );

					m4a = _mm_shuffle_ps(*rC_j, *rC_k, _MM_SHUFFLE(1, 2, 1, 2) );
					m4b = _mm_shuffle_ps(*rC_k, *rC_l, _MM_SHUFFLE(1, 2, 1, 2) );
					m4c = _mm_shuffle_ps(*rC_l, *rC_j, _MM_SHUFFLE(1, 2, 1, 2) );

					for(size_t j = 0; j < 8; j++)
					{
						mA1 = (j & 1) ? _mm_mul_pX(m1a,m2a) : _mm_mul_pX(m3a,m4a);
						mA2 = (j & 2) ? _mm_mul_pX(m1b,m2b) : _mm_mul_pX(m3b,m4b);
						mA3 = (j & 4) ? _mm_mul_pX(m1c,m2c) : _mm_mul_pX(m3c,m4c);
						mA =  _mm_mul_pX(mA1,mA2);
						mA =  _mm_mul_pX(mA,mA3);
						*rC_i = _mm_add_pX(*rC_i,mA);
					}
				}

				if(node.isRoot())
				{
					*uC_i = _mm_mul_pX(*uC_i,mrf);
					*uI_i = _mm_mul_pX(*uI_i,mrf);
				}
#else
				RealVector::iterator               u_i = p_node   + offset;
				RealVector::const_iterator         u_j = p_left   + offset;
				RealVector::const_iterator         u_k = p_right  + offset;
				RealVector::const_iterator         u_l = p_middle + offset;

				for(size_t ci = 0; ci < 2; ci++)
				{
					RealVector::iterator         uC_i = u_i  + ci*2;
					RealVector::iterator         uI_i = uC_i + 4;

					for(size_t c = 0; c < 2; c++)
					{

						uC_i[c] = 0.0;
						uI_i[c] = 0.0;

						for(size_t cj = 0; cj < 2; cj++)
						{
							RealVector::const_iterator         uC_j = u_j  + cj*2;
							RealVector::const_iterator         uI_j = uC_j + 4;

							for(size_t ck = 0; ck < 2; ck++)
							{
								RealVector::const_iterator         uC_k = u_k  + ck*2;
								RealVector::const_iterator         uI_k = uC_k + 4;

								for(size_t cl = 0; cl < 2; cl++)
								{
									RealVector::const_iterator         uC_l = u_l  + cl*2;
									RealVector::const_iterator         uI_l = uC_l + 4;

									RealNumber Pij =   t_left[2*ci + cj];
									RealNumber Pik =  t_right[2*ci + ck];
									RealNumber Pil = t_middle[2*ci + cl];

									// probability of constant state c descending from this node
									// when the state at this node is ci, with children states cj, ck, and cl
									uC_i[c] += Pij*uC_j[c] * Pik*uC_k[c] * Pil*uC_l[c];

									// probability of invert singleton state c descending from
									// when the state at this node is ci, with children states cj, ck, and cl
									uI_i[c] += Pij*uI_j[c] * Pik*uC_k[c] * Pil*uC_l[c]
											 + Pij*uC_j[c] * Pik*uI_k[c] * Pil*uC_l[c]
											 + Pij*uC_j[c] * Pik*uC_k[c] * Pil*uI_l[c];

								}
							}
						}

						if(node.isRoot())
						{
							uC_i[c] *= f[ci];
							uI_i[c] *= f[ci];
						}
					}
				}
#endif
			}
		}
    }
}



void BinarySubstitutionModel::computeNodeCorrection(const TopologyNode &node, size_t nodeIndex, size_t left, size_t right)
{
    // get the pointers to the partial likelihoods for this node and the two descendant subtrees
    RealVector::const_iterator   p_left  = correctionLikelihoods.begin() + activeLikelihood[left]*activeCorrectionOffset + left*correctionNodeOffset;
    RealVector::const_iterator   p_right = correctionLikelihoods.begin() + activeLikelihood[right]*activeCorrectionOffset + right*correctionNodeOffset;
    RealVector::iterator         p_node  = correctionLikelihoods.begin() + activeLikelihood[nodeIndex]*activeCorrectionOffset + nodeIndex*correctionNodeOffset;

    RealVector::iterator    pi_left   = transitionProbabilities.begin() + activeProbability[left]*tActiveOffset + left*tNodeOffset;
    RealVector::iterator    pi_right  = transitionProbabilities.begin() + activeProbability[right]*tActiveOffset + right*tNodeOffset;

    // iterate over all mixture categories
    for (size_t freq = 0; freq < numSiteFrequencies; ++freq)
    {
#ifdef SIMD_CORRECTION_ENABLED
		RealNumber rf1 = getStationaryFrequency(root);
		RealNumber rf0 = 1.0 - rf1;

		__m128 mrf = _mm_setr_ps(rf0, rf0, rf1, rf1);
#else
		std::vector<RealNumber> f(2, 1.0);
		f[1] = getStationaryFrequency(nodeIndex, freq);
		f[0] = 1.0 - f[1];
#endif

		for (size_t rate = 0; rate < numSiteRates; ++rate)
		{
			size_t tOffset = rate*tRateOffset + freq*tMixtureOffset;

#ifdef SIMD_CORRECTION_ENABLED
			__m128*    t_left   = (__m128*)&*(pi_left   + tOffset);
			__m128*    t_right  = (__m128*)&*(pi_right  + tOffset);
#else
			RealVector::iterator    t_left   = pi_left   + tOffset;
			RealVector::iterator    t_right  = pi_right  + tOffset;
#endif

			for(size_t mask = 0; mask < numCorrectionMasks; mask++){

				size_t offset = rate*correctionRateOffset + freq*correctionMixtureOffset + mask*8;

#ifdef SIMD_CORRECTION_ENABLED
				__m128* uC_i = (__m128*)&*(p_node  + offset);
				__m128* uC_j = (__m128*)&*(p_left  + offset);
				__m128* uC_k = (__m128*)&*(p_right + offset);

				__m128* uI_i = (__m128*)&*(p_node  + offset + 4);
				__m128* uI_j = (__m128*)&*(p_left  + offset + 4);
				__m128* uI_k = (__m128*)&*(p_right + offset + 4);

				*uC_i = _mm_setzero_ps();
				*uI_i = _mm_setzero_ps();

				__m128 m1a,m1b,m2a,m2b,m3a,m3b,m4a,m4b,mA1,mA2,mA;

				// compute for constant sites
				for(size_t i = 0; i < 3; i++)
				{
					__m128* rC_j = i == 0 ? uI_j : uC_j;
					__m128* rC_k = i == 1 ? uI_k : uC_k;
					__m128* rC_i = i == 2 ? uC_i : uI_i;

					// compute for constant sites
					m1a = _mm_shuffle_ps(*t_left,   *t_right,  _MM_SHUFFLE(3, 2, 1, 0) );
					m1b = _mm_shuffle_ps(*t_right, *t_left,   _MM_SHUFFLE(3, 2, 1, 0) );

					m2a = _mm_shuffle_ps(*rC_j, *rC_k, _MM_SHUFFLE(3, 0, 3, 0) );
					m2b = _mm_shuffle_ps(*rC_k, *rC_j, _MM_SHUFFLE(3, 0, 3, 0) );

					m3a = _mm_shuffle_ps(*t_left,   *t_right,  _MM_SHUFFLE(2, 3, 0, 1) );
					m3b = _mm_shuffle_ps(*t_right, *t_left,   _MM_SHUFFLE(2, 3, 0, 1) );

					m4a = _mm_shuffle_ps(*rC_j, *rC_k, _MM_SHUFFLE(1, 2, 1, 2) );
					m4b = _mm_shuffle_ps(*rC_k, *rC_j, _MM_SHUFFLE(1, 2, 1, 2) );

					for(size_t j = 0; j < 4; j++)
					{
						mA1 = (j & 1) ? _mm_mul_pX(m1a,m2a) : _mm_mul_pX(m3a,m4a);
						mA2 = (j & 2) ? _mm_mul_pX(m1b,m2b) : _mm_mul_pX(m3b,m4b);
						mA =  _mm_mul_pX(mA1,mA2);
						*rC_i = _mm_add_pX(*rC_i,mA);
					}
				}

				if(node.isRoot())
				{
					*uC_i = _mm_mul_pX(*uC_i,mrf);
					*uI_i = _mm_mul_pX(*uI_i,mrf);
				}
#else
				RealVector::iterator               u_i = p_node  + offset;
				RealVector::const_iterator         u_j = p_left  + offset;
				RealVector::const_iterator         u_k = p_right + offset;

				for(size_t ci = 0; ci < 2; ci++)
				{
					RealVector::iterator         uC_i = u_i  + ci*2;
					RealVector::iterator         uI_i = uC_i + 4;

					for(size_t c = 0; c < 2; c++)
					{

						uC_i[c] = 0.0;
						uI_i[c] = 0.0;

						for(size_t cj = 0; cj < 2; cj++)
						{
							RealVector::const_iterator         uC_j = u_j  + cj*2;
							RealVector::const_iterator         uI_j = uC_j + 4;

							for(size_t ck = 0; ck < 2; ck++)
							{
								RealVector::const_iterator         uC_k = u_k  + ck*2;
								RealVector::const_iterator         uI_k = uC_k + 4;

								RealNumber Pij =  t_left[2*ci + cj];
								RealNumber Pik = t_right[2*ci + ck];

								// probability of constant state c descending from this node
								// when the state at this node is ci, with children states cj and ck
								uC_i[c] += Pij*uC_j[c] * Pik*uC_k[c];

								// probability of invert singleton state c descending from this node
								// when the state at this node is ci, with children states cj and ck
								uI_i[c] += Pij*uC_j[c] * Pik*uI_k[c] + Pij*uI_j[c] * Pik*uC_k[c];
							}
						}

						if(node.isRoot())
						{
							uC_i[c] *= f[ci];
							uI_i[c] *= f[ci];
						}
					}
				}
#endif
			}
		}
    }
}


RealNumber BinarySubstitutionModel::sumRootLikelihood( void )
{
    RealNumber sumPartialProbs = sumUncorrectedRootLikelihood();
    
    if(coding == AscertainmentBias::ALL)
        return sumPartialProbs;
                         
    // get the root node
    const TopologyNode &root = this->tau->getValue().getRoot();

    // get the index of the root node
    size_t nodeIndex = root.getIndex();
    
    double sampling = 1.0 - getSamplingRate();

    RealVector::const_iterator p_node = correctionLikelihoods.begin() + activeLikelihood[nodeIndex] * activeCorrectionOffset  + nodeIndex*correctionNodeOffset;
    
    std::fill(perMaskCorrections.begin(), perMaskCorrections.end(), 0.0);
    perCodingProbs = std::vector<RealNumber>(4, 0.0);
    
    // iterate over each correction mask
    for(size_t mask = 0; mask < numCorrectionMasks; mask++)
    {   
        RealNumber logScalingFactor = useScaling ? perNodeCorrectionLogScalingFactors[activeLikelihood[nodeIndex]*activeCorrectionScalingOffset + nodeIndex*numCorrectionMasks + mask] : 0.0;
        // iterate over all mixture categories
        for (size_t freq = 0; freq < numSiteFrequencies; ++freq)
        {

			for (size_t rate = 0; rate < numSiteRates; ++rate)
			{
				size_t offset = rate*correctionRateOffset + freq*correctionMixtureOffset + mask*8;

				RealVector::const_iterator         u_i = p_node   + offset;

				RealNumber max = -std::numeric_limits<RealNumber>::infinity();

				RealNumber prob = 0.0;

				std::vector<RealNumber> logCorrections;

				for(size_t ci = 0; ci < 2; ci++)
				{
					// constant site pattern likelihoods
					RealVector::const_iterator         uC_i = u_i  + ci*2;
					// invert singleton likelihoods
					RealVector::const_iterator         uI_i = uC_i + 4;

					for(size_t c = 0; c < 2; c++)
					{
						RealNumber tmp = 0.0;

						// c is the character state of the correction pattern
						if(c == 0)
						{
							if(coding & AscertainmentBias::NOABSENCESITES)
								tmp += uC_i[c];

							if(coding & AscertainmentBias::NOSINGLETONPRESENCE)
								tmp += uI_i[c];
							else
								tmp += uI_i[c] * sampling;

							if(mask == 0)
							{
								perCodingProbs[0] += uC_i[c];
								perCodingProbs[1] += uI_i[c];
							}
						}

						if(c == 1)
						{
							// if there is only one observed tip, then don't double-count singleton gains
							if((coding & AscertainmentBias::NOPRESENCESITES) && maskObservationCounts[mask] > 1)
								tmp += uC_i[c];

							// if there are only two observed tips, then don't double-count singleton gains
							// if there is only one observed tip, then don't double-count absence sites
							if((coding & AscertainmentBias::NOSINGLETONABSENCE) && maskObservationCounts[mask] > 2)
								tmp += uI_i[c];

							if(mask == 0)
							{
								perCodingProbs[2] += uI_i[c];
								perCodingProbs[3] += uC_i[c];
							}
						}

						if(tmp == 0.0)
							continue;

						if(useScaling)
						{
							tmp = log(tmp) + logScalingFactor;

							max = std::max(tmp, max);

							logCorrections.push_back(tmp);
						}
						else
						{
							prob += tmp;
						}
					}
				}

				if(useScaling)
				{
					//std::cerr << logScalingFactor << std::endl;
					//max = std::max(logScalingFactor, max);

					//logCorrections.push_back(logScalingFactor);
					// use the log-exp-sum to get the sum of the corrections
					prob = exp(Math::log_sum_exp(logCorrections, max));
				}

				perMaskCorrections[mask] += prob;

				RealNumber mixprob = 1.0 - prob;

				// impose a boundary
				if(mixprob <= 0)
					mixprob = 0;

				if(mask == 0)
					perMixtureCorrections[freq*numSiteRates*numCorrectionMasks + rate*numCorrectionMasks] = mixprob;
			}
        }

        if(mask == 0)
        	for(size_t i = 0; i < 4; i++)
        		perCodingProbs[i] /= numSiteRates*numSiteFrequencies;

        // normalize and invert the probability
        perMaskCorrections[mask] = 1.0 - perMaskCorrections[mask]/(numSiteRates*numSiteFrequencies);

        if(perMaskCorrections[mask] <= 0)
        	perMaskCorrections[mask] = std::numeric_limits<RealNumber>::infinity();

        // log transform
        perMaskCorrections[mask] = log(perMaskCorrections[mask]);
        
        //std::cerr << perMaskCorrections[mask]*correctionMaskCounts[mask] << "\t" << perCodingProbs << std::endl;

        // apply the correction for this correction mask
        sumPartialProbs -= perMaskCorrections[mask]*correctionMaskCounts[mask];
    }

    return sumPartialProbs;
}


RealNumber BinarySubstitutionModel::sumUncorrectedRootLikelihood( void )
{
    // get the root node
    const TopologyNode &root = tau->getValue().getRoot();

    // get the index of the root node
    size_t rootIndex = root.getIndex();

	double logSampling = log(getSamplingRate());

    // get the pointers to the partial likelihoods of the left and right subtree
    RealVector::iterator p_node  = partialLikelihoods.begin() + activeLikelihood[rootIndex] * activeLikelihoodOffset  + rootIndex*nodeOffset;
   
    //reset the per site likelihood
    std::fill(per_site_Likelihoods.begin(), per_site_Likelihoods.end(), 0.0);
    // iterate over all mixture categories
	for (size_t freq = 0; freq < numSiteFrequencies; ++freq)
    {
		for (size_t rate = 0; rate < numSiteRates; ++rate)
		{

#ifdef SIMD_ENABLED
			SIMDRegister * mTotals = (SIMDRegister *)&per_site_Likelihoods[0];

			for (size_t pattern = 0; pattern < numSIMDBlocks; pattern++)
			{
				// get the pointers to the likelihoods for this site and mixture category
				SIMDRegister *   p_site_mixture     = (SIMDRegister *)&*(p_node + pattern*siteOffset + rate*rateOffset + freq*mixtureOffset);
#ifdef AVX_ENABLED
				*mTotals = _mm256_add_ps(*mTotals, _mm256_add_ps(p_site_mixture[0], p_site_mixture[1]));
#else
				*mTotals = _mm_add_pX(*mTotals, _mm_add_pX(p_site_mixture[0], p_site_mixture[1]));
#endif

				mTotals++;

			} // end-for over all mixtures (=rate categories)
		}
    }

	// sum the log-likelihoods for all sites together
	RealNumber sumPartialProbs = 0.0;

	for (size_t block = 0; block < numSIMDBlocks; block++)
	{
		for (size_t ss = 0; ss < REALS_PER_SIMD_REGISTER; ss++)
		{
			size_t pattern = block*REALS_PER_SIMD_REGISTER + ss;

			if(pattern >= numPatterns)
				continue;

			per_site_Likelihoods[pattern] = log( per_site_Likelihoods[pattern] / (numSiteRates * numSiteFrequencies) );

			if ( useScaling == true )
			{
				per_site_Likelihoods[pattern] += perNodeSiteLogScalingFactors[activeLikelihood[rootIndex]*activeScalingOffset + rootIndex*numAllocatedPatterns + pattern];
			}

			per_site_Likelihoods[pattern] += (pattern2numPresent[pattern] == 1) ? logSampling : 0.0;

			sumPartialProbs += per_site_Likelihoods[pattern]*patternCounts[pattern];
		}
	}
#else
			for (size_t pattern = 0; pattern < numPatterns; ++pattern)
			{
				// get the pointers to the likelihoods for this site and mixture category
				RealVector::iterator   p_site_mixture     = p_node + pattern*siteOffset + rate*rateOffset + freq*mixtureOffset;

				per_site_Likelihoods[pattern] += p_site_mixture[0] + p_site_mixture[1];

			}

		} // end-for over all mixtures (=rate categories)
	}

	// sum the log-likelihoods for all sites together
	double sumPartialProbs = 0.0;

	for (size_t pattern = 0; pattern < numPatterns; ++pattern)
	{
		per_site_Likelihoods[pattern] = log( per_site_Likelihoods[pattern] / (numSiteRates * numSiteFrequencies) );

		if ( useScaling == true )
		{

			per_site_Likelihoods[pattern] += perNodeSiteLogScalingFactors[activeLikelihood[rootIndex]*activeScalingOffset + rootIndex*numPatterns + pattern] * patternCounts[pattern];
		}

		per_site_Likelihoods[pattern] += (pattern2numPresent[pattern] == 1) ? logSampling : 0.0;

		sumPartialProbs += per_site_Likelihoods[pattern]*patternCounts[pattern];
	}
#endif
	//std::cerr << sumPartialProbs << std::endl;
    
	return sumPartialProbs;
}



void BinarySubstitutionModel::scale( size_t nodeIndex, size_t left, size_t right )
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
            for (size_t freq = 0; freq < numSiteFrequencies; ++freq)
            {
				for (size_t rate = 0; rate < numSiteRates; ++rate)
				{
					// get the pointers to the likelihood for this mixture category
					size_t offset = rate*rateOffset + freq*mixtureOffset + site*siteOffset;

					SIMDRegister *          p_site_mixture = (SIMDRegister*)&*(p_node + offset);
#ifdef AVX_ENABLED
					SIMDRegister tmp = _mm256_max_ps(p_site_mixture[0],p_site_mixture[1]);
					if(freq == 0 && rate == 0)
						max = tmp;
					else
						max = _mm256_max_ps(max, tmp);
#else
					SIMDRegister tmp = _mm_max_pX(p_site_mixture[0],p_site_mixture[1]);
					if(freq == 0 && rate == 0)
						max = tmp;
					else
						max = _mm_max_pX(max, tmp);
#endif
				}
            }

#ifdef AVX_ENABLED
            p_scaler[site] = _mm256_add_ps(p_scaler_left[site],p_scaler_right[site]);
#ifdef DOUBLE_PRECISION
            SIMDRegister tmp = max;
            double* ptmp = (double*)&tmp;
            double* pmax = (double*)&max;
            for(size_t i = 0; i < REALS_PER_SIMD_REGISTER; i++)
            	ptmp[i] = log(pmax[i]);

            p_scaler[site] = _mm256_add_pX(p_scaler[site], tmp);
#else
            p_scaler[site] = _mm256_add_pX(p_scaler[site], log256_ps(max));
#endif
#else
            p_scaler[site] = _mm_add_pX(p_scaler_left[site],p_scaler_right[site]);
#ifdef DOUBLE_PRECISION
            SIMDRegister tmp = max;
            double* ptmp = (double*)&tmp;
            double* pmax = (double*)&max;
            for(size_t i = 0; i < REALS_PER_SIMD_REGISTER; i++)
            	ptmp[i] = log(pmax[i]);

            p_scaler[site] = _mm_add_pX(p_scaler[site], tmp);
#else
            p_scaler[site] = _mm_add_pX(p_scaler[site], log_ps(max));
#endif
#endif

            // compute the per site probabilities
			for (size_t freq = 0; freq < numSiteFrequencies; ++freq)
			{
				// compute the per site probabilities
				for (size_t rate = 0; rate < numSiteRates; ++rate)
				{
					// get the pointers to the likelihood for this mixture category
					size_t offset = rate*rateOffset + freq*mixtureOffset + site*siteOffset;

					SIMDRegister *          p_site_mixture = (SIMDRegister*)&*(p_node + offset);
#ifdef AVX_ENABLED
					p_site_mixture[0] = _mm256_div_ps(p_site_mixture[0], max);
					p_site_mixture[1] = _mm256_div_ps(p_site_mixture[1], max);
#else
					p_site_mixture[0] = _mm_div_pX(p_site_mixture[0], max);
					p_site_mixture[1] = _mm_div_pX(p_site_mixture[1], max);
#endif
				}
			}
#else
        // iterate over all mixture categories
        for (size_t site = 0; site < numPatterns ; ++site)
        {   
            // the max probability
            RealNumber max = 0.0;

            // compute the per site probabilities
            for (size_t freq = 0; freq < numSiteFrequencies; ++freq)
            {
				for (size_t rate = 0; rate < numSiteRates; ++rate)
				{
					// get the pointers to the likelihood for this mixture category
					size_t offset = rate*rateOffset + freq*mixtureOffset + site*siteOffset;

					RealVector::iterator          p_site_mixture          = p_node + offset;

					for ( size_t i=0; i<2; ++i)
					{
						if ( p_site_mixture[i] > max )
						{
							max = p_site_mixture[i];
						}
					}

				}
            }

            p_scaler[site] = p_scaler_left[site] + p_scaler_right[site] + log(max);


            // compute the per site probabilities
            for (size_t freq = 0; freq < numSiteFrequencies; ++freq)
            {
				for (size_t rate = 0; rate < numSiteRates; ++rate)
				{
					// get the pointers to the likelihood for this mixture category
					size_t offset = rate*rateOffset + freq*mixtureOffset + site*siteOffset;

					RealVector::iterator          p_site_mixture          = p_node + offset;

					for ( size_t i=0; i<2; ++i)
					{
						p_site_mixture[i] /= max;
					}

				}
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
#ifdef AVX_ENABLED
            p_scaler[site] = _mm256_add_ps(p_scaler_left[site],p_scaler_right[site]);
#else
            p_scaler[site] = _mm_add_pX(p_scaler_left[site],p_scaler_right[site]);
#endif
        }
#else
        for (size_t site = 0; site < numPatterns ; ++site)
        {
            p_scaler[site] = p_scaler_left[site] + p_scaler_right[site];
        }
#endif

    }
}


void BinarySubstitutionModel::scale( size_t nodeIndex, size_t left, size_t right, size_t middle )
{
    RealVector::iterator p_node = partialLikelihoods.begin() + activeLikelihood[nodeIndex]*activeLikelihoodOffset + nodeIndex*nodeOffset;
    
#ifdef SIMD_ENABLED
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
            for (size_t freq = 0; freq < numSiteFrequencies; ++freq)
            {
				for (size_t rate = 0; rate < numSiteRates; ++rate)
				{
					// get the pointers to the likelihood for this mixture category
					size_t offset = rate*rateOffset + freq*mixtureOffset + site*siteOffset;

					SIMDRegister *          p_site_mixture = (SIMDRegister*)&*(p_node + offset);
#ifdef AVX_ENABLED
					SIMDRegister tmp = _mm256_max_ps(p_site_mixture[0],p_site_mixture[1]);
					if(freq == 0 && rate == 0)
						max = tmp;
					else
						max = _mm256_max_ps(max, tmp);
#else
					SIMDRegister tmp = _mm_max_pX(p_site_mixture[0],p_site_mixture[1]);
					if(freq == 0 && rate == 0)
						max = tmp;
					else
						max = _mm_max_pX(max, tmp);
#endif
				}
            }

#ifdef AVX_ENABLED
            p_scaler[site] = _mm256_add_ps(p_scaler_left[site], p_scaler_right[site]);
            p_scaler[site] = _mm256_add_ps(p_scaler[site], p_scaler_middle[site]);
#ifdef DOUBLE_PRECISION
            SIMDRegister tmp = max;
            double* ptmp = (double*)&tmp;
            double* pmax = (double*)&max;
            for(size_t i = 0; i < REALS_PER_SIMD_REGISTER; i++)
            	ptmp[i] = log(pmax[i]);

            p_scaler[site] = _mm256_add_pX(p_scaler[site], tmp);
#else
            p_scaler[site] = _mm256_add_pX(p_scaler[site], log256_ps(max));
#endif
#else
            p_scaler[site] = _mm_add_pX(p_scaler_left[site], p_scaler_right[site]);
            p_scaler[site] = _mm_add_pX(p_scaler[site], p_scaler_middle[site]);
#ifdef DOUBLE_PRECISION
            SIMDRegister tmp = max;
            double* ptmp = (double*)&tmp;
            double* pmax = (double*)&max;
            for(size_t i = 0; i < REALS_PER_SIMD_REGISTER; i++)
            	ptmp[i] = log(pmax[i]);

            p_scaler[site] = _mm_add_pX(p_scaler[site], tmp);
#else
            p_scaler[site] = _mm_add_pX(p_scaler[site], log_ps(max));
#endif
#endif

            // compute the per site probabilities
            for (size_t freq = 0; freq < numSiteFrequencies; ++freq)
            {
				for (size_t rate = 0; rate < numSiteRates; ++rate)
				{
					// get the pointers to the likelihood for this mixture category
					size_t offset = rate*rateOffset + freq*mixtureOffset + site*siteOffset;

					SIMDRegister *          p_site_mixture = (SIMDRegister*)&*(p_node + offset);
#ifdef AVX_ENABLED
					p_site_mixture[0] = _mm256_div_ps(p_site_mixture[0], max);
					p_site_mixture[1] = _mm256_div_ps(p_site_mixture[1], max);
#else
					p_site_mixture[0] = _mm_div_pX(p_site_mixture[0], max);
					p_site_mixture[1] = _mm_div_pX(p_site_mixture[1], max);
#endif
				}
            }
#else
        // iterate over all mixture categories
        for (size_t site = 0; site < numPatterns ; ++site)
        {   
            // the max probability
            RealNumber max = 0.0;

            // compute the per site probabilities
            for (size_t freq = 0; freq < numSiteFrequencies; ++freq)
            {
				for (size_t rate = 0; rate < numSiteRates; ++rate)
				{
					// get the pointers to the likelihood for this mixture category
					size_t offset = rate*rateOffset + freq*mixtureOffset + site*siteOffset;

					RealVector::iterator          p_site_mixture          = p_node + offset;

					for ( size_t i=0; i<2; ++i)
					{
						if ( p_site_mixture[i] > max )
						{
							max = p_site_mixture[i];
						}
					}

				}
            }

            p_scaler[site] = p_scaler_left[site] + p_scaler_right[site] + p_scaler_middle[site] + log(max);


            // compute the per site probabilities
            for (size_t freq = 0; freq < numSiteFrequencies; ++freq)
            {
				for (size_t rate = 0; rate < numSiteRates; ++rate)
				{
					// get the pointers to the likelihood for this mixture category
					size_t offset = rate*rateOffset + freq*mixtureOffset + site*siteOffset;

					RealVector::iterator          p_site_mixture          = p_node + offset;

					for ( size_t i=0; i<2; ++i)
					{
						p_site_mixture[i] /= max;
					}

				}
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
#ifdef AVX_ENABLED
            p_scaler[site] = _mm256_add_ps(p_scaler_left[site], p_scaler_right[site]);
            p_scaler[site] = _mm256_add_ps(p_scaler[site], p_scaler_middle[site]);
#else
            p_scaler[site] = _mm_add_pX(p_scaler_left[site], p_scaler_right[site]);
            p_scaler[site] = _mm_add_pX(p_scaler[site], p_scaler_middle[site]);
#endif
        }
#else
        for (size_t site = 0; site < numPatterns ; ++site)
        {
            p_scaler[site] = p_scaler_left[site] + p_scaler_right[site] + p_scaler_middle[site];
        }
#endif
    }
}



void BinarySubstitutionModel::scaleCorrection( size_t nodeIndex, size_t left, size_t right )
{

    RealVector::iterator p_node   = correctionLikelihoods.begin() + activeLikelihood[nodeIndex]*activeCorrectionOffset + nodeIndex*correctionNodeOffset;

               RealVector::iterator p_scaler  = perNodeCorrectionLogScalingFactors.begin() + activeLikelihood[nodeIndex]*activeCorrectionScalingOffset + nodeIndex*numCorrectionMasks;
    RealVector::const_iterator p_scaler_left  = perNodeCorrectionLogScalingFactors.begin() + activeLikelihood[left]*activeCorrectionScalingOffset      + left*numCorrectionMasks;
    RealVector::const_iterator p_scaler_right = perNodeCorrectionLogScalingFactors.begin() + activeLikelihood[right]*activeCorrectionScalingOffset     + right*numCorrectionMasks;
    
    if ( useScaling == true && nodeIndex % scalingDensity == 0 )
    {   
/*#ifdef SIMD_ENABLED
        for (size_t mask = 0; mask < numCorrectionMasks; ++mask)
        {
            
            SIMDRegister max;
            
            for (size_t freq = 0; freq < numSiteFrequencies; ++freq)
            {
				for (size_t rate = 0; rate < numSiteRates; ++rate)
				{
					size_t offset = rate*correctionRateOffset + freq*correctionMixtureOffset + mask*correctionMaskOffset;

					SIMDRegister *   u_i  = (SIMDRegister *)&*(p_node  + offset);
#ifdef AVX_ENABLED
					if(freq == 0 && rate == 0)
						max = *u_i;
					else
						max = _mm256_max_ps(max, *u_i);
#else
					SIMDRegister tmp = _mm_max_pX(u_i[0],u_i[1]);
					if(freq == 0 && rate == 0)
						max = tmp;
					else
						max = _mm_max_pX(max, tmp);
#endif
				}
            }
            
            RealNumber maximum = 0.0;
            
            RealNumber* tmp = (RealNumber*)&max;

            for(size_t i = 0; i < REALS_PER_SIMD_REGISTER; i++)
                maximum = std::max(maximum, tmp[i]);

            p_scaler[mask] = p_scaler_left[mask] + p_scaler_right[mask] + log(maximum);
            
#ifdef AVX_ENABLED
            max = _mm256_broadcast_ss(&maximum);
#else
            max = _mm_load1_pX(&maximum);
#endif
            
            // compute the per site probabilities
            for (size_t freq = 0; freq < numSiteFrequencies; ++freq)
            {
				for (size_t rate = 0; rate < numSiteRates; ++rate)
				{
					size_t offset = rate*correctionRateOffset + freq*correctionMixtureOffset + mask*correctionMaskOffset;

					SIMDRegister *   u_i  = (SIMDRegister *)&*(p_node  + offset);

#ifdef AVX_ENABLED
					*u_i = _mm256_div_ps(*u_i, max);
#else
					u_i[0] = _mm_div_pX(u_i[0], max);
					u_i[1] = _mm_div_pX(u_i[1], max);
#endif
				}
            }
        }
#else*/
        for (size_t mask = 0; mask < numCorrectionMasks; ++mask)
		{

			RealNumber max = 0.0;

			for (size_t freq = 0; freq < numSiteFrequencies; ++freq)
			{
				for (size_t rate = 0; rate < numSiteRates; ++rate)
				{
					size_t offset = rate*correctionRateOffset + freq*correctionMixtureOffset + mask*correctionMaskOffset;

					RealVector::const_iterator   u_i  = p_node  + offset;

					for(size_t i = 0; i < correctionMaskOffset; i++)
						max = std::max(max, u_i[i]);
				}
			}

			p_scaler[mask] = p_scaler_left[mask] + p_scaler_right[mask] + log(max);


			// compute the per site probabilities
			for (size_t freq = 0; freq < numSiteFrequencies; ++freq)
			{
				for (size_t rate = 0; rate < numSiteRates; ++rate)
				{
					size_t offset = rate*correctionRateOffset + freq*correctionMixtureOffset + mask*correctionMaskOffset;

					RealVector::iterator   u_i  = p_node  + offset;

					for(size_t i = 0; i < correctionMaskOffset; i++)
						u_i[i] /= max;
				}
			}
		}
//#endif
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


void BinarySubstitutionModel::scaleCorrection( size_t nodeIndex, size_t left, size_t right, size_t middle )
{
    RealVector::iterator p_node   = correctionLikelihoods.begin() + activeLikelihood[nodeIndex]*activeCorrectionOffset + nodeIndex*correctionNodeOffset;

               RealVector::iterator p_scaler  = perNodeCorrectionLogScalingFactors.begin() + activeLikelihood[nodeIndex]*activeCorrectionScalingOffset + nodeIndex*numCorrectionMasks;
    RealVector::const_iterator p_scaler_left  = perNodeCorrectionLogScalingFactors.begin() + activeLikelihood[left]*activeCorrectionScalingOffset      + left*numCorrectionMasks;
    RealVector::const_iterator p_scaler_right = perNodeCorrectionLogScalingFactors.begin() + activeLikelihood[right]*activeCorrectionScalingOffset     + right*numCorrectionMasks;
   RealVector::const_iterator p_scaler_middle = perNodeCorrectionLogScalingFactors.begin() + activeLikelihood[middle]*activeCorrectionScalingOffset    + middle*numCorrectionMasks;
    
    if ( useScaling == true && nodeIndex % scalingDensity == 0 )
    {
    	/*
#ifdef SIMD_ENABLED
        for (size_t mask = 0; mask < numCorrectionMasks; ++mask)
        {
            
            SIMDRegister max;
            
            for (size_t freq = 0; freq < numSiteFrequencies; ++freq)
            {
				for (size_t rate = 0; rate < numSiteRates; ++rate)
				{
					size_t offset = rate*correctionRateOffset + freq*correctionMixtureOffset + mask*8;

					SIMDRegister *   u_i  = (SIMDRegister *)&*(p_node  + offset);
#ifdef AVX_ENABLED
					if(freq == 0 && rate == 0)
						max = *u_i;
					else
						max = _mm256_max_ps(max, *u_i);
#else
					SIMDRegister tmp = _mm_max_pX(u_i[0],u_i[1]);
					if(freq == 0 && rate == 0)
						max = tmp;
					else
						max = _mm_max_pX(max, tmp);
#endif
				}
            }
            
            RealNumber maximum = 0.0;
            
            RealNumber* tmp = (RealNumber*)&max;

            for(size_t i = 0; i < REALS_PER_SIMD_REGISTER; i++)
                maximum = std::max(maximum, tmp[i]);

            p_scaler[mask] = p_scaler_left[mask] + p_scaler_right[mask] + p_scaler_middle[mask] + log(maximum);
            
#ifdef AVX_ENABLED
            max = _mm256_broadcast_ss(&maximum);
#else
            max = _mm_load1_pX(&maximum);
#endif
            
            // compute the per site probabilities
            for (size_t freq = 0; freq < numSiteFrequencies; ++freq)
            {
				for (size_t rate = 0; rate < numSiteRates; ++rate)
				{
					size_t offset = rate*correctionRateOffset + freq*correctionMixtureOffset + mask*8;

					SIMDRegister *   u_i  = (SIMDRegister *)&*(p_node  + offset);

#ifdef AVX_ENABLED
					*u_i = _mm256_div_ps(*u_i, max);
#else
					u_i[0] = _mm_div_pX(u_i[0], max);
					u_i[1] = _mm_div_pX(u_i[1], max);
#endif
				}
            }
        }
#else*/
        for (size_t mask = 0; mask < numCorrectionMasks; ++mask)
		{

			RealNumber max = 0.0;

			for (size_t freq = 0; freq < numSiteFrequencies; ++freq)
			{
				for (size_t rate = 0; rate < numSiteRates; ++rate)
				{
					size_t offset = rate*correctionRateOffset + freq*correctionMixtureOffset + mask*correctionMaskOffset;

					RealVector::const_iterator   u_i  = p_node  + offset;

					for(size_t i = 0; i < correctionMaskOffset; i++)
						max = std::max(max, u_i[i]);
				}
			}

			p_scaler[mask] = p_scaler_left[mask] + p_scaler_right[mask] + p_scaler_middle[mask] + log(max);


			// compute the per site probabilities
			for (size_t freq = 0; freq < numSiteFrequencies; ++freq)
			{
				for (size_t rate = 0; rate < numSiteRates; ++rate)
				{
					size_t offset = rate*correctionRateOffset + freq*correctionMixtureOffset + mask*correctionMaskOffset;

					RealVector::iterator   u_i  = p_node  + offset;

					for(size_t i = 0; i < correctionMaskOffset; i++)
						u_i[i] /= max;
				}
			}
		}
//#endif
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


size_t BinarySubstitutionModel::getNumSites() const {
    
    return numSites;
}


size_t BinarySubstitutionModel::getNumPatterns() const {
    
    return numPatterns;
}


RealNumber BinarySubstitutionModel::getLnCorrection() const {
    
    return perMaskCorrections[0];
}

std::vector<RealNumber> BinarySubstitutionModel::getLnCorrections() const {

    return perCodingProbs;
}


RealVector BinarySubstitutionModel::getPerSiteLnProbs() const {
    
    RealVector ret;
    
    for(size_t i = 0; i < numPatterns; i++)
    {
        if(coding != AscertainmentBias::ALL)
            ret.push_back((per_site_Likelihoods[i] - perMaskCorrections[pattern2site[site2mask[i]]])*patternCounts[i]);
        else
            ret.push_back(per_site_Likelihoods[i]*patternCounts[i]);
    }
    
    return ret;
}


void BinarySubstitutionModel::setVerbose(bool v) {
    
    verbose = v;
}

void BinarySubstitutionModel::setUseScaling(bool s) {
    
    if(s != useScaling)
    {
        useScaling = s;
        this->resizeLikelihoodVectors();
    }
}


void BinarySubstitutionModel::updateTransitionProbabilities() {

    for(std::set<size_t>::iterator it = touchedNodes.begin(); it != touchedNodes.end(); it++)
    {
        size_t nodeIndex = *it;
        
        // skip root
        if(nodeIndex == numNodes - 1)
            continue;
        
        RealNumber rate = getClockRate(nodeIndex);
        
        RealNumber brlen = getBranchLength(nodeIndex);
    
        RealVector::iterator p_node = transitionProbabilities.begin() + activeProbability[nodeIndex]*tActiveOffset + nodeIndex*tNodeOffset;
            
        for (size_t freq = 0; freq < numSiteFrequencies; ++freq)
        {
        	RealNumber pi1 = getStationaryFrequency(nodeIndex, freq);
			RealNumber pi0 = 1.0 - pi1;

			RealNumber mu = 1.0/(2.0*pi1*pi0);

			for (size_t rate = 0; rate < numSiteRates; ++rate)
			{
				RealNumber r = 1.0;
				if(rateVariationAcrossSites)
					r = siteRates->getValue()[rate];

				RealVector::iterator p_node_mixture = p_node + rate*tRateOffset + freq*tMixtureOffset;

    			RealNumber expPart = exp( - mu * rate * brlen * r );

        		p_node_mixture[0] = pi0 + pi1 * expPart;
				p_node_mixture[1] = pi1 - pi1 * expPart;
				p_node_mixture[2] = pi0 - pi0 * expPart;
				p_node_mixture[3] = pi1 + pi0 * expPart;
        	}
        }
    }

}




void BinarySubstitutionModel::keepSpecialization( DagNode* affecter ) {
    // reset all flags
    
    touchedNodes.clear();
    changedNodes.clear();
    
    dirtyNodes = std::vector<bool>(numNodes, false);
    
}


void BinarySubstitutionModel::restoreSpecialization( DagNode* affecter ) {
    
    for(std::set<size_t>::iterator it = touchedNodes.begin(); it != touchedNodes.end(); it++)
        activeProbability[*it] = (activeProbability[*it] ? 0 : 1);
    
    for(std::set<size_t>::iterator it = changedNodes.begin(); it != changedNodes.end(); it++)
    {
        // set all flags to false
        dirtyNodes[*it] = false;
        activeLikelihood[*it] = (activeLikelihood[*it] ? 0 : 1);
    }
    
    touchedNodes.clear();
    changedNodes.clear();
    
}


void BinarySubstitutionModel::touchSpecialization( DagNode* affecter) {
    
    // if the topology wasn't the culprit for the touch, then we just flag everything as dirty
    std::set<size_t> indices;

    if ( affecter == heterogeneousClockRates )
    {
        indices = heterogeneousClockRates->getTouchedElementIndices();
    }
    else if ( affecter == heterogeneousFrequencies )
    { 
        indices = heterogeneousFrequencies->getTouchedElementIndices();
    }
    else if ( affecter == branchLengths )
    { 
        indices = branchLengths->getTouchedElementIndices();
        
        const std::vector<TopologyNode*> nodes = tau->getValue().getNodes();
        
        if(indices.empty())
        {
        	for (size_t index = 0; index < dirtyNodes.size(); ++index)
			{
				if(!nodes[index]->isRoot())
					nodes[index]->setBranchLength(branchLengths->getValue()[index]);

				activeProbability[index] = (activeProbability[index] ? 0 : 1);
				touchedNodes.insert(index);
			}
        }
        else
        {
        	for (std::set<size_t>::iterator it = indices.begin(); it != indices.end(); ++it)
			{
				nodes[*it]->setBranchLength(branchLengths->getValue()[*it]);

				if ( touchedNodes.find(*it) == touchedNodes.end() )
				{
					activeProbability[*it] = (activeProbability[*it] ? 0 : 1);
					touchedNodes.insert(*it);
				}
			}
            
            indices.clear();
        }
        
    }
    
    if ( indices.size() != 0 )
    {
    	const std::vector<TopologyNode*> nodes = tau->getValue().getNodes();

        // flag recomputation only for the nodes
        for (std::set<size_t>::iterator it = indices.begin(); it != indices.end(); ++it)
        {
        	if ( touchedNodes.find(*it) == touchedNodes.end() )
            {
                activeProbability[*it] = (activeProbability[*it] ? 0 : 1);
                touchedNodes.insert(*it);
            }
            
            recursivelyFlagNodeDirty( *nodes[*it] );
        }
    }
    else if ( affecter != tau && affecter != branchLengths)
    {
        // flip the active likelihood pointers
        for (size_t index = 0; index < dirtyNodes.size(); ++index) 
        {
            dirtyNodes[index] = true;
            
            if ( touchedNodes.find(index) == touchedNodes.end() ) 
            {
                activeProbability[index] = (activeProbability[index] ? 0 : 1);
                touchedNodes.insert(index);
            }
            
            if ( changedNodes.find(index) == changedNodes.end() ) 
            {
                activeLikelihood[index] = (activeLikelihood[index] ? 0 : 1);
                changedNodes.insert(index);
            }
        }
    }
    
    if(this->dagNode != NULL)
        this->dagNode->touchAffected();
}


void BinarySubstitutionModel::fireTreeChangeEvent( const TopologyNode &n ) {
    
    // call a recursive flagging of all node above (closer to the root) and including this node
    recursivelyFlagNodeDirty( n );
}



void BinarySubstitutionModel::recursivelyFlagNodeDirty( const TopologyNode &n ) {
    
    // we need to flag this node and all ancestral nodes for recomputation
    size_t index = n.getIndex();

    // if this node is already dirty, the also all the ancestral nodes must have been flagged as dirty
    if ( !dirtyNodes[index] ) 
    {
        // the root doesn't have an ancestor
        if ( !n.isRoot() ) 
        {
            recursivelyFlagNodeDirty( n.getParent() );
        }
        
        // set the flag
        dirtyNodes[index] = true;
        
        // if we previously haven't touched this node, then we need to change the active likelihood pointer
        if ( changedNodes.find(index) == changedNodes.end() ) 
        {
            activeLikelihood[index] = (activeLikelihood[index] ? 0 : 1);
            changedNodes.insert(index);
        }
        
    }
    
}



void BinarySubstitutionModel::setValue(BinaryCharacterData *v) {
    
    // delegate to the parent class
    TypedDistribution< BinaryCharacterData >::setValue(v);
    
    continuous = v->getDatatype() == "Continuous";
    
    compress();
}




void BinarySubstitutionModel::swapParameter(const DagNode *oldP, const DagNode *newP) {
    
    if (oldP == homogeneousClockRate)
    {
        homogeneousClockRate = static_cast<const TypedDagNode< double >* >( newP );
    }
    else if (oldP == heterogeneousClockRates)
    {
        heterogeneousClockRates = static_cast<const TypedDagNode< std::vector< double > >* >( newP );
    }
    else if (oldP == branchLengths)
    {
        branchLengths = static_cast<const TypedDagNode< std::vector< double > >* >( newP );
    }
    else if (oldP == heterogeneousFrequencies)
    {
        heterogeneousFrequencies = static_cast<const TypedDagNode< std::vector<double> >* >( newP );
    }
    else if (oldP == homogeneousFrequency)
    {
        homogeneousFrequency = static_cast<const TypedDagNode< double >* >( newP );
    }
    else if (oldP == siteRates)
    {
        siteRates = static_cast<const TypedDagNode< std::vector< double > >* >( newP );
    }
    else if (oldP == samplingRate)
	{
    	samplingRate = static_cast<const TypedDagNode< double >* >( newP );
	}
    else if (oldP == tau)
    {
        tau->getValue().getTreeChangeEventHandler().removeListener( this );

        tau = static_cast<const TypedDagNode<Tree>* >( newP );

        tau->getValue().getTreeChangeEventHandler().addListener( this );

        numNodes = tau->getValue().getNumberOfNodes();
    }
    
}


void BinarySubstitutionModel::setClockRate(const TypedDagNode< double > *r)
{

    // remove the old parameter first
    if ( homogeneousClockRate != NULL )
    {
        this->removeParameter( homogeneousClockRate );
        homogeneousClockRate = NULL;
    }
    else if ( heterogeneousClockRates != NULL )
    {
        this->removeParameter( heterogeneousClockRates );
        heterogeneousClockRates = NULL;
    }

    // set the value
    branchHeterogeneousClockRates = false;
    homogeneousClockRate = r;

    // add the new parameter
    this->addParameter( homogeneousClockRate );

    // redraw the current value
    if ( this->dagNode == NULL || !this->dagNode->isClamped() )
    {
        this->redrawValue();
    }

}


void BinarySubstitutionModel::setClockRate(const TypedDagNode< std::vector< double > > *r)
{

    // remove the old parameter first
    if ( homogeneousClockRate != NULL )
    {
        this->removeParameter( homogeneousClockRate );
        homogeneousClockRate = NULL;
    }
    else if ( heterogeneousClockRates != NULL )
    {
        this->removeParameter( heterogeneousClockRates );
        heterogeneousClockRates = NULL;
    }

    // set the value
    branchHeterogeneousClockRates = true;
    heterogeneousClockRates = r;

    // add the new parameter
    this->addParameter( heterogeneousClockRates );

    // redraw the current value
    if ( this->dagNode == NULL || !this->dagNode->isClamped() )
    {
        this->redrawValue();
    }

}


void BinarySubstitutionModel::setBranchLengths(const TypedDagNode< std::vector< double > > *r)
{

    // remove the old parameter first
    if ( branchLengths != NULL )
    {
        this->removeParameter( branchLengths );
        branchLengths = NULL;
    }

    // set the value
    branchLengths = r;

    // add the new parameter
    this->addParameter( branchLengths );

    // redraw the current value
    if ( this->dagNode == NULL || !this->dagNode->isClamped() )
    {
        this->redrawValue();
    }

}


void BinarySubstitutionModel::setSiteRates(const TypedDagNode< std::vector< double > > *r) {
    
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
        numSiteRates = r->getValue().size();
        resizeLikelihoodVectors();
    }
    else
    {
        // set the value
        rateVariationAcrossSites = false;
        siteRates = NULL;
        numSiteRates = 1;
        resizeLikelihoodVectors();
        
    }
    
    // add the parameter
    this->addParameter( siteRates );

    // redraw the current value
    if ( this->dagNode != NULL && !this->dagNode->isClamped() )
    {
        redrawValue();
    }
}


double BinarySubstitutionModel::getStationaryFrequency( size_t nodeIndex, size_t mixture ) {
    
	if ( frequencyVariationAcrossSites)
	{
		return heterogeneousFrequencies->getValue()[mixture];
	}
	else if ( branchHeterogeneousFrequencies )
    {
		return heterogeneousFrequencies->getValue()[nodeIndex];
    } 
    else if ( homogeneousFrequency != NULL ) 
    {
        return homogeneousFrequency->getValue();
    }
    else
    {
        return 0.5;
    }
}

double BinarySubstitutionModel::getSamplingRate(void ) {

	if ( samplingRate)
	{
		return samplingRate->getValue();
	}
    else
    {
        return 1.0;
    }
}


double BinarySubstitutionModel::getClockRate( size_t nodeIndex ) {
    
    // second, get the clock rate for the branch
    double rate = 1.0;
    
    if ( branchHeterogeneousClockRates == true )
    {
    	rate = heterogeneousClockRates->getValue()[nodeIndex];
    }
    else if(homogeneousClockRate != NULL)
    {
        rate = homogeneousClockRate->getValue();
    }
    
    return rate;
}


double BinarySubstitutionModel::getBranchLength( size_t nodeIndex ) {
    
	// second, get the clock rate for the branch
    double brlen = 1.0;
    
    if ( branchLengths != NULL )
    {
        brlen = branchLengths->getValue()[nodeIndex];
    }
    else
    {
        brlen = tau->getValue().getNode(nodeIndex).getBranchLength();
    }
    
    return brlen;
}


void BinarySubstitutionModel::setStationaryFrequency(const TypedDagNode< double > *r)
{

    // remove the old parameter first
    if ( homogeneousFrequency != NULL )
    {
        this->removeParameter( homogeneousFrequency );
        homogeneousFrequency = NULL;
    }
    else if ( heterogeneousFrequencies != NULL )
    {
        this->removeParameter( heterogeneousFrequencies );
        heterogeneousFrequencies = NULL;
    }

    // set the value
    branchHeterogeneousFrequencies = false;
    frequencyVariationAcrossSites = false;
    homogeneousFrequency = r;

    // add the new parameter
    this->addParameter( homogeneousFrequency );

    // redraw the current value
    if ( this->dagNode == NULL || !this->dagNode->isClamped() )
    {
        this->redrawValue();
    }

}


void BinarySubstitutionModel::setStationaryFrequency(const TypedDagNode< std::vector< double > > *f) {

    // remove the old parameter first
    if ( homogeneousFrequency != NULL )
    {
        this->removeParameter( homogeneousFrequency );
        homogeneousFrequency = NULL;
    }
    else if ( heterogeneousFrequencies != NULL )
    {
        this->removeParameter( heterogeneousFrequencies );
        heterogeneousFrequencies = NULL;
    }

    // set the value
    branchHeterogeneousFrequencies = true;
    frequencyVariationAcrossSites = false;
    heterogeneousFrequencies = f;

    // add the new parameter
    this->addParameter( heterogeneousFrequencies );

    // redraw the current value
    if ( this->dagNode == NULL || !this->dagNode->isClamped() )
    {
        this->redrawValue();
    }

}

void BinarySubstitutionModel::setSiteFrequencies(const TypedDagNode< std::vector< double > > *r) {

    // remove the old parameter first
    if ( heterogeneousFrequencies != NULL )
    {
        this->removeParameter( heterogeneousFrequencies );
        heterogeneousFrequencies = NULL;
    }

    if ( r != NULL )
    {
        // set the value
        frequencyVariationAcrossSites = true;
        heterogeneousFrequencies = r;
        numSiteFrequencies = r->getValue().size();
        resizeLikelihoodVectors();
    }
    else
    {
        // set the value
    	frequencyVariationAcrossSites = false;
    	heterogeneousFrequencies = NULL;
    	numSiteFrequencies = 1;
        resizeLikelihoodVectors();

    }

    // add the parameter
    this->addParameter( heterogeneousFrequencies );

    // redraw the current value
    if ( this->dagNode != NULL && !this->dagNode->isClamped() )
    {
        redrawValue();
    }
}

void BinarySubstitutionModel::setSamplingRate(const TypedDagNode< double > *r)
{

    // remove the old parameter first
    if ( samplingRate != NULL )
    {
        this->removeParameter( samplingRate );
        samplingRate = NULL;
    }

    // set the value
    samplingRate = r;

    // add the new parameter
    this->addParameter( samplingRate );

    // redraw the current value
    if ( this->dagNode == NULL || !this->dagNode->isClamped() )
    {
        this->redrawValue();
    }

}

std::vector<int> BinarySubstitutionModel::getCountDistribution( void ) const
{
	return countDistribution;
}

void BinarySubstitutionModel::redrawValue( void ) {

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

    const TopologyNode &root = tau->getValue().getRoot();
    size_t rootIndex = tau->getValue().getRoot().getIndex();
    
    std::vector< BinaryTaxonData* > taxa;

    for(size_t i = 0; i < numTaxa; i++)
    {
		if(continuous)
			taxa.push_back( new ContinuousBinaryTaxonData("") );
		else
			taxa.push_back( new DiscreteBinaryTaxonData("") );
    }
    
    /*double total = 1.0;

    if(coding != AscertainmentBias::ALL && numSites > 0)
    {
        // first sample a total number of characters (M) from the marginal posterior: 
        // M - N | N ~ NegBinomial(N, exp(lnCorrection) )
        //double M_minus_N = Statistics::NegativeBinomial::rv(N, exp(perMaskCorrections[0]), *rng);
    
        // then sample the observed number of characters (numSites) from the likelihood:
        // numSites | M ~ Binomial(M, exp(lnCorrection) )
        //numSites = Statistics::Binomial::rv( M_minus_N + N, exp(perMaskCorrections[0]), *rng);
    
        // sample the rate categories in proportion to the total probability (correction) for each mixture.
        for ( size_t i = 0; i < numSiteFrequencies; ++i )
        	for ( size_t j = 0; j < numSiteRates; ++j )
        		total += perMixtureCorrections[i*numSiteRates*numCorrectionMasks + j*numCorrectionMasks];

        // draw the state
			double u = rng->uniform01()*total;

			double tmp = 0.0;
			while(tmp < u){
				tmp += perMixtureCorrections[freqIndex*numSiteRates*numCorrectionMasks + rateIndex*numCorrectionMasks];
				if(tmp < u)
				{
					rateIndex++;
					if(rateIndex == numSiteRates)
					{
						freqIndex++;
						rateIndex = 0;
					}
				}
			}
    }*/

    double sampling = getSamplingRate();
    std::fill(countDistribution.begin(), countDistribution.end(), 0);

    // then sample site-patterns using rejection sampling,
    // rejecting those that match the unobservable ones.
    for ( size_t i = 0; i < numSites; i++ )
    {
    	size_t rateIndex = 0;
    	size_t freqIndex = 0;

    	// draw the state
		double u = rng->uniform01();

		rateIndex = (int)(u*numSiteRates);

		u = rng->uniform01();
		freqIndex = (int)(u*numSiteFrequencies);


        std::vector<RealNumber> siteData(numNodes, 0.0);

        RealNumber rf = getStationaryFrequency(rootIndex, freqIndex);
        // draw the root state
        u = rng->uniform01();
        if(u < rf)
            siteData[rootIndex] = 1.0;

        // recursively simulate the sequences
        std::pair<size_t, size_t> charCounts;
        simulate( root, siteData, rateIndex, freqIndex, charCounts);

        countDistribution[charCounts.second]++;

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


void BinarySubstitutionModel::simulate( const TopologyNode &node, std::vector<RealNumber> &data, size_t rateIndex, size_t freqIndex, std::pair<size_t, size_t>& charCounts) {

    // get the children of the node
    const std::vector<TopologyNode*>& children = node.getChildren();

    // get the sequence of this node
    size_t nodeIndex = node.getIndex();
    RealNumber parentState = data[ nodeIndex ];

    // simulate the sequence for each child
    RandomNumberGenerator* rng = GLOBAL_RNG;
    for (std::vector< TopologyNode* >::const_iterator it = children.begin(); it != children.end(); ++it)
    {
        const TopologyNode &child = *(*it);

        RealVector::iterator pi = transitionProbabilities.begin() + activeProbability[nodeIndex]*tActiveOffset + nodeIndex*tNodeOffset + rateIndex*tRateOffset + freqIndex*tMixtureOffset;

        if(parentState == 1.0)
            pi += 2;
        
        // create the character
        data[ child.getIndex() ] = 0.0;
        // draw the state
        double u = rng->uniform01()*(pi[0] + pi[1]);
        if(u < pi[1])
            data[ child.getIndex() ] = 1.0;

        RealNumber& childState = data[ child.getIndex() ];
        if(child.isTip())
        {
            if(continuous)
                data[ child.getIndex() ] = pi[1]/(pi[0] + pi[1]);
            
            if(childState == 1.0)
                charCounts.second++;
            else if(childState == 0.0)
                charCounts.first++;
        }
        else
            simulate( child, data, rateIndex, freqIndex, charCounts);
    }

}



const std::vector< DiscreteBinaryTaxonData >& BinarySubstitutionModel::getMapping(void) {

    computeLnProbability();

    mapping = std::vector< DiscreteBinaryTaxonData >( tau->getValue().getNumberOfNodes(), DiscreteBinaryTaxonData() );

    // first, simulate the per site rates
    RandomNumberGenerator* rng = GLOBAL_RNG;

    std::vector<size_t> perSiteRates;
    std::vector<size_t> perSiteFrequencies;
    if(coding != AscertainmentBias::ALL)
	{
		double total = 0.0;
		for ( size_t i = 0; i < numSiteFrequencies; ++i )
			for ( size_t j = 0; j < numSiteRates; ++j )
				total += perMixtureCorrections[i*numSiteRates*numCorrectionMasks + j*numCorrectionMasks];

		for ( size_t i = 0; i < numSites; ++i )
		{
			// draw the state
			double u = rng->uniform01()*total;
			size_t rateIndex = 0;
			size_t freqIndex = 0;

			double tmp = 0.0;
			while(tmp < u){
				tmp += perMixtureCorrections[freqIndex*numSiteRates*numCorrectionMasks + rateIndex*numCorrectionMasks];
				if(tmp < u)
				{
					rateIndex++;
					if(rateIndex == numSiteRates)
					{
						freqIndex++;
						rateIndex = 0;
					}
				}
			}

			perSiteRates.push_back( rateIndex );
			perSiteFrequencies.push_back( freqIndex );
		}
	}
	else
	{
		for ( size_t i = 0; i < numSites; ++i )
		{
			// draw the state
			double u = rng->uniform01();

			size_t rateIndex = (int)(u*numSiteRates);
			perSiteRates.push_back( rateIndex );

			u = rng->uniform01();
			size_t freqIndex = (int)(u*numSiteFrequencies);

			perSiteFrequencies.push_back( rateIndex );
		}
	}
    
    // compute the ln probability by recursively calling the probability calculation for each node
    const TopologyNode &root = tau->getValue().getRoot();
    
    // we start with the root and then traverse down the tree
    size_t rootIndex = root.getIndex();

    DiscreteBinaryTaxonData &rootData = mapping[ rootIndex ];
#ifdef SIMD_ENABLED
    for ( size_t site = 0; site < numSites; ++site )
    {
        size_t offset = perSiteRates[site]*rateOffset + perSiteFrequencies[site]*mixtureOffset + (site2pattern[site] % REALS_PER_SIMD_REGISTER) + size_t(site2pattern[site]/REALS_PER_SIMD_REGISTER)*siteOffset;
                
        RealVector::const_iterator   p_node  = partialLikelihoods.begin() + activeLikelihood[rootIndex]*activeLikelihoodOffset + rootIndex*nodeOffset + offset;
        
        RealNumber p0 = p_node[0];
        RealNumber p1 = p_node[REALS_PER_SIMD_REGISTER];
#else
    for ( size_t site = 0; site < numSites; ++site )
    {
        size_t offset = perSiteRates[site]*rateOffset + perSiteFrequencies[site]*mixtureOffset + site2pattern[site]*siteOffset;

        RealVector::const_iterator   p_root  = partialLikelihoods.begin() + activeLikelihood[rootIndex]*activeLikelihoodOffset + rootIndex*nodeOffset + offset;

        RealNumber p0 = p_root[0];
        RealNumber p1 = p_root[1];
#endif    
        // create the character
        bool c = false;
        // draw the state
        double u = rng->uniform01()*(p0 + p1);

        if(u < p1)
            c = true;

        rootData.addCharacter(c);
    }

    // we start with the root and then traverse down the tree
    const std::vector<TopologyNode*> children = root.getChildren();
    
    for(std::vector<TopologyNode*>::const_iterator it = children.begin(); it != children.end(); it++)
    {
        simulateMapping(root, *(*it), perSiteRates, perSiteFrequencies);
    }
        
    return mapping;
}



void BinarySubstitutionModel::simulateMapping(const TopologyNode &parent, const TopologyNode &node, std::vector<size_t>& perSiteRates, std::vector<size_t>& perSiteFrequencies)
{
    // first, simulate the per site rates
    RandomNumberGenerator* rng = GLOBAL_RNG;

    size_t nodeIndex = node.getIndex();
    size_t parentIndex = parent.getIndex();
    
    DiscreteBinaryTaxonData &nodeData   = mapping[ nodeIndex ];
    DiscreteBinaryTaxonData &parentData = mapping[ parentIndex ];
    
#ifdef SIMD_ENABLED
    for ( size_t site = 0; site < numSites; ++site )
    {
        size_t offset = perSiteRates[site]*rateOffset + perSiteFrequencies[site]*mixtureOffset + (site2pattern[site] % REALS_PER_SIMD_REGISTER) + size_t(site2pattern[site]/REALS_PER_SIMD_REGISTER)*siteOffset;
        
        RealVector::const_iterator   p_node  = partialLikelihoods.begin() + activeLikelihood[nodeIndex]*activeLikelihoodOffset + nodeIndex*nodeOffset + offset;
        
        RealVector::iterator         pi_node = transitionProbabilities.begin() + activeProbability[nodeIndex]*tActiveOffset + nodeIndex*tNodeOffset;
        
        if(parentData[site])
            pi_node += 2;
        
        RealNumber p0 = pi_node[0]*p_node[0];
        RealNumber p1 = pi_node[1]*p_node[REALS_PER_SIMD_REGISTER];
#else
    for ( size_t site = 0; site < numSites; ++site )
    {
        size_t offset = perSiteRates[site]*rateOffset + perSiteFrequencies[site]*mixtureOffset + site2pattern[site]*siteOffset;

        RealVector::const_iterator   p_node  = partialLikelihoods.begin() + activeLikelihood[nodeIndex]*activeLikelihoodOffset + nodeIndex*nodeOffset + offset;
        
        RealVector::iterator         pi_node = transitionProbabilities.begin() + activeProbability[nodeIndex]*tActiveOffset + nodeIndex*tNodeOffset;
        
        if(parentData[site])
            pi_node += 2;

        RealNumber p0 = pi_node[0]*p_node[0];
        RealNumber p1 = pi_node[1]*p_node[1];
#endif
        // create the character
        bool c = false;
        // draw the state
        double u = rng->uniform01()*(p0 + p1);

        if(u < p1)
            c = true;

        nodeData.addCharacter(c);
    }

    // we start with the root and then traverse down the tree
    const std::vector<TopologyNode*> children = node.getChildren();
    
    for(std::vector<TopologyNode*>::const_iterator it = children.begin(); it != children.end(); it++)
    {
        simulateMapping(node, *(*it), perSiteRates, perSiteFrequencies);
    }
}

set(BIPHY_FILES
	${BIPHY_DIR}/Biphy.cpp
	${BIPHY_DIR}/BiphyFreeTopology.cpp
	${BIPHY_DIR}/MultiBiphy.cpp
	${BIPHY_DIR}/MultiBiphyRhoPi.cpp
	${BIPHY_DIR}/MultiBiphyLambdaMu.cpp
	${BIPHY_DIR}/CladeReader.cpp
	${BIPHY_DIR}/CrossValidationScoreMonitor.cpp
	${BIPHY_DIR}/DPPBetaSimplexMove.cpp
	${BIPHY_DIR}/FreeBinaryRateMatrixVectorFunction.cpp
	${BIPHY_DIR}/MultiMove.cpp
	${BIPHY_DIR}/NewickTreeMonitor.cpp
	${BIPHY_DIR}/NexusTreeMonitor.cpp
	${BIPHY_DIR}/PosteriorPredictiveStateFrequencyMonitor.cpp
	${BIPHY_DIR}/TruncatedDistributionUnnormalized.cpp
	${BIPHY_DIR}/UniformRootedTopologyDistribution.cpp
	${BIPHY_DIR}/NegativeBinomialDistribution.cpp
	${BIPHY_DIR}/FreeBinaryRateMatrixVectorFunction.cpp
	${BIPHY_DIR}/NHXTreeReader.cpp
	${BIPHY_DIR}/NHXConverter.cpp
	${BIPHY_DIR}/NewickTreeMonitor.cpp
	${BIPHY_DIR}/NexusTreeMonitor.cpp
)
add_library(biphy-lib ${BIPHY_FILES})
target_link_libraries(biphy-lib rb-core)

set(BIPHY_FILES
${BIPHY_DIR}/BinarySubstitutionModel.cpp
${BIPHY_DIR}/BinaryDolloSubstitutionModel.cpp
${BIPHY_DIR}/Biphy.cpp
${BIPHY_DIR}/CladeReader.cpp
${BIPHY_DIR}/MultiMove.cpp
${BIPHY_DIR}/NewickTreeMonitor.cpp
${BIPHY_DIR}/NexusTreeMonitor.cpp
${BIPHY_DIR}/TruncatedDistributionUnnormalized.cpp
${BIPHY_DIR}/UniformTopologyDistribution.cpp
${BIPHY_DIR}/NegativeBinomialDistribution.cpp
${BIPHY_DIR}/NewickTreeMonitor.cpp
${BIPHY_DIR}/NexusTreeMonitor.cpp
${BIPHY_DIR}/CrossValidationScoreMonitor.cpp
${BIPHY_DIR}/PosteriorPredictiveStateFrequencyMonitor.cpp
${BIPHY_DIR}/MappingMonitor.cpp
${BIPHY_DIR}/LnCorrectionFunction.cpp
)
add_library(biphy-lib ${BIPHY_FILES})
target_link_libraries(biphy-lib rb-core)

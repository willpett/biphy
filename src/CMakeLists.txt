set(SRC_FILES
${PROJECT_SRC_DIR}/BranchHeterogeneousBinaryModel.cpp
${PROJECT_SRC_DIR}/CladeReader.cpp
${PROJECT_SRC_DIR}/CrossValidationScoreMonitor.cpp
${PROJECT_SRC_DIR}/DPPBetaSimplexMove.cpp
${PROJECT_SRC_DIR}/FreeBinaryRateMatrixVectorFunction.cpp
${PROJECT_SRC_DIR}/ModelStreamMonitor.cpp
${PROJECT_SRC_DIR}/ModelStreamReader.cpp
${PROJECT_SRC_DIR}/MultiMove.cpp
${PROJECT_SRC_DIR}/NewickTreeMonitor.cpp
${PROJECT_SRC_DIR}/NexusTreeMonitor.cpp
${PROJECT_SRC_DIR}/PosteriorPredictiveStateFrequencyMonitor.cpp
${PROJECT_SRC_DIR}/TruncatedDistributionUnnormalized.cpp
${PROJECT_SRC_DIR}/UniformRootedTopologyDistribution.cpp
${PROJECT_SRC_DIR}/UniformConstrainedTimeTreeDistribution.cpp
${PROJECT_SRC_DIR}/VectorScaleFunction.cpp
)
add_library(wp-lib ${SRC_FILES})
target_link_libraries(wp-lib rb-core)

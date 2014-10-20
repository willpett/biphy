cmake_minimum_required(VERSION 2.6)
project(biphy)

find_package(OpenMP REQUIRED)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -g -msse -msse2 -msse3 ${OpenMP_CXX_FLAGS}")
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O3 -g ${OpenMP_C_FLAGS}")


set(REVBAYES_DIR ${CMAKE_SOURCE_DIR}/revbayes)
set(PROJECT_SRC_DIR ${CMAKE_SOURCE_DIR}/src)
# Split into much smaller libraries

add_subdirectory(revbayes)
include_directories(${RBCORE_INCLUDES})

add_subdirectory(src)
############# executables #################
# basic rev-bayes binary

set(TARGET_NAME "biphy")
add_executable(${TARGET_NAME} ${PROJECT_SRC_DIR}/main.cpp)
set_target_properties(${TARGET_NAME} PROPERTIES PREFIX "../")
target_link_libraries(${TARGET_NAME} rb-core wp-lib gomp)

#add_executable(convert ${PROJECT_SRC_DIR}/convert.cpp)
#set_target_properties(convert PROPERTIES PREFIX "../")
#target_link_libraries(convert rb-core wp-lib)

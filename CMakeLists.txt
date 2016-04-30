cmake_minimum_required(VERSION 2.6)
project(biphy)

find_package(Git)
if(GIT_FOUND)
    execute_process(
            COMMAND ${GIT_EXECUTABLE} describe --abbrev=4 --dirty --always --tags
            OUTPUT_VARIABLE DESCRIPTION
            OUTPUT_STRIP_TRAILING_WHITESPACE)
endif()
add_definitions(-DVERSION="${DESCRIPTION}" -DSINGLE_PRECISION)

find_package(OpenMP REQUIRED)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -msse3 -O3 ${OpenMP_CXX_FLAGS}")
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -msse3 -O3 ${OpenMP_C_FLAGS}")

#-Wa,-q -mavx 

set(SRC_DIR ${CMAKE_SOURCE_DIR}/src)

add_subdirectory(src)
include_directories(${SRC_INCLUDES})

set(TARGET_NAME "biphy")
add_executable(${TARGET_NAME} ${SRC_DIR}/main.cpp)
set_target_properties(${TARGET_NAME} PROPERTIES PREFIX "../")
target_link_libraries(${TARGET_NAME} biphy-core gomp)

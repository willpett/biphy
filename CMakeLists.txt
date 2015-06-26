cmake_minimum_required(VERSION 2.6)
project(biphy)

find_package(Git)
if(GIT_FOUND)
    execute_process(
            COMMAND ${GIT_EXECUTABLE} describe --abbrev=4 --dirty --always --tags
            OUTPUT_VARIABLE DESCRIPTION
            OUTPUT_STRIP_TRAILING_WHITESPACE)
endif()
add_definitions(-DVERSION="${DESCRIPTION}")

find_package(OpenMP REQUIRED)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -msse -msse2 -msse3 ${OpenMP_CXX_FLAGS}")
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O3 ${OpenMP_C_FLAGS}")


set(REVBAYES_DIR ${CMAKE_SOURCE_DIR}/revbayes)
set(BIPHY_DIR ${CMAKE_SOURCE_DIR}/src)

add_subdirectory(revbayes)
include_directories(${RBCORE_INCLUDES})

add_subdirectory(src)
include_directories(src)

set(TARGET_NAME "biphy")
add_executable(${TARGET_NAME} ${BIPHY_DIR}/main.cpp)
set_target_properties(${TARGET_NAME} PROPERTIES PREFIX "../")
target_link_libraries(${TARGET_NAME} rb-core biphy-lib gomp)

set(TARGET_NAME_MULTI "multibiphy")
add_executable(${TARGET_NAME_MULTI} ${BIPHY_DIR}/main_multi.cpp)
set_target_properties(${TARGET_NAME_MULTI} PROPERTIES PREFIX "../")
target_link_libraries(${TARGET_NAME_MULTI} rb-core biphy-lib gomp)

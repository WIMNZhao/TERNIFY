cmake_minimum_required(VERSION 3.30)
cmake_policy(SET CMP0167 NEW)
project(Ternify VERSION 1.0)


set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# set building type as Debug
#set(CMAKE_BUILD_TYPE Debug)

# set building type as Release
set(CMAKE_BUILD_TYPE Release)

# for Debug  with special tags
#set(CMAKE_CXX_FLAGS_DEBUG "-g -O0")

#if(NOT CMAKE_BUILD_TYPE)
#    set(CMAKE_BUILD_TYPE Release)
#endif()

# compile options
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra")
#set(CMAKE_CXX_FLAGS_DEBUG "-g -O0")
set(CMAKE_CXX_FLAGS_RELEASE "-O3")

# set rdkit, nlopt lib !!!!
set(CONDA_ENV "/opt/anaconda3/envs/ternify/")
link_directories(${CONDA_ENV}/lib)
set(RDKIT_DIR ${CONDA_ENV})
set(NLOPT_DIR ${CONDA_ENV})
find_package(RDKit REQUIRED)
find_package(Eigen3 REQUIRED)
find_package(Threads REQUIRED)
find_package(NLopt REQUIRED)

# source list
set(SOURCES
    src/align.cpp
    src/getgriden.cpp
    src/grid.cpp
    src/main.cpp
    src/minimize_h.cpp
    src/parameters.cpp
    src/protac.cpp
    src/protein.cpp
)

# header list
set(HEADERS
    src/align.hpp
    src/getgriden.hpp
    src/grid.hpp
    src/minimize_h.hpp
    src/parameters.hpp
    src/protac.hpp
    src/protein.hpp
    src/ternify.hpp
    src/common.hpp
    src/lowEnerDihe.hpp
)

# building executable file
add_executable(ternify ${SOURCES} ${HEADERS})

# including directories
target_include_directories(ternify
    PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/src
        ${RDKit_INCLUDE_DIRS}
        ${EIGEN3_INCLUDE_DIR}
        ${NLOPT_INCLUDE_DIRS} 
)

# link libraries
target_link_libraries(ternify
    PRIVATE
        ${RDKit_LIBRARIES}
        RDKitFileParsers
        RDKitForceField
        RDKitMolTransforms
        RDKitSubstructMatch
        RDKitForceFieldHelpers
        RDKitForceField
        RDKitSmilesParse
        RDKitGraphMol
        RDKitRDGeneral
        RDKitMolAlign
        RDKitPartialCharges
        RDKitRDGeometryLib
        RDKitDistGeomHelpers
        Eigen3::Eigen
        Threads::Threads
        nlopt
)

# installation rules
install(TARGETS ternify
    RUNTIME DESTINATION bin
)

# adding testing
enable_testing()
# add_subdirectory(tests)  # if there is any test

# print config
message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")
message(STATUS "RDKit version: ${RDKit_VERSION}")
message(STATUS "Eigen3 version: ${EIGEN3_VERSION}")

# set output directory
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

# if in debug type
if(CMAKE_BUILD_TYPE STREQUAL "Debug")
    target_compile_definitions(ternify PRIVATE DEBUG)
endif()

# check os
if(UNIX AND NOT APPLE)
    target_link_libraries(ternify PRIVATE stdc++fs)
endif()

# check verssion
configure_file(
    ${CMAKE_CURRENT_SOURCE_DIR}/version.h.in
    ${CMAKE_CURRENT_BINARY_DIR}/version.h
)

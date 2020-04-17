cmake_minimum_required(VERSION 3.15)

list(APPEND CMAKE_MODULE_PATH
    "${CMAKE_SOURCE_DIR}/cmake/"
    "${CMAKE_SOURCE_DIR}/cmake/Modules/"
)

if("${CMAKE_SOURCE_DIR}" STREQUAL "${CMAKE_BINARY_DIR}")
    message(FATAL_ERROR "In-source builds are not supported; \
generation should be done in a local subfolder, e.g. '${CMAKE_SOURCE_DIR}/build'.")
endif()

include(FrameworkUtils)
include(ExternalProject)
include(DownloadProject)
include(GetGitRevisionDescription)
include(CheckFunctionExists)

project(covid19_p2p_geo C CXX)

include(GNUInstallDirs)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/${CMAKE_INSTALL_LIBDIR})
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/${CMAKE_INSTALL_LIBDIR})
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/${CMAKE_INSTALL_BINDIR})
set(CMAKE_INSTALL_RPATH "$ORIGIN/../lib:$ORIGIN/")
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
set_property(GLOBAL PROPERTY USE_FOLDERS ON)

set(COVID19_P2P_GEO_VERSION_MAJOR 0)
set(COVID19_P2P_GEO_VERSION_MINOR 1)
set(COVID19_P2P_GEO_VERSION_PATCH 0)
set(COVID19_P2P_GEO_VERSION
    "${COVID19_P2P_GEO_VERSION_MAJOR}.${COVID19_P2P_GEO_VERSION_MINOR}.${COVID19_P2P_GEO_VERSION_PATCH}")
set(COVID19_P2P_GEO_VERSION_PLAIN
    "${COVID19_P2P_GEO_VERSION_MAJOR}${COVID19_P2P_GEO_VERSION_MINOR}${COVID19_P2P_GEO_VERSION_PATCH}")
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")

add_subdirectory(3rdparty)

add_definitions(-Wfatal-errors)

if(CMAKE_CXX_COMPILER_ID MATCHES GNU)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fprofile-arcs -ftest-coverage")
endif()

add_library(covid19_p2p_geo_api STATIC)
target_sources(
    covid19_p2p_geo_api
        PRIVATE
            "api/src/api.cpp"
            "api/include/generic_utils.hpp"
            "api/include/geo_utils.hpp"
            "api/include/hdf5_utils.hpp"
            "api/include/regions.hpp"
        PUBLIC
            "api/include/api.hpp"
)
add_dependencies(
    covid19_p2p_geo_api
        hdf5_local
)
target_include_directories(
    covid19_p2p_geo_api
        PUBLIC
            "${CMAKE_SOURCE_DIR}/api/include"
            "${_hdf5_include_dirs}"
)
target_link_libraries(
    covid19_p2p_geo_api
        PRIVATE
            geos
            "$<BUILD_INTERFACE:${_hdf5_libs}>"
            lz4
)

add_executable(
    covid19_p2p_geo_tester
        "apps/src/random_tester.cpp"
)
add_dependencies(
    covid19_p2p_geo_tester
        hdf5_local
)
target_include_directories(
    covid19_p2p_geo_tester
        PUBLIC
        "$<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/api/include>"
        "$<BUILD_INTERFACE:${_hdf5_include_dirs}>"
)
target_link_libraries(
    covid19_p2p_geo_tester
        PUBLIC
            geos
            "$<BUILD_INTERFACE:${_hdf5_libs}>"
            lz4
            covid19_p2p_geo_api
)

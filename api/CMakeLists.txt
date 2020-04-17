
add_library(covid19_p2p_geo_api STATIC)

target_sources(
    covid19_p2p_geo_api
        PRIVATE
            "src/api.cpp"
            "include/generic_utils.hpp"
            "include/geo_utils.hpp"
            "include/hdf5_utils.hpp"
            "include/regions.hpp"
        PUBLIC
            "include/api.hpp"
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


add_library(covid19_p2p_geo_api STATIC)

target_sources(
    covid19_p2p_geo_api
        PRIVATE
            "src/api.cpp"
            "include/private/generic_utils.hpp"
            "include/private/geo_utils.hpp"
            "include/private/hdf5_utils.hpp"
            "include/private/regions.hpp"
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

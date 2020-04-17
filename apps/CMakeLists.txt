
add_executable(
    covid19_p2p_geo_tester
        "src/random_tester.cpp"
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

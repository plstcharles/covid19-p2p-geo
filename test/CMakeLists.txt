
add_executable(covid19_p2p_geo_tests)

target_sources(
    covid19_p2p_geo_tests
        PUBLIC
            "src/test_api.cpp"
            "src/test_geo.cpp"
            "src/test_uids.cpp"
)

target_link_libraries(
    covid19_p2p_geo_tests
        PUBLIC
            geos
            "$<BUILD_INTERFACE:${_hdf5_libs}>"
            lz4
            covid19_p2p_geo_api
            gtest_main
)

target_compile_definitions(
    covid19_p2p_geo_tests
        PUBLIC
            -DDATA_ROOT="${DATA_ROOT}"
)

add_test(
    NAME
        covid19_p2p_geo
    COMMAND
        "${CMAKE_BINARY_DIR}/${CMAKE_INSTALL_BINDIR}/covid19_p2p_geo_tests"
)

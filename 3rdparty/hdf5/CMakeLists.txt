
# derived from https://github.com/bast/hdf5-cmake-example/blob/master/CMakeLists.txt

include(ExternalProject)
set(ExternalProjectCMakeArgs -DHDF5_BUILD_CPP_LIB=ON)
set(_hdf5_version "1.12.0")
set(_hdf5_pkg_url "https://support.hdfgroup.org/ftp/HDF5/releases/hdf5-1.12/hdf5-1.12.0/src/hdf5-1.12.0.tar.gz")
file(MAKE_DIRECTORY ${PROJECT_BINARY_DIR}/hdf5_local/src)
ExternalProject_Add(hdf5_local
    DOWNLOAD_COMMAND
        curl ${_hdf5_pkg_url} | tar xvz -C ${PROJECT_BINARY_DIR}/hdf5_local/src
    CMAKE_ARGS
        ${ExternalProjectCMakeArgs}
    SOURCE_DIR
        ${PROJECT_BINARY_DIR}/hdf5_local/src/hdf5-${_hdf5_version}
    BINARY_DIR
        ${PROJECT_BINARY_DIR}/hdf5_local/build
    STAMP_DIR
        ${PROJECT_BINARY_DIR}/hdf5_local/stamp
    TMP_DIR
        ${PROJECT_BINARY_DIR}/hdf5_local/tmp
    INSTALL_COMMAND
        ""
)

set(_hdf5_include_dirs
    "${PROJECT_BINARY_DIR}/hdf5_local/src/hdf5-${_hdf5_version}/src"
    "${PROJECT_BINARY_DIR}/hdf5_local/src/hdf5-${_hdf5_version}/c++/src"
    "${PROJECT_BINARY_DIR}/hdf5_local/build"
)
set(_hdf5_include_dirs ${_hdf5_include_dirs} CACHE INTERNAL "Local HDF5 include directory paths")

set(_hdf5_libs
    "${PROJECT_BINARY_DIR}/hdf5_local/build/bin/libhdf5_cpp.a"
    "${PROJECT_BINARY_DIR}/hdf5_local/build/bin/libhdf5.a"
    "-ldl"
)
set(_hdf5_libs ${_hdf5_libs} CACHE INTERNAL "Local HDF5 libraries")

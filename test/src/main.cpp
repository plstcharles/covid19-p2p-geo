#include <sys/stat.h>

#include "gtest/gtest.h"
#include "api.hpp"

inline bool file_exists(const std::string& sFilePath) {
    struct stat oBuffer;
    return (stat(sFilePath.c_str(), &oBuffer) == 0);
}

const std::string g_sDataRootPath = DATA_ROOT;
const std::string g_sExpectedFullHDF5Path = g_sDataRootPath + "/full.hdf5";

TEST(covid19_p2p_geo_api, data_exists) {
    std::cout << "g_sExpectedFullHDF5Path = " << g_sExpectedFullHDF5Path << std::endl;
    ASSERT_TRUE(file_exists(g_sExpectedFullHDF5Path));
}

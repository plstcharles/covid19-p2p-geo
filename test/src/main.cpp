#include <sys/stat.h>

#include "gtest/gtest.h"
#include "api.hpp"

inline bool checkPathExists(const std::string& sFilePath) {
    struct stat oBuffer;
    return (stat(sFilePath.c_str(), &oBuffer) == 0);
}

const std::string g_sDataRootPath = DATA_ROOT;

TEST(covid19_p2p_geo_api, data_exists) {
    const std::string g_sExpectedFullHDF5Path = g_sDataRootPath + "/full.hdf5";
    ASSERT_TRUE(checkPathExists(g_sExpectedFullHDF5Path));
    const std::string g_sExpectedCDHDF5Path = g_sDataRootPath + "/cd.hdf5";
    ASSERT_TRUE(checkPathExists(g_sExpectedCDHDF5Path));
    const std::string g_sExpectedDAHDF5Path = g_sDataRootPath + "/divisions";
    ASSERT_TRUE(checkPathExists(g_sExpectedDAHDF5Path));
}

TEST(covid19_p2p_geo_api, random_build_and_queries) {
    testRandomBuildAndQueries(g_sDataRootPath);
}

TEST(covid19_p2p_geo_api, mila_coord_check) {
    // check if we can match the data obtained manually from censusmapper.ca
    const std::string sExpectedMilaUID = "24661626";
    const double dLat = 45.530637, dLon = -73.613989;
    const std::string sUID = fetchRegionUID(dLat, dLon, g_sDataRootPath);
    ASSERT_EQ(sUID, sExpectedMilaUID);
}

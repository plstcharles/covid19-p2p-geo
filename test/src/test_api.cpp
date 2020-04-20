#include <sys/stat.h>

#include "gtest/gtest.h"
#include "private/generic_utils.hpp"
#include "api.hpp"

const std::string g_sDataRootPath = DATA_ROOT;

// check if we can match the data obtained manually from censusmapper.ca
const uint32_t g_nExpectedMilaUID = 24661626;
const double g_dMilaLatitude = 45.530637, g_dMilaLongitude = -73.613989;

TEST(covid19_p2p_geo_api, data_exists) {
    // assume all data is already available and no need to download externally
    const std::string g_sExpectedFullHDF5Path = g_sDataRootPath + "/full.hdf5";
    ASSERT_TRUE(checkPathExists(g_sExpectedFullHDF5Path));
    const std::string g_sExpectedCDHDF5Path = g_sDataRootPath + "/cd.hdf5";
    ASSERT_TRUE(checkPathExists(g_sExpectedCDHDF5Path));
    const std::string g_sExpectedDAHDF5Path = g_sDataRootPath + "/divisions";
    ASSERT_TRUE(checkPathExists(g_sExpectedDAHDF5Path));
}

TEST(covid19_p2p_geo_api, mila_coord_check) {
    const uint32_t nUID = fetchUID(g_dMilaLatitude, g_dMilaLongitude, g_sDataRootPath.c_str());
    ASSERT_EQ(nUID, g_nExpectedMilaUID);
}

TEST(covid19_p2p_geo_api, mila_prepared_coord_check) {
    prepareNear(g_dMilaLatitude, g_dMilaLongitude, g_sDataRootPath.c_str());
    ASSERT_GT(getCachedRegionTreeCount(), 0u);
    const uint32_t nUID = fetchUID(g_dMilaLatitude, g_dMilaLongitude, g_sDataRootPath.c_str());
    ASSERT_EQ(nUID, g_nExpectedMilaUID);
}

TEST(covid19_p2p_geo_api, prep_and_release) {
    releaseUnusedCache(0);
    ASSERT_EQ(getCachedRegionTreeCount(), 0u);
    prepareNear(g_dMilaLatitude, g_dMilaLongitude, g_sDataRootPath.c_str());
    const uint32_t nInitCachedRegionCount = getCachedRegionTreeCount();
    ASSERT_GT(nInitCachedRegionCount, 0u);
    sleep(1);
    ASSERT_EQ(getCachedRegionTreeCount(), nInitCachedRegionCount);
    releaseUnusedCache(1000); // ...should not fail, right?
    ASSERT_EQ(getCachedRegionTreeCount(), nInitCachedRegionCount);
    releaseUnusedCache(0);
    ASSERT_EQ(getCachedRegionTreeCount(), 0u);
}

TEST(covid19_p2p_geo_api, random_build_and_queries) {
    testRandomBuildAndQueries(g_sDataRootPath.c_str());
}

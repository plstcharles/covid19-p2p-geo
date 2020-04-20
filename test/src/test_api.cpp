#include <sys/stat.h>

#ifdef __linux__
#include <iostream>
#include <fstream>
#include <sys/types.h>
#endif //def(__linux__)

#include "gtest/gtest.h"
#include "private/generic_utils.hpp"
#include "api.hpp"

const std::string g_sDataRootPath = DATA_ROOT;

const std::string g_sEmptyDataRootPath = "/tmp/covid19-p2p-geo-test";

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

#ifdef __linux__

inline void copyFile(const std::string& sInput, const std::string& sOutput) {
    std::ifstream hInput(sInput, std::ios::binary);
    std::ofstream hOutput(sOutput, std::ios::binary);
    hOutput << hInput.rdbuf();
}

TEST(covid19_p2p_geo_api, mila_full_coord_check_from_nothing) {
    ASSERT_TRUE(checkPathExists("/tmp")); // can't run test if /tmp does not exist...
    if(!checkPathExists(g_sEmptyDataRootPath)) {
        ASSERT_EQ(mkdir(g_sEmptyDataRootPath.c_str(), 0777), 0);
        ASSERT_EQ(mkdir((g_sEmptyDataRootPath + "/divisions").c_str(), 0777), 0);
        copyFile(g_sDataRootPath + "/cd.hdf5", g_sEmptyDataRootPath + "/cd.hdf5");
    }
    // remove old divisions if the test was already executed
    if(checkPathExists(g_sEmptyDataRootPath + "/divisions/2465.hdf5"))
        remove((g_sEmptyDataRootPath + "/divisions/2465.hdf5").c_str());
    if(checkPathExists(g_sEmptyDataRootPath + "/divisions/2466.hdf5"))
        remove((g_sEmptyDataRootPath + "/divisions/2466.hdf5").c_str());
    const char* acTargetFilenames =
            prepareNear(g_dMilaLatitude, g_dMilaLongitude, g_sEmptyDataRootPath.c_str());
    ASSERT_NE(acTargetFilenames, nullptr);
    ASSERT_NE(acTargetFilenames[0], 0);
    std::string sTargetFilenames = acTargetFilenames;
    std::vector<std::string> vsTargetFilenames;
    size_t nStrPos = 0u;
    while((nStrPos = sTargetFilenames.find(',')) != std::string::npos) {
        vsTargetFilenames.push_back(sTargetFilenames.substr(0, nStrPos));
        sTargetFilenames.erase(0, nStrPos + 1u);
    }
    vsTargetFilenames.push_back(sTargetFilenames);
    for(const auto& sTargetFilename: vsTargetFilenames) {
        const std::string sInput = g_sDataRootPath + "/divisions/" + sTargetFilename + ".hdf5";
        const std::string sOutput = g_sEmptyDataRootPath + "/divisions/" + sTargetFilename + ".hdf5";
        copyFile(sInput, sOutput);
    }
    acTargetFilenames = prepareNear(g_dMilaLatitude, g_dMilaLongitude, g_sEmptyDataRootPath.c_str());
    ASSERT_NE(acTargetFilenames, nullptr);
    ASSERT_EQ(acTargetFilenames[0], 0);
    ASSERT_GT(getCachedRegionTreeCount(), 0u);
    const uint32_t nUID = fetchUID(g_dMilaLatitude, g_dMilaLongitude, g_sDataRootPath.c_str());
    ASSERT_EQ(nUID, g_nExpectedMilaUID);
}

#endif //def(__linux__)


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

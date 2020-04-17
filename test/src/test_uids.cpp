
#include "gtest/gtest.h"
#include "private/uid_utils.hpp"


TEST(covid19_p2p_geo_uids, province_and_territory) {
    const GeoRegionUID anPositiveTestVals[] = {10, 42, 99};
    for(GeoRegionUID nPosVal : anPositiveTestVals) {
        ASSERT_TRUE(uid::isProvinceOrTerritory(nPosVal));
        ASSERT_TRUE(uid::isProvinceOrTerritory(std::to_string(nPosVal)));
    }
    const GeoRegionUID anNegativeTestVals[] = {0, 9, 100, 101};
    for(GeoRegionUID nNegVal : anNegativeTestVals) {
        ASSERT_FALSE(uid::isProvinceOrTerritory(nNegVal));
        ASSERT_FALSE(uid::isProvinceOrTerritory(std::to_string(nNegVal)));
    }
}

TEST(covid19_p2p_geo_uids, census_division) {
    const GeoRegionUID anPositiveTestVals[] = {1000, 1001, 1234, 4242, 9999};
    for(GeoRegionUID nPosVal : anPositiveTestVals) {
        ASSERT_TRUE(uid::isCensusDivision(nPosVal));
        ASSERT_TRUE(uid::isCensusDivision(std::to_string(nPosVal)));
    }
    const GeoRegionUID anNegativeTestVals[] = {0, 9, 42, 99, 10000, 10001, 125133};
    for(GeoRegionUID nNegVal : anNegativeTestVals) {
        ASSERT_FALSE(uid::isCensusDivision(nNegVal));
        ASSERT_FALSE(uid::isCensusDivision(std::to_string(nNegVal)));
    }
}

TEST(covid19_p2p_geo_uids, dissemination_areas) {
    const GeoRegionUID anPositiveTestVals[] = {10000000, 42536475, 99999999};
    for(GeoRegionUID nPosVal : anPositiveTestVals) {
        ASSERT_TRUE(uid::isDisseminationArea(nPosVal));
        ASSERT_TRUE(uid::isDisseminationArea(std::to_string(nPosVal)));
    }
    const GeoRegionUID anNegativeTestVals[] = {
            0, 10, 1000, 9999999, std::numeric_limits<GeoRegionUID>::max()};
    for(GeoRegionUID nNegVal : anNegativeTestVals) {
        ASSERT_FALSE(uid::isDisseminationArea(nNegVal));
        ASSERT_FALSE(uid::isDisseminationArea(std::to_string(nNegVal)));
    }
}

TEST(covid19_p2p_geo_uids, parents_uint32) {
    ASSERT_EQ(uid::getParentUID(0u), 0u);
    ASSERT_EQ(uid::getParentUID(1u), 0u);
    ASSERT_EQ(uid::getParentUID(10u), 0u);
    ASSERT_EQ(uid::getParentUID(1000u), 10u);
    ASSERT_EQ(uid::getParentUID(4253u), 42u);
    ASSERT_EQ(uid::getParentUID(42536475u), 4253u);
}

TEST(covid19_p2p_geo_uids, parents_str) {
    ASSERT_EQ(uid::getParentUID("0"), "0");
    ASSERT_EQ(uid::getParentUID("1"), "0");
    ASSERT_EQ(uid::getParentUID("10"), "0");
    ASSERT_EQ(uid::getParentUID("1000"), "10");
    ASSERT_EQ(uid::getParentUID("4253"), "42");
    ASSERT_EQ(uid::getParentUID("42536475"), "4253");
}

TEST(covid19_p2p_geo_uids, parents_iter) {
    const size_t nUIDs = 100u;
    std::vector<GeoRegionUID> vMixedCDUIDs, vCDUIDs, vDAUIDs;
    ASSERT_EQ(uid::getParentUID(vCDUIDs.begin(), vCDUIDs.end()), 0u);
    for(size_t nUIDIdx = 0u; nUIDIdx < nUIDs; ++nUIDIdx) {
        if(std::rand()%2)
            vMixedCDUIDs.push_back(4200u + (uint32_t)std::rand() % 100u);
        else
            vMixedCDUIDs.push_back(8400u + (uint32_t)std::rand() % 100u);
        vCDUIDs.push_back(4200u + (uint32_t)std::rand() % 100u);
        vDAUIDs.push_back(42420000u + (uint32_t)std::rand() % 10000u);
    }
    ASSERT_EQ(uid::getParentUID(vMixedCDUIDs.begin(), vMixedCDUIDs.end()), 0u);
    ASSERT_EQ(uid::getParentUID(vCDUIDs.begin(), vCDUIDs.end()), 42u);
    ASSERT_EQ(uid::getParentUID(vDAUIDs.begin(), vDAUIDs.end()), 4242u);
}

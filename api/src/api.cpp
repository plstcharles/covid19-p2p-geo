
#include "api.hpp"

#include "regions.hpp"

void prepareRegionsNear(double dLatitude, double dLongitude, const std::string& sDataRootPath) {
    GeoRegionTreePtr pHighLevelTree = prepareHighLevelRegionTree(sDataRootPath);
    assert(pHighLevelTree);
    const std::vector<GeoRegionPtr>& vHits = fetchRegionHits(dLatitude, dLongitude, pHighLevelTree);
    for(GeoRegionPtr pHitRegion : vHits) {
        assert(pHitRegion && pHitRegion->nParentUID == GLOBAL_REGION_UID);
        const std::string& sTargetSubRegion = std::to_string(pHitRegion->nUID);
        prepareRegionTree(sTargetSubRegion, sDataRootPath);
    }
}

std::string fetchRegionsUID(double dLatitude, double dLongitude, const std::string& sDataRootPath) {
    GeoRegionTreePtr pHighLevelTree = prepareHighLevelRegionTree(sDataRootPath);
    assert(pHighLevelTree);
    const std::vector<GeoRegionPtr>& vHits = fetchRegionHits(dLatitude, dLongitude, pHighLevelTree);
    for(GeoRegionPtr pHitRegion : vHits) {
        assert(pHitRegion && pHitRegion->nParentUID == GLOBAL_REGION_UID);
        const std::string& sTargetSubRegion = std::to_string(pHitRegion->nUID);
        GeoRegionTreePtr pRegionTree = prepareRegionTree(sTargetSubRegion, sDataRootPath);
        assert(pRegionTree && pRegionTree->nParentUID != GLOBAL_REGION_UID);
        GeoRegionPtr pResult = fetchRegion(dLatitude, dLongitude, pRegionTree);
        if(pResult)
            return std::to_string(pResult->nUID);
    }
    return GLOBAL_REGION_STR; // no clean intersection found, return global ID (0)
}

void testRandomBuildAndQueries(const std::string& sDataRootPath) {
    GeoRegionTreePtr pHighLevelTree = prepareHighLevelRegionTree(sDataRootPath);
    const std::string sHDF5DivisonsFolderPath = sDataRootPath + "/divisions/";
    const size_t nRandomBuildAndQueryCount = 1000u;
    std::vector<double> vdBuildTotalTime_msec, vdQueryTotalTime_msec;
    for(size_t nTrialIdx = 0u; nTrialIdx < nRandomBuildAndQueryCount; ++nTrialIdx) {
        // pick a random census division for the test
        const size_t nTargetRegionIdx = (size_t)std::rand() % pHighLevelTree->mRegions.size();
        auto pRegionIter = pHighLevelTree->mRegions.begin();
        std::advance(pRegionIter, nTargetRegionIdx);
        const GeoRegionPtr pTargetRegion = pRegionIter->second;
        const GeoRegionUID nTargetCDUID = pTargetRegion->nUID;
        const std::string sTargetCDUID = std::to_string(nTargetCDUID);
        assert(sTargetCDUID.size() == 4u); // all CDUIDs should be 4-digit
        assert(!pTargetRegion->pGeometry);
        GeoRegionTreePtr pSubTree = GeoRegionTreeCacher::getGeoRegionTree(sTargetCDUID);
        if(!pSubTree) {
            const std::string sHDF5FilePath = sHDF5DivisonsFolderPath + sTargetCDUID + ".hdf5";
            auto tPreBuildTimestamp = std::chrono::high_resolution_clock::now();
            const GeoRegionMap& mSubRegions = createDisseminationAreaMap(sHDF5FilePath);
            if(std::rand() % 2) {
                // 50% of the time, ask with a non-default name
                const std::string sSessionName = sTargetCDUID + "_test";
                auto oInserted = GeoRegionTreeCacher::createGeoRegionTree(mSubRegions, sSessionName);
                pSubTree = oInserted.second;
                assert(oInserted.first == sSessionName); // we requested that name, we should have it
            }
            else {
                // 50% of the time, ask with the default name (should resolve to sTargetCDUID)
                auto oInserted = GeoRegionTreeCacher::createGeoRegionTree(mSubRegions);
                pSubTree = oInserted.second;
                assert(oInserted.first == sTargetCDUID); // parent of all DAs in a CD = the CD itself
            }
            auto tPostBuildTimestamp = std::chrono::high_resolution_clock::now();
            const std::chrono::duration<double, std::milli> tBuildTime =
                    tPostBuildTimestamp - tPreBuildTimestamp;
            vdBuildTotalTime_msec.push_back(tBuildTime.count());
        }
        assert(pSubTree);
        // within the census division region tree, pick a random dissemination area as a query
        const size_t nTargetSubRegionIdx = (size_t)std::rand() % pSubTree->mRegions.size();
        auto pSubRegionIter = pSubTree->mRegions.begin();
        std::advance(pSubRegionIter, nTargetSubRegionIdx);
        const GeoRegionPtr pTargetSubRegion = pSubRegionIter->second;
        const GeoRegionUID nTargetDAUID = pTargetSubRegion->nUID;
        const std::string sTargetDAUID = std::to_string(nTargetDAUID);
        assert(sTargetDAUID.size() == 8u); // all DAUIDs should be 8-digit
        assert(pTargetSubRegion->pGeometry);
        const auto tPreQueryTimestamp = std::chrono::high_resolution_clock::now();
        Geometry pQueryPt = pTargetSubRegion->pGeometry->getCentroid();
        const GeoRegionPtr pResult = pSubTree->findGeoRegion(pQueryPt.get());
        const auto tPostQueryTimestamp = std::chrono::high_resolution_clock::now();
        const std::chrono::duration<double, std::milli> tQueryTime =
                tPostQueryTimestamp - tPreQueryTimestamp;
        vdQueryTotalTime_msec.push_back(tQueryTime.count());
        if(pTargetSubRegion->pGeometry->contains(pQueryPt.get())) {
            assert(pResult->nUID == nTargetDAUID);
        }
        else {
            assert(!pResult || pResult->nUID != nTargetDAUID);
        }
    }
    std::cout << std::endl;
    std::cout << "census division tree build time (msec):\n" << writeStatsToString(vdBuildTotalTime_msec) << std::endl;
    std::cout << "census division tree query time (msec):\n" << writeStatsToString(vdQueryTotalTime_msec) << std::endl;
}

/*
# other impl ideas:  @@@@@
#  - query w/ buffer around GPS coord, return average of whatev-values stored in geometries?
#      (compute intersection w/ bbox instead of just a point)
#      buffer size could be determined based on GPS coord accuracy
#  - skeletize all 400-700 pop geometries down to 100-200? (hyperparam)
#     (use voronoi diagram builder inside geos, increase size/overlap with buffer)
#  - provide 2d map of min/max distances between all cd-level region pairs?
*/
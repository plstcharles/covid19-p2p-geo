
#include "api.hpp"

#include "hdf5_utils.hpp"

/// Creates a high-level region map (based on census division envelopes) for spatial querying
GeoRegionMap createHighLevelRegionMap(const std::string& sHDF5FilePath) {
    const H5::H5File oH5Archive(sHDF5FilePath, H5F_ACC_RDONLY);
    const hsize_t nCDCount = readHDF5Int64Attrib<hsize_t>(oH5Archive, "cd_count");
    assert(nCDCount > 0);
    const GeoRegionStatsArray vGeoRegionStats = readHDF5GeoRegionStats(oH5Archive, "cd_stats");
    assert(nCDCount == vGeoRegionStats.size());
    const GeoRegionUIDArray vnUIDs = readHDF5GeoRegionUIDs(oH5Archive, "cduid");
    assert(nCDCount == vnUIDs.size());
    const std::vector<GeomEnvelope> vEnvelopes = readHDF5GeomEnvelopes(oH5Archive, "cd_envelope");
    assert(nCDCount == vEnvelopes.size());
    // note: high-level geometries don't come with actual geometries, just envelopes
    GeoRegionMap mRegions;
    for(size_t nCDIdx = 0u; nCDIdx < nCDCount; ++nCDIdx) {
        auto pRegion = std::make_shared<GeoRegion>(vGeoRegionStats[nCDIdx], vnUIDs[nCDIdx], vEnvelopes[nCDIdx]);
        mRegions.insert(std::make_pair(vnUIDs[nCDIdx], pRegion));
    }
    return mRegions;
}

/// Creates a region map (based on dissemination area geometries) for spatial querying
GeoRegionMap createDisseminationAreaMap(const std::string& sHDF5FilePath) {
    const H5::H5File oH5Archive(sHDF5FilePath, H5F_ACC_RDONLY);
    const hsize_t nDACount = readHDF5Int64Attrib<hsize_t>(oH5Archive, "da_count");
    assert(nDACount > 0);
    const GeoRegionStatsArray vGeoRegionStats = readHDF5GeoRegionStats(oH5Archive, "da_stats");
    assert(nDACount == vGeoRegionStats.size());
    const GeoRegionUIDArray vnUIDs = readHDF5GeoRegionUIDs(oH5Archive, "dauid");
    assert(nDACount == vnUIDs.size());
    std::vector<Geometry> vDAGeometries = readHDF5DAWKBGeometries(oH5Archive);
    assert(nDACount == vDAGeometries.size());
    GeoRegionMap mRegions;
    for(size_t nDAIdx = 0u; nDAIdx < nDACount; ++nDAIdx) {
        const GeomEnvelope* pRegionEnvelope = vDAGeometries[nDAIdx]->getEnvelopeInternal();
        auto pRegion = std::make_shared<GeoRegion>(
                vGeoRegionStats[nDAIdx], vnUIDs[nDAIdx], *pRegionEnvelope, std::move(vDAGeometries[nDAIdx]));
        mRegions.insert(std::make_pair(vnUIDs[nDAIdx], pRegion));
    }
    return mRegions;
}

void testRandomBuildAndQueries(const std::string& sDataRootPath) {
    const std::string sHighLevelHDF5FilePath = sDataRootPath + "/cd.hdf5";
    const std::string sHDF5DivisonsFolderPath = sDataRootPath + "/divisions/";
    const GeoRegionMap& mHighLevelRegions = createHighLevelRegionMap(sHighLevelHDF5FilePath);
    std::pair<SessionNameType, GeoRegionTreePtr> oHighLevelTree =
            GeoRegionTreeCacher::createGeoRegionTree(mHighLevelRegions);
    assert(oHighLevelTree.first == "0"); // high-level stuff parent should be global region (0)
    const size_t nRandomBuildAndQueryCount = 1000u;
    std::vector<double> vdBuildTotalTime_msec, vdQueryTotalTime_msec;
    for(size_t nTrialIdx = 0u; nTrialIdx < nRandomBuildAndQueryCount; ++nTrialIdx) {
        // pick a random census division for the test
        const size_t nTargetRegionIdx = (size_t)std::rand() % mHighLevelRegions.size();
        auto pRegionIter = mHighLevelRegions.begin();
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
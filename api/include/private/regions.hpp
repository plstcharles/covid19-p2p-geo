#pragma once

#include "private/hdf5_utils.hpp"

/// creates a high-level region map (based on census division envelopes) for spatial querying
inline GeoRegionMap createHighLevelRegionMap(const std::string& sHDF5FilePath) {
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

/// creates a region map (based on dissemination area geometries) for spatial querying
inline GeoRegionMap createDisseminationAreaMap(const std::string& sHDF5FilePath) {
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

/// prepares the high-level region tree (based on census division envelopes)
inline GeoRegionTreePtr prepareHighLevelRegionTree(const std::string& sDataRootPath) {
    GeoRegionTreePtr pHighLevelTree = GeoRegionTreeCacher::getGeoRegionTree(GLOBAL_REGION_STR);
    if(!pHighLevelTree) {
        const std::string sHighLevelHDF5FilePath = sDataRootPath + "/cd.hdf5";
        if(!checkPathExists(sHighLevelHDF5FilePath))
            throw std::runtime_error(std::string("invalid archive path: ") + sHighLevelHDF5FilePath);
        const GeoRegionMap& mHighLevelRegions = createHighLevelRegionMap(sHighLevelHDF5FilePath);
        std::pair<SessionNameType, GeoRegionTreePtr> oHighLevelTree =
                GeoRegionTreeCacher::createGeoRegionTree(mHighLevelRegions);
        assert(oHighLevelTree.first == GLOBAL_REGION_STR); // top parent should be global region
        pHighLevelTree = oHighLevelTree.second;
    }
    return pHighLevelTree;
}

/// prepares a region tree (based on dissemination area geometries)
inline GeoRegionTreePtr prepareRegionTree(const std::string& sUID, const std::string& sDataRootPath) {
    assert(sUID.size() == 4u); // should always be querying census divisions
    GeoRegionTreePtr pTree = GeoRegionTreeCacher::getGeoRegionTree(sUID);
    if(!pTree) {
        const std::string sHDF5FilePath = sDataRootPath + "/divisions/" + sUID +".hdf5";
        if(!checkPathExists(sHDF5FilePath))
            throw std::runtime_error(std::string("invalid archive path: ") + sHDF5FilePath);
        const GeoRegionMap& mRegions = createDisseminationAreaMap(sHDF5FilePath);
        std::pair<SessionNameType, GeoRegionTreePtr> oTree =
                GeoRegionTreeCacher::createGeoRegionTree(mRegions);
        assert(oTree.first == sUID); // should provide queried key, always
        pTree = oTree.second;
    }
    return pTree;
}

/// returns the list of tree regions that intersect the given lat/lon coordinates
inline std::vector<GeoRegionPtr> fetchRegionHits(
        double dLatitude,
        double dLongitude,
        GeoRegionTreePtr pTree) {
    geos::geom::GeometryFactory::Ptr pGeomFact = geos::geom::GeometryFactory::create();
    const geos::geom::Coordinate oTargetCoords(dLatitude, dLongitude);
    std::unique_ptr<geos::geom::Point> pPoint(pGeomFact->createPoint(oTargetCoords));
    const geos::geom::Envelope* pGeomEnv = pPoint->getEnvelopeInternal();
    return pTree->findGeoRegionHits(pGeomEnv);
}

/// returns the first tree region that intersects the given lat/lon coordinates
inline GeoRegionPtr fetchRegion(
        double dLatitude,
        double dLongitude,
        GeoRegionTreePtr pTree) {
    geos::geom::GeometryFactory::Ptr pGeomFact = geos::geom::GeometryFactory::create();
    const geos::geom::Coordinate oTargetCoords(dLatitude, dLongitude);
    std::unique_ptr<geos::geom::Point> pPoint(pGeomFact->createPoint(oTargetCoords));
    return pTree->findGeoRegion(pPoint); // will return nullptr if no hit is found
}

/// returns the duration (in seconds) since the last access of a specific region's cache
inline double getCacheLastAccess(const SessionNameType& sUID) {
    if(uid::isDisseminationArea(sUID))
        return getCacheLastAccess(uid::getParentUID(sUID));
    const GeoRegionTreePtr pTree = GeoRegionTreeCacher::getGeoRegionTree(sUID);
    if(!pTree)
        return std::numeric_limits<double>::quiet_NaN();
    const std::chrono::duration<double> tTimeDiff =
            std::chrono::high_resolution_clock::now() - pTree->getLastAccessedTimestamp();
    return tTimeDiff.count();
}

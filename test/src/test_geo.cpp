
#include "gtest/gtest.h"
#include "private/geo_utils.hpp"

inline GeoRegionArray getDummyRegions() {
    std::vector<GeoRegionPtr> vRegions;
    vRegions.push_back(
            std::make_shared<GeoRegion>(
                    GeoRegionStats(), 1u, GeomEnvelope(0, 1, 0, -1)));
    vRegions.push_back(
            std::make_shared<GeoRegion>(
                    GeoRegionStats(), 2u, GeomEnvelope(1, 2, 0, -1)));
    vRegions.push_back(
            std::make_shared<GeoRegion>(
                    GeoRegionStats(), 3u, GeomEnvelope(0, 1, -1, -2)));
    vRegions.push_back(
            std::make_shared<GeoRegion>(
                    GeoRegionStats(), 4u, GeomEnvelope(1, 2, -1, -2)));
    return vRegions;
}

TEST(covid19_p2p_geo_utils, region_map_getter) {
    const auto vRegions = getDummyRegions();
    auto mRegions = getGeoRegionMapFromArray(vRegions);
    ASSERT_EQ(mRegions.find(0u), mRegions.end());
    ASSERT_NE(mRegions.find(1u), mRegions.end());
    ASSERT_NE(mRegions.find(2u), mRegions.end());
    ASSERT_NE(mRegions.find(3u), mRegions.end());
    ASSERT_NE(mRegions.find(4u), mRegions.end());
    ASSERT_EQ(mRegions.find(1u)->second.get(), vRegions[0u].get());
}

inline GeomArray getDummyGeoms(const GeoRegionArray& vRegions) {
    std::vector<Geometry> vGeoms;
    geos::geom::GeometryFactory::Ptr pGeomFact = geos::geom::GeometryFactory::create();
    for(size_t nRegionIdx = 0u; nRegionIdx < vRegions.size(); ++nRegionIdx)
        vGeoms.push_back(pGeomFact->toGeometry(&vRegions[nRegionIdx]->oEnvelope));
    return vGeoms;
}

TEST(covid19_p2p_geo_utils, env_getter) {
    const auto vRegions = getDummyRegions();
    auto vGeoms = getDummyGeoms(vRegions);
    auto pGeomEnv = getGeomArrayEnvelope(vGeoms);
    auto pEnv = pGeomEnv->getEnvelopeInternal();
    ASSERT_DOUBLE_EQ(pEnv->getMinX(), 0.0);
    ASSERT_DOUBLE_EQ(pEnv->getMaxX(), 2.0);
    ASSERT_DOUBLE_EQ(pEnv->getMinY(), -2.0);
    ASSERT_DOUBLE_EQ(pEnv->getMaxY(), 0.0);
    auto pGeomEnv2 = getGeomArrayEnvelope(vRegions);
    auto pEnv2 = pGeomEnv2->getEnvelopeInternal();
    ASSERT_DOUBLE_EQ(pEnv->getMinX(), pEnv2->getMinX());
    ASSERT_DOUBLE_EQ(pEnv->getMaxX(), pEnv2->getMaxX());
    ASSERT_DOUBLE_EQ(pEnv->getMinY(), pEnv2->getMinY());
    ASSERT_DOUBLE_EQ(pEnv->getMaxY(), pEnv2->getMaxY());
}

TEST(covid19_p2p_geo_utils, geo_region_utils) {
    const auto vRegions = getDummyRegions();
    GeoRegionTree oTree(vRegions);
    geos::geom::GeometryFactory::Ptr pGeomFact = geos::geom::GeometryFactory::create();
    auto tCurrTime = std::chrono::high_resolution_clock::now();
    auto pGeomEnv = pGeomFact->toGeometry(&vRegions[0]->oEnvelope);
    auto pCentroid = pGeomEnv->getCentroid();
    auto pHitRegion = oTree.findGeoRegion(pCentroid.get());
    ASSERT_TRUE(pHitRegion != nullptr);
    ASSERT_EQ(pHitRegion->oEnvelope, vRegions[0]->oEnvelope);
    ASSERT_LT(tCurrTime, oTree.getLastAccessedTimestamp());
    auto vRegionHits =  oTree.findGeoRegionHits(pCentroid->getEnvelopeInternal());
    ASSERT_EQ(vRegionHits.size(), 1u);
    ASSERT_EQ(vRegionHits[0u].get(), vRegionHits[0u].get());
}

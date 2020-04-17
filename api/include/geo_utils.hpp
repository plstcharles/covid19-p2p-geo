
/// contains GEOS and geo-related classes and utilities

#pragma once

#include <geos/geom/PrecisionModel.h>
#include <geos/geom/GeometryFactory.h>
#include <geos/geom/Geometry.h>
#include <geos/geom/Point.h>
#include <geos/geom/GeometryCollection.h>
#include <geos/geom/Coordinate.h>
#include <geos/geom/CoordinateSequence.h>
#include <geos/geom/IntersectionMatrix.h>
#include <geos/io/WKBReader.h>
#include <geos/index/strtree/STRtree.h>
#include <geos/operation/linemerge/LineMerger.h>
#include <geos/constants.h>

#include "generic_utils.hpp"

struct GeoRegionStats;
struct GeoRegionTree;
struct GeoRegion;

using Geometry = geos::geom::Geometry::Ptr;
using GeomEnvelope = geos::geom::Envelope;
using GeomArray = std::vector<Geometry>;

using GeoRegionUID = uint32_t;
using GeoRegionPtr = std::shared_ptr<GeoRegion>;
using GeoRegionUIDArray = std::vector<GeoRegionUID>;
using GeoRegionArray = std::vector<GeoRegionPtr>;
using GeoRegionStatsArray = std::vector<GeoRegionStats>;
using GeoRegionMap = std::map<GeoRegionUID, GeoRegionPtr>;
using GeoRegionTreePtr = std::shared_ptr<GeoRegionTree>;

using SessionNameType = std::string;

#define GLOBAL_REGION_UID GeoRegionUID(0u)
#define GLOBAL_REGION_STR std::to_string(GLOBAL_REGION_UID)
#define AUTO_SESSION_NAME std::string()

/// returns the unique identifier (UID) of the parent region specified via UID
GeoRegionUID getParentUID(GeoRegionUID nUID) {
    // we use GEO CODES for provincies/territories, census divisions, and dissemination areas
    // (this means they can only be 2-digit, 4-digit, or 8-digit)
    if(nUID < 100u) // we are already at the province/territory level; no parent
        return GLOBAL_REGION_UID;
    else if(nUID < 10000000u) { // we are at the census division level; get first two digits
        assert(nUID >= 1000u && nUID < 10000u); // we should be in the 4-digit range
        return nUID / 100u;
    }
    else { // we are at the dissemination area level; get first four digits
        assert(nUID > 10000000u && nUID < 100000000u);
        return nUID / 10000u;
    }
}

/// returns the unique identifier (UID) of the parent region specified via a UID/value pair
template<typename TKey, typename TVal>
GeoRegionUID getParentUID(const std::pair<TKey, TVal>& oPair) {
    return getParentUID(oPair.first); // assume first key is UID-related
}

/// returns the unique identifier (UID) of the parent region specified via GeoRegion
GeoRegionUID getParentUID(const GeoRegion& oRegion);

/// returns the unique identifier (UID) of the parent region specified via GeoRegionPtr
GeoRegionUID getParentUID(const GeoRegionPtr& pRegion);

/// returns the unique identifier (UID) of the parent region covering an input region range
template<typename TIter>
GeoRegionUID getParentUID(TIter iBegin, TIter iEnd) {
    if(iBegin == iEnd)
        return GLOBAL_REGION_UID;
    std::set<GeoRegionUID> mParentUIDs;
    for(; iBegin != iEnd; ++iBegin)
        mParentUIDs.insert(getParentUID(*iBegin));
    if(mParentUIDs.size() == 1u)
        return *mParentUIDs.begin();
    return getParentUID(mParentUIDs.begin(), mParentUIDs.end());
}

/// geographic region statistics data block; may be enlarged/shrunk based on app needs
struct GeoRegionStats {
    // this is a standalone class not merged with 'GeoRegion' to simplify parsing

    /// constructs a NaN-filled region statistics data block
    GeoRegionStats():
            fPopulation(std::numeric_limits<float>::quiet_NaN()),
            fDwellings(std::numeric_limits<float>::quiet_NaN()),
            fArea(std::numeric_limits<float>::quiet_NaN()) {}

    /// constructs a region statistics data block filled with provided values
    GeoRegionStats(float fPopulation_, float fDwellings_, float fArea_):
            fPopulation(fPopulation_),
            fDwellings(fDwellings_),
            fArea(fArea_) {}

    /// constructs a region statistics data block with an array of values
    explicit GeoRegionStats(const float* afStatsArray):
            fPopulation(afStatsArray[0]),
            fDwellings(afStatsArray[1]),
            fArea(afStatsArray[2]) {}

    /// for now, we only expect HDF5 files to contain three statistics
    static constexpr size_t s_nExpectedStatCount = 3u;
    /// for now, we only expect HDF5 files to contain three statistics
    static std::array<std::string, s_nExpectedStatCount> getExpectedStatNames() {
        return std::array<std::string, s_nExpectedStatCount>{"pop", "dwellings", "area"};
    }

    /// the population contained in the geographic region
    float fPopulation;
    /// the number of private residential dwellings in the geographic region
    float fDwellings;
    /// the area (in square kilometers) of the geographic region
    float fArea;
};

/// geographic region data block instantiated as part of search trees user data
struct GeoRegion {

    /// constructs a geographic region; the class members should never change once built
    GeoRegion(
            const GeoRegionStats& oStats_, GeoRegionUID nUID_, const GeomEnvelope& oEnvelope_,
            Geometry&& pGeometry_ = nullptr):
            nUID(nUID_), nParentUID(getParentUID(nUID_)), oEnvelope(oEnvelope_), oStats(oStats_),
            pGeometry(std::move(pGeometry_)) {}

    /// the unique identifier (UID) of this geographic region
    const GeoRegionUID nUID;
    /// the unique identifier (UID) of this geographic region's parent
    const GeoRegionUID nParentUID;
    /// the envelope (or bounding box) of this geographic region
    const GeomEnvelope oEnvelope;
    /// the array of statistics associated with this geographic region
    const GeoRegionStats oStats;
    /// the detailed geometry associated with this geographic region (may be null if too high level)
    Geometry pGeometry;
};

/// returns the unique identifier (UID) of the parent region specified via GeoRegion
GeoRegionUID getParentUID(const GeoRegion& oRegion) {
    return getParentUID(oRegion.nUID);
}

/// returns the unique identifier (UID) of the parent region specified via GeoRegionPtr
GeoRegionUID getParentUID(const GeoRegionPtr& pRegion) {
    if(!pRegion)
        return GLOBAL_REGION_UID; // if the region is empty/unavailable, return global ID
    return getParentUID(pRegion->nUID);
}

/// returns a geographic region map (UID,REGION) given an array of regions
GeoRegionMap getGeoRegionMapFromArray(const GeoRegionArray& vRegions) {
    GeoRegionMap mRegions;
    for(const auto& pRegion : vRegions)
        mRegions[pRegion->nUID] = pRegion;
    return mRegions;
}

/// returns the bounding box (polygon) that envelopes all geometries in the given array
Geometry getGeomArrayEnvelope(const GeomArray& vGeoms) {
    geos::geom::GeometryFactory::Ptr pGeomFact = geos::geom::GeometryFactory::create();
    GeomArray vGeomClones;
    vGeomClones.reserve(vGeoms.size());
    for(const auto& oGeom : vGeoms)
        vGeomClones.push_back(oGeom->clone());
    geos::geom::GeometryCollection::Ptr pGeomCollection =
            pGeomFact->createGeometryCollection(std::move(vGeomClones));
    return pGeomCollection->getEnvelope();
}

/// returns the bounding box (polygon) that envelopes all geometries/envelopes in the given array
Geometry getGeomArrayEnvelope(const GeoRegionArray& vRegions) {
    // note: if a region does not possess an actual geometry, we will use its envelope instead
    geos::geom::GeometryFactory::Ptr pGeomFact = geos::geom::GeometryFactory::create();
    GeomArray vGeomClones;
    vGeomClones.reserve(vRegions.size());
    for(const auto& pRegion : vRegions) {
        if(pRegion->pGeometry)
            vGeomClones.push_back(pRegion->pGeometry->clone());
        else
            vGeomClones.push_back(pGeomFact->toGeometry(&pRegion->oEnvelope));
    }
    geos::geom::GeometryCollection::Ptr pGeomCollection =
            pGeomFact->createGeometryCollection(std::move(vGeomClones));
    return pGeomCollection->getEnvelope();
}

/// returns the bounding box (polygon) that envelopes all geometries in the given map
Geometry getGeomArrayEnvelope(const GeoRegionMap& mRegions) {
    return getGeomArrayEnvelope(getValArrayFromMap(mRegions));
}

/// spatial tree used to accelerate the lookup of geographic regions
struct GeoRegionTree {

    /// constructs a geographic region tree for a given array of subregions
    explicit GeoRegionTree(
            const GeoRegionArray& vRegions,
            size_t nDefaultNodeCapacity = 10u):
            tBuildTimestamp(std::chrono::high_resolution_clock::now()),
            mRegions(getGeoRegionMapFromArray(vRegions)),
            nParentUID(getParentUID(vRegions.begin(), vRegions.end())),
            pEnvelope(getGeomArrayEnvelope(vRegions)),
            oTree(nDefaultNodeCapacity) {
        buildTree();
    }

    /// constructs a geographic region tree for a given map of subregions
    explicit GeoRegionTree(
            const GeoRegionMap& mRegions_,
            size_t nDefaultNodeCapacity = 10u):
            tBuildTimestamp(std::chrono::high_resolution_clock::now()),
            mRegions(mRegions_),
            nParentUID(getParentUID(mRegions.begin(), mRegions.end())),
            pEnvelope(getGeomArrayEnvelope(mRegions)), oTree(nDefaultNodeCapacity) {
        buildTree();
    }

    /// returns the first region of the tree that contains the provided geometry
    template<typename TGeometry>
    GeoRegionPtr findGeoRegion(const std::unique_ptr<TGeometry>& pGeometry) {
        return findGeoRegion(pGeometry.get());
    }

    /// returns the first region of the tree that contains the provided geometry
    template<typename TGeometry>
    GeoRegionPtr findGeoRegion(const std::shared_ptr<TGeometry>& pGeometry) {
        return findGeoRegion(pGeometry.get());
    }

    /// returns the first region of the tree that contains the provided geometry
    template<typename TGeometry>
    GeoRegionPtr findGeoRegion(const TGeometry* pGeometry) {
        const geos::geom::Envelope* pGeomEnv = pGeometry->getEnvelopeInternal();
        const std::vector<GeoRegionPtr>& vQueryHits = findGeoRegionHits(pGeomEnv);
        for(GeoRegionPtr pHitRegion : vQueryHits) {
            if(pHitRegion->pGeometry) {
                if(pHitRegion->pGeometry->contains(pGeometry))
                    return pHitRegion;
            }
            else {
                // will always fall back to envelope checks if full geometry is not available
                if(pHitRegion->oEnvelope.contains(pGeomEnv))
                    return pHitRegion;
            }
        }
        // if we did not hit a region, the point is assumed to be "out of bounds"
        return nullptr;
    }

    /// returns the list of tree regions that contain the query envelope
    std::vector<GeoRegionPtr> findGeoRegionHits(const geos::geom::Envelope* pQueryEnv) {
        // abstract tree query func is not constant, better use a mutex for thread safety...
        std::vector<GeoRegionPtr> vHitRegions;
        const std::lock_guard<std::mutex> oLock(m_oTreeMutex);
        std::vector<void*> vQueryHits;
        oTree.query(pQueryEnv, vQueryHits);
        vHitRegions.reserve(vQueryHits.size());
        for(const auto& pHitData : vQueryHits)
            vHitRegions.push_back(*static_cast<GeoRegionPtr*>(pHitData));
        return vHitRegions;
    }

    /// timestamp at which the tree was built (for memory cleanup purposes)
    const std::chrono::high_resolution_clock::time_point tBuildTimestamp;
    /// map of the geographic regions covered by the tree
    const GeoRegionMap mRegions;
    /// unique idenfitier (UID) of the parent of all tree regions
    const GeoRegionUID nParentUID;
    /// envelope (bounding box) geometry covered by this tree
    const Geometry pEnvelope;

private:

    void buildTree() {
        // calling this function mutliple times is undefined behavior
        assert(!mRegions.empty());
        for(auto& pRegionPair : mRegions)
            oTree.insert(&pRegionPair.second->oEnvelope, (void*)(&pRegionPair.second));
        oTree.build();
    }

    geos::index::strtree::STRtree oTree;
    std::mutex m_oTreeMutex;
};

/// tree cache manager used to create/release trees based on user query patterns
struct GeoRegionTreeCacher {

    /// creates a new tree (or overwrites an existing one) for a set of geographic regions
    static std::pair<SessionNameType, GeoRegionTreePtr> createGeoRegionTree(
            const GeoRegionArray& vRegions,
            SessionNameType oSessionName = AUTO_SESSION_NAME,
            size_t nDefaultNodeCapacity = 10u) {
        // if using the automatic session name deduction, the parent region ID will be deducted
        // warning: name deduction + insertion might overwrite an existing parent tree!
        assert(!vRegions.empty());
        if(oSessionName == AUTO_SESSION_NAME)
            oSessionName = suggestSessionName(vRegions);
        auto pRegionTree = std::make_shared<GeoRegionTree>(
                vRegions, nDefaultNodeCapacity);
        const std::lock_guard<std::mutex> oLock(s_oGeoTreeMapMutex);
        return *s_mGeoTrees.insert(std::make_pair(oSessionName, pRegionTree)).first;
    }

    /// creates a new tree (or overwrites an existing one) for a set of geographic regions
    static std::pair<SessionNameType, GeoRegionTreePtr> createGeoRegionTree(
            const GeoRegionMap& mRegions,
            SessionNameType oSessionName = AUTO_SESSION_NAME,
            size_t nDefaultNodeCapacity = 10u) {
        assert(!mRegions.empty());
        if(oSessionName == AUTO_SESSION_NAME)
            oSessionName = suggestSessionName(mRegions);
        auto pRegionTree = std::make_shared<GeoRegionTree>(
                mRegions, nDefaultNodeCapacity);
        const std::lock_guard<std::mutex> oLock(s_oGeoTreeMapMutex);
        return *s_mGeoTrees.insert(std::make_pair(oSessionName, pRegionTree)).first;
    }

    /// returns a tree based on its session name (will return a nullptr if it does not exist)
    static GeoRegionTreePtr getGeoRegionTree(const SessionNameType& sSessionName) {
        assert(sSessionName != AUTO_SESSION_NAME); // not supported here
        const std::lock_guard<std::mutex> oLock(s_oGeoTreeMapMutex);
        auto iMatchedGeoRegionTree = s_mGeoTrees.find(sSessionName);
        if(iMatchedGeoRegionTree == s_mGeoTrees.end())
            return GeoRegionTreePtr();
        return iMatchedGeoRegionTree->second;
    }

    /// returns a proper session name for a set of regions that corresponds to their common parent UID
    static SessionNameType suggestSessionName(const GeoRegionMap& mRegions) {
        const GeoRegionUIDArray vGeoRegionUIds = getKeyArrayFromMap(mRegions);
        const GeoRegionUID nID = getParentUID(mRegions.begin(), mRegions.end());
        return std::to_string(nID);
    }

    /// returns a proper session name for a set of regions that corresponds to their common parent UID
    static SessionNameType suggestSessionName(const GeoRegionArray& vRegions) {
        GeoRegionUIDArray vGeoRegionUIDs;
        vGeoRegionUIDs.reserve(vRegions.size());
        for(const auto& pRegion : vRegions)
            vGeoRegionUIDs.push_back(pRegion->nUID);
        const GeoRegionUID nID = getParentUID(vRegions.begin(), vRegions.end());
        return std::to_string(nID);
    }

private:
    static std::map<SessionNameType, GeoRegionTreePtr> s_mGeoTrees;
    static std::mutex s_oGeoTreeMapMutex;
};

#pragma once

#include "generic_utils.hpp"
#include "geos_utils.hpp"


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

#define AUTO_SESSION_NAME std::string()

#define GLOBAL_REGION_UID GeoRegionUID(0u)


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

template<typename TKey, typename TVal>
GeoRegionUID getParentUID(const std::pair<TKey, TVal>& oPair) {
    return getParentUID(oPair.first); // assume first key is UID-related
}

GeoRegionUID getParentUID(const GeoRegion& oRegion);
GeoRegionUID getParentUID(const GeoRegionPtr& pRegion);

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

/// standalone class not merged with 'GeoRegion' to simplify parsing
struct GeoRegionStats {

    GeoRegionStats():
            fPopulation(std::numeric_limits<float>::quiet_NaN()),
            fDwellings(std::numeric_limits<float>::quiet_NaN()),
            fArea(std::numeric_limits<float>::quiet_NaN()) {}

    GeoRegionStats(float fPopulation_, float fDwellings_, float fArea_):
            fPopulation(fPopulation_),
            fDwellings(fDwellings_),
            fArea(fArea_) {}

    explicit GeoRegionStats(const float* afStatsArray):
            fPopulation(afStatsArray[0]),
            fDwellings(afStatsArray[1]),
            fArea(afStatsArray[2]) {}

    static constexpr size_t s_nExpectedStatCount = 3u;
    static std::array<std::string, s_nExpectedStatCount> getExpectedStatNames() {
        return std::array<std::string, s_nExpectedStatCount>{"pop", "dwellings", "area"};
    }

    float fPopulation;
    float fDwellings;
    float fArea;
};

struct GeoRegion {

    GeoRegion():
            nUID(0), nParentUID(0), oEnvelope(), oStats(), pGeometry(nullptr) {}

    GeoRegion(
            const GeoRegionStats& oStats_, GeoRegionUID nUID_, const GeomEnvelope& oEnvelope_,
            Geometry&& pGeometry_ = nullptr):
            nUID(nUID_), nParentUID(getParentUID(nUID_)), oEnvelope(oEnvelope_), oStats(oStats_),
            pGeometry(std::move(pGeometry_)) {}

    const GeoRegionUID nUID, nParentUID;
    const GeomEnvelope oEnvelope;
    const GeoRegionStats oStats;
    Geometry pGeometry;
};

GeoRegionUID getParentUID(const GeoRegion& oRegion) {
    return getParentUID(oRegion.nUID);
}

GeoRegionUID getParentUID(const GeoRegionPtr& pRegion) {
    if(!pRegion)
        return GLOBAL_REGION_UID;
    return getParentUID(pRegion->nUID);
}

GeoRegionMap getGeoRegionMapFromArray(const GeoRegionArray& vRegions) {
    GeoRegionMap mRegions;
    for(const auto& pRegion : vRegions)
        mRegions[pRegion->nUID] = pRegion;
    return mRegions;
}

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

/// Will use envelope if no actual geometry is available
Geometry getGeomArrayEnvelope(const GeoRegionArray& vRegions) {
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

Geometry getGeomArrayEnvelope(const GeoRegionMap& mRegions) {
    return getGeomArrayEnvelope(getValArrayFromMap(mRegions));
}

struct GeoRegionTree {

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

    explicit GeoRegionTree(
            const GeoRegionMap& mRegions_,
            size_t nDefaultNodeCapacity = 10u):
            tBuildTimestamp(std::chrono::high_resolution_clock::now()),
            mRegions(mRegions_),
            nParentUID(getParentUID(mRegions.begin(), mRegions.end())),
            pEnvelope(getGeomArrayEnvelope(mRegions)), oTree(nDefaultNodeCapacity) {
        buildTree();
    }

    template<typename TGeometry>
    GeoRegionPtr findGeoRegion(const std::unique_ptr<TGeometry>& pGeometry) {
        return findGeoRegion(pGeometry.get());
    }

    template<typename TGeometry>
    GeoRegionPtr findGeoRegion(const std::shared_ptr<TGeometry>& pGeometry) {
        return findGeoRegion(pGeometry.get());
    }

    template<typename TGeometry>
    GeoRegionPtr findGeoRegion(const TGeometry* pGeometry) {
        // note: output pointer will remain valid as long as tree is valid
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

    /// Returns the first-hit tree region that fully contains the query envelope
    std::vector<GeoRegionPtr> findGeoRegionHits(const geos::geom::Envelope* pQueryEnv) {
        // note: output pointers will remain valid as long as tree is valid
        // abstract tree query func is not constant, better use a mutex for thread safety...
        std::vector<void*> vQueryHits;
        {
            const std::lock_guard<std::mutex> oLock(m_oTreeMutex);
            oTree.query(pQueryEnv, vQueryHits);
        }
        std::vector<GeoRegionPtr> vHitRegions;
        vHitRegions.reserve(vQueryHits.size());
        for(const auto& pHitData : vQueryHits)
            vHitRegions.push_back(*static_cast<GeoRegionPtr*>(pHitData));
        return vHitRegions;
    }

    const std::chrono::high_resolution_clock::time_point tBuildTimestamp;
    const GeoRegionMap mRegions;
    const GeoRegionUID nParentUID;
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

struct GeoRegionTreeCacher {

    /// creating a tree with an already-in-use session ID will overwrite that session
    /// if using the automatic session name deduction, the parent region ID will be deducted
    /// automatic session name deduction + insertion might overwrite an existing parent tree
    static std::pair<SessionNameType, GeoRegionTreePtr> createGeoRegionTree(
            const GeoRegionArray& vRegions,
            SessionNameType oSessionName = AUTO_SESSION_NAME,
            size_t nDefaultNodeCapacity = 10u) {
        assert(!vRegions.empty());
        if(oSessionName == AUTO_SESSION_NAME)
            oSessionName = suggestSessionName(vRegions);
        auto pRegionTree = std::make_shared<GeoRegionTree>(
                vRegions, nDefaultNodeCapacity);
        const std::lock_guard<std::mutex> oLock(s_oGeoTreeMapMutex);
        return *s_mGeoTrees.insert(std::make_pair(oSessionName, pRegionTree)).first;
    }

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

    static GeoRegionTreePtr getGeoRegionTree(const SessionNameType& sSessionName) {
        assert(sSessionName != AUTO_SESSION_NAME); // not supported here
        const std::lock_guard<std::mutex> oLock(s_oGeoTreeMapMutex);
        auto iMatchedGeoRegionTree = s_mGeoTrees.find(sSessionName);
        if(iMatchedGeoRegionTree == s_mGeoTrees.end())
            return GeoRegionTreePtr();
        return iMatchedGeoRegionTree->second;
    }

    static SessionNameType suggestSessionName(const GeoRegionMap& mRegions) {
        const GeoRegionUIDArray vGeoRegionUIds = getKeyArrayFromMap(mRegions);
        const GeoRegionUID nID = getParentUID(mRegions.begin(), mRegions.end());
        return std::to_string(nID);
    }

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

std::map<SessionNameType, GeoRegionTreePtr> GeoRegionTreeCacher::s_mGeoTrees;
std::mutex GeoRegionTreeCacher::s_oGeoTreeMapMutex;

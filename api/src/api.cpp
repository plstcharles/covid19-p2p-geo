
#include "api.hpp"

#include <chrono>
#include <cmath>
#include <cstdint>
#include <mutex>
#include <numeric>
#include <sstream>

#include <geos/geom/PrecisionModel.h>
#include <geos/geom/GeometryFactory.h>
#include <geos/geom/Geometry.h>
#include <geos/geom/Point.h>
#include <geos/geom/LinearRing.h>
#include <geos/geom/LineString.h>
#include <geos/geom/Polygon.h>
#include <geos/geom/GeometryCollection.h>
#include <geos/geom/Coordinate.h>
#include <geos/geom/CoordinateSequence.h>
#include <geos/geom/CoordinateArraySequence.h>
#include <geos/geom/IntersectionMatrix.h>
#include <geos/io/WKBReader.h>
#include <geos/io/WKBWriter.h>
#include <geos/io/WKTWriter.h>
#include <geos/index/strtree/STRtree.h>
#include <geos/util/GeometricShapeFactory.h>
#include <geos/geom/util/SineStarFactory.h>
#include <geos/util/GEOSException.h>
#include <geos/util/IllegalArgumentException.h>
#include <geos/operation/linemerge/LineMerger.h>
#include <geos/operation/polygonize/Polygonizer.h>
#include <geos/constants.h>
#include <vector>
#include <sstream>
#include <iomanip>
#include <iostream>
#include <cstdlib>
#include <lz4.h>

#include "H5Cpp.h"


#define DEFAULT_VERBOSE false


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


struct BufferStream: std::streambuf {
    BufferStream(uint8_t* pBegin, uint8_t* pEnd) {
        this->setg((char*)pBegin, (char*)pBegin, (char*)pEnd);
    }
};

template<typename TKey, typename TVal>
std::vector<TKey> getKeyArrayFromMap(const std::map<TKey, TVal>& mMap) {
    std::vector<TKey> vKeys;
    vKeys.reserve(mMap.size());
    for(const auto& oPair : mMap)
        vKeys.push_back(oPair.first);
    return vKeys;
}

template<typename TKey, typename TVal>
std::vector<TVal> getValArrayFromMap(const std::map<TKey, TVal>& mMap) {
    std::vector<TVal> vVals;
    vVals.reserve(mMap.size());
    for(const auto& oPair : mMap)
        vVals.push_back(oPair.second);
    return vVals;
}

template<typename TVal>
std::string writeStatsToString(const std::vector<TVal>& vVals) {
    static_assert(std::is_arithmetic<TVal>::value, "unexpected type for stats");
    assert(!vVals.empty());
    std::stringstream ssStr;
    const double dMin = (double)*std::min_element(vVals.begin(), vVals.end());
    ssStr << "\tmin = " << dMin << "\n";
    const double dMax = (double)*std::max_element(vVals.begin(), vVals.end());
    ssStr << "\tmax = " << dMax << "\n";
    const double dSum = std::accumulate(vVals.begin(), vVals.end(), 0.0);
    const double dMean = dSum / (double)vVals.size();
    ssStr << "\tmean = " << dMean << "\n";
    const double dSqSum = std::inner_product(vVals.begin(), vVals.end(), vVals.begin(), 0.0);
    const double dStdDev = std::sqrt(dSqSum / (double)vVals.size() - dMean * dMean);
    ssStr << "\tstddev = " << dStdDev << "\n";
    return ssStr.str();
}

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

    GeoRegion(const GeoRegionStats& oStats_, GeoRegionUID nUID_, const GeomEnvelope& oEnvelope_,
            Geometry&& pGeometry_=nullptr):
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

    explicit GeoRegionTree(const GeoRegionArray& vRegions,
            size_t nDefaultNodeCapacity=10u, bool bVerboseBuild=DEFAULT_VERBOSE):
            tBuildTimestamp(std::chrono::high_resolution_clock::now()),
            mRegions(getGeoRegionMapFromArray(vRegions)),
            nParentUID(getParentUID(vRegions.begin(), vRegions.end())),
            pEnvelope(getGeomArrayEnvelope(vRegions)),
            oTree(nDefaultNodeCapacity) {
        buildTree(bVerboseBuild);
    }

    explicit GeoRegionTree(
            const GeoRegionMap& mRegions_,
            size_t nDefaultNodeCapacity = 10u, bool bVerboseBuild = DEFAULT_VERBOSE):
            tBuildTimestamp(std::chrono::high_resolution_clock::now()),
            mRegions(mRegions_),
            nParentUID(getParentUID(mRegions.begin(), mRegions.end())),
            pEnvelope(getGeomArrayEnvelope(mRegions)), oTree(nDefaultNodeCapacity) {
        buildTree(bVerboseBuild);
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

    void buildTree(bool bVerboseBuild) {
        // calling this function mutliple times is undefined behavior
        assert(!mRegions.empty());
        for(auto& pRegionPair : mRegions)
            oTree.insert(&pRegionPair.second->oEnvelope, (void*)(&pRegionPair.second));
        oTree.build();
        if(bVerboseBuild) {
            std::chrono::duration<double, std::milli> tBuildTime =
                    std::chrono::high_resolution_clock::now() - tBuildTimestamp;
            std::cout << " tree build time: " << tBuildTime.count() << " msec" << std::endl;
        }
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
            size_t nDefaultNodeCapacity = 10u,
            bool bVerboseBuild = DEFAULT_VERBOSE) {
        assert(!vRegions.empty());
        if(oSessionName == AUTO_SESSION_NAME)
            oSessionName = suggestSessionName(vRegions);
        auto pRegionTree = std::make_shared<GeoRegionTree>(
                vRegions, nDefaultNodeCapacity, bVerboseBuild);
        const std::lock_guard<std::mutex> oLock(s_oGeoTreeMapMutex);
        return *s_mGeoTrees.insert(std::make_pair(oSessionName, pRegionTree)).first;
    }

    static std::pair<SessionNameType, GeoRegionTreePtr> createGeoRegionTree(
            const GeoRegionMap& mRegions,
            SessionNameType oSessionName=AUTO_SESSION_NAME,
            size_t nDefaultNodeCapacity = 10u,
            bool bVerboseBuild = DEFAULT_VERBOSE) {
        assert(!mRegions.empty());
        if(oSessionName == AUTO_SESSION_NAME)
            oSessionName = suggestSessionName(mRegions);
        auto pRegionTree = std::make_shared<GeoRegionTree>(
                mRegions, nDefaultNodeCapacity, bVerboseBuild);
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


template<typename TOut, typename TDataset>
TOut readHDF5Int64Attrib(const TDataset& oH5, const std::string& sAttribName) {
    const H5::Attribute oCountAttrib = oH5.openAttribute(sAttribName);
    const H5::DataType tCountType(H5::PredType::NATIVE_INT64);
    assert(oCountAttrib.getDataType() == tCountType);
    int64_t nCount_int64 = 0;
    oCountAttrib.read(tCountType, &nCount_int64);
    return (TOut)nCount_int64;
}

template<typename TDataset>
std::vector<std::string> readHDF5StrVecAttrib(const TDataset& oH5, const std::string& sAttribName) {
    const H5::Attribute oStrVecAttrib = oH5.openAttribute(sAttribName);
    const H5::DataType oDataType = oStrVecAttrib.getDataType();
    assert(oDataType == H5::StrType(H5::PredType::C_S1, H5T_VARIABLE));
    const H5::DataSpace oDataSpace = oStrVecAttrib.getSpace();
    assert(oDataSpace.getSimpleExtentNdims() == 1);
    hsize_t nStrCount = 0u;
    oDataSpace.getSimpleExtentDims(&nStrCount);
    std::vector<char*> vsRawStrings(nStrCount);  // will need to dealloc manually
    oStrVecAttrib.read(oDataType, (void*)vsRawStrings.data());
    std::vector<std::string> vOutputStrings;
    for(size_t nStrIdx = 0u; nStrIdx < nStrCount; ++nStrIdx) {
        vOutputStrings.emplace_back(vsRawStrings[nStrIdx]);
        delete[] vsRawStrings[nStrIdx];
    }
    return vOutputStrings;
}

GeomArray readWKBGeometries(uint8_t* aBuffer, size_t nBufferSize, size_t nGeomCount) {
    GeomArray vResult;
    vResult.reserve(nGeomCount);
    BufferStream ssBuffer(aBuffer, aBuffer + nBufferSize);
    std::istream ssInputStream(&ssBuffer);
    geos::io::WKBReader oWKBReader;
    for(size_t nGeomIdx = 0u; nGeomIdx < nGeomCount; ++nGeomIdx)
        vResult.push_back(oWKBReader.read(ssInputStream));
    return vResult;
}

GeomArray readHDF5DAWKBGeometries(const H5::H5File& oH5Archive) {
    const hsize_t nDACount = readHDF5Int64Attrib<hsize_t>(oH5Archive, "da_count");
    const H5::DataSet oWKBDataset = oH5Archive.openDataSet("da_wkb");
    const H5::Attribute oWKBByteCountAttrib = oWKBDataset.openAttribute("byte_count");
    const H5::DataType tCountType(H5::PredType::NATIVE_INT64);
    assert(oWKBByteCountAttrib.getDataType() == tCountType);
    int64_t nWKBByteCount_int64 = 0;
    oWKBByteCountAttrib.read(tCountType, &nWKBByteCount_int64);
    const auto nWKBByteCount = (hsize_t)nWKBByteCount_int64;
    H5::DataSpace oWKBSpace = oWKBDataset.getSpace();
    assert(oWKBSpace.getSimpleExtentNdims() == 1);
    hsize_t nWKBDims;
    oWKBSpace.getSimpleExtentDims(&nWKBDims);
    assert(nWKBDims == (hsize_t)nDACount);
    const H5::DataType oWKBType = oWKBDataset.getDataType();
    assert(oWKBType == H5::VarLenType(H5::PredType::NATIVE_UINT8));
    std::unique_ptr<uint8_t[]> abWKBBuffer(new uint8_t[nWKBByteCount]);
    const hsize_t nDAReadCount = 1u;
    hsize_t nLookupIdx = 0u, nWKBOffset = 0u;
    H5::DataSpace oWKBBufferSpace(1, &nDAReadCount);
    oWKBBufferSpace.selectHyperslab(H5S_SELECT_SET, &nDAReadCount, &nLookupIdx);
    for(; nLookupIdx < nDACount; ++nLookupIdx) {
        oWKBSpace.selectHyperslab(H5S_SELECT_SET, &nDAReadCount, &nLookupIdx);
        const hsize_t nWKBMemSize = oWKBDataset.getVlenBufSize(oWKBType, oWKBSpace);
        assert(nWKBMemSize > 0u && nWKBOffset + nWKBMemSize <= nWKBByteCount);
        hvl_t oTempVLenBuffer;
        oWKBDataset.read(
                &oTempVLenBuffer,
                H5::VarLenType(&H5::PredType::NATIVE_UINT8), oWKBBufferSpace, oWKBSpace);
        assert(oTempVLenBuffer.len == nWKBMemSize);
        std::copy(
                (uint8_t*)oTempVLenBuffer.p,
                (uint8_t*)oTempVLenBuffer.p + oTempVLenBuffer.len,
                abWKBBuffer.get() + nWKBOffset);
        H5::DataSet::vlenReclaim(
                &oTempVLenBuffer,
                H5::VarLenType(&H5::PredType::NATIVE_UINT8), oWKBBufferSpace);
        nWKBOffset += nWKBMemSize;
    }
    return readWKBGeometries(abWKBBuffer.get(), nWKBByteCount, nDACount);
}

template<typename TDataset>
GeoRegionStatsArray readHDF5GeoRegionStats(const TDataset& oH5, const std::string& sDatasetName) {
    const std::vector<std::string> vStatNames = readHDF5StrVecAttrib(oH5, "stat_cols");
    assert(vStatNames.size() == GeoRegionStats::s_nExpectedStatCount);
    const auto aExpectedStatNames = GeoRegionStats::getExpectedStatNames();
    assert(std::equal(vStatNames.begin(), vStatNames.end(), aExpectedStatNames.begin()));
    const H5::DataSet oDataset = oH5.openDataSet(sDatasetName);
    const H5::DataType oDataType = oDataset.getDataType();
    assert(oDataType == H5::PredType::NATIVE_FLOAT);
    const H5::DataSpace oDataSpace = oDataset.getSpace();
    assert(oDataSpace.getSimpleExtentNdims() == 2);
    hsize_t anDatasetShape[2] = {0u, 0u};
    oDataSpace.getSimpleExtentDims(anDatasetShape);
    assert(anDatasetShape[0] > 0u && anDatasetShape[1] == GeoRegionStats::s_nExpectedStatCount);
    std::vector<float> vfRawStats(anDatasetShape[0] * anDatasetShape[1]);
    oDataset.read((void*)vfRawStats.data(), oDataType);
    GeoRegionStatsArray vGeoRegionStats;
    vGeoRegionStats.reserve(anDatasetShape[0]);
    for(size_t nStatsIdx = 0u; nStatsIdx < anDatasetShape[0]; ++nStatsIdx) {
        const float* aStatsArray = vfRawStats.data() + nStatsIdx * anDatasetShape[1];
        vGeoRegionStats.emplace_back(aStatsArray);
    }
    return vGeoRegionStats;
}

template<typename TDataset>
GeoRegionUIDArray readHDF5GeoRegionUIDs(const TDataset& oH5, const std::string& sDatasetName) {
    const H5::DataSet oDataset = oH5.openDataSet(sDatasetName);
    const H5::DataType oDataType = oDataset.getDataType();
    assert(oDataType == H5::PredType::NATIVE_UINT32);
    const H5::DataSpace oDataSpace = oDataset.getSpace();
    assert(oDataSpace.getSimpleExtentNdims() == 1);
    hsize_t nDatasetShape = 0u;
    oDataSpace.getSimpleExtentDims(&nDatasetShape);
    assert(nDatasetShape > 0u);
    GeoRegionUIDArray vnUIDs(nDatasetShape);
    oDataset.read((void*)vnUIDs.data(), oDataType);
    return vnUIDs;
}

template<typename TDataset>
std::vector<GeomEnvelope> readHDF5GeomEnvelopes(const TDataset& oH5, const std::string& sDatasetName) {
    const H5::DataSet oDataset = oH5.openDataSet(sDatasetName);
    const H5::DataType oDataType = oDataset.getDataType();
    assert(oDataType == H5::PredType::NATIVE_FLOAT);
    const H5::DataSpace oDataSpace = oDataset.getSpace();
    assert(oDataSpace.getSimpleExtentNdims() == 2);
    hsize_t anDatasetShape[2] = {0u, 0u};
    oDataSpace.getSimpleExtentDims(anDatasetShape);
    assert(anDatasetShape[0] > 0u && anDatasetShape[1] == 4u);
    std::vector<float> vfRawCoords(anDatasetShape[0] * anDatasetShape[1]);
    oDataset.read((void*)vfRawCoords.data(), oDataType);
    std::vector<GeomEnvelope> vEnvelopes;
    vEnvelopes.reserve(anDatasetShape[0]);
    for(size_t nEnvIdx = 0u; nEnvIdx < anDatasetShape[0]; ++nEnvIdx) {
        const float* aEnvArray = vfRawCoords.data() + nEnvIdx * anDatasetShape[1];
        vEnvelopes.emplace_back(aEnvArray[0], aEnvArray[1], aEnvArray[2], aEnvArray[3]);
    }
    return vEnvelopes;
}

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
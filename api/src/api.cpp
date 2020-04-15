//
// Created by perf6 on 4/13/20.
//

#include "api.hpp"

#include <chrono>
#include <cmath>
#include <cstdint>

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

struct BufferStream: std::streambuf {
    BufferStream(uint8_t* pBegin, uint8_t* pEnd) {
        this->setg((char*)pBegin, (char*)pBegin, (char*)pEnd);
    }
};

struct GeomStats {

    GeomStats() :
        fPopulation(std::numeric_limits<float>::quiet_NaN()),
        fDwellings(std::numeric_limits<float>::quiet_NaN()),
        fArea(std::numeric_limits<float>::quiet_NaN()) {}

    GeomStats(float fPopulation_, float fDwellings_, float fArea_):
            fPopulation(fPopulation_),
            fDwellings(fDwellings_),
            fArea(fArea_) {}

    float fPopulation;
    float fDwellings;
    float fArea;
};

using GeomArray = std::vector<geos::geom::Geometry::Ptr>;
using GeomStatsArray = std::vector<GeomStats>;

hsize_t getHDF5GeomCount(const H5::H5File& oH5Archive) {
    const H5::Attribute oDACountAttrib = oH5Archive.openAttribute("da_count");
    const H5::DataType tCountType(H5::PredType::NATIVE_INT64);
    assert(oDACountAttrib.getDataType() == tCountType);
    int64_t nDACount_int64 = 0;
    oDACountAttrib.read(tCountType, &nDACount_int64);
    return (hsize_t)nDACount_int64;
}

GeomArray readWKBGeometries(uint8_t* aBuffer, size_t nBufferSize, size_t nGeomCount) {
    std::vector<std::unique_ptr<geos::geom::Geometry>> vResult;
    vResult.reserve(nGeomCount);
    BufferStream ssBuffer(aBuffer, aBuffer + nBufferSize);
    std::istream ssInputStream(&ssBuffer);
    geos::io::WKBReader oWKBReader;
    for(size_t nGeomIdx = 0u; nGeomIdx < nGeomCount; ++nGeomIdx)
        vResult.push_back(oWKBReader.read(ssInputStream));
    return vResult;
}

GeomArray readHDF5Geometries(const H5::H5File& oH5Archive) {
    const hsize_t nDACount = getHDF5GeomCount(oH5Archive);
    const H5::DataSet oWKBDataset = oH5Archive.openDataSet("wkb");
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

GeomStatsArray readHDF5Stats(const H5::H5File& oH5Archive) {
    const hsize_t nDACount = getHDF5GeomCount(oH5Archive);
    const H5::DataSet oStatsDataset = oH5Archive.openDataSet("stats");
    H5::DataSpace oStatsSpace = oStatsDataset.getSpace();
    assert(oStatsSpace.getSimpleExtentNdims() == 2);
    hsize_t anStatsDims[2];
    oStatsSpace.getSimpleExtentDims(anStatsDims);
    assert(anStatsDims[0] == nDACount && anStatsDims[1] == 3);
    const H5::DataType oStatsType = oStatsDataset.getDataType();
    assert(oStatsType == H5::PredType::NATIVE_FLOAT);
    std::vector<GeomStats> vGeomStats;
    vGeomStats.reserve(nDACount);
    const hsize_t anDAReadCount[2] = {1u, 3u};
    hsize_t anLookupIdxs[2] = {0u, 0u};
    H5::DataSpace oStatsBufferSpace(2, anDAReadCount);
    oStatsBufferSpace.selectHyperslab(H5S_SELECT_SET, anDAReadCount, anLookupIdxs);
    float afStats[3];
    for(anLookupIdxs[0] = 0u; anLookupIdxs[0] < nDACount; ++anLookupIdxs[0]) {
        oStatsSpace.selectHyperslab(H5S_SELECT_SET, anDAReadCount, anLookupIdxs);
        oStatsDataset.read(afStats, oStatsType, oStatsBufferSpace, oStatsSpace);
        // by definition, stats should encode 0=pop, 1=dwellings, 2=area
        vGeomStats.emplace_back(afStats[0], afStats[1], afStats[2]);
    }
    return vGeomStats;
}

double stopwatch_tree_random_query(
        const geos::geom::GeometryCollection::Ptr& pGeomCollection,
        geos::index::strtree::STRtree& oSTRTree
        ) {
    const size_t nTargetGeomIdx = (size_t)std::rand() % pGeomCollection->getNumGeometries();
    const geos::geom::Geometry* pTargetGeom = pGeomCollection->getGeometryN(nTargetGeomIdx);
    const geos::geom::Point::Ptr pTargetGeomCentroid = pTargetGeom->getCentroid();
    const geos::geom::Envelope* pQueryEnvelope = pTargetGeomCentroid->getEnvelopeInternal();
    std::vector<void*> vQueryHits;
    const auto tPreQueryTimestamp = std::chrono::high_resolution_clock::now();
    oSTRTree.query(pQueryEnvelope, vQueryHits);
    const auto tPostQueryTimestamp = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::milli> tQueryTime = tPostQueryTimestamp - tPreQueryTimestamp;
    bool bFoundTargetGeom = false;
    for(auto& pHitData : vQueryHits)
        bFoundTargetGeom |= (*(size_t*)pHitData) == nTargetGeomIdx;
    assert(bFoundTargetGeom);
    return tQueryTime.count();
}

void print_debug_info() {

    // add proper exception handling w/o rethrowing to top level?

    // @@@@ DIRTY CHEAT
    const std::string sHDF5FilePath = "/home/perf6/dev/covid19-p2p-geo/data/da_pop_geometries.hdf5";


    const H5::H5File oH5Archive(sHDF5FilePath, H5F_ACC_RDONLY);
    const GeomStatsArray vStats = readHDF5Stats(oH5Archive);
    double tot_area = 0.0, tot_pop = 0.0;
    for(size_t nIdx = 0u; nIdx < vStats.size(); ++nIdx) {
        if(!std::isnan(vStats[nIdx].fPopulation)) {
            tot_pop += (double)vStats[nIdx].fPopulation;
            tot_area += (double)vStats[nIdx].fArea;
        }
    }
    std::cout << " tot pop: " << tot_pop << std::endl;
    std::cout << " tot area: " << tot_area << std::endl;
    std::cout << " average density: " << tot_pop / tot_area << std::endl;

    geos::geom::GeometryFactory::Ptr pGeomFact = geos::geom::GeometryFactory::create();
    geos::geom::GeometryCollection::Ptr pGeomCollection =
            pGeomFact->createGeometryCollection(readHDF5Geometries(oH5Archive));

    geos::io::WKTWriter oWKTWriter;
    geos::geom::Geometry::Ptr pGlobalROI = pGeomCollection->getEnvelope();
    std::cout << " global roi: " << oWKTWriter.writeFormatted(pGlobalROI.get()) << std::endl;
    std::cout << " geometry count: " << pGeomCollection->getNumGeometries() << std::endl;

    const size_t nDefaultNodeCapacity = 10u;
    geos::index::strtree::STRtree oGlobalTree(nDefaultNodeCapacity);
    std::vector<size_t> vnGeomIdxsLookupBuffer(pGeomCollection->getNumGeometries());
    for(size_t nGeomIdx = 0u; nGeomIdx < pGeomCollection->getNumGeometries(); ++nGeomIdx) {
        vnGeomIdxsLookupBuffer[nGeomIdx] = nGeomIdx;
        oGlobalTree.insert(
                pGeomCollection->getGeometryN(nGeomIdx)->getEnvelopeInternal(),
                &vnGeomIdxsLookupBuffer[nGeomIdx]);
    }
    const auto tPreBuildTimestamp = std::chrono::high_resolution_clock::now();
    oGlobalTree.build();
    const auto tPostBuildTimestamp = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::milli> tBuildTime = tPostBuildTimestamp - tPreBuildTimestamp;
    std::cout << " build time: " << tBuildTime.count() << " msec" << std::endl;

    const size_t nRandomQueryCount = 10000u;
    double dTotQueryTime = 0.0;
    for(size_t nQueryIdx = 0u; nQueryIdx < nRandomQueryCount; ++nQueryIdx)
        dTotQueryTime += stopwatch_tree_random_query(pGeomCollection, oGlobalTree);
    const double dAvgQueryTime = dTotQueryTime / nRandomQueryCount;
    std::cout << " average query time: " << dAvgQueryTime << " msec " <<
            "(" << nRandomQueryCount << " trials" << ")" << std::endl;
}

// multi-level tree idea w/ overlapping recursive region lookup
//   - get global envelope, split into grid of blocks (with some overlap)
//       - size of blocks should be based on "how wide" a person can move without new data query
//       - idea for smart blocks: use VoronoiDiagramBuilder, get diagram polygons, buffer them up
//   - for each block, get the list of all geometries (DAs/sub-blocks) that intersect the block
//       - some geometries will be present in multiple blocks
//       - the geometries could be sub-blocks from another gridding level (recursive if needed)
//   - compute envelope for all geometries, build STRTree, prepare for spatial lookups
//       - user data tied to each envelope is geometry (polygon) + DAUID if available
//       - after making a query, look at each resulting geometry and check for fine-grained inters
//       - last tree level should only contain DAs, and guarantee a DAUID no matter what
//   - ultimate result of an API query should be one ID per geometry that was hit
//       - e.g.:  BLOCK#2365236; PT#30; CD#3058; DAUID#30583099;
//

// @@@@ VERSION 0: use STRtree/quadtree from geos
//   - will use a lot of memory, hard to serialize in a smart way
// @@@@ VERSION 1: use a custom quad tree
//   - maybe not as efficient for lookups as STRtree, but can be serialized well
// @@@@ VERSION 2: reimplement a STRtree (based on geos impl) for proper serialization
//   - will take a lot of dev/test time

/*

#  quad tree impl notes:
#  - use internal vector(s) for geoms overlapping with neighbors (should have overlap w/ 2nd internal from children; add ptr to neighb)
#  - use 2nd internal vector for geoms fully contained in node, but covering multiple children
#  - test by randomly casting in digital boundary --- should always find a geometry
#  - query w/ buffer around GPS coord, return average of whatev-values stored in geometries?
#      (compute intersection w/ bbox instead of just a point)
#      buffer size could be determined based on GPS coord accuracy
#  - skeletize all 400-700 pop geometries down to 100-200? (hyperparam)

*/
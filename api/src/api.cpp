//
// Created by perf6 on 4/13/20.
//

#include "../include/api.hpp"

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

void api_hello() {
    // add proper exception handling w/o rethrowing to top level?

    const std::string sHDF5FilePath = "/home/perf6/dev/covid19-p2p-geo/data/da_pop_geometries.hdf5";
    H5::H5File oH5Archive(sHDF5FilePath, H5F_ACC_RDONLY);
    const H5::Attribute oDACountAttrib = oH5Archive.openAttribute("da_count");
    const H5::DataType tCountType(H5::PredType::NATIVE_INT64);
    assert(oDACountAttrib.getDataType() == tCountType);
    int64_t nDACount_int64 = 0;
    oDACountAttrib.read(tCountType, &nDACount_int64);
    const auto nDACount = (hsize_t)nDACount_int64;
    H5::DataSet oWKBDataset = oH5Archive.openDataSet("wkb");
    const H5::Attribute oWKBByteCountAttrib = oWKBDataset.openAttribute("byte_count");
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

    H5::DataSet oStatsDataset = oH5Archive.openDataSet("stats");
    H5::DataSpace oStatsSpace = oStatsDataset.getSpace();
    assert(oStatsSpace.getSimpleExtentNdims() == 2);
    hsize_t anStatsDims[2];
    oStatsSpace.getSimpleExtentDims(anStatsDims);
    assert(anStatsDims[0] == nDACount && anStatsDims[1] == 3);
    const H5::DataType oStatsType = oStatsDataset.getDataType();
    assert(oStatsType == H5::PredType::NATIVE_FLOAT);
    // by definition, stats should encode 0=pop, 1=dwellings, 2=area

    std::vector<hsize_t> vnWKBMemSizes((size_t)nDACount);
    std::vector<float> vfArea((size_t)nDACount), vfPopCount((size_t)nDACount);
    // weird alloc below skips buffer zero-init
    std::unique_ptr<uint8_t[]> abWKBBuffer(new uint8_t[nWKBByteCount]);
    const hsize_t anDAReadCount[2] = {1u, 3u};
    hsize_t anLookupIdxs[2] = {0u, 0u};
    hsize_t nWKBOffset = 0u;
    H5::DataSpace oWKBBufferSpace(1, anDAReadCount);
    H5::DataSpace oStatsBufferSpace(2, anDAReadCount);
    oWKBBufferSpace.selectHyperslab(H5S_SELECT_SET, anDAReadCount, anLookupIdxs);
    oStatsBufferSpace.selectHyperslab(H5S_SELECT_SET, anDAReadCount, anLookupIdxs);
    float afStats[3];
    size_t nTotMemSize = 0u;
    for(; anLookupIdxs[0] < nDACount; ++anLookupIdxs[0]) {
        oWKBSpace.selectHyperslab(H5S_SELECT_SET, anDAReadCount, anLookupIdxs);
        oStatsSpace.selectHyperslab(H5S_SELECT_SET, anDAReadCount, anLookupIdxs);
        const hsize_t nWKBMemSize = oWKBDataset.getVlenBufSize(oWKBType, oWKBSpace);
        auto something = oWKBDataset.getVarLenType();
        nTotMemSize += nWKBMemSize;
        assert(nWKBMemSize > 0u && nWKBOffset + nWKBMemSize <= nWKBByteCount);
        //std::cout << "#" << anLookupIdxs[0] << "/" << nDACount <<
        //        "  (size=" << nWKBMemSize << ")  " << nWKBOffset <<
        //        "/" << nWKBByteCount << std::endl;
        hvl_t tmp;
        oWKBDataset.read(&tmp,
                H5::VarLenType(&H5::PredType::NATIVE_UINT8), oWKBBufferSpace, oWKBSpace);
        assert(tmp.len == nWKBMemSize);
        std::copy((uint8_t*)tmp.p, (uint8_t*)tmp.p + tmp.len, abWKBBuffer.get() + nWKBOffset);
        vnWKBMemSizes[anLookupIdxs[0]] = nWKBMemSize;
        nWKBOffset += nWKBMemSize;
        oStatsDataset.read(afStats, oStatsType, oStatsBufferSpace, oStatsSpace);
        vfPopCount[anLookupIdxs[0]] = afStats[0];
        vfArea[anLookupIdxs[0]] = afStats[2];
    }

    double tot_area = 0.0, tot_pop = 0.0;
    for(size_t nIdx = 0u; nIdx < vfPopCount.size(); ++nIdx) {
        if(!std::isnan(vfPopCount[nIdx])) {
            tot_pop += (double)vfPopCount[nIdx];
            tot_area += (double)vfArea[nIdx];
        }
    }

    std::cout << " tot pop = " << tot_pop << std::endl;
    std::cout << " tot area = " << tot_area << std::endl;
    std::cout << " average density = " << tot_pop / tot_area << std::endl;
}

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
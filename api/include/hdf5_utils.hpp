#pragma once

#include "H5Cpp.h"

#include "regions.hpp"


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

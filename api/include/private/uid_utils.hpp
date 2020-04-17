
/// contains region unique identifier (UID)-related utilities

#pragma once

#include "private/generic_utils.hpp"

// NOTES:
//   we use GEO CODES for provincies/territories, census divisions, and dissemination areas
//   (this means they can only be 2-digit, 4-digit, or 8-digit)
//   (UID = 0 is reserved as the "global region")

using GeoRegionUID = uint32_t;
using GeoRegionUIDArray = std::vector<GeoRegionUID>;
using SessionNameType = std::string;

#define AUTO_SESSION_NAME std::string()
#define GLOBAL_REGION_UID GeoRegionUID(0u)
#define GLOBAL_REGION_STR std::to_string(GLOBAL_REGION_UID)

/// returns whether the UID corresponds to a province/territory
bool isProvinceOrTerritory(GeoRegionUID nUID) {
    return nUID >= 10u && nUID < 100u;
}

/// returns whether the UID corresponds to a province/territory
bool isProvinceOrTerritory(const SessionNameType& sUID) {
    return sUID.size() == 2u;
}

/// returns whether the UID corresponds to a census division
bool isCensusDivision(GeoRegionUID nUID) {
    return nUID >= 1000u && nUID < 10000u;
}

/// returns whether the UID corresponds to a census division
bool isCensusDivision(const SessionNameType& sUID) {
    return sUID.size() == 4;
}

/// returns whether the UID corresponds to a dissemination area
bool isDisseminationArea(GeoRegionUID nUID) {
    return nUID >= 10000000u && nUID < 100000000u;
}

/// returns whether the UID corresponds to a dissemination area
bool isDisseminationArea(const SessionNameType& sUID) {
    return sUID.size() == 8;
}

/// returns the unique identifier (UID) of the parent region specified via UID
GeoRegionUID getParentUID(GeoRegionUID nUID) {
    if(isDisseminationArea(nUID))
        return nUID / 10000u;
    if(isCensusDivision(nUID))
        return nUID / 100u;
    if(isProvinceOrTerritory(nUID))
        return GLOBAL_REGION_UID;
    return GLOBAL_REGION_UID;
}

/// returns the unique identifier (UID) of the parent region specified via UID
SessionNameType getParentUID(const SessionNameType& sUID) {
    if(isDisseminationArea(sUID))
        return sUID.substr(0u, 4u);
    if(isCensusDivision(sUID))
        return sUID.substr(0u, 2u);
    if(isProvinceOrTerritory(sUID))
        return GLOBAL_REGION_STR;
    return GLOBAL_REGION_STR;
}

/// returns the unique identifier (UID) of the parent region specified via a UID/value pair
template<typename TKey, typename TVal>
GeoRegionUID getParentUID(const std::pair<TKey, TVal>& oPair) {
    return getParentUID(oPair.first); // assume first key is UID-related
}

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

/// returns the unique identifier (UID) of the parent region specified via GeoRegion
template<typename TGeoRegion>
GeoRegionUID getParentUID(const TGeoRegion& oRegion) {
    return getParentUID(oRegion.nUID);
}

/// returns the unique identifier (UID) of the parent region specified via GeoRegionPtr
template<typename TGeoRegion>
GeoRegionUID getParentUID(const std::shared_ptr<TGeoRegion>& pRegion) {
    if(!pRegion)
        return GLOBAL_REGION_UID; // if the region is empty/unavailable, return global ID
    return getParentUID(pRegion->nUID);
}

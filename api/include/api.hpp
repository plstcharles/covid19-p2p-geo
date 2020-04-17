#pragma once

#include <string>

// add namespace?

/// prepares the region lookup datastructures near the given GPS coordinates
void prepareRegionsNear(double dLatitude, double dLongitude, const std::string& sDataRootPath);

/// returns the dissemination area UID associated with the given GPS coordinates
std::string fetchRegionUID(double dLatitude, double dLongitude, const std::string& sDataRootPath);

/// runs a 1000-trial random building/querying routine on all regions
void testRandomBuildAndQueries(const std::string& sDataRootPath);

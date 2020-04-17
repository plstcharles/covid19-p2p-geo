#pragma once

#include <string>

#define TIMEOUT_12_HOURS (((60) * 60) * 12)

// add namespace?

/// prepares the region lookup datastructures near the given GPS coordinates
void prepareNear(double dLatitude, double dLongitude, const std::string& sDataRootPath);

/// returns the region UID associated with the given GPS coordinates
std::string fetchUID(double dLatitude, double dLongitude, const std::string& sDataRootPath);

/// returns the current number of cached region subtrees
uint32_t getCachedRegionTreeCount();

/// releases the cache associated with regions that has not been accessed in some time
void releaseUnusedCache(double dTimeout_seconds=TIMEOUT_12_HOURS);

/// runs a 1000-trial random building/querying routine on all regions
void testRandomBuildAndQueries(const std::string& sDataRootPath);

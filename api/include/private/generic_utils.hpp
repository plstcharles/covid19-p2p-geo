
/// contains generic C++ utility functions/defines

#pragma once

#include <algorithm>
#include <cassert>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <cstdint>
#include <map>
#include <memory>
#include <mutex>
#include <numeric>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <sys/stat.h>
#include <vector>

//#include <lz4.h> // might add compression utils later


struct BufferStream: std::streambuf {
    BufferStream(uint8_t* pBegin, uint8_t* pEnd) {
        this->setg((char*)pBegin, (char*)pBegin, (char*)pEnd);
    }
};

inline bool checkPathExists(const std::string& sFilePath) {
    struct stat oBuffer;
    return (stat(sFilePath.c_str(), &oBuffer) == 0);
}

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

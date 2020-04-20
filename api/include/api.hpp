#pragma once

#include <string>

#define ONE_WEEK_IN_SECONDS (((60) * 60) * 24 * 7)

/*
 * Note on usage:
 *   All functions below *should* be thread-safe EXCEPT for the "prepareNear" function
 *   since it returns a pointer to an internal buffer. That buffer should remain valid
 *   until the next call to "prepareNear".
 *
 *   The "data root path" (acDataRootPath) argument expected by several functions should
 *   be the path to a folder that contains "cd.hdf5", i.e. the preprocessed census division
 *   boundaries archive, and a subdirectory named "divisions". This subdirectory should
 *   contain the subdivided region archives named as "XXXX.hdf5", where XXXX is the UID(s)
 *   of the subdivision(s) as returned by "prepareNear".
 *
 *   If the census division archive ("cd.hdf5") is missing or if a required subdivision
 *   archive ("divisions/XXXX.hdf5") is missing during a query, an exception will be thrown.
 *
 *   A typical usage sequence for this API is detailed below:
 *     1) initialize a data root directory: download "cd.hdf5" & create a "divisions" subdirectory
 *     2) given a set of lat/lon coordinates, determine what subdivision(s) must be downloaded
 *        by calling the "prepareNear" function (it will return a string of file names).
 *     3) decompose the returned (comma-separated) string into filenames, download each
 *        archive ("XXXX.hdf5", "YYYY.hdf5", ...), and place them in the "divisions" folder
 *     4) call the "prepareNear" function with the same arguments as before to confirm that
 *        all the necessary regions are available (it should return an empty string)
 *     5) call "fetchUID" to get the final UID of the discretized region in which the lat/long
 *        coordinates fall.
 *     x) repeat steps 2-5 for all new GPS coordinates.
 *     y) once in a while: call "releaseUnusedCache" to cleanup memory allocations.
 *
 */

/// Prepares the region lookup data structures near the given GPS coordinates, returning
/// a c-style string of comma-separated divisions UIDs to download externally, or an empty
/// string if all regions are already downloaded/available. The provided data root path MUST
/// contain the census division archive ("cd.hdf5"), otherwise this function will throw.
const char* prepareNear(double dLatitude, double dLongitude, const char* acDataRootPath);

/// Returns the fine-grained (discretized) region UID (in integer format) for the given GPS
/// location. This UID should be a 10-digit number that corresponds to a census dissemination
/// area. This function assumes that the HDF5 archives with filenames returned by calling the
/// "prepareNear" function with the same coordinates are already downloaded and placed in
/// the data root's "divisions" subdirectory, otherwise this function will throw.
uint32_t fetchUID(double dLatitude, double dLongitude, const char* acDataRootPath);

/// Returns the current number of cached region subtrees (for memory tracking purposes).
uint32_t getCachedRegionTreeCount();

/// Releases the cache associated with regions that has not been accessed in some time.
void releaseUnusedCache(double dTimeout_seconds=ONE_WEEK_IN_SECONDS);

/// Runs a 1000-trial random building/querying routine on all regions (for testing purposes).
/// This function will fail if all subdivision archives are not already available in the
/// provided data root directory.
void testRandomBuildAndQueries(const char* acDataRootPath);

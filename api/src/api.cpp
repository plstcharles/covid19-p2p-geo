//
// Created by perf6 on 4/13/20.
//

#include "../include/api.hpp"

void api_hello() {
    std::cout << "Hello, World!" << std::endl;
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
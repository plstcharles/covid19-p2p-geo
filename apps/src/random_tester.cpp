
/// dummy testing app that will run random region tree builds & queries via the API

#include <iostream>

#include "api.hpp"

int main(int argc, char* argv[]) {
    // assumes the first and only argument passed in is the data root directory path
    if(argc < 2) {
        std::cout << " ERROR: missing command line argument (data root path)" << std::endl;
        return 1;
    }
    testRandomBuildAndQueries(argv[1]);
    return 0;
}

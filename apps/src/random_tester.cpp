
#include "api.hpp"

int main(int argc, char* argv[]) {
    // assumes the first and only argument passed in is the data root directory path
    testRandomBuildAndQueries(argv[1]);
    return 0;
}

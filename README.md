# Discrete Geographic Region Projection Engine

This is a C++ library based on GEOS. It relies on data from the 2016 census to project
GPS coordinates across Canada and returns the corresponding dissemination area. These
areas typically contain 400-700 people, but may stretch over large regions in sparsely
populated areas.

The build system should be fully self-contained. In other words, CMake will automatically
download and build all dependencies of the library. A typical build procedure is:
```
  <PROJECT_ROOT> $  mkdir build
  <PROJECT_ROOT> $  cd build
       .../build $  cmake ../
       .../build $  cmake --build .
```
To use mutiple workers in the build process, you can use:
```
       .../build $  cmake --build . -j X
```
where `X` is roughly the number of workers.

Once built, you should have `bin` and `lib` subdirectories in the output directory. These
should contain CLI apps and API libraries, respectively.

To download and preprocess the census data, a set of scripts is provided in the `preproc`
directory.

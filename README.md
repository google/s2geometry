# S2 Geometry Library

## Overview

This is a package for manipulating geometric shapes. Unlike many geometry
libraries, S2 is primarily designed to work with _spherical geometry_, i.e.,
shapes drawn on a sphere rather than on a planar 2D map. This makes it
especially suitable for working with geographic data.

If you want to learn more about the library, start by reading the
[overview](http://s2geometry.io/about/overview) and [quick start
document](http://s2geometry.io/devguide/cpp/quickstart), then read the
introduction to the [basic types](http://s2geometry.io/devguide/basic_types).

S2 documentation can be found on [s2geometry.io](http://s2geometry.io).

## API/ABI Stability

Note that all [releases](https://github.com/google/s2geometry/releases) are
version 0.x, so there are
[no API or ABI stability guarantees](https://semver.org/#spec-item-4). Starting
with 1.0 we will adhere to [SemVer](https://semver.org/) and follow the
[Google OSS breaking change policy](https://opensource.google/documentation/policies/library-breaking-change)

The Python API is particularly unstable, and it is planned that the SWIGged
API will be replaced by a pybind11 version with more Pythonic names and more
complete functionality.

## Requirements for End Users

*   We aim to support all platforms supported by the
    [Google foundational C++ support policy](https://opensource.google/documentation/policies/cplusplus-support)
*   [CMake](http://www.cmake.org/) >= 3.22
*   A C++ compiler with C++17 support, such as
    [g++ >= 7.5](https://gcc.gnu.org/) or
    [clang >= 14.0.0](https://clang.llvm.org/)
*   [Abseil](https://github.com/abseil/abseil-cpp) LTS
    [`20250814`](https://github.com/abseil/abseil-cpp/releases/tag/20250814.1)
    (standard library extensions). This exact version must be used.
*   [OpenSSL](https://github.com/openssl/openssl) (for its bignum library)

On Ubuntu, all of these other than abseil can be installed via apt-get:

```
sudo apt-get install cmake libssl-dev
```

Otherwise, you may need to install some from source.

Currently, Abseil must always be installed from source.  See the use of
`-DCMAKE_PREFIX_PATH` in the [build instructions below](#building).
This is likely to change.

On macOS, use [MacPorts](http://www.macports.org/) or
[Homebrew](http://brew.sh/).  For MacPorts:

```
sudo port install cmake openssl
```

## Build and Install

You may either download the source as a ZIP archive, or [clone the git
repository](https://help.github.com/articles/cloning-a-repository/).

### Via ZIP archive

Download [ZIP file](https://github.com/google/s2geometry/archive/master.zip)

```
cd [parent of directory where you want to put S2]
unzip [path to ZIP file]/s2geometry-master.zip
cd s2geometry-master
```

### Via `git clone`

```
cd [parent of directory where you want to put S2]
git clone https://github.com/google/s2geometry.git
cd s2geometry
```

### Building

First, [install Abseil](https://github.com/abseil/abseil-cpp/blob/master/CMake/README.md#traditional-cmake-set-up).
It must be configured with `-DCMAKE_POSITION_INDEPENDENT_CODE=ON`.
s2geometry must be configured to use the came C++ version that
abseil uses.  The easiest way to achieve this is to pass
`-DCMAKE_CXX_STANDARD=17` to `cmake` when compiling both abseil and
s2geometry.

From the appropriate directory depending on how you got the source:

```
mkdir build
cd build
# You can use -DBUILD_TESTS=no to skip tests.
# Use the same CMAKE_CXX_STANDARD value that was used with absl.
cmake -DBUILD_TESTS=yes -DCMAKE_PREFIX_PATH=/path/to/absl/install -DCMAKE_CXX_STANDARD=17 ..
make -j $(nproc)
make test ARGS="-j$(nproc)"  # If -DBUILD_TESTS=yes was used above.
sudo make install
```

On macOS, `sysctl -n hw.logicalcpu` is the equivalent of `nproc`.

Disable building of shared libraries with `-DBUILD_SHARED_LIBS=OFF`.

## Python

If you want the Python interface, you will also need:

* [SWIG 4](https://github.com/swig/swig) (for Python support, optional)

which can be installed via

```
sudo apt-get install swig
```

or on macOS:

```
sudo port install swig
```
Version 4.0 is required, but it should be easy to make it work 3.0 or probably
even 2.0.

## Other S2 implementations

* [Go](https://github.com/golang/geo) (Approximately 40% complete.)
* [Java](https://github.com/google/s2-geometry-library-java)

## Disclaimer

This is not an official Google product.

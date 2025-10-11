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

## Requirements for End Users using Bazel

## Build

This bazel build requires 7.0.0 – bzlmod default. Builds were tested using 
C++20 as set in .bazelrc. This setup relies on abseil-cpp, boringssl, and 
googletest from the bazel central repository as set in MODULE.bazel.

To build and test using bazel, from within s2geometry/src, run:

`bazel test "//:*"`

To build the libary without testing, from within s2geometry/src, run:

`bazel build //:s2`

## Status

All tests enumerated in the BUILD.bazel file pass on x86 machines. On 
Apple M1, s2loop_measures_test fails due to a 6% excess accumlated error. 
This issue may require revision of boringssl or exactfloat.

## Requirements for End Users using CMake

*   We aim to support all platforms supported by the
    [Google foundational C++ support policy](https://opensource.google/documentation/policies/cplusplus-support)
*   [CMake](http://www.cmake.org/) >= 3.22
*   A C++ compiler with C++17 support, such as
    [g++ >= 7.5](https://gcc.gnu.org/) or
    [clang >= 14.0.0](https://clang.llvm.org/)
*   [Abseil](https://github.com/abseil/abseil-cpp) >= LTS
    [`20250814`](https://github.com/abseil/abseil-cpp/releases/tag/20250814.1)
    (standard library extensions)
*   [googletest testing framework >= 1.10](https://github.com/google/googletest)
    (to build tests and example programs, optional)

On Ubuntu, all of these other than abseil can be installed via apt-get:

```
sudo apt-get install cmake googletest libssl-dev
```

abseil-cpp may need to be installed from source if an LTS release is not
packaged for the platform.  See the use of `-DCMAKE_PREFIX_PATH` in the
[build instructions below](#building).

On macOS, use [MacPorts](http://www.macports.org/) or
[Homebrew](http://brew.sh/).  For MacPorts:

```
sudo port install cmake abseil gtest openssl
```

then use

```
cmake -DGOOGLETEST_ROOT=/opt/local/src -DCMAKE_PREFIX_PATH=/opt/local ..
```

in the build instructions below.

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
s2geometry must be configured to use the same C++ version that
abseil uses.  The easiest way to achieve this is to pass
`-DCMAKE_CXX_STANDARD=17` to `cmake` when compiling both abseil and
s2geometry.

From the appropriate directory depending on how you got the source:

```
mkdir build
cd build
# You can omit -DGOOGLETEST_ROOT to skip tests; see above for macOS.
# Use the same CMAKE_CXX_STANDARD value that was used with absl.
cmake -DGOOGLETEST_ROOT=/usr/src/googletest -DCMAKE_PREFIX_PATH=/path/to/absl/install -DCMAKE_CXX_STANDARD=17 ..
make -j $(nproc)
make test ARGS="-j$(nproc)"  # If GOOGLETEST_ROOT specified above.
sudo make install
```

On macOS, `sysctl -n hw.logicalcpu` is the equivalent of `nproc`.

Disable building of shared libraries with `-DBUILD_SHARED_LIBS=OFF`.

Enable the python interface with `-DWITH_PYTHON=ON`.

# For Testing

If BUILD_TESTS is 'on' (the default), then OpenSSL must be available to build
some tests:

*   [OpenSSL](https://github.com/openssl/openssl) (for its bignum library)

If OpenSSL is installed in a non-standard location set `OPENSSL_ROOT_DIR`
before running configure, for example on macOS:
```
OPENSSL_ROOT_DIR=/opt/homebrew/Cellar/openssl@3/3.1.0 cmake -DCMAKE_PREFIX_PATH=/opt/homebrew -DCMAKE_CXX_STANDARD=17
```

## Installing

From `build` subdirectory:

```
make install
```

Prefix it with `sudo` if needed:

```
sudo make install
```

_NOTE_: There is not `uninstall` target but `install_manifest.txt` may be helpful.

All files will be installed at location specified in `CMAKE_INSTALL_PREFIX` variable.

Several suffix variables used for some file groups:

Variable | Default | Description
-------- | ------- | -----------
`CMAKE_INSTALL_INCLUDEDIR` | `include` | For header files
`CMAKE_INSTALL_BINDIR`     | `bin`     | For executables and `*.dll` files on `DLL`-based platforms
`CMAKE_INSTALL_LIBDIR`     | `lib`     | For library files (`*.so`, `*.a`, `*.lib` etc)

If needed set this variables on command line as `cmake` arguments with `-D` prefix or edit from `build` subdirectory:

```
make edit_cache
```

For more info read: [The CMake Cache](https://cmake.org/cmake/help/latest/guide/user-interaction/index.html#the-cmake-cache).

## Python

If you want the Python interface, you need to run cmake using
`-DWITH_PYTHON=ON`. You will also need to install the following dependencies:

* [SWIG 4](https://github.com/swig/swig) (for Python support, optional)
* python3-dev (for Python support, optional)

which can be installed via

```
sudo apt-get install swig python3-dev
```

or on macOS:

```
sudo port install swig
```
Version 4.0 is required, but it should be easy to make it work 3.0 or probably
even 2.0.

Python 3 is required.

### Creating wheels
First, make a virtual environment and install `build` into it:
```
python3 -m venv venv
source venv/bin/activate
pip install build
```

Then build the wheel:
```
python -m build
```

The resulting wheel will be in the `dist` directory.


## Other S2 implementations

* [Go](https://github.com/golang/geo) (Approximately 40% complete.)
* [Java](https://github.com/google/s2-geometry-library-java)
* [Kotlin](https://github.com/Enovea/s2-geometry-kotlin) (Complete except binary serialization)

## Disclaimer

This is not an official Google product.

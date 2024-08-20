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
[no API or ABI stability guarantees](https://semver.org/#spec-item-4).
Starting with 1.0 we will adhere to [SemVer](https://semver.org/).

The Python API is particularly unstable, and it is planned that the SWIGged
API will be replaced by a pybind11 version with more Pythonic names and more
complete functionality.

## Requirements for End Users

* [CMake](http://www.cmake.org/)
* A C++ compiler with C++14 support, such as [g++ >= 5](https://gcc.gnu.org/)
* [Abseil](https://github.com/abseil/abseil-cpp) >= LTS
  [`20240116`](https://github.com/abseil/abseil-cpp/releases/tag/20240116.1)
  (standard library extensions)
* [OpenSSL](https://github.com/openssl/openssl) (for its bignum library)
* [googletest testing framework >= 1.10](https://github.com/google/googletest)
  (to build tests and example programs, optional)

On Ubuntu, all of these other than abseil can be installed via apt-get:

```
sudo apt-get install cmake googletest libssl-dev
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

Do not install `gtest` from MacPorts; instead download [release
1.10.0](https://github.com/google/googletest/releases/tag/release-1.10.0), unpack,
and substitute

```
cmake -DGOOGLETEST_ROOT=/...absolute path to.../googletest-release-1.10.0 ..
```

in the build instructions below.

Thorough testing has only been done on Ubuntu 14.04.3 and macOS 10.12.

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
`-DCMAKE_CXX_STANDARD=14` (or `-DCMAKE_CXX_STANDARD=17`) to `cmake`
when compiling both abseil and s2geometry.

From the appropriate directory depending on how you got the source:

```
mkdir build
cd build
# You can omit -DGOOGLETEST_ROOT to skip tests; see above for macOS.
# Use the same CMAKE_CXX_STANDARD value that was used with absl.
cmake -DGOOGLETEST_ROOT=/usr/src/googletest -DCMAKE_PREFIX_PATH=/path/to/absl/install -DCMAKE_CXX_STANDARD=14 ..
make -j $(nproc)
make test ARGS="-j$(nproc)"  # If GOOGLETEST_ROOT specified above.
sudo make install
```

On macOS, `sysctl -n hw.logicalcpu` is the equivalent of `nproc`.

Disable building of shared libraries with `-DBUILD_SHARED_LIBS=OFF`.

Enable the python interface with `-DWITH_PYTHON=ON`.

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
First, make a virtual environment and install `cmake_build_extension` and `wheel`
into it:
```
python3 -m venv venv
source venv/bin/activate
pip install cmake_build_extension wheel
```

Then build the wheel:
```
python setup.py bdist_wheel
```

The resulting wheel will be in the `dist` directory.

> If OpenSSL is in a non-standard location make sure to set `OPENSSL_ROOT_DIR` 
> when calling `setup.py`, see above for more information.

## Other S2 implementations

* [Go](https://github.com/golang/geo) (Approximately 40% complete.)
* [Java](https://github.com/google/s2-geometry-library-java)
* [Kotlin](https://github.com/Enovea/s2-geometry-kotlin) (Complete except binary serialization)

## Disclaimer

This is not an official Google product.

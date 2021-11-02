# S2 Geometry Library

[![Build Status](https://travis-ci.org/google/s2geometry.svg?branch=master)](https://travis-ci.org/google/s2geometry)

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

## Requirements for End Users

* [CMake](http://www.cmake.org/)
* A C++ compiler with C++11 support, such as [g++](https://gcc.gnu.org/)
  \>= 4.7.
* [OpenSSL](https://github.com/openssl/openssl) (for its bignum library)
* [gflags command line flags](https://github.com/gflags/gflags), optional
* [glog logging module](https://github.com/google/glog), optional
* [googletest testing framework](https://github.com/google/googletest)
  (to build tests and example programs, optional)

On Ubuntu, all of these can be installed via apt-get:

```
sudo apt-get install cmake libgflags-dev libgoogle-glog-dev libgtest-dev libssl-dev
```
## Build in Alpha docker image

To avoit the absl conflicts, this branch version can only build in Alpha image using bazel.
We have removed the cmake.

```
Make sure to build it in Alpha env since boringssl
baezel build src:all
```

## Other S2 implementations

* [Go](https://github.com/golang/geo) (Approximately 40% complete.)
* [Java](https://github.com/google/s2-geometry-library-java) (Older version;
  last updated in 2011.)

## Disclaimer

This is not an official Google product.

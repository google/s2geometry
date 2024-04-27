
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

## Requirements for End Users

## Build

This bazel build version uses version 7.1.1 as noted in .bazelversion. Builds 
were tested using C++20 as set in .bazelrc. This setup relies on abseil-cpp, boringssl, and googletest from the bazel central repository as set in MODULE.bazel.

To build and test using bazel, run:

`bazel test "//:*"`

To build the libary without testing, run:

`bazel build //:s2lib`

## Status

All tests enumerated in the BUILD.bazel file pass using an x86 machine. On Apple M1, s2loop_measures_test fails due to a 6% excess accumlated error. It's unclear if this test would pass if boringssl were replaced with openssl.

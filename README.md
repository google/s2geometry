S2 Geometry Library
===================

Overview
--------

This is a package for manipulating geometric shapes.  It is rather
heavily weighted towards spherical geometry because currently it is
mainly used for geographic data.  Other types of geometric primitives
will be added as necessary.

Please see [the documentation](doc/intro.md) for more information and examples.

Requirements for End Users
--------------------------

* [CMake](http://www.cmake.org/)
* A C++ compiler with minimal C++11 support, such as [gcc](https://gcc.gnu.org/)
* [gflags command line flags](https://github.com/gflags/gflags)
* [glog logging module](https://github.com/google/glog)
* [OpenSSL](https://github.com/openssl/openssl) (for its bignums)
* [googletest testing framework](https://github.com/google/googletest)
  (to build tests, optional)

On Ubuntu, all of these can be installed via apt-get.  Otherwise, you may need
to install some from source.

Thorough testing has only been done on Ubuntu 14.04.3.  Several tests are
known to fail on OS X (patches appreciated!).

Build and Install
-----------------

```
cd $S2_DIR/geometry
mkdir build
cd build
cmake -DGTEST_ROOT=/usr/src/gtest ..  # Omit -DGTEST_ROOT to skip tests.
make
make test  # If GTEST_ROOT specified above.
sudo make install
```

Who Is Using S2?
----------------

* [Google Maps](https://www.google.com/maps)

Send us a PR to add yourself here.

Disclaimer
----------
This is not an official Google product.

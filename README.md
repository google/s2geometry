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
* A C++ compiler with C++11 support, such as [g++](https://gcc.gnu.org/)
  \>= 4.7.
* [gflags command line flags](https://github.com/gflags/gflags)
* [glog logging module](https://github.com/google/glog)
* [OpenSSL](https://github.com/openssl/openssl) (for its bignums)
* [googletest testing framework](https://github.com/google/googletest)
  (to build tests, optional)
* A POSIX system (for pthreads and getrusage).

On Ubuntu, all of these can be installed via apt-get.  Otherwise, you may need
to install some from source.

Thorough testing has only been done on Ubuntu 14.04.3.

Build and Install
-----------------

```
cd $S2_DIR
mkdir build
cd build
cmake -DGTEST_ROOT=/usr/src/gtest ..  # Omit -DGTEST_ROOT to skip tests.
make
make test  # If GTEST_ROOT specified above.
sudo make install
```

Disclaimer
----------
This is not an official Google product.

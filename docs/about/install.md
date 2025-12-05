---
title: S2 Installation
---

This document lists the platform requirements and installation
instructions for working with the S2 Geometry Library in both
C++ and Python. (The S2 Python interfaces uses SWIG to interact
with the C++ code.)

## Requirements

Running the S2 code requires the following:

*   A MacOSX or Linux platform. (Windows is not supported at this time.)
*   POSIX support (for `getrusage()`).
*   A compatible C++ compiler *supporting at least C++11*, such as
    [g++](https://gcc.gnu.org/){:target="_blank"} >= 4.7.
*   Installation of the following libraries:
    *   [OpenSSL](https://github.com/openssl/openssl){:target="_blank"} (for its
        bignum library)
    *   [gflags command line flags](https://github.com/gflags/gflags)
        {:target="_blank"} (optional, disabled by default)
    *   [glog logging module](https://github.com/google/glog) {:target="_blank"}
        (optional, disabled by default)
    *   [googletest testing framework](https://github.com/google/googletest)
        {:target="_blank"}
    *   [SWIG](https://github.com/swig/swig) (optional, for Python interface)
*   [Git](https://git-scm.com/) for interacting with the S2 source code
    repository, which is contained on [GitHub](http://github.com). To install
    Git, consult the [Set Up Git](https://help.github.com/articles/set-up-git/)
    guide on GitHub.

Installation instructions for the above components are [noted below](#Install).

<p class="note">
Note: this installation guide uses CMake as the official build system for
S2, which is supported on most major platforms and compilers. The S2
source code assumes you are using CMake and contains
<code>CMakeList.txt</code> files for that purpose.
</p>

Although you are free to use your own build system, most of the documentation
within this guide will assume you are using
[CMake](https://cmake.org/){:target="_blank"}.

## Installation Instructions {#Install}

Note: thorough testing has only been done on Ubuntu 14.04.3
and MacOSX 10.12.

### Setting Up Your Development Environment (Linux)

Note: we recommend you use a Linux package manager such as
`apt-get`. Alternatively, you may build the dependent
libraries from source.

1\. Install the additional libraries S2 code requires:

<pre>
<b>$ sudo apt-get install libgflags-dev libgoogle-glog-dev libgtest-dev libssl-dev</b>
Reading package lists... Done
Building dependency tree
Reading state information... Done
...
After this operation, 3,090 kB of additional disk space will be used.
Do you want to continue? [Y/n] Y
...
Setting up libgtest-dev (1.6.0-1ubuntu6) ...
Processing triggers for libc-bin (2.19-0ubuntu6.13) ...
<b>$</b>
</pre>

2\. Install CMake

<pre>
<b>$ sudo apt-get install cmake</b>
Reading package lists... Done
Building dependency tree
Reading state information... Done
...
After this operation, 16.6 MB of additional disk space will be used.
Do you want to continue? [Y/n] Y
...
Setting up cmake (2.8.12.2-0ubuntu3) ...
<b>$</b>
</pre>

### Setting Up Your Development Environment (MacOSX)

Note: we recommend you use a MacOSX package manager such as
[MacPorts](https://macports.org) or [Homebrew](https://brew.sh).
Alternatively, you may build the dependent libraries from source.

1\. Install the additional libraries S2 code requires:

<pre>
# MacPorts
<b>$ sudo port install gflags google-glog openssl</b>

# Homebrew
<b>$ sudo brew install gflags glog openssl</b>
</pre>

2\. Download Googletest
[release 1.8.0](https://github.com/google/googletest/releases/tag/release-1.8.0){:target="_blank"}
and unpack it in a directory of your choosing. (Take note of this
directory as you will need to point CMake to it.)

3\. Install CMake

<pre>
# Note: XCode requires command-line tools, which can be installed with:
<b>$ xcode-select --install</b>

# MacPorts
<b>$ sudo port install cmake</b>

# Homebrew
<b>$ sudo brew install cmake</b>
</pre>

### Getting the S2 Code

The [Reference implementation of S2](https://github.com/google/s2geometry)
is written in C++. Once you have CMake and Git installed, you can obtain
the S2 code from its repository on GitHub:

<pre>
# Change to the directory where you want to create the code repository
<b>$ cd ~
$ mkdir Source; cd Source</b>
</pre>

Clone the S2 code into your development directory:

<pre>
<b>$ git clone https://github.com/google/s2geometry.git</b>
Cloning into 's2geometry'...
remote: Counting objects: 5672, done.
remote: Compressing objects: 100% (52/52), done.
remote: Total 5672 (delta 21), reused 36 (delta 17), pack-reused 5601
Receiving objects: 100% (5672/5672), 7.15 MiB | 27.21 MiB/s, done.
Resolving deltas: 100% (4503/4503), done.
<b>$</b>
</pre>

Git will create the repository within a directory named `s2geometry`.
Navigate into this directory.

<p class="note">
Alternatively, you can download an archive of the S2 library as a
<a href="https://github.com/google/s2geometry/archive/master.zip">ZIP file</a>
and unzip it within your development directory.

<pre>
<b>$ unzip <i>path_to_ZIP_file</i>/s2geometry-master.zip</b>
</pre>
</p>

### Compiling S2

Once you have installed S2's requirements and the S2 library, you
are ready to build S2.

1\. Configure the make files using CMake:

<pre>
<b>$ cd s2geometry  # or geometry-master if you installed from the ZIP file
$ mkdir build
$ cd build</b>
# See Notes below.
<b>$ cmake -DWITH_GFLAGS=ON -WITH_GTEST=ON \
-DGTEST_ROOT=<i>googletest_root_dir</i> \
-DOPENSSL_INCLUDE_DIR=<i>openssl_include_dir</i> ..</b>
-- Found OpenSSL: /usr/lib/libcrypto.dylib (found version "1.0.2m")
-- Looking for pthread.h
-- Looking for pthread.h - found
-- Looking for pthread_create
-- Looking for pthread_create - found
-- Found Threads: TRUE
-- Could NOT find SWIG (missing: SWIG_EXECUTABLE SWIG_DIR)
-- Found PythonInterp: /Library/Frameworks/Python.framework/Versions/2.7/bin/python (found version "2.7.11")
-- Found PythonLibs: /Library/Frameworks/Python.framework/Versions/2.7/lib/libpython2.7.dylib (found version "2.7.11")
GTEST_ROOT: /Users/shreck/SOURCE/googletest
-- Configuring done
-- Generating done
-- Build files have been written to: /Users/shreck/SOURCE/s2geometry-master/build
<b>$</b>
</pre>

Notes:

*   `-DWITH_GFLAGS` and `-DWITH_GTEST` may be omitted to avoid depending
    on those packages.
*   `-DGTEST_ROOT` and `-DOPENSSL_INCLUDE_DIR` must be absolute paths.
*   `-DGTEST_ROOT` is usually `/usr/src/gtest` on Linux systems; on
    MacOSX use the directory where you installed GoogleTest.
*   MacOSX/Homebrew users may need to set `-DOPENSSL_INCLUDE_DIR` to
    enable CMake to find your openssl `include` directory.
*   You can omit the `DGTEST_ROOT` flag to skip tests.

2\. Make the S2 binary (and associated tests if `GTEST_ROOT` was
    specified in CMake:

<pre>
<b>$ make</b>
Scanning dependencies of target gtest
...
Scanning dependencies of target s2
...
[ 35%] Built target s2
...
[100%] Built target point_index
<b>$</b>
</pre>

3\. Run the S2 tests (if `GTEST_ROOT` was specified in CMake):

<pre>
<b>$ make test</b>
Running tests...
Test project /.../s2geometry-master/build
      Start  1: id_set_lexicon_test
...
83/83 Test #83: value_lexicon_test .............................   Passed    0.01 sec
100% tests passed, 0 tests failed out of 83
Total Test time (real) =  78.67 sec
<b>$</b>
</pre>


4\. Install the built S2 library:

<pre>
<b>$ make install</b>
[  1%] Built target gtest
...
[100%] Built target point_index
Install the project...
-- Install configuration: ""
...
<b>$</b>
</pre>

Congratulations! You've successfully installed S2, built the library, and run
all of the tests! You're now ready to move on to the
[C++ Quickstart](/devguide/cpp/quickstart).

cmake_minimum_required(VERSION 2.8)
project(s2-geometry)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
include(CheckCXXCompilerFlag)
include(FindPackageHandleStandardArgs)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/third_party/cmake")

find_package(GFlags REQUIRED)
find_package(Glog REQUIRED)
find_package(OpenSSL REQUIRED)
find_package(Threads REQUIRED)

add_definitions(-Wno-deprecated-declarations -std=c++11)
# Controls HASH_NAMESPACE in port.h.
# Use std::hash instead of __gnu_cxx::hash.
add_definitions(-DUSE_STD_HASH)

include_directories(${GFLAGS_INCLUDE_DIRS} ${GLOG_INCLUDE_DIRS})
include_directories(src)

# Explicitly list S2 files until we move test files to test dir.
# TODO: Move test files to another dir.
file(GLOB S2FilesGlob
     src/base/*.cc
     src/r2rect.cc
     src/s1angle.cc
     src/s1chordangle.cc
     src/s1interval.cc
     src/s2.cc
     src/s2cap.cc
     src/s2cell.cc
     src/s2cellid.cc
     src/s2cellunion.cc
     src/s2closestedgequery.cc
     src/s2edgequery.cc
     src/s2edgeutil.cc
     src/s2error.cc
     src/s2latlng.cc
     src/s2latlngrect.cc
     src/s2loop.cc
     src/s2paddedcell.cc
     src/s2pointcompression.cc
     src/s2pointregion.cc
     src/s2polygon.cc
     src/s2polygonbuilder.cc
     src/s2polyline.cc
     src/s2r2rect.cc
     src/s2region.cc
     src/s2regioncoverer.cc
     src/s2regionintersection.cc
     src/s2regionunion.cc
     src/s2shapeindex.cc
     src/s2shapeutil.cc
     src/s2textformat.cc
     src/strings/*.cc
     src/util/bits/*.cc
     src/util/coding/*.cc
     src/util/hash/*.cc
     src/util/math/exactfloat/*.cc
     src/util/math/*.cc)

add_library(s2geometry STATIC ${S2FilesGlob})

file(GLOB S2Headers src/*.h)

# We don't need to install all headers, only those
# transitively included by s2 headers we are exporting.
install(FILES ${S2Headers}
        DESTINATION include/s2geometry)
install(FILES src/base/atomicops.h
              src/base/casts.h
              src/base/int128.h
              src/base/integral_types.h
              src/base/macros.h
              src/base/mutex.h
              src/base/port.h
              src/base/port_hash.h
              src/base/spinlock.h
              src/base/stl_decl_msvc.h
              src/base/template_util.h
              src/base/type_traits.h
        DESTINATION include/s2geometry/base)
install(FILES src/util/bits/bits.h
              src/util/bits/bits-internal-unknown.h
              src/util/bits/bits-internal-windows.h
        DESTINATION include/s2geometry/util/bits)
install(FILES src/util/gtl/inlined_vector.h
              src/util/gtl/manual_constructor.h
        DESTINATION include/s2geometry/util/gtl)
install(FILES src/util/hash/builtin_type_hash.h
              src/util/hash/city.h
              src/util/hash/hash.h
              src/util/hash/hash128to64.h
              src/util/hash/jenkins.h
              src/util/hash/jenkins_lookup2.h
              src/util/hash/string_hash.h
        DESTINATION include/s2geometry/util/hash)
install(FILES src/util/math/mathlimits.h
              src/util/math/mathutil.h
              src/util/math/matrix3x3.h
              src/util/math/vector2.h
              src/util/math/vector3.h
              src/util/math/vector3_hash.h
              src/util/math/vector4.h
        DESTINATION include/s2geometry/util/math)
install(FILES src/util/units/length-units.h
              src/util/units/physical-units.h
        DESTINATION include/s2geometry/util/units)
install(TARGETS s2geometry DESTINATION lib)

message("GTEST_ROOT: ${GTEST_ROOT}")
if (GTEST_ROOT)
  add_library(gtest_main STATIC ${GTEST_ROOT}/src/gtest-all.cc ${GTEST_ROOT}/src/gtest_main.cc)
  add_library(s2testing STATIC src/s2testing.cc)

  file(GLOB S2TestFilesGlob src/*_test.cc)

  enable_testing()
  include_directories(${GTEST_ROOT}/include)
  include_directories(${GTEST_ROOT})

  foreach (test_cc ${S2TestFilesGlob})
    get_filename_component(test ${test_cc} NAME_WE)
    add_executable(${test} ${test_cc})
    target_link_libraries(
        ${test}
        s2testing s2geometry gtest_main
        ${GFLAGS_LIBRARIES} ${GLOG_LIBRARIES} ${OPENSSL_LIBRARIES}
       ${CMAKE_THREAD_LIBS_INIT})
    add_test(${test} ${test})
  endforeach()
endif()


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
find_package(SWIG)
find_package(PythonInterp)
find_package(PythonLibs)

add_definitions(-Wno-deprecated-declarations -std=c++11)

# If OpenSSL is installed in a non-standard location, configure with
# something like:
# OPENSSL_ROOT_DIR=/usr/local/opt/openssl cmake ..
include_directories(
    ${GFLAGS_INCLUDE_DIRS} ${GLOG_INCLUDE_DIRS} ${OPENSSL_INCLUDE_DIR}
    ${PYTHON_INCLUDE_DIRS})
include_directories(src)

# Explicitly list S2 files until we move test files to test dir.
# TODO: Move test files to another dir.
file(GLOB S2FilesGlob
     src/s2/base/*.cc
     src/s2/id_set_lexicon.cc
     src/s2/r2rect.cc
     src/s2/s1angle.cc
     src/s2/s1chordangle.cc
     src/s2/s1interval.cc
     src/s2/s2boundary_operation.cc
     src/s2/s2builder.cc
     src/s2/s2builder_graph.cc
     src/s2/s2builderutil_layers.cc
     src/s2/s2builderutil_snap_functions.cc
     src/s2/s2cap.cc
     src/s2/s2cell.cc
     src/s2/s2cellid.cc
     src/s2/s2cellunion.cc
     src/s2/s2centroids.cc
     src/s2/s2closestedgequery.cc
     src/s2/s2convexhullquery.cc
     src/s2/s2coords.cc
     src/s2/s2crossingedgequery.cc
     src/s2/s2debug.cc
     src/s2/s2earth.cc
     src/s2/s2edgeutil.cc
     src/s2/s2error.cc
     src/s2/s2latlng.cc
     src/s2/s2latlngrect.cc
     src/s2/s2loop.cc
     src/s2/s2measures.cc
     src/s2/s2metrics.cc
     src/s2/s2paddedcell.cc
     src/s2/s2pointcompression.cc
     src/s2/s2pointregion.cc
     src/s2/s2pointutil.cc
     src/s2/s2polygon.cc
     src/s2/s2polyline.cc
     src/s2/s2polyline_simplifier.cc
     src/s2/s2predicates.cc
     src/s2/s2r2rect.cc
     src/s2/s2region.cc
     src/s2/s2regioncoverer.cc
     src/s2/s2regionintersection.cc
     src/s2/s2regionunion.cc
     src/s2/s2shapeindex.cc
     src/s2/s2shapeutil.cc
     src/s2/s2textformat.cc
     src/s2/sequence_lexicon.cc
     src/s2/strings/*.cc
     src/s2/third_party/absl/port/internal/*.cc
     src/s2/util/bits/*.cc
     src/s2/util/coding/*.cc
     src/s2/util/hash/*.cc
     src/s2/util/math/exactfloat/*.cc
     src/s2/util/math/*.cc
     src/s2/value_lexicon.cc)

add_library(s2 SHARED ${S2FilesGlob})
target_link_libraries(
    s2
    ${GFLAGS_LIBRARIES} ${GLOG_LIBRARIES} ${OPENSSL_LIBRARIES}
    ${CMAKE_THREAD_LIBS_INIT})

file(GLOB S2Headers src/s2/*.h)

# We don't need to install all headers, only those
# transitively included by s2 headers we are exporting.
install(FILES ${S2Headers}
        DESTINATION include/s2)
install(FILES src/s2/base/atomicops.h
              src/s2/base/casts.h
              src/s2/base/mutex.h
              src/s2/base/spinlock.h
              src/s2/base/template_util.h
              src/s2/base/type_traits.h
        DESTINATION include/s2/base)
install(FILES src/s2/third_party/absl/base/casts.h
              src/s2/third_party/absl/base/config.h
              src/s2/third_party/absl/base/integral_types.h
              src/s2/third_party/absl/base/macros.h
              src/s2/third_party/absl/base/port.h
        DESTINATION include/s2/third_party/absl/base)
install(FILES src/s2/third_party/absl/base/internal/port_hash.inc
        DESTINATION include/s2/third_party/absl/base/internal)
install(FILES src/s2/third_party/absl/container/inlined_vector.h
        DESTINATION include/s2/third_party/absl/container)
install(FILES src/s2/third_party/absl/meta/type_traits.h
        DESTINATION include/s2/third_party/absl/meta)
install(FILES src/s2/third_party/absl/numeric/int128.h
        DESTINATION include/s2/third_party/absl/numeric)
install(FILES src/s2/util/bits/bits.h
        DESTINATION include/s2/util/bits)
install(FILES src/s2/util/btree/btree.h
              src/s2/util/btree/btree_container.h
              src/s2/util/btree/btree_map.h
              src/s2/util/btree/btree_set.h
        DESTINATION include/s2/util/btree)
install(FILES src/s2/util/coding/coder.h
              src/s2/util/coding/varint.h
        DESTINATION include/s2/util/coding)
install(FILES src/s2/util/endian/endian.h
        DESTINATION include/s2/util/endian)
install(FILES src/s2/util/gtl/compact_array.h
              src/s2/util/gtl/container_logging.h
              src/s2/util/gtl/dense_hash_set.h
              src/s2/util/gtl/densehashtable.h
              src/s2/util/gtl/hashtable_common.h
              src/s2/util/gtl/libc_allocator_with_realloc.h
              src/s2/util/gtl/manual_constructor.h
        DESTINATION include/s2/util/gtl)
install(FILES src/s2/util/hash/mix.h
        DESTINATION include/s2/util/hash)
install(FILES src/s2/util/math/mathutil.h
              src/s2/util/math/matrix3x3.h
              src/s2/util/math/vector.h
              src/s2/util/math/vector3_hash.h
        DESTINATION include/s2/util/math)
install(FILES src/s2/util/units/length-units.h
              src/s2/util/units/physical-units.h
        DESTINATION include/s2/util/units)
install(TARGETS s2 DESTINATION lib)

message("GTEST_ROOT: ${GTEST_ROOT}")
if (GTEST_ROOT)
  add_library(gtest_main STATIC ${GTEST_ROOT}/src/gtest-all.cc ${GTEST_ROOT}/src/gtest_main.cc)
  add_library(s2testing STATIC src/s2/s2testing.cc)

  file(GLOB S2TestFilesGlob src/s2/*_test.cc)

  enable_testing()
  include_directories(${GTEST_ROOT}/include)
  include_directories(${GTEST_ROOT})

  foreach (test_cc ${S2TestFilesGlob})
    get_filename_component(test ${test_cc} NAME_WE)
    add_executable(${test} ${test_cc})
    target_link_libraries(
        ${test}
        s2testing s2 gtest_main)
    add_test(${test} ${test})
  endforeach()
endif()

if (${SWIG_FOUND} AND ${PYTHONLIBS_FOUND})
  include(${SWIG_USE_FILE})
  include_directories(${PYTHON_INCLUDE_PATH})
  set(CMAKE_SWIG_FLAGS "")
  set_property(SOURCE src/s2/s2.i PROPERTY SWIG_FLAGS "-module" "pywraps2")
  set_property(SOURCE src/s2/s2.i PROPERTY CPLUSPLUS ON)
  swig_add_module(pywraps2 python src/s2/s2.i)
  swig_link_libraries(pywraps2 ${PYTHON_LIBRARIES} s2)
  enable_testing()
  add_test(NAME pywraps2_test COMMAND
           ${PYTHON_EXECUTABLE} "${PROJECT_SOURCE_DIR}/src/s2/pywraps2_test.py")
  set_property(TEST pywraps2_test PROPERTY ENVIRONMENT
               "PYTHONPATH=$ENV{PYTHONPATH}:${PROJECT_BINARY_DIR}")
endif()

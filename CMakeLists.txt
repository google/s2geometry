cmake_minimum_required(VERSION 3.18)

project(s2-geometry
        VERSION 0.12.0)

include(CMakeDependentOption)
include(CheckCXXCompilerFlag)
include(FeatureSummary)
include(FindPackageHandleStandardArgs)
include(GNUInstallDirs)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# Avoid "Policy CMP0042 is not set" warning on macOS.
if (APPLE)
    set(CMAKE_MACOSX_RPATH TRUE)
endif()

# s2geometry needs to use the same C++ standard that absl used to avoid
# undefined symbol errors since ABSL_HAVE_STD_STRING_VIEW etc will
# end up defined differently.  There is probably a better way to achieve
# this than assuming what absl used.  We require at least C++17.
# Using CACHE allows the user to override the default.
set(CMAKE_CXX_STANDARD 17 CACHE STRING "The C++ standard to build with")
set(CMAKE_CXX_STANDARD_REQUIRED ON)
# No compiler-specific extensions, i.e. -std=c++17, not -std=gnu++17.
set(CMAKE_CXX_EXTENSIONS OFF)

option(BUILD_SHARED_LIBS "Build shared libraries instead of static." ON)
add_feature_info(SHARED_LIBS BUILD_SHARED_LIBS
                 "builds shared libraries instead of static.")

option(BUILD_EXAMPLES "Build s2 documentation examples." ON)
option(BUILD_TESTS "Build s2 unittests." ON)

option(WITH_PYTHON "Add python interface" OFF)
add_feature_info(PYTHON WITH_PYTHON "provides python interface to S2")

option(S2_USE_SYSTEM_INCLUDES
    "Silence warnings in s2 headers by marking them as SYSTEM includes." OFF)
add_feature_info(S2_USE_SYSTEM_INCLUDES S2_USE_SYSTEM_INCLUDES
    "Silence warnings in s2 headers by marking them as SYSTEM includes.")

# when s2geometry is included as subproject (i.e. using add_subdirectory(s2geometry))
# in the source tree of a project that uses it, install rules are disabled.
if(NOT CMAKE_SOURCE_DIR STREQUAL PROJECT_SOURCE_DIR)
  option(S2_ENABLE_INSTALL "Enable install rule" OFF)
else()
  option(S2_ENABLE_INSTALL "Enable install rule" ON)
endif()

feature_summary(WHAT ALL)

# Don't output anything for LOG(INFO).
add_definitions(-DABSL_MIN_LOG_LEVEL=1)

# If we already have abseil we don't need to find it, for an example:
# add_subdirectory(absl-submbodule)
# add_subdirectory(s2-submodule)
if (NOT TARGET absl::base)
    find_package(absl REQUIRED)
endif()
# pthreads isn't used directly, but this is still required for std::thread.
find_package(Threads REQUIRED)

if (NOT TARGET absl::vlog_is_on)
    message(FATAL_ERROR "Could not find absl vlog module. Are you using an older version?")
endif()

if (WITH_PYTHON)
    # Should be easy to make it work with swig3, but some args to %pythonprepend
    # seem to be different and were changed.
    find_package(SWIG 4.0)
    # Use Python3_ROOT_DIR to help find python3, if the correct location is not
    # being found by default.
    find_package(Python3 COMPONENTS Interpreter Development.Module)
endif()

if (MSVC)
    # Use unsigned characters
    add_definitions(-J)
    # Make sure cmath header defines things like M_PI
    add_definitions(-D_USE_MATH_DEFINES)
    # Make sure Windows doesn't define min/max macros that interfere with STL
    add_definitions(-DNOMINMAX)
else()
    # Avoid megabytes of warnings like:
    # util/math/vector.h:178:16: warning: optimization attribute on
    # ‘double sqrt(double)’ follows definition but the attribute doesn’t
    # match [-Wattributes]
    add_definitions(-Wno-attributes)
    add_definitions(-Wno-deprecated-declarations)
    # `base/port.h` uses sized deallocation, which should be enabled by
    # default for C++14 and later.  There appears to be a bug with clang
    # < 19 causing this not to be enabled.
    # https://github.com/google/s2geometry/issues/411#issuecomment-2726949607
    # This can be removed when clang19 is the minimum supported version.
    add_compile_options(-fsized-deallocation)
endif()


if (WITH_PYTHON)
    include_directories(${Python3_INCLUDE_DIRS})
endif()

include_directories(src)

add_library(s2
            src/s2/base/malloc_extension.cc
            src/s2/encoded_s2cell_id_vector.cc
            src/s2/encoded_s2point_vector.cc
            src/s2/encoded_s2shape_index.cc
            src/s2/encoded_string_vector.cc
            src/s2/id_set_lexicon.cc
            src/s2/internal/s2incident_edge_tracker.cc
            src/s2/internal/s2index_cell_data.cc
            src/s2/mutable_s2shape_index.cc
            src/s2/r2rect.cc
            src/s2/s1angle.cc
            src/s2/s1chord_angle.cc
            src/s2/s1interval.cc
            src/s2/s2boolean_operation.cc
            src/s2/s2buffer_operation.cc
            src/s2/s2builder.cc
            src/s2/s2builder_graph.cc
            src/s2/s2builderutil_closed_set_normalizer.cc
            src/s2/s2builderutil_find_polygon_degeneracies.cc
            src/s2/s2builderutil_get_snapped_winding_delta.cc
            src/s2/s2builderutil_lax_polygon_layer.cc
            src/s2/s2builderutil_lax_polyline_layer.cc
            src/s2/s2builderutil_s2point_vector_layer.cc
            src/s2/s2builderutil_s2polygon_layer.cc
            src/s2/s2builderutil_s2polyline_layer.cc
            src/s2/s2builderutil_s2polyline_vector_layer.cc
            src/s2/s2builderutil_snap_functions.cc
            src/s2/s2cap.cc
            src/s2/s2cell.cc
            src/s2/s2cell_id.cc
            src/s2/s2cell_index.cc
            src/s2/s2cell_union.cc
            src/s2/s2centroids.cc
            src/s2/s2chain_interpolation_query.cc
            src/s2/s2closest_cell_query.cc
            src/s2/s2closest_edge_query.cc
            src/s2/s2closest_point_query.cc
            src/s2/s2contains_vertex_query.cc
            src/s2/s2convex_hull_query.cc
            src/s2/s2coords.cc
            src/s2/s2crossing_edge_query.cc
            src/s2/s2debug.cc
            src/s2/s2density_tree.cc
            src/s2/s2earth.cc
            src/s2/s2edge_clipping.cc
            src/s2/s2edge_crosser.cc
            src/s2/s2edge_crossings.cc
            src/s2/s2edge_distances.cc
            src/s2/s2edge_tessellator.cc
            src/s2/s2error.cc
            src/s2/s2fractal.cc
            src/s2/s2furthest_edge_query.cc
            src/s2/s2hausdorff_distance_query.cc
            src/s2/s2latlng.cc
            src/s2/s2latlng_rect.cc
            src/s2/s2latlng_rect_bounder.cc
            src/s2/s2lax_loop_shape.cc
            src/s2/s2lax_polygon_shape.cc
            src/s2/s2lax_polyline_shape.cc
            src/s2/s2loop.cc
            src/s2/s2loop_measures.cc
            src/s2/s2max_distance_targets.cc
            src/s2/s2measures.cc
            src/s2/s2memory_tracker.cc
            src/s2/s2metrics.cc
            src/s2/s2min_distance_targets.cc
            src/s2/s2padded_cell.cc
            src/s2/s2point_compression.cc
            src/s2/s2point_region.cc
            src/s2/s2pointutil.cc
            src/s2/s2polygon.cc
            src/s2/s2polyline.cc
            src/s2/s2polyline_alignment.cc
            src/s2/s2polyline_measures.cc
            src/s2/s2polyline_simplifier.cc
            src/s2/s2predicates.cc
            src/s2/s2projections.cc
            src/s2/s2r2rect.cc
            src/s2/s2random.cc
            src/s2/s2region_coverer.cc
            src/s2/s2region_intersection.cc
            src/s2/s2region_sharder.cc
            src/s2/s2region_term_indexer.cc
            src/s2/s2region_union.cc
            src/s2/s2shape_index.cc
            src/s2/s2shape_index_buffered_region.cc
            src/s2/s2shape_index_measures.cc
            src/s2/s2shape_measures.cc
            src/s2/s2shape_nesting_query.cc
            src/s2/s2shapeutil_build_polygon_boundaries.cc
            src/s2/s2shapeutil_coding.cc
            src/s2/s2shapeutil_contains_brute_force.cc
            src/s2/s2shapeutil_conversion.cc
            src/s2/s2shapeutil_count_vertices.cc
            src/s2/s2shapeutil_edge_iterator.cc
            src/s2/s2shapeutil_edge_wrap.cc
            src/s2/s2shapeutil_get_reference_point.cc
            src/s2/s2shapeutil_visit_crossing_edge_pairs.cc
            src/s2/s2text_format.cc
            src/s2/s2wedge_relations.cc
            src/s2/s2winding_operation.cc
            src/s2/util/bits/bit-interleave.cc
            src/s2/util/coding/coder.cc
            src/s2/util/coding/varint.cc
            src/s2/util/math/exactfloat/bignum.cc
            src/s2/util/math/exactfloat/exactfloat.cc
            src/s2/util/math/mathutil.cc
            src/s2/util/units/length-units.cc)

if (GOOGLETEST_ROOT)
  add_library(s2testing STATIC
              src/s2/gmock_matchers.cc
              src/s2/s2builderutil_testing.cc
              src/s2/s2shapeutil_testing.cc
              src/s2/s2testing.cc
              src/s2/thread_testing.cc)
endif()


target_link_libraries(
    s2
    absl::absl_vlog_is_on
    absl::base
    absl::btree
    absl::check
    absl::config
    absl::core_headers
    absl::dynamic_annotations
    absl::endian
    absl::fixed_array
    absl::flags
    absl::flat_hash_map
    absl::flat_hash_set
    absl::hash
    absl::inlined_vector
    absl::int128
    absl::log
    absl::log_severity
    absl::memory
    absl::span
    absl::status
    absl::str_format
    absl::strings
    absl::type_traits
    absl::utility
    absl::vlog_is_on
    ${CMAKE_THREAD_LIBS_INIT})

if (GOOGLETEST_ROOT)
  target_link_libraries(
      s2testing
      absl::check
      absl::flags
      absl::log
      absl::memory
      absl::strings)
endif()

# Allow other CMake projects to use this one with:
# list(APPEND CMAKE_MODULE_PATH "<path_to_s2geometry_dir>/third_party/cmake")
# add_subdirectory(<path_to_s2geometry_dir> s2geometry)
# target_link_libraries(<target_name> s2)
if (S2_USE_SYSTEM_INCLUDES)
    target_include_directories(s2 SYSTEM PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/src)
else ()
    target_include_directories(s2 PUBLIC
      $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src>
      $<INSTALL_INTERFACE:include>)
endif ()

# Add version information to the target
set_target_properties(s2 PROPERTIES
    SOVERSION ${PROJECT_VERSION_MAJOR}
    VERSION ${PROJECT_VERSION})

if(S2_ENABLE_INSTALL)
  # We don't need to install all headers, only those
  # transitively included by s2 headers we are exporting.
  install(FILES src/s2/_fp_contract_off.h
                src/s2/encoded_s2cell_id_vector.h
                src/s2/encoded_s2point_vector.h
                src/s2/encoded_s2shape_index.h
                src/s2/encoded_string_vector.h
                src/s2/encoded_uint_vector.h
                src/s2/gmock_matchers.h
                src/s2/id_set_lexicon.h
                src/s2/mutable_s2shape_index.h
                src/s2/r1interval.h
                src/s2/r2.h
                src/s2/r2rect.h
                src/s2/s1angle.h
                src/s2/s1chord_angle.h
                src/s2/s1interval.h
                src/s2/s2boolean_operation.h
                src/s2/s2buffer_operation.h
                src/s2/s2builder.h
                src/s2/s2builder_graph.h
                src/s2/s2builder_layer.h
                src/s2/s2builderutil_closed_set_normalizer.h
                src/s2/s2builderutil_find_polygon_degeneracies.h
                src/s2/s2builderutil_get_snapped_winding_delta.h
                src/s2/s2builderutil_graph_shape.h
                src/s2/s2builderutil_lax_polygon_layer.h
                src/s2/s2builderutil_lax_polyline_layer.h
                src/s2/s2builderutil_s2point_vector_layer.h
                src/s2/s2builderutil_s2polygon_layer.h
                src/s2/s2builderutil_s2polyline_layer.h
                src/s2/s2builderutil_s2polyline_vector_layer.h
                src/s2/s2builderutil_snap_functions.h
                src/s2/s2builderutil_testing.h
                src/s2/s2cap.h
                src/s2/s2cell.h
                src/s2/s2cell_id.h
                src/s2/s2cell_index.h
                src/s2/s2cell_iterator.h
                src/s2/s2cell_iterator_join.h
                src/s2/s2cell_range_iterator.h
                src/s2/s2cell_union.h
                src/s2/s2centroids.h
                src/s2/s2chain_interpolation_query.h
                src/s2/s2closest_cell_query.h
                src/s2/s2closest_cell_query_base.h
                src/s2/s2closest_edge_query.h
                src/s2/s2closest_edge_query_base.h
                src/s2/s2closest_point_query.h
                src/s2/s2closest_point_query_base.h
                src/s2/s2coder.h
                src/s2/s2contains_point_query.h
                src/s2/s2contains_vertex_query.h
                src/s2/s2convex_hull_query.h
                src/s2/s2coords.h
                src/s2/s2coords_internal.h
                src/s2/s2crossing_edge_query.h
                src/s2/s2debug.h
                src/s2/s2density_tree.h
                src/s2/s2density_tree_internal.h
                src/s2/s2distance_target.h
                src/s2/s2earth.h
                src/s2/s2edge_clipping.h
                src/s2/s2edge_crosser.h
                src/s2/s2edge_crossings.h
                src/s2/s2edge_crossings_internal.h
                src/s2/s2edge_distances.h
                src/s2/s2edge_tessellator.h
                src/s2/s2edge_vector_shape.h
                src/s2/s2error.h
                src/s2/s2fractal.h
                src/s2/s2furthest_edge_query.h
                src/s2/s2hausdorff_distance_query.h
                src/s2/s2latlng.h
                src/s2/s2latlng_rect.h
                src/s2/s2latlng_rect_bounder.h
                src/s2/s2lax_loop_shape.h
                src/s2/s2lax_polygon_shape.h
                src/s2/s2lax_polyline_shape.h
                src/s2/s2loop.h
                src/s2/s2loop_measures.h
                src/s2/s2max_distance_targets.h
                src/s2/s2measures.h
                src/s2/s2memory_tracker.h
                src/s2/s2metrics.h
                src/s2/s2min_distance_targets.h
                src/s2/s2padded_cell.h
                src/s2/s2point.h
                src/s2/s2point_compression.h
                src/s2/s2point_index.h
                src/s2/s2point_region.h
                src/s2/s2point_span.h
                src/s2/s2point_vector_shape.h
                src/s2/s2pointutil.h
                src/s2/s2polygon.h
                src/s2/s2polyline.h
                src/s2/s2polyline_alignment.h
                src/s2/s2polyline_measures.h
                src/s2/s2polyline_simplifier.h
                src/s2/s2predicates.h
                src/s2/s2predicates_internal.h
                src/s2/s2projections.h
                src/s2/s2r2rect.h
                src/s2/s2random.h
                src/s2/s2region.h
                src/s2/s2region_coverer.h
                src/s2/s2region_intersection.h
                src/s2/s2region_sharder.h
                src/s2/s2region_term_indexer.h
                src/s2/s2region_union.h
                src/s2/s2shape.h
                src/s2/s2shape.h
                src/s2/s2shape_index.h
                src/s2/s2shape_index_buffered_region.h
                src/s2/s2shape_index_region.h
                src/s2/s2shape_measures.h
                src/s2/s2shape_nesting_query.h
                src/s2/s2shapeutil_build_polygon_boundaries.h
                src/s2/s2shapeutil_coding.h
                src/s2/s2shapeutil_contains_brute_force.h
                src/s2/s2shapeutil_conversion.h
                src/s2/s2shapeutil_count_edges.h
                src/s2/s2shapeutil_count_vertices.h
                src/s2/s2shapeutil_edge_iterator.h
                src/s2/s2shapeutil_edge_wrap.h
                src/s2/s2shapeutil_get_reference_point.h
                src/s2/s2shapeutil_shape_edge.h
                src/s2/s2shapeutil_shape_edge_id.h
                src/s2/s2shapeutil_testing.h
                src/s2/s2shapeutil_visit_crossing_edge_pairs.h
                src/s2/s2testing.h
                src/s2/s2text_format.h
                src/s2/s2validation_query.h
                src/s2/s2wedge_relations.h
                src/s2/s2winding_operation.h
                src/s2/s2wrapped_shape.h
                src/s2/sequence_lexicon.h
                src/s2/thread_testing.h
                src/s2/value_lexicon.h
          DESTINATION include/s2)
  install(FILES src/s2/internal/s2disjoint_set.h
                src/s2/internal/s2incident_edge_tracker.h
                src/s2/internal/s2index_cell_data.h
                src/s2/internal/s2meta.h
          DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/s2/internal")
  install(FILES src/s2/base/casts.h
                src/s2/base/commandlineflags.h
                src/s2/base/commandlineflags_declare.h
                src/s2/base/log_severity.h
                src/s2/base/malloc_extension.h
                src/s2/base/port.h
                src/s2/base/spinlock.h
          DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/s2/base")
  install(FILES src/s2/testing/gtest_prod.h
          DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/s2/testing")
  install(FILES src/s2/util/bitmap/bitmap.h
          DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/s2/util/bitmap")
  install(FILES src/s2/util/coding/coder.h
                src/s2/util/coding/varint.h
          DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/s2/util/coding")
  install(FILES src/s2/util/endian/endian.h
          DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/s2/util/endian")
  install(FILES src/s2/util/gtl/compact_array.h
                src/s2/util/gtl/container_logging.h
                src/s2/util/gtl/dense_hash_set.h
                src/s2/util/gtl/densehashtable.h
                src/s2/util/gtl/hashtable_common.h
                src/s2/util/gtl/requires.h
          DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/s2/util/gtl")
  install(FILES src/s2/util/hash/mix.h
          DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/s2/util/hash")
  install(FILES src/s2/util/math/mathutil.h
                src/s2/util/math/matrix3x3.h
                src/s2/util/math/vector.h
          DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/s2/util/math")
  install(FILES src/s2/util/math/exactfloat/bignum.h
                src/s2/util/math/exactfloat/exactfloat.h
          DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/s2/util/math/exactfloat")
  install(FILES src/s2/util/units/length-units.h
                src/s2/util/units/physical-units.h
          DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/s2/util/units")

  if (GOOGLETEST_ROOT)
    set(S2_TARGETS s2 s2testing)
  else()
    set(S2_TARGETS s2)
  endif()

  install(TARGETS ${S2_TARGETS}
          RUNTIME DESTINATION "${CMAKE_INSTALL_BINDIR}"
          ARCHIVE DESTINATION "${CMAKE_INSTALL_LIBDIR}"
          LIBRARY DESTINATION "${CMAKE_INSTALL_LIBDIR}")

  # Create an export type "s2Targets" detailing the targetable artifacts created
  # by this project.
  install(TARGETS s2 EXPORT s2Targets)
  # Install the export targets as a CMake config file in the share/s2 folder so
  # that they can referenced  by downstream projects as `s2::s2` after a
  # successful `find_package` call.
  install(EXPORT s2Targets
          NAMESPACE s2::
          FILE s2Targets.cmake
          DESTINATION share/s2/)
endif() # S2_ENABLE_INSTALL

if (BUILD_TESTS)
  if (NOT GOOGLETEST_ROOT)
    message(FATAL_ERROR "BUILD_TESTS requires GOOGLETEST_ROOT")
  endif()
  message("GOOGLETEST_ROOT: ${GOOGLETEST_ROOT}")

  add_subdirectory(${GOOGLETEST_ROOT}/googlemock build_gmock)
  include_directories(${GOOGLETEST_ROOT}/googlemock/include)
  include_directories(${GOOGLETEST_ROOT}/googletest/include)

  set(S2TestFiles
      src/s2/encoded_s2cell_id_vector_test.cc
      src/s2/encoded_s2point_vector_test.cc
      src/s2/encoded_s2shape_index_test.cc
      src/s2/encoded_string_vector_test.cc
      src/s2/encoded_uint_vector_test.cc
      src/s2/gmock_matchers_test.cc
      src/s2/id_set_lexicon_test.cc
      src/s2/internal/s2disjoint_set_test.cc
      src/s2/internal/s2index_cell_data_test.cc
      src/s2/mutable_s2shape_index_test.cc
      src/s2/r1interval_test.cc
      src/s2/r2rect_test.cc
      src/s2/s1angle_test.cc
      src/s2/s1chord_angle_test.cc
      src/s2/s1interval_test.cc
      src/s2/s2boolean_operation_test.cc
      src/s2/s2buffer_operation_test.cc
      src/s2/s2builder_graph_test.cc
      src/s2/s2builder_test.cc
      src/s2/s2builderutil_closed_set_normalizer_test.cc
      src/s2/s2builderutil_find_polygon_degeneracies_test.cc
      src/s2/s2builderutil_get_snapped_winding_delta_test.cc
      src/s2/s2builderutil_lax_polygon_layer_test.cc
      src/s2/s2builderutil_lax_polyline_layer_test.cc
      src/s2/s2builderutil_s2point_vector_layer_test.cc
      src/s2/s2builderutil_s2polygon_layer_test.cc
      src/s2/s2builderutil_s2polyline_layer_test.cc
      src/s2/s2builderutil_s2polyline_vector_layer_test.cc
      src/s2/s2builderutil_snap_functions_test.cc
      src/s2/s2builderutil_testing_test.cc
      src/s2/s2cap_test.cc
      src/s2/s2cell_id_test.cc
      src/s2/s2cell_index_test.cc
      src/s2/s2cell_iterator_join_test.cc
      src/s2/s2cell_iterator_testing_test.cc
      src/s2/s2cell_range_iterator_test.cc
      src/s2/s2cell_test.cc
      src/s2/s2cell_union_test.cc
      src/s2/s2centroids_test.cc
      src/s2/s2chain_interpolation_query_test.cc
      src/s2/s2closest_cell_query_base_test.cc
      src/s2/s2closest_cell_query_test.cc
      src/s2/s2closest_edge_query_base_test.cc
      src/s2/s2closest_edge_query_test.cc
      src/s2/s2closest_point_query_base_test.cc
      src/s2/s2closest_point_query_test.cc
      src/s2/s2contains_point_query_test.cc
      src/s2/s2contains_vertex_query_test.cc
      src/s2/s2convex_hull_query_test.cc
      src/s2/s2coords_test.cc
      src/s2/s2crossing_edge_query_test.cc
      src/s2/s2density_tree_test.cc
      src/s2/s2earth_test.cc
      src/s2/s2edge_clipping_test.cc
      src/s2/s2edge_crosser_test.cc
      src/s2/s2edge_crossings_test.cc
      src/s2/s2edge_distances_test.cc
      src/s2/s2edge_tessellator_test.cc
      src/s2/s2edge_vector_shape_test.cc
      src/s2/s2error_test.cc
      src/s2/s2fractal_test.cc
      src/s2/s2furthest_edge_query_test.cc
      src/s2/s2hausdorff_distance_query_test.cc
      src/s2/s2latlng_rect_bounder_test.cc
      src/s2/s2latlng_rect_test.cc
      src/s2/s2latlng_test.cc
      src/s2/s2lax_loop_shape_test.cc
      src/s2/s2lax_polygon_shape_test.cc
      src/s2/s2lax_polyline_shape_test.cc
      src/s2/s2loop_measures_test.cc
      src/s2/s2loop_test.cc
      src/s2/s2max_distance_targets_test.cc
      src/s2/s2measures_test.cc
      src/s2/s2memory_tracker_test.cc
      src/s2/s2metrics_test.cc
      src/s2/s2min_distance_targets_test.cc
      src/s2/s2padded_cell_test.cc
      src/s2/s2point_compression_test.cc
      src/s2/s2point_index_test.cc
      src/s2/s2point_region_test.cc
      src/s2/s2point_test.cc
      src/s2/s2point_vector_shape_test.cc
      src/s2/s2pointutil_test.cc
      src/s2/s2polygon_test.cc
      src/s2/s2polyline_alignment_test.cc
      src/s2/s2polyline_measures_test.cc
      src/s2/s2polyline_simplifier_test.cc
      src/s2/s2polyline_test.cc
      src/s2/s2predicates_test.cc
      src/s2/s2projections_test.cc
      src/s2/s2r2rect_test.cc
      src/s2/s2random_test.cc
      src/s2/s2region_coverer_test.cc
      src/s2/s2region_sharder_test.cc
      src/s2/s2region_term_indexer_test.cc
      src/s2/s2region_test.cc
      src/s2/s2region_union_test.cc
      src/s2/s2shape_index_buffered_region_test.cc
      src/s2/s2shape_index_measures_test.cc
      src/s2/s2shape_index_region_test.cc
      src/s2/s2shape_index_test.cc
      src/s2/s2shape_measures_test.cc
      src/s2/s2shape_nesting_query_test.cc
      src/s2/s2shapeutil_build_polygon_boundaries_test.cc
      src/s2/s2shapeutil_coding_test.cc
      src/s2/s2shapeutil_contains_brute_force_test.cc
      src/s2/s2shapeutil_conversion_test.cc
      src/s2/s2shapeutil_count_edges_test.cc
      src/s2/s2shapeutil_count_vertices_test.cc
      src/s2/s2shapeutil_edge_iterator_test.cc
      src/s2/s2shapeutil_edge_wrap_test.cc
      src/s2/s2shapeutil_get_reference_point_test.cc
      src/s2/s2shapeutil_shape_edge_id_test.cc
      src/s2/s2shapeutil_visit_crossing_edge_pairs_test.cc
      src/s2/s2text_format_test.cc
      src/s2/s2validation_query_test.cc
      src/s2/s2wedge_relations_test.cc
      src/s2/s2winding_operation_test.cc
      src/s2/s2wrapped_shape_test.cc
      src/s2/sequence_lexicon_test.cc
      src/s2/util/math/exactfloat/bignum_test.cc
      src/s2/util/math/exactfloat/exactfloat_test.cc
      src/s2/value_lexicon_test.cc)

  enable_testing()

  foreach (test_cc ${S2TestFiles})
    get_filename_component(test ${test_cc} NAME_WE)
    add_executable(${test} ${test_cc})
    target_link_libraries(
        ${test}
        s2testing s2
        absl::base
        absl::btree
        absl::check
        absl::core_headers
        absl::flags
        absl::flags_reflection
        absl::log
        absl::memory
        absl::random_random
        absl::span
        absl::status
        absl::strings
        absl::synchronization
        gmock_main
        crypto)
    add_test(${test} ${test})
  endforeach()
endif()

if (BUILD_EXAMPLES AND TARGET s2testing)
  add_subdirectory("doc/examples" examples)
endif()

if (SWIG_FOUND AND Python3_FOUND)
  add_subdirectory("src/python" python)
endif()

include(CMakePackageConfigHelpers)

# Generate the config file that includes the exports.
configure_package_config_file(${CMAKE_CURRENT_SOURCE_DIR}/Config.cmake.in
  "${CMAKE_CURRENT_BINARY_DIR}/s2Config.cmake"
  INSTALL_DESTINATION "share/s2/"
  NO_SET_AND_CHECK_MACRO
  NO_CHECK_REQUIRED_COMPONENTS_MACRO)

write_basic_package_version_file(
  "${CMAKE_CURRENT_BINARY_DIR}/s2ConfigVersion.cmake"
  VERSION ${PROJECT_VERSION}
  COMPATIBILITY SameMinorVersion)

install(FILES ${CMAKE_CURRENT_BINARY_DIR}/s2Config.cmake
              ${CMAKE_CURRENT_BINARY_DIR}/s2ConfigVersion.cmake
        DESTINATION "share/s2/")

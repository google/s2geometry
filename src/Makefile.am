# Make sure that when we re-make ./configure, we get the macros we need
ACLOCAL_AMFLAGS = -I m4

AM_CXXFLAGS = -std=c++11 -Wno-deprecated
AM_CPPFLAGS = -DGOOGLE_NAMESPACE=base \
              -D_START_GOOGLE_NAMESPACE_='namespace base {' \
              -D_END_GOOGLE_NAMESPACE_='}' -DNDEBUG

lib_LTLIBRARIES = libs2.la
check_LIBRARIES = libgtest.a

s2_includedir = $(includedir)/s2
s2_base_includedir = $(includedir)/s2/base
s2_strings_includedir = $(includedir)/s2/strings
s2_util_bits_includedir = $(includedir)/s2/util/bits
s2_util_coding_includedir = $(includedir)/s2/util/coding
s2_util_endian_includedir = $(includedir)/s2/util/endian
s2_util_gtl_includedir = $(includedir)/s2/util/gtl
s2_util_hash_includedir = $(includedir)/s2/util/hash
s2_util_math_includedir = $(includedir)/s2/util/math

s2_base_include_HEADERS = \
  base/atomicops.h \
  base/casts.h \
  base/int128.h \
  base/integral_types.h \
  base/macros.h \
  base/mutex.h \
  base/port.h \
  base/port_hash.h \
  base/scoped_ptr.h \
  base/scoped_ptr_internals.h \
  base/spinlock.h \
  base/stringprintf.h \
  base/template_util.h \
  base/type_traits.h

s2_strings_include_HEADERS = \
  strings/ascii_ctype.h \
  strings/case.h \
  strings/charset.h \
  strings/fastmem.h \
  strings/join_internal.h \
  strings/numbers.h \
  strings/split_internal.h \
  strings/strcat.h \
  strings/stringpiece.h \
  strings/stringpiece_utils.h \
  strings/strip.h

s2_util_bits_include_HEADERS = \
  util/bits/bits-internal-unknown.h \
  util/bits/bits-internal-windows.h \
  util/bits/bits.h

s2_util_coding_include_HEADERS = \
  util/coding/varint.h

s2_util_endian_include_HEADERS = \
  util/endian/endian.h

s2_util_gtl_include_HEADERS = \
  util/gtl/algorithm.h \
  util/gtl/container_logging.h \
  util/gtl/inlined_vector.h \
  util/gtl/manual_constructor.h

s2_util_hash_include_HEADERS = \
  util/hash/builtin_type_hash.h \
  util/hash/city.h \
  util/hash/farmhash.h \
  util/hash/hash.h \
  util/hash/hash128to64.h \
  util/hash/jenkins.h \
  util/hash/jenkins_lookup2.h \
  util/hash/murmur.h \
  util/hash/string_hash.h

s2_util_math_include_HEADERS = \
  util/math/mathlimits.h \
  util/math/mathutil.h \
  util/math/matrix3x3.h \
  util/math/vector2.h \
  util/math/vector3.h \
  util/math/vector3_hash.h \
  util/math/vector4.h

s2_include_HEADERS = \
  r1interval.h \
  r2.h \
  r2rect.h \
  s1chordangle.h \
  s1angle.h \
  s1interval.h \
  s2.h \
  s2cap.h \
  s2cell.h \
  s2cellid.h \
  s2closestedgequery.h \
  s2cellunion.h \
  s2edgequery.h \
  s2edgeutil.h \
  s2error.h \
  s2latlng.h \
  s2latlngrect.h \
  s2loop.h \
  s2paddedcell.h \
  s2pointcompression.h \
  s2pointregion.h \
  s2polygon.h \
  s2polygonbuilder.h \
  s2polyline.h \
  s2r2rect.h \
  s2region.h \
  s2regioncoverer.h \
  s2regionintersection.h \
  s2regionunion.h \
  s2shapeindex.h \
  s2shapeutil.h \
  s2testing.h

# TODO(jrosenstock): split out non-s2 libraries / testing library
libs2_la_SOURCES = \
  base/int128.cc \
  base/stringprintf.cc \
  base/strtoint.cc \
  strings/ascii_ctype.cc \
  strings/case.cc \
  strings/charset.cc \
  strings/escaping.cc \
  strings/join.cc \
  strings/memutil.cc \
  strings/numbers.cc \
  strings/serialize.cc \
  strings/split.cc \
  strings/strcat.cc \
  strings/stringpiece.cc \
  strings/stringpiece_utils.cc \
  strings/strip.cc \
  strings/util.cc \
  util/bits/bit-interleave.cc \
  util/bits/bits.cc \
  util/coding/coder.cc \
  util/coding/varint.cc \
  util/hash/city.cc \
  util/hash/hash.cc \
  util/hash/jenkins.cc \
  util/hash/murmur.cc \
  util/math/mathutil.cc \
  util/math/mathlimits.cc \
  util/math/exactfloat/exactfloat.cc \
  s1angle.cc \
  s1chordangle.cc \
  s2.cc \
  s2cellid.cc \
  s2latlng.cc \
  r2rect.cc \
  s1interval.cc \
  s2cap.cc \
  s2cell.cc \
  s2cellunion.cc \
  s2closestedgequery.cc \
  s2edgequery.cc \
  s2edgeutil.cc \
  s2error.cc \
  s2latlngrect.cc \
  s2loop.cc \
  s2paddedcell.cc \
  s2pointcompression.cc \
  s2pointregion.cc \
  s2polygon.cc \
  s2polygonbuilder.cc \
  s2polyline.cc \
  s2r2rect.cc \
  s2region.cc \
  s2regioncoverer.cc \
  s2regionintersection.cc \
  s2regionunion.cc \
  s2shapeindex.cc \
  s2shapeutil.cc \
  s2testing.cc

TEST_CXXFLAGS = \
  -I$(top_srcdir)/gtest/include -I$(top_srcdir)/gtest \
  -std=c++11 -Wno-deprecated
TEST_LDADD = $(lib_LTLIBRARIES) $(check_LIBRARIES) -lpthread

libgtest_a_SOURCES = \
  gtest/src/gtest-all.cc \
  gtest/src/gtest_main.cc
libgtest_a_CXXFLAGS = $(TEST_CXXFLAGS)

TESTS = \
  r1interval_test \
  r2rect_test \
  s1angle_test \
  s1interval_test \
  s2_test \
  s2cap_test \
  s2cell_test \
  s2cellid_test \
  s2cellunion_test \
  s2edgequery_test \
  s2edgeutil_test \
  s2latlng_test \
  s2latlngrect_test \
  s2loop_test \
  s2paddedcell_test \
  s2pointcompression_test \
  s2pointregion_test \
  s2polygon_test \
  s2polygonbuilder_test \
  s2polyline_test \
  s2r2rect_test \
  s2regioncoverer_test \
  s2regionunion_test \
  s2shapeindex_test \
  s2shapeutil_test \
  s2testing_test

check_PROGRAMS = $(TESTS)

r1interval_test_SOURCES = r1interval_test.cc
r1interval_test_CXXFLAGS = $(TEST_CXXFLAGS)
r1interval_test_LDADD = $(TEST_LDADD)

r2rect_test_SOURCES = r2rect_test.cc
r2rect_test_CXXFLAGS = $(TEST_CXXFLAGS)
r2rect_test_LDADD = $(TEST_LDADD)

s1angle_test_SOURCES = s1angle_test.cc
s1angle_test_CXXFLAGS = $(TEST_CXXFLAGS)
s1angle_test_LDADD = $(TEST_LDADD)

s1interval_test_SOURCES = s1interval_test.cc
s1interval_test_CXXFLAGS = $(TEST_CXXFLAGS)
s1interval_test_LDADD = $(TEST_LDADD)

s2_test_SOURCES = s2_test.cc
s2_test_CXXFLAGS = $(TEST_CXXFLAGS)
s2_test_LDADD = $(TEST_LDADD)

s2cap_test_SOURCES = s2cap_test.cc
s2cap_test_CXXFLAGS = $(TEST_CXXFLAGS)
s2cap_test_LDADD = $(TEST_LDADD)

s2cell_test_SOURCES = s2cell_test.cc
s2cell_test_CXXFLAGS = $(TEST_CXXFLAGS)
s2cell_test_LDADD = $(TEST_LDADD)

s2cellid_test_SOURCES = s2cellid_test.cc
s2cellid_test_CXXFLAGS = $(TEST_CXXFLAGS)
s2cellid_test_LDADD = $(TEST_LDADD)

s2cellunion_test_SOURCES = s2cellunion_test.cc
s2cellunion_test_CXXFLAGS = $(TEST_CXXFLAGS)
s2cellunion_test_LDADD = $(TEST_LDADD)

s2edgequery_test_SOURCES = s2edgequery_test.cc
s2edgequery_test_CXXFLAGS = $(TEST_CXXFLAGS)
s2edgequery_test_LDADD = $(TEST_LDADD)

s2edgeutil_test_SOURCES = s2edgeutil_test.cc
s2edgeutil_test_CXXFLAGS = $(TEST_CXXFLAGS)
s2edgeutil_test_LDADD = $(TEST_LDADD)

s2latlng_test_SOURCES = s2latlng_test.cc
s2latlng_test_CXXFLAGS = $(TEST_CXXFLAGS)
s2latlng_test_LDADD = $(TEST_LDADD)

s2latlngrect_test_SOURCES = s2latlngrect_test.cc
s2latlngrect_test_CXXFLAGS = $(TEST_CXXFLAGS)
s2latlngrect_test_LDADD = $(TEST_LDADD)

s2loop_test_SOURCES = s2loop_test.cc
s2loop_test_CXXFLAGS = $(TEST_CXXFLAGS)
s2loop_test_LDADD = $(TEST_LDADD)

s2paddedcell_test_SOURCES = s2paddedcell_test.cc
s2paddedcell_test_CXXFLAGS = $(TEST_CXXFLAGS)
s2paddedcell_test_LDADD = $(TEST_LDADD)

s2pointcompression_test_SOURCES = s2pointcompression_test.cc
s2pointcompression_test_CXXFLAGS = $(TEST_CXXFLAGS)
s2pointcompression_test_LDADD = $(TEST_LDADD)

s2pointregion_test_SOURCES = s2pointregion_test.cc
s2pointregion_test_CXXFLAGS = $(TEST_CXXFLAGS)
s2pointregion_test_LDADD = $(TEST_LDADD)

s2polygon_test_SOURCES = s2polygon_test.cc
s2polygon_test_CXXFLAGS = $(TEST_CXXFLAGS)
s2polygon_test_LDADD = $(TEST_LDADD)

s2polygonbuilder_test_SOURCES = s2polygonbuilder_test.cc
s2polygonbuilder_test_CXXFLAGS = $(TEST_CXXFLAGS)
s2polygonbuilder_test_LDADD = $(TEST_LDADD)

s2polyline_test_SOURCES = s2polyline_test.cc
s2polyline_test_CXXFLAGS = $(TEST_CXXFLAGS)
s2polyline_test_LDADD = $(TEST_LDADD)

s2r2rect_test_SOURCES = s2r2rect_test.cc
s2r2rect_test_CXXFLAGS = $(TEST_CXXFLAGS)
s2r2rect_test_LDADD = $(TEST_LDADD)

s2regioncoverer_test_SOURCES = s2regioncoverer_test.cc
s2regioncoverer_test_CXXFLAGS = $(TEST_CXXFLAGS)
s2regioncoverer_test_LDADD = $(TEST_LDADD)

s2regionunion_test_SOURCES = s2regionunion_test.cc
s2regionunion_test_CXXFLAGS = $(TEST_CXXFLAGS)
s2regionunion_test_LDADD = $(TEST_LDADD)

s2shapeindex_test_SOURCES = s2shapeindex_test.cc
s2shapeindex_test_CXXFLAGS = $(TEST_CXXFLAGS)
s2shapeindex_test_LDADD = $(TEST_LDADD)

s2shapeutil_test_SOURCES = s2shapeutil_test.cc
s2shapeutil_test_CXXFLAGS = $(TEST_CXXFLAGS)
s2shapeutil_test_LDADD = $(TEST_LDADD)

s2testing_test_SOURCES = s2testing_test.cc
s2testing_test_CXXFLAGS = $(TEST_CXXFLAGS)
s2testing_test_LDADD = $(TEST_LDADD)

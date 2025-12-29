// Copyright 2018 Google Inc. All Rights Reserved.
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//     http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS-IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.
//

// Author: ericv@google.com (Eric Veach)

#include "s2/encoded_s2shape_index.h"

#include <algorithm>
#include <cstddef>
#include <cstring>
#include <memory>
#include <random>
#include <utility>
#include <vector>

#include <benchmark/benchmark.h>
#include <gtest/gtest.h>
#include "absl/base/call_once.h"
#include "absl/flags/flag.h"
#include "absl/log/absl_check.h"
#include "absl/log/absl_log.h"
#include "absl/log/log_streamer.h"
#include "absl/random/random.h"
#include "absl/strings/escaping.h"
#include "absl/strings/str_cat.h"
#include "s2/util/coding/coder.h"
#include "s2/mutable_s2shape_index.h"
#include "s2/s1angle.h"
#include "s2/s2builder.h"
#include "s2/s2builder_layer.h"
#include "s2/s2builderutil_s2polyline_layer.h"
#include "s2/s2builderutil_snap_functions.h"
#include "s2/s2cap.h"
#include "s2/s2cell_id.h"
#include "s2/s2closest_edge_query.h"
#include "s2/s2coder.h"
#include "s2/s2contains_point_query.h"
#include "s2/s2edge_distances.h"
#include "s2/s2error.h"
#include "s2/s2fractal.h"
#include "s2/s2latlng.h"
#include "s2/s2lax_polygon_shape.h"
#include "s2/s2lax_polyline_shape.h"
#include "s2/s2loop.h"
#include "s2/s2point.h"
#include "s2/s2point_vector_shape.h"
#include "s2/s2pointutil.h"
#include "s2/s2polygon.h"
#include "s2/s2polyline.h"
#include "s2/s2random.h"
#include "s2/s2shape.h"
#include "s2/s2shapeutil_coding.h"
#include "s2/s2shapeutil_testing.h"
#include "s2/s2testing.h"
#include "s2/s2text_format.h"
#include "s2/thread_testing.h"
#include "s2/util/math/matrix3x3.h"
#include "s2/util/random/shared_bit_gen.h"

using absl::StrCat;
using s2builderutil::S2CellIdSnapFunction;
using s2builderutil::S2PolylineLayer;
using std::make_unique;
using std::max;
using std::string;
using std::unique_ptr;
using std::vector;

template <class Shape>
bool DecodeHomegeneousShapeIndex(EncodedS2ShapeIndex* index, Decoder* decoder) {
  return index->Init(decoder,
                     s2shapeutil::HomogeneousShapeFactory<Shape>(decoder));
}

template <class InShape, class OutShape>
void TestEncodedS2ShapeIndex(const MutableS2ShapeIndex& expected,
                             size_t expected_bytes) {
  Encoder encoder;
  s2shapeutil::EncodeHomogeneousShapes<InShape>(expected, &encoder);
  size_t shapes_bytes = encoder.length();
  expected.Encode(&encoder);
  EXPECT_EQ(expected_bytes, encoder.length() - shapes_bytes);
  Decoder decoder(encoder.base(), encoder.length());
  EncodedS2ShapeIndex actual;
  ASSERT_TRUE(DecodeHomegeneousShapeIndex<OutShape>(&actual, &decoder));
  EXPECT_EQ(expected.options().max_edges_per_cell(),
            actual.options().max_edges_per_cell());
  s2testing::ExpectEqual(expected, actual);

  // Make sure that re-encoding the index gives us back the original bytes.
  Encoder new_encoder;
  actual.Encode(&new_encoder);
  EXPECT_EQ(encoder.length() - shapes_bytes, new_encoder.length());

  const char* index_base = encoder.base() + shapes_bytes;
  EXPECT_EQ(memcmp(index_base, new_encoder.base(), new_encoder.length()), 0);
}

TEST(EncodedS2ShapeIndex, Empty) {
  MutableS2ShapeIndex index;
  TestEncodedS2ShapeIndex<S2LaxPolylineShape, S2LaxPolylineShape>(index, 4);
}

TEST(EncodedS2ShapeIndex, OneEdge) {
  MutableS2ShapeIndex index;
  index.Add(s2textformat::MakeLaxPolylineOrDie("1:1, 2:2"));
  TestEncodedS2ShapeIndex<S2LaxPolylineShape, S2LaxPolylineShape>(index, 8);
}

TEST(EncodedS2ShapeIndex, RegularLoops) {
  struct TestCase {
    int num_edges;
    size_t expected_bytes;
  };
  vector<TestCase> test_cases = {
    {4, 8},
    {8, 8},
    {16, 16},
    {64, 77},
    {256, 327},
    {4096, 8813},
    {65536, 168291},
  };
  for (const auto& test_case : test_cases) {
    MutableS2ShapeIndex index;
    SCOPED_TRACE(StrCat("num_edges = ", test_case.num_edges));
    S2Polygon polygon(S2Loop::MakeRegularLoop(S2Point(3, 2, 1).Normalize(),
                                              S1Angle::Degrees(0.1),
                                              test_case.num_edges));
    index.Add(make_unique<S2LaxPolygonShape>(polygon));
    TestEncodedS2ShapeIndex<S2LaxPolygonShape, EncodedS2LaxPolygonShape>(
        index, test_case.expected_bytes);
  }
}

#if !defined( __EMSCRIPTEN__) && !(defined(__ANDROID__) && defined(__i386__))
// TODO(b/232496949): This test relies on `mt19937_64` return values because
// it tests an exact encoded byte size.  Either change it to accept a range
// of sizes, or decode and check either the number of shapes, or possibly
// the points themselves by resetting the RNG state.
TEST(EncodedS2ShapeIndex, OverlappingPointClouds) {
  std::mt19937_64 bitgen;
  struct TestCase {
    int num_shapes, num_points_per_shape;
    size_t expected_bytes;
  };
  vector<TestCase> test_cases = {
      {1, 50, 85},
      {2, 100, 600},
      {4, 100, 1426},
  };
  S2Cap cap(S2Point(0.1, -0.4, 0.3).Normalize(), S1Angle::Degrees(1));
  for (const auto& test_case : test_cases) {
    MutableS2ShapeIndex index;
    SCOPED_TRACE(StrCat("num_shapes = ", test_case.num_shapes));
    for (int i = 0; i < test_case.num_shapes; ++i) {
      vector<S2Point> points;
      for (int j = 0; j < test_case.num_points_per_shape; ++j) {
        points.push_back(s2random::SamplePoint(bitgen, cap));
      }
      index.Add(make_unique<S2PointVectorShape>(points));
    }
    TestEncodedS2ShapeIndex<S2PointVectorShape, EncodedS2PointVectorShape>(
        index, test_case.expected_bytes);
  }
}

// TODO(b/232496949): This test relies on `mt19937_64` return values.
TEST(EncodedS2ShapeIndex, OverlappingPolylines) {
  std::mt19937_64 bitgen;
  struct TestCase {
    int num_shapes, num_shape_edges;
    size_t expected_bytes;
  };
  vector<TestCase> test_cases = {
      {2, 50, 128},
      {10, 50, 808},
      {20, 50, 2499},
  };
  S2Cap cap(S2Point(-0.2, -0.3, 0.4).Normalize(), S1Angle::Degrees(0.1));
  for (const auto& test_case : test_cases) {
    S1Angle edge_len = 2 * cap.GetRadius() / test_case.num_shape_edges;
    MutableS2ShapeIndex index;
    SCOPED_TRACE(StrCat("num_shapes = ", test_case.num_shapes));
    for (int i = 0; i < test_case.num_shapes; ++i) {
      S2Point a = s2random::SamplePoint(bitgen, cap);
      S2Point b = s2random::Point(bitgen);
      vector<S2Point> vertices;
      int n = test_case.num_shape_edges;
      for (int j = 0; j <= n; ++j) {
        vertices.push_back(S2::GetPointOnLine(a, b, j * edge_len));
      }
      index.Add(make_unique<S2LaxPolylineShape>(vertices));
    }
    TestEncodedS2ShapeIndex<S2LaxPolylineShape, EncodedS2LaxPolylineShape>(
        index, test_case.expected_bytes);
  }
}

// TODO(b/232496949): This test relies on `mt19937_64` return values.
TEST(EncodedS2ShapeIndex, OverlappingLoops) {
  std::mt19937_64 bitgen;
  struct TestCase {
    int num_shapes, max_edges_per_loop;
    size_t expected_bytes;
  };
  vector<TestCase> test_cases = {
      {2, 250, 412},
      {5, 250, 1208},
      {25, 50, 3162},
  };
  S2Cap cap(S2Point(-0.1, 0.25, 0.2).Normalize(), S1Angle::Degrees(3));
  for (const auto& test_case : test_cases) {
    MutableS2ShapeIndex index;
    SCOPED_TRACE(StrCat("num_shapes = ", test_case.num_shapes));
    for (int i = 0; i < test_case.num_shapes; ++i) {
      S2Point center = s2random::SamplePoint(bitgen, cap);
      double radius_fraction = absl::Uniform(bitgen, 0.0, 1.0);
      // Scale the number of edges so that they are all about the same length
      // (similar to modeling all geometry at a similar resolution).
      int num_edges = max(3.0, test_case.max_edges_per_loop * radius_fraction);
      S2Polygon polygon(S2Loop::MakeRegularLoop(
          center, cap.GetRadius() * radius_fraction, num_edges));
      index.Add(make_unique<S2LaxPolygonShape>(polygon));
    }
    TestEncodedS2ShapeIndex<S2LaxPolygonShape, EncodedS2LaxPolygonShape>(
        index, test_case.expected_bytes);
  }
}
#endif  // defined(__EMSCRIPTEN__)

// Like S2PolylineLayer, but converts the polyline to an S2LaxPolylineShape
// and adds it to an S2ShapeIndex (if the polyline is non-empty).
class IndexedLaxPolylineLayer : public S2Builder::Layer {
 public:
  using Options = S2PolylineLayer::Options;
  explicit IndexedLaxPolylineLayer(MutableS2ShapeIndex* index,
                                   const Options& options = Options())
      : index_(index), polyline_(make_unique<S2Polyline>()),
        layer_(polyline_.get(), options) {}

  GraphOptions graph_options() const override {
    return layer_.graph_options();
  }

  void Build(const Graph& g, S2Error* error) override {
    layer_.Build(g, error);
    if (error->ok() && polyline_->num_vertices() > 0) {
      index_->Add(make_unique<S2LaxPolylineShape>(*polyline_));
    }
  }

 private:
  MutableS2ShapeIndex* index_;
  unique_ptr<S2Polyline> polyline_;
  S2PolylineLayer layer_;
};

TEST(EncodedS2ShapeIndex, SnappedFractalPolylines) {
  absl::BitGen bitgen(S2Testing::MakeTaggedSeedSeq(
      "SNAPPED_FRACTAL_POLYLINES", absl::LogInfoStreamer(__FILE__, __LINE__).stream()));
  MutableS2ShapeIndex index;
  S2Builder builder{S2Builder::Options{S2CellIdSnapFunction()}};
  for (int i = 0; i < 5; ++i) {
    builder.StartLayer(make_unique<IndexedLaxPolylineLayer>(&index));
    S2Fractal fractal(bitgen);
    fractal.SetLevelForApproxMaxEdges(3 * 256);
    auto frame = S2::GetFrame(S2LatLng::FromDegrees(10, i).ToPoint());
    auto loop = fractal.MakeLoop(frame, S1Angle::Degrees(0.1));
    vector<S2Point> vertices;
    S2Testing::AppendLoopVertices(*loop, &vertices);
    S2Polyline polyline(vertices);
    builder.AddPolyline(polyline);
  }
  S2Error error;
  ASSERT_TRUE(builder.Build(&error)) << error.message();
  TestEncodedS2ShapeIndex<S2LaxPolylineShape, EncodedS2LaxPolylineShape>(
      index, 8698);
}

// A test that repeatedly minimizes "index_" in one thread and then reads the
// index_ concurrently from several other threads.  When all threads have
// finished reading, the first thread minimizes the index again.
//
// Note that Minimize() is non-const and therefore does not need to be tested
// concurrently with the const methods.
class LazyDecodeTest : public s2testing::ReaderWriterTest {
 public:
  LazyDecodeTest() {
    // We generate one shape per dimension.  Each shape has vertices uniformly
    // distributed across the sphere, and the vertices for each dimension are
    // different.  Having fewer cells in the index is more likely to trigger
    // race conditions, and so shape 0 has 384 points, shape 1 is a polyline
    // with 96 vertices, and shape 2 is a polygon with 24 vertices.
    MutableS2ShapeIndex input;
    for (int dim = 0; dim < 3; ++dim) {
      int level = 3 - dim;  // See comments above.
      vector<S2Point> vertices;
      for (auto id = S2CellId::Begin(level);
           id != S2CellId::End(level); id = id.next()) {
        vertices.push_back(id.ToPoint());
      }
      switch (dim) {
        case 0: input.Add(make_unique<S2PointVectorShape>(vertices)); break;
        case 1: input.Add(make_unique<S2LaxPolylineShape>(vertices)); break;
        default:
          input.Add(make_unique<S2LaxPolygonShape>(
              vector<vector<S2Point>>{std::move(vertices)}));
          break;
      }
    }
    Encoder encoder;
    ABSL_CHECK(s2shapeutil::CompactEncodeTaggedShapes(input, &encoder));
    input.Encode(&encoder);
    encoded_.assign(encoder.base(), encoder.length());

    Decoder decoder(encoded_.data(), encoded_.size());
    ABSL_CHECK(
        index_.Init(&decoder, s2shapeutil::LazyDecodeShapeFactory(&decoder)));
  }

  void WriteOp() override {
    index_.Minimize();
  }

  void ReadOp() override {
    util_random::SharedBitGen bitgen;
    S2ClosestEdgeQuery query(&index_);
    for (int iter = 0; iter < 10; ++iter) {
      S2ClosestEdgeQuery::PointTarget target(s2random::Point(bitgen));
      query.FindClosestEdge(&target);
    }
  }

 protected:
  string encoded_;
  EncodedS2ShapeIndex index_;
};

TEST(EncodedS2ShapeIndex, LazyDecode) {
  // Ensure that lazy decoding is thread-safe.  In other words, make sure that
  // nothing bad happens when multiple threads call "const" methods that cause
  // index and/or shape data to be decoded.
  LazyDecodeTest test;

  // The number of readers should be large enough so that it is likely that
  // several readers will be running at once (with a multiple-core CPU).
  constexpr int kNumReaders = 8;
  constexpr int kIters = 1000;
  test.Run(kNumReaders, kIters);
}

TEST(EncodedS2ShapeIndex, JavaByteCompatibility) {
  MutableS2ShapeIndex expected;
  expected.Add(make_unique<S2Polyline::OwningShape>(
      s2textformat::MakePolylineOrDie("0:0, 1:1")));
  expected.Add(make_unique<S2Polyline::OwningShape>(
      s2textformat::MakePolylineOrDie("1:1, 2:2")));
  expected.Release(0);

  // bytes is the encoded data of an S2ShapeIndex with a null shape and a
  // polyline with one edge. It was derived by base-16 encoding the buffer of
  // an encoder to which expected was encoded.
  string bytes = absl::HexStringToBytes(
      "100036020102000000B4825F3C81FDEF3F27DCF7C958DE913F1EDD892B0BDF913FFC7FB8"
      "B805F6EF3F28516A6D8FDBA13F27DCF7C958DEA13F28C809010408020010");
  Decoder decoder(bytes.data(), bytes.length());
  MutableS2ShapeIndex actual;
  ASSERT_TRUE(
      actual.Init(&decoder, s2shapeutil::FullDecodeShapeFactory(&decoder)));

  s2testing::ExpectEqual(expected, actual);
}

// The number of different polygons used for benchmarking purposes.
static constexpr int kNumBenchmarkPolygons = 3;

// Returns a pointer to a cached polygon of the given type and snap state.
//
// REQUIRES: type < kNumBenchmarkPolygons
const S2Polygon* GetBenchmarkPolygon(int type, bool snapped) {
  // Store an array of pointers to cache test polygons.  We'll never call
  // destructors so this doesn't violate the trivial destruction rule.
  struct CachedPolygon {
    absl::once_flag once;
    const S2Polygon* polygon;
    const S2Polygon* snapped;
  };

  static CachedPolygon cache[kNumBenchmarkPolygons] = {};

  ABSL_CHECK_GE(type, 0);
  ABSL_CHECK_LT(type, kNumBenchmarkPolygons);

  // Initialize cache entry if it hasn't been done yet.
  absl::call_once(cache[type].once, [type]() {
    const string seed_str =
        StrCat("BENCHMARK_POLYGON", absl::GetFlag(FLAGS_s2_random_seed), type);
    std::seed_seq seed(seed_str.begin(), seed_str.end());
    std::mt19937_64 bitgen(seed);

    S2Polygon* polygon;
    switch (type) {
      case 0: {
        // A small regular loop.
        polygon = new S2Polygon(S2Loop::MakeRegularLoop(
            s2random::Point(bitgen), S2Testing::KmToAngle(10), 8));
        break;
      }
      case 1: {
        // A medium-sized fractal loop, e.g. a small natural feature.
        S2Fractal fractal(bitgen);
        fractal.SetLevelForApproxMaxEdges(3 * 256);
        Matrix3x3_d frame = s2random::Frame(bitgen);
        auto loop = fractal.MakeLoop(frame, S2Testing::KmToAngle(10));
        polygon = new S2Polygon(std::move(loop));
        break;
      }
      case 2: {
        // A large number of small loops.
        polygon = new S2Polygon();
        S2Testing::ConcentricLoopsPolygon(S2Point(1, 0, 0), 100, 10, polygon);
        break;
      }
      default:
        ABSL_LOG(FATAL) << "Unsupported benchmark polygon index: " << type;
    }

    // Build snapped variant of polygon
    S2Polygon* snapped = new S2Polygon();
    snapped->InitToSnapped(*polygon, S2CellId::kMaxLevel);

    cache[type].polygon = polygon;
    cache[type].snapped = snapped;
  });

  // Cache will always be initialized at this point, safe to read values.
  return snapped ? cache[type].snapped : cache[type].polygon;
}

// Returns a pointer to a cached polygon specified by the benchmark state.
//
// REQUIRES: state.range(0) < kNumBenchmarkPolygons
const S2Polygon* GetBenchmarkPolygon(const benchmark::State& state) {
  return GetBenchmarkPolygon(state.range(0), state.range(1));
}

// Returns an owning pointer to a polygon specified by the benchmark state.
// Since the returned pointer is mutable, polygons are cloned from the cached
// values returned by `GetBenchmarkPolygon`
//
// REQUIRES: state.range(0) < kNumBenchmarkPolygons
unique_ptr<S2Polygon> GetMutableBenchmarkPolygon(
    const benchmark::State& state) {
  return unique_ptr<S2Polygon>(GetBenchmarkPolygon(state)->Clone());
}

// Returns a pointer to a lax polygon specified by the benchmark state.
//
// REQUIRES: state.range(0) < kNumBenchmarkPolygons
const S2LaxPolygonShape* GetLaxBenchmarkPolygon(const benchmark::State& state) {
  // Store an array of pointers to cache test polygons.  We'll never call
  // destructors so this doesn't violate the trivial destruction rule.
  static absl::once_flag once;
  static const S2LaxPolygonShape* cache[kNumBenchmarkPolygons][2] = {};

  // Construct the lax polygon cache once.
  absl::call_once(once, []() {
    for (int type = 0; type < kNumBenchmarkPolygons; type++) {
      for (int snapped = 0; snapped <= 1; snapped++) {
        const S2Polygon* polygon = GetBenchmarkPolygon(type, snapped);
        cache[type][snapped] = new S2LaxPolygonShape(*polygon);
      }
    }
  });

  int type = state.range(0);
  bool snapped = state.range(1);

  ABSL_CHECK_GE(type, 0);
  ABSL_CHECK_LT(type, kNumBenchmarkPolygons);

  // Cache will always be initialized at this point, safe to read values.
  return cache[type][snapped];
}

// Runs the given benchmark for all benchmark polygons using the given value
// of the "snapped" parameter.
template <bool snap>
void BenchmarkArgsSnapped(benchmark::internal::Benchmark* benchmark) {
  for (int i = 0; i < kNumBenchmarkPolygons; ++i) {
    benchmark->ArgPair(i, snap);
  }
}

// Runs the given benchmark for all benchmark polygons, with and without
// snapping (which controls whether the FAST or COMPACT encoding is used).
void BenchmarkArgs(benchmark::internal::Benchmark* benchmark) {
  BenchmarkArgsSnapped<false>(benchmark);
  BenchmarkArgsSnapped<true>(benchmark);
}

// Measure the time to build and destruct a MutableS2ShapeIndex.
static void BM_MutableIndexBuildDestruct(benchmark::State& state) {
  auto polygon = GetBenchmarkPolygon(state);
  state.SetLabel(StrCat("vertices = ", polygon->num_vertices(),
                        ", loops = ", polygon->num_loops()));
  for (auto _ : state) {
    MutableS2ShapeIndex index;
    index.Add(make_unique<S2Polygon::Shape>(polygon));
    index.ForceBuild();
  }
}
BENCHMARK(BM_MutableIndexBuildDestruct)->Apply(BenchmarkArgsSnapped<false>);

// Measure the time to decode and destruct a MutableS2ShapeIndex.
static void BM_MutableIndexDecodeDestruct(benchmark::State& state) {
  Encoder encoder;
  {
    MutableS2ShapeIndex index;
    index.Add(make_unique<S2LaxPolygonShape>(*GetBenchmarkPolygon(state)));
    ABSL_CHECK(s2shapeutil::FastEncodeTaggedShapes(index, &encoder));
    index.Encode(&encoder);
  }
  for (auto _ : state) {
    Decoder decoder(encoder.base(), encoder.length());
    MutableS2ShapeIndex index;
    ABSL_CHECK(
        index.Init(&decoder, s2shapeutil::LazyDecodeShapeFactory(&decoder)));
  }
}
BENCHMARK(BM_MutableIndexDecodeDestruct)->Apply(BenchmarkArgsSnapped<false>);

// Measure the time to initialize and destruct an EncodedS2ShapeIndex.
static void BM_EncodedIndexInitDestruct(benchmark::State& state) {
  Encoder encoder;
  {
    MutableS2ShapeIndex index;
    index.Add(make_unique<S2LaxPolygonShape>(*GetBenchmarkPolygon(state)));
    ABSL_CHECK(s2shapeutil::FastEncodeTaggedShapes(index, &encoder));
    int index_start = encoder.length();
    index.Encode(&encoder);
    Decoder decoder(encoder.base(), encoder.length());
    state.SetLabel(StrCat("index_bytes = ",
                          encoder.length() - index_start));
  }
  for (auto _ : state) {
    Decoder decoder(encoder.base(), encoder.length());
    EncodedS2ShapeIndex index;
    ABSL_CHECK(
        index.Init(&decoder, s2shapeutil::LazyDecodeShapeFactory(&decoder)));
  }
}
BENCHMARK(BM_EncodedIndexInitDestruct)->Apply(BenchmarkArgsSnapped<false>);

// Measure the time to decode and destruct a snapped/unsnapped S2Polygon.
static void BM_S2PolygonDecodeDestructSnapped(benchmark::State& state) {
  Encoder encoder;
  GetBenchmarkPolygon(state)->Encode(&encoder);
  state.SetLabel(StrCat("bytes = ", encoder.length()));
  for (auto _ : state) {
    Decoder decoder(encoder.base(), encoder.length());
    S2Polygon polygon;
    ABSL_CHECK(polygon.Decode(&decoder));
  }
}
BENCHMARK(BM_S2PolygonDecodeDestructSnapped)->Apply(BenchmarkArgs);

// Measure the time to fully decode and destruct a snapped/unsnapped
// S2LaxPolygonShape.
static void BM_LaxPolygonDecodeDestructSnapped(
    benchmark::State& state) {
  Encoder encoder;
  GetLaxBenchmarkPolygon(state)->Encode(&encoder,
                                        s2coding::CodingHint::COMPACT);
  state.SetLabel(StrCat("bytes = ", encoder.length()));
  for (auto _ : state) {
    Decoder decoder(encoder.base(), encoder.length());
    S2LaxPolygonShape polygon;
    ABSL_CHECK(polygon.Init(&decoder));
  }
}
BENCHMARK(BM_LaxPolygonDecodeDestructSnapped)->Apply(BenchmarkArgs);

// Measure the time to initialize and destruct a snapped/unsnapped
// EncodedS2LaxPolygonShape.
static void BM_EncodedPolygonInitDestructSnapped(
    benchmark::State& state) {
  Encoder encoder;
  GetLaxBenchmarkPolygon(state)->Encode(&encoder,
                                        s2coding::CodingHint::COMPACT);
  state.SetLabel(StrCat("bytes = ", encoder.length()));
  for (auto _ : state) {
    Decoder decoder(encoder.base(), encoder.length());
    EncodedS2LaxPolygonShape polygon;
    ABSL_CHECK(polygon.Init(&decoder));
  }
}
BENCHMARK(BM_EncodedPolygonInitDestructSnapped)->Apply(BenchmarkArgs);

// Measure the time to access an edge of an S2Polygon::Shape.
static void BM_S2PolygonGetEdge(benchmark::State& state) {
  S2Polygon::OwningShape shape(GetMutableBenchmarkPolygon(state));
  int i = 0, num_edges = shape.num_edges();
  for (auto _ : state) {
    (void) shape.edge(i);
    if (++i == num_edges) i = 0;
  }
}
BENCHMARK(BM_S2PolygonGetEdge)->Apply(BenchmarkArgsSnapped<false>);

static void BM_LaxPolygonGetEdge(benchmark::State& state) {
  const S2LaxPolygonShape& shape = *GetLaxBenchmarkPolygon(state);
  int i = 0, num_edges = shape.num_edges();
  for (auto _ : state) {
    benchmark::DoNotOptimize(shape.edge(i));
    if (++i == num_edges) i = 0;
  }
}
BENCHMARK(BM_LaxPolygonGetEdge)->Apply(BenchmarkArgsSnapped<false>);

// Measure the time to access an edge of a snapped/unsnapped
// EncodedS2LaxPolygonShape.
static void BM_EncodedPolygonGetEdge(benchmark::State& state) {
  Encoder encoder;
  GetLaxBenchmarkPolygon(state)->Encode(&encoder,
                                        s2coding::CodingHint::COMPACT);
  state.SetLabel(StrCat("bytes = ", encoder.length()));
  Decoder decoder(encoder.base(), encoder.length());
  EncodedS2LaxPolygonShape shape;
  ABSL_CHECK(shape.Init(&decoder));
  int i = 0, num_edges = shape.num_edges();
  for (auto _ : state) {
    (void) shape.edge(i);
    if (++i == num_edges) i = 0;
  }
}
BENCHMARK(BM_EncodedPolygonGetEdge)->Apply(BenchmarkArgs);

// Measure the time to check whether the polygon contains one of its vertices,
// using a MutableS2ShapeIndex containing an S2Polygon::Shape.
static void BM_ContainsPointMutableIndexS2Polygon(benchmark::State& state) {
  auto polygon = GetBenchmarkPolygon(state);
  MutableS2ShapeIndex index;
  index.Add(make_unique<S2Polygon::Shape>(polygon));
  index.ForceBuild();
  int i = 0, num_edges = index.shape(0)->num_edges();
  for (auto _ : state) {
    S2Point test_point = index.shape(0)->edge(i).v0;
    if (S2ContainsPointQuery(&index).Contains(test_point)) ++i;
    if ((i += 10) >= num_edges) i = 0;
  }
}
BENCHMARK(BM_ContainsPointMutableIndexS2Polygon)
->Apply(BenchmarkArgsSnapped<false>);

// Measure the time to check whether the polygon contains one of its vertices,
// using an EncodedS2ShapeIndex containing a snapped/unsnapped
// EncodedS2LaxPolygonShape.
static void BM_ContainsPointEncodedIndexPolygon(benchmark::State& state) {
  Encoder encoder;
  {
    MutableS2ShapeIndex index;
    index.Add(make_unique<S2LaxPolygonShape>(*GetBenchmarkPolygon(state)));
    ABSL_CHECK(s2shapeutil::CompactEncodeTaggedShapes(index, &encoder));
    index.Encode(&encoder);
  }
  Decoder decoder(encoder.base(), encoder.length());
  EncodedS2ShapeIndex index;
  ABSL_CHECK(
      index.Init(&decoder, s2shapeutil::LazyDecodeShapeFactory(&decoder)));
  int i = 0, num_edges = index.shape(0)->num_edges();
  for (auto _ : state) {
    S2Point test_point = index.shape(0)->edge(i).v0;
    if (S2ContainsPointQuery(&index).Contains(test_point)) ++i;
    if ((i += 10) >= num_edges) {
      // In theory we should stop benchmark timing for this, but ironically
      // that throws off the benchmark results.
      i = 0;
      index.Minimize();
    }
  }
}
BENCHMARK(BM_ContainsPointEncodedIndexPolygon)->Apply(BenchmarkArgs);

// Measure the time to compute distance from the polygon to one of its
// vertices, using a MutableS2ShapeIndex containing an S2Polygon::Shape.
static void BM_VertexDistanceMutableIndexS2Polygon(benchmark::State& state) {
  auto polygon = GetBenchmarkPolygon(state);
  MutableS2ShapeIndex index;
  index.Add(make_unique<S2Polygon::Shape>(polygon));
  index.ForceBuild();
  int i = 0, num_edges = index.shape(0)->num_edges();
  for (auto _ : state) {
    S2ClosestEdgeQuery query(&index);
    S2Point v0 = index.shape(0)->edge(i).v0;
    S2Point test_point = (v0 + 1e-4 * S2::Ortho(v0)).Normalize();
    S2ClosestEdgeQuery::PointTarget target(test_point);
    query.GetDistance(&target);
    if ((i += 10) >= num_edges) i = 0;
  }
}
BENCHMARK(BM_VertexDistanceMutableIndexS2Polygon)
->Apply(BenchmarkArgsSnapped<false>);

// Measure the time to compute the distance from the polygon to one of its
// vertices, using an EncodedS2ShapeIndex containing a snapped/unsnapped
// EncodedS2LaxPolygonShape.
static void BM_VertexDistanceEncodedIndexPolygon(benchmark::State& state) {
  Encoder encoder;
  {
    MutableS2ShapeIndex index;
    index.Add(make_unique<S2LaxPolygonShape>(*GetBenchmarkPolygon(state)));
    ABSL_CHECK(s2shapeutil::CompactEncodeTaggedShapes(index, &encoder));
    index.Encode(&encoder);
  }
  Decoder decoder(encoder.base(), encoder.length());
  EncodedS2ShapeIndex index;
  ABSL_CHECK(
      index.Init(&decoder, s2shapeutil::LazyDecodeShapeFactory(&decoder)));
  int i = 0, num_edges = index.shape(0)->num_edges();
  for (auto _ : state) {
    S2ClosestEdgeQuery query(&index);
    S2Point v0 = index.shape(0)->edge(i).v0;
    S2Point test_point = (v0 + 1e-4 * S2::Ortho(v0)).Normalize();
    S2ClosestEdgeQuery::PointTarget target(test_point);
    query.GetDistance(&target);
    if ((i += 10) >= num_edges) {
      // In theory we should stop benchmark timing for this, but ironically
      // that throws off the benchmark results.
      i = 0;
      index.Minimize();
    }
  }
}
BENCHMARK(BM_VertexDistanceEncodedIndexPolygon)->Apply(BenchmarkArgs);

// Measure the combined time to decode the snapped/unsnapped polygon, build an
// S2ShapeIndex for it, and check whether the polygon contains one of its
// vertices.
static void BM_DecodePolygonBuildIndexContainsPointSnapped(
    benchmark::State& state) {
  Encoder encoder;
  GetBenchmarkPolygon(state)->Encode(&encoder);
  for (auto _ : state) {
    Decoder decoder(encoder.base(), encoder.length());
    S2Polygon polygon;
    ABSL_CHECK(polygon.Decode(&decoder));
    MutableS2ShapeIndex index;
    index.Add(make_unique<S2Polygon::Shape>(&polygon));
    S2ContainsPointQuery query(&index);
    benchmark::DoNotOptimize(query.Contains(polygon.loop(0)->vertex(0)));
  }
}
BENCHMARK(BM_DecodePolygonBuildIndexContainsPointSnapped)->Apply(BenchmarkArgs);

// Measure the combined time to decode the snapped/unsnapped polygon and its
// S2ShapeIndex, and check whether the polygon contains one of its vertices.
static void BM_DecodeIndexAndPolygonContainsPointSnapped(
    benchmark::State& state) {
  Encoder encoder;
  {
    auto polygon = GetBenchmarkPolygon(state);
    MutableS2ShapeIndex index;
    index.Add(make_unique<S2Polygon::Shape>(polygon));
    ABSL_CHECK(s2shapeutil::CompactEncodeTaggedShapes(index, &encoder));
    index.Encode(&encoder);
  }
  for (auto _ : state) {
    Decoder decoder(encoder.base(), encoder.length());
    MutableS2ShapeIndex index;
    ABSL_CHECK(
        index.Init(&decoder, s2shapeutil::FullDecodeShapeFactory(&decoder)));
    S2ContainsPointQuery query(&index);
    benchmark::DoNotOptimize(query.Contains(index.shape(0)->edge(0).v0));
  }
}
BENCHMARK(BM_DecodeIndexAndPolygonContainsPointSnapped)->Apply(BenchmarkArgs);

// Measure the combined time to initialize encoded versions of the
// snapped/unsnapped polygon and its index, and check whether the polygon
// contains one of its vertices.
static void BM_EncodedIndexAndPolygonContainsPointSnapped(
    benchmark::State& state) {
  Encoder encoder;
  {
    MutableS2ShapeIndex index;
    index.Add(make_unique<S2LaxPolygonShape>(*GetBenchmarkPolygon(state)));
    ABSL_CHECK(s2shapeutil::CompactEncodeTaggedShapes(index, &encoder));
    index.Encode(&encoder);
  }
  for (auto _ : state) {
    Decoder decoder(encoder.base(), encoder.length());
    EncodedS2ShapeIndex index;
    ABSL_CHECK(
        index.Init(&decoder, s2shapeutil::LazyDecodeShapeFactory(&decoder)));
    S2ContainsPointQuery query(&index);
    benchmark::DoNotOptimize(query.Contains(index.shape(0)->edge(0).v0));
  }
}
BENCHMARK(BM_EncodedIndexAndPolygonContainsPointSnapped)->Apply(BenchmarkArgs);

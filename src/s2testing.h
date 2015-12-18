// Copyright 2005 Google Inc. All Rights Reserved.
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

#ifndef S2GEOMETRY_S2TESTING_H_
#define S2GEOMETRY_S2TESTING_H_

#include <string>
#include <vector>

#include <gflags/gflags.h>
#include "base/integral_types.h"
#include "base/macros.h"
#include "fpcontractoff.h"
#include "r2.h"
#include "s1angle.h"
#include "s2.h"
#include "s2cellid.h"
#include "util/math/matrix3x3.h"

class S1Angle;
class S2Cap;
class S2CellUnion;
class S2LatLng;
class S2LatLngRect;
class S2Loop;
class S2Polygon;
class S2Polyline;
class S2Region;

// You can optionally call S2Testing::rnd.Reset(FLAGS_s2_random_seed) at the
// start of a test or benchmark to ensure that its results do not depend on
// which other tests of benchmarks have run previously.  This can help with
// debugging.
//
// This flag currently does *not* affect the initial seed value for
// S2Testing::rnd.  TODO(user): Fix this.
DECLARE_int32(s2_random_seed);

// This class defines various static functions that are useful for writing
// unit tests.
class S2Testing {
 public:
  // Returns a vector of points shaped as a regular polygon with
  // num_vertices vertices, all on a circle of the specified angular
  // radius around the center.  The radius is the actual distance from
  // the center to the circle along the sphere.
  //
  // If you want to construct a regular polygon, try this:
  //   S2Polygon polygon(S2Loop::MakeRegularLoop(center, radius, num_vertices));
  static std::vector<S2Point> MakeRegularPoints(S2Point const& center,
                                           S1Angle radius,
                                           int num_vertices);

  // Append the vertices of "loop" to "vertices".
  static void AppendLoopVertices(S2Loop const& loop,
                                 std::vector<S2Point>* vertices);

  // A simple class that generates "Koch snowflake" fractals (see Wikipedia
  // for an introduction).  There is an option to control the fractal
  // dimension (between 1.0 and 2.0); values between 1.02 and 1.50 are
  // reasonable simulations of various coastlines.  The default dimension
  // (about 1.26) corresponds to the standard Koch snowflake.  (The west coast
  // of Britain has a fractal dimension of approximately 1.25.)
  //
  // The fractal is obtained by starting with an equilateral triangle and
  // recursively subdividing each edge into four segments of equal length.
  // Therefore the shape at level "n" consists of 3*(4**n) edges.  Multi-level
  // fractals are also supported: if you set min_level() to a non-negative
  // value, then the recursive subdivision has an equal probability of
  // stopping at any of the levels between the given min and max (inclusive).
  // This yields a fractal where the perimeter of the original triangle is
  // approximately equally divided between fractals at the various possible
  // levels.  If there are k distinct levels {min,..,max}, the expected number
  // of edges at each level "i" is approximately 3*(4**i)/k.
  class Fractal {
   public:
    // You must call set_max_level() or SetLevelForApproxMaxEdges() before
    // calling MakeLoop().
    Fractal();

    // Set the maximum subdivision level for the fractal (see above).
    // REQUIRES: max_level >= 0
    void set_max_level(int max_level);
    int max_level() const { return max_level_; }

    // Set the minimum subdivision level for the fractal (see above).  The
    // default value of -1 causes the min and max levels to be the same.  A
    // min_level of 0 should be avoided since this creates a significant
    // chance that none of the three original edges will be subdivided at all.
    void set_min_level(int min_level_arg);
    int min_level() const { return min_level_arg_; }

    // Set the min and/or max level to produce approximately the given number
    // of edges.  (The values are rounded to a nearby value of 3*(4**n).)
    void SetLevelForApproxMinEdges(int min_edges);
    void SetLevelForApproxMaxEdges(int max_edges);

    // Set the fractal dimension.  The default value of approximately 1.26
    // corresponds to the stardard Koch curve.  The value must lie in the
    // range [1.0, 2.0).
    void set_fractal_dimension(double dimension);
    double fractal_dimension() const { return dimension_; }

    // Return a lower bound on ratio (Rmin / R), where "R" is the radius
    // passed to MakeLoop() and "Rmin" is the minimum distance from the
    // fractal boundary to its center.  This can be used to inscribe another
    // geometric figure within the fractal without intersection.
    double min_radius_factor() const;

    // Return the ratio (Rmax / R), where "R" is the radius passed
    // to MakeLoop() and "Rmax" is the maximum distance from the fractal
    // boundary to its center.  This can be used to inscribe the fractal
    // within some other geometric figure without intersection.
    double max_radius_factor() const;

    // Return a fractal loop centered around the a-axis of the given
    // coordinate frame, with the first vertex in the direction of the
    // positive x-axis, and the given nominal radius.
    S2Loop* MakeLoop(Matrix3x3_d const& frame, S1Angle nominal_radius) const;

   private:
    void ComputeMinLevel();
    void ComputeOffsets();
    void GetR2Vertices(std::vector<R2Point>* vertices) const;
    void GetR2VerticesHelper(R2Point const& v0, R2Point const& v4, int level,
                             std::vector<R2Point>* vertices) const;

    int max_level_;
    int min_level_arg_;  // Value set by user
    int min_level_;      // Actual min level (depends on max_level_)
    double dimension_;

    // The ratio of the sub-edge length to the original edge length at each
    // subdivision step.
    double edge_fraction_;

    // The distance from the original edge to the middle vertex at each
    // subdivision step, as a fraction of the original edge length.
    double offset_fraction_;

    DISALLOW_COPY_AND_ASSIGN(Fractal);
  };

  // Convert a distance on the Earth's surface to an angle.
  static S1Angle KmToAngle(double km);
  static S1Angle MetersToAngle(double meters);

  // The Earth's mean radius in kilometers (according to NASA).
  static double const kEarthRadiusKm;

  // A deterministically-seeded random number generator.
  class Random;

  static Random rnd;

  // Return a random unit-length vector.
  static S2Point RandomPoint();

  // Return a right-handed coordinate frame (three orthonormal vectors).
  static void GetRandomFrame(S2Point* x, S2Point* y, S2Point* z);
  static Matrix3x3_d GetRandomFrame();

  // Given a unit-length z-axis, compute x- and y-axes such that (x,y,z) is a
  // right-handed coordinate frame (three orthonormal vectors).
  static void GetRandomFrameAt(S2Point const& z, S2Point* x, S2Point *y);
  static Matrix3x3_d GetRandomFrameAt(S2Point const& z);

  // Return a cap with a random axis such that the log of its area is
  // uniformly distributed between the logs of the two given values.
  // (The log of the cap angle is also approximately uniformly distributed.)
  static S2Cap GetRandomCap(double min_area, double max_area);

  // Return a point chosen uniformly at random (with respect to area)
  // from the given cap.
  static S2Point SamplePoint(S2Cap const& cap);

  // Return a point chosen uniformly at random (with respect to area on the
  // sphere) from the given latitude-longitude rectangle.
  static S2Point SamplePoint(S2LatLngRect const& rect);

  // Return a random cell id at the given level or at a randomly chosen
  // level.  The distribution is uniform over the space of cell ids,
  // but only approximately uniform over the surface of the sphere.
  static S2CellId GetRandomCellId(int level);
  static S2CellId GetRandomCellId();

  // Return a polygon with the specified center, number of concentric loops
  // and vertices per loop.
  static void ConcentricLoopsPolygon(S2Point const& center,
                                     int num_loops,
                                     int num_vertices_per_loop,
                                     S2Polygon* polygon);

  // Checks that "covering" completely covers the given region.  If
  // "check_tight" is true, also checks that it does not contain any cells
  // that do not intersect the given region.  ("id" is only used internally.)
  static void CheckCovering(S2Region const& region,
                            S2CellUnion const& covering,
                            bool check_tight,
                            S2CellId id = S2CellId());

  // Returns the user time consumed by this process, in seconds.
  static double GetCpuTime();

  // Delete all loops in the vector, and clear the vector.
  static void DeleteLoops(std::vector<S2Loop*>* loops);


 private:
  DISALLOW_IMPLICIT_CONSTRUCTORS(S2Testing);  // Contains static methods
};

// Functions in this class return random numbers that are as good as random()
// is.  The results are reproducible since the seed is deterministic.  This
// class is *NOT* thread-safe; it is only intended for testing purposes.
class S2Testing::Random {
 public:
  // Initialize using a deterministic seed.
  Random();

  // Reset the generator state using the given seed.
  void Reset(int32 seed);

  // Return a uniformly distributed 64-bit unsigned integer.
  uint64 Rand64();

  // Return a uniformly distributed 32-bit unsigned integer.
  uint32 Rand32();

  // Return a uniformly distributed "double" in the range [0,1).  Note that
  // the values returned are all multiples of 2**-53, which means that not all
  // possible values in this range are returned.
  double RandDouble();

  // Return a uniformly distributed integer in the range [0,n).
  int32 Uniform(int32 n);

  // Return a uniformly distributed "double" in the range [min, limit).
  double UniformDouble(double min, double limit);

  // A functor-style version of Uniform, so that this class can be used with
  // STL functions that require a RandomNumberGenerator concept.
  int32 operator() (int32 n) {
    return Uniform(n);
  }

  // Return true with probability 1 in n.
  bool OneIn(int32 n);

  // Skewed: pick "base" uniformly from range [0,max_log] and then
  // return "base" random bits.  The effect is to pick a number in the
  // range [0,2^max_log-1] with bias towards smaller numbers.
  int32 Skewed(int max_log);

 private:
  // Currently this class is based on random(), therefore it makes no sense to
  // make a copy.
  DISALLOW_COPY_AND_ASSIGN(Random);
};

#endif  // S2GEOMETRY_S2TESTING_H_

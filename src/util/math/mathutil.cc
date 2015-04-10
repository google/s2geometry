// Copyright 2008 Google Inc. All Rights Reserved.
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

#include "util/math/mathutil.h"

#include <stdlib.h>
#include <vector>

#include "base/integral_types.h"
#include <glog/logging.h>

MathUtil::QuadraticRootType MathUtil::DegenerateQuadraticRoots(
    long double b,
    long double c,
    long double *r1,
    long double *r2) {
  // This degenerate quadratic is really a linear equation b * x = -c.
  if (b == 0.0) {
    // The equation is constant, c == 0.
    if (c == 0.0) {
      // Quadratic equation is 0==0; treat as ambiguous, as if a==epsilon.
      *r1 = *r2 = 0.0;
      return kAmbiguous;
    }
    return kNoRealRoots;
  }
  // The linear equation has a single root at x = -c / b, not a double
  // one.  Respond as if a==epsilon: The other root is at "infinity",
  // which we signal with HUGE_VAL so that the behavior stays consistent
  // as a->0.
  *r1 = -c / b;
  *r2 = HUGE_VAL;
  return kTwoRealRoots;
}

bool MathUtil::RealRootsForCubic(long double const a,
                                 long double const b,
                                 long double const c,
                                 long double *const r1,
                                 long double *const r2,
                                 long double *const r3) {
  // According to Numerical Recipes (pp. 184-5), what
  // follows is an arrangement of computations to
  // compute the roots of a cubic that minimizes
  // roundoff error (as pointed out by A.J. Glassman).

  long double const a_squared = a * a, a_third = a / 3.0, b_tripled = 3.0 * b;
  long double const Q = (a_squared - b_tripled) / 9.0;
  long double const R =
      (2.0 * a_squared * a - 3.0 * a * b_tripled + 27.0 * c) / 54.0;

  long double const R_squared = R * R;
  long double const Q_cubed = Q * Q * Q;

  if (R_squared < Q_cubed) {
    long double const root_Q = sqrt(Q);
    long double const two_pi_third = 2.0 * M_PI / 3.0;
    long double const theta_third = acos(R / sqrt(Q_cubed)) / 3.0;
    long double const minus_two_root_Q = -2.0 * root_Q;

    *r1 = minus_two_root_Q * cos(theta_third) - a_third;
    *r2 = minus_two_root_Q * cos(theta_third + two_pi_third) - a_third;
    *r3 = minus_two_root_Q * cos(theta_third - two_pi_third) - a_third;

    return true;
  }

  long double const A =
    -sgn(R) * pow(std::abs(R) + sqrt(R_squared - Q_cubed), 1.0 / 3.0L);

  if (A != 0.0) {  // in which case, B from NR is zero
    *r1 = A + Q / A - a_third;
    return false;
  }

  *r1 = *r2 = *r3 = -a_third;
  return true;
}

namespace {

inline bool IsCloseToZero(long double x) {
  static const long double kEqnEPS = 1e-9;
  return x > -kEqnEPS && x < kEqnEPS;
}

}  // namespace

int MathUtil::RealRootsForQuartic(long double a,
                                  long double b,
                                  long double c,
                                  long double d,
                                  long double *roots) {
  // This is converted from code in Graphics Gems 1,
  // Schwarze, Jochen, Cubic and Quartic Roots, p. 404-407, code: p. 738-786.
  int number_of_roots = 0;

  // Converts to normal form: x^4 + Ax^3 + Bx^2 + Cx + D = 0.

  // Substitutes x = y - a/4 to eliminate cubic term: x^4 + px^2 + qx + r = 0.

  const long double sq_A = a * a;
  const long double p = -0.375 * sq_A + b;
  const long double q = 0.125 * sq_A * a - 0.5 * a * b + c;
  const long double r = -0.01171875 * sq_A * sq_A +
      0.0625 * sq_A * b - 0.25 * a * c + d;

  if (IsCloseToZero(r)) {
    // Handles case where there is no absolute term: y(y^3 + py + q) = 0.
    if (RealRootsForCubic(0, p, q, roots, roots + 1, roots + 2)) {
      number_of_roots = 3;
    } else {
      number_of_roots = 1;
    }
  } else {
    RealRootsForCubic(-0.5 * p,
                      -r,
                      0.5 * r * p - 0.125 * q * q,
                      roots,
                      roots + 1,
                      roots + 2);

    // Take one of the real roots...if there are multiple roots it doesn't seem
    // to matter which one is taken.
    const long double z = roots[0];

    long double u = z * z - r;

    if (IsCloseToZero(u))
      u = 0.0;
    else if (u > 0)
      u = sqrt(u);
    else
      return 0;

    long double v = 2.0 * z - p;

    if (IsCloseToZero(v))
      v = 0.0;
    else if (v > 0)
      v = sqrt(v);
    else
      return 0;

    number_of_roots = static_cast<int>(RealRootsForQuadratic(1.0,
                                                             q < 0 ? -v : v,
                                                             z - u,
                                                             roots,
                                                             roots + 1));

    number_of_roots +=
        static_cast<int>(RealRootsForQuadratic(1.0,
                                               q < 0 ? v : -v,
                                               z + u,
                                               roots + number_of_roots,
                                               roots + number_of_roots + 1));
  }

  // Resubstitutes into original equation.
  const long double sub = 0.25 * a;

  for (int i = 0; i < number_of_roots; ++i) {
    roots[i] -= sub;
  }
  return number_of_roots;
}

// Returns the greatest common divisor of two unsigned integers x and y,
// and assigns a, and b such that a*x + b*y = gcd(x, y).
unsigned int MathUtil::ExtendedGCD(unsigned int x, unsigned int y,
                                   int* a, int* b) {
  *a = 1;
  *b = 0;
  int c = 0;
  int d = 1;
  // before and after each loop:
  // current_x == a * original_x + b * original_y
  // current_y == c * original_x + d * original_y
  while (y != 0) {
    // div() takes int parameters; there is no version that takes unsigned int
    div_t r = div(static_cast<int>(x), static_cast<int>(y));
    x = y;
    y = r.rem;

    int tmp = c;
    c = *a - r.quot * c;
    *a = tmp;

    tmp = d;
    d = *b - r.quot * d;
    *b = tmp;
  }
  return x;
}


void MathUtil::ShardsToRead(const std::vector<bool>& shards_to_write,
                            std::vector<bool>* shards_to_read) {
  const int N = shards_to_read->size();
  const int M = shards_to_write.size();
  CHECK(N > 0 || M == 0) << ": have shards to write but not to read";

  // Input shard n of N can contribute to output shard m of M if there
  // exists a record with sharding hash x s.t. n = x % N and m = x % M.
  // Equivalently, there must exist s and t s.t. x = tN + n = sM + m,
  // i.e., tN - sM = m - n.  Since G = gcd(N, M) evenly divides tN - sM,
  // G must also evenly divide m - n.  Proof in the other direction is
  // left as an exercise.
  // Given output shard m, we should, therefore, read input shards n
  // that satisfy (n - m) = kG, i.e., n = m + kG.  Let 0 <= n < N.
  // Then, 0 <= m + kG < N and, finally, -m / G <= k < (N - m) / G.

  const int G = GCD(N, M);
  shards_to_read->assign(N, false);
  for (int m = 0; m < M; m++) {
    if (!shards_to_write[m]) continue;
    const int k_min = -m / G;
    const int k_max = k_min + N / G;
    for (int k = k_min; k < k_max; k++) {
      (*shards_to_read)[m + k * G] = true;
    }
  }
}

double MathUtil::Harmonic(int64 const n, double *const e) {
  CHECK_GT(n, 0);

  //   Hn ~ ln(n) + 0.5772156649 +
  //        + 1/(2n) - 1/(12n^2) + 1/(120n^4) - error,
  //   with 0 < error < 1/(256*n^4).

  double const
    d = static_cast<double>(n),
    d2 = d * d,
    d4 = d2 * d2;

  return (log(d) + 0.5772156649)  // ln + Gamma constant
    + 1 / (2 * d) - 1 / (12 * d2) + 1 / (120 * d4)
    - (*e = 1 / (256 * d4));
}

// The formula is extracted from the following page
// http://en.wikipedia.org/w/index.php?title=Stirling%27s_approximation
double MathUtil::Stirling(double n) {
  static const double kLog2Pi = log(2 * M_PI);
  const double logN = log(n);
  return (n * logN
          - n
          + 0.5 * (kLog2Pi + logN)      // 0.5 * log(2 * M_PI * n)
          + 1 / (12 * n)
          - 1 / (360 * n * n * n));
}

double MathUtil::LogCombinations(int n, int k) {
  CHECK_GE(n, k);
  CHECK_GT(n, 0);
  CHECK_GE(k, 0);

  // use symmetry to pick the shorter calculation
  if (k > n / 2) {
    k = n - k;
  }

  // If we have more than 30 logarithms to calculate, we'll use
  // Stirling's approximation for log(n!).
  if (k > 15) {
    return Stirling(n) - Stirling(k) - Stirling(n - k);
  } else {
    double result = 0;
    for (int i = 1; i <= k; i++) {
      result += log(static_cast<double>(n - k + i)) -
          log(static_cast<double>(i));
    }
    return result;
  }
}
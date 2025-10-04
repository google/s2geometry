// Copyright 2009 Google Inc. All Rights Reserved.
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

#include "s2/util/math/exactfloat/exactfloat.h"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <vector>

#include <gtest/gtest.h>
#include "absl/base/casts.h"
#include "absl/base/macros.h"
#include "absl/log/absl_log.h"
#include "absl/random/random.h"
#include "s2/util/math/vector.h"

namespace {

// To check the handling of boundary cases (NaN, infinities, out of range,
// overflow, etc), we check the ExactFloat versions of all the operators and
// functions against the glibc and/or hardware versions.
//
// The results are not always identical for several reasons.  The main reasons
// are that ExactFloat has a much higher precision, a much larger exponent
// range, and does not round.  But some differences are due to bugs in glibc.
// The worst of these are fixed below.

#ifndef DONT_FIX_GLIBC_MATH_FUNCTIONS

// There are a few cases where the glibc math functions return an incorrect
// result.  We handle this by defining our own corrected versions here.
//
// Note that these corrections are only known to be sufficient on exactly one
// platform, namely "--march=k8 --mfpmath=sse" with glibc 2.3.6.  Glibc
// undoubtedly has other bugs in other platforms and versions which would need
// to be corrected separately.

double fdim(double a, double b) {
  // ::fdim(-inf, -inf) incorrectly returns NaN.  The result should be zero
  // because fdim is defined to return 0 when a <= b, and IEEE 754 defines
  // -inf to be equal to itself (even though (-inf) - (-inf) is NaN).
  //
  // Note that ::fdim(inf, inf) does return 0 (the correct result).
  if (a == b) return 0;
  return ::fdim(a, b);
}

double fmax(double a, double b) {
  // fmax(0, -0) returns -0 rather than +0.  This is technically not a bug,
  // but it is different from (and inferior to) the MPFloat behavior.
  if (a == 0 && b == 0) return a + b;
  return ::fmax(a, b);
}

double fmin(double a, double b) {
  // fmin(-0, 0) returns +0 rather than -0.  This is technically not a bug,
  // but it is different from (and inferior to) the MPFloat behavior.
  if (a == 0 && b == 0 && (std::signbit(a) || std::signbit(b))) {
    return copysign(0, -1);
  }
  return ::fmin(a, b);
}

double logb(double a) {
  // If "a" is denormalized, logb() is supposed to return the exponent that
  // "a" would have if it were normalized.  (But it doesn't.)
  if (a != 0 && fabs(a) < std::numeric_limits<double>::min()) {
    return ::logb(scalbn(a, 100)) - 100;
  }
  return ::logb(a);
}

double ldexp(double a, int exp) {
  // ldexp() incorrectly returns infinity rather than zero when the
  // "a" is finite and "exp" is a very large negative value.
  double r = ::ldexp(a, exp);
  if (std::isinf(r) && !std::isinf(a) && exp < 0) {
    return copysign(0, r);
  }
  return r;
}

double scalbn(double a, int exp) {
  return ldexp(a, exp);  // See ldexp().
}

double scalbln(double a, long exp) {
  // Clamp the exponent to the range of "int" in order to avoid truncation.
  exp = std::max(static_cast<long>(INT_MIN),
                 std::min(static_cast<long>(INT_MAX), exp));
  return ldexp(a, exp);  // See ldexp().
}

// The C standard does not specify the result of converting a floating-point
// number to an integer if the argument is NaN or out of range.  This applies
// to both static_cast<int_type> and functions such as lrint, lround, etc.
//
// With the glibc/Intel platform tested here, these conversions return the
// minimum possible value of the result type.  (Except if the argument is a
// constant, in which case gcc converts NaN to zero and clamps out of range
// arguments to the minimum or maximum possible value.)
//
// Here we fix the rounding functions to match MPFloat, which clamps out of
// range values and returns the maximum possible value for NaN.

#define FIX_INT_ROUNDING(T, fname)                      \
  T fname(double a) {                                   \
    if (std::isnan(a)) return std::numeric_limits<T>::max(); \
    if (a <= std::numeric_limits<T>::min())             \
      return std::numeric_limits<T>::min();             \
    if (a >= std::numeric_limits<T>::max())             \
      return std::numeric_limits<T>::max();             \
    return ::fname(a);                                  \
  }

FIX_INT_ROUNDING(long, lrint)
FIX_INT_ROUNDING(long, lround)
FIX_INT_ROUNDING(long long, llrint)
FIX_INT_ROUNDING(long long, llround)

#endif  // DONT_FIX_GLIBC_MATH_FUNCTIONS

// A list of double-precision constants to use as arguments when evaluating
// math intrinsics.  The negated values of these constants are used as well,
// so only one constant with a given absolute value needs to be listed.

const double kSpecialUnsignedDoubleValues[] = {
    std::numeric_limits<double>::quiet_NaN(),

    // Minimum and maximum values of various relevant C++ types.
    std::numeric_limits<double>::infinity(),
    std::numeric_limits<double>::max(),
    std::numeric_limits<double>::min(),
    std::numeric_limits<double>::denorm_min(),
    static_cast<double>(std::numeric_limits<int>::max()),
    static_cast<double>(std::numeric_limits<int>::min()),
    static_cast<double>(std::numeric_limits<long long>::max()),
    static_cast<double>(std::numeric_limits<long long>::min()),
    static_cast<double>(std::numeric_limits<unsigned long long>::max()),
    0,

    // Small and large values that are not quite minimums or maximums.
    1.23e-300,
    1e-20,
    1e20,
    1.23e300,

    // Interesting values for trigonometric, exponential, and logarithm
    // functions.
    2 * M_PI,
    M_PI,
    M_PI_2,
    M_PI_4,
    M_E,
    M_LN2,

    // Positive and negative powers of two.
    1 << 30,
    256,
    16,
    2,
    1,
    0.5,
    1.0 / 256,
    1.0 / (1 << 30),

    // Interesting values for integer rounding functions.
    1.5,
    2.5,
    3.5,

    // Miscellaneous integer and non-integer "ordinary" values.
    42,
    7,
    5,
    0.3,
    0.17,
};

class ExactFloatTest : public ::testing::Test {
 public:
  // Initialize the list of constants to be used for testing intrinsics.
  static void SetUpTestSuite() {
    for (int i = 0; i < ABSL_ARRAYSIZE(kSpecialUnsignedDoubleValues); ++i) {
      double d = kSpecialUnsignedDoubleValues[i];
      kSpecialDoubleValues.push_back(d);
      // Glibc and MPFloat handle negative NaN values differently.  To avoid
      // discrepancies, we only test positively-signed NaN values.
      if (!std::isnan(d)) {
        kSpecialDoubleValues.push_back(-d);
      }
    }
  }

  static void TearDownTestSuite() { kSpecialDoubleValues.clear(); }

  // Return the difference measured in ulps (units in the last place) between
  // two floating-point values.  Return 0 if the values are equal or both are
  // NaN, and return the largest possible uint64_t if exactly one value is NaN.
  // Note that +0 and -0 are equal (i.e., they differ by 0 ulps), and that the
  // smallest positive and negative values differ by 2 ulps.  Infinity is one
  // ulp larger than the largest finite number.
  static uint64_t GetErrorUlps(double a, double b) {
    if (std::isnan(a) && std::isnan(b)) return 0;
    if (std::isnan(a) || std::isnan(b)) {
      return std::numeric_limits<uint64_t>::max();
    }

    // Floating-point numbers are arranged so that for numbers of the same
    // sign, the difference in ulps is just the difference between the two
    // numbers viewed as 64-bit unsigned integers.
    uint64_t a_bits = absl::bit_cast<uint64_t>(a);
    uint64_t b_bits = absl::bit_cast<uint64_t>(b);
    if (std::signbit(a) == std::signbit(b)) {
      return (a_bits > b_bits) ? (a_bits - b_bits) : (b_bits - a_bits);
    }
    // For numbers of opposite sign, we take the difference in ulps between
    // each number and the zero of the same sign, and add them together.
    a_bits ^= static_cast<uint64_t>(1) << 63;
    return a_bits + b_bits;
  }

  // Return true if the difference between "expected" and "actual" is at most
  // the given number of ulps (units in the last place).  The two values are
  // also required to have the same sign bit unless they are both NaN.  (So for
  // example, +0 and -0 are not equivalent.)
  static bool IsExpected(double expected, double actual, uint64_t ulps) {
    // We require the sign bit to match unless the values are NaN.
    if (!std::isnan(expected) &&
        std::signbit(expected) != std::signbit(actual)) {
      return false;
    }
    return GetErrorUlps(expected, actual) <= ulps;
  }

  // Expect "actual" to have the given value when converted to a "double".
  // Two values are considered equivalent if they have the same bit pattern or
  // they are both NaN.  (So for example, +0 and -0 are not equivalent.)
  void ExpectSame(double expected, const ExactFloat &xf_actual) {
    double actual = xf_actual.ToDouble();
    if (std::isnan(expected)) {
      EXPECT_TRUE(std::isnan(actual));
    } else {
      // Keep the ugly std::signbit() macro out of the error messages.
      bool expected_sign = std::signbit(expected);
      bool actual_sign = std::signbit(actual);
      EXPECT_EQ(expected_sign, actual_sign);
      EXPECT_EQ(expected, actual);
    }
  }

  // Like ExpectSame() but also check that "actual" has the expected precision.
  void ExpectSameWithPrec(double expected_value, int expected_prec,
                          const ExactFloat &xf_actual) {
    ExpectSame(expected_value, xf_actual);
    EXPECT_EQ(expected_prec, xf_actual.prec());
  }

  // Log an error when a math intrinsic does not return the expected result.
  static void AddMathcallFailure(const testing::Message &call_msg,
                                 double expected, double actual) {
    ADD_FAILURE() << call_msg << "\nExpected (glibc): " << ExactFloat(expected)
                  << "\nActual (ExactFloat): " << ExactFloat(actual)
                  << "\nError: " << GetErrorUlps(expected, actual) << " ulps";
  }

 protected:
  // Given two versions "f" and "mp_f" of the unary function called "fname",
  // check that their results agree to within the given number of ulps for a
  // range of test arguments.
  void TestMathcall1(const char *fname, double f(double),
                     ExactFloat mp_f(const ExactFloat &), uint64_t ulps) {
    for (int i = 0; i < kSpecialDoubleValues.size(); ++i) {
      double a = kSpecialDoubleValues[i];
      double expected = f(a);
      double actual = mp_f(ExactFloat(a)).ToDouble();
      if (!IsExpected(expected, actual, ulps)) {
        AddMathcallFailure(testing::Message()
                               << fname << "(" << ExactFloat(a) << ")",
                           expected, actual);
      }
    }
  }

  // Given two versions "f" and "mp_f" of the binary function called "fname",
  // check that their results agree to within the given number of ulps for a
  // range of test arguments.
  void TestMathcall2(const char *fname, double f(double, double),
                     ExactFloat mp_f(const ExactFloat &, const ExactFloat &),
                     uint64_t ulps) {
    for (int i = 0; i < kSpecialDoubleValues.size(); ++i) {
      double a = kSpecialDoubleValues[i];
      for (int j = 0; j < kSpecialDoubleValues.size(); ++j) {
        double b = kSpecialDoubleValues[j];
        double expected = f(a, b);
        double actual = mp_f(ExactFloat(a), ExactFloat(b)).ToDouble();
        if (!IsExpected(expected, actual, ulps)) {
          AddMathcallFailure(testing::Message() << fname << "(" << ExactFloat(a)
                                                << ", " << ExactFloat(b) << ")",
                             expected, actual);
        }
      }
    }
  }

  // Given a one-argument function "f" and a zero-argument ExactFloat member
  // function "mp_f", check that they return the same result on a range of
  // test arguments.
  template <typename ResultType>
  void TestMethod0(const char *fname, ResultType f(double),
                   ResultType (ExactFloat::*mp_f)() const) {
    for (int i = 0; i < kSpecialDoubleValues.size(); ++i) {
      double a = kSpecialDoubleValues[i];
      ResultType expected = f(a);
      ResultType actual = (ExactFloat(a).*mp_f)();
      if (expected != actual) {
        AddMathcallFailure(testing::Message()
                               << ExactFloat(a) << "." << fname << "()",
                           expected, actual);
      }
    }
  }

  // Given two versions "f" and "mp_f" of a unary function called "fname" that
  // returns the integer type ResultType, check that they return the same
  // value for a range of test arguments.
  template <typename ResultType>
  void TestIntMathcall1(const char *fname, ResultType f(double),
                        ResultType mp_f(const ExactFloat &)) {
    for (int i = 0; i < kSpecialDoubleValues.size(); ++i) {
      double a = kSpecialDoubleValues[i];
      ResultType expected = f(a);
      ResultType actual = mp_f(ExactFloat(a));
      if (actual != expected) {
        AddMathcallFailure(
            testing::Message() << fname << "(" << ExactFloat(a) << ")",
            static_cast<double>(expected), static_cast<double>(actual));
      }
    }
  }

  // Given two versions "f" and "mp_f" of a function called "fname" that has
  // an interface similar to ldexp() (i.e., with one "double" argument and one
  // integer argument), check that they return the same result for a range of
  // test arguments.  "ExpType" is the type of the integer argument.
  template <typename ExpType>
  void TestLdexpCall(const char *fname, double f(double, ExpType),
                     ExactFloat mp_f(const ExactFloat &, ExpType)) {
    static const ExpType kUnsignedExpValues[] = {
        // Doesn't test with numeric_limits<ExpType>::min() because it's
        // undefined
        // to negate the min value of a signed number.
        std::numeric_limits<ExpType>::max(),
        0,
        1,
        2,
        3,
        10,
        20,
        50,
        100,
        200,
        500,
        999,
        1022,
        1023,
        1024,
        1073,
        1074,
        2046,
        2047,
        2048,
        4096,
        1000000,
    };
    for (int i = 0; i < kSpecialDoubleValues.size(); ++i) {
      double a = kSpecialDoubleValues[i];
      for (int j = 0; j < ABSL_ARRAYSIZE(kUnsignedExpValues); ++j) {
        for (int sign = -1; sign <= +1; sign += 2) {
          ExpType exp = sign * kUnsignedExpValues[j];
          double expected = f(a, exp);
          double actual = mp_f(ExactFloat(a), exp).ToDouble();
          if (!IsExpected(expected, actual, 0)) {
            AddMathcallFailure(testing::Message()
                                   << fname << "(" << ExactFloat(a) << ", "
                                   << exp << ")",
                               expected, actual);
          }
        }
      }
    }
  }

 private:
  static std::vector<double> kSpecialDoubleValues;
};

std::vector<double> ExactFloatTest::kSpecialDoubleValues;

TEST_F(ExactFloatTest, GetErrorUlps) {
  // Verify some of the assertions made by GetErrorUlps().
  EXPECT_EQ(0, GetErrorUlps(NAN, -NAN));
  EXPECT_EQ(0, GetErrorUlps(0, copysign(0, -1)));
  EXPECT_EQ(2, GetErrorUlps(std::numeric_limits<double>::denorm_min(),
                            -std::numeric_limits<double>::denorm_min()));
  EXPECT_EQ(1, GetErrorUlps(std::numeric_limits<double>::max(), INFINITY));
  EXPECT_EQ(std::numeric_limits<uint64_t>::max(), GetErrorUlps(0, NAN));
}

TEST_F(ExactFloatTest, Constructors) {
  // Default constructor.
  ExactFloat a;
  ExpectSameWithPrec(0.0, 0, a);

  // Construct an ExactFloat from the literal constant "0".  It is necessary to
  // have an implicit constructor that takes an "int" argument for this to
  // work; otherwise the compiler thinks that "0" is ambiguous (it could be a
  // "double" or a "const char*").
  ExactFloat a0 = 0;
  ExpectSameWithPrec(0.0, 0, a);

  // Implicit constructor from "double".
  ExactFloat b = 2.5;
  ExpectSameWithPrec(2.5, 3, b);

  // Implicit constructor from "int".
  ExactFloat c = -125;
  ExpectSameWithPrec(-125, 7, c);

  // Copy constructor.
  ExactFloat e = c;
  ExpectSameWithPrec(-125, 7, e);

  // Ensure that construction with INT_MIN works properly.
  ExactFloat f = INT_MIN;
  ExpectSame(INT_MIN, f);
}

TEST_F(ExactFloatTest, Constants) {
  EXPECT_TRUE(ExactFloat::SignedZero(+1).is_zero());
  EXPECT_EQ(false, ExactFloat::SignedZero(+1).sign_bit());
  EXPECT_TRUE(ExactFloat::SignedZero(-1).is_zero());
  EXPECT_EQ(true, ExactFloat::SignedZero(-1).sign_bit());

  EXPECT_TRUE(ExactFloat::Infinity(+1).is_inf());
  EXPECT_EQ(false, ExactFloat::Infinity(+1).sign_bit());
  EXPECT_TRUE(ExactFloat::Infinity(-1).is_inf());
  EXPECT_EQ(true, ExactFloat::Infinity(-1).sign_bit());

  EXPECT_TRUE(ExactFloat::NaN().is_nan());
}

TEST_F(ExactFloatTest, Accessors) {
  // A few quick checks of prec(), max_prec(), and exp().
  ExactFloat x = 7.25;
  ExpectSameWithPrec(7.25, 5, x);
  EXPECT_EQ(3, x.exp());
  EXPECT_EQ(ExactFloat::kMaxPrec, x.max_prec());

  x = 255;
  ExpectSameWithPrec(255, 8, x);
  EXPECT_EQ(8, x.exp());

  x = 256;
  ExpectSameWithPrec(256, 1, x);
  EXPECT_EQ(9, x.exp());
}

TEST_F(ExactFloatTest, set_zero) {
  ExactFloat x = 5.0;
  EXPECT_FALSE(x.is_zero());

  x.set_zero(+1);
  EXPECT_TRUE(x.is_zero());
  ExpectSameWithPrec(0.0, 0, x);

  x.set_zero(-1);
  EXPECT_TRUE(x.is_zero());
  ExpectSameWithPrec(copysign(0.0, -1), 0, x);
}

TEST_F(ExactFloatTest, set_inf) {
  ExactFloat x;
  EXPECT_FALSE(x.is_inf());

  x.set_inf(+1);
  EXPECT_TRUE(x.is_inf());
  ExpectSameWithPrec(INFINITY, 0, x);

  x.set_inf(-1);
  EXPECT_TRUE(x.is_inf());
  ExpectSameWithPrec(-INFINITY, 0, x);

  x = NAN;
  EXPECT_FALSE(x.is_inf());

  x = INFINITY;
  EXPECT_TRUE(x.is_inf());
  ExpectSameWithPrec(INFINITY, 0, x);

  x = -INFINITY;
  EXPECT_TRUE(x.is_inf());
  ExpectSameWithPrec(-INFINITY, 0, x);
}

TEST_F(ExactFloatTest, set_nan) {
  ExactFloat x;
  EXPECT_FALSE(x.is_nan());
  ExpectSameWithPrec(0.0, 0, x);

  x.set_nan();
  EXPECT_TRUE(x.is_nan());
  ExpectSameWithPrec(NAN, 0, x);

  x = 2;
  EXPECT_FALSE(x.is_nan());

  x = NAN;
  EXPECT_TRUE(x.is_nan());
  ExpectSameWithPrec(NAN, 0, x);
}

TEST_F(ExactFloatTest, ToDouble) {
  EXPECT_EQ(0, ExactFloat(0).ToDouble());
  ExpectSame(copysign(0.0, -1), -ExactFloat(0));
  EXPECT_EQ(std::numeric_limits<double>::max(),
            ExactFloat(std::numeric_limits<double>::max()).ToDouble());
  EXPECT_EQ(-std::numeric_limits<double>::min(),
            ExactFloat(-std::numeric_limits<double>::min()).ToDouble());
  EXPECT_EQ(std::numeric_limits<double>::denorm_min(),
            ExactFloat(std::numeric_limits<double>::denorm_min()).ToDouble());
  EXPECT_EQ(-12.7, ExactFloat(-12.7).ToDouble());
  EXPECT_EQ(M_PI, ExactFloat(M_PI).ToDouble());
  EXPECT_EQ(INFINITY, ExactFloat(INFINITY).ToDouble());
  EXPECT_EQ(-INFINITY, ExactFloat(-INFINITY).ToDouble());
  EXPECT_TRUE(std::isnan(ExactFloat(NAN).ToDouble()));
}

TEST_F(ExactFloatTest, ToString) {
  EXPECT_EQ("0.001", ExactFloat(0.001).ToString());
  EXPECT_EQ("0.10000000000000001", ExactFloat(0.1).ToString());
  EXPECT_EQ("0.001953125", ExactFloat(1. / 512).ToString());
  EXPECT_EQ("1e-50", ExactFloat(1e-50).ToString());
  EXPECT_EQ("1e-50", ExactFloat(1e-50).ToString());
  EXPECT_EQ("0", ExactFloat(0).ToString());
  EXPECT_EQ("nan", ExactFloat(NAN).ToString());
  EXPECT_EQ("-inf", ExactFloat(-INFINITY).ToString());
  EXPECT_EQ("65536", ExactFloat(65536).ToString());
  EXPECT_EQ("1048576", ExactFloat(1 << 20).ToString());
  EXPECT_EQ("1048575", ExactFloat((1 << 20) - 1).ToString());
  EXPECT_EQ("1073741824", ExactFloat(1 << 30).ToString());
  EXPECT_EQ("1073741823", ExactFloat((1 << 30) - 1).ToString());
  EXPECT_EQ("1.099511628e+12", ExactFloat(pow(2, 40)).ToString());
  EXPECT_EQ("1099511627775", ExactFloat(pow(2, 40) - 1).ToString());
  EXPECT_EQ("1.606938044e+60", ExactFloat(pow(2, 200)).ToString());
  EXPECT_EQ("1606938044258990275541962092341162602522202993782792835301375",
            (ExactFloat(pow(2, 200)) - 1).ToString());
  EXPECT_EQ("3.1415926535897931", ExactFloat(M_PI).ToString());

  // Verify that ToString implements kRoundTiesToEven (the same as glibc).
  EXPECT_EQ("0.62", ExactFloat(0.625).ToStringWithMaxDigits(2));
  EXPECT_EQ("0.62", ExactFloat(0.625).ToStringWithMaxDigits(2));
  EXPECT_EQ("1.88", ExactFloat(1.875).ToStringWithMaxDigits(3));

  // Three digits are sufficient to distinguish any two numbers with 6-bit
  // mantissas, but up to four digits are required to distinguish two numbers
  // with 7-bit mantissas (e.g. 1016 and 1024).
  EXPECT_EQ(3, ExactFloat::NumSignificantDigitsForPrec(6));
  EXPECT_EQ(4, ExactFloat::NumSignificantDigitsForPrec(7));

  ABSL_LOG(INFO) << "Test of ExactFloat stream output: " << ExactFloat(1.23456);
}

TEST_F(ExactFloatTest, ToUniqueString) {
  // The following test demonstrates that ExactFloats with different values and
  // different precisions can generate the same ToString() value.
  ExactFloat value = 0.01;
  ExactFloat v1 = value.RoundToMaxPrec(33, ExactFloat::kRoundTiesToEven);
  ExactFloat v2 = value.RoundToMaxPrec(46, ExactFloat::kRoundTiesToEven);
  EXPECT_NE(v1, v2);
  EXPECT_EQ("0.01", v1.ToString());
  EXPECT_EQ("0.01", v2.ToString());
  EXPECT_EQ("0.01<33>", v1.ToUniqueString());
  EXPECT_EQ("0.01<46>", v2.ToUniqueString());
}

TEST_F(ExactFloatTest, RoundToMaxPrec) {
  // This method is used by the ToDouble() implementation and also the
  // rounding function such as trunc(), floor(), etc, so it is already pretty
  // well tested.  We do a few other tests here.
  EXPECT_EQ(ExactFloat(pow(2, 500)),
            (ExactFloat(pow(2, 500)) - 1)
                .RoundToMaxPrec(2, ExactFloat::kRoundTiesToEven));
  EXPECT_EQ(ExactFloat(pow(2, 500)),
            (ExactFloat(pow(2, 500)) - 1)
                .RoundToMaxPrec(499, ExactFloat::kRoundTiesToEven));
  EXPECT_EQ(ExactFloat(pow(2, 500)) - 4,
            (ExactFloat(pow(2, 500)) - 3)
                .RoundToMaxPrec(499, ExactFloat::kRoundTiesToEven));
  EXPECT_EQ(ExactFloat(pow(2, 500)) + 16,
            (ExactFloat(pow(2, 500)) + 1)
                .RoundToMaxPrec(497, ExactFloat::kRoundTowardPositive));
  EXPECT_EQ(2.7, ExactFloat(2.7).RoundToMaxPrec(1000000,
                                                ExactFloat::kRoundTiesToEven));
}

// Test the ExactFloat version of the unary operator "op" against the
// corresponding C++ operator.
#define TEST_MATHOP1(op_name, op)                                \
  double op_name(double a) { return op(a); }                     \
  ExactFloat mp_##op_name(const ExactFloat &a) { return op(a); } \
  TEST_F(ExactFloatTest, op_name) {                              \
    TestMathcall1(#op_name, op_name, mp_##op_name, 0);           \
  }

TEST_MATHOP1(plus, +)
TEST_MATHOP1(minus, -)

// Test the ExactFloat version of the binary operator "op" against the
// corresponding C++ operator.
#define TEST_MATHOP2(op_name, op)                                     \
  double op_name(double a, double b) { return (a)op(b); }             \
  ExactFloat mp_##op_name(const ExactFloat &a, const ExactFloat &b) { \
    return (a)op(b);                                                  \
  }                                                                   \
  TEST_F(ExactFloatTest, op_name) {                                   \
    TestMathcall2(#op_name, op_name, mp_##op_name, 0);                \
  }

TEST_MATHOP2(add, +);
TEST_MATHOP2(subtract, -);
TEST_MATHOP2(multiply, *);
// Division is not defined.
TEST_MATHOP2(equal, ==);
TEST_MATHOP2(less, <);
TEST_MATHOP2(greater, >);
TEST_MATHOP2(not_equal, !=);
TEST_MATHOP2(not_less, >=);
TEST_MATHOP2(not_greater, <=);

// Test the various assignment operators.
#define TEST_ASSIGNOP(op_name, op)                                    \
  double op_name(double a, double b) {                                \
    (a) op(b);                                                        \
    return a;                                                         \
  }                                                                   \
  ExactFloat mp_##op_name(const ExactFloat &a, const ExactFloat &b) { \
    ExactFloat x = a;                                                 \
    x op(b);                                                          \
    return x;                                                         \
  }                                                                   \
  TEST_F(ExactFloatTest, op_name) {                                   \
    TestMathcall2(#op_name, op_name, mp_##op_name, 0);                \
  }
TEST_ASSIGNOP(assignment, =);
TEST_ASSIGNOP(plus_equals, +=);
TEST_ASSIGNOP(minus_equals, -=);
TEST_ASSIGNOP(times_equals, *=);
// Division is not defined.

// Check that the ExactFloat and glibc versions of "func" always return the same
// value to within the given number of ulps.
#define TEST_MATHCALL1(func, ulps)                                       \
  /* We define a wrapper function around ExactFloat version of "func" */ \
  /* so that we can take its address (gcc can't find it otherwise). */   \
  ExactFloat mp_##func(const ExactFloat &a) { return func(a); }          \
  TEST_F(ExactFloatTest, func) { TestMathcall1(#func, func, mp_##func, ulps); }

// Test all the unary math instrinsics (in the same order as the .h file).

TEST_MATHCALL1(fabs, 0)
TEST_MATHCALL1(ceil, 0)
TEST_MATHCALL1(floor, 0)
TEST_MATHCALL1(trunc, 0)
TEST_MATHCALL1(rint, 0)
TEST_MATHCALL1(nearbyint, 0)
TEST_MATHCALL1(logb, 0)
TEST_MATHCALL1(round, 0)

// Check that the ExactFloat and glibc versions of "func" always return the
// same value to within the given number of ulps.
#define TEST_MATHCALL2(func, ulps)                                 \
  ExactFloat mp_##func(const ExactFloat &a, const ExactFloat &b) { \
    return func(a, b);                                             \
  }                                                                \
  TEST_F(ExactFloatTest, func) { TestMathcall2(#func, func, mp_##func, ulps); }

// Test all the binary math instrinsics (in the same order as the .h file).

TEST_MATHCALL2(fmax, 0)
TEST_MATHCALL2(fmin, 0)
TEST_MATHCALL2(fdim, 0)
TEST_MATHCALL2(copysign, 0)

// Check that the ExactFloat and glibc versions of "func" return the same
// integer value.

#define TEST_INTEGER_MATHCALL1(ResultType, func)                \
  ResultType mp_##func(const ExactFloat &a) { return func(a); } \
  TEST_F(ExactFloatTest, func) { TestIntMathcall1(#func, func, mp_##func); }

TEST_INTEGER_MATHCALL1(long, lrint);
TEST_INTEGER_MATHCALL1(long long, llrint);
TEST_INTEGER_MATHCALL1(long, lround);
TEST_INTEGER_MATHCALL1(long long, llround);
TEST_INTEGER_MATHCALL1(int, ilogb);

////////////////////////////////////////////////////////////////////////////
// frexp(): Test using wrapper functions to check each result separately.
double frexp_frac(double a) {
  int exp_part;
  return frexp(a, &exp_part);
}
int frexp_exp(double a) {
  int exp_part;
  (void)frexp(a, &exp_part);
  return exp_part;
}
ExactFloat mp_frexp_frac(const ExactFloat &a) {
  int exp_part;
  return frexp(a, &exp_part);
}
int mp_frexp_exp(const ExactFloat &a) {
  int exp_part;
  (void)frexp(a, &exp_part);
  return exp_part;
}

#if !defined(__wasm__)
TEST_F(ExactFloatTest, frexp) {
  TestMathcall1("frexp_frac", frexp_frac, mp_frexp_frac, 0);
  TestIntMathcall1("frexp_exp", frexp_exp, mp_frexp_exp);
}
#endif
////////////////////////////////////////////////////////////////////////////
// ldexp(), scalbn(), scalbln()

#define TEST_LDEXP_CALL(ExpType, func)                     \
  ExactFloat mp_##func(const ExactFloat &a, ExpType exp) { \
    return func(a, exp);                                   \
  }                                                        \
  TEST_F(ExactFloatTest, func) { TestLdexpCall(#func, func, mp_##func); }
// Results are different on Apple platforms; it's not clear if it's an Apple
// or glibc problem.
#if !defined(__APPLE__)
TEST_LDEXP_CALL(int, ldexp);
TEST_LDEXP_CALL(int, scalbn);
TEST_LDEXP_CALL(long, scalbln);
#endif
// Test a zero-argument ExactFloat member function against a corresponding
// one-argument reference function.
#define TEST_METHOD0_VS_FUNCTION(ResultType, method, function)       \
  TEST_F(ExactFloatTest, method) {                                   \
    TestMethod0<ResultType>(#method, function, &ExactFloat::method); \
  }

// "Reference versions" of various zero-argument methods.
bool ref_is_zero(double a) { return a == 0; }
bool ref_is_normal(double a) {
  return (a != 0) && !std::isinf(a) && !std::isnan(a);
}
bool ref_sign_bit(double a) { return std::signbit(a); }
int ref_sgn(double a) { return (a > 0) ? 1 : (a < 0) ? -1 : 0; }

TEST_METHOD0_VS_FUNCTION(bool, is_zero, ref_is_zero)
TEST_METHOD0_VS_FUNCTION(bool, is_normal, ref_is_normal)
TEST_METHOD0_VS_FUNCTION(bool, sign_bit, ref_sign_bit)
TEST_METHOD0_VS_FUNCTION(int, sgn, ref_sgn)

// Test a zero-argument ExactFloat member function against a corresponding
// one-argument std:: function.
#define TEST_METHOD0_VS_STD_FN(ResultType, method, fn) \
  ResultType ref_##fn(double a) { return std::fn(a); }  \
  TEST_METHOD0_VS_FUNCTION(ResultType, method, ref_##fn)

TEST_METHOD0_VS_STD_FN(bool, is_inf, isinf)
TEST_METHOD0_VS_STD_FN(bool, is_nan, isnan)
TEST_METHOD0_VS_STD_FN(bool, is_finite, isfinite)

typedef Vector2<ExactFloat> MyVec2;
typedef Vector3<ExactFloat> MyVec3;

TEST_F(ExactFloatTest, Vector3_Part1) {
  // Check that ExactFloat can be used with templatized classes such as Vector2
  // and Vector3.  We exercise all the methods of Vector3 to check that they
  // compile and do something reasonable.
  //
  // IsNaN() is not tested because it cannot be used with vectors of ExactFloats
  // (see the header file).
  //
  // This test is split into parts because of the stack frame size limit.

  MyVec3 zero;
  MyVec3 zero2(0, 0, 0);
  EXPECT_EQ(zero, zero2);
  MyVec3 zero3(MyVec2(0, 0), 0);
  EXPECT_EQ(zero, zero3);
  MyVec3 almost_zero(1e-200, 1e-100, 1e-50);
  EXPECT_NE(zero, almost_zero);
  EXPECT_TRUE(almost_zero.aequal(zero, 1e-15));
  EXPECT_LT(zero, almost_zero);  // Total order, for STL purposes.
  MyVec3 almost_zero2(ExactFloat(1e-320), 0, 0);
  EXPECT_NE(zero, almost_zero2);
  EXPECT_TRUE(almost_zero2.aequal(zero, 1e-300));
  EXPECT_LT(almost_zero2, almost_zero);
  EXPECT_LT(zero, almost_zero2);

  EXPECT_EQ(MyVec3(1, -2, 3), MyVec3::Cast(Vector3_d(1, -2, 3)));
  EXPECT_EQ(MyVec3(1, -2, 3), MyVec3::Cast(Vector3<ExactFloat>(1, -2, 3)));

  MyVec3 x(1.5, 2.5, 3.5);
  // Vector3_d y = Vector3_d::Cast(x);  // NOT ALLOWED
  Vector3_d y = Vector3_d(x[0].ToDouble(), x[1].ToDouble(), x[2].ToDouble());
  static_cast<void>(y);
}

TEST_F(ExactFloatTest, Vector3_Part2) {
  MyVec3 x;
  x.Set(1, 2, 3.5);
  EXPECT_EQ(MyVec3(1, 2, 3.5), x);
  x = MyVec3(4, 5, 6);
  EXPECT_EQ(MyVec3(4, 5, 6), x);
  x += MyVec3(1, 2, 3);
  EXPECT_EQ(MyVec3(5, 7, 9), x);
  x -= MyVec3(3, 2, 1);
  EXPECT_EQ(MyVec3(2, 5, 8), x);
  x *= 3;
  EXPECT_EQ(MyVec3(6, 15, 24), x);
  x = MyVec3(2, 5, 8);
  MyVec3 y(1, 2, 3);
  MyVec3 z = x.MulComponents(y);
  EXPECT_EQ(MyVec3(2, 10, 24), z);
  z = MyVec3(1, 2, 3);
  EXPECT_EQ(MyVec3(5, 7, 9), z + MyVec3(4, 5, 6));
  EXPECT_EQ(MyVec3(-5, -3, -1), z - MyVec3(6, 5, 4));
  EXPECT_EQ(MyVec3(-1, -2, -3), -z);
  EXPECT_EQ(z.DotProd(z), 14);
  EXPECT_EQ(MyVec3(2, 4, 6), 2 * z);
}

TEST_F(ExactFloatTest, Vector3_Part3) {
  MyVec3 y(3, 4, 1);
  MyVec3 z(2, -3, 1);
  EXPECT_EQ(MyVec3(7, -1, -17), y.CrossProd(z));
  EXPECT_EQ(3, y[0]);
  y[0] = 5;
  EXPECT_EQ(5, y.x());
  EXPECT_EQ(4, y.y());
  EXPECT_EQ(1, y.z());
  y.x(3);
  EXPECT_EQ(MyVec3(3, 4, 1), y);
  ExactFloat *y_data = y.Data();
  y_data[2] = 12;
  EXPECT_EQ(MyVec3(3, 4, 12), y);
  EXPECT_EQ(y.Data()[1], ExactFloat(4));
  EXPECT_EQ(y.Norm2(), 169);
  EXPECT_EQ(6, z.DotProd(y));
  EXPECT_EQ(2, y.LargestAbsComponent());
  EXPECT_EQ(Vector3<int>(1, 2, 0), z.ComponentOrder());
  EXPECT_EQ(MyVec3(4, -6, 2), 2 * z);
  EXPECT_EQ(MyVec3(4, 9, 1), z.MulComponents(z));
  EXPECT_EQ(MyVec3(2, 3, 1), z.Abs());
  EXPECT_EQ(MyVec3(2, 3, 1), z.Fabs());
  z.Set(-0.5, 1.5, 2.5);
  EXPECT_EQ(MyVec3(-1, 1, 2), z.Floor());
  EXPECT_EQ(MyVec3(0, 2, 3), z.Ceil());
  EXPECT_EQ(MyVec3(0, 2, 2), z.FRound());
  EXPECT_EQ(Vector3_i(0, 2, 2), z.IRound());
  y.Set(2.0, -2.0, 0.0);
  EXPECT_EQ(MyVec3(2.0, 1.5, 2.5), Max(z, y));
  EXPECT_EQ(MyVec3(-0.5, -2.0, 0.0), Min(z, y));
  y.Clear();
  EXPECT_EQ(MyVec3(0, 0, 0), y);
  ABSL_LOG(INFO) << "Test of Vector3 stream output: " << y;
}

TEST_F(ExactFloatTest, ImplicitConversions) {
  // Check that constants are implicitly converted to ExactFloats when they are
  // used in arithmetic expressions.

  // Check that ints/doubles are converted to ExactFloats in comparisons.
  EXPECT_EQ(1, ExactFloat(1.0));
  EXPECT_LE(1.5, ExactFloat(2.5));

  // Check that ints/doubles can be used as the first argument in arithmetic
  // operations.
  ExactFloat x = 11;
  EXPECT_EQ(24, 13 + x);
  EXPECT_EQ(10.5, 21.5 - x);
  EXPECT_EQ(-22, 2.0 * -x);

  // Intrinsic where one argument is an ExactFloat and the other is an int.
  EXPECT_EQ(3.0, fmax(ExactFloat(1), 3));
}

TEST_F(ExactFloatTest, ExactCalculations) {
  // Try a calculation that requires high precision.
  ExactFloat v1 = pow(2, 1000);
  ExactFloat v2 = v1 - pow(2, -1000);
  ExactFloat v3 = v1 + pow(2, -1000);
  EXPECT_EQ(v3 * v3 - v2 * v2, 4 * v1 * pow(2, -1000));
  EXPECT_EQ((v3 + v2) * (v3 - v2), 4 * v1 * pow(2, -1000));

  // Check that precision is increased as necessary to avoid rounding.
  typedef Vector3<ExactFloat> Vector3_xf;
  Vector3_xf x(1.1e300, 1.1e-300, 1.1);
  Vector3_xf y(1.1, 1.1e300, 1.1e-300);
  Vector3_xf z(1.1e-300, 1.1, 1.1e300);

  ExactFloat det = x.CrossProd(y).DotProd(z);
  EXPECT_GE(det.prec(), 6000);
  EXPECT_LT(det.prec(), 6300);
}

}  // namespace

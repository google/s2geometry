// Copyright 2003 Google Inc. All Rights Reserved.
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

// All Rights Reserved.
//
//
// A simple class to handle vectors in 2D
// The aim of this class is to be able to manipulate vectors in 2D
// as naturally as possible and make calculations readable.
// For that reason, the operators +, -, * are overloaded.
// (Reading a = a + b*2 - c is much easier to read than
// a = Sub(Add(a, Mul(b,2)),c)   )
// The code generated using this vector class is easily optimized by
// the compiler and does not generate overhead compared to manually
// writing the operations component by component
// (e.g a.x = b.x + c.x; a.y = b.y + c.y...)
//
// Operator overload is not usually allowed, but in this case an
// exemption has been granted by the C++ style committee.
//
// Please be careful about overflows when using those vectors with integer types
// The calculations are carried with the same type as the vector's components
// type. eg : if you are using uint8 as the base type, all values will be modulo
// 256.
// This feature is necessary to use the class in a more general framework with
// VType != plain old data type.

#ifndef UTIL_MATH_VECTOR2_H__
#define UTIL_MATH_VECTOR2_H__

#include <cstdlib>
#include <algorithm>
#include <cmath>
#include <iosfwd>
#include <iostream>  // NOLINT(readability/streams)
#include <limits>
#include <type_traits>

#include <glog/logging.h>

#include "s2/base/integral_types.h"
#include "s2/base/macros.h"
#include "s2/base/template_util.h"
#include "s2/util/math/vector.h"

// TODO(user): Look into creating conversion operators to remove the
// need to forward-declare.
template <typename VType> class Vector3;
template <typename VType> class Vector4;

// Template class for 2D vectors.
template <typename VType>
class Vector2 {
 private:
  VType c_[2];

 public:
  // FloatType is the type returned by Norm() and Angle().  These methods are
  // special because they return floating-point values even when VType is an
  // integer.
  typedef typename std::conditional<std::is_integral<VType>::value,
                             double, VType>::type FloatType;
  typedef Vector2<VType> Self;
  typedef VType BaseType;

  // The size, made accessible at compile time.
  enum { SIZE = 2 };

  // Create a new vector (0,0)
  Vector2();
  // Create a new vector (x,y)
  Vector2(const VType x, const VType y);
  // Keep only the two first coordinates of the vector vb
  explicit Vector2(const Vector3<VType> &vb);
  // Keep only the two first coordinates of the vector vb
  explicit Vector2(const Vector4<VType> &vb);
  // Convert from another vector type
  template <typename VType2>
  static Vector2 Cast(const Vector2<VType2> &vb);
  // Return the size of the vector
  static int Size() { return SIZE; }
  // Modify the coordinates of the current vector
  void Set(const VType x, const VType y);
  // Add two vectors, component by component
  Vector2& operator+=(const Vector2 &vb);
  // Subtract two vectors, component by component
  Vector2& operator-=(const Vector2 &vb);
  // Multiply a vector by a scalar
  Vector2& operator*=(const VType k);
  // Divide a vector by a scalar
  Vector2& operator/=(const VType k);
  // Multiply two vectors component by component
  Vector2 MulComponents(const Vector2 &vb) const;
  // Divide two vectors component by component
  Vector2 DivComponents(const Vector2 &vb) const;
  // Add two vectors, component by component
  Vector2 operator+(const Vector2 &vb) const;
  // Subtract two vectors, component by component
  Vector2 operator-(const Vector2 &vb) const;
  // Change the sign of the components of a vector
  Vector2 operator-() const;
  // Dot product.  Be aware that if VType is an integer type, the high bits of
  // the result are silently discarded.
  VType DotProd(const Vector2 &vb) const;
  // Multiplication by a scalar
  Vector2 operator*(const VType k) const;
  // Division by a scalar
  Vector2 operator/(const VType k) const;
  // Cross product.  Be aware that if VType is an integer type, the high bits
  // of the result are silently discarded.
  VType CrossProd(const Vector2 &vb) const;
  // Access component #b for read/write operations
  VType& operator[](const int b);
  // Access component #b for read only operations
  VType operator[](const int b) const;
  // Labeled Accessor methods.
  void x(const VType &v);
  VType x() const;
  void y(const VType &v);
  VType y() const;

  // return a pointer to the data array for interface with other libraries
  // like opencv
  VType* Data();
  const VType* Data() const;
  // Return the squared Euclidean norm of the vector.  Be aware that if VType
  // is an integer type, the high bits of the result are silently discarded.
  VType Norm2(void) const;
  // Return the Euclidean norm of the vector.  Note that if VType is an
  // integer type, the return value is correct only if the *squared* norm does
  // not overflow VType.
  FloatType Norm(void) const;
  // return the angle between "this" and v in radians
  FloatType Angle(const Vector2 &v) const;
  // Return a normalized version of the vector if the norm of the
  // vector is not 0.  Not to be used with integer types.
  Vector2 Normalize() const;
  // Compare two vectors, return true if all their components are equal
  // this operator is mostly useful for  integer types
  // for floating point types prefer "aequal"
  bool operator==(const Vector2 &vb) const;
  bool operator!=(const Vector2 &vb) const;
  // Compare two vectors, return true if all their components are within
  // a difference of margin.
  bool aequal(const Vector2 &vb, FloatType margin) const;

  // Compare two vectors, these comparisons are mostly for interaction
  // with STL.
  bool operator<(const Vector2 &vb) const;
  bool operator>(const Vector2 &vb) const;
  bool operator<=(const Vector2 &vb) const;
  bool operator>=(const Vector2 &vb) const;

  // return a vector orthogonal to the current one
  // with the same norm and counterclockwise to it
  Vector2 Ortho() const;
  // take the sqrt of each component and return a vector containing those values
  Vector2 Sqrt() const;
  // Take the fabs of each component and return a vector containing
  // those values.
  Vector2 Fabs() const;
  // Take the absolute value of each component and return a vector containing
  // those values.  This method should only be used when VType is a signed
  // integer type that is not wider than "int".
  Vector2 Abs() const;
  // take the floor of each component and return a vector containing
  // those values
  Vector2 Floor() const;
  // Take the ceil of each component and return a vector containing
  // those values.
  Vector2 Ceil() const;
  // take the round of each component and return a vector containing those
  // values
  Vector2 FRound() const;
  // take the round of each component and return an integer vector containing
  // those values
  Vector2<int> IRound() const;
  // Reset all the coordinates of the vector to 0
  void Clear();

  // return true if one of the components is not a number
  bool IsNaN() const;

  // return an invalid floating point vector
  static Vector2 NaN();
};

//
// Inline definitions of member functions
//

template <typename VType>
Vector2<VType>::Vector2() {
  Clear();
}
template <typename VType>
Vector2<VType>::Vector2(const VType x, const VType y) {
  c_[0] = x;
  c_[1] = y;
}
template <typename VType>
Vector2<VType>::Vector2(const Vector3<VType> &vb) {
  c_[0] = vb.x();
  c_[1] = vb.y();
}
template <typename VType>
Vector2<VType>::Vector2(const Vector4<VType> &vb) {
  c_[0] = vb.x();
  c_[1] = vb.y();
}

template <typename VType> template <typename VType2>
Vector2<VType> Vector2<VType>::Cast(const Vector2<VType2> &vb) {
  return Vector2(static_cast<VType>(vb[0]),
              static_cast<VType>(vb[1]));
}

template <typename VType>
void Vector2<VType>::Set(const VType x, const VType y) {
  c_[0] = x;
  c_[1] = y;
}

template <typename VType>
Vector2<VType>& Vector2<VType>::operator+=(const Vector2 &vb) {
  c_[0] += vb.c_[0];
  c_[1] += vb.c_[1];
  return (*this);
}

template <typename VType>
Vector2<VType>& Vector2<VType>::operator-=(const Vector2 &vb) {
  c_[0] -= vb.c_[0];
  c_[1] -= vb.c_[1];
  return (*this);
}

template <typename VType>
Vector2<VType>& Vector2<VType>::operator*=(const VType k) {
  c_[0] *= k;
  c_[1] *= k;
  return (*this);
}

template <typename VType>
Vector2<VType>& Vector2<VType>::operator/=(const VType k) {
  c_[0] /= k;
  c_[1] /= k;
  return (*this);
}

template <typename VType>
Vector2<VType> Vector2<VType>::MulComponents(const Vector2 &vb) const {
  return Vector2(c_[0] * vb.c_[0], c_[1] * vb.c_[1]);
}

template <typename VType>
Vector2<VType> Vector2<VType>::DivComponents(const Vector2 &vb) const {
  return Vector2(c_[0] / vb.c_[0], c_[1] / vb.c_[1]);
}

template <typename VType>
Vector2<VType> Vector2<VType>::operator+(const Vector2 &vb) const {
  return Vector2(*this) += vb;
}

template <typename VType>
Vector2<VType> Vector2<VType>::operator-(const Vector2 &vb) const {
  return Vector2(*this) -= vb;
}

template <typename VType>
Vector2<VType> Vector2<VType>::operator-() const {
  return Vector2(-c_[0], -c_[1]);
}

template <typename VType>
VType Vector2<VType>::DotProd(const Vector2 &vb) const {
  return c_[0] * vb.c_[0] + c_[1] * vb.c_[1];
}

template <typename VType>
Vector2<VType> Vector2<VType>::operator*(const VType k) const {
  return Vector2(*this) *= k;
}

template <typename VType>
Vector2<VType> Vector2<VType>::operator/(const VType k) const {
  return Vector2(*this) /= k;
}

template <typename VType>
VType Vector2<VType>::CrossProd(const Vector2 &vb) const {
  return c_[0] * vb.c_[1] - c_[1] * vb.c_[0];
}

template <typename VType>
VType& Vector2<VType>::operator[](const int b) {
  DCHECK_GE(b, 0);
  DCHECK_LE(b, 1);
  return c_[b];
}

template <typename VType>
VType Vector2<VType>::operator[](const int b) const {
  DCHECK_GE(b, 0);
  DCHECK_LE(b, 1);
  return c_[b];
}

template <typename VType>
void Vector2<VType>::x(const VType &v) {
  c_[0] = v;
}

template <typename VType>
VType Vector2<VType>::x() const {
  return c_[0];
}

template <typename VType>
void Vector2<VType>::y(const VType &v) {
  c_[1] = v;
}

template <typename VType>
VType Vector2<VType>::y() const {
  return c_[1];
}

template <typename VType>
VType* Vector2<VType>::Data() {
  return reinterpret_cast<VType*>(c_);
}

template <typename VType>
const VType* Vector2<VType>::Data() const {
  return reinterpret_cast<const VType*>(c_);
}

template <typename VType>
VType Vector2<VType>::Norm2(void) const {
  return c_[0]*c_[0] + c_[1]*c_[1];
}

template <typename VType>
typename Vector2<VType>::FloatType Vector2<VType>::Norm(void) const {
  return sqrt(Norm2());
}

template <typename VType>
typename Vector2<VType>::FloatType
Vector2<VType>::Angle(const Vector2 &v) const {
  return atan2(this->CrossProd(v), this->DotProd(v));
}

template <typename VType>
Vector2<VType> Vector2<VType>::Normalize() const {
  static_assert(!std::is_integral<VType>::value, "must be floating point");
  VType n = Norm();
  if (n != VType(0)) {
    n = VType(1.0) / n;
  }
  return Vector2(*this) *= n;
}

template <typename VType>
bool Vector2<VType>::operator==(const Vector2 &vb) const {
  return  (c_[0] == vb.c_[0]) && (c_[1] == vb.c_[1]);
}

template <typename VType>
bool Vector2<VType>::operator!=(const Vector2 &vb) const {
  return  (c_[0] != vb.c_[0]) || (c_[1] != vb.c_[1]);
}

template <typename VType>
bool Vector2<VType>::aequal(const Vector2 &vb, FloatType margin) const {
  return (fabs(c_[0]-vb.c_[0]) < margin) && (fabs(c_[1]-vb.c_[1]) < margin);
}

template <typename VType>
bool Vector2<VType>::operator<(const Vector2 &vb) const {
  if ( c_[0] < vb.c_[0] ) return true;
  if ( vb.c_[0] < c_[0] ) return false;
  if ( c_[1] < vb.c_[1] ) return true;
  return false;
}

template <typename VType>
bool Vector2<VType>::operator>(const Vector2 &vb) const {
  return vb.operator<(*this);
}

template <typename VType>
bool Vector2<VType>::operator<=(const Vector2 &vb) const {
  return !operator>(vb);
}

template <typename VType>
bool Vector2<VType>::operator>=(const Vector2 &vb) const {
  return !operator<(vb);
}

template <typename VType>
Vector2<VType> Vector2<VType>::Ortho() const {
  return Vector2(-c_[1], c_[0]);
}

template <typename VType>
Vector2<VType> Vector2<VType>::Sqrt() const {
  return Vector2(sqrt(c_[0]), sqrt(c_[1]));
}

template <typename VType>
Vector2<VType> Vector2<VType>::Fabs() const {
  return Vector2(fabs(c_[0]), fabs(c_[1]));
}

template <typename VType>
Vector2<VType> Vector2<VType>::Abs() const {
  static_assert(std::is_integral<VType>::value, "use Fabs for float_types");
  static_assert(static_cast<VType>(-1) == -1, "type must be signed");
  static_assert(sizeof(c_[0]) <= sizeof(int), "Abs truncates to int");
  return Vector2(abs(c_[0]), abs(c_[1]));
}

template <typename VType>
Vector2<VType> Vector2<VType>::Floor() const {
  return Vector2(floor(c_[0]), floor(c_[1]));
}

template <typename VType>
Vector2<VType> Vector2<VType>::Ceil() const {
  return Vector2(ceil(c_[0]), ceil(c_[1]));
}

template <typename VType>
Vector2<VType> Vector2<VType>::FRound() const {
  return Vector2(rint(c_[0]), rint(c_[1]));
}

template <typename VType>
Vector2<int> Vector2<VType>::IRound() const {
  return Vector2<int>(lrint(c_[0]), lrint(c_[1]));
}

template <typename VType>
void Vector2<VType>::Clear() {
  c_[1] = c_[0] = VType();
}

template <typename VType>
bool Vector2<VType>::IsNaN() const {
  return isnan(c_[0]) || isnan(c_[1]);
}

template <typename VType>
Vector2<VType> Vector2<VType>::NaN() {
  using Lim = std::numeric_limits<VType>;
  return Vector2(Lim::quiet_NaN(), Lim::quiet_NaN());
}

template <typename ScalarType, typename VType2>
Vector2<VType2> operator*(const ScalarType k, const Vector2<VType2> v) {
  return Vector2<VType2>( k * v[0], k * v[1]);
}

template <typename ScalarType, typename VType2>
Vector2<VType2> operator/(const ScalarType k, const Vector2<VType2> v) {
  return Vector2<VType2>(k / v[0], k / v[1]);
}

template <typename VType>
Vector2<VType> Max(const Vector2<VType> &v1, const Vector2<VType> &v2) {
  return Vector2<VType>(max(v1[0], v2[0]), max(v1[1], v2[1]));
}

template <typename VType>
Vector2<VType> Min(const Vector2<VType> &v1, const Vector2<VType> &v2) {
  return Vector2<VType>(min(v1[0], v2[0]), min(v1[1], v2[1]));
}

template <typename VType>
std::ostream &operator <<(std::ostream &out, const Vector2<VType> &va) {
  out << "["
      << va[0] << ", "
      << va[1] << "]";
  return out;
}

template <>
inline std::ostream &operator <<(std::ostream &out, const Vector2<uint8> &va) {
  // ostream << uint8 prints the ASCII character, which is not useful.
  // Cast to int so that numbers will be printed instead.
  out << Vector2<int>::Cast(va);
  return out;
}

// TODO(user): Declare extern templates for these types.
typedef Vector2<uint8>  Vector2_b;
typedef Vector2<int>    Vector2_i;
typedef Vector2<float>  Vector2_f;
typedef Vector2<double> Vector2_d;

// TODO(user): Vector2<T> does not actually satisfy the definition of a POD
// type even when T is a POD. Pretending that Vector2<T> is a POD probably
// won't cause any immediate problems, but eventually this should be fixed.


#endif  // UTIL_MATH_VECTOR2_H__

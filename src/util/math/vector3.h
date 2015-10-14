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
// A simple class to handle vectors in 3D
// The aim of this class is to be able to manipulate vectors in 3D
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

#ifndef UTIL_MATH_VECTOR3_H__
#define UTIL_MATH_VECTOR3_H__

#include <stdlib.h>
#include <algorithm>
#include <cmath>
#include <iosfwd>
#include <iostream>  // NOLINT(readability/streams)

#include <glog/logging.h>

#include "base/integral_types.h"
#include "base/macros.h"
#include "base/template_util.h"
#include "base/type_traits.h"
#include "util/math/mathutil.h"
#include "util/math/vector2.h"
#include "util/math/vector4.h"

// TODO(user): Look into creating conversion operators to remove the
// need to forward-declare.
template <typename VType> class Vector2;
template <typename VType> class Vector4;

// Template class for 3D vectors.
template <typename VType>
class Vector3 {
 private:
  VType c_[3];

 public:
  // FloatType is the type returned by Norm() and Angle().  These methods are
  // special because they return floating-point values even when VType is an
  // integer.
  typedef typename std::conditional<std::is_integral<VType>::value,
                             double, VType>::type FloatType;
  typedef Vector3<VType> Self;
  typedef VType BaseType;

  // The size, made accessible at compile time.
  enum { SIZE = 3 };

  // Create a new vector (0,0,0)
  Vector3();
  // Create a new vector (x,y,z)
  Vector3(const VType x, const VType y, const VType z);
  // Create a new 3D vector using the two first coordinates of a 2D vectors
  // and an additional z argument.
  Vector3(const Vector2<VType> &vb, VType z);
  // Keep only the three first coordinates of the 4D vector vb
  explicit Vector3(const Vector4<VType> &vb);
  // Convert from another vector type
  template <typename VType2>
  static Self Cast(const Vector3<VType2> &vb);
  // Compare two vectors, return true if all their components are equal
  bool operator==(const Self& vb) const;
  bool operator!=(const Self& vb) const;
  // Compare two vectors, return true if all their components are within
  // a difference of margin.
  bool aequal(const Self &vb, FloatType margin) const;
  // Compare two vectors, these comparisons are mostly for interaction
  // with STL.
  bool operator<(const Self &vb) const;
  bool operator>(const Self &vb) const;
  bool operator<=(const Self &vb) const;
  bool operator>=(const Self &vb) const;

  // Return the size of the vector
  static int Size() { return SIZE; }
  // Modify the coordinates of the current vector
  void Set(const VType x, const VType y, const VType z);
  // Add two vectors, component by component
  Self& operator+=(const Self &vb);
  // Subtract two vectors, component by component
  Self& operator-=(const Self &vb);
  // Multiply a vector by a scalar
  Self& operator*=(const VType k);
  // Divide a vector by a scalar
  Self& operator/=(const VType k);
  // Multiply two vectors component by component
  Self MulComponents(const Self &vb) const;
  // Divide two vectors component by component
  Self DivComponents(const Self &vb) const;
  // Add two vectors, component by component
  Self operator+(const Self &vb) const;
  // Subtract two vectors, component by component
  Self operator-(const Self &vb) const;
  // Dot product.  Be aware that if VType is an integer type, the high bits of
  // the result are silently discarded.
  VType DotProd(const Self &vb) const;
  // Multiplication by a scalar
  Self operator*(const VType k) const;
  // Divide by a scalar
  Self operator/(const VType k) const;
  // Cross product.  Be aware that if VType is an integer type, the high bits
  // of the result are silently discarded.
  Self CrossProd(const Self& vb) const;
  // Access component #b for read/write operations
  VType& operator[](const int b);
  // Access component #b for read only operations
  VType operator[](const int b) const;
  // Labeled Accessor methods.
  void x(const VType &v);
  VType x() const;
  void y(const VType &v);
  VType y() const;
  void z(const VType &v);
  VType z() const;
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
  // Return a normalized version of the vector if the norm of the
  // vector is not 0.  Not to be used with integer types.
  Self Normalize() const;
  // return a vector orthogonal to this one
  Self Ortho() const;
  // return the index of the largest component (fabs)
  int LargestAbsComponent() const;
  // return the index of the smallest, median ,largest component of the vector
  Vector3<int> ComponentOrder() const;
  // return the angle between two vectors in radians
  FloatType Angle(const Self &va) const;
  // take the sqrt of each component and return a vector containing those values
  Self Sqrt() const;
  // take the fabs of each component and return a vector containing those values
  Self Fabs() const;
  // Take the absolute value of each component and return a vector containing
  // those values.
  Self Abs() const;
  // take the floor of each component and return a vector containing
  // those values
  Self Floor() const;
  // take the ceil of each component and return a vector containing those values
  Self Ceil() const;
  // take the round of each component and return a vector containing those
  // values
  Self FRound() const;
  // take the round of each component and return an integer vector containing
  // those values
  Vector3<int> IRound() const;
  // Reset all the coordinates of the vector to 0
  void Clear();

  // return true if one of the components is not a number
  bool IsNaN() const;

  // return an invalid floating point vector
  static Self NaN();
};

//
// Inline definitions of member functions.
//

template <typename VType>
Vector3<VType>::Vector3() {
  Clear();
}

template <typename VType>
Vector3<VType>::Vector3(const VType x, const VType y, const VType z) {
  c_[0] = x;
  c_[1] = y;
  c_[2] = z;
}

template <typename VType>
Vector3<VType>::Vector3(const Vector2<VType> &vb, VType z) {
  c_[0] = vb.x();
  c_[1] = vb.y();
  c_[2] = z;
}

template <typename VType>
Vector3<VType>::Vector3(const Vector4<VType> &vb) {
  c_[0] = vb.x();
  c_[1] = vb.y();
  c_[2] = vb.z();
}

template <typename VType> template <typename VType2>
Vector3<VType> Vector3<VType>::Cast(const Vector3<VType2> &vb) {
  return Self(VType(vb[0]),
              VType(vb[1]),
              VType(vb[2]));
}

template <typename VType>
bool Vector3<VType>::operator==(const Self& vb) const {
  return  (c_[0] == vb.c_[0]) && (c_[1] == vb.c_[1]) && (c_[2] == vb.c_[2]);
}

template <typename VType>
bool Vector3<VType>::operator!=(const Self& vb) const {
  return  (c_[0] != vb.c_[0]) || (c_[1] != vb.c_[1]) || (c_[2] != vb.c_[2]);
}

template <typename VType>
bool Vector3<VType>::aequal(const Self &vb, FloatType margin) const {
  using std::abs;
  return (abs(c_[0] - vb.c_[0]) < margin)
      && (abs(c_[1] - vb.c_[1]) < margin)
      && (abs(c_[2] - vb.c_[2]) < margin);
}

template <typename VType>
bool Vector3<VType>::operator<(const Self &vb) const {
  if ( c_[0] < vb.c_[0] ) return true;
  if ( vb.c_[0] < c_[0] ) return false;
  if ( c_[1] < vb.c_[1] ) return true;
  if ( vb.c_[1] < c_[1] ) return false;
  if ( c_[2] < vb.c_[2] ) return true;
  return false;
}

template <typename VType>
bool Vector3<VType>::operator>(const Self &vb) const {
  return vb.operator<(*this);
}

template <typename VType>
bool Vector3<VType>::operator<=(const Self &vb) const {
  return !operator>(vb);
}

template <typename VType>
bool Vector3<VType>::operator>=(const Self &vb) const {
  return !operator<(vb);
}

template <typename VType>
void Vector3<VType>::Set(const VType x, const VType y, const VType z) {
  c_[0] = x;
  c_[1] = y;
  c_[2] = z;
}

template <typename VType>
Vector3<VType>& Vector3<VType>::operator+=(const Self &vb) {
  c_[0] += vb.c_[0];
  c_[1] += vb.c_[1];
  c_[2] += vb.c_[2];
  return (*this);
}

template <typename VType>
Vector3<VType>& Vector3<VType>::operator-=(const Self &vb) {
  c_[0] -= vb.c_[0];
  c_[1] -= vb.c_[1];
  c_[2] -= vb.c_[2];
  return (*this);
}

template <typename VType>
Vector3<VType>& Vector3<VType>::operator*=(const VType k) {
  c_[0] *= k;
  c_[1] *= k;
  c_[2] *= k;
  return (*this);
}

template <typename VType>
Vector3<VType>& Vector3<VType>::operator/=(const VType k) {
  c_[0] /= k;
  c_[1] /= k;
  c_[2] /= k;
  return (*this);
}

template <typename VType>
Vector3<VType> Vector3<VType>::MulComponents(const Self &vb) const {
  return Self(c_[0] * vb.c_[0], c_[1] * vb.c_[1], c_[2] * vb.c_[2]);
}

template <typename VType>
Vector3<VType> Vector3<VType>::DivComponents(const Self &vb) const {
  return Self(c_[0] / vb.c_[0], c_[1] / vb.c_[1], c_[2] / vb.c_[2]);
}

template <typename VType>
Vector3<VType> Vector3<VType>::operator+(const Self &vb) const {
  return Self(*this) += vb;
}

template <typename VType>
Vector3<VType> Vector3<VType>::operator-(const Self &vb) const {
  return Self(*this) -= vb;
}

template <typename VType>
VType Vector3<VType>::DotProd(const Self &vb) const {
  return c_[0]*vb.c_[0] + c_[1]*vb.c_[1] + c_[2]*vb.c_[2];
}

template <typename VType>
Vector3<VType> Vector3<VType>::operator*(const VType k) const {
  return Self(*this) *= k;
}

template <typename VType>
Vector3<VType> Vector3<VType>::operator/(const VType k) const {
  return Self(*this) /= k;
}

template <typename VType>
Vector3<VType> Vector3<VType>::CrossProd(const Self& vb) const {
  return Self( c_[1] * vb.c_[2] - c_[2] * vb.c_[1],
               c_[2] * vb.c_[0] - c_[0] * vb.c_[2],
               c_[0] * vb.c_[1] - c_[1] * vb.c_[0]);
}

template <typename VType>
VType& Vector3<VType>::operator[](const int b) {
  DCHECK_GE(b, 0);
  DCHECK_LE(b, 2);
  return c_[b];
}

template <typename VType>
VType Vector3<VType>::operator[](const int b) const {
  DCHECK_GE(b, 0);
  DCHECK_LE(b, 2);
  return c_[b];
}

template <typename VType>
void Vector3<VType>::x(const VType &v) {
  c_[0] = v;
}

template <typename VType>
VType Vector3<VType>::x() const {
  return c_[0];
}

template <typename VType>
void Vector3<VType>::y(const VType &v) {
  c_[1] = v;
}

template <typename VType>
VType Vector3<VType>::y() const {
  return c_[1];
}

template <typename VType>
void Vector3<VType>::z(const VType &v) {
  c_[2] = v;
}

template <typename VType>
VType Vector3<VType>::z() const {
  return c_[2];
}

template <typename VType>
VType* Vector3<VType>::Data() {
  return reinterpret_cast<VType*>(c_);
}

template <typename VType>
const VType* Vector3<VType>::Data() const {
  return reinterpret_cast<const VType*>(c_);
}

template <typename VType>
VType Vector3<VType>::Norm2(void) const {
  return c_[0]*c_[0] + c_[1]*c_[1] + c_[2]*c_[2];
}

template <typename VType>
typename Vector3<VType>::FloatType Vector3<VType>::Norm(void) const {
  return sqrt(Norm2());
}

template <typename VType>
Vector3<VType> Vector3<VType>::Normalize() const {
  COMPILE_ASSERT(!std::is_integral<VType>::value, must_be_floating_point);
  VType n = Norm();
  if (n != VType(0.0)) {
    n = VType(1.0) / n;
  }
  return Self(*this) *= n;
}

template <typename VType>
Vector3<VType> Vector3<VType>::Ortho() const {
  int k = LargestAbsComponent() - 1;
  if (k < 0) k = 2;
  Self temp;
  temp[k] = VType(1);
  return (this->CrossProd(temp)).Normalize();
}

template <typename VType>
int Vector3<VType>::LargestAbsComponent() const {
  Self temp = Abs();
  if (temp[0] > temp[1]) {
    if (temp[0] > temp[2]) {
      return 0;
    } else {
      return 2;
    }
  } else {
    if (temp[1] > temp[2]) {
      return 1;
    } else {
      return 2;
    }
  }
}

template <typename VType>
Vector3<int> Vector3<VType>::ComponentOrder() const {
  using std::swap;
  Vector3<int> temp(0, 1, 2);
  if (c_[temp[0]] > c_[temp[1]]) swap(temp[0], temp[1]);
  if (c_[temp[1]] > c_[temp[2]]) swap(temp[1], temp[2]);
  if (c_[temp[0]] > c_[temp[1]]) swap(temp[0], temp[1]);
  return temp;
}

template <typename VType>
typename Vector3<VType>::FloatType Vector3<VType>::Angle(const Self &va) const {
  return atan2(this->CrossProd(va).Norm(), this->DotProd(va));
}

template <typename VType>
Vector3<VType> Vector3<VType>::Sqrt() const {
  return Self(sqrt(c_[0]), sqrt(c_[1]), sqrt(c_[2]));
}

template <typename VType>
Vector3<VType> Vector3<VType>::Fabs() const {
  return Abs();
}

template <typename VType>
Vector3<VType> Vector3<VType>::Abs() const {
  COMPILE_ASSERT(
      !std::is_integral<VType>::value || static_cast<VType>(-1) == -1,
      type_must_be_signed);
  using std::abs;
  return Self(abs(c_[0]), abs(c_[1]), abs(c_[2]));
}

template <typename VType>
Vector3<VType> Vector3<VType>::Floor() const {
  return Self(floor(c_[0]), floor(c_[1]), floor(c_[2]));
}

template <typename VType>
Vector3<VType> Vector3<VType>::Ceil() const {
  return Self(ceil(c_[0]), ceil(c_[1]), ceil(c_[2]));
}

template <typename VType>
Vector3<VType> Vector3<VType>::FRound() const {
  return Self(rint(c_[0]), rint(c_[1]), rint(c_[2]));
}

template <typename VType>
Vector3<int> Vector3<VType>::IRound() const {
  return Vector3<int>(lrint(c_[0]), lrint(c_[1]), lrint(c_[2]));
}

template <typename VType>
void Vector3<VType>::Clear() {
  c_[2] = c_[1] = c_[0] = VType();
}

template <typename VType>
bool Vector3<VType>::IsNaN() const {
  return isnan(c_[0]) || isnan(c_[1]) || isnan(c_[2]);
}

template <typename VType>
Vector3<VType> Vector3<VType>::NaN() {
  return Self(MathUtil::NaN(), MathUtil::NaN(), MathUtil::NaN());
}

template <typename VType>
Vector3<VType> operator-(const Vector3<VType> &vb) {
  return Vector3<VType>(-vb[0], -vb[1], -vb[2]);
}

template <typename ScalarType, typename VType>
Vector3<VType> operator*(const ScalarType k, const Vector3<VType> &v) {
  return Vector3<VType>(k*v[0], k*v[1], k*v[2]);
}

template <typename ScalarType, typename VType>
Vector3<VType> operator/(const ScalarType k, const Vector3<VType> &v) {
  return Vector3<VType>(k/v[0], k/v[1], k/v[2]);
}

template <typename VType>
Vector3<VType> Max(const Vector3<VType> &v1, const Vector3<VType> &v2) {
  return Vector3<VType>(std::max(v1[0], v2[0]), std::max(v1[1], v2[1]),
                        std::max(v1[2], v2[2]));
}

template <typename VType>
Vector3<VType> Min(const Vector3<VType> &v1, const Vector3<VType> &v2) {
  return Vector3<VType>(std::min(v1[0], v2[0]), std::min(v1[1], v2[1]),
                        std::min(v1[2], v2[2]));
}

template <typename VType>
std::ostream &operator <<(std::ostream &out, const Vector3<VType> &va) {
  out << "["
      << va[0] << ", "
      << va[1] << ", "
      << va[2] << "]";
  return out;
}

template <>
inline std::ostream &operator <<(std::ostream &out, const Vector3<uint8> &va) {
  // ostream << uint8 prints the ASCII character, which is not useful.
  // Cast to int so that numbers will be printed instead.
  out << Vector3<int>::Cast(va);
  return out;
}

// TODO(user): Declare extern templates for these types.
typedef Vector3<uint8>  Vector3_b;
typedef Vector3<int>    Vector3_i;
typedef Vector3<float>  Vector3_f;
typedef Vector3<double> Vector3_d;

// TODO(user): Vector3<T> does not actually satisfy the definition of a POD
// type even when T is a POD. Pretending that Vector3<T> is a POD probably
// won't cause any immediate problems, but eventually this should be fixed.


#endif  // UTIL_MATH_VECTOR3_H__

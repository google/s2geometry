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
// A simple class to handle vectors in 4D
// The aim of this class is to be able to manipulate vectors in 4D
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

#ifndef UTIL_MATH_VECTOR4_H__
#define UTIL_MATH_VECTOR4_H__

#include <cstdlib>
#include <algorithm>
#include <cmath>
#include <iosfwd>
#include <iostream>  // NOLINT(readability/streams)

#include <glog/logging.h>

#include "base/integral_types.h"
#include "base/macros.h"
#include "base/template_util.h"
#include <type_traits>
#include "util/math/mathutil.h"
#include "util/math/vector2.h"
#include "util/math/vector3.h"

// TODO(user): Look into creating conversion operators to remove the
// need to forward-declare.
template <typename VType> class Vector2;
template <typename VType> class Vector3;

// Template class for 4D vectors
template <typename VType>
class Vector4 {
 private:
  VType c_[4];

 public:
  // FloatType is the type returned by Norm().  This method is special because
  // it returns floating-point values even when VType is an integer.
  typedef typename std::conditional<std::is_integral<VType>::value,
                             double, VType>::type FloatType;
  typedef Vector4<VType> Self;
  typedef VType BaseType;

  // The size, made accessible at compile time.
  enum { SIZE = 4 };

  // Create a new vector (0,0)
  Vector4();
  // Create a new vector (x,y,z,w)
  Vector4(const VType x, const VType y, const VType z, const VType w);
  // Create a new 4D vector from 2D vector and two scalars
  // (vb.x,vb.y,z,w)
  Vector4(const Vector2<VType> &vb, const VType z, const VType w);
  // Create a 4D vector from two 2D vectors (vb1.x,vb1.y,vb2.x,vb2.y)
  Vector4(const Vector2<VType> &vb1, const Vector2<VType> &vb2);
  // Create a new 4D vector from 3D vector and scalar
  // (vb.x,vb.y,vb.z,w)
  Vector4(const Vector3<VType> &vb, const VType w);
  // Convert from another vector type
  template <typename VType2>
  static Vector4 Cast(const Vector4<VType2> &vb);
  // Compare two vectors, return true if all their components are equal
  bool operator==(const Vector4& vb) const;
  bool operator!=(const Vector4& vb) const;
  // Compare two vectors, return true if all their components are within
  // a difference of margin.
  bool aequal(const Vector4 &vb, FloatType margin) const;
  // Compare two vectors, these comparisons are mostly for interaction
  // with STL.
  bool operator<(const Vector4 &vb) const;
  bool operator>(const Vector4 &vb) const;
  bool operator<=(const Vector4 &vb) const;
  bool operator>=(const Vector4 &vb) const;

  // Return the size of the vector
  static int Size() { return SIZE; }
  // Modify the coordinates of the current vector
  void Set(const VType x, const VType y, const VType z, const VType w);
  // add two vectors, component by component
  Vector4& operator+=(const Vector4& vb);
  // subtract two vectors, component by component
  Vector4& operator-=(const Vector4& vb);
  // multiply a vector by a scalar
  Vector4& operator*=(const VType k);
  // divide a vector by a scalar : implemented that way for integer vectors
  Vector4& operator/=(const VType k);
  // multiply two vectors component by component
  Vector4 MulComponents(const Vector4 &vb) const;
  // divide two vectors component by component
  Vector4 DivComponents(const Vector4 &vb) const;
  // add two vectors, component by component
  Vector4 operator+(const Vector4 &vb) const;
  // subtract two vectors, component by component
  Vector4 operator-(const Vector4 &vb) const;
  // Dot product.  Be aware that if VType is an integer type, the high bits of
  // the result are silently discarded.
  VType DotProd(const Vector4 &vb) const;
  // Multiplication by a scalar
  Vector4 operator*(const VType k) const;
  // Division by a scalar
  Vector4 operator/(const VType k) const;
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
  void w(const VType &v);
  VType w() const;
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
  Vector4 Normalize() const;
  // take the sqrt of each component and return a vector containing those values
  Vector4 Sqrt() const;
  // take the fabs of each component and return a vector containing those values
  Vector4 Fabs() const;
  // Take the absolute value of each component and return a vector containing
  // those values.  This method should only be used when VType is a signed
  // integer type that is not wider than "int".
  Vector4 Abs() const;
  // take the floor of each component and return a vector containing
  // those values
  Vector4 Floor() const;
  // take the ceil of each component and return a vector containing those values
  Vector4 Ceil() const;
  // take the round of each component and return a vector containing those
  // values
  Vector4 FRound() const;
  // take the round of each component and return an integer vector containing
  // those values
  Vector4<int> IRound() const;
  // Reset all the coordinates of the vector to 0
  void Clear();

  // return true if one of the components is not a number
  bool IsNaN() const;

  // return an invalid floating point vector
  static Vector4 NaN();
};

//
// Inline definitions of member functions.
//

template <typename VType>
Vector4<VType>::Vector4() {
  Clear();
}

template <typename VType>
Vector4<VType>::Vector4(const VType x, const VType y, const VType z,
                        const VType w) {
  c_[0] = x;
  c_[1] = y;
  c_[2] = z;
  c_[3] = w;
}

template <typename VType>
Vector4<VType>::Vector4(const Vector2<VType> &vb, const VType z,
                        const VType w) {
  c_[0] = vb.x();
  c_[1] = vb.y();
  c_[2] = z;
  c_[3] = w;
}

template <typename VType>
Vector4<VType>::Vector4(const Vector2<VType> &vb1, const Vector2<VType> &vb2) {
  c_[0] = vb1.x();
  c_[1] = vb1.y();
  c_[2] = vb2.x();
  c_[3] = vb2.y();
}

template <typename VType>
Vector4<VType>::Vector4(const Vector3<VType> &vb, const VType w) {
  c_[0] = vb.x();
  c_[1] = vb.y();
  c_[2] = vb.z();
  c_[3] = w;
}

template <typename VType> template <typename VType2>
Vector4<VType> Vector4<VType>::Cast(const Vector4<VType2> &vb) {
  return Vector4(static_cast<VType>(vb[0]),
                 static_cast<VType>(vb[1]),
                 static_cast<VType>(vb[2]),
                 static_cast<VType>(vb[3]));
}

template <typename VType>
bool Vector4<VType>::operator==(const Vector4& vb) const {
  return  (c_[0] == vb.c_[0])
      && (c_[1] == vb.c_[1])
      && (c_[2] == vb.c_[2])
      && (c_[3] == vb.c_[3]);
}

template <typename VType>
bool Vector4<VType>::operator!=(const Vector4& vb) const {
  return  (c_[0] != vb.c_[0])
      || (c_[1] != vb.c_[1])
      || (c_[2] != vb.c_[2])
      || (c_[3] != vb.c_[3]);
}

template <typename VType>
bool Vector4<VType>::aequal(const Vector4 &vb, FloatType margin) const {
  return (fabs(c_[0] - vb.c_[0]) < margin)
      && (fabs(c_[1] - vb.c_[1]) < margin)
      && (fabs(c_[2] - vb.c_[2]) < margin)
      && (fabs(c_[3] - vb.c_[3]) < margin);
}

template <typename VType>
bool Vector4<VType>::operator<(const Vector4 &vb) const {
  if ( c_[0] < vb.c_[0] ) return true;
  if ( vb.c_[0] < c_[0] ) return false;
  if ( c_[1] < vb.c_[1] ) return true;
  if ( vb.c_[1] < c_[1] ) return false;
  if ( c_[2] < vb.c_[2] ) return true;
  if ( vb.c_[2] < c_[2] ) return false;
  if ( c_[3] < vb.c_[3] ) return true;
  return false;
}

template <typename VType>
bool Vector4<VType>::operator>(const Vector4 &vb) const {
  return vb.operator<(*this);
}

template <typename VType>
bool Vector4<VType>::operator<=(const Vector4 &vb) const {
  return !operator>(vb);
}

template <typename VType>
bool Vector4<VType>::operator>=(const Vector4 &vb) const {
  return !operator<(vb);
}

template <typename VType>
void Vector4<VType>::Set(const VType x, const VType y, const VType z,
                         const VType w) {
  c_[0] = x;
  c_[1] = y;
  c_[2] = z;
  c_[3] = w;
}

template <typename VType>
Vector4<VType>& Vector4<VType>::operator+=(const Vector4& vb) {
  c_[0] += vb.c_[0];
  c_[1] += vb.c_[1];
  c_[2] += vb.c_[2];
  c_[3] += vb.c_[3];
  return (*this);
}

template <typename VType>
Vector4<VType>& Vector4<VType>::operator-=(const Vector4& vb) {
  c_[0] -= vb.c_[0];
  c_[1] -= vb.c_[1];
  c_[2] -= vb.c_[2];
  c_[3] -= vb.c_[3];
  return (*this);
}

template <typename VType>
Vector4<VType>& Vector4<VType>::operator*=(const VType k) {
  c_[0] *= k;
  c_[1] *= k;
  c_[2] *= k;
  c_[3] *= k;
  return (*this);
}

template <typename VType>
Vector4<VType>& Vector4<VType>::operator/=(const VType k) {
  c_[0] /= k;
  c_[1] /= k;
  c_[2] /= k;
  c_[3] /= k;
  return (*this);
}

template <typename VType>
Vector4<VType> Vector4<VType>::MulComponents(const Vector4 &vb) const {
  return Vector4(c_[0] * vb.c_[0], c_[1] * vb.c_[1],
                 c_[2] * vb.c_[2], c_[3] * vb.c_[3]);
}

template <typename VType>
Vector4<VType> Vector4<VType>::DivComponents(const Vector4 &vb) const {
  return Vector4(c_[0] / vb.c_[0], c_[1] / vb.c_[1],
                 c_[2] / vb.c_[2], c_[3] / vb.c_[3]);
}

template <typename VType>
Vector4<VType> Vector4<VType>::operator+(const Vector4 &vb) const {
  return Vector4(*this) += vb;
}

template <typename VType>
Vector4<VType> Vector4<VType>::operator-(const Vector4 &vb) const {
  return Vector4(*this) -= vb;
}

template <typename VType>
VType Vector4<VType>::DotProd(const Vector4 &vb) const {
  return c_[0]*vb.c_[0] + c_[1]*vb.c_[1] + c_[2]*vb.c_[2] + c_[3]*vb.c_[3];
}

template <typename VType>
Vector4<VType> Vector4<VType>::operator*(const VType k) const {
  return Vector4(*this) *= k;
}

template <typename VType>
Vector4<VType> Vector4<VType>::operator/(const VType k) const {
  return Vector4(*this) /= k;
}

template <typename VType>
VType& Vector4<VType>::operator[](const int b) {
  DCHECK_GE(b, 0);
  DCHECK_LE(b, 3);
  return c_[b];
}

template <typename VType>
VType Vector4<VType>::operator[](const int b) const {
  DCHECK_GE(b, 0);
  DCHECK_LE(b, 3);
  return c_[b];
}

template <typename VType>
void Vector4<VType>::x(const VType &v) {
  c_[0] = v;
}

template <typename VType>
VType Vector4<VType>::x() const {
  return c_[0];
}

template <typename VType>
void Vector4<VType>::y(const VType &v) {
  c_[1] = v;
}

template <typename VType>
VType Vector4<VType>::y() const {
  return c_[1];
}

template <typename VType>
void Vector4<VType>::z(const VType &v) {
  c_[2] = v;
}

template <typename VType>
VType Vector4<VType>::z() const {
  return c_[2];
}

template <typename VType>
void Vector4<VType>::w(const VType &v) {
  c_[3] = v;
}

template <typename VType>
VType Vector4<VType>::w() const {
  return c_[3];
}

template <typename VType>
VType* Vector4<VType>::Data() {
  return reinterpret_cast<VType*>(c_);
}

template <typename VType>
const VType* Vector4<VType>::Data() const {
  return reinterpret_cast<const VType*>(c_);
}

template <typename VType>
VType Vector4<VType>::Norm2(void) const {
  return c_[0]*c_[0] + c_[1]*c_[1] + c_[2]*c_[2] + c_[3]*c_[3];
}

template <typename VType>
typename Vector4<VType>::FloatType Vector4<VType>::Norm(void) const {
  return sqrt(Norm2());
}

template <typename VType>
Vector4<VType> Vector4<VType>::Normalize() const {
  COMPILE_ASSERT(!std::is_integral<VType>::value, must_be_floating_point);
  VType n = Norm();
  if (n != VType(0.0)) {
    n = VType(1.0) / n;
  }
  return Vector4(*this) *= n;
}

template <typename VType>
Vector4<VType> Vector4<VType>::Sqrt() const {
  return Vector4(sqrt(c_[0]), sqrt(c_[1]), sqrt(c_[2]), sqrt(c_[3]));
}

template <typename VType>
Vector4<VType> Vector4<VType>::Fabs() const {
  return Vector4(fabs(c_[0]), fabs(c_[1]), fabs(c_[2]), fabs(c_[3]));
}

template <typename VType>
Vector4<VType> Vector4<VType>::Abs() const {
  COMPILE_ASSERT(std::is_integral<VType>::value, use_Fabs_for_float_types);
  COMPILE_ASSERT(static_cast<VType>(-1) == -1, type_must_be_signed);
  COMPILE_ASSERT(sizeof(c_[0]) <= sizeof(int), Abs_truncates_to_int);
  return Vector4(abs(c_[0]), abs(c_[1]), abs(c_[2]), abs(c_[3]));
}

template <typename VType>
Vector4<VType> Vector4<VType>::Floor() const {
  return Vector4(floor(c_[0]),
                 floor(c_[1]),
                 floor(c_[2]),
                 floor(c_[3]));
}

template <typename VType>
Vector4<VType> Vector4<VType>::Ceil() const {
  return Vector4(ceil(c_[0]), ceil(c_[1]), ceil(c_[2]), ceil(c_[3]));
}

template <typename VType>
Vector4<VType> Vector4<VType>::FRound() const {
  return Vector4(rint(c_[0]), rint(c_[1]),
                 rint(c_[2]), rint(c_[3]));
}

template <typename VType>
Vector4<int> Vector4<VType>::IRound() const {
  return Vector4<int>(lrint(c_[0]), lrint(c_[1]),
                      lrint(c_[2]), lrint(c_[3]));
}

template <typename VType>
void Vector4<VType>::Clear() {
  c_[3] = c_[2] = c_[1] = c_[0] = VType();
}

template <typename VType>
bool Vector4<VType>::IsNaN() const {
  return isnan(c_[0]) || isnan(c_[1]) || isnan(c_[2]) || isnan(c_[3]);
}

template <typename VType>
Vector4<VType> Vector4<VType>::NaN() {
  return Vector4(MathUtil::NaN(), MathUtil::NaN(),
              MathUtil::NaN(), MathUtil::NaN());
}

template <typename VType>
Vector4<VType> operator-(const Vector4<VType> &vb) {
  return Vector4<VType>( -vb[0], -vb[1], -vb[2], -vb[3]);
}

template <typename ScalarType, typename VType>
Vector4<VType> operator*(const ScalarType k, const Vector4<VType> &v) {
  return Vector4<VType>(k*v[0], k*v[1], k*v[2], k*v[3]);
}

template <typename ScalarType, typename VType>
Vector4<VType> operator/(const ScalarType k, const Vector4<VType> &v) {
  return Vector4<VType>(k/v[0], k/v[1], k/v[2], k/v[3]);
}

template <typename VType>
Vector4<VType> Max(const Vector4<VType> &v1, const Vector4<VType> &v2) {
  return Vector4<VType>(std::max(v1[0], v2[0]), std::max(v1[1], v2[1]),
                        std::max(v1[2], v2[2]), std::max(v1[3], v2[3]));
}

template <typename VType>
Vector4<VType> Min(const Vector4<VType> &v1, const Vector4<VType> &v2) {
  return Vector4<VType>(std::min(v1[0], v2[0]), std::min(v1[1], v2[1]),
                        std::min(v1[2], v2[2]), std::min(v1[3], v2[3]));
}

template <typename VType>
std::ostream &operator <<(std::ostream &out, const Vector4<VType> &va) {
  out << "["
      << va[0] << ", "
      << va[1] << ", "
      << va[2] << ", "
      << va[3] << "]";
  return out;
}

template <>
inline std::ostream &operator <<(std::ostream &out, const Vector4<uint8> &va) {
  // ostream << uint8 prints the ASCII character, which is not useful.
  // Cast to int so that numbers will be printed instead.
  out << Vector4<int>::Cast(va);
  return out;
}

typedef Vector4<uint8>  Vector4_b;
typedef Vector4<int>    Vector4_i;
typedef Vector4<float>  Vector4_f;
typedef Vector4<double> Vector4_d;

// TODO(user): Vector4<T> does not actually satisfy the definition of a POD
// type even when T is a POD. Pretending that Vector4<T> is a POD probably
// won't cause any immediate problems, but eventually this should be fixed.


#endif  // UTIL_MATH_VECTOR4_H__

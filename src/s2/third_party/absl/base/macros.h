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

//
// Various Google-specific macros.
//
// This code is compiled directly on many platforms, including client
// platforms like Windows, Mac, and embedded systems.  Before making
// any changes here, make sure that you're not breaking any platforms.

#ifndef S2_THIRD_PARTY_ABSL_BASE_MACROS_H_
#define S2_THIRD_PARTY_ABSL_BASE_MACROS_H_

#include <cstddef>         // For size_t
#include "s2/third_party/absl/base/port.h"

// ABSTRACT is equivalent to "= 0" in most builds.  Expands to {} in SWIG
// builds only.  Avoid this macro in APIs that are not SWIG wrapped.
#ifndef SWIG
#define ABSTRACT = 0
#endif

// COMPILE_ASSERT is provided for backwards compatability.  New code
// should use static_assert instead.
#define COMPILE_ASSERT(expr, msg) static_assert((expr), #msg)

// DISALLOW_COPY_AND_ASSIGN disallows the copy constructor and copy assignment
// operator. DISALLOW_IMPLICIT_CONSTRUCTORS is like DISALLOW_COPY_AND_ASSIGN,
// but also disallows the default constructor, intended to help make a
// class uninstantiable.
//
// Unless building with a pre-C++11 compiler, prefer the language-supported
// "= delete" syntax in the public: section of the class.
//
// These must be placed in the private: declarations for a class,
// as a technical requirement of the non-C++11 branch.
#if LANG_CXX11
#define DISALLOW_COPY_AND_ASSIGN(TypeName) \
  TypeName(const TypeName&) = delete;      \
  TypeName& operator=(const TypeName&) = delete
#define DISALLOW_IMPLICIT_CONSTRUCTORS(TypeName) \
  TypeName() = delete;                           \
  DISALLOW_COPY_AND_ASSIGN(TypeName)
#else  // C++98 case follows
// Note that these C++98 implementations cannot completely disallow copying,
// as members and friends can still accidentally make elided copies without
// triggering a linker error.
#define DISALLOW_COPY_AND_ASSIGN(TypeName) \
  TypeName(const TypeName&);               \
  TypeName& operator=(const TypeName&)
#define DISALLOW_IMPLICIT_CONSTRUCTORS(TypeName) \
  TypeName();                                    \
  DISALLOW_COPY_AND_ASSIGN(TypeName)
#endif  // LANG_CXX11

// The arraysize(arr) macro returns the # of elements in an array arr.
// The expression is a compile-time constant, and therefore can be
// used in defining new arrays, for example.  If you use arraysize on
// a pointer by mistake, you will get a compile-time error.
//
// This template function declaration is used in defining arraysize.
// Note that the function doesn't need an implementation, as we only
// use its type.
template <typename T, size_t N>
char (&ArraySizeHelper(T (&array)[N]))[N];

// That gcc wants both of these prototypes seems mysterious. VC, for
// its part, can't decide which to use (another mystery). Matching of
// template overloads: the final frontier.
#ifndef _MSC_VER
template <typename T, size_t N>
char (&ArraySizeHelper(const T (&array)[N]))[N];
#endif

#define arraysize(array) (sizeof(ArraySizeHelper(array)))

// DEPRECATED: arraysize() is prefered over ARRAYSIZE().
//
// ARRAYSIZE performs essentially the same calculation as arraysize,
// but can be used on anonymous types or types defined inside
// functions.  It's less safe than arraysize as it accepts some
// (although not all) pointers.  Therefore, you should use arraysize
// whenever possible.
//
// The expression ARRAYSIZE(a) is a compile-time constant of type
// size_t.
//
// ARRAYSIZE catches a few type errors.  If you see a compiler error
//
//   "warning: division by zero in ..."
//
// when using ARRAYSIZE, you are (wrongfully) giving it a pointer.
// You should only use ARRAYSIZE on statically allocated arrays.
//
// The following comments are on the implementation details, and can
// be ignored by the users.
//
// ARRAYSIZE(arr) works by inspecting sizeof(arr) (the # of bytes in
// the array) and sizeof(*(arr)) (the # of bytes in one array
// element).  If the former is divisible by the latter, perhaps arr is
// indeed an array, in which case the division result is the # of
// elements in the array.  Otherwise, arr cannot possibly be an array,
// and we generate a compiler error to prevent the code from
// compiling.
//
// Since the size of bool is implementation-defined, we need to cast
// !(sizeof(a) & sizeof(*(a))) to size_t in order to ensure the final
// result has type size_t.
//
// This macro is not perfect as it wrongfully accepts certain
// pointers, namely where the pointer size is divisible by the pointee
// size.  For code that goes through a 32-bit compiler, where a pointer
// is 4 bytes, this means all pointers to a type whose size is 3 or
// greater than 4 will be (righteously) rejected.
//
// Kudos to Jorg Brown for this simple and elegant implementation.
//
// - wan 2005-11-16
//
// Starting with Visual C++ 2005, WinNT.h includes ARRAYSIZE.
#if !defined(_MSC_VER)
#define ARRAYSIZE(a) \
  ((sizeof(a) / sizeof(*(a))) / \
   static_cast<size_t>(!(sizeof(a) % sizeof(*(a)))))
#endif

// A macro to turn a symbol into a string
// These should be used judiciously in modern code. Prefer raw-string literal
// should you simply need to specify a large string.
#define AS_STRING(x)   AS_STRING_INTERNAL(x)
#define AS_STRING_INTERNAL(x)   #x

// The following enum should be used only as a constructor argument to indicate
// that the variable has static storage class, and that the constructor should
// do nothing to its state.  It indicates to the reader that it is legal to
// declare a static instance of the class, provided the constructor is given
// the base::LINKER_INITIALIZED argument.  Normally, it is unsafe to declare a
// static variable that has a constructor or a destructor because invocation
// order is undefined.  However, IF the type can be initialized by filling with
// zeroes (which the loader does for static variables), AND the type's
// destructor does nothing to the storage, then a constructor for static
// initialization can be declared as
//       explicit MyClass(base::LinkerInitialized x) {}
// and invoked as
//       static MyClass my_variable_name(base::LINKER_INITIALIZED);
namespace base {
enum LinkerInitialized { LINKER_INITIALIZED };
}

// The FALLTHROUGH_INTENDED macro can be used to annotate implicit fall-through
// between switch labels:
//  switch (x) {
//    case 40:
//    case 41:
//      if (truth_is_out_there) {
//        ++x;
//        FALLTHROUGH_INTENDED;  // Use instead of/along with annotations in
//                               // comments.
//      } else {
//        return x;
//      }
//    case 42:
//      ...
//
//  As shown in the example above, the FALLTHROUGH_INTENDED macro should be
//  followed by a semicolon. It is designed to mimic control-flow statements
//  like 'break;', so it can be placed in most places where 'break;' can, but
//  only if there are no statements on the execution path between it and the
//  next switch label.
//
//  When compiled with clang in C++11 mode, the FALLTHROUGH_INTENDED macro is
//  expanded to [[clang::fallthrough]] attribute, which is analysed when
//  performing switch labels fall-through diagnostic ('-Wimplicit-fallthrough').
//  See clang documentation on language extensions for details:
//  http://clang.llvm.org/docs/AttributeReference.html#fallthrough-clang-fallthrough
//
//  When used with unsupported compilers, the FALLTHROUGH_INTENDED macro has no
//  effect on diagnostics.
//
//  In either case this macro has no effect on runtime behavior and performance
//  of code.
#if defined(__clang__) && defined(LANG_CXX11) && defined(__has_warning)
#if __has_feature(cxx_attributes) && __has_warning("-Wimplicit-fallthrough")
#define FALLTHROUGH_INTENDED [[clang::fallthrough]]  // NOLINT
#endif
#endif

#ifndef FALLTHROUGH_INTENDED
#define FALLTHROUGH_INTENDED do { } while (0)
#endif

// The ABSL_DEPRECATED(...) macro can be used to mark deprecated class,
// struct, enum, function, method and variable declarations. The macro argument
// is used as a custom diagnostic message (e.g. suggestion of a better
// alternative):
//
//   class ABSL_DEPRECATED("Use Bar instead") Foo {...};
//   ABSL_DEPRECATED("Use Baz instead") void Bar() {...}
//
// Every usage of a deprecated entity will trigger a warning when compiled with
// clang's -Wdeprecated-declarations option. This option is turned off by
// default, but the warnings will be reported by go/clang-tidy.
#if defined(__clang__) && __cplusplus >= 201103L && defined(__has_warning)
#define ABSL_DEPRECATED(message) __attribute__((deprecated(message)))  // NOLINT
#endif

#ifndef ABSL_DEPRECATED
#define ABSL_DEPRECATED(message)
#endif

// ABCL_DEPRECATED will be replaced by ABSL_DEPRECATED once all usage are
// replaced.
#if defined(__clang__) && __cplusplus >= 201103L && defined(__has_warning)
#define ABCL_DEPRECATED(message) __attribute__((deprecated(message)))  // NOLINT
#endif

#ifndef ABCL_DEPRECATED
#define ABCL_DEPRECATED(message)
#endif

// GOOGLE_DEPRECATED will be replaced by ABSL_DEPRECATED once all usage are
// replaced.
#if defined(__clang__) && defined(LANG_CXX11) && defined(__has_warning)
#define GOOGLE_DEPRECATED(message) \
  __attribute__((deprecated(message)))  // NOLINT
#endif

#ifndef GOOGLE_DEPRECATED
#define GOOGLE_DEPRECATED(message)
#endif

// The CLANG_WARN_UNUSED_RESULT macro can be used on a class or struct to mark
// the class or struct as one that should never be ignored when it is a return
// value. Which is to say, any function or method that returns an instance of
// the marked class/struct and is not a member function of the class/struct
// will be treated as if it has the "warn_unused_result" attribute. The macro
// is only expanded into an attribute annotation when compiling with clang, as
// the use of the "warn_unused_result" attribute on a class or struct is a
// clang-specific extension of the eponymous function attribute.
//
// For example, in:
//   class CLANG_WARN_UNUSED_RESULT Status {
//     ...
//     void CheckSuccess();
//     Status StripMessage() const;
//   };
//
//   Status CreateResource();
//
//   void DoIt() {
//     Status s = CreateResource();
//     s.StripMessage();
//     CreateResource();
//     CreateResource().CheckSuccess();
//   }
//
// The first call to CreateResource in DoIt will not trigger a warning because
// the returned Status object was assigned to a variable. The call to
// Status::StripMessage also won't raise a warning despite the returned Status
// object being unused because the method is a member of the Status class.
// The second call to CreateResource will raise a warning because CreateResource
// returns a Status object and that object is unused (even though CreateResource
// was not explicitly declared with the "warn_unused_result" attribute). The
// final call to CreateResource is fine since the CheckSuccess method is called
// for the returned Status object.
#if defined(__clang__)
# if defined(LANG_CXX11) && __has_feature(cxx_attributes)
#  define CLANG_WARN_UNUSED_RESULT [[clang::warn_unused_result]]  // NOLINT
# else
#  define CLANG_WARN_UNUSED_RESULT __attribute__((warn_unused_result))  // NOLINT
# endif
#else
# define CLANG_WARN_UNUSED_RESULT
#endif

// The CLANG_BAD_CALL_IF macro can be used on a function overload to trap
// bad calls: any call that matches the overload will cause a compile-time
// error.  This uses a clang-specific "enable_if" attribute, as described at
// http://clang.llvm.org/docs/AttributeReference.html#enable-if
//
// Overloads which use this macro should be surrounded by
// "#ifdef CLANG_BAD_CALL_IF".  For example:
//
// int isdigit(int c);
// #ifdef CLANG_BAD_CALL_IF
// int isdigit(int c)
//     CLANG_BAD_CALL_IF(c <= -1 || c > 255,
//                       "'c' must have the value of an unsigned char or EOF");
// #endif // CLANG_BAD_CALL_IF
#if defined(__clang__)
# if __has_attribute(enable_if)
#  define CLANG_BAD_CALL_IF(expr, msg) \
    __attribute__((enable_if(expr, "Bad call trap"), unavailable(msg)))
# endif
#endif

#endif  // S2_THIRD_PARTY_ABSL_BASE_MACROS_H_

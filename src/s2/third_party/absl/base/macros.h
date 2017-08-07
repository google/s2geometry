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

#include <cstddef>

#include "s2/third_party/absl/base/port.h"

// The ABSL_ARRAYSIZE(arr) macro returns the # of elements in an array arr.
// The expression is a compile-time constant, and therefore can be
// used in defining new arrays, for example.  If you use arraysize on
// a pointer by mistake, you will get a compile-time error.
//
// This template function declaration is used in defining arraysize.
// Note that the function doesn't need an implementation, as we only
// use its type.
template <typename T, size_t N>
char (&ArraySizeHelper(T (&array)[N]))[N];

// TODO(b/62370839): Replace arraysize() with ABSL_ARRAYSIZE().
#define arraysize(array) (sizeof(ArraySizeHelper(array)))
#define ABSL_ARRAYSIZE(array) (sizeof(ArraySizeHelper(array)))

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

// The ABSL_FALLTHROUGH_INTENDED macro can be used to annotate implicit
// fall-through between switch labels:
//  switch (x) {
//    case 40:
//    case 41:
//      if (truth_is_out_there) {
//        ++x;
//        ABSL_FALLTHROUGH_INTENDED;  // Use instead of/along with annotations
//                                    // in comments
//      } else {
//        return x;
//      }
//    case 42:
//      ...
//
//  As shown in the example above, the ABSL_FALLTHROUGH_INTENDED macro should be
//  followed by a semicolon. It is designed to mimic control-flow statements
//  like 'break;', so it can be placed in most places where 'break;' can, but
//  only if there are no statements on the execution path between it and the
//  next switch label.
//
//  When compiled with clang in C++11 mode, the ABSL_FALLTHROUGH_INTENDED macro
//  is expanded to [[clang::fallthrough]] attribute, which is analysed when
//  performing switch labels fall-through diagnostic ('-Wimplicit-fallthrough').
//  See clang documentation on language extensions for details:
//  http://clang.llvm.org/docs/AttributeReference.html#fallthrough-clang-fallthrough
//
//  When used with unsupported compilers, the ABSL_FALLTHROUGH_INTENDED macro
//  has no effect on diagnostics.
//
//  In either case this macro has no effect on runtime behavior and performance
//  of code.
// TODO(b/62370839): Replace FALLTHROUGH_INTENDED with
// ABSL_FALLTHROUGH_INTENDED.
#if defined(__clang__) && defined(__has_warning)
#if __has_feature(cxx_attributes) && __has_warning("-Wimplicit-fallthrough")
#define FALLTHROUGH_INTENDED [[clang::fallthrough]]  // NOLINT
#endif
#endif

#ifndef FALLTHROUGH_INTENDED
#define FALLTHROUGH_INTENDED do { } while (0)
#endif
#ifdef ABSL_FALLTHROUGH_INTENDED
#error "ABSL_FALLTHROUGH_INTENDED should not be defined."
#endif

#if defined(__clang__) && defined(__has_warning)
#if __has_feature(cxx_attributes) && __has_warning("-Wimplicit-fallthrough")
#define ABSL_FALLTHROUGH_INTENDED [[clang::fallthrough]]  // NOLINT
#endif
#endif

#ifndef ABSL_FALLTHROUGH_INTENDED
#define ABSL_FALLTHROUGH_INTENDED \
  do {                            \
  } while (0)
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

// The ABSL_BAD_CALL_IF macro can be used on a function overload to trap
// bad calls: any call that matches the overload will cause a compile-time
// error.  This uses a clang-specific "enable_if" attribute, as described at
// http://clang.llvm.org/docs/AttributeReference.html#enable-if
//
// Overloads which use this macro should be surrounded by
// "#ifdef ABSL_BAD_CALL_IF".  For example:
//
// int isdigit(int c);
// #ifdef ABSL_BAD_CALL_IF
// int isdigit(int c)
//     ABSL_BAD_CALL_IF(c <= -1 || c > 255,
//                       "'c' must have the value of an unsigned char or EOF");
// #endif // ABSL_BAD_CALL_IF

#if defined(__clang__)
# if __has_attribute(enable_if)
#  define ABSL_BAD_CALL_IF(expr, msg) \
    __attribute__((enable_if(expr, "Bad call trap"), unavailable(msg)))
# endif
#endif

#endif  // S2_THIRD_PARTY_ABSL_BASE_MACROS_H_

// Copyright 2016 Google Inc. All Rights Reserved.
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

// Temporary header providing transition support for the gtl namespace

#ifndef S2GEOMETRY_UTIL_GTL_GTL_NAMESPACE_H_
#define S2GEOMETRY_UTIL_GTL_GTL_NAMESPACE_H_

namespace util {
namespace gtl {
}
}

// Namespace alias to permit user code to begin referring to the ::gtl
// namespace.
namespace gtl = ::util::gtl;

// Macros for declaring the gtl namespace. Code can use these macros instead
// of explicitly naming the namespace in order to be oblivious to the eventual
// switch to make ::gtl the canonical name of the namespace.
#define GTL_NAMESPACE_BEGIN \
  namespace util {          \
  namespace gtl {

#define GTL_NAMESPACE_END }}

#define GTL_NAMESPACE util::gtl

#endif  // S2GEOMETRY_UTIL_GTL_GTL_NAMESPACE_H_

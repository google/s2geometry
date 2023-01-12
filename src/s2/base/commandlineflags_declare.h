// Copyright Google Inc. All Rights Reserved.
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

#ifndef S2_BASE_COMMANDLINEFLAGS_DECLARE_H_
#define S2_BASE_COMMANDLINEFLAGS_DECLARE_H_

#include <string>

#include "s2/base/integral_types.h"

#ifdef S2_USE_GFLAGS

#include <gflags/gflags.h>

#define S2_DECLARE_bool(name) DECLARE_bool(name)
#define S2_DECLARE_double(name) DECLARE_double(name)
#define S2_DECLARE_int32(name) DECLARE_int32(name)
#define S2_DECLARE_int64(name) DECLARE_int64(name)
#define S2_DECLARE_string(name) DECLARE_string(name)

#else  // !defined(S2_USE_GFLAGS)

#define S2_DECLARE_bool(name) \
  extern bool FLAGS_##name

#define S2_DECLARE_double(name) \
  extern double FLAGS_##name

#define S2_DECLARE_int32(name) \
  extern int32 FLAGS_##name

#define S2_DECLARE_int64(name) \
  extern int64 FLAGS_##name

#define S2_DECLARE_string(name) \
  extern std::string FLAGS_##name

#endif  // !defined(S2_USE_GFLAGS)

#endif  // S2_BASE_COMMANDLINEFLAGS_DECLARE_H_

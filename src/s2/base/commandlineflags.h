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

#ifndef S2_BASE_COMMANDLINEFLAGS_H_
#define S2_BASE_COMMANDLINEFLAGS_H_

#include <string>

#include "absl/flags/flag.h"

#include "s2/base/commandlineflags_declare.h"
#include "s2/base/integral_types.h"

#define S2_DEFINE_bool(name, default_value, description) \
  ABSL_FLAG(bool, name, default_value, description)

#define S2_DEFINE_double(name, default_value, description) \
  ABSL_FLAG(double, name, default_value, description)

#define S2_DEFINE_int32(name, default_value, description) \
  ABSL_FLAG(int32, name, default_value, description)

#define S2_DEFINE_int64(name, default_value, description) \
  ABSL_FLAG(int64, name, default_value, description)

#define S2_DEFINE_string(name, default_value, description) \
  ABSL_FLAG(std::string, name, default_value, description)

#endif  // S2_BASE_COMMANDLINEFLAGS_H_

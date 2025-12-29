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

#ifndef S2_UTIL_TASK_STATUS_MACROS_H_
#define S2_UTIL_TASK_STATUS_MACROS_H_

#include <utility>

// Executes an expression that returns a StatusOr, extracting its value
// into the variable defined by lhs (or returning the status on error).
//
// Example:
//   absl::StatusOr<uint64_t> ParseId(absl::string_view id_string);
//   ASSIGN_OR_RETURN(uint64_t id, ParseId(id_string));
#define ASSIGN_OR_RETURN(lhs, rhs)                           \
  auto status_or = (rhs);                                    \
  if (!status_or.ok()) return std::move(status_or).status(); \
  lhs = *std::move(status_or)

#endif  // S2_UTIL_TASK_STATUS_MACROS_H_

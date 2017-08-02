// Copyright 2017 Google Inc. All Rights Reserved.
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

#include "s2/third_party/absl/strings/match.h"

#include "s2/third_party/absl/strings/internal/memutil.h"

namespace absl {

namespace {
// Forked from string_view_utils. Should re-unify with numbers.cc at
// some point.
bool CaseEqual(absl::string_view piece1, absl::string_view piece2) {
  return (piece1.size() == piece2.size() &&
          0 == absl::strings_internal::memcasecmp(piece1.data(), piece2.data(),
                                                  piece1.size()));
  // memcasecmp uses ascii_tolower().
}
}  // namespace

bool StartsWithIgnoreCase(absl::string_view text,
                          absl::string_view starts_with) {
  if (text.size() < starts_with.size()) return false;
  return CaseEqual(text.substr(0, starts_with.size()), starts_with);
}

bool EndsWithIgnoreCase(absl::string_view text, absl::string_view ends_with) {
  if (text.size() < ends_with.size()) return false;
  return CaseEqual(text.substr(text.size() - ends_with.size()), ends_with);
}

}  // namespace absl

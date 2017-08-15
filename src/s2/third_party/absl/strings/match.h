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

#ifndef S2_THIRD_PARTY_ABSL_STRINGS_MATCH_H_
#define S2_THIRD_PARTY_ABSL_STRINGS_MATCH_H_

// String matching utilities.

#include "s2/third_party/absl/strings/string_view.h"

namespace absl {

// Returns whether s contains x.
inline bool StrContains(absl::string_view s, absl::string_view x) {
  return static_cast<absl::string_view::size_type>(s.find(x, 0)) != s.npos;
}

// Returns whether s begins with x.
inline bool StartsWith(absl::string_view s, absl::string_view x) {
  return x.empty() ||
         (s.size() >= x.size() &&
          absl::strings_internal::memeq(s.data(), x.data(), x.size()));
  /* absl:oss-replace-with
  return x.empty() ||
         (s.size() >= x.size() && memcmp(s.data(), x.data(), x.size()) == 0);
  absl:oss-replace-end */
}

// Returns whether s ends with x.
inline bool EndsWith(absl::string_view s, absl::string_view x) {
  return x.empty() ||
         (s.size() >= x.size() &&
          absl::strings_internal::memeq(s.data() + (s.size() - x.size()),
                                        x.data(), x.size()));
  /* absl:oss-replace-with
  return x.empty() ||
         (s.size() >= x.size() &&
          memcmp(s.data() + (s.size() - x.size()), x.data(), x.size()) == 0);
  absl:oss-replace-end */
}

// Returns true if "text" starts with "starts_with". The comparison ignores
// case.
bool StartsWithIgnoreCase(absl::string_view text,
                          absl::string_view starts_with);

// Returns true if "text" ends with "ends_with". The comparison ignores
// case.
bool EndsWithIgnoreCase(absl::string_view text, absl::string_view ends_with);

}  // namespace absl

#endif  // S2_THIRD_PARTY_ABSL_STRINGS_MATCH_H_

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

#include "s2/third_party/absl/strings/string_view_utils.h"

namespace absl {

bool StartsWithIgnoreCase(absl::string_view text,
                          absl::string_view starts_with) {
  if (text.size() < starts_with.size()) return false;
  return strings::EqualIgnoreCase(text.substr(0, starts_with.size()),
                                  starts_with);
}

bool EndsWithIgnoreCase(absl::string_view text, absl::string_view ends_with) {
  if (text.size() < ends_with.size()) return false;
  return strings::EqualIgnoreCase(text.substr(text.size() - ends_with.size()),
                                  ends_with);
}

}  // namespace absl

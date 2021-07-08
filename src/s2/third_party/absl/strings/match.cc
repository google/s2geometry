// Copyright 2017 The Abseil Authors.
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//      http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

#include "s2/third_party/absl/strings/match.h"

#include "s2/third_party/absl/strings/internal/memutil.h"

namespace s2::absl {

bool EqualsIgnoreCase(::s2::absl::string_view piece1, ::s2::absl::string_view piece2) {
  return (piece1.size() == piece2.size() &&
          0 == ::s2::absl::strings_internal::memcasecmp(piece1.data(), piece2.data(),
                                                  piece1.size()));
  // memcasecmp uses ::s2::absl::ascii_tolower().
}

bool StartsWithIgnoreCase(::s2::absl::string_view text, ::s2::absl::string_view prefix) {
  return (text.size() >= prefix.size()) &&
         EqualsIgnoreCase(text.substr(0, prefix.size()), prefix);
}

bool EndsWithIgnoreCase(::s2::absl::string_view text, ::s2::absl::string_view suffix) {
  return (text.size() >= suffix.size()) &&
         EqualsIgnoreCase(text.substr(text.size() - suffix.size()), suffix);
}

}  // namespace s2::absl

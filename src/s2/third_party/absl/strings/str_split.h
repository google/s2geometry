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

#ifndef S2_THIRD_PARTY_ABSL_STRINGS_STR_SPLIT_H_
#define S2_THIRD_PARTY_ABSL_STRINGS_STR_SPLIT_H_

#include <functional>
#include <string>
#include <vector>

#include "s2/third_party/absl/strings/strip.h"

namespace strings {

std::vector<std::string> Split(
    std::string const& text, char delim,
    std::function<bool(std::string const&)> predicate);
std::vector<std::string> Split(std::string const& text, char delim);

// Returns false if the given StringPiece is empty, indicating that the
// strings::Split() API should omit the empty string.
//
// std::vector<string> v = Split(" a , ,,b,", ',', SkipEmpty());
// EXPECT_THAT(v, ElementsAre(" a ", " ", "b"));
struct SkipEmpty {
  bool operator()(std::string const& s) const { return !s.empty(); }
};

// Returns false if the given string is empty or contains only whitespace,
// indicating that the strings::Split() API should omit the string.
//
// std::vector<string> v = Split(" a , ,,b,", ',', SkipWhitespace());
// EXPECT_THAT(v, ElementsAre(" a ", "b"));
struct SkipWhitespace {
  bool operator()(std::string s) const {
    StripWhitespace(&s);
    return !s.empty();
  }
};

}  // namespace strings

#endif  // S2_THIRD_PARTY_ABSL_STRINGS_STR_SPLIT_H_

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

#include "strings/split.h"

#include <functional>
#include <sstream>
#include <string>
#include <vector>

namespace strings {

std::vector<std::string> Split(
    std::string const& text, char const delim,
    std::function<bool(std::string const&)> predicate) {
  std::stringstream ss(text);
  std::string item;
  std::vector<std::string> elems;
  while (std::getline(ss, item, delim)) {
    if (predicate(item))
      elems.push_back(std::move(item));
  }
  return elems;
}

std::vector<std::string> Split(std::string const& text, char const delim) {
  return Split(text, delim, [](const std::string&) { return true; });
}

}  // namespace strings

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

#include "s2/third_party/absl/strings/str_split.h"

#include <functional>
#include <string>
#include <vector>

#include "s2/third_party/absl/strings/string_view.h"

using absl::string_view;

namespace strings {

template <typename String>
std::vector<String> Split(
    String const& text, char const delim,
    std::function<bool(string_view)> predicate) {
  std::vector<String> elems;
  typename String::size_type begin = 0;
  typename String::size_type end;
  while ((end = text.find(delim, begin)) != String::npos) {
    string_view view(text.data() + begin, end - begin);
    if (predicate(view))
      elems.emplace_back(view);
    begin = end + 1;
  }
  // Try to add the portion after the last delim.
  string_view view(text.data() + begin, text.size() - begin);
  if (predicate(view))
    elems.emplace_back(view);
  return elems;
}
template std::vector<std::string> Split(
    std::string const& text, char const delim,
    std::function<bool(string_view)> predicate);
template std::vector<string_view> Split(
    string_view const& text, char const delim,
    std::function<bool(string_view)> predicate);

template <typename String>
std::vector<String> Split(String const& text, char const delim) {
  return Split(text, delim, [](string_view) { return true; });
}
template std::vector<std::string> Split(std::string const& text,
                                        char const delim);
template std::vector<string_view> Split(string_view const& text,
                                        char const delim);

}  // namespace strings

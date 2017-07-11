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

#include "s2/strings/serialize.h"

#include <string>
#include <vector>

#include "s2/third_party/absl/strings/str_split.h"
#include "s2/third_party/absl/strings/string_view.h"

using absl::string_view;

bool DictionaryParse(string_view encoded_str,
                     std::vector<std::pair<std::string, std::string>>* items) {
  if (encoded_str.empty())
    return true;
  std::vector<string_view> const entries  = strings::Split(encoded_str, ',');
  for (int i = 0; i < entries.size(); ++i) {
    std::vector<string_view> const fields = strings::Split(entries[i], ':');
    if (fields.size() != 2)  // parsing error
      return false;
    items->push_back(std::make_pair(std::string(fields[0]),
                                    std::string(fields[1])));
  }
  return true;
}

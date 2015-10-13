// Copyright 2015 Google Inc. All Rights Reserved.
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

#include "strings/serialize.h"

#include <string>
#include <vector>

#include "strings/split.h"

bool DictionaryParse(std::string const& encoded_str,
                     std::vector<std::pair<std::string, std::string>>* items) {
  if (encoded_str.empty())
    return true;
  std::vector<std::string> const entries  = strings::Split(encoded_str, ',');
  for (int i = 0; i < entries.size(); ++i) {
    std::vector<std::string> const fields = strings::Split(entries[i], ':');
    if (fields.size() != 2)  // parsing error
      return false;
    items->push_back(std::make_pair(fields[0], fields[1]));
  }
  return true;
}

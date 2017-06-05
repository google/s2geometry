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

#ifndef S2_THIRD_PARTY_ABSL_STRINGS_STRIP_H_
#define S2_THIRD_PARTY_ABSL_STRINGS_STRIP_H_

#include <string>

#include "s2/third_party/absl/base/port.h"

// Removes whitespace from both ends of the given string.
void StripWhitespace(string* str);

namespace strings {

// Calls StripWhitespace() on each element in the given collection.
//
// Note: this implementation is conceptually similar to
//
//   std::for_each(c.begin(), c.end(), StripWhitespace);
//
// except that StripWhitespace requires a *pointer* to the element, so the above
// std::for_each solution wouldn't work.
template <typename Collection>
inline void StripWhitespaceInCollection(Collection* collection) {
  for (typename Collection::iterator it = collection->begin();
       it != collection->end(); ++it)
    StripWhitespace(&(*it));
}

}  // namespace strings

#endif  // S2_THIRD_PARTY_ABSL_STRINGS_STRIP_H_

#ifndef S2_THIRD_PARTY_ABSL_STRINGS_STR_JOIN_H_
#define S2_THIRD_PARTY_ABSL_STRINGS_STR_JOIN_H_

#include <string>

namespace strings {

template <typename Range>
std::string Join(const Range& strs, const char* delim) {
  std::string joined;
  bool first = true;
  for (const auto& s : strs) {
    if (!first) joined.append(delim);
    joined.append(s);
    first = false;
  }
  return joined;
}

}  // namespace strings

#endif  // S2_THIRD_PARTY_ABSL_STRINGS_STR_JOIN_H_

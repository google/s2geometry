#include "s2/third_party/absl/strings/str_split.h"

#include <functional>
#include <string>
#include <vector>

#include "s2/third_party/absl/strings/string_view.h"

using absl::string_view;
using std::function;
using std::string;
using std::vector;

namespace strings {

template <typename String>
vector<String> Split(String const& text, char const delim,
                     function<bool(string_view)> predicate) {
  vector<String> elems;
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
template vector<string> Split(string const& text, char const delim,
                              function<bool(string_view)> predicate);
template vector<string_view> Split(string_view const& text, char const delim,
                                   function<bool(string_view)> predicate);

template <typename String>
vector<String> Split(String const& text, char const delim) {
  return Split(text, delim, [](string_view) { return true; });
}
template vector<string> Split(string const& text, char const delim);
template vector<string_view> Split(string_view const& text, char const delim);

}  // namespace strings

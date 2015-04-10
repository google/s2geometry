// Copyright 2008 Google Inc. All Rights Reserved.
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

#include "strings/join.h"

#include <glog/logging.h>

#include "base/scoped_ptr.h"
#include "strings/ascii_ctype.h"
#include "strings/escaping.h"

using std::vector;

void JoinCSVLineWithDelimiter(const vector<string>& cols, char delimiter,
                              string* output) {
  CHECK(output);
  CHECK(output->empty());
  vector<string> quoted_cols;

  const string delimiter_str(1, delimiter);
  const string escape_chars = delimiter_str + "\"";

  // If the string contains the delimiter or " anywhere, or begins or ends with
  // whitespace (ie ascii_isspace() returns true), escape all double-quotes and
  // bracket the string in double quotes. string.rbegin() evaluates to the last
  // character of the string.
  for (int i = 0; i < cols.size(); ++i) {
    if ((cols[i].find_first_of(escape_chars) != string::npos) ||
        (!cols[i].empty() && (ascii_isspace(*cols[i].begin()) ||
                              ascii_isspace(*cols[i].rbegin())))) {
      // Double the original size, for escaping, plus two bytes for
      // the bracketing double-quotes, and one byte for the closing \0.
      int size = 2 * cols[i].size() + 3;
      scoped_ptr<char[]> buf(new char[size]);

      // Leave space at beginning and end for bracketing double-quotes.
      int escaped_size = strings::EscapeStrForCSV(cols[i].c_str(),
                                                  buf.get() + 1, size - 2);
      CHECK_GE(escaped_size, 0) << "Buffer somehow wasn't large enough.";
      CHECK_GE(size, escaped_size + 3)
        << "Buffer should have one space at the beginning for a "
        << "double-quote, one at the end for a double-quote, and "
        << "one at the end for a closing '\0'";
      *buf.get() = '"';
      *((buf.get() + 1) + escaped_size) = '"';
      *((buf.get() + 1) + escaped_size + 1) = '\0';
      quoted_cols.push_back(string(buf.get(), buf.get() + escaped_size + 2));
    } else {
      quoted_cols.push_back(cols[i]);
    }
  }
  *output = strings::Join(quoted_cols, delimiter_str);
}

void JoinCSVLine(const vector<string>& cols, string* output) {
  JoinCSVLineWithDelimiter(cols, ',', output);
}

string JoinCSVLine(const vector<string>& cols) {
  string output;
  JoinCSVLine(cols, &output);
  return output;
}

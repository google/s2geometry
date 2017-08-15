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

#include "s2/third_party/absl/base/log_severity.h"

namespace absl {
#ifdef NDEBUG
const absl::LogSeverity kLogDebugFatal = absl::LogSeverity::kError;
#else
const absl::LogSeverity kLogDebugFatal = absl::LogSeverity::kFatal;
#endif
}  // namespace absl

namespace base_logging {

const char* const LogSeverityNames[NUM_SEVERITIES] = {
  "INFO", "WARNING", "ERROR", "FATAL"
};

LogSeverity NormalizeSeverity(LogSeverity s) {
  if (s < INFO) {
    return INFO;
  }
  if (s > FATAL) {
    return ERROR;
  }
  return s;
}

}  // namespace base_logging

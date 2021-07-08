// Copyright 2017 The Abseil Authors.
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//      http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.
//
//
// NOTE FOR GOOGLERS:
//
// IWYU pragma: private, include "base/log_severity.h"

#ifndef S2_THIRD_PARTY_ABSL_BASE_INTERNAL_LOG_SEVERITY_H_
#define S2_THIRD_PARTY_ABSL_BASE_INTERNAL_LOG_SEVERITY_H_

#include <array>

#include "s2/third_party/absl/base/attributes.h"

namespace s2::absl {

// Four severity levels are defined.  Logging APIs should terminate the program
// when a message is logged at severity `kFatal`; the other levels have no
// special semantics.
enum class LogSeverity : int {
  kInfo = 0,
  kWarning = 1,
  kError = 2,
  kFatal = 3,
};

// Returns an iterable of all standard `::s2::absl::LogSeverity` values, ordered from
// least to most severe.
constexpr std::array<::s2::absl::LogSeverity, 4> LogSeverities() {
  return {{::s2::absl::LogSeverity::kInfo, ::s2::absl::LogSeverity::kWarning,
           ::s2::absl::LogSeverity::kError, ::s2::absl::LogSeverity::kFatal}};
}

// `::s2::absl::kLogDebugFatal` equals `::s2::absl::LogSeverity::kFatal` in debug builds
// (i.e. when `NDEBUG` is not defined) and `::s2::absl::LogSeverity::kError`
// otherwise.  It is extern to prevent ODR violations when compilation units
// with different build settings are linked together.
ABSL_CONST_INIT extern const ::s2::absl::LogSeverity kLogDebugFatal;

// Returns the all-caps string representation (e.g. "INFO") of the specified
// severity level if it is one of the normal levels and "UNKNOWN" otherwise.
constexpr const char* LogSeverityName(::s2::absl::LogSeverity s) {
  return s == ::s2::absl::LogSeverity::kInfo
             ? "INFO"
             : s == ::s2::absl::LogSeverity::kWarning
                   ? "WARNING"
                   : s == ::s2::absl::LogSeverity::kError
                         ? "ERROR"
                         : s == ::s2::absl::LogSeverity::kFatal ? "FATAL" : "UNKNOWN";
}

// Values less than `kInfo` normalize to `kInfo`; values greater than `kFatal`
// normalize to `kError` (**NOT** `kFatal`).
constexpr ::s2::absl::LogSeverity NormalizeLogSeverity(::s2::absl::LogSeverity s) {
  return s < ::s2::absl::LogSeverity::kInfo
             ? ::s2::absl::LogSeverity::kInfo
             : s > ::s2::absl::LogSeverity::kFatal ? ::s2::absl::LogSeverity::kError : s;
}
constexpr ::s2::absl::LogSeverity NormalizeLogSeverity(int s) {
  return NormalizeLogSeverity(static_cast<::s2::absl::LogSeverity>(s));
}

}  // namespace s2::absl

#endif  // S2_THIRD_PARTY_ABSL_BASE_INTERNAL_LOG_SEVERITY_H_

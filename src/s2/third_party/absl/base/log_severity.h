// Copyright 2007 Google Inc. All Rights Reserved.
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


#ifndef S2_THIRD_PARTY_ABSL_BASE_LOG_SEVERITY_H_
#define S2_THIRD_PARTY_ABSL_BASE_LOG_SEVERITY_H_

#include "s2/third_party/absl/base/attributes.h"
#include "s2/third_party/absl/base/port.h"

// WinGDI.h defines ERROR, undef to avoid conflict naming.
#ifdef _WIN32
#undef ERROR
#endif

namespace absl {

enum class LogSeverity : int {
  kInfo = 0,
  kWarning = 1,
  kError = 2,
  kFatal = 3,
};

// absl::kLogDebugFatal equals absl::LogSeverity::kFatal in debug builds (i.e.
// when NDEBUG is not defined) and absl::LogSeverity::kError otherwise.
// It is extern to prevent ODR violations when compilation units with different
// build settings are linked together.
ABSL_CONST_INIT extern const absl::LogSeverity kLogDebugFatal;

constexpr const char* LogSeverityName(absl::LogSeverity s) {
  return s == absl::LogSeverity::kInfo
             ? "INFO"
             : s == absl::LogSeverity::kWarning
                   ? "WARNING"
                   : s == absl::LogSeverity::kError
                         ? "ERROR"
                         : s == absl::LogSeverity::kFatal ? "FATAL" : "UNKNOWN";
}

// Note that out-of-range large severities normalize to kError, not kFatal.
constexpr absl::LogSeverity NormalizeLogSeverity(absl::LogSeverity s) {
  return s < absl::LogSeverity::kInfo
             ? absl::LogSeverity::kInfo
             : s > absl::LogSeverity::kFatal ? absl::LogSeverity::kError : s;
}
constexpr absl::LogSeverity NormalizeLogSeverity(int s) {
  return NormalizeLogSeverity(static_cast<absl::LogSeverity>(s));
}

}  // namespace absl

namespace base_logging {

// Variables of type LogSeverity are widely taken to lie in the range
// [0, NUM_SEVERITIES-1].  Be careful to preserve this assumption if
// you ever need to change their values or add a new severity.
typedef int LogSeverity;

// Severity-level definitions.  These are C++ constants representing
// the four severity levels INFO through FATAL.
//
// These severity constants can be used outside of a LOG statement, or
// in a variable-severity statement such as:
//
//   LOG(LEVEL(is_error ? base_logging::ERROR : base_logging::WARNING)) << ...;
//
// Another common use of these symbols is:
//
//   using base_logging::INFO;
//   testing::ScopedMockLog log;
//   EXPECT_CALL(log, Log(INFO, testing::_, "..."));
const LogSeverity INFO = 0;
const LogSeverity WARNING = 1;
const LogSeverity ERROR = 2;
const LogSeverity FATAL = 3;
const int NUM_SEVERITIES = 4;

// An array containing "INFO", "WARNING", "ERROR", "FATAL".
extern const char* const LogSeverityNames[NUM_SEVERITIES];

}  // namespace base_logging

namespace base_logging {

// Use base_logging::DFATAL to refer to the LogSeverity level
// "DFATAL" outside of a log statement.
#ifdef NDEBUG
const LogSeverity DFATAL = ERROR;
#else
const LogSeverity DFATAL = FATAL;
#endif

// If "s" is less than base_logging::INFO, returns base_logging::INFO.
// If "s" is greater than base_logging::FATAL, returns
// base_logging::ERROR.  Otherwise, returns "s".
LogSeverity NormalizeSeverity(LogSeverity s);

}  // namespace base_logging

// Special-purpose macros for use OUTSIDE //base/logging.h.  These are
// to support hygienic forms of other logging macros (e.g.,
// raw_logging.h).
#define BASE_LOG_SEVERITY_INFO (::base_logging::INFO)
#define BASE_LOG_SEVERITY_WARNING (::base_logging::WARNING)
#define BASE_LOG_SEVERITY_ERROR (::base_logging::ERROR)
#define BASE_LOG_SEVERITY_FATAL (::base_logging::FATAL)
#define BASE_LOG_SEVERITY_DFATAL (::base_logging::DFATAL)

#define BASE_LOG_SEVERITY_LEVEL(severity) \
  (::base_logging::NormalizeSeverity(severity))

#endif  // S2_THIRD_PARTY_ABSL_BASE_LOG_SEVERITY_H_

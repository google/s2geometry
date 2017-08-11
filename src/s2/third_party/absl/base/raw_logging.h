// Copyright 2006 Google Inc. All Rights Reserved.
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

//
// Moved into its own module by sanjay@google.com (Sanjay Ghemawat)
//
// Thread-safe logging routines that do not allocate any memory or
// acquire any locks, and can therefore be used by low-level memory
// allocation, synchronization, and signal-handling code.

#ifndef S2_THIRD_PARTY_ABSL_BASE_RAW_LOGGING_H_
#define S2_THIRD_PARTY_ABSL_BASE_RAW_LOGGING_H_

#include "s2/third_party/absl/base/log_severity.h"
#include "s2/third_party/absl/base/port.h"

// This is similar to LOG(severity) << format..., but
// * it is to be used ONLY by low-level modules that can't use normal LOG()
// * it is designed to be a low-level logger that does not allocate any
//   memory and does not need any locks, hence:
// * it logs straight and ONLY to STDERR w/o buffering
// * it uses an explicit printf-format and arguments list
// * it will silently chop off really long message strings
// Usage example:
//   RAW_LOG(ERROR, "Failed foo with %i: %s", status, error);
// This will print an almost standard log line like this to stderr only:
//   E0821 211317 file.cc:123] RAW: Failed foo with 22: bad_file

#if !defined(STRIP_LOG) || STRIP_LOG == 0
#define RAW_LOG(severity, ...)                                         \
  do {                                                                 \
    ::base_raw_logging::RawLog(BASE_LOG_SEVERITY_##severity, __FILE__, \
                               __LINE__, __VA_ARGS__);                 \
  } while (0)
#else
#define RAW_LOG(severity, ...)                                                \
  do {                                                                        \
    if (BASE_LOG_SEVERITY_##severity == BASE_LOG_SEVERITY_FATAL)              \
      ::base_raw_logging::RawLog(BASE_LOG_SEVERITY_FATAL, __FILE__, __LINE__, \
                                 __VA_ARGS__);                                \
  } while (0)
#endif

// Similar to CHECK(condition) << message, but for low-level modules:
// we use only RAW_LOG that does not allocate memory.
// We do not want to provide args list here to encourage this usage:
//   if (!cond)  RAW_LOG(FATAL, "foo ...", hard_to_compute_args);
// so that the args are not computed when not needed.
#define RAW_CHECK(condition, message) \
  do { \
    if (ABSL_PREDICT_FALSE(!(condition))) {                            \
      RAW_LOG(FATAL, "Check %s failed: %s", #condition, message);      \
    } \
  } while (0)

// Debug versions of RAW_LOG and RAW_CHECK
#ifndef NDEBUG

#define RAW_DLOG(severity, ...) RAW_LOG(severity, __VA_ARGS__)
#define RAW_DCHECK(condition, message) RAW_CHECK(condition, message)

#else  // NDEBUG

#define RAW_DLOG(severity, ...) \
  while (false) \
    RAW_LOG(severity, __VA_ARGS__)
#define RAW_DCHECK(condition, message) \
  while (false) \
    RAW_CHECK(condition, message)

#endif  // NDEBUG

namespace base_raw_logging {

// Helper function to implement RAW_LOG
// Logs format... at "severity" level, reporting it
// as called from file:line.
// This does not allocate memory or acquire locks.
void RawLog(base_logging::LogSeverity severity, const char* file, int line,
            const char* format, ...) ABSL_PRINTF_ATTRIBUTE(4, 5);

// Function type for a raw_logging customization hook for suppressing messages
// by severity, and for writing custom prefixes on non-suppressed messages.
//
// The installed hook is called for every raw log invocation.  The message will
// be logged to stderr only if the hook returns true.  FATAL errors will cause
// the process to abort, even if writing to stderr is suppressed.  The hook is
// also provided with an output buffer, where it can write a custom log message
// prefix.
//
// The raw_logging system does not allocate memory or grab locks.  User-provided
// hooks must avoid these operations, and must not throw exceptions.
//
// 'severity' is the severity level of the message being written.
// 'file' and 'line' are the file and line number where the RAW_LOG macro was
// located.
// 'buffer' and 'buf_size' are pointers to the buffer and buffer size.  If the
// hook writes a prefix, it must increment *buffer and decrement *buf_size
// accordingly.
using LogPrefixHook = bool (*)(base_logging::LogSeverity severity,
                               const char* file, int line, char** buffer,
                               int* buf_size);

// Function type for a raw_logging customization hook called to abort a process
// when a FATAL message is logged.  If the provided AbortHook() returns, the
// logging system will call abort().
//
// 'file' and 'line' are the file and line number where the RAW_LOG macro was
// located.
// The null-terminated logged message lives in the buffer between 'buf_start'
// and 'buf_end'.  'prefix_end' points to the first non-prefix character of the
// buffer (as written by the LogPrefixHook.)
using AbortHook = void (*)(const char* file, int line, const char* buf_start,
                           const char* prefix_end, const char* buf_end);

// Registers hooks of the above types.  Only a single hook of each type may be
// registered.  It is an error to call these functions multiple times with
// different input arguments.
//
// These functions are safe to call at any point during initialization; they do
// not block or malloc, and are async-signal safe.
void RegisterLogPrefixHook(LogPrefixHook fn);
void RegisterAbortHook(AbortHook fn);

namespace internal {  // For testing.
// Returns true if raw logging is fully supported. When it is not
// fully supported, no messages will be emitted, but a log at FATAL
// severity will cause an abort.
//
// TODO(user): Come up with a better name for this method.
bool RawLoggingFullySupported();
}  // namespace internal

}  // namespace base_raw_logging

// TODO(user): Note: There are 10s of "RAW_LOG_xxxx" that prevent us
// from removing these.  Update these, then remove: DO NOT USE.
//
// DO NOT USE.
#define RAW_LOG_INFO(...) RAW_LOG(INFO, __VA_ARGS__)
#define RAW_LOG_WARNING(...) RAW_LOG(WARNING, __VA_ARGS__)
#define RAW_LOG_ERROR(...) RAW_LOG(ERROR, __VA_ARGS__)
#define RAW_LOG_FATAL(...) RAW_LOG(FATAL, __VA_ARGS__)

#endif  // S2_THIRD_PARTY_ABSL_BASE_RAW_LOGGING_H_

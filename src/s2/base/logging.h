// Copyright Google Inc. All Rights Reserved.
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

#ifndef S2_BASE_LOGGING_H_
#define S2_BASE_LOGGING_H_

// TODO(user): Get rid of `base/logging.h` includes and
// include the relevant absl file directly instead.
#include "absl/log/check.h"
#include "absl/log/log.h"
#include "s2/base/log_severity.h"

// The names CHECK, etc. are too common and may conflict with other
// packages.  We use S2_CHECK to make it easier to switch to
// something other than abseil-cpp for logging.
// TODO(user): Remove these or make absl::log optional.

#define S2_LOG LOG
#define S2_LOG_IF LOG_IF
#define S2_DLOG DLOG
#define S2_DLOG_IF DLOG_IF

#define S2_CHECK CHECK
#define S2_CHECK_EQ CHECK_EQ
#define S2_CHECK_NE CHECK_NE
#define S2_CHECK_LT CHECK_LT
#define S2_CHECK_LE CHECK_LE
#define S2_CHECK_GT CHECK_GT
#define S2_CHECK_GE CHECK_GE

#define S2_DCHECK DCHECK
#define S2_DCHECK_EQ DCHECK_EQ
#define S2_DCHECK_NE DCHECK_NE
#define S2_DCHECK_LT DCHECK_LT
#define S2_DCHECK_LE DCHECK_LE
#define S2_DCHECK_GT DCHECK_GT
#define S2_DCHECK_GE DCHECK_GE

// Logging stream that does nothing.
struct S2NullStream {
  template <typename T>
  S2NullStream& operator<<(const T& v) { return *this; }
};

// Abseil-cpp doesn't support VLOG yet.  Make VLOG a no-op.
#define S2_VLOG(verbose_level) S2NullStream()
#define S2_VLOG_IS_ON(verbose_level) (false)

#endif  // S2_BASE_LOGGING_H_

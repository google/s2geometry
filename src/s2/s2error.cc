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

#include "s2/s2error.h"

#include <string>

#include "absl/base/optimization.h"
#include "absl/status/status.h"
#include "absl/strings/str_cat.h"
#include "absl/strings/string_view.h"

using absl::string_view;

S2Error ToS2Error(const absl::Status& status) {
  S2Error error;

  if (status.ok()) {
    return error;
  }

  const string_view message = status.message();

  switch (status.code()) {
    case absl::StatusCode::kCancelled:
      error = S2Error::Cancelled(message);
      break;

    case absl::StatusCode::kInvalidArgument:
      error = S2Error::InvalidArgument(message);
      break;

    case absl::StatusCode::kDataLoss:
      error = S2Error::DataLoss(message);
      break;

    case absl::StatusCode::kResourceExhausted:
      error = S2Error::ResourceExhausted(message);
      break;

    case absl::StatusCode::kFailedPrecondition:
      error = S2Error::FailedPrecondition(message);
      break;

    case absl::StatusCode::kOutOfRange:
      error = S2Error::OutOfRange(message);
      break;

    case absl::StatusCode::kUnimplemented:
      error = S2Error::Unimplemented(message);
      break;

    case absl::StatusCode::kInternal:
      error = S2Error::Internal(message);
      break;

    default:
      error = S2Error::Unknown(message);
  }

  return error;
}

absl::Status ToStatus(const S2Error& error) {
  // Map user-defined error codes to `UnknownError`.  We do this outside of
  // the `switch` so we can handle all cases explicitly without `default`.
  if (S2Error::USER_DEFINED_START <= error.code() &&
      error.code() <= S2Error::USER_DEFINED_END) {
    return absl::UnknownError(
        absl::StrCat(error.message(), " (", error.code(), ")"));
  }

  switch (static_cast<S2Error::Code>(error.code())) {
    case S2Error::OK:
      return absl::OkStatus();

    case S2Error::CANCELLED:
      return absl::CancelledError(error.message());

    case S2Error::INVALID_ARGUMENT:
      return absl::InvalidArgumentError(error.message());

    case S2Error::DATA_LOSS:
      return absl::DataLossError(error.message());

    case S2Error::RESOURCE_EXHAUSTED:
      return absl::ResourceExhaustedError(error.message());

    case S2Error::FAILED_PRECONDITION:
      return absl::FailedPreconditionError(error.message());

    // Custom error space are deprecated and not open-sourced, so we map
    // these to invalid argument errors.
    case S2Error::NOT_CONTINUOUS:
    case S2Error::INVALID_VERTEX:
    case S2Error::INVALID_DIMENSION:
    case S2Error::SPLIT_INTERIOR:
    case S2Error::OVERLAPPING_GEOMETRY:
    case S2Error::NOT_UNIT_LENGTH:
    case S2Error::DUPLICATE_VERTICES:
    case S2Error::ANTIPODAL_VERTICES:
    case S2Error::LOOP_NOT_ENOUGH_VERTICES:
    case S2Error::LOOP_SELF_INTERSECTION:
    case S2Error::POLYGON_LOOPS_SHARE_EDGE:
    case S2Error::POLYGON_LOOPS_CROSS:
    case S2Error::POLYGON_EMPTY_LOOP:
    case S2Error::POLYGON_EXCESS_FULL_LOOP:
    case S2Error::POLYGON_INCONSISTENT_LOOP_ORIENTATIONS:
    case S2Error::POLYGON_INVALID_LOOP_DEPTH:
    case S2Error::POLYGON_INVALID_LOOP_NESTING:
    case S2Error::BUILDER_SNAP_RADIUS_TOO_SMALL:
    case S2Error::BUILDER_MISSING_EXPECTED_SIBLING_EDGES:
    case S2Error::BUILDER_UNEXPECTED_DEGENERATE_EDGE:
    case S2Error::BUILDER_EDGES_DO_NOT_FORM_LOOPS:
    case S2Error::BUILDER_EDGES_DO_NOT_FORM_POLYLINE:
    case S2Error::BUILDER_IS_FULL_PREDICATE_NOT_SPECIFIED:
      return absl::InvalidArgumentError(
          absl::StrCat(error.message(), " (", error.code(), ")"));

    case S2Error::OUT_OF_RANGE:
      return absl::OutOfRangeError(error.message());

    case S2Error::UNIMPLEMENTED:
      return absl::UnimplementedError(error.message());

    case S2Error::INTERNAL:
      return absl::InternalError(error.message());

    case S2Error::UNKNOWN:
      return absl::UnknownError(
          absl::StrCat(error.message(), " (", error.code(), ")"));

    case S2Error::USER_DEFINED_START:
    case S2Error::USER_DEFINED_END:
      // USER_DEFINED_START/END (or anything in between) can't happen here;
      // they're handled above.
      ABSL_UNREACHABLE();
  }
  ABSL_UNREACHABLE();
}

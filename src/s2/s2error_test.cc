// Copyright 2013 Google Inc. All Rights Reserved.
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

// Author: ericv@google.com (Eric Veach)

#include "s2/s2error.h"

#include <string>

#include <gmock/gmock.h>
#include <gtest/gtest.h>
#include "absl/status/status.h"
#include "absl/strings/str_format.h"

TEST(S2Error, Basic) {
  S2Error error;
  error = S2Error(S2Error::DUPLICATE_VERTICES,
                  "Vertex 23 is the same as vertex 47");
  // Prepend additional context to the message.
  error =
      S2Error(error.code(), absl::StrFormat("Loop %d: %s", 5, error.message()));
  EXPECT_EQ(error.code(), S2Error::DUPLICATE_VERTICES);
  EXPECT_EQ(error.message(), "Loop 5: Vertex 23 is the same as vertex 47");
}

TEST(S2Error, Constructor) {
  S2Error error(S2Error::RESOURCE_EXHAUSTED,
                "Memory limit exceeded (100 vs 50)");
  EXPECT_EQ(error.code(), S2Error::RESOURCE_EXHAUSTED);
  EXPECT_EQ(error.message(), "Memory limit exceeded (100 vs 50)");
}

TEST(S2Error, OkIsOk) {
  EXPECT_TRUE(S2Error::Ok().ok());
  EXPECT_EQ(S2Error::Ok().code(), S2Error::OK);
}

TEST(S2Error, UnknownCode) {
  EXPECT_EQ(S2Error::Unknown("").code(), S2Error::UNKNOWN);
}

TEST(S2Error, UnimplementedCode) {
  EXPECT_EQ(S2Error::Unimplemented("").code(), S2Error::UNIMPLEMENTED);
}

TEST(S2Error, OutOfRangeCode) {
  EXPECT_EQ(S2Error::OutOfRange("").code(), S2Error::OUT_OF_RANGE);
}

TEST(S2Error, InvalidArgumentCode) {
  EXPECT_EQ(S2Error::InvalidArgument("").code(), S2Error::INVALID_ARGUMENT);
}

TEST(S2Error, FailedPreconditionCode) {
  EXPECT_EQ(S2Error::FailedPrecondition("").code(),
            S2Error::FAILED_PRECONDITION);
}

TEST(S2Error, InternalCode) {
  EXPECT_EQ(S2Error::Internal("").code(), S2Error::INTERNAL);
}

TEST(S2Error, ResourceExhaustedCode) {
  EXPECT_EQ(S2Error::ResourceExhausted("").code(), S2Error::RESOURCE_EXHAUSTED);
}

TEST(S2Error, CancelledCode) {
  EXPECT_EQ(S2Error::Cancelled("").code(), S2Error::CANCELLED);
}

TEST(S2Error, ToS2Error) {
  const S2Error ok = ToS2Error(absl::Status());
  EXPECT_EQ(ok.code(), S2Error::OK);

  const S2Error cancelled = ToS2Error(absl::CancelledError("cancelled"));
  EXPECT_EQ(cancelled.code(), S2Error::CANCELLED);
  EXPECT_EQ(cancelled.message(), "cancelled");

  const S2Error invalid_argument =
      ToS2Error(absl::InvalidArgumentError("invalid_argument"));
  EXPECT_EQ(invalid_argument.code(), S2Error::INVALID_ARGUMENT);
  EXPECT_EQ(invalid_argument.message(), "invalid_argument");

  const S2Error data_loss = ToS2Error(absl::DataLossError("data_loss"));
  EXPECT_EQ(data_loss.code(), S2Error::DATA_LOSS);
  EXPECT_EQ(data_loss.message(), "data_loss");

  const S2Error resource_exhausted =
      ToS2Error(absl::ResourceExhaustedError("resource_exhausted"));
  EXPECT_EQ(resource_exhausted.code(), S2Error::RESOURCE_EXHAUSTED);
  EXPECT_EQ(resource_exhausted.message(), "resource_exhausted");

  const S2Error failed_precondition =
      ToS2Error(absl::FailedPreconditionError("failed_precondition"));
  EXPECT_EQ(failed_precondition.code(), S2Error::FAILED_PRECONDITION);
  EXPECT_EQ(failed_precondition.message(), "failed_precondition");

  const S2Error out_of_range = ToS2Error(absl::OutOfRangeError("out_of_range"));
  EXPECT_EQ(out_of_range.code(), S2Error::OUT_OF_RANGE);
  EXPECT_EQ(out_of_range.message(), "out_of_range");

  const S2Error unimplemented =
      ToS2Error(absl::UnimplementedError("unimplemented"));
  EXPECT_EQ(unimplemented.code(), S2Error::UNIMPLEMENTED);
  EXPECT_EQ(unimplemented.message(), "unimplemented");

  const S2Error internal = ToS2Error(absl::InternalError("internal"));
  EXPECT_EQ(internal.code(), S2Error::INTERNAL);
  EXPECT_EQ(internal.message(), "internal");

  const S2Error unknown = ToS2Error(absl::UnknownError("unknown"));
  EXPECT_EQ(unknown.code(), S2Error::UNKNOWN);
  EXPECT_THAT(unknown.message(), "unknown");

  // All absl::Status codes that don't have an exact mapping to S2Error codes
  // should be mapped to S2Error::UNKNOWN.
  for (const absl::Status& status :
       {absl::AbortedError("other"), absl::AlreadyExistsError("other"),
        absl::DeadlineExceededError("other"), absl::NotFoundError("other"),
        absl::PermissionDeniedError("other"),
        absl::UnauthenticatedError("other"), absl::UnavailableError("other")}) {
    const S2Error other = ToS2Error(status);
    EXPECT_EQ(other.code(), S2Error::UNKNOWN);
    EXPECT_EQ(other.message(), "other");
  }
}

TEST(S2Error, ToStatus) {
  S2Error error;

  const absl::Status ok = ToStatus(error);
  EXPECT_EQ(ok.code(), absl::StatusCode::kOk);

  error = S2Error::Cancelled("cancelled");
  const absl::Status cancelled = ToStatus(error);
  EXPECT_EQ(cancelled.code(), absl::StatusCode::kCancelled);
  EXPECT_EQ(cancelled.message(), "cancelled");

  error = S2Error::InvalidArgument("invalid_argument");
  const absl::Status invalid_argument = ToStatus(error);
  EXPECT_EQ(invalid_argument.code(), absl::StatusCode::kInvalidArgument);
  EXPECT_EQ(invalid_argument.message(), "invalid_argument");

  error = S2Error::DataLoss("data_loss");
  const absl::Status data_loss = ToStatus(error);
  EXPECT_EQ(data_loss.code(), absl::StatusCode::kDataLoss);
  EXPECT_EQ(data_loss.message(), "data_loss");

  error = S2Error::ResourceExhausted("resource_exhausted");
  const absl::Status resource_exhausted = ToStatus(error);
  EXPECT_EQ(resource_exhausted.code(), absl::StatusCode::kResourceExhausted);
  EXPECT_EQ(resource_exhausted.message(), "resource_exhausted");

  error = S2Error::FailedPrecondition("failed_precondition");
  const absl::Status failed_precondition = ToStatus(error);
  EXPECT_EQ(failed_precondition.code(), absl::StatusCode::kFailedPrecondition);
  EXPECT_EQ(failed_precondition.message(), "failed_precondition");

  for (const S2Error::Code code :
       {S2Error::NOT_UNIT_LENGTH,
        S2Error::DUPLICATE_VERTICES,
        S2Error::ANTIPODAL_VERTICES,
        S2Error::LOOP_NOT_ENOUGH_VERTICES,
        S2Error::LOOP_SELF_INTERSECTION,
        S2Error::POLYGON_LOOPS_SHARE_EDGE,
        S2Error::POLYGON_LOOPS_CROSS,
        S2Error::POLYGON_EMPTY_LOOP,
        S2Error::POLYGON_EXCESS_FULL_LOOP,
        S2Error::POLYGON_INCONSISTENT_LOOP_ORIENTATIONS,
        S2Error::POLYGON_INVALID_LOOP_DEPTH,
        S2Error::POLYGON_INVALID_LOOP_NESTING,
        S2Error::BUILDER_SNAP_RADIUS_TOO_SMALL,
        S2Error::BUILDER_MISSING_EXPECTED_SIBLING_EDGES,
        S2Error::BUILDER_UNEXPECTED_DEGENERATE_EDGE,
        S2Error::BUILDER_EDGES_DO_NOT_FORM_LOOPS,
        S2Error::BUILDER_EDGES_DO_NOT_FORM_POLYLINE,
        S2Error::BUILDER_IS_FULL_PREDICATE_NOT_SPECIFIED}) {
    error = S2Error(code, "other_invalid_argument");
    const absl::Status other = ToStatus(error);
    EXPECT_EQ(other.code(), absl::StatusCode::kInvalidArgument);
    EXPECT_THAT(other.message(),
                testing::HasSubstr("other_invalid_argument"));
  }

  error = S2Error::OutOfRange("out_of_range");
  const absl::Status out_of_range = ToStatus(error);
  EXPECT_EQ(out_of_range.code(), absl::StatusCode::kOutOfRange);
  EXPECT_EQ(out_of_range.message(), "out_of_range");

  error = S2Error::Unimplemented("unimplemented");
  const absl::Status unimplemented = ToStatus(error);
  EXPECT_EQ(unimplemented.code(), absl::StatusCode::kUnimplemented);
  EXPECT_EQ(unimplemented.message(), "unimplemented");

  error = S2Error::Internal("internal");
  const absl::Status internal = ToStatus(error);
  EXPECT_EQ(internal.code(), absl::StatusCode::kInternal);
  EXPECT_EQ(internal.message(), "internal");

  error = S2Error::Unknown("unknown");
  const absl::Status unknown = ToStatus(error);
  EXPECT_EQ(unknown.code(), absl::StatusCode::kUnknown);
  EXPECT_THAT(unknown.message(), testing::HasSubstr("unknown"));

  // All S2Error codes that don't have an exact mapping to absl::Status codes
  // should be mapped to absl::StatusCode::kUnknown.
  for (const S2Error::Code code :
       {S2Error::USER_DEFINED_START,
        S2Error::USER_DEFINED_END}) {
    error = S2Error(code, "other");
    const absl::Status other = ToStatus(error);
    EXPECT_EQ(other.code(), absl::StatusCode::kUnknown);
    EXPECT_THAT(other.message(), testing::HasSubstr("other"));
  }
}

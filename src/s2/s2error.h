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
//
// S2Error is a simple class consisting of an error code and a human-readable
// error message.

#ifndef S2_S2ERROR_H_
#define S2_S2ERROR_H_

#include <cstdarg>
#include <ostream>
#include <string>

#include "absl/base/attributes.h"
#include "absl/status/status.h"
#include "absl/strings/str_format.h"
#include "absl/strings/string_view.h"

#include "s2/base/port.h"

// This class is intended to be copied by value as desired.  It uses
// the default copy constructor and assignment operator.
class S2Error {
 public:
  enum Code {
    OK = 0,  // No error.

    ////////////////////////////////////////////////////////////////////
    // Generic errors, not specific to geometric objects:

    UNKNOWN = 1000,              // Unknown error.
    UNIMPLEMENTED = 1001,        // Operation is not implemented.
    OUT_OF_RANGE = 1002,         // Argument is out of range.
    INVALID_ARGUMENT = 1003,     // Invalid argument (other than a range error).
    FAILED_PRECONDITION = 1004,  // Object is not in the required state.
    INTERNAL = 1005,             // An internal invariant has failed.
    DATA_LOSS = 1006,            // Data loss or corruption.
    RESOURCE_EXHAUSTED = 1007,   // A resource has been exhausted.
    CANCELLED = 1008,            // Operation was cancelled.

    ////////////////////////////////////////////////////////////////////
    // Error codes in the following range can be defined by clients:

    USER_DEFINED_START = 1000000,
    USER_DEFINED_END = 9999999,

    ////////////////////////////////////////////////////////////////////
    // Errors that apply to more than one type of geometry:

    NOT_UNIT_LENGTH = 1,     // Vertex is not unit length.
    DUPLICATE_VERTICES = 2,  // There are two identical vertices.
    ANTIPODAL_VERTICES = 3,  // There are two antipodal vertices.
    NOT_CONTINUOUS = 4,      // Edges of a chain aren't continuous.
    INVALID_VERTEX = 5,      // Vertex has value that's inf or NaN.

    ////////////////////////////////////////////////////////////////////
    // S2Loop errors:

    LOOP_NOT_ENOUGH_VERTICES = 100,  // Loop with fewer than 3 vertices.
    LOOP_SELF_INTERSECTION = 101,    // Loop has a self-intersection.

    ////////////////////////////////////////////////////////////////////
    // S2Polygon/S2Shape errors:

    POLYGON_LOOPS_SHARE_EDGE = 200,  // Two polygon loops share an edge.
    POLYGON_LOOPS_CROSS = 201,       // Two polygon loops cross.
    POLYGON_EMPTY_LOOP = 202,        // Polygon has an empty loop.
    POLYGON_EXCESS_FULL_LOOP = 203,  // Non-full polygon has a full loop.

    // Inconsistent loop orientations were detected, indicating that the
    // interior is not on the left of all edges.
    POLYGON_INCONSISTENT_LOOP_ORIENTATIONS = 204,

    // Loop depths don't correspond to any valid nesting hierarchy.
    POLYGON_INVALID_LOOP_DEPTH = 205,

    // Actual polygon nesting does not correspond to the nesting hierarchy
    // encoded by the loop depths.
    POLYGON_INVALID_LOOP_NESTING = 206,

    INVALID_DIMENSION = 207,     // Shape dimension isn't valid.
    SPLIT_INTERIOR = 208,        // Interior split by holes.
    OVERLAPPING_GEOMETRY = 209,  // Geometry overlaps where it shouldn't

    ////////////////////////////////////////////////////////////////////
    // S2Builder errors:

    // The S2Builder snap function moved a vertex by more than the specified
    // snap radius.
    BUILDER_SNAP_RADIUS_TOO_SMALL = 300,

    // S2Builder expected all edges to have siblings (as specified by
    // S2Builder::GraphOptions::SiblingPairs::REQUIRE), but some were missing.
    BUILDER_MISSING_EXPECTED_SIBLING_EDGES = 301,

    // S2Builder found an unexpected degenerate edge.  For example,
    // Graph::GetLeftTurnMap() does not support degenerate edges.
    BUILDER_UNEXPECTED_DEGENERATE_EDGE = 302,

    // S2Builder found a vertex with (indegree != outdegree), which means
    // that the given edges cannot be assembled into loops.
    BUILDER_EDGES_DO_NOT_FORM_LOOPS = 303,

    // The edges provided to S2Builder cannot be assembled into a polyline.
    BUILDER_EDGES_DO_NOT_FORM_POLYLINE = 304,

    // There was an attempt to assemble a polygon from degenerate geometry
    // without having specified a predicate to decide whether the output is
    // the empty polygon (containing no points) or the full polygon
    // (containing all points).
    BUILDER_IS_FULL_PREDICATE_NOT_SPECIFIED = 305,
  };

  // Creates error with OK status and empty message.
  S2Error() = default;

  // Create an error with the specified code and message.
  //
  // Note that you can prepend text to an existing error by reinitializeing
  // like this:
  //
  //   error = S2Error(error.code(),
  //                   absl::StrFormat("Loop %d: %s", j, error.message()));
  S2Error(Code code, absl::string_view message) : code_(code), text_(message) {}

  ABSL_MUST_USE_RESULT bool ok() const { return code_ == OK; }
  Code code() const { return code_; }
  absl::string_view message() const ABSL_ATTRIBUTE_LIFETIME_BOUND {
    return text_;
  }
  [[deprecated("Use message() instead.")]]
  std::string text() const {
    return text_;
  }

  // Returns a value with an `OK` code.  Prefer this to `S2Error()`.
  static S2Error Ok() {
    // This will be changed to `return absl::OkStatus();` when we switch over.
    return S2Error();
  }

  static S2Error Unknown(absl::string_view message) {
    // This will be changed to `absl::UnknownError()` when we switch over.
    return S2Error(UNKNOWN, message);
  }

  static S2Error Unimplemented(absl::string_view message) {
    return S2Error(UNIMPLEMENTED, message);
  }

  static S2Error OutOfRange(absl::string_view message) {
    return S2Error(OUT_OF_RANGE, message);
  }

  static S2Error InvalidArgument(absl::string_view message) {
    return S2Error(INVALID_ARGUMENT, message);
  }

  static S2Error FailedPrecondition(absl::string_view message) {
    return S2Error(FAILED_PRECONDITION, message);
  }

  static S2Error Internal(absl::string_view message) {
    return S2Error(INTERNAL, message);
  }

  static S2Error DataLoss(absl::string_view message) {
    return S2Error(DATA_LOSS, message);
  }

  static S2Error ResourceExhausted(absl::string_view message) {
    return S2Error(RESOURCE_EXHAUSTED, message);
  }

  static S2Error Cancelled(absl::string_view message) {
    return S2Error(CANCELLED, message);
  }

 private:
  Code code_ = OK;
  std::string text_;
};

// Converts an absl::Status object into an S2Error.
//
// absl::Status codes with no exact equivalent in S2Error::Code are converted to
// S2Error::UNKNOWN.
ABSL_MUST_USE_RESULT S2Error ToS2Error(const absl::Status& status);

// Converts an S2Error into an absl::Status.
//
// Custom error space are deprecated,
// so S2Error codes with no exact equivalent in absl::StatusCode are converted
// as follows:
//
//   - errors related to S2Loop, S2Polygon, and S2Builder, as well as those that
//     apply to more than one type of geometry are mapped to
//     absl::StatusCode::kInvalidArgument.
//
//   - client-defined errors (those in the [USER_DEFINED_START,USER_DEFINED_END]
//     range), UNKNOWN errors, and all other errors not covered above are mapped
//     to absl::StatusCode::kUnknown
//
// The above mapping significantly reduces the set of expressible errors.
// Clients reacting to specific S2 errors (e.g. polygon validation errors)
// should continue to use S2Error for its richer error API.
absl::Status ToStatus(const S2Error& error);


//////////////////   Implementation details follow   ////////////////////


inline std::ostream& operator<<(std::ostream& os, const S2Error& error) {
  return os << error.message();
}

#endif  // S2_S2ERROR_H_

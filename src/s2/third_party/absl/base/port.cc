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

// Platform-specific code that doesn't need to go in a header file.

// a given platform, so we include the body of the .cc file
// for the appropriate platform, using macros.

#include "s2/third_party/absl/base/port.h"

// Preprocessor token used to ensure that the header file mentioned below is not
// included by clients
#define THIRD_PARTY_ABSL_BASE_PORT_CC_

// Add other platforms here

#ifdef _MSC_VER
#include "s2/third_party/absl/base/internal/port_msvc_cc.inc"
#endif

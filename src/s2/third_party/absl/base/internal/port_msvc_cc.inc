#ifndef THIRD_PARTY_ABSL_BASE_INTERNAL_PORT_MSVC_CC_INC_
#define THIRD_PARTY_ABSL_BASE_INTERNAL_PORT_MSVC_CC_INC_

// MSVC Platform-specific code that doesn't need to be in port.h

// This header may be included only from base/port.cc
#ifndef THIRD_PARTY_ABSL_BASE_PORT_CC_
#error third_party/absl/base/internal/port_msvc_cc.inc \
may not be included directly
#endif  // THIRD_PARTY_ABSL_BASE_PORT_CC_

#include <cstdio>
#include <cstdarg>

#include "third_party/absl/base/port.h"

// Microsoft's _vsnprintf() and _snprintf() do not nul-terminate on buffer
// overflow, and do not return the correct length when buffer overflow occurs.
// These wrappers fix both issues.

int base_port_MSVC_vsnprintf(char *str, size_t size,
                             const char *format, va_list ap) {
  int count = _vsnprintf(str, size, format, ap);
  if (count < 0) {
    count = _vscprintf(format, ap);  // Yields character count.
  }
  if (size > 0 && count >= size) {
    str[size-1] = '\0';
  }
  return count;
}

int base_port_MSVC_snprintf(char *str, size_t size, const char *format, ...) {
  va_list ap;
  va_start(ap, format);
  int count = base_port_MSVC_vsnprintf(str, size, format, ap);
  va_end(ap);
  return count;
}

#endif  // THIRD_PARTY_ABSL_BASE_INTERNAL_PORT_MSVC_CC_INC_

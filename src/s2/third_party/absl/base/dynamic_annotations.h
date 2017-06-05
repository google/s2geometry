// This file includes the corresponding file from third_party.
// In debug mode (NDEBUG not defined) the dynamic annotations are enabled;
// in opt mode they are disabled.

#ifndef S2_THIRD_PARTY_ABSL_BASE_DYNAMIC_ANNOTATIONS_H_
#define S2_THIRD_PARTY_ABSL_BASE_DYNAMIC_ANNOTATIONS_H_

// Compiler-based ThreadSanitizer may be used with -c opt,
// but still requires dynamic annotations.
#ifndef DYNAMIC_ANNOTATIONS_ENABLED
#define BASE_DYNAMIC_ANNOTATIONS_UNDEF_WHEN_DONE
#if !defined(__APPLE__) && !defined(__asmjs__) &&              \
    (defined(MEMORY_SANITIZER) || defined(THREAD_SANITIZER) || \
     !defined(NDEBUG))
#define DYNAMIC_ANNOTATIONS_ENABLED 1
#else
#define DYNAMIC_ANNOTATIONS_ENABLED 0
#endif
#endif

#if !defined(__native_client__)
# include "s2/third_party/dynamic_annotations/dynamic_annotations.h"
#else
# include "s2/nacl/dynamic_annotations.h"
// Stub out the macros missing from the NaCl version.
# ifndef ANNOTATE_CONTIGUOUS_CONTAINER
#  define ANNOTATE_CONTIGUOUS_CONTAINER(beg, end, old_mid, new_mid)
# endif
# ifndef ANNOTATE_RWLOCK_CREATE_STATIC
#  define ANNOTATE_RWLOCK_CREATE_STATIC(lock)
# endif
# ifndef ADDRESS_SANITIZER_REDZONE
#  define ADDRESS_SANITIZER_REDZONE(name)
# endif
# ifndef ANNOTATE_MEMORY_IS_UNINITIALIZED
#  define ANNOTATE_MEMORY_IS_UNINITIALIZED(address, size)
# endif
#endif

#ifdef BASE_DYNAMIC_ANNOTATIONS_UNDEF_WHEN_DONE
# undef DYNAMIC_ANNOTATIONS_ENABLED
# undef BASE_DYNAMIC_ANNOTATIONS_UNDEF_WHEN_DONE
#endif

#endif  // S2_THIRD_PARTY_ABSL_BASE_DYNAMIC_ANNOTATIONS_H_

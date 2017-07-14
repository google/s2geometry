#ifdef __cplusplus
# error "This file should be built as pure C to avoid name mangling"
#endif

#include <cstdlib>
#include <cstring>

#include "third_party/absl/base/dynamic_annotations.h"

#ifndef __has_feature
#define __has_feature(x) 0
#endif

/* Compiler-based ThreadSanitizer defines
   DYNAMIC_ANNOTATIONS_EXTERNAL_IMPL = 1
   and provides its own definitions of the functions. */

#ifndef DYNAMIC_ANNOTATIONS_EXTERNAL_IMPL
# define DYNAMIC_ANNOTATIONS_EXTERNAL_IMPL 0
#endif

/* Each function is empty and called (via a macro) only in debug mode.
   The arguments are captured by dynamic tools at runtime. */

#if DYNAMIC_ANNOTATIONS_EXTERNAL_IMPL == 0 && !defined(__native_client__)

#if __has_feature(memory_sanitizer)
#include <sanitizer/msan_interface.h>
#endif

void AnnotateRWLockCreate(const char *file, int line,
                          const volatile void *lock){}
void AnnotateRWLockDestroy(const char *file, int line,
                           const volatile void *lock){}
void AnnotateRWLockAcquired(const char *file, int line,
                            const volatile void *lock, long is_w){}
void AnnotateRWLockReleased(const char *file, int line,
                            const volatile void *lock, long is_w){}
void AnnotateBenignRace(const char *file, int line,
                        const volatile void *address,
                        const char *description){}
void AnnotateBenignRaceSized(const char *file, int line,
                             const volatile void *address,
                             size_t size,
                             const char *description) {}
void AnnotateThreadName(const char *file, int line,
                        const char *name){}
void AnnotateIgnoreReadsBegin(const char *file, int line){}
void AnnotateIgnoreReadsEnd(const char *file, int line){}
void AnnotateIgnoreWritesBegin(const char *file, int line){}
void AnnotateIgnoreWritesEnd(const char *file, int line){}
void AnnotateEnableRaceDetection(const char *file, int line, int enable){}
void AnnotateMemoryIsInitialized(const char *file, int line,
                                 const volatile void *mem, size_t size) {
#if __has_feature(memory_sanitizer)
  __msan_unpoison(mem, size);
#endif
}

void AnnotateMemoryIsUninitialized(const char *file, int line,
                                   const volatile void *mem, size_t size) {
#if __has_feature(memory_sanitizer)
  __msan_allocated_memory(mem, size);
#endif
}

static int GetRunningOnValgrind(void) {
#ifdef RUNNING_ON_VALGRIND
  if (RUNNING_ON_VALGRIND) return 1;
#endif
  char *running_on_valgrind_str = getenv("RUNNING_ON_VALGRIND");
  if (running_on_valgrind_str) {
    return strcmp(running_on_valgrind_str, "0") != 0;
  }
  return 0;
}

/* See the comments in dynamic_annotations.h */
int RunningOnValgrind(void) {
  static volatile int running_on_valgrind = -1;
  int local_running_on_valgrind = running_on_valgrind;
  /* C doesn't have thread-safe initialization of statics, and we
     don't want to depend on pthread_once here, so hack it. */
  ANNOTATE_BENIGN_RACE(&running_on_valgrind, "safe hack");
  if (local_running_on_valgrind == -1)
    running_on_valgrind = local_running_on_valgrind = GetRunningOnValgrind();
  return local_running_on_valgrind;
}

/* See the comments in dynamic_annotations.h */
double ValgrindSlowdown(void) {
  /* Same initialization hack as in RunningOnValgrind(). */
  static volatile double slowdown = 0.0;
  double local_slowdown = slowdown;
  ANNOTATE_BENIGN_RACE(&slowdown, "safe hack");
  if (RunningOnValgrind() == 0) {
    return 1.0;
  }
  if (local_slowdown == 0.0) {
    char *env = getenv("VALGRIND_SLOWDOWN");
    slowdown = local_slowdown = env ? atof(env) : 50.0;
  }
  return local_slowdown;
}

#endif  /* DYNAMIC_ANNOTATIONS_EXTERNAL_IMPL == 0 */

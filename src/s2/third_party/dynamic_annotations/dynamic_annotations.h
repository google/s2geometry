/* Copyright (c) 2008-2009, Google Inc.
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are
 * met:
 *
 *     * Redistributions of source code must retain the above copyright
 * notice, this list of conditions and the following disclaimer.
 *     * Redistributions in binary form must reproduce the above
 * copyright notice, this list of conditions and the following disclaimer
 * in the documentation and/or other materials provided with the
 * distribution.
 *     * Neither the name of Google Inc. nor the names of its
 * contributors may be used to endorse or promote products derived from
 * this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 * A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
 * OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 * LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 * DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
 * THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 * ---
 * Author: Kostya Serebryany
 */

/* This file defines dynamic annotations for use with dynamic analysis
   tool such as valgrind, PIN, etc.

   Dynamic annotation is a source code annotation that affects
   the generated code (that is, the annotation is not a comment).
   Each such annotation is attached to a particular
   instruction and/or to a particular object (address) in the program.

   The annotations that should be used by users are macros in all upper-case
   (e.g., ANNOTATE_THREAD_NAME).

   Actual implementation of these macros may differ depending on the
   dynamic analysis tool being used.

   See http://code.google.com/p/data-race-test/  for more information.

   This file supports the following dynamic analysis tools:
   - None (DYNAMIC_ANNOTATIONS_ENABLED is not defined or zero).
      Macros are defined empty.
   - ThreadSanitizer and AddressSanitizer (DYNAMIC_ANNOTATIONS_ENABLED is 1).
      Macros are defined as calls to non-inlinable empty functions
      that are intercepted by the tool. */

#ifndef __DYNAMIC_ANNOTATIONS_H__
#define __DYNAMIC_ANNOTATIONS_H__

#ifndef DYNAMIC_ANNOTATIONS_ENABLED
# define DYNAMIC_ANNOTATIONS_ENABLED 0
#endif

#if DYNAMIC_ANNOTATIONS_ENABLED != 0

  /* -------------------------------------------------------------
     Annotations that suppress errors.  It is usually better to express the
     program's synchronization using the other annotations, but these can
     be used when all else fails. */

  /* Report that we may have a benign race at "pointer", with size
     "sizeof(*(pointer))". "pointer" must be a non-void* pointer.  Insert at the
     point where "pointer" has been allocated, preferably close to the point
     where the race happens.  See also ANNOTATE_BENIGN_RACE_STATIC. */
  #define ANNOTATE_BENIGN_RACE(pointer, description) \
    AnnotateBenignRaceSized(__FILE__, __LINE__, pointer, \
                            sizeof(*(pointer)), description)

  /* Same as ANNOTATE_BENIGN_RACE(address, description), but applies to
     the memory range [address, address+size). */
  #define ANNOTATE_BENIGN_RACE_SIZED(address, size, description) \
    AnnotateBenignRaceSized(__FILE__, __LINE__, address, size, description)

  /* Request the analysis tool to ignore all reads in the current thread
     until ANNOTATE_IGNORE_READS_END is called.
     Useful to ignore intentional racey reads, while still checking
     other reads and all writes.
     See also ANNOTATE_UNPROTECTED_READ. */
  #define ANNOTATE_IGNORE_READS_BEGIN() \
    AnnotateIgnoreReadsBegin(__FILE__, __LINE__)

  /* Stop ignoring reads. */
  #define ANNOTATE_IGNORE_READS_END() \
    AnnotateIgnoreReadsEnd(__FILE__, __LINE__)

  /* Similar to ANNOTATE_IGNORE_READS_BEGIN, but ignore writes. */
  #define ANNOTATE_IGNORE_WRITES_BEGIN() \
    AnnotateIgnoreWritesBegin(__FILE__, __LINE__)

  /* Stop ignoring writes. */
  #define ANNOTATE_IGNORE_WRITES_END() \
    AnnotateIgnoreWritesEnd(__FILE__, __LINE__)

  /* Start ignoring all memory accesses (reads and writes). */
  #define ANNOTATE_IGNORE_READS_AND_WRITES_BEGIN() \
    do {\
      ANNOTATE_IGNORE_READS_BEGIN();\
      ANNOTATE_IGNORE_WRITES_BEGIN();\
    }while(0)\

  /* Stop ignoring all memory accesses. */
  #define ANNOTATE_IGNORE_READS_AND_WRITES_END() \
    do {\
      ANNOTATE_IGNORE_WRITES_END();\
      ANNOTATE_IGNORE_READS_END();\
    }while(0)\

  /* Enable (enable!=0) or disable (enable==0) race detection for all threads.
     This annotation could be useful if you want to skip expensive race analysis
     during some period of program execution, e.g. during initialization. */
  #define ANNOTATE_ENABLE_RACE_DETECTION(enable) \
    AnnotateEnableRaceDetection(__FILE__, __LINE__, enable)

  /* -------------------------------------------------------------
     Annotations useful for debugging. */

  /* Request to trace every access to "address". */
  #define ANNOTATE_TRACE_MEMORY(address) \
    AnnotateTraceMemory(__FILE__, __LINE__, address)

  /* Report the current thread name to a race detector. */
  #define ANNOTATE_THREAD_NAME(name) \
    AnnotateThreadName(__FILE__, __LINE__, name)

  /* -------------------------------------------------------------
     Annotations useful when implementing locks.  They are not
     normally needed by modules that merely use locks.
     The "lock" argument is a pointer to the lock object. */

  /* Report that a lock has been created at address "lock". */
  #define ANNOTATE_RWLOCK_CREATE(lock) \
    AnnotateRWLockCreate(__FILE__, __LINE__, lock)

  /* Report that a linker initialized lock has been created at address "lock".
   */
#ifdef THREAD_SANITIZER
  #define ANNOTATE_RWLOCK_CREATE_STATIC(lock) \
    AnnotateRWLockCreateStatic(__FILE__, __LINE__, lock)
#else
  #define ANNOTATE_RWLOCK_CREATE_STATIC(lock) ANNOTATE_RWLOCK_CREATE(lock)
#endif

  /* Report that the lock at address "lock" is about to be destroyed. */
  #define ANNOTATE_RWLOCK_DESTROY(lock) \
    AnnotateRWLockDestroy(__FILE__, __LINE__, lock)

  /* Report that the lock at address "lock" has been acquired.
     is_w=1 for writer lock, is_w=0 for reader lock. */
  #define ANNOTATE_RWLOCK_ACQUIRED(lock, is_w) \
    AnnotateRWLockAcquired(__FILE__, __LINE__, lock, is_w)

  /* Report that the lock at address "lock" is about to be released. */
  #define ANNOTATE_RWLOCK_RELEASED(lock, is_w) \
    AnnotateRWLockReleased(__FILE__, __LINE__, lock, is_w)

  /* -------------------------------------------------------------
     Annotations useful when implementing barriers.  They are not
     normally needed by modules that merely use barriers.
     The "barrier" argument is a pointer to the barrier object. */

  /* Report that the "barrier" has been initialized with initial "count".
   If 'reinitialization_allowed' is true, initialization is allowed to happen
   multiple times w/o calling barrier_destroy() */
  #define ANNOTATE_BARRIER_INIT(barrier, count, reinitialization_allowed) \
    AnnotateBarrierInit(__FILE__, __LINE__, barrier, count, \
                        reinitialization_allowed)

  /* Report that we are about to enter barrier_wait("barrier"). */
  #define ANNOTATE_BARRIER_WAIT_BEFORE(barrier) \
    AnnotateBarrierWaitBefore(__FILE__, __LINE__, barrier)

  /* Report that we just exited barrier_wait("barrier"). */
  #define ANNOTATE_BARRIER_WAIT_AFTER(barrier) \
    AnnotateBarrierWaitAfter(__FILE__, __LINE__, barrier)

  /* Report that the "barrier" has been destroyed. */
  #define ANNOTATE_BARRIER_DESTROY(barrier) \
    AnnotateBarrierDestroy(__FILE__, __LINE__, barrier)

  /* -------------------------------------------------------------
     Annotations useful for testing race detectors. */

  /* Report that we expect a race on the variable at "address".
     Use only in unit tests for a race detector. */
  #define ANNOTATE_EXPECT_RACE(address, description) \
    AnnotateExpectRace(__FILE__, __LINE__, address, description)

  /* A no-op. Insert where you like to test the interceptors. */
  #define ANNOTATE_NO_OP(arg) \
    AnnotateNoOp(__FILE__, __LINE__, arg)

  #define ANNOTATE_MEMORY_IS_INITIALIZED(address, size) \
    AnnotateMemoryIsInitialized(__FILE__, __LINE__, address, size)

  #define ANNOTATE_MEMORY_IS_UNINITIALIZED(address, size) \
    AnnotateMemoryIsUninitialized(__FILE__, __LINE__, address, size)

#else  /* DYNAMIC_ANNOTATIONS_ENABLED == 0 */

  #define ANNOTATE_RWLOCK_CREATE(lock) /* empty */
  #define ANNOTATE_RWLOCK_CREATE_STATIC(lock) /* empty */
  #define ANNOTATE_RWLOCK_DESTROY(lock) /* empty */
  #define ANNOTATE_RWLOCK_ACQUIRED(lock, is_w) /* empty */
  #define ANNOTATE_RWLOCK_RELEASED(lock, is_w) /* empty */
  #define ANNOTATE_BARRIER_INIT(barrier, count, reinitialization_allowed) /* */
  #define ANNOTATE_BARRIER_WAIT_BEFORE(barrier) /* empty */
  #define ANNOTATE_BARRIER_WAIT_AFTER(barrier) /* empty */
  #define ANNOTATE_BARRIER_DESTROY(barrier) /* empty */
  #define ANNOTATE_EXPECT_RACE(address, description) /* empty */
  #define ANNOTATE_BENIGN_RACE(address, description) /* empty */
  #define ANNOTATE_BENIGN_RACE_SIZED(address, size, description) /* empty */
  #define ANNOTATE_TRACE_MEMORY(arg) /* empty */
  #define ANNOTATE_THREAD_NAME(name) /* empty */
  #define ANNOTATE_IGNORE_READS_BEGIN() /* empty */
  #define ANNOTATE_IGNORE_READS_END() /* empty */
  #define ANNOTATE_IGNORE_WRITES_BEGIN() /* empty */
  #define ANNOTATE_IGNORE_WRITES_END() /* empty */
  #define ANNOTATE_IGNORE_READS_AND_WRITES_BEGIN() /* empty */
  #define ANNOTATE_IGNORE_READS_AND_WRITES_END() /* empty */
  #define ANNOTATE_ENABLE_RACE_DETECTION(enable) /* empty */
  #define ANNOTATE_NO_OP(arg) /* empty */
  #define ANNOTATE_MEMORY_IS_INITIALIZED(address, size) /* empty */
  #define ANNOTATE_MEMORY_IS_UNINITIALIZED(address, size) /* empty */

#endif  /* DYNAMIC_ANNOTATIONS_ENABLED */

/* Clang provides limited support for static thread-safety analysis
   through a feature called Annotalysis. We configure macro-definitions
   according to whether (and which) Annotalysis support is available. */

/* TODO(user) -- Replace __CLANG_SUPPORT_DYN_ANNOTATION__ with the
   appropriate feature ID. */
#if defined(__clang__) && (!defined(SWIG)) \
    && defined(__CLANG_SUPPORT_DYN_ANNOTATION__)
#define ANNOTALYSIS_SUPPORT
#endif

/* TODO(user) -- The exclusive lock here ignores writes as well, but
   allows INGORE_READS_AND_WRITES to work properly. */
#ifdef ANNOTALYSIS_SUPPORT
  #define ANNOTALYSIS_IGNORE_READS_BEGIN \
      __attribute__((exclusive_lock_function("*")))
  #define ANNOTALYSIS_IGNORE_READS_END \
      __attribute__((unlock_function("*")))
#else
  #define ANNOTALYSIS_IGNORE_READS_BEGIN
  #define ANNOTALYSIS_IGNORE_READS_END
#endif  /* ANNOTALYSIS_SUPPORT */

#if DYNAMIC_ANNOTATIONS_ENABLED == 0 && defined(ANNOTALYSIS_SUPPORT)
/* Turn on certain macros for static analysis, even if dynamic annotations are
   not enabled. */
  #define ANNOTALYSIS_STATIC_INLINE static inline
  #define ANNOTALYSIS_SEMICOLON_OR_EMPTY_BODY { (void)file; (void)line; }
#else
  #define ANNOTALYSIS_STATIC_INLINE
  #define ANNOTALYSIS_SEMICOLON_OR_EMPTY_BODY ;
#endif  /* DYNAMIC_ANNOTATIONS_ENABLED == 0 && defined(ANNOTALYSIS_SUPPORT) */

/* Use the macros above rather than using these functions directly. */
#include <cstddef>
#ifdef __cplusplus
extern "C" {
#endif
void AnnotateRWLockCreate(const char *file, int line,
                          const volatile void *lock);
void AnnotateRWLockCreateStatic(const char *file, int line,
                          const volatile void *lock);
void AnnotateRWLockDestroy(const char *file, int line,
                           const volatile void *lock);
void AnnotateRWLockAcquired(const char *file, int line,
                            const volatile void *lock, long is_w);
void AnnotateRWLockReleased(const char *file, int line,
                            const volatile void *lock, long is_w);
void AnnotateBarrierInit(const char *file, int line,
                         const volatile void *barrier, long count,
                         long reinitialization_allowed);
void AnnotateBarrierWaitBefore(const char *file, int line,
                               const volatile void *barrier);
void AnnotateBarrierWaitAfter(const char *file, int line,
                              const volatile void *barrier);
void AnnotateBarrierDestroy(const char *file, int line,
                            const volatile void *barrier);
void AnnotateExpectRace(const char *file, int line,
                        const volatile void *address,
                        const char *description);
void AnnotateBenignRace(const char *file, int line,
                        const volatile void *address,
                        const char *description);
void AnnotateBenignRaceSized(const char *file, int line,
                        const volatile void *address,
                        size_t size,
                        const char *description);
void AnnotateTraceMemory(const char *file, int line,
                         const volatile void *arg);
void AnnotateThreadName(const char *file, int line,
                        const char *name);
ANNOTALYSIS_STATIC_INLINE
void AnnotateIgnoreReadsBegin(const char *file, int line)
    ANNOTALYSIS_IGNORE_READS_BEGIN ANNOTALYSIS_SEMICOLON_OR_EMPTY_BODY
ANNOTALYSIS_STATIC_INLINE
void AnnotateIgnoreReadsEnd(const char *file, int line)
    ANNOTALYSIS_IGNORE_READS_END ANNOTALYSIS_SEMICOLON_OR_EMPTY_BODY
ANNOTALYSIS_STATIC_INLINE
void AnnotateIgnoreWritesBegin(const char *file, int line)
    ANNOTALYSIS_SEMICOLON_OR_EMPTY_BODY
ANNOTALYSIS_STATIC_INLINE
void AnnotateIgnoreWritesEnd(const char *file, int line)
    ANNOTALYSIS_SEMICOLON_OR_EMPTY_BODY
void AnnotateEnableRaceDetection(const char *file, int line, int enable);
void AnnotateNoOp(const char *file, int line,
                  const volatile void *arg);
void AnnotateMemoryIsInitialized(const char *file, int line,
                                 const volatile void *mem, size_t size);
void AnnotateMemoryIsUninitialized(const char *file, int line,
                                   const volatile void *mem, size_t size);

/* Return non-zero value if running under valgrind.

  If "valgrind.h" is included into dynamic_annotations.c,
  the regular valgrind mechanism will be used.
  See http://valgrind.org/docs/manual/manual-core-adv.html about
  RUNNING_ON_VALGRIND and other valgrind "client requests".
  The file "valgrind.h" may be obtained by doing
     svn co svn://svn.valgrind.org/valgrind/trunk/include

  If for some reason you can't use "valgrind.h" or want to fake valgrind,
  there are two ways to make this function return non-zero:
    - Use environment variable: export RUNNING_ON_VALGRIND=1
    - Make your tool intercept the function RunningOnValgrind() and
      change its return value.
 */
int RunningOnValgrind(void);

/* ValgrindSlowdown returns:
    * 1.0, if (RunningOnValgrind() == 0)
    * 50.0, if (RunningOnValgrind() != 0 && getenv("VALGRIND_SLOWDOWN") == nullptr)
    * atof(getenv("VALGRIND_SLOWDOWN")) otherwise
   This function can be used to scale timeout values:
   EXAMPLE:
   for (;;) {
     DoExpensiveBackgroundTask();
     SleepForSeconds(5 * ValgrindSlowdown());
   }
 */
double ValgrindSlowdown(void);

#ifdef __cplusplus
}
#endif

#if DYNAMIC_ANNOTATIONS_ENABLED != 0 && defined(__cplusplus)

  /* ANNOTATE_UNPROTECTED_READ is the preferred way to annotate racey reads.

     Instead of doing
        ANNOTATE_IGNORE_READS_BEGIN();
        ... = x;
        ANNOTATE_IGNORE_READS_END();
     one can use
        ... = ANNOTATE_UNPROTECTED_READ(x); */
  template <class T>
  inline T ANNOTATE_UNPROTECTED_READ(const volatile T &x) {
    ANNOTATE_IGNORE_READS_BEGIN();
    T res = x;
    ANNOTATE_IGNORE_READS_END();
    return res;
  }
  /* Apply ANNOTATE_BENIGN_RACE_SIZED to a static variable. */
  #define ANNOTATE_BENIGN_RACE_STATIC(static_var, description)        \
    namespace {                                                       \
      class static_var ## _annotator {                                \
       public:                                                        \
        static_var ## _annotator() {                                  \
          ANNOTATE_BENIGN_RACE_SIZED(&static_var,                     \
                                      sizeof(static_var),             \
            # static_var ": " description);                           \
        }                                                             \
      };                                                              \
      static static_var ## _annotator the ## static_var ## _annotator;\
    }
#else /* DYNAMIC_ANNOTATIONS_ENABLED == 0 */

  #define ANNOTATE_UNPROTECTED_READ(x) (x)
  #define ANNOTATE_BENIGN_RACE_STATIC(static_var, description)  /* empty */

#endif /* DYNAMIC_ANNOTATIONS_ENABLED */

#ifdef ADDRESS_SANITIZER
/* Describe the current state of a contiguous container such as e.g.
 * std::vector or std::string. For more details see
 * sanitizer/common_interface_defs.h, which is provided by the compiler. */
#include <sanitizer/common_interface_defs.h>
#define ANNOTATE_CONTIGUOUS_CONTAINER(beg, end, old_mid, new_mid) \
  __sanitizer_annotate_contiguous_container(beg, end, old_mid, new_mid)
#define ADDRESS_SANITIZER_REDZONE(name)         \
  struct { char x[8] __attribute__ ((aligned (8))); } name
#else
#define ANNOTATE_CONTIGUOUS_CONTAINER(beg, end, old_mid, new_mid)
#define ADDRESS_SANITIZER_REDZONE(name)
#endif  // ADDRESS_SANITIZER

#if DYNAMIC_ANNOTATIONS_ENABLED==0 && defined(ANNOTALYSIS_SUPPORT)

/* Turn on macros that the static analyzer understands.  These should be on
 * even if dynamic annotations are off. */

  #undef ANNOTATE_IGNORE_READS_BEGIN
  #define ANNOTATE_IGNORE_READS_BEGIN() \
    AnnotateIgnoreReadsBegin(__FILE__, __LINE__)

  #undef ANNOTATE_IGNORE_READS_END
  #define ANNOTATE_IGNORE_READS_END() \
    AnnotateIgnoreReadsEnd(__FILE__, __LINE__)

  #undef ANNOTATE_IGNORE_READS_AND_WRITES_BEGIN
  #define ANNOTATE_IGNORE_READS_AND_WRITES_BEGIN() \
    do {                                           \
      ANNOTATE_IGNORE_READS_BEGIN();               \
      ANNOTATE_IGNORE_WRITES_BEGIN();              \
    } while (0)                                    \

  #undef ANNOTATE_IGNORE_READS_AND_WRITES_END
  #define ANNOTATE_IGNORE_READS_AND_WRITES_END()   \
    do {                                           \
      ANNOTATE_IGNORE_WRITES_END();                \
      ANNOTATE_IGNORE_READS_END();                 \
    } while (0)                                    \

  #if defined(__cplusplus)
  #undef ANNOTATE_UNPROTECTED_READ
  template <class T>
  inline T ANNOTATE_UNPROTECTED_READ(const volatile T &x) {
    ANNOTATE_IGNORE_READS_BEGIN();
    T res = x;
    ANNOTATE_IGNORE_READS_END();
    return res;
  }
  #endif

#endif  /* DYNAMIC_ANNOTATIONS_ENABLED==0 && defined(ANNOTALYSIS_SUPPORT) */


/* Undefine the macros intended only in this file. */
#undef ANNOTALYSIS_STATIC_INLINE
#undef ANNOTALYSIS_SEMICOLON_OR_EMPTY_BODY
#undef ANNOTALYSIS_SUPPORT

#endif  /* __DYNAMIC_ANNOTATIONS_H__ */

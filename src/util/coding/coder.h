// Copyright 2000 Google Inc. All Rights Reserved.
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
//
// This holds the encoding/decoding routines that used to live in netutil

#ifndef UTIL_CODING_CODER_H__
#define UTIL_CODING_CODER_H__

#include <string.h>

// Avoid adding expensive includes here.
#include <glog/logging.h>

#include "base/integral_types.h"
#include "base/macros.h"
#include "base/port.h"
#include "base/type_traits.h"
#include "util/coding/varint.h"
#include "util/endian/endian.h"

/* Class for encoding data into a memory buffer */
class Encoder {
 public:
  // Creates an empty Encoder with no room that is enlarged
  // (if necessary) when "Encoder::Ensure(N)" is called.
  Encoder();
  ~Encoder();

  // Initialize encoder to encode into "buf"
  Encoder(void* buf, size_t maxn);
  void reset(void* buf, size_t maxn);
  void clear();

  // Encoding routines.  Note that these do not check bounds
  void put8(unsigned char v);
  void put16(uint16 v);
  void put32(uint32 v);
  void put64(uint64 v);
  void putword(uword_t v);
  void putn(const void* mem, size_t n);

  // Put no more than n bytes, stopping when c is put.
  void putcn(const void* mem, int c, size_t n);

  void puts(const void* mem);                // put a c-string including \0
  void puts_without_null(const char* mem);   // put a c-string without \0
  void putfloat(float f);
  void putdouble(double d);

  // Support for variable length encoding with 7 bits per byte
  // (these are just simple wrappers around the Varint module)
  static const int kVarintMax32 = Varint::kMax32;
  static const int kVarintMax64 = Varint::kMax64;

  void put_varint32(uint32 v);
  void put_varint64(uint64 v);
  static int varint32_length(uint32 v);  // Length of var encoding of "v"
  static int varint64_length(uint64 v);  // Length of var encoding of "v"

  // DEPRECATED
  //
  // For new code use put_varint32(ZigZagEncode(signed_value));
  // ZigZag coding is defined in utils/coding/transforms.h
  void put_varsigned32(int32 v);


  // Return number of bytes encoded so far
  size_t length() const;

  // Return number of bytes of space remaining in buffer
  size_t avail() const;

  // REQUIRES: Encoder was created with the 0-argument constructor interface.
  //
  // This interface ensures that at least "N" more bytes are available
  // in the underlying buffer by resizing the buffer (if necessary).
  //
  // Note that no bounds checking is done on any of the put routines,
  // so it is the client's responsibility to call Ensure() at
  // appropriate intervals to ensure that enough space is available
  // for the data being added.
  void Ensure(size_t N);

  // Returns true if Ensure is allowed to be called on "this"
  bool ensure_allowed() const { return underlying_buffer_ != NULL; }

  // Return ptr to start of encoded data.  This pointer remains valid
  // until reset or Ensure is called.
  const char* base() const { return (const char*)orig_; }

  // Advances the write pointer by "N" bytes.
  void skip(size_t N) { buf_ += N; }

  // REQUIRES: length() >= N
  // Removes the last N bytes out of the encoded buffer
  void RemoveLast(size_t N);

  // REQUIRES: length() >= N
  // Removes the last length()-N bytes to make the encoded buffer have length N
  void Resize(size_t N);

 private:
  void EnsureSlowPath(size_t N);

  unsigned char* orig_;
  unsigned char* buf_;
  unsigned char* limit_;

  // If constructed with the zero-argument constructor, we're allowed
  // to use Ensure; otherwise we're not.  If Ensure is allowed,
  // underlying_buffer_ is non-NULL; otherwise it is set to NULL.
  unsigned char* underlying_buffer_;

  static unsigned char kEmptyBuffer;

  DISALLOW_COPY_AND_ASSIGN(Encoder);
};

/* Class for decoding data from a memory buffer */
class Decoder {
 public:
  // Empty constructor to create uninitialized decoder
  inline Decoder() { }

  // NOTE: for efficiency reasons, this is not virtual.  so don't add
  // any members that really need to be destructed, and be careful about
  // inheritance.
  ~Decoder() { }

  // Initialize decoder to decode from "buf"
  Decoder(const void* buf, size_t maxn);
  void reset(const void* buf, size_t maxn);

  // Decoding routines.  Note that these do not check bounds
  unsigned char  get8();
  uint16 get16();
  uint32 get32();
  uint64 get64();
  uword_t getword();
  float  getfloat();
  double getdouble();
  void   getn(void* mem, size_t n);
  void   getcn(void* mem, int c, size_t n);    // get no more than n bytes,
                                               // stopping after c is got
  void   gets(void* mem, size_t n);            // get a c-string no more than
                                               // n bytes. always appends '\0'
  void   skip(size_t n);
  unsigned char const* ptr();       // Return ptr to current position in buffer

  // "get_varint" actually checks bounds
  bool get_varint32(uint32* v);
  bool get_varint64(uint64* v);

  // DEPRECATED
  //
  // For new code use
  //   get_varint32(&unsigned_temp);
  //   signed_value = ZigZagDecode(unsigned_temp);
  // ZigZag coding is defined in utils/coding/transforms.h
  bool get_varsigned32(int32* v);


  size_t pos() const;
  // Return number of bytes decoded so far

  size_t avail() const;
  // Return number of available bytes to read

 private:
  friend class IndexBlockDecoder;
  const unsigned char* orig_;
  const unsigned char* buf_;
  const unsigned char* limit_;
};

/***** Implementation details.  Clients should ignore them. *****/

inline Encoder::Encoder(void* b, size_t maxn) {
  orig_ = buf_ = reinterpret_cast<unsigned char*>(b);
  limit_ = orig_ + maxn;
  underlying_buffer_ = NULL;
}

inline void Encoder::reset(void* b, size_t maxn) {
  orig_ = buf_ = reinterpret_cast<unsigned char*>(b);
  limit_ = orig_ + maxn;
  // Can't use the underlying buffer anymore
  if (underlying_buffer_ != &kEmptyBuffer) {
    delete[] underlying_buffer_;
  }
  underlying_buffer_ = NULL;
}

inline void Encoder::clear() {
  buf_ = orig_;
}

inline void Encoder::Ensure(size_t N) {
  DCHECK(ensure_allowed());
  if (avail() < N) {
    EnsureSlowPath(N);
  }
}

inline size_t Encoder::length() const {
  DCHECK_GE(buf_, orig_);
  return buf_ - orig_;
}

inline size_t Encoder::avail() const {
  DCHECK_GE(limit_, buf_);
  return limit_ - buf_;
}

inline void Encoder::putn(const void* src, size_t n) {
  memcpy(buf_, src, n);
  buf_ += n;
}

inline void Encoder::putcn(const void* src, int c, size_t n) {
  unsigned char *old = buf_;
  buf_ = static_cast<unsigned char *>(memccpy(buf_, src, c, n));
  if (buf_ == NULL)
    buf_ = old + n;
}

inline void Encoder::puts(const void* src) {
  putcn(src, '\0', avail());
}

inline void Encoder::puts_without_null(const char* mem) {
  while (*mem != '\0' && buf_ < limit_) {
    *buf_++ = *mem++;
  }
}

inline void Encoder::put_varint32(uint32 v) {
  buf_ = reinterpret_cast<unsigned char*>
         (Varint::Encode32(reinterpret_cast<char*>(buf_), v));
}

inline void Encoder::put_varint64(uint64 v) {
  buf_ = reinterpret_cast<unsigned char*>
         (Varint::Encode64(reinterpret_cast<char*>(buf_), v));
}

// DEPRECATED
//
// For new code use put_varint32(ZigZagEncode(signed_value));
// ZigZag coding is defined in utils/coding/transforms.h
inline void Encoder::put_varsigned32(int32 n) {
  // Encode sign in low-bit
  int sign = (n < 0) ? 1 : 0;
  uint32 mag = (n < 0) ? -n : n;
  put_varint32((mag << 1) | sign);
}

inline Decoder::Decoder(const void* b, size_t maxn) {
  orig_ = buf_ = reinterpret_cast<const unsigned char*>(b);
  limit_ = orig_ + maxn;
}

inline void Decoder::reset(const void* b, size_t maxn) {
  orig_ = buf_ = reinterpret_cast<const unsigned char*>(b);
  limit_ = orig_ + maxn;
}

inline size_t Decoder::pos() const {
  DCHECK_GE(buf_, orig_);
  return buf_ - orig_;
}

inline size_t Decoder::avail() const {
  DCHECK_GE(limit_, buf_);
  return limit_ - buf_;
}

inline void Decoder::getn(void* dst, size_t n) {
  memcpy(dst, buf_, n);
  buf_ += n;
}

inline void Decoder::getcn(void* dst, int c, size_t n) {
  void *ptr;
  ptr = memccpy(dst, buf_, c, n);
  if (ptr == NULL)
    buf_ = buf_ + n;
  else
    buf_ = buf_ + (reinterpret_cast<unsigned char *>(ptr) -
                   reinterpret_cast<unsigned char *>(dst));
}

inline void Decoder::gets(void* dst, size_t n) {
  size_t len = n - 1;
  DCHECK_GE(limit_, buf_);
  if (n > 1 + limit_ - buf_) {
    len = limit_ - buf_;
  }
  (reinterpret_cast<char *>(dst))[len] = '\0';
  getcn(dst, '\0', len);
}

inline void Decoder::skip(size_t n) {
  buf_ += n;
}

inline unsigned char const* Decoder::ptr() {
  return buf_;
}


// DEPRECATED
//
// For new code use
//   get_varint32(&unsigned_temp);
//   signed_value = ZigZagDecode(unsigned_temp);
// ZigZag coding is defined in utils/coding/transforms.h
inline bool Decoder::get_varsigned32(int32* v) {
  uint32 coding;
  if (get_varint32(&coding)) {
    int sign = coding & 1;
    int32 mag = coding >> 1;
    if (sign) {
      // Special handling for encoding of kint32min
      *v = (mag == 0) ? kint32min : -mag;
    } else {
      *v = mag;
    }
    return true;
  } else {
    return false;
  }
}

inline void Encoder::put8(unsigned char v) {
  DCHECK_GE(avail(), sizeof(v));
  *buf_ = v;
  buf_ += sizeof(v);
}

inline void Encoder::put16(uint16 v) {
  DCHECK_GE(avail(), sizeof(v));
  LittleEndian::Store16(buf_, v);
  buf_ += sizeof(v);
}

inline void Encoder::put32(uint32 v) {
  DCHECK_GE(avail(), sizeof(v));
  LittleEndian::Store32(buf_, v);
  buf_ += sizeof(v);
}

inline void Encoder::put64(uint64 v) {
  DCHECK_GE(avail(), sizeof(v));
  LittleEndian::Store64(buf_, v);
  buf_ += sizeof(v);
}

inline void Encoder::putword(uword_t v) {
#ifdef _LP64
  LittleEndian::Store64(buf_, v);
#else
  LittleEndian::Store32(buf_, v);
#endif /* _LP64 */
  buf_ += sizeof(v);
}


inline void Encoder::putfloat(float f) {
  uint32 v;
  typedef char VerifySizesAreEqual[sizeof(f) == sizeof(v) ? 1 : -1]
    ATTRIBUTE_UNUSED;
  memcpy(&v, &f, sizeof(f));
  put32(v);
}

inline void Encoder::putdouble(double d) {
  uint64 v;
  typedef char VerifySizesAreEqual[sizeof(d) == sizeof(v) ? 1 : -1]
    ATTRIBUTE_UNUSED;
  memcpy(&v, &d, sizeof(d));
  put64(v);
}

inline unsigned char Decoder::get8() {
  const unsigned char v = *buf_;
  buf_ += sizeof(v);
  return v;
}

inline uint16 Decoder::get16() {
  const uint16 v = LittleEndian::Load16(buf_);
  buf_ += sizeof(v);
  return v;
}

inline uint32 Decoder::get32() {
  const uint32 v = LittleEndian::Load32(buf_);
  buf_ += sizeof(v);
  return v;
}

inline uint64 Decoder::get64() {
  const uint64 v = LittleEndian::Load64(buf_);
  buf_ += sizeof(v);
  return v;
}

inline uword_t Decoder::getword() {
#ifdef _LP64
  const uword_t v = LittleEndian::Load64(buf_);
#else
  const uword_t v = LittleEndian::Load32(buf_);
#endif /* _LP64 */
  buf_ += sizeof(v);
  return v;
}


inline float Decoder::getfloat() {
  uint32 v = get32();
  float f;
  typedef char VerifySizesAreEqual[sizeof(f) == sizeof(v) ? 1 : -1]
    ATTRIBUTE_UNUSED;
  memcpy(&f, &v, sizeof(f));
  return f;
}

inline double Decoder::getdouble() {
  uint64 v = get64();
  double d;
  typedef char VerifySizesAreEqual[sizeof(d) == sizeof(v) ? 1 : -1]
    ATTRIBUTE_UNUSED;
  memcpy(&d, &v, sizeof(d));
  return d;
}

#endif  // UTIL_CODING_CODER_H__

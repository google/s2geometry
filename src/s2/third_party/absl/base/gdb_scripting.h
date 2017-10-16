// Copyright 2017 The Abseil Authors.
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//      http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.
//
// DEFINE_GDB_AUTO_SCRIPT is intended to be added to source files (or headers)
// of data structures that have accompanying pretty-printers or other related
// support for gdb.
// It adds a reference to the script to the .debug_gdb_scripts section of the
// resulting executable or shared library.  GDB >= 7.2 automagically loads
// scripts mentioned in this section.
//
// NOTE: The definition is done in a way that makes the linker remove
// duplicates, so there is no harm to the resulting binary/library when
// using this macro in a header.
//
// Example:
// #include "third_party/absl/base/gdb_scripting.h"
// DEFINE_GDB_AUTO_SCRIPT("mylib/mylib_gdb.py")
//
// NOTE FOR GOOGLERS:
//
// Overview docs can be found at go/gdbscripts
//
// IWYU pragma: private, include "base/gdb-scripting.h"

#ifndef S2_THIRD_PARTY_ABSL_BASE_GDB_SCRIPTING_H_
#define S2_THIRD_PARTY_ABSL_BASE_GDB_SCRIPTING_H_

// The assembler syntax has only been tested with linux on the
// architectures below.  So only support those for now.
#if defined(__GNUC__) && defined(__linux__) &&      \
  (defined(__i386__) || defined(__amd64__) ||       \
   defined(__arm__) || defined(__aarch64__) ||      \
   defined(__powerpc__) || defined(__powerpc64__))
#define DEFINE_GDB_AUTO_SCRIPT(script_name) \
  asm(".pushsection \".debug_gdb_scripts\", \"MS\",%progbits,1\n" \
      ".byte 1\n" \
      ".asciz \"" script_name "\"\n" \
      ".popsection \n");
#else
#define DEFINE_GDB_AUTO_SCRIPT(script_name)
#endif

#endif  // S2_THIRD_PARTY_ABSL_BASE_GDB_SCRIPTING_H_

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

// main() for tests, running benchmarks if specified.
// Use --benchmark_filter=. to run all benchmarks.
// https://github.com/google/benchmark/blob/main/docs/user_guide.md#running-a-subset-of-benchmarks

#include <algorithm>
#include <cstdlib>
#include <cstring>

// cpplint thinks these are C headers and should go before the C++ ones.
// NOLINTBEGIN(build/include_order)
#include <benchmark/benchmark.h>
#include <gtest/gtest.h>
// NOLINTEND(build/include_order)

int main(int argc, char *argv[]) {
  // Disable ASLR for benchmarks, but not tests.  `MaybeReenterWithoutASLR`
  // calls `execve`, so we need to call it before we change `argv`.
  const bool has_benchmark_filter =
      std::any_of(&argv[1], &argv[argc], [](const char* arg) {
        // `string_view::starts_with()` is C++20, so use `strncmp`.
        const char prefix[] = "--benchmark_filter";
        return strncmp(arg, prefix, sizeof(prefix) - 1) == 0;
      });
  if (has_benchmark_filter)
    benchmark::MaybeReenterWithoutASLR(argc, argv);

  benchmark::Initialize(&argc, argv);
  testing::InitGoogleTest(&argc, argv);

  // If benchmarks were requested, just run them, then exit.
  if (!benchmark::GetBenchmarkFilter().empty()) {
    benchmark::RunSpecifiedBenchmarks();
    std::exit(EXIT_SUCCESS);
  }

  return RUN_ALL_TESTS();
}

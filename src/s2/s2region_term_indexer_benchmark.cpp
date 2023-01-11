#include "s2region_term_indexer.h"
#include "s2testing.h"
#include <benchmark/benchmark.h>

static void BM_OldApiGetIndexTerms(benchmark::State &state) {
  S2RegionTermIndexer indexer;
  benchmark::DoNotOptimize(indexer);
  auto const count = state.range(0);
  for (auto _ : state) {
    for (int64_t i = 0; i < count; ++i) {
      auto value = indexer.GetIndexTerms(S2Testing::RandomPoint(), {});
      benchmark::DoNotOptimize(value);
    }
  }
}

BENCHMARK(BM_OldApiGetIndexTerms)->RangeMultiplier(4)->Range(1, 65536);

static void BM_NewApiGetIndexTerms(benchmark::State &state) {
  S2RegionTermIndexer indexer;
  benchmark::DoNotOptimize(indexer);
  auto const count = state.range(0);
  for (auto _ : state) {
    std::vector<std::string> terms;
    benchmark::DoNotOptimize(terms);
    for (int64_t i = 0; i < count; ++i) {
      terms.clear();
      indexer.GetIndexTerms(S2Testing::RandomPoint(), {}, &terms);
      benchmark::DoNotOptimize(terms);
    }
  }
}

BENCHMARK(BM_NewApiGetIndexTerms)->RangeMultiplier(4)->Range(1, 65536);

static void BM_OldApiGetQueryTerms(benchmark::State &state) {
  S2RegionTermIndexer indexer;
  benchmark::DoNotOptimize(indexer);
  auto const count = state.range(0);
  for (auto _ : state) {
    for (int64_t i = 0; i < count; ++i) {
      auto value = indexer.GetQueryTerms(S2Testing::RandomPoint(), {});
      benchmark::DoNotOptimize(value);
    }
  }
}

BENCHMARK(BM_OldApiGetQueryTerms)->RangeMultiplier(4)->Range(1, 65536);

static void BM_NewApiGetQueryTerms(benchmark::State &state) {
  S2RegionTermIndexer indexer;
  benchmark::DoNotOptimize(indexer);
  auto const count = state.range(0);
  for (auto _ : state) {
    std::vector<std::string> terms;
    benchmark::DoNotOptimize(terms);
    for (int64_t i = 0; i < count; ++i) {
      terms.clear();
      indexer.GetQueryTerms(S2Testing::RandomPoint(), {}, &terms);
      benchmark::DoNotOptimize(terms);
    }
  }
}

BENCHMARK(BM_NewApiGetQueryTerms)->RangeMultiplier(4)->Range(1, 65536);

BENCHMARK_MAIN();
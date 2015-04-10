// Copyright 2009 Google Inc. All Rights Reserved.
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
// Author: jyrki@google.com (Jyrki Alakuijala)
//
// Tests for the fast hashing algorithm based on Austin Appleby's
// MurmurHash 2.0 algorithm. See http://murmurhash.googlepages.com/

#include "util/hash/murmur.h"

#include <string.h>
#include <string>

#include <glog/logging.h>

#include "base/integral_types.h"
#include "gtest/gtest.h"
#include "util/hash/fingerprint96.h"
#include "util/hash/fingerprinting.h"
#include "util/hash/hasheval/string_hash_eval.h"
#include "util/random/acmrandom.h"


namespace util_hash {

TEST(Murmur, EmptyData64) {
  EXPECT_EQ(0ULL, MurmurHash64(NULL, 0ULL));
  EXPECT_EQ(0ULL, MurmurHash64WithSeed(NULL, 0ULL, 0));
}

TEST(MurmurCat, EmptyData64) {
  MurmurCat cat;
  cat.Init(0ULL, 0ULL);
  cat.Append(NULL, 0ULL);
  EXPECT_EQ(0ULL, cat.GetHash());
  EXPECT_EQ(0ULL, MurmurHash64WithSeed(NULL, 0ULL, 0));
}

TEST(Murmur, VaryWithDifferentSeeds) {
  // While in theory different seeds could return the same
  // hash for the same data this is unlikely.
  char data1 = 'x';
  EXPECT_NE(MurmurHash64WithSeed(&data1, 1, 100),
            MurmurHash64WithSeed(&data1, 1, 101));
}

TEST(MurmurCat, VaryWithDifferentSeeds) {
  MurmurCat cat1;
  cat1.Init(100ULL, 1ULL);
  cat1.Append("x", 1);

  MurmurCat cat2;
  cat2.Init(101ULL, 1ULL);
  cat2.Append("x", 1);

  EXPECT_NE(cat1.GetHash(), cat2.GetHash());
}

// Hashes don't change.
TEST(Murmur, Idempotence) {
  const char data[] = "deadbeef";
  const size_t dlen = strlen(data);
  const uint64 orig64 = MurmurHash64(data, dlen);

  EXPECT_EQ(MurmurHash64(data, dlen), MurmurHash64(data, dlen));

  for (int i = 0; i < 10; i++) {
    EXPECT_EQ(MurmurHash64WithSeed(data, dlen, i),
              MurmurHash64WithSeed(data, dlen, i));
  }

  const char next_data[] = "deadbeef000---";
  const size_t next_dlen = strlen(next_data);

  EXPECT_EQ(MurmurHash64(next_data, next_dlen),
            MurmurHash64(next_data, next_dlen));

  for (int i = 0; i < 10; i++) {
    EXPECT_EQ(MurmurHash64WithSeed(next_data, next_dlen, i),
              MurmurHash64WithSeed(next_data, next_dlen, i));
  }

  // Go back to the first test data ('data') and make sure it hasn't changed.
  EXPECT_EQ(MurmurHash64(data, dlen), orig64);
}

TEST(MurmurCat, CompatibilityWithMurmurHash) {
  ACMRandom rnd(ACMRandom::DeterministicSeed());
  for (int i = 0; i < 1000; ++i) {
    string data(i + 1, 0);
    InitializeRandomString(data.size(), &data[0], &rnd);
    uint64 expected = MurmurHash64(data.data(), data.size());

    MurmurCat cat;
    cat.Init(0ULL, data.size());
    cat.Append(data.data(), data.size());
    EXPECT_EQ(expected, cat.GetHash());
  }
}

TEST(MurmurCat, CompatibilityWithMurmurHashSplit2) {
  ACMRandom rnd(ACMRandom::DeterministicSeed());
  for (int i = 0; i < 1000; ++i) {
    string data(i + 9, 0);
    InitializeRandomString(data.size(), &data[0], &rnd);
    uint64 expected = MurmurHash64(data.data(), data.size());

    MurmurCat split_cat;
    split_cat.Init(0ULL, data.size());
    int split = 1 + (i % (data.size() - 2));
    split_cat.Append(data.data(), split);
    split_cat.Append(data.data() + split, data.size() - split);
    EXPECT_EQ(expected, split_cat.GetHash());
  }
}

TEST(MurmurCat, CompatibilityWithMurmurHashSplit3) {
  ACMRandom rnd(ACMRandom::DeterministicSeed());
  for (int i = 0; i < 1000; ++i) {
    string data(i + 9, 0);
    InitializeRandomString(data.size(), &data[0], &rnd);
    uint64 expected = MurmurHash64(data.data(), data.size());

    MurmurCat triple_cat;
    triple_cat.Init(0ULL, data.size());
    int split0 = i & 7;
    int split1 = 8 + (i % (data.size() - 2));
    triple_cat.Append(data.data(), split0);
    triple_cat.Append(data.data() + split0, split1 - split0);
    triple_cat.Append(data.data() + split1, data.size() - split1);
    EXPECT_EQ(expected, triple_cat.GetHash());
  }
}

}  // namespace util_hash

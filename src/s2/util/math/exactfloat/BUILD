load("@rules_cc//cc:defs.bzl", "cc_library", "cc_test")

package(default_visibility = ["//visibility:public"])

cc_library(
    name = "bignum",
    srcs = ["bignum.cc"],
    hdrs = ["bignum.h"],
    visibility = "//visibility:private",
    deps = [
        "@abseil-cpp//absl/algorithm:container",
        "@abseil-cpp//absl/container:inlined_vector",
        "@abseil-cpp//absl/log:absl_check",
    ],
)

cc_library(
    name = "exactfloat",
    srcs = [
        "exactfloat.cc",
    ],
    hdrs = [
        "exactfloat.h",
    ],
    deps = [
        ":bignum",
        "//s2/base:logging",
        "@abseil-cpp//absl/container:inlined_vector",
        "@abseil-cpp//absl/log:absl_check",
        "@abseil-cpp//absl/log:absl_log",
        "@abseil-cpp//absl/numeric:bits",
        "@abseil-cpp//absl/numeric:int128",
        "@abseil-cpp//absl/random:random",
        "@abseil-cpp//absl/strings:ascii",
        "@abseil-cpp//absl/strings:str_cat",
        "@abseil-cpp//absl/strings:str_format",
    ] + select({
        "//:use_openssl_setting": ["@openssl//:crypto"],
        "//conditions:default": ["@boringssl//:crypto"],
    }),
)

cc_test(
    name = "bignum_test",
    srcs = ["bignum_test.cc"],
    deps = [
        ":bignum",
        "@abseil-cpp//absl/random:bit_gen_ref",
        "@boringssl//:crypto",
        "@googletest//:gtest_main",
    ],
)

cc_test(
    name = "exactfloat_test",
    srcs = ["exactfloat_test.cc"],
    deps = [
        ":exactfloat",
        "//s2/util/math:vector",
        "@abseil-cpp//absl/log:absl_check",
        "@abseil-cpp//absl/log:absl_log",
        "@abseil-cpp//absl/random",
        "@googletest//:gtest_main",
    ],
)

cc_test(
    name = "exactfloat_underflow_test",
    srcs = ["exactfloat_underflow_test.cc"],
    deps = [
        ":exactfloat",
        "@googletest//:gtest_main",
    ],
)

load("@rules_cc//cc:defs.bzl", "cc_library", "cc_test")

package(default_visibility = ["//visibility:public"])

cc_library(
    name = "exactfloat",
    srcs = [
        "exactfloat.cc",
        "exactfloat_internal.cc",
    ],
    hdrs = [
        "exactfloat.h",
        "exactfloat_internal.h",
    ],
    deps = [
        "//s2/base:logging",
        "@abseil-cpp//absl/container:inlined_vector",
        "@abseil-cpp//absl/log:absl_check",
        "@abseil-cpp//absl/log:log",
        "@abseil-cpp//absl/numeric:bits",
        "@abseil-cpp//absl/numeric:int128",
        "@abseil-cpp//absl/strings:ascii",
        "@abseil-cpp//absl/strings:str_format",
    ] + select({
        "//:use_openssl_setting": ["@openssl//:crypto"],
        "//conditions:default": ["@boringssl//:crypto"],
    }),
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

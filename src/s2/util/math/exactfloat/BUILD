package(default_visibility = ["//visibility:public"])

cc_library(
    name = "exactfloat",
    srcs = ["exactfloat.cc"],
    hdrs = ["exactfloat.h"],
    deps = [
        "//s2/base:logging",
        "@abseil-cpp//absl/log:log",
        "@abseil-cpp//absl/log:absl_check",
        "@boringssl//:crypto",
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

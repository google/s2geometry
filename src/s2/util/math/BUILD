load("@rules_cc//cc:defs.bzl", "cc_library")

package(default_visibility = ["//visibility:public"])

cc_library(
    name = "vector",
    hdrs = ["vector.h"],
)

cc_library(
    name = "mathutil",
    srcs = ["mathutil.cc"],
    hdrs = ["mathutil.h"],
    deps = [
        "//:s2-includes",
        "//s2/base:integral_types",
        "//s2/base:logging",
        "//s2/util/bits",
        "@abseil-cpp//absl/log",
        "@abseil-cpp//absl/log:check",
    ],
)

cc_library(
    name = "matrix3x3",
    hdrs = ["matrix3x3.h"],
)

load("@rules_cc//cc:defs.bzl", "cc_library")

package(default_visibility = ["//visibility:public"])

cc_library(
    name = "exactfloat",
    srcs = ["exactfloat.cc"],
    hdrs = ["exactfloat.h"],
    deps = [
        "//:s2-includes",
        "//s2/base:integral_types",
        "//s2/base:logging",
        "//s2/base:port",
        "@abseil-cpp//absl/log",
        "@boringssl//:crypto",
    ],
)

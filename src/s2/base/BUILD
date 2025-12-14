load("@rules_cc//cc:defs.bzl", "cc_library")

package(default_visibility = ["//visibility:public"])

cc_library(
    name = "logging",
    hdrs = [
        "log_severity.h",
    ],
    deps = [
        "@abseil-cpp//absl/log",
        "@abseil-cpp//absl/log:absl_log",
        "@abseil-cpp//absl/log:check",
    ],
)

cc_library(
    name = "malloc_extension",
    srcs = ["malloc_extension.cc"],
    hdrs = ["malloc_extension.h"],
    deps = [
        "@abseil-cpp//absl/base:core_headers",
    ],
)

cc_library(
    name = "casts",
    hdrs = ["casts.h"],
    deps = [
        "@abseil-cpp//absl/log",
    ],
)

cc_library(
    name = "timer",
    hdrs = [
        "timer.h",
    ],
)

cc_library(
    name = "spinlock",
    hdrs = [
        "spinlock.h",
    ],
    deps = [
        "@abseil-cpp//absl/base:core_headers",
    ],
)

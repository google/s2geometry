load("@rules_cc//cc:defs.bzl", "cc_library")

package(default_visibility = ["//visibility:public"])

cc_library(
    name = "logging",
    hdrs = [
        "log_severity.h",
        "logging.h",
    ],
    deps = [
        "@abseil-cpp//absl/log",
        "@abseil-cpp//absl/log:absl_log",
        "@abseil-cpp//absl/log:check",
    ],
)

cc_library(
    name = "types",
    hdrs = ["types.h"],
)

cc_library(
    name = "casts",
    hdrs = ["casts.h"],
    deps = [
        "@abseil-cpp//absl/log",
    ],
)

cc_library(
    name = "port",
    hdrs = ["port.h"],
)

cc_library(
    name = "commandlineflags",
    hdrs = [
        "commandlineflags.h",
        "commandlineflags_declare.h",
    ],
    deps = [
        ":types",
        "@abseil-cpp//absl/flags:flag",
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
)

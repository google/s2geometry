package(default_visibility = ["//visibility:public"])

cc_library(
    name = "exactfloat",
    srcs = ["exactfloat.cc"],
    hdrs = ["exactfloat.h"],
    deps = [
        "//s2/base:types",
        "//s2/base:port",
        "//s2/base:logging",
        "@abseil-cpp//absl/log:log",
        "@abseil-cpp//absl/log:absl_check",
        "@boringssl//:crypto",
    ],
)

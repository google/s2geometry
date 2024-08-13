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
        "//s2/base:types",
        "//s2/util/bits",
        "@abseil-cpp//absl/log",
        "@abseil-cpp//absl/log:check",
        "@abseil-cpp//absl/log:absl_check",
    ],
)

cc_library(
    name = "matrix3x3",
    hdrs = ["matrix3x3.h"],
)

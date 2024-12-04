package(default_visibility = ["//visibility:public"])

cc_library(
    name = "coding",
    srcs = [
        "coder.cc",
        "varint.cc",
    ],
    hdrs = [
        "coder.h",
        "nth-derivative.h",
        "transforms.h",
        "varint.h",
    ],
    deps = [
        "//s2/base:casts",
        "//s2/base:types",
        "//s2/base:logging",
        "//s2/base:malloc_extension",
        "//s2/base:port",
        "//s2/util/bits",
        "//s2/util/endian",
        "@abseil-cpp//absl/base",
        "@abseil-cpp//absl/container:btree",
        "@abseil-cpp//absl/container:flat_hash_set",
        "@abseil-cpp//absl/log:absl_check",
        "@abseil-cpp//absl/strings:cord",
    ],
)

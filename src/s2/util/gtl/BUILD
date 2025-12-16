load("@rules_cc//cc:defs.bzl", "cc_library")

package(default_visibility = ["//visibility:public"])

cc_library(
    name = "gtl",
    hdrs = glob(["*.h"]),
    deps = [
        "//s2/base:malloc_extension",
    ],
)

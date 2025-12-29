load("@rules_cc//cc:defs.bzl", "cc_library")

package(default_visibility = ["//visibility:public"])

cc_library(
    name = "test_main",
    srcs = ["test_main.cc"],
    deps = [
        "@google_benchmark//:benchmark",
        "@googletest//:gtest",
    ],
)

cc_library(
    name = "testing",
    hdrs = ["gtest_prod.h"],
)

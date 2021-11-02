load("@bazel_tools//tools/build_defs/repo:http.bzl", "http_archive")
load("@bazel_tools//tools/build_defs/repo:git.bzl", "git_repository")


git_repository(
    name = "dm_core_cpp",
    commit = "282bc312d06aec790e98e4495af9b43478c42ce2",
    remote = "https://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@github.com/deepmirrorinc/CoreCpp.git",
)

load("@dm_core_cpp//bazel:corecpp_deps.bzl", "corecpp_deps")

corecpp_deps()

# grpc
# Change the boringssl url from googlesource.com used by grpc to github.com.
http_archive(
    name = "boringssl",
    sha256 = "6f17e41a24f3327ccd174ce3435a5a5b3f36fc9999ddcc0a5b616d64bc2b214c",
    strip_prefix = "boringssl-2c9ef59078c8ad527158fa73b02d517d1b4044a7",
    url = "file:///home/tmp/boringssl-2c9ef59078c8ad527158fa73b02d517d1b4044a7.zip",
)

#
# same as the absl version in grpc
# https://github.com/grpc/grpc/blob/v1.38.1/third_party/upb/bazel/workspace_deps.bzl
git_repository(
    name = "com_google_absl",
    commit = "df3ea785d8c30a9503321a3d35ee7d35808f190d",  # LTS 2020-02-25
    remote = "https://github.com/abseil/abseil-cpp.git",
    shallow_since = "1583355457 -0500",
)

# latest absl version
#http_archive(
#  name = "com_google_absl",
#  urls = ["https://github.com/abseil/abseil-cpp/archive/98eb410c93ad059f9bba1bf43f5bb916fc92a5ea.zip"],
#  strip_prefix = "abseil-cpp-98eb410c93ad059f9bba1bf43f5bb916fc92a5ea",
#)

http_archive(
  name = "rules_cc",
  urls = ["https://github.com/bazelbuild/rules_cc/archive/262ebec3c2296296526740db4aefce68c80de7fa.zip"],
  strip_prefix = "rules_cc-262ebec3c2296296526740db4aefce68c80de7fa",
)

http_archive(
  name = "com_google_googletest",
  urls = ["https://github.com/google/googletest/archive/011959aafddcd30611003de96cfd8d7a7685c700.zip"],
  strip_prefix = "googletest-011959aafddcd30611003de96cfd8d7a7685c700",
)

http_archive(
    name = "com_github_google_benchmark",
    urls = ["https://github.com/google/benchmark/archive/bf585a2789e30585b4e3ce6baf11ef2750b54677.zip"],
    strip_prefix = "benchmark-bf585a2789e30585b4e3ce6baf11ef2750b54677",
    sha256 = "2a778d821997df7d8646c9c59b8edb9a573a6e04c534c01892a40aa524a7b68c",
)

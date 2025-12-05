---
title: Java Benchmarks for S2
---

<!--*
# Document freshness: For more information, see go/fresh-source.
freshness: { owner: 'torrey' reviewed: '2022-02-07' }
*-->



## Overview

The Java implementation of S2 has a collection of benchmarks that measure
performance of many important parts of the library. Their intent is to provide
accurate data to guide work on performance improvements, and also to detect
unwanted performance regressions.

<!-- TODO(torrey): Document Perfgate integration.
Benchmark results are tracked over time on Perfgate, and can be viewed at ...

-->

The benchmark code can be found in
google3/javatests/com/google/common/geometry/benchmarks .


## Benchmark Organization

The Java S2 benchmarks are written using JMH, the Java Microbenchmark Harness
framework. Documentation on JMH is available from the [official site](https://github.com/openjdk/jmh).

The S2 benchmarks are in multiple files, one for each S2 class. For example,
benchmarks for S2CellId are in
[S2CellIdBenchmark.java](https://source.corp.google.com/piper///depot/google3/javatests/com/google/common/geometry/benchmarks/S2CellIdBenchmark.java).
Each benchmark has a corresponding Blaze build target, with a
`Benchmark_jmh` suffix, e.g. `S2CellIdBenchmark_jmh`.

Within each file, there are many individual JMH benchmarks, each of which has
the **@Benchmark** annotation. Most of these measures the performance of a
single method, while some may measure a sequence of operations on the class.
For example, the S2CellId Benchmark
[getAllNeighbors](https://source.corp.google.com/search?q=S2CellIdBenchmark%20getAllNeighbors&sq=package:piper%20file:%2F%2Fdepot%2Fgoogle3%20-file:google3%2Fexperimental)
measures the time to compute all neighbors of leaf cells.

In a normal JMH run, the framework forks a new process which starts a new Java
virtual machine. Within that process, it calls each Benchmark method as many
times as possible until a target time has elapsed. Each call to the method is
called an **Invocation**. The whole sequence of calls until the target time
elapses is called a benchmark **Iteration**.

JMH performs multiple Iterations for each Benchmark. The set of iterations is
called a **Trial**. When a trial completes, the JVM running it shuts down.
Multiple trials can be run for a benchmark (or all benchmarks) by increasing the
number of **Forks**. Every fork starts a new JVM and runs a complete Trial
within it. Statistics are computed across all Iterations in all Trials, or
forks, for the average time per invocation with a confidence interval.

(If the number of forks is set to zero, everything runs in the same single JVM,
with one trial for each benchmark. This approach is not recommended for getting
accurate results, but it is faster, and is used in benchmark unit tests.)

Before JMH runs the timed iterations in a Trial, it will run several **Warmup**
iterations. These run in the same way as the timed ones, but the results
are not used in computing the resulting statistics. Warmup iterations allow the
new JVM to fully optimize the code for maximum performance before the timed
iterations, which is important to get consistent, low-variance results. When you
look at the output of the run, you'll usually see something like the following
example, with the first few (warmup) iterations having somewhat lower
performance. This effect can be *very* dramatic for some code.

```
# Warmup Iteration   1: 34.030 ns/op
# Warmup Iteration   2: 29.849 ns/op
# Warmup Iteration   3: 27.266 ns/op
Iteration   1: 27.690 ns/op
Iteration   2: 27.213 ns/op
Iteration   3: 27.292 ns/op
Iteration   4: 27.351 ns/op
Iteration   5: 27.533 ns/op
```

If the first timed iteration of a benchmark is noticeably slower (a longer
average time per operation) than following ones, you should use more or longer
warmup iterations; otherwise, the final statistics will have lower accuracy and
a wider confidence interval.

Within the code, Benchmark methods are usually contained within a **@State**
object. Benchmarks typically iterate through a variety of test data to provided
by the State to ensure the code is run with a representative range of values.

The `State` also manages **@Parameters**. By default, JMH will run every
@Benchmark in a @State with every combination of @Parameters, and report the
results separately. Sometimes, you may not want to benchmark every possible
combination, as this can result in an explosion of cases. There are several
possible solutions to trim these number of cases. One such way is to specify
parameter values on the command line, as [described below](#commandline-parameters)

However, if you *never* want to explore the entire parameter space,
you can define a Java enum for the combinations of parameters you want to
benchmark, where the enum instances each contain a set of values which are the
actual parameters. A simpler solution is to throw an exception in @Setup for
unwanted combinations of Parameters. This is effective if there is a constraint
like "Param A must always be less than Param B".

## Running Benchmarks

Running benchmarks on a local workstation or CloudTop is convenient and
reasonably accurate during development, but you should avoid running other
significant tasks such as Chrome or an IDE at the same time. Run benchmarks
for one class like:

```
export B=javatests/com/google/common/geometry/benchmarks
blaze run -c opt $B:S2LoopBenchmark_jmh
```

[Various flags](#flags) can be passed to JMH to control how the benchmark runs.

For greater accuracy and repeatability, benchmarks can be run on dedicated
**Perflab** machines in Borg. To run on Perflab, you'll first need to complete
the quickstart instructions in the
[Perflab Tutorial](https://g3doc.corp.google.com/devtools/crosstool/perflab/g3doc/tutorial.md)
page.

Then, you can run a benchmark on Perflab as follows:

```
export B=javatests/com/google/common/geometry/benchmarks
blaze run -c opt --run_under='perflab --constraints=arch=x86_64' \
  $B:S2CellIdBenchmark_jmh
```

When the job starts, URLs will be printed where you can track the progress of
the run on Borg, and find the final results in CNS when the job completes.
The output files have a TTL, which is currently 60 days.

Note: you can ignore the message: `*Unable to locate a valid CitC path for CL
None, CitC workspace "None" and snapshot ID "None"*`. Under the covers,
`run_under` is using [mcd_perflab](go/perflab), which is designed for A/B
experiments. This message is just mcd_perflab saying it doesn't have a B side
for an A/B experiment.)

Add the `--perfgate` flag if you want to upload benchmark results to PerfGate,
like this:

```
export B=javatests/com/google/common/geometry/benchmarks
blaze run -c opt --run_under='perflab --constraints=arch=x86_64' \
  $B:S2CellIdBenchmark_jmh -- --perfgate
```

## More Examples

#### Run only one or more specific @Benchmarks for a class

This technique is useful if you're working on improving or understanding
performance of a single S2 method or something similar, and just want to
repeatedly test the effect of changes you are making with the most applicable
benchmark.

Add the name(s) of the @Benchmark methods to the end of the command line. The
strings you provide are actually matched as regular expressions to the fully
qualified method name. So, depending on the organization and naming of your
State and Benchmarks, you can select what you want in various ways.

For example, these three commands all do the same thing:

```
export B=javatests/com/google/common/geometry/benchmarks

blaze run -c opt $B:S2CellIdBenchmark_jmh com.google.common.geometry.benchmarks.S2CellIdBenchmark.EvenlySpacedLeafCellIdsState.toPointRaw
blaze run -c opt $B:S2CellIdBenchmark_jmh EvenlySpacedLeafCellIdsState.toPointRaw
blaze run -c opt $B:S2CellIdBenchmark_jmh toPointRaw
```

To run both the `toPointRaw` and `toPoint` benchmarks using a regular expression
match:

```
export B=javatests/com/google/common/geometry/benchmarks

blaze run -c opt $B:S2CellIdBenchmark_jmh toPoint*
```

Providing just the name of a State object will match all the benchmarks within
it, like

```
export B=javatests/com/google/common/geometry/benchmarks
blaze run -c opt $B:S2CellIdBenchmark_jmh EvenlySpacedLeafCellIdsState
```

You can also specify multiple benchmark names:

```
export B=javatests/com/google/common/geometry/benchmarks
blaze run -c opt $B:S2EdgeUtilBenchmark_jmh robustCrossing edgeCrosser100RobustCrossings
```

#### Run a specific @Benchmark with specific parameter values{#commandline-parameters}

This technique is useful (for example) if you are working to improve performance
for a specific scenario, perhaps large overlapping polygons, and don't need to
get results for other scenarios.

Use the `-p` flag, repeatedly if needed, and with comma-separated values for
each parameter. This example sets values for two parameters:

```
blaze run -c opt $B:S2ClosestEdgeQueryBenchmark_jmh \
  FindClosestToSmallAbuttingIndex.findClosestEdge -- -p numIndexEdges=1024,65536 -p factory=FRACTAL_LOOP
```

## JMH Benchmark Flags{#flags}

JMH has many command flags that can specify values for many of the settings
described above. It is convenient to provide defaults for these values using
annotations like @Warmup and @Measurement in the code, at the level of
individual @Benchmarks or @States, but those can be overridden from the command
line with the following flags:

| Flag Example     | Meaning                                                   |
| ---------------- | ----------------------------------------------------------|
| -wi=3            | Set the number of warmup iterations per trial.            |
| -i=5             | Set the number of timed benchmark iterations per trial.   |
| -f=2             | Set the number of forks, i.e. how many times JMH will fork
:                  : a new JVM process to run a set of benchmark trials.       :
| -w=2s            | Set the target time per warmup iteration to 2 seconds     |
| -r=10s           | Set the target time per benchmark iteration to 10 seconds |
| -to=10           | Set the timeout for iterations to 10 seconds. (note, no   |
:                  : trailing 's')                                             :
| -p arg=val1,val2 | Override the parameter values used for parameter "arg" to :
:                  : only the specific provided values. Can be repeated.       :

**Flags Example**

Note that blaze differentiates flags for the blaze binary itself from flags for
the run target using a double-hyphen.

To run the `S2ClosestEdgeQuery.FindClosestToSmallAbuttingIndex` benchmark with
five timed iterations, 10 seconds per measured iteration, the
`numIndexedEdges` @Param set to 64K, and other values using the values set in
the benchmark code:

```
blaze run -c opt $B:S2ClosestEdgeQueryBenchmark_jmh FindClosestToSmallAbuttingIndex -- -i=5 -r=10s -p numIndexEdges=65536
```

## Developing Benchmarks

**Test the benchmarks for a single class:**

Testing a benchmark with "blaze test" invokes each @Benchmark one time, for
every combination of parameters. The blaze target names are like
`test_S2LoopBenchmark_jmh`. For example:

```
export B=javatests/com/google/common/geometry/benchmarks
blaze test -c opt $B:test_S2ClosestEdgeQueryBenchmark_jmh
```

"Blaze test" is supposed to be fast, but if some combinations of parameters take
a long time to run even once, test timeouts can be a problem. The "isUnitTest()"
method can be used in @State `setup()` methods to replace expensive parameter
values with fast ones. See the
[S2PolygonBenchmark](https://source.corp.google.com/piper///depot/google3/javatests/com/google/common/geometry/benchmarks/S2PolygonBenchmark.java)
for some examples.

**Test all the benchmarks:**

```
export B=javatests/com/google/common/geometry/benchmarks
blaze test -c opt $B/...
```

## Comparing Java to C++ benchmarks

Many (but not all) of the Java S2 benchmarks were designed to do exactly the
same operations as corresponding C++ benchmarks. This allows for comparing the
performance of the Java and C++ implementations.

The C++ benchmarks are part of the unit tests. You can build one like:

```
blaze build -c opt --dynamic_mode=off --copt=-gmlt util/geometry:s2closest_edge_query_test
```

And run all the benchmarks in a unit test:

```
blaze-bin/util/geometry/s2closest_edge_query_test --benchmark_filter=all
```

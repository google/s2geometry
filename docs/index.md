---
title: S2 Geometry
---


<center><b>Welcome to the S2 Geometry Library!</b></center>
<br/>

![](devguide/img/s2hierarchy.gif)

A unique feature of the S2 library is that unlike traditional geographic
information systems, which represent data as flat two-dimensional projections
(similar to an atlas), the S2 library represents all data on a three-dimensional
sphere (similar to a globe). This makes it possible to build a worldwide
geographic database with no seams or singularities, using a single coordinate
system, and with low distortion everywhere compared to the true shape of the
Earth. While the Earth is not quite spherical, it is much closer to being a
sphere than it is to being flat!

If you want to learn more about the library, start by reading the
[Overview](about/overview) and
[QuickStart document](devguide/cpp/quickstart), then read the
introduction to the [basic types](devguide/basic_types)

**<a href="https://github.com/google/s2geometry" target="_blank">
Get S2 on GitHub</a>**

## S2 Features

Notable features of the library include:

* Flexible support for spatial indexing, including the ability
  to approximate arbitrary regions as collections of discrete
  S2 cells. This feature makes it easy to build large distributed
  spatial indexes.
* Fast in-memory spatial indexing of collections of points,
  polylines, and polygons.
* Robust constructive operations (such as intersection, union, and
  simplification) and boolean predicates (such as testing for
  containment).
* Efficient query operations for finding nearby objects, measuring
  distances, computing centroids, etc.
* A flexible and robust implementation of snap rounding (a geometric
  technique that allows operations to be implemented 100% robustly
  while using small and fast coordinate representations).
* A collection of efficient yet exact mathematical predicates for
  testing relationships among geometric primitives.
* Extensive testing on Google's vast collection of geographic data.
* Flexible Apache 2.0 license.

## Table of Contents

*   About S2

    *   [Overview](about/overview)
    *   [Platforms Guide](about/platforms)

*   Tutorials

    *   [Quick Start](devguide/cpp/quickstart)

*   Developer Guides

    *   [Basic Types](devguide/basic_types)
    *   [Finding Nearby Edges](devguide/s2closestedgequery)
    *   [S2Cell Hierarchy](devguide/s2cell_hierarchy)
    *   [Covering Examples](devguide/examples/coverings)

*   Resources

    *   The [S2 Earthcube](resources/earthcube)
    *   [S2Cell Statistics](resources/s2cell_statistics)

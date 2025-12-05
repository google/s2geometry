---
title: Degenerate Geometry in S2
---

## Overview

Unlike many libraries, S2 fully supports degenerate geometry.  Degeneracies not
only do not cause errors, they have well-defined semantics including support for
open, semi-open, and closed boundaries.  Supported degeneracy types include the
following:

 - *Point polylines* consisting of a single degenerate edge *AA* (where *A*
   represents a vertex).

 - *Point loops* consisting of a single vertex *A*.  Such loops may represent
   either *point shells* or *point holes* according to whether the loop adds to
   or subtracts from the surrounding region of the polygon.

 - *Sibling edge pairs* of the form \{*AB*, *BA*\}.  Such sibling pairs may
   represent either degenerate shells or degenerate holes according to whether
   they add to or subtract from the surrounding region.  The edges of a sibling
   pair may belong to the same polygon loop (e.g. a loop *AB*) or to different
   polygon loops (e.g. the polygon consisting of the loops *ABC* and *CBD*).

Supporting such degeneracies has one big advantage, namely that geometry can be
simplified without completely losing small details.  For example consider a
polygon representing a land region with many small lakes.  Each lake would be
represented as a hole.  When this polygon is simplified, some or even all of the
lakes may collapse to single points.  Without support for degeneracies, such
lakes would disappear from the simplified result.  If we allow degeneracies, on
the other hand, we can guarantee that for every point in the original region
there is a nearby corresponding point in the simplified region and vice versa.
In other words every lake in the input is still near some lake in the output.

## Dimension Reduction: A Poor Substitute {#dimension-reduction}

An alternative technique for handling degeneracies is *dimension reduction*,
whereby degenerate geometry of dimension 2 or 1 is replaced by non-degenerate
geometry of dimension 1 or 0.  For example, the sibling pair \{*AB*, *BA*\}
would be replaced by a polyline edge (either *AB* or *BA*, chosen arbitrarily),
and similarly the degenerate polyline *AA* would be replaced by a single point
*A*.

This technique has four very significant drawbacks compared to true degeneracy
support as implemented by the S2 library:

 - **Dimension reduction is only capable of handling** ***positive
   degeneracies*** (e.g., polygon shells) as opposed to *negative degeneracies*
   (e.g., polygon holes).  Consider again the example of simplifying a land
   region with many small lakes that collapse to single points.  Since these
   lakes are represented as holes in the polygon rather than shells, they cannot
   be replaced by points but instead must be simply discarded.  Similarly,
   dimension reduction is not capable of representing negative degeneracies of
   dimension 1 (e.g., a land polygon containing a river that when simplified has
   its width collapse to zero).

 - **Dimension reduction changes the boundary of the affected geometry.** In
   general the boundary of a polygon is defined to consist of its edges, the
   boundary of a polyline is defined to consist of its two endpoints, and points
   are defined to have no boundary at all.  This implies that when the
   degenerate sibling pair \{*AB*, *BA*\} is replaced by the polyline edge *AB*,
   for example, the definition of its boundary changes from the two edges
   \{*AB*, *BA*\} to the two points \{*A*, *B*\}.  This is a very significant
   difference, since it can change the result of virtually all the spatial
   relationship predicates defined by OGC Simple Features model (which are based
   on evaluating whether the interior, boundary, and exterior of one region
   intersect the interior, boundary, and exterior of the other as summarized by
   the `DE-9IM` matrix).

 - **Dimension reduction forces clients to deal with heterogeneous geometry as a
   possible result of virtually any geometric operation.** For example,
   intersecting two polygons may yield a collection of geometry consisting of
   polygons, polylines, and points.  Similarly, intersecting two polylines may
   yield a collection of polylines and points.  In contrast, the true degeneracy
   support offered by S2 means that any boolean operation on two polygons is
   guaranteed to yield a polygon, and any boolean operation on two polylines is
   guaranteed to yield a polyline.  This is not only simpler for clients but
   also reduces the opportunities for bugs (since the cases resulting in
   mixed-dimension geometry may be unexpected and rare).

 - **Dimension reduction makes it hard for clients to say how degeneracies
   should be handled.** Some clients may wish to discard degeneracies, generate
   errors, or convert degeneracies back into non-degenerate forms (e.g.
   replacing point shells with tiny triangles).  True degeneracy support makes
   it easy to implement any of these models, and even to delegate this decision
   to the functions that convert geometry to non-degeneracy-supporting formats
   (such as geoJSON).  Dimension reduction, on the other hand, makes it
   difficult to distinguish degeneracies from legitimate lower-dimensional
   objects.  For example, imagine intersecting two collections of polylines and
   polygons.  Some output polylines will correspond to portions of input
   polylines, while others may correspond to intersecting polygons that abut
   along an edge.

Other alternatives, such as extracting small details into separate objects to
protect them from simplification and then merging them back in afterwards, can
best be described as hacks.  Not only are such algorithms virtually impossible
to make robust, they cannot offer the mathematical guarantees that can be
achieved by simplifying with a uniform error tolerance using degeneracies (see
below).

## Degeneracies and Hausdorff Distance

The fact that S2 supports both positive and negative degeneracies allows us to
make a very powerful guarantee when simplifying geometry.  Recall that the
*Hausdorff distance* between two geometries *A* and *B* is defined as

$$H(A, B) = {\rm max} \{ h(A, B), h(B, A) \}$$

where $$h(A, B)$$ is the *directed Hausdorff distance*

$$h(A, B) = {\rm max}_{a \in A} \big\{ {\rm min}_{b \in B} \{ d(a, b) \} \big\}
. $$

**The Hausdorff distance is a measure of how different two geometries are;** in
particular it is an upper bound on how far away a point in one geometry can be
from the other geometry.

The most important property of Hausdorff distance is that it can be used to
bound the distance measurement errors due to simplifying geometry.  Suppose that
geometry *B* is a simplified version of geometry *A*, and that we wish to
measure the distance to some third geometry *X*.  Then the error due to
simplification, i.e. the difference between $$d(A, X)$$ and $$d(B, X)$$, is
bounded by the Hausdorff distance between *A* and *B*:

$$| d(A, X) - d(B, X) | \leq H(A, B) \, .$$

We can now state the guarantee that the S2 library makes when simplifying
geometry.  Let geometry *A* be an arbitrary collection of points, polylines, and
polygons, and suppose that it is simplified using `S2Builder` with a snap radius
of *r* to yield a new geometry collection *B*.  Then if degeneracies are allowed
when constructing *B*, we can guarantee that

$$H(A, B) \leq r_{\rm e} \, .$$

The quantity $$r_{\rm e}$$ is called the *maximum edge deviation*, which
`S2Builder` ensures is at most 10% greater than the given snap radius *r* (i.e.,
$$r_{\rm e} \leq 1.1 r$$).  **In other words, the maximum possible distance
measurement error due to simplification is only slightly greater than the snap
radius used during simplification.** (Note that the additional 10% is a
guaranteed bound rather than a heuristic; it is only necessary because S2 works
with spherical geometry and would not be needed in a planar version of the
algorithm.)

In fact, under the view that positive and negative degeneracies have
infinitesimal areas, the bound above applies not only to the whole geometries *A*
and *B*, but also to the interiors, boundaries, and exteriors of these
geometries:

$$H({\rm Int}\, A, {\rm Int}\, B) \leq r_{\rm e}$$

$$H(\partial A, \partial B) \leq r_{\rm e}$$

$$H({\rm Ext}\, A, {\rm Ext}\, B) \leq r_{\rm e} \, .$$

These properties are the holy grail of simplification algorithms; it is a
mathematical way of saying that no details larger than the maximum edge devation
$$r_{\rm e}$$ have been lost.

## Representation and Validity

Degenerate geometry is not supported by the legacy classes `S2Polyline` and
`S2Polygon`.  Instead the more modern, efficient, lightweight classes
`S2LaxPolylineShape` and `S2LaxPolygonShape` must be used.  (The `Lax` in the
names refers to the fact that these classes allow degeneracies.)  For example,
while `S2Polygon` requires that all loops have at least three vertices,
`S2LaxPolygonShape` supports loops with two, one, or even zero vertices.  (The
zero-vertex case is called the *full loop* and represents the full sphere.)

The `Lax` classes themselves do not have any validity requirements, but the
geometry they represent must satisfy certain conditions before it can be used in
S2 operations.  Some operations are more restrictive than others; for example,
distance measurement has fewer requirements than `S2BooleanOperation`.  You will
generally want to construct your geometry to satisfy the most restrictive rules,
outlined below.  This can easily be accomplished by using `S2Builder` with an
appropriate output layer type (e.g. `S2LaxPolygonShape`).

### Polylines

Polylines must consist of either a single degenerate edge *AA*, or a sequence of
non-degenerate edges (e.g., *ABC*).  Polylines with both degenerate and
non-degenerate edges are not allowed (e.g., *AABC* or *ABBC*).  On the other
hand, polylines may self-intersect or have duplicate edges.

Note that a point polyline is represented as two vertices (*AA*) whereas a point
polygon loop is represented as one vertex (*A*).  Each of these objects thus
consists of a single degenerate edge AA, corresponding to the general rule that
an *n*-vertex loop has *n* edges, whereas an *n*-vertex polyline has *n-1*
edges.

Here are a few examples of valid and invalid polylines:

 - *AB*: a valid polyline defining a single directed edge.

 - *ABCD*: a valid polyline consisting of three directed edges \{*AB*, *BC*,
   *CD*\}

 - \[\] (the empty vertex sequence): a valid polyline defining no edges.

 - *AA*: a valid polyline defining a single degenerate edge at point *A*.

 - *ABBC*: invalid because degenerate and non-degenerate edges are used in the
   same polyline.

 - *A*: invalid because it defines a vertex but no edges.  (A degenerate edge at
   *A* is represented as *AA*, while an empty polyline is represented as \[\].)

### Polygons

Polygons must consist of a set of oriented loops such that the interior of the
polygon is always on the left.  Polygon edges may intersect only at their
endpoint vertices.  A polygon may contain both an edge and its reverse edge
(i.e., *AB* and *BA*), but duplicate edges are not allowed (*AB* and *AB*).
Loops also must not contain repeated adjacent vertices (e.g., *ABBC*).

Point loops consisting of a single vertex (*A*) are allowed, but such loops
cannot be incident to any other edges.  The *full loop* (i.e., the loop with
zero vertices mentioned above) is allowed only if all the remaining loops (taken
together) define only degeneracies.

Here are a few examples of valid and invalid polygons:

 - {} (the empty set of loops): a valid polygon containing no points (the empty
   polygon).

 - \{*A*, *BC*\}: a valid polygon consisting of one point shell and one sibling
   pair shell.

 - \{*AB*, *BAC*\}: invalid because the edge *BA* occurs twice (note that loop
   *AB* consists of the two edges *AB* and *BA*).

 - \{*ABC*, *A*\}: invalid because the point shell *A* is incident to other
   edges.

 - {full}: a valid polygon containing the entire sphere (the full polygon).

 - {full, *ABCB*\}: a valid polygon consisting of the entire sphere except for
   two sibling pair holes (*AB* and *BC*).

 - {full, *ABC*\}: invalid because the full loop is used together with
   non-degenerate geometry.

 - {full, *ABC*, *ACB*\}: a valid polygon consisting of the entire sphere except
   for three sibling pair holes (*AB*, *BC*, and *CA*\}.  Note that even though
   the loops *ABC* and *ACB* are themselves non-degenerate, together they define
   only degenerate geometry and therefore can be used with the full loop.

### Mixed-Dimensional Geometry {#mixed-dimensional}

The `S2ShapeIndex` class represents an arbitrary collection of points,
polylines, and polygons.  Most S2 operations (e.g. `S2ClosestEdgeQuery`,
`S2RegionCoverer`, `S2ContainsPointQuery`) do not impose any additional validity
requirements on `S2ShapeIndex` beyond requiring that its component shapes are
valid.  The notable exception to this rule is `S2BooleanOperation`, which
requires that **polygon interiors must be disjoint from all other geometry**
(including other polygon interiors).  So for example, polygon interiors must not
intersect any point or polyline geometry.

In order to support degeneracies, `S2BooleanOperation` requires in addition that
**duplicate polygon edges are not allowed.** Note that this rule is only
necessary when an `S2ShapeIndex` contains more than one polygon, since duplicate
edges are already disallowed within individual polygons.  This additional rule
essentially requires that the collection of all polygons in each `S2ShapeIndex`
must satisfy the same validity requirements as a single polygon.  (Also note
that more than one polygon is never needed, because polygons in S2 can have more
than one outer shell and therefore can represent an arbitrary subset of the
sphere.)

Note that although duplicate polygon edges are not allowed, duplicate points and
polyline edges are permitted.  So for example, an `S2ShapeIndex` containing the
points \{*A*, *A*\}, the polyline ABCABC, and the degenerate loop *AB* is a
valid input to `S2BooleanOperation`.  This is true whether the boundary model is
open, semi-open, or closed.  However an `S2ShapeIndex` input containing two
polygons each consisting of the degenerate loop *AB* would be invalid.

Supporting points and polyline edges as multisets makes it much easier to
support time series such as GPS tracks.  For example, a car driving around a
track might easily generate a polyline such as ABCABCABC.  Similarly, this makes
it easier to reconstruct geometry that has been snapped and/or simplified, which
can cause distinct vertices to merge together.  Note that if duplicate values
are not desired, they can easily be removed by choosing the appropriate
`S2Builder` output layer options.

## Semantics of Degeneracies {#semantics}

In general **degeneracies are modeled as infinitesimal regions of the
appropriate dimension.** So for example, a point shell *A* should be thought of
as a tiny closed loop in the vicinity of the vertex *A*, and a degenerate shell
loop *AB* (consisting of the two edges \{*AB*, *BA*\}) should be thought of as a
narrow sliver in the vicinity of edge *AB*.  Similarly, the polyline *AA* should
be thought of as a tiny polyline in the vicinity of vertex *A*.  We say "in the
vicinity of" because under certain boundary models, degeneracies may not contain
their nominal vertices and edges.  For example, under the open boundary model
the point shell *A* does not contain the point *A*.  Instead it should be
thought of as a tiny loop near the point *A*, so tiny that it does not contain
any representable points along its edges or in its interior.

This model implies that **lower-dimensional degeneracies are infinitesimal
subsets of higher-dimensional degeneracies.** For example, the point *A* is
vanishingly small compared to the degenerate polyline *AA*.  This is because the
point *A* is truly zero-dimensional, whereas the polyline *AA* represents a very
small one-dimensional set.  As another example, consider a point *A*, a closed
polyline *AA*, and a closed point shell *A*.  The intersection of the point *A*
and the closed polyline *AA* consists only of the point *A*, because it is an
infinitesimal subset of the polyline, and similarly the closed polyline *AA*
intersected with the closed polygon *A* yields only the polyline *AA*.

It is also worth recalling that in the S2 library, **polyline and polygon edges
do not contain any interior points.**  For example, the closed polygon shell
loop *AB* (consisting of the sibling edge pair \{*AB*, *BA*\}) does not contain
any points except *A* and *B*.  This is because the S2 library models all
geometry as being infinitesimally perturbed such that edges do not pass through
any representable points on the sphere.  This technique is known as *simulation
of simplicity* (Edelsbrunner and Muecke, 1990) and greatly reduces the number of
special cases that need to be handled when implementing robust geometric
algorithms.

## Boundary Models

The S2 library supports modeling polylines and polygons as having open,
semi-open, or closed boundaries.  Here we briefly review these models and
discuss the treatment of degeneracies within them.

Note that the polyline and polygon boundary models are specified independently.
The `PolylineModel` class specifies whether polylines contain their start and/or
end vertex, whereas the `PolygonModel` class specifies whether polygons contain
their vertices and edges.

Also recall that the boundary model is specified for each operation rather than
being a property of the geometry itself.  This allows more flexibility; for
example, if two geometries *A* and *B* intersect under the CLOSED boundary model
but not under the OPEN boundary model, this would mean that the boundaries of
the two geometries intersect but not their interiors.

### Closed Boundary Model

**A closed polyline contains all of its vertices.** For example, the closed
polyline *ABC* intersected with the three points \{*A*, *B*, *C*\} yields \{*A*,
*B*, *C*\}.

Similarly, **a closed polygon contains all of its vertices, edges, and reversed
edges.** So for example, the closed polygon *ABC* intersected with the points
\{*A*, *B*, *C*\} yields \{*A*, *B*, *C*\}, and the closed polygon *ABC*
intersected with the six polylines \{*AB*, *BA*, *BC*, *CB*, *AC*, *CA*\} yields
\{*AB*, *BA*, *BC*, *CB*, *AC*, *CA*\}.

These rules also apply to degeneracies.  So for example, a closed point polyline
*AA* contains its vertex *A*, and a closed sibling pair shell *AB* contains its
vertices \{*A*, *B*\} and its edges \{*AB*, *BA*\}.

Note that **certain degeneracies have no effect on point or edge containment; we
call such degeneracies *redundant*.** In the case of the closed boundary model,
**degenerate holes are redundant** because they do not exclude any points or
edges.  For example, consider a polygon *ABC* containing a point hole *D*.
Since the boundary is closed, the polygon contains the point *D* whether or not
the hole is present.  (The degenerate hole is still considered to exclude an
infinitesimal region from the polygon, it's simply that this region does not
contain any representable points.)

Redundant degeneracies are considered to be part of the polygon's boundary and
therefore still affect distance measurement to the boundary.  In fact, they may
be viewed as a method of expanding the polygon's boundary without affecting its
interior.  They are treated in exactly the same way as non-degenerate regions in
boolean operations; for example, a redundant degeneracy intersected or unioned
with itself yields the original degeneracy, whereas a redundant degeneracy
subtracted from itself yields the empty set.

### Open Boundary Model

**An open polyline contain all of its vertices except the first and last.** For
example, the open polyline *ABC* intersected with the three points \{*A*, *B*,
*C*\} yields \{*B*\}.  The only exception is if the
`polyline_loops_have_boundaries()` option is false and the first and last
vertices of a polyline are the same (e.g. *ABCDA*); in this case the first/last
vertex of the polyline loop is defined to be contained and its boundary is
defined to be empty.  So for example, with this option the open polyline *ABCA*
intersected with \{*A*, *B*, *C*\} yields \{*A*, *B*, *C*\}.

**An open polygon contains none of its vertices, edges, or reversed edges.** For
example, the open polygon *ABC* intersected with the points \{*A*, *B*, *C*\}
yields the empty set, and similarly the open polygon *ABC* intersected with the
six polylines \{*AB*, *BA*, *BC*, *CB*, *AC*, *CA*\} yields the empty set.

These rules also apply to degeneracies.  So for example, an open point polyline
*AA* contains no points at all, and an open sibling pair shell *AB* contains no
points or edges.  In other words, just as degenerate polygon holes are redundant
in the closed boundary model, **degenerate polylines and degenerate polygon
shells are redundant in the open boundary model.**  These types of degeneracies
do not affect point or edge containment, but simply add to the boundary of the
affected geometry.  The geometry contains exactly the same set of points and
edges whether these degeneracies are present or not, however since the
degeneracies change the boundary they can still affect distance measurement.

Also note that **geometric objects have exactly the same boundary under all
three boundary models.** So for example, the boundary of an open degenerate
polyline *AA* is the point *A*, even though the polyline is defined not to
contain that point.

### Semi-Open Boundary Model

The purpose of the semi-open boundary model is to allow geometry to **cover a
region without gaps or overlaps.**

In order to achieve this, **a semi-open polyline contains all of its vertices
except the last.** For example, the semi-open polyline *ABC* intersected with
the points \{*A*, *B*, *C*\} yields \{*A*, *B*\}.  The major advantage of this
model is that polylines can be split into pieces without affecting vertex
containment.  For example, the single polyline *ABCDEFG* contains exactly the
same set of vertices as the three polylines \{*ABC*, *CDE*, *EFG*\}, and
furthermore each vertex is contained by exactly one of these polylines.

Similarly, semi-open polygon point containment is defined such that **if several
semi-open polygons tile the region around a vertex, then exactly one of those
polygons contains that vertex.** This implies that a triangle *ABC* may or may
not contain each of its vertices \{*A*, *B*, *C*\}.  However it ensures the
important property that **if a set of polygons tiles the entire S2 sphere, then
every point on the sphere is contained by exactly one polygon.** This property
can be very useful for certain algorithms.

**Semi-open polygons contain all of their edges, but none of their reversed
edges.** For example, consider two abutting triangles *ABC* and *CBD*.  The
edge *BC* is contained by *ABC* but not *CBD*, and the edge *CB* is contained
by *CBD* but not *ABC*.  This property ensures **when a set of polygons tiles
the sphere, each polygon edge is contained by exactly one polygon.**

This rule also ensures that if a polyline is intersected with a set of polygons
that tile the sphere, the resulting polylines can be unioned to obtain the
original polyline without gaps or overlaps.  For example, let *ABCD* be a small
square polygon and consider its complement *DCBA*.  These two polygons do not
intersect and their union is the entire S2 sphere.  Now consider the zig-zag
polyline *ABDC*: its intersection with the polygon *ABCD* is the polyline *ABCD*
(, whereas its intersection with the polygon *DCBA* is *DC*.  Together these two
polylines \{*ABD*, *DC*\} can be unioned to obtain the original polyline *ABDC*
without gaps or overlaps.

**All types of degeneracies are redundant in the semi-open model.** In other
words, degeneracies never affect point or edge containment but simply add to the
boundary of the geometry.  For example, a semi-open point polyline *AA* does not
contain the point *A* (just like the open model), semi-open degenerate polygon
shells do not contain any points or edges (just like the open model), and
semi-open degenerate polygon holes do not exclude any points or edges (just like
the closed model).  All of these degeneracies simply add to the boundary of the
geometry without affecting its interior or exterior.  The semi-open model is
purest of the three models in that shells and holes are treated symmetrically.

## Boolean Operations

Many properties of boolean operations involving degeneracies follow naturally
from the definitions above.  Here we summarize the additional principles that
define the results of such operations (intersection, union, difference, and
symmetric difference).

**Operations on geometry of a given dimension yields geometry of the same
dimension.** So for example, a polyline *AB* intersected with a polyline *BC*
under the closed model yields the degenerate polyline *BB* rather than a point.
Similarly a polygon *ABC* intersected with an abutting polygon *CBD* under the
closed model yields the degenerate polygon shell *BC*, not a polyline.  And
symmetrically, a polygon hole *ACB* unioned with an abutting polygon hole *BCD*
under the open model yields a degenerate polygon hole *BC* (which has no
alternative representation).  This behavior is exactly what allows S2 to avoid
the serious problems of [dimension reduction](#dimension-reduction).

**In all boundary models, subtracting an object from itself yields the empty
set.** This is true even for redundant degeneracies that contain no points or
edges, such as subtracting a point polyline *AA* from itself in the open
boundary model.

**In all boundary models, intersecting or unioning an object with itself yields
the original object(s).** The reason for the plural is that points and polyline
edges are treated as multisets.  For example, if a point *A* is intersected with
another point *A*, the result is two points \{*A*, *A*\}.  (Note that such
duplicates can easily be filtered away by choosing appropriate options in the
`S2Builder` output layer.)  Again, this is true even for redundant degeneracies
that contain no points or edges, so for example a point shell *A* intersected
with itself in the open boundary model yields the same point shell *A*, and a
point polyline *AA* intesected with itself in the open boundary model yields two
point polylines \{*AA*, *AA*\} (neither of which contains the point *A*).

**Intersecting geometry of different dimensions yields only the
lower-dimensional geometry.** This is consistent with the principles described
under [Semantics of Degeneracies](#semantics) above.  For example, intersecting
a point *A* with a closed polygon *ABC* yields only the point *A*, not a point
polygon shell *A* as well.  Of course, if the input geometries do not intersect
then the result is empty.  For example intersecting a point *A* with an open
polygon *ABC* yields nothing at all.

**Unioning geometry of different dimensions yields only the higher-dimensional
geometry.** So for example, unioning a point *A* with a closed polyline *AA*
yields only the polyline *AA*.  Of course if the input geometries do not
intersect then the result is simply the collection of both inputs.  For example
the union of a point *A* with an open or semi-open polyline *AA* is the
collection of both objects \{*A*, *AA*\} since they do not intersect.

**Unioning points with other points, or polylines with other polylines, yields a
multiset.** For example, the union of two identical points *A* is the multiset
\{*A*, *A*\}, and the union of two identical polylines *ABC* is the multiset
\{*ABC*, *ABC*\}.  As [noted earlier](#mixed-dimensional) this behavior makes it
possible to reconstruct time series such as GPS tracks accurately.  Unwanted
duplicates can be filtered by choosing appropriate options in the `S2Builder`
output layer.

### Semi-Open Semantics

The semi-open boundary model has a few properties that deserve further
explanation.

**In the semi-open model, complementary degeneracies do not intersect.**  So for
example, a point shell *A* does not intersect a point hole *A*, and a sibling
pair shell *AB* does not intersect a sibling pair hole *AB*.

**In the semi-open model, a point polyline *AA* or point shell *A* intersects
another different geometry if and only if that other geometry contains the point
*A*.** So for example, a point shell *A* intersects a semi-open polygon *ABC* if
and only if *ABC* contains its vertex *A* under the usual semi-open rules.  In
other words, positive point degeneracies are treated exactly like points (except
when a degeneracy is intersected with itself---see above).

The purpose of this rule is to ensure that the defining property of the
semi-open model (i.e. the ability to cover a region without overlaps or gaps)
applies to point degeneracies as well as points.  For example, if a degenerate
polyline *AA* or degenerate shell *A* is intersected with a set of semi-open
polygons that tile the region around a shared vertex *A*, then exactly one of
those polygons contains that degeneracy.  Similarly, if a semi-open polyline
*ABC* is expressed as the disjoint union of sub-polylines \{*AB*, *BC*\} then
exactly one of those polylines contains the degenerate polyline *BB* (in
particular, *BC* contains *BB* while *AB* does not).

As noted above, the only way that point degeneracies are handled differently
 than points is with respect to self-intersection.  Specifically, in the
 semi-open model:

 - A point *A* intersected with a polyline *AA* or point shell *A* is empty
   (since the latter objects do not contain any points).
 - A polyline *AA* intersected with a point shell *A* is empty (since the point
   shell contains no points or polylines).
 - And yet, a point *A*, polyline *AA*, or point shell *A* intersected with
   itself yields the original object.

These statements are not contradictory; there is in fact a valid representation
of these objects as infinitesimal regions that yields this behavior.

The complementary property is that **the union of a point hole *A* with another
different geometry eliminates the hole if and only if that other geometry
contains the point *A*.**  So for example, a point hole *A* unioned with the
semi-open polygon *ABC* eliminates the hole from the result if and only if *ABC*
contains its vertex *A*.

A similar property for edge degeneracies is that **sibling pair shells do not
intersect abutting polygons.** So for example, a shell *BC* does not intersect
the polygon *ABC* or *CBD*.  The best way to think about this is that the
intersection region is only one-half of the sibling pair degeneracy region and
therefore cannot be represented.  The complementary property is that **the union
of a sibling pair hole with an abutting polygon eliminates the hole**.  So for
example, the union of the polygon *ABC* with the full sphere except for the hole
*BC* is the full sphere.

### Limitations

It is important to realize that in some cases the true mathematical result of an
operation cannot be represented.  This is true even when no degeneracies are
present.

For example, consider a triangle *ABC* that contains a nested triangle *DEF*,
and suppose that *DEF* is subtracted from *ABC* in the CLOSED boundary model.
Ideally the result would contain the edges \{*AB*, *BC*, *CA*\} but not the
edges \{*DE*, *EF*, *FD*\} (since the latter three edges belong to *DEF*, which
has been subtracted).  However it is not possible to represent this in the
CLOSED model since all of these edges are on the boundary of the result.

Similarly, results involving degeneracies sometimes cannot be exactly
represented.  In such cases the preference is always to correctly represent the
**interior** of the result (including infinitesimal interiors) even when the
boundary of the result cannot be perfectly represented.

For example, consider subtracting a point shell *D* from the interior of a
closed polygon *ABC*.  The result is defined to be the polygon \{*ABC*, *D*\}
where *D* is now a point hole.  This may seem strange since (1) the result
still contains the point *D* because the boundary is closed, and (2) there is
no obvious reason to include *D* in the result because point holes exclude no
points under the closed model (i.e., they are redundant).  And yet this is the
**only** correct result under the S2 degeneracy model, because the original
point shell *D* is considered to include an infinitesimal two-dimensional area
near the point *D* in addition to *D* itself.  The only way to exclude this
infinitesimal two-dimensional area from the result is to include *D* as a point
hole.

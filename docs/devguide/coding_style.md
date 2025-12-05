---
title: S2 C++ Style Guide
---



## Purpose

This guide documents the C++ coding conventions used by the S2 Geometry Library.
It does not attempt to justify why these conventions are used, and indeed in
some cases they could be improved upon.

In general, the S2 C++ coding style is built on the [Google C++ Style
Guide](https://google.github.io/styleguide/cppguide.html).  However it's
important to remember that the S2 library has a long history, and many parts of
the library were developed before modern C++ features even existed.  The S2
library has many clients both within and outside of Google (through the open
source release), which can make it very difficult to adopt new conventions.

**The overriding objective when adding new code should be to maintain
consistency within the library.** In cases where the S2 conventions conflict
with the current Google style guide, the existing S2 conventions should take
precedence.  Changes in coding style should only be made when it is practical to
adopt them across the entire library, including the important step of updating
all affected clients so that any deprecated classes or methods can be removed.

The first step when adding new classes or methods should be to find similar
classes or methods that already exist in the library, and then follow their
conventions as closely as possible.

## External Dependencies

The S2 library is built on many platforms, and has an open source version with
many external clients.  New dependencies should be avoided.  In particular, note
that S2 does not currently depend on the following:

 - GMock (except for a few files that not currently part of the open-source
   release).  Instead please limit yourself to the more basic features of GUnit.

 - Protocol buffers.

 - Any Google maps or geo library.

 - Many features of Abseil, including `StatusOr` (see "Error Handling" below).

## Filenames

The S2 library has a flat directory structure.  Generally each header file
contains either a single class (plus any associated helper classes) or a
collection of functions in a namespace.

Filenames should be in lowercase with words separated by underscores.  Note that
"s2", "s1", etc. are not considered to be words.  Examples:

 - `s2builder.h`
 - `s2cell_id.h`
 - `encoded_s2shape_index.h`

For historical reasons, `latlng` is treated as one word even though the
corresponding classes are camel-cased.  For example, `S2LatLngRect` is defined
in `s2latlng_rect.h`.

Files that define functions or classes within the `s2shapeutil` and
`s2builderutil` namespaces include that namespace as a prefix:

 - `s2shapeutil_conversion.h`
 - `s2shapeutil_count_edges.h`
 - `s2builderutil_lax_polyon_layer.h`
 - `s2builderutil_snap_functions.h`

Classes or functions that are needed in multiple places within the S2 library
but are not part of the public API can be put in header files with the
`_internal.h` suffix.  Examples:

 - `s2coords_internal.h`
 - `s2predicates_internal.h`

This can help remove clutter from the public API and also facilitate testing.
Note that the internal content can be implemented in the same `.cc` file as the
public part of the API.

"Util" header files are strongly discouraged.  Instead try to break your API
into meaningful components.  The only remaining "util"-style header file in the
S2 library is `s2pointutil.h` which contains about ten functions for
manipulating points.

## Documentation

### Comment Style

Please put effort into your comments; they are the most important part of your
code.

Comments should be full sentences, starting with a capital letter and ending
with a period.  This rule applies even to partial-line comments when reasonable:

```c++
  enum IndexStatus {
    STALE,     // There are pending updates.
    UPDATING,  // Updates are currently being applied.
    FRESH,     // There are no pending updates.
  };
```

For bonus points, put two spaces after your periods: contrary to popular belief,
this is not some archaic convention from the days of typewriters, but rather an
effort to mimic professional typography (take a close look at a book sometime).

### Header Files

Public classes and methods should be thoroughly documented in the header file.
Important topics to cover include:

 - Description of the problem being solved
 - Code samples showing how to use the API
 - Limitations and known problems

Requirements of classes and methods may be documented via a "REQUIRES" comment.
For example:

```c++
// Encodes an unsigned integer in little-endian format using "length" bytes.
// (The client must ensure that the encoder's buffer is large enough.)
//
// REQUIRES: T is an unsigned integer type.
// REQUIRES: 2 <= sizeof(T) <= 8
// REQUIRES: 0 <= length <= sizeof(T)
// REQUIRES: value < 256 ** length
// REQUIRES: encoder->avail() >= length
template <class T>
void EncodeUintWithLength(T value, int length, Encoder* encoder);
```

The following callouts may also be useful:

 - `ENSURES:` properties guaranteed by the class/method.
 - `DEFAULT:` documents the default value of an option.
 - `CAVEAT:` an unexpected property or potential pitfall.
 - `NOTE:` to draw special attention to a comment.

### Implementations

Implementation details should be documented where the corresponding types or
methods are defined, not where they are declared.  This means that types and
member variables declared in the header file should be documented there, whereas
private non-inline methods should be documented in the .cc file.  This
convention makes it easier to find documentation when looking at code, and also
makes it more likely that the documentation will be updated when implementations
are changed.

Here is an example of documenting a private method in a `.cc` file:

```c++
// Outputs the current buffered path (which is assumed to be a loop), and
// resets the state to prepare for buffering a new loop.
void S2BufferOperation::OutputPath() {
  op_.AddLoop(path_);
  path_.clear();  // Does not change capacity.
  have_input_start_ = false;
  have_offset_start_ = false;
}
```

### External Documentation

In addition to header file comments, significant new classes should have
documentation added in the "g3doc" directory.  Typically this will expand
significantly on what is in the header file, giving more background information
and examples.  Documentation added here will be published on the `s2geometry.io`
website.

Ideally documentation should be added and/or updated in the same changelist that
adds the relevant code.

## Namespaces

The S2 library as a whole is not defined within a namespace.  (Namespaces did
not even exist when it was first designed.)  Instead most of the global symbols
include `S2` in their names (except for a few classes that use `S1`, `R2`,
etc. instead.)

Nevertheless, some parts of the library do use namespaces.  In general the names
of namespaces should be in lowercase with words separated by underscores
(e.g. `s2polyline_alignment`).  Note that the "util" suffix is not considered a
word, so we have:

 - `namespace s2shapeutil { ... }`
 - `namespace s2builderutil { ... }`

There is also a namespace called `S2` that contains many of the low-level
functions defined by the library (e.g. `S2::TurnAngle`, `S2::RobustCrossProd`,
etc.)  Note that while class names often include `S2` as a prefix rather
than being placed in a namespace, this convention is not typically used for
global functions or constants and so these entities should be put in a namespace
instead.

In general S2 does not use nested namespaces, except for internal content which
may be placed in a nested namespace called `internal`.  This technique can be
used alone to protect internal content (e.g. the file `s2testing.h` defines some
content in the `S2::internal` namespace) and/or with an `_internal.h` header
file.

## Using Declarations

`using` declarations are used extensively in `.cc` files in order to reduce code
verbosity.  This includes types and functions from the `std` and `absl`
namespaces as well as aliases for nested types.  For example,
`s2buffer_operation_test.cc` has the following declarations:

```c++
using absl::make_unique;
using std::max;
using std::string;
using std::unique_ptr;
using std::vector;

using EndCapStyle = S2BufferOperation::EndCapStyle;
using PolylineSide = S2BufferOperation::PolylineSide;
```

Note that this is allowable only in `.cc` files, not header files, and that
`using namespace` declarations are prohibited.

## Classes

Class names in S2 should be nouns with the first letter of each word
capitalized.  Examples:

 - `S2Shape`
 - `S2LatLngRect`
 - `S1Angle`

Short, descriptive names are highly preferred.  It is reasonable to put as much
thought into the names of your classes, methods, variables, etc. as you put into
the code itself.  Good names are essential to making APIs easier to understand.

The names of subtypes should be obviously related to their parent class.  For
example, `S2LaxPolygonShape` is a subtype of `S2Shape`, and
`s2builderutil::S2CellIdSnapFunction` is a subtype of `S2Builder::SnapFunction`.

### Index, Query, Operation

Classes in S2 should be designed so as to separate the data representation from
the actual operations as much as possible.A class that contains many geometric
objects for the purpose of performing operations on them should preferentially
be called an `Index`.  Such classes should define only the minimum number of
functions necessary to access the data in the index.  Examples:

 - `S2ShapeIndex`
 - `S2CellIndex`
 - `S2PointIndex`

Rather than having a single class with a wide API, prefer to have several
classes or standalone functions with narrow APIs.  For example, rather than
having a method in `S2ShapeIndex` to count the total number of edges, instead
there is a separate function to do this defined in `s2shapeutil_count_edges.h`.

A class that does something useful with an index should preferentially be called
a `Query`, unless it can be viewed as function that maps an index to another
similar index in which case it may be called an `Operation`.  Examples:

 - `S2ClosestPointQuery`
 - `S2ConvexHullQuery`
 - `S2CrossingEdgeQuery`
 - `S2BooleanOperation`
 - `S2BufferOperation`

`Index` classes should be thread-safe but `Query` and `Operation` classes
generally are not.  Instead the expectation is that each thread should create
its own `Query` objects (typically on the stack) as needed.  The fact that
`Query` classes are not thread-safe allows them to maintain internal state while
executing the query and/or to cache information between queries.

`Query` and `Operation` classes typically follow one of two design patterns.
One is the builder pattern, where data can be added incrementally using several
method calls, and then a different method is called to compute the result.
Examples following this pattern include `S2WindingOperation`,
`S2BufferOperation`, `S2ConvexHullQuery`, etc.

The other common pattern is a pure function pattern where the arguments are
supplied and the result is returned during a single function call.  Examples of
this pattern include `S2ClosestEdgeQuery`, `S2ClosestPointQuery`,
`S2ContainsPointQuery`, etc.  Note that even though this pattern in theory
allows the query method to be stateless and therefore thread-safe, S2 classes
avoid this in favor of rather allowing query objects to maintain whatever state
they may need in order to make queries as efficient as possible.  This means
that **query methods should not be const** even if the current implementation of
the method would allow this, since this might prevent future optimizations.

When query methods need to return several values of different types, the query
method should return a result object.  For example, `S2ClosestEdgeQuery` returns
a `Result` object with methods such as `distance()`, `shape_id()`, and
`edge_id()`.

This technique is strongly preferred over the alternative technique of having
the query object store the latest result in its own state and having methods to
retrieve the result field-by-field.  For example, a Result object makes it much
easier to combine the results of several queries:

```
  auto a_result = query.GetResult(a);
  auto b_result = query.GetResult(b);
  ... do something with both results ...
```

Note that as long as the optimizer does its job correctly, returning a result by
value does not involve any actual copying.  Instead the caller passes a hidden
pointer and the result is constructed directly on the caller's stack.

### Options

Classes with options should have a nested `Options` class.  For
example, `S2Builder` has an `S2Builder::Options` class.  Options
classes generally consist of a set of fields with "get" and "set" methods.  For
example:

```c++
class S2RegionCoverer {
 public:
  class Options {
   public:
    Options();

    int max_cells() const;
    void set_max_cells(int max_cells);
    ...
  };
  explicit S2RegionCoverer(const Options& options);
  const Options& options() const;
  ...
};
```

The set methods generally return `void`, but it is also acceptable to return a
reference to the Options object in order to allow chained initialization:

```
  S2Foo foo{S2Foo::Options().set_bar(x).set_baz(y)};
```

Options are usually read-only once the relevant object has been constructed.
Classes that allow options to be modified after construction do so via a
`mutable_options()` method:

```c++
  Options* mutable_options();
```

Options classes are copyable and movable.

## Methods

Accessor method names are nouns, written in lowercase with words separated by
underscores.  Examples:

 - `S1Angle::degrees()`
 - `S2Builder::options()`

Non-accessor method names should be verbs, written in uppercase with the first
letter of each word capitalized.  Examples:

 - `S2Builder::StartLayer(...)`
 - `S2Builder::Build(...)`
 - `S2ShapeIndex::Minimize()`

Methods that return something but that don't have an obvious verb should
typically start with `Get`.  For example:

 - `S2ShapeIndex::GetShape(int id)`
 - `S2::GetPointOnLine(...)`

If the method may return an empty set of results, the verb `Find` can also be
used.  For example:

 - `S2ClosestEdgeQuery::FindClosestEdges(...)`

There are exceptions for well-known geometric operations that yield objects of
the same type.  For example, `S2Cap::Union` computes the union of two `S2Cap`
objects.  This naming convention is generally reserved for low-level geometric
objects such as `S1Interval`, `S2LatLngRect`, `S2Cap`, etc.

Methods that return boolean values (predicates) should be active verbs
(e.g. `S2LatLngRect::Contains()`) or start with the word "is"
(e.g. `S2Shape::is_empty()`, `S2::IsUnitLength()`).  (Note that this convention
differs from `std::vector<T>::empty()` for example.)

### Inline Methods

Except for accessor methods, inline methods should be declared and defined
separately (just like non-inline methods).  Note that only the definition needs
to be declared inline, not the declaration.  Inline definitions should be placed
at the bottom of the corresponding header file after a comment saying
"Implementation details follow".  For example:

```c++
class S2PaddedCell {
 public:
  // Construct an S2PaddedCell for the given cell id and padding.
  S2PaddedCell(S2CellId id, double padding);

  S2CellId id() const { return id_; }
  double padding() const { return padding_; }
  int level() const { return level_; }

  // Return the bound for this cell (including padding).
  const R2Rect& bound() const { return bound_; }

  // Return the (i,j) coordinates for the child cell at the given traversal
  // position.  The traversal position corresponds to the order in which child
  // cells are visited by the Hilbert curve.
  void GetChildIJ(int pos, int* i, int* j) const;

  ...
};


//////////////////   Implementation details follow   ////////////////////


inline void S2PaddedCell::GetChildIJ(int pos, int* i, int* j) const {
  int ij = S2::internal::kPosToIJ[orientation_][pos];
  *i = ij >> 1;
  *j = ij & 1;
}
```

### Constructors, Init, Factory Methods

Virtually all classes should have constructors.  If the constructor parameters
may be ambiguous, static factory methods should be used instead.  For example,
`S2LatLng` has the following factory methods:

```c++
  static constexpr S2LatLng FromRadians(double lat_radians, double lng_radians);
  static constexpr S2LatLng FromDegrees(double lat_degrees, double lng_degrees);
  static constexpr S2LatLng FromE5(int32_t lat_e5, int32_t lng_e5);
  static constexpr S2LatLng FromE6(int32_t lat_e6, int32_t lng_e6);
  static constexpr S2LatLng FromE7(int32_t lat_e7, int32_t lnt_e7);
```

Factory method names should begin with "From", except for simple and frequently
used classes where brevity takes precedence.  For example, `S1Angle` has the
factory methods `S1Angle::Radians(double radians)` and `S1Angle::Degrees(double
degrees)`.

Factory methods should also be used when the parameter order may be ambiguous.
Examples:

```c++
  static R2Rect FromCenterSize(const R2Point& center, const R2Point& size);
  static S2CellId FromFacePosLevel(int face, uint64 pos, int level);
```

Complex classes that allocate memory should also have a default constructor and
an `Init` method.  This allows clients to declare instances of these classes as
member variables and defer their initialization as necessary.  This allows
clients to avoid allocating objects on the heap and accessing them through
`unique_ptr<>` just because the arguments needed to initialize them are not
readily available in their own constructor.  This not only reduces memory
allocation but also increases locality of reference.

So for example, `S2Builder` has the following:

```c++
class S2Builder {
  ...
  // Default constructor; requires Init() to be called.
  S2Builder();

  // Convenience constructor that calls Init().
  explicit S2Builder(const Options& options);

  // Initializes an S2Builder with the given options.
  void Init(const Options& options);

  const Options& options() const { return options_; }
  ...
};
```

Factory methods that return `unique_ptr<T>` should be avoided, since this
forces clients to allocate memory and also reduces locality of reference.

## Types

### Templates

Templated types and functions are discouraged in the public parts of the API.
Templated code is necessarily harder to document and understand, and the
generalization that it offers is frequently overrated and underused.

Templates may be used internally when this substantially reduces code
duplication or otherwise makes the library easier to maintain.  For example, the
S2 classes `S2ClosestEdgeQuery` and `S2FurtherEdgeQuery` are each wrappers
around a templated internal class called `S2ClosestEdgeQueryBase`.

### Distances and Angles

Distances and angles should be represented using the classses `S1Angle` or
`S1ChordAngle`.  S2 APIs should always use these classes in preference to
accepting distances in radians, degrees, meters, etc. since it avoids the need
to specify units.  (Note that the S2 library operates on the unit sphere, not
the Earth's surface, and therefore the natural unit of distance is the radian.)

Whenever possible, classes should support both `S1Angle` and `S1ChordAngle`.
`S1ChordAngle` is frequently used internally, since (1) it is much faster to
compute and (2) it supports the exact distance predicates that are needed to
implement algorithms robustly (see `s2predicates.h`).  Methods that accept an
`S1Angle` often convert it to `S1ChordAngle` internally.

Note that `S1ChordAngle` is only capable of representing distances up to 180
degrees, and that its accuracy declines for distances approaching 180 degrees.
Distances which may be greater than 180 degrees, such as polyline lengths or
polygon perimeters, should always be represented using `S1Angle`.

### Parameters

Small types (16 bytes or less) should generally be passed by value.  This
includes the following commonly used S2 types:

 - `S2CellId`
 - `S2LatLng`
 - `S1Angle`
 - `S1ChordAngle`

Types larger than 16 bytes should be passed as follows:

 - If ownership is being transferred, use `unique_ptr<T>`.

 - If the callee keeps a pointer to the object that outlives the function call,
   use `const T*`.  For example, `S2ClosestEdgeQuery` keeps a pointer to its
   `S2ShapeIndex` constructor argument, so it accepts that argument by const
   pointer:

   ```c++
   explicit S2ClosestEdgeQuery(const S2ShapeIndex* index,
                               const Options& options = Options());
   ```

 - Otherwise, the parameter should be passed by const reference.  Note in
   particular that `S2Point` is passed by const reference.  Const pointers
   should not be used unless the callee keeps a pointer as outlined above.

 - The only exception to the rule above is when the callee needs to make a copy
   of the argument **and*** the given type has an efficient move operator (i.e.,
   one that transfers ownership of data rather than simply copying it).  In this
   case the parameter may be passed by value rather than by const reference.
   For example:

   ```c++
   class S2CellUnion final : public S2Region {
   ...

   // Constructs a cell union with the given S2CellIds, then calls Normalize()
   // to sort them, remove duplicates, and merge cells when possible.
   //
   // The argument is passed by value, so if you are passing a named variable
   // and have no further use for it, consider using std::move().
   //
   // A cell union containing a single S2CellId may be constructed like this:
   //
   //     S2CellUnion example({cell_id});
   explicit S2CellUnion(std::vector<S2CellId> cell_ids);
   ...
   }

   inline S2CellUnion::S2CellUnion(std::vector<S2CellId> cell_ids)
   : cell_ids_(std::move(cell_ids)) {
   Normalize();
   }
   ```

   Note that this guidance is considerably more strict than advocated by the
   current Google style guide.  Do not pass parameters larger than 16 bytes
   unless (1) the type has an efficient `std::move` operator and (2) the callee
   makes a copy of the argument using `std::move`.

#### Const Value Parameters

Parameters passed by value should never be "const" in method declarations.  They
should only be declared const (if desired) in method **definitions**.  For
example:

```c++
class S2CellId {
  ...
  // Return true if the given cell is contained within this one.
  bool contains(S2CellId other) const;
  ...
};


//////////////////   Implementation details follow   ////////////////////


inline bool S2CellId::contains(const S2CellId other) const {
  DCHECK(is_valid());
  DCHECK(other.is_valid());
  return other >= range_min() && other <= range_max();
}
```

Note that most S2 code does not bother declaring local variables and parameters
passed by value as "const" even when those values are not modified, since the
additional documentation benefit often does not outweigh the additional code
clutter.

### Geometry Output

Classes that construct geometry and return it to the client should do so via an
`S2Builder::Layer` parameter.  For example, the `S2BufferOperation` constructor
looks like this:

```c++
  explicit S2BufferOperation(std::unique_ptr<S2Builder::Layer> result_layer,
                             const Options& options = Options{});
```

When the client eventually calls the `S2BufferOperation::Build()` method, the
output geometry is sent to the given `result_layer`.  This allows the client to
have complete control over how the output geometry is represented.  There are
various possible layer types that build different representations such as
`S2Polyon` or `S2LaxPolygonShape`.

### Obsolete Types

Clients are discouraged from using certain S2 classes for new code.  Classes
that have officially been deprecated are marked as such using the
`ABSL_DEPRECATED` macro.  Currently the only such class is `S2PolygonBuilder`
which has been superseded by the far more robust and flexible `S2Builder`.

However it is important to note that the classes `S2Polygon` and `S2Polyline`
are also obsolete and should ideally be deprecated at some point in the future.
New code should prefer to use their `S2Shape`-based replacements,
`S2LaxPolygonShape` and `S2LaxPolylineShape`.  These new classes can be
constructed and decoded far more efficiently than `S2Polygon` and `S2Polyline`,
and are also capable of representing degenerate geometry if desired (which is
what the `Lax` in their names refers to).

## Error Handling

In general S2 code should try to avoid generating runtime errors.

### Parameter Range Checking

If a parameter has a limited valid range, the preferred technique in S2 is to
`DCHECK` that the parameter is valid and then clamp it to the
allowed range.  For example:

```c++
  class S2BufferOperation::Options {
    ...
    // Specifies the allowable error when buffering, expressed as a fraction
    // of buffer_radius().
    //
    // REQUIRES: error_fraction() >= kMinErrorFraction
    // REQUIRES: error_fraction() <= 1.0
    //
    // DEFAULT: 0.01  (i.e., maximum error of 1%)
    static constexpr double kMinErrorFraction = 1e-6;
    double error_fraction() const;
    void set_error_fraction(double error_fraction);
    ...
  };

void S2BufferOperation::Options::set_error_fraction(double error_fraction) {
  DCHECK_GE(error_fraction, kMinErrorFraction);
  DCHECK_LE(error_fraction, 1.0);
  error_fraction_ = max(kMinErrorFraction, min(1.0, error_fraction));
}
```

### Invariants

Internal invariants should be enforced via `DCHECK`.  This is an important
aspect of both documentation and testing.  For example, the following is a
private method that should only be called when certain conditions are true:

```c++
void S2BufferOperation::BufferEdgeAndVertex(const S2Point& a, const S2Point& b,
                                            const S2Point& c) {
  DCHECK_NE(a, b);
  DCHECK_NE(b, c);
  DCHECK_NE(buffer_sign_, 0);
  if (!tracker_.ok()) return;
  ...
}
```

Programming errors within the S2 library should not be converted to runtime
errors if at all possible.


### S2Error

For historical reasons, S2 defines its own error class `S2Error` consisting of
an error code and a status message.  This class is used exclusively within the
library in preference to other error representations such as `absl::Status`.
While it may be desirable to eventually convert the library to use a different
error representation, the overriding principle is to maintain consistency within
the library so that users (including external users of the open source version)
do not need to deal with two different error representations at once.

Here is a partial definition of the `S2Error` class:

```c++
class S2Error {
 public:
  enum Code {
    OK = 0,                  // No error.

    ////////////////////////////////////////////////////////////////////
    // Generic errors, not specific to geometric objects:

    UNKNOWN = 1000,              // Unknown error.
    UNIMPLEMENTED = 1001,        // Operation is not implemented.
    OUT_OF_RANGE = 1002,         // Argument is out of range.
    INVALID_ARGUMENT = 1003,     // Invalid argument (other than a range error).
    FAILED_PRECONDITION = 1004,  // Object is not in the required state.
    INTERNAL = 1005,             // An internal invariant has failed.
    DATA_LOSS = 1006,            // Data loss or corruption.
    RESOURCE_EXHAUSTED = 1007,   // A resource has been exhausted.
    CANCELLED = 1008,            // Operation was cancelled.

    ...
  };
  S2Error() : code_(OK), text_() {}

  // Convenience constructor that calls Init().
  template <typename... Args>
  S2Error(Code code, const absl::FormatSpec<Args...>& format,
          const Args&... args);
  }

  // Set the error to the given code and printf-style message.  Note that you
  // can prepend text to an existing error by calling Init() more than once:
  //
  //   error->Init(error->code(), "Loop %d: %s", j, error->text().c_str());
  template <typename... Args>
  void Init(Code code, const absl::FormatSpec<Args...>& format,
            const Args&... args);

  bool ok() const { return code_ == OK; }
  Code code() const { return code_; }
  std::string text() const { return text_; }

  // Clear the error to contain the OK code and no error message.
  void Clear();

 private:
  Code code_;
  std::string text_;
};
```

Errors are typically initialized using `S2Error::Init()`, although they can also
be constructed and copied.  For example:

```c++
  error->Init(S2Error::RESOURCE_EXHAUSTED,
              "Memory limit exceeded (tracked usage %d bytes, limit %d bytes)",
              usage_bytes_, limit_bytes_);
```

`S2Error` is typically passed as a pointer parameter and is required to be
non-`nullptr`.  (Note that the current Google style recommendation is to pass such
parameters by reference, which precludes `nullptr` arguments at the expense of
giving up the visible documentation hint that the argument may be modified.)
Public methods typically call `error->Clear()` to remove any existing error,
whereas private methods either set the error or leave its state unchanged.

Methods that would otherwise return `void` often return `error->ok()` instead in
order to simplify client code.  For example:

```c++
class S2BooleanOperation {
  ...
  // Executes the given operation.  Returns true on success, and otherwise
  // sets "error" appropriately.  (This class does not generate any errors
  // itself, but the S2Builder::Layer might.)
  bool Build(const S2ShapeIndex& a, const S2ShapeIndex& b, S2Error* error);
}
```

Client code might then look like this:

```c++
  bool DoSomething(..., S2Error* error) {
    ...
    if (!op.Build(a, b, error)) return false;
    ...
  }
```

If a function has a void or other non-boolean result, the caller must check
`error->ok()` explicitly:

```c++
  void DoSomething(S2Builder::Layer* layer, const S2Builder::Graph& graph,
                   S2Error* error) {
    layer->Build(graph, error);
    if (!error->ok()) return;
    ...
  }
```

## Testing

Tests should attempt to achieve excellent coverage with a minimum of code.  This
will usually mean writing helper functions and/or test fixtures in order to
avoid code duplication.  Do not be afraid to write test code with loops,
multiple levels of helpers, etc., if this will result in better testing with
less code duplication.

S2 tests should always be run **in a debug build** (not the default "fast
build") since this will give full line number information on all your stack
traces.  In other words, build S2 with `-c dbg` when you are testing new code.

Non-optimized builds also enable various internal consistency checks.  This
includes not only the `DCHECK` macros but also the `--s2debug` flag which
enables certain internal validity checking.  (This flag can also be enabled by
hand if you are suspect a problem with geometry that is too large or slow to
process using a debug build.)

Tests should focus on the functionality provided by your class, not the
functionality of the classes that your class depends on.  For example, if your
class uses `S2Builder` then you should not write tests whose purpose is simply
to verify the guarantees that `S2Builder` makes (which are thoroughly tested
elsewhere).

### Allow For Errors

Tests should generally allow for small amounts of error in their results rather
than expecting an exact match.  This helps to make tests less fragile, avoiding
the need to update them whenever small changes to the underlying algorithms are
made.

Useful Gunit helper functions include `EXPECT_DOUBLE_EQ` and `EXPECT_NEAR`.
It can also be useful to snap results to a fixed resolution using
`s2builderutil::IntLatLngSnapFunction()`.

"Golden" files (which compare the results of a series of operations to a
"golden" result deemed to be correct) are strictly prohibited.` Such tests are
inherently fragile and will break whenever the smallest details of S2 algorithms
are changed.

### S2Testing

The header file `s2testing.h` provides many utility functions and classes that
are useful when writing tests.  In addition the distance-related functions
mentioned above, there is extensive support for generating random geometry of
various kinds.

Support code for writing tests can also be found in the following files:

 - `thread_testing.h`
 - `s2shapeutil_testing.h`
 - `s2closest_edge_query_testing.h`

### S2Earth

S2 tests are not allowed to depend on `s2earth.h` (which provides basic
functionality for converting distances on the Earth's surface to S2-compatible
types).  Instead tests should use the even more basic functionality provided in
`s2testing.h`.  For example:

```c++
  S1Angle tolerance = S2Testing::MetersToAngle(1);
```

### Randomization and Robustness

S2 makes heavy use of randomized tests.  This is particularly useful as a way to
test robustness guarantees.  Typically, a specialized function is written to
generate test cases that are very likely to trigger problems, and then a large
number of cases are checked for correctness.  Here is an example:

```c++
// Chooses a random S2Point that is often near the intersection of one of the
// coodinates planes or coordinate axes with the unit sphere.  (It is possible
// to represent very small perturbations near such points.)
S2Point ChoosePoint() {
  S2Point x = S2Testing::RandomPoint();
  for (int i = 0; i < 3; ++i) {
    if (S2Testing::rnd.OneIn(3)) {
      x[i] *= pow(1e-50, S2Testing::rnd.RandDouble());
    }
  }
  return x.Normalize();
}

TEST(S2, ProjectError) {
  for (int iter = 0; iter < 1000; ++iter) {
    SCOPED_TRACE(StrCat("Iteration ", iter));
    S2Testing::rnd.Reset(iter);  // Easier to reproduce a specific case.
    S2Point a = ChoosePoint();
    S2Point b = ChoosePoint();
    S2Point n = S2::RobustCrossProd(a, b).Normalize();
    S2Point x = S2Testing::SamplePoint(S2Cap(n, S1Angle::Radians(1e-15)));
    S2Point p = S2::Project(x, a, b);
    EXPECT_LT(s2pred::CompareEdgeDistance(
        p, a, b, S1ChordAngle(S2::kProjectPerpendicularError)), 0);
  }
}
```

Note the use of `SCOPED_TRACE` so that any failing iterations can be identified.
This is obviously only useful if the random number seed is deterministic, which
is one reason that `S2Testing` defines its own random number generator available
as `S2Testing::rnd`.  Tests **should not** use random number generators with
non-deterministic seeds, since this makes problems impossible to debug.

Also note the use of `S2Testing::rnd.Reset(iter)` which allows a specific
failing test case to be reproduced more easily by simply replacing 'iter' by the
failing iteration number.  (Note that most S2 code actually calls
`rnd.Reset(iter + 1)` because it turns out that `rnd.Reset(0)` and
`rnd.Reset(1)` have the same effect.)

Observe the use of `S2Testing` helper methods such as `RandomPoint` and
`SamplePoint`.  Also note the use of `pow(10**x, S2Testing::rnd.RandDouble())`,
which produces a value whose *exponent* is uniformly distributed between "x"
and 0.  This is a useful technique for generating numbers that are
well-distributed over a wide range of magnitudes.

# S2 Geometry Pybind11 Bindings

*UNDER DEVELOPMENT*

This directory contains SWIG and pybind11-based Python SDK for S2 Geometry.

## Migration Strategy

The S2 Geometry library is transitioning from SWIG-based bindings to pybind11-based bindings. During this migration:

- **SWIG bindings** (`s2geometry`): The current production bindings, built with CMake. Use `import s2geometry` to access these.
- **pybind11 bindings** (`s2geometry_pybind`): The new bindings under development, built with Bazel. Use `import s2geometry_pybind` to access these.

Once the pybind11 bindings are feature-complete and stable, the SWIG bindings will be deprecated and the pybind11 package will be renamed to `s2geometry` to become the primary Python API.

## User Guide

### Usage Example

```python
import s2geometry_pybind as s2

p1 = s2.S2Point(1.0, 0.0, 0.0)
p2 = s2.S2Point(0.0, 1.0, 0.0)
sum_point = p1 + p2
print(sum_point)
```

### Interface Notes

The Python bindings follow the C++ API closely but with Pythonic conventions:

**Naming Conventions:**
- Core classes exist within the top-level module; we may define submodules for utility classes.
- Class names remain unchanged (e.g., `S2Point`, `S1Angle`, `R1Interval`)
- Method names are converted to snake_case (converted from UpperCamelCase C++ function names)

**Properties vs. Methods:**
- Simple coordinate accessors are properties: `point.x`, `point.y`, `interval.lo`, `interval.hi`
- Properties are always read-only. To create a modified object, use a constructor or factory method.
- Other functions are not properties: `angle.radians()`, `angle.degrees()`, `interval.length()`

**Invalid Values:**
- Invalid inputs to constructions or functions raises `ValueError`.
- Example: `S1Interval(0.0, 4.0)` raises `ValueError` because `4.0 > π`.
- Note: In C++, these conditions trigger `ABSL_DCHECK` assertions. The bindings prevent these assertions from firing by pre-validating inputs.
- Note: Python bindings check for invalid inputs and throw C++ exceptions which are caught by
  pybind and converted to Python exceptions. Exceptions are normally prohibited by the C++
  style guide, but this is the preferred approach for pybind.

**Documentation:**
- Python docstrings provide essential information about parameters, return values, and key behaviors
- For comprehensive documentation including edge cases and algorithmic details, refer to the C++ header files
- The C++ documentation is the authoritative source of truth

**Operators:**
- Standard Python operators work as expected: `+`, `-`, `*`, `==`, `!=`, `<`, `>` (for C++ classes that implement those operators)

**String Representations:**
- `repr()` prefixes the class name and delegates to C++ `operator<<` for the value
- `str()` delegates to C++ `operator<<` for a cleaner output
- Example: `repr(S1Interval(0.0, 2.0))` returns `'S1Interval([0, 2])'` while `str()` returns `'[0, 2]'`

**Vector Inheritance:**
- In C++, various geometry classes inherit from or expose vector types (e.g., `S2Point` inherits from `Vector3_d`, `R2Point` is a type alias for `Vector2_d`, `R1Interval` returns bounds as `Vector2_d`)
- The Python bindings do **not** expose this inheritance hierarchy; it is treated as an implementation detail
- Instead, classes that inherit from a vector expose key functions from the `BasicVector` interface (e.g., `norm()`, `dot_prod()`, `cross_prod()`)
- C++ functions that accept or return a vector object use a Python tuple (of length matching the vector dimension)
- Array indexing operators (e.g., `point[0]`) are not currently supported

**Serialization:**
- The C++ Encoder/Decoder serialization functions are not currently supported

## Development

### Directory Structure

```
python/
├── module.cc                     # Binding module entry point
├── s2point_bindings.cc           # Bindings for S2Point (add more *_bindings.cc as needed)
├── s2geometry_pybind/            # Dir for Python package
│   └── __init__.py               # Package initialization
├── s2point_test.py               # Tests for S2Point (add more *_test.py as needed)
└── BUILD.bazel                   # Build rules for bindings, library, and tests
```

### Building with Bazel (pybind11 bindings)

Bazel can be used for development and testing of the new pybind11-based bindings.

To run all tests:
```bash
cd src
bazel test //python/...
```

### Building with CMake (SWIG bindings)

CMake currently builds **only the SWIG-based bindings** (the legacy `s2geometry` package). The pybind11 bindings are not yet integrated into the CMake build system.

For detailed CMake build instructions, dependency installation, and Python wheel creation, see the [parent directory README](../README.md#build-and-install). Key points:

- Install dependencies: `sudo apt-get install cmake libssl-dev swig python3-dev`
- Enable Python bindings with `-DWITH_PYTHON=ON` when running cmake
- This will build the SWIG-based `s2geometry` package only

Example:
```bash
mkdir build && cd build
cmake -DBUILD_TESTS=yes -DWITH_PYTHON=ON -DCMAKE_PREFIX_PATH=/path/to/absl/install -DCMAKE_CXX_STANDARD=17 ..
make -j $(nproc)
make test ARGS="-j$(nproc)"
sudo make install
```

To run the SWIG tests directly without cmake:
```bash
cd build/python
PYTHONPATH=. python3 ../../src/python/s2geometry_test.py
```

**Note:** Once the pybind11 bindings are complete, they will be integrated into the CMake build system as a replacement for the SWIG bindings.

## Extending the Bindings

To add bindings for a new class:

1. Create `<classname>_bindings.cc` with pybind11 bindings
2. Update `BUILD.bazel` to add a new `pybind_library` target
3. Update `module.cc` to call your binding function
4. Create tests in `<classname>_test.py`

### Binding File Organization

Use the following sections to organize functions within the bindings files and tests. Secondarily, follow the order in which functions are declared in the C++ headers.

1. **Constructors** - Default constructors and constructors with parameters
2. **Factory methods** - Static factory methods (e.g., `from_degrees`, `from_radians`, `zero`, `invalid`)
3. **Properties** - Mutable and read-only properties (e.g., coordinate accessors like `x`, `y`, `lo`, `hi`)
4. **Predicates** - Simple boolean state checks (e.g., `is_empty`, `is_valid`, `is_full`)
5. **Geometric operations** - All other methods including conversions, computations, containment checks, set operations, normalization, and distance calculations
6. **Vector operations** - Methods from the Vector base class (e.g., `norm`, `norm2`, `normalize`, `dot_prod`, `cross_prod`, `angle`). Only applicable to classes that inherit from `util/math/vector.h`
7. **Operators** - Operator overloads (e.g., `==`, `+`, `*`, comparison operators)
8. **String representation** - `__repr__` (which also provides `__str__`), and string conversion methods like `to_string_in_degrees`
9. **Module-level functions** - Standalone functions (e.g., trigonometric functions for S1Angle)

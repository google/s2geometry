# S2 Geometry Pybind11 Bindings

*UNDER DEVELOPMENT*

This directory contains pybind11-based Python SDK for S2 Geometry.

## Migration Strategy

The S2 Geometry library is transitioning from SWIG-based bindings to pybind11-based bindings. During this migration:

- **SWIG bindings** (`s2geometry`): The current production bindings, built with CMake. Use `import s2geometry` to access these.
- **pybind11 bindings** (`s2geometry_pybind`): The new bindings under development, built with Bazel. Use `import s2geometry_pybind` to access these.

Once the pybind11 bindings are feature-complete and stable, the SWIG bindings will be deprecated and the pybind11 package will be renamed to `s2geometry` to become the primary Python API.

## Directory Structure

```
python/
├── bindings/                     # Dir for C++ pybind11 bindings
│   ├── BUILD.bazel               # Binding build rules
│   ├── module.cc                 # Binding module entry point
│   ├── s2point_bindings.cc       # Bindings for S2Point
│   └── ...                       # Bindings for additional classes
├── s2geometry_pybind/            # Dir for Python package
│   └── __init__.py               # Package initialization
├── tests/                        # Dir for Python unit tests
│   ├── BUILD.bazel               # Test build rules
│   ├── s2point_test.py           # Tests for S2Point
│   └── ...                       # Tests for additional classes
└── BUILD.bazel                   # Top-level build rules
```

## Usage Example

```python
import s2geometry_pybind as s2

p1 = s2.S2Point(1.0, 0.0, 0.0)
p2 = s2.S2Point(0.0, 1.0, 0.0)
sum_point = p1 + p2
print(sum_point)
```

## Development

Bazel can be used for development and testing.

To run all tests:
```bash
cd src
bazel test //python/tests:all
```

## Extending the Bindings

To add bindings for a new class:

1. Create `bindings/<classname>_bindings.cc` with pybind11 bindings
2. Update `bindings/BUILD.bazel` to add a new `pybind_library` target
3. Update `bindings/module.cc` to call your binding function
4. Create tests in `tests/<classname>_test.py`
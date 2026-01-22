# S2 Geometry Pybind11 Bindings

*UNDER DEVELOPMENT*

This directory contains SWIG and pybind11-based Python SDK for S2 Geometry.

## Migration Strategy

The S2 Geometry library is transitioning from SWIG-based bindings to pybind11-based bindings. During this migration:

- **SWIG bindings** (`s2geometry`): The current production bindings, built with CMake. Use `import s2geometry` to access these.
- **pybind11 bindings** (`s2geometry_pybind`): The new bindings under development, built with Bazel. Use `import s2geometry_pybind` to access these.

Once the pybind11 bindings are feature-complete and stable, the SWIG bindings will be deprecated and the pybind11 package will be renamed to `s2geometry` to become the primary Python API.

## Directory Structure

```
python/
├── module.cc                     # Binding module entry point
├── s2point_bindings.cc           # Bindings for S2Point (add more *_bindings.cc as needed)
├── s2geometry_pybind/            # Dir for Python package
│   └── __init__.py               # Package initialization
├── s2point_test.py               # Tests for S2Point (add more *_test.py as needed)
└── BUILD.bazel                   # Build rules for bindings, library, and tests
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
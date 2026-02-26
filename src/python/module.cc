#include <pybind11/pybind11.h>

namespace py = pybind11;

// Forward declarations for binding functions
void bind_r1interval(py::module& m);
void bind_r2point(py::module& m);
void bind_r2rect(py::module& m);
void bind_s1interval(py::module& m);
void bind_s2point(py::module& m);

PYBIND11_MODULE(s2geometry_bindings, m) {
  m.doc() = "S2 Geometry Python bindings using pybind11";

  // Bind core geometry classes in dependency order
  bind_r1interval(m);
  bind_r2point(m);
  bind_r2rect(m);
  bind_s1interval(m);
  bind_s2point(m);
}

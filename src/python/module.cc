#include <pybind11/pybind11.h>

namespace py = pybind11;

// Forward declarations for binding functions
void bind_r1interval(py::module& m);
void bind_r2point(py::module& m);
void bind_r2rect(py::module& m);
void bind_s1angle(py::module& m);
void bind_s1interval(py::module& m);
void bind_s2point(py::module& m);
void bind_s2latlng(py::module& m);
void bind_s2cell_id(py::module& m);

PYBIND11_MODULE(s2geometry_bindings, m) {
  m.doc() = "S2 Geometry Python bindings using pybind11";

  // Bind core geometry classes in dependency order.
  // Each class must be registered before classes that use it as a
  // parameter or return type.

  // No dependencies
  bind_r1interval(m);
  bind_r2point(m);
  bind_s1interval(m);
  bind_s2point(m);

  // Deps: r1interval, r2point
  bind_r2rect(m);

  // Deps: s2point
  bind_s1angle(m);

  // Deps: s1angle, s2point
  bind_s2latlng(m);
  bind_s2cell_id(m);
}

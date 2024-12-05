# Generate standard target names.
cmake_policy(SET CMP0078 NEW)
# Honor SWIG_MODULE_NAME via -module flag.
cmake_policy(SET CMP0086 NEW)

# Handle where to install the resulting Python package
if (CALL_FROM_SETUP_PY)
    # The CMakeExtension will set CMAKE_INSTALL_PREFIX to the root
    # of the resulting wheel archive
    set(S2GEOMETRY_INSTALL_PREFIX ${CMAKE_INSTALL_PREFIX})
else()
    # The Python package is installed directly in the folder of the
    # detected interpreter (system, user, or virtualenv)
    set(S2GEOMETRY_INSTALL_PREFIX ${Python3_SITELIB})
endif()

include(${SWIG_USE_FILE})
include_directories(${Python3_INCLUDE_DIRS})

# SWIG does not understand `uint64_t`.  `SWIGWORDSIZE64` controls whether
# it maps `uint64_t` to `unsigned long` or `unsigned long long`.
# (See `stdint.i`.)  `SWIGWORDSIZE64` should really be called
# `INT64_IS_LONG_INT`.  Unices other that macOS use `long`.  Windows uses
# `long long`.  If you get error involving `vector<long long unsigned int>`
# and `vector<long unsigned int>` while compiling `s2PYTHON_wrap.cxx`, then
# this is probably defined incorrectly.
# https://github.com/swig/swig/issues/2923
if (CMAKE_SIZEOF_VOID_P EQUAL 8 AND UNIX AND NOT APPLE)
    set(CMAKE_SWIG_FLAGS "-DSWIGWORDSIZE64")
else()
    set(CMAKE_SWIG_FLAGS "")
endif()
set_property(SOURCE s2.i PROPERTY SWIG_FLAGS "-module" "s2geometry")
set_property(SOURCE s2.i PROPERTY CPLUSPLUS ON)

swig_add_library(s2geometry LANGUAGE python SOURCES s2.i)

swig_link_libraries(s2geometry Python3::Module s2)
enable_testing()
add_test(NAME s2geometry_test COMMAND
         ${Python3_EXECUTABLE}
         "${PROJECT_SOURCE_DIR}/src/python/s2geometry_test.py")
set_property(TEST s2geometry_test PROPERTY ENVIRONMENT
             "PYTHONPATH=$ENV{PYTHONPATH}:${PROJECT_BINARY_DIR}/python")

# Install the wrapper.
install(TARGETS s2geometry DESTINATION ${S2GEOMETRY_INSTALL_PREFIX})

# Install swig-generated Python file (we rename it to __init__.py as it will
# ultimately end up in a directory called s2geometry in site-packages, which will
# serve as the module directory.
install(FILES "${CMAKE_CURRENT_BINARY_DIR}/s2geometry.py"
        DESTINATION ${S2GEOMETRY_INSTALL_PREFIX}
        RENAME __init__.py
        COMPONENT s2geometry)

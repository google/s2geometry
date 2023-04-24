# Set SWIG policies
#cmake_policy(SET CMP0078 NEW)
#cmake_policy(SET CMP0086 NEW)

# Handle where to install the resulting Python package
if(CALL_FROM_SETUP_PY)
    # The CMakeExtension will set CMAKE_INSTALL_PREFIX to the root
    # of the resulting wheel archive
    set(PYWRAPS2_INSTALL_PREFIX ${CMAKE_INSTALL_PREFIX})
else()
    # The Python package is installed directly in the folder of the
    # detected interpreter (system, user, or virtualenv)
    set(PYWRAPS2_INSTALL_PREFIX ${Python3_SITELIB})
endif()

include(${SWIG_USE_FILE})
include_directories(${Python3_INCLUDE_DIRS})

set(CMAKE_SWIG_FLAGS "")
set_property(SOURCE s2.i PROPERTY SWIG_FLAGS "-module" "pywraps2")
set_property(SOURCE s2.i PROPERTY CPLUSPLUS ON)

swig_add_library(pywraps2 LANGUAGE python SOURCES s2.i)

swig_link_libraries(pywraps2 ${Python3_LIBRARIES} s2)
enable_testing()
add_test(NAME pywraps2_test COMMAND
         ${Python3_EXECUTABLE}
         "${PROJECT_SOURCE_DIR}/src/python/pywraps2_test.py")
set_property(TEST pywraps2_test PROPERTY ENVIRONMENT
             "PYTHONPATH=$ENV{PYTHONPATH}:${PROJECT_BINARY_DIR}/python")

# Install the wrapper.
install(TARGETS _pywraps2 DESTINATION ${PYWRAPS2_INSTALL_PREFIX})

# Install swig-generated Python file (we rename it to __init__.py as it will
# ultimately end up in a directory called pywraps2 in site-packages, which will
# serve as the module directory.
install(FILES "${CMAKE_CURRENT_BINARY_DIR}/pywraps2.py"
        DESTINATION ${PYWRAPS2_INSTALL_PREFIX}
        RENAME __init__.py
        COMPONENT pywraps2)

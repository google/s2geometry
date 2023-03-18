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

execute_process(COMMAND ${Python3_EXECUTABLE} -c
                "from distutils import sysconfig; print(sysconfig.get_python_lib(prefix=''))"
                OUTPUT_VARIABLE CMAKE_INSTALL_PYLIBDIR_DEFAULT
                OUTPUT_STRIP_TRAILING_WHITESPACE)

set(CMAKE_INSTALL_PYLIBDIR "${CMAKE_INSTALL_PYLIBDIR_DEFAULT}"
    CACHE STRING "Python site-packages installation directory")

# Install the wrapper.
install(TARGETS _pywraps2 DESTINATION "${CMAKE_INSTALL_PYLIBDIR}")
install(FILES "${PROJECT_BINARY_DIR}/python/pywraps2.py"
        DESTINATION "${CMAKE_INSTALL_PYLIBDIR}")

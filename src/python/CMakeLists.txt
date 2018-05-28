find_package(SWIG)
find_package(PythonInterp)
find_package(PythonLibs)

if (${SWIG_FOUND} AND ${PYTHONLIBS_FOUND})

include(${SWIG_USE_FILE})

#set include for swig
get_property(S2_INCLUDES TARGET s2 PROPERTY INTERFACE_INCLUDE_DIRECTORIES) #FIXME
include_directories(${PYTHON_INCLUDE_PATH} ${S2_INCLUDES})

#lib definition
set(CMAKE_SWIG_FLAGS "")
set_property(SOURCE s2.i PROPERTY SWIG_FLAGS "-module" "pywraps2")
set_property(SOURCE s2.i PROPERTY CPLUSPLUS ON)
swig_add_module(pywraps2 python s2.i)
swig_link_libraries(pywraps2 ${PYTHON_LIBRARIES} s2)

#testing
enable_testing()
add_test(NAME pywraps2_test COMMAND
         ${PYTHON_EXECUTABLE}
	 "${PROJECT_SOURCE_DIR}/src/python/pywraps2_test.py")
set_property(TEST pywraps2_test PROPERTY ENVIRONMENT
             "PYTHONPATH=$ENV{PYTHONPATH}:${PROJECT_BINARY_DIR}/python")

#install
install(TARGETS _pywraps2 DESTINATION share/python)
install(FILES "${PROJECT_BINARY_DIR}/python/pywraps2.py"
	DESTINATION share/python)

endif()

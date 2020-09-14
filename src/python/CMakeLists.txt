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

# Package pip wheel
set(PYWRAPS2_PKG "${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}")
set(PYWRAPS2_DIST "${CMAKE_CURRENT_BINARY_DIR}/dist")
set(PYWRAPS2_SETUP "${CMAKE_CURRENT_BINARY_DIR}/setup.py")

configure_file("${CMAKE_CURRENT_SOURCE_DIR}/setup.py" ${PYWRAPS2_SETUP})
add_custom_command(OUTPUT ${PYWRAPS2_PKG}
                   COMMAND ${CMAKE_COMMAND} -E
                     make_directory ${PYWRAPS2_PKG}
                   COMMAND ${CMAKE_COMMAND} -E
                     copy $<TARGET_FILE:s2> ${PYWRAPS2_PKG}
                   COMMAND ${CMAKE_COMMAND} -E
                     copy "$<TARGET_FILE_DIR:_pywraps2>/*pywraps2.*"
                          ${PYWRAPS2_PKG}
                   COMMAND ${CMAKE_COMMAND} -E
                     copy "${CMAKE_CURRENT_SOURCE_DIR}/__init__.py"
                          ${PYWRAPS2_PKG}
                   DEPENDS s2 _pywraps2)
add_custom_command(OUTPUT ${PYWRAPS2_DIST}
                   COMMAND ${Python3_EXECUTABLE} ${PYWRAPS2_SETUP} bdist_wheel
                   DEPENDS ${PYWRAPS2_PKG} ${PYWRAPS2_SETUP})
add_custom_target(bdist_wheel ALL DEPENDS ${PYWRAPS2_DIST})

# Install the package
install(CODE "execute_process(
                COMMAND ${Python3_EXECUTABLE} ${PYWRAPS2_SETUP} install
                WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})")

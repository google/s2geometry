import os
import sys
from pathlib import Path

import cmake_build_extension
import setuptools

# Extra options passed to the CI/CD pipeline that uses cibuildwheel
CIBW_CMAKE_OPTIONS = []
if "CIBUILDWHEEL" in os.environ and os.environ["CIBUILDWHEEL"] == "1":
    # The manylinux variant runs in Debian Stretch and it uses lib64 folder
    if sys.platform == "linux":
        CIBW_CMAKE_OPTIONS += ["-DCMAKE_INSTALL_LIBDIR=lib"]

setuptools.setup(
    ext_modules=[
        cmake_build_extension.CMakeExtension(
            name="s2geometry",
            install_prefix="s2geometry",
            source_dir=str(Path(__file__).parent.absolute()),
            cmake_configure_options=[
                # This option points CMake to the right Python interpreter, and helps
                # the logic of FindPython3.cmake to find the active version
                # f"-DPython3_ROOT_DIR={Path(sys.prefix)}",
                "-DCALL_FROM_SETUP_PY:BOOL=ON",
                "-DBUILD_SHARED_LIBS:BOOL=OFF",
                "-DCMAKE_POSITION_INDEPENDENT_CODE:BOOL=ON",
                "-DWITH_PYTHON:BOOL=ON",
                "-DBUILD_TESTS:BOOL=OFF",
            ]
            + CIBW_CMAKE_OPTIONS,
        )
    ],
    cmdclass=dict(
        build_ext=cmake_build_extension.BuildExtension,
        sdist=cmake_build_extension.GitSdistFolder,
    ),
)

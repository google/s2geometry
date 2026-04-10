import sys
from pathlib import Path

import skbuild

skbuild.setup(
    cmake_args=[
        # This option points CMake to the right Python interpreter, and helps
        # the logic of FindPython3.cmake to find the active version
        f"-DPython3_ROOT_DIR={Path(sys.prefix)}",
        "-DCALL_FROM_SETUP_PY:BOOL=ON",
        "-DBUILD_SHARED_LIBS:BOOL=OFF",
        "-DCMAKE_POSITION_INDEPENDENT_CODE=ON",
        "-DWITH_PYTHON=ON",
        "-DBUILD_TESTS:BOOL=OFF",
    ],
    cmake_languages=("CXX",),
)

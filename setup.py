import sys
from pathlib import Path

import cmake_build_extension
import setuptools

# This example is compliant with PEP517 and PEP518. It uses the setup.cfg file to store
# most of the package metadata. However, build extensions are not supported and must be
# configured in the setup.py.
setuptools.setup(
    # The resulting "mymath" archive contains two packages: mymath_swig and mymath_pybind.
    # This approach separates the two bindings types, typically just one of them is used.
    ext_modules=[
        cmake_build_extension.CMakeExtension(
            # This could be anything you like, it is used to create build folders
            name="SwigBindings",
            # Name of the resulting package name (import mymath_swig)
            install_prefix="s2geometry",
            # Selects the folder where the main CMakeLists.txt is stored
            # (it could be a subfolder)
            source_dir=str(Path(__file__).parent.absolute()),
            cmake_configure_options=[
                                        # This option points CMake to the right Python interpreter, and helps
                                        # the logic of FindPython3.cmake to find the active version
                                        f"-DPython3_ROOT_DIR={Path(sys.prefix)}",
                                        '-DCALL_FROM_SETUP_PY:BOOL=ON',
                                        '-DBUILD_SHARED_LIBS:BOOL=OFF',
                                        '-DCMAKE_CXX_STANDARD=17',
                                        '-DWITH_PYTHON=ON'
                                    ]
        )
    ],
    cmdclass=dict(
        # Enable the CMakeExtension entries defined above
        build_ext=cmake_build_extension.BuildExtension,
    ),
)

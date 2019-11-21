from setuptools import setup, Distribution


class BinaryDistribution(Distribution):
    def has_ext_modules(self):
        return True


setup(
    name="@PROJECT_NAME@",
    version="@PROJECT_VERSION@",
    packages=["@PROJECT_NAME@"],
    package_data={"": ["*.so", "*.dylib", "*.dll"]},
    distclass=BinaryDistribution
)

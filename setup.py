from setuptools import setup, Distribution


class BinaryDistribution(Distribution):
    def has_ext_modules(foo):
        return True


setup(
    name='s2geometry',
    version="0.0.0",
    packages=['s2geometry'],
    package_dir={'s2geometry': 'build/python'},
    package_data={'': ['_pywraps2.so']},
    distclass=BinaryDistribution
)

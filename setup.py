import os
from setuptools import setup, Distribution


class BinaryDistribution(Distribution):
    def has_ext_modules(foo):
        return True


def get_git_version():
    commit = os.environ.get("TRAVIS_COMMIT")

    if commit is None:
        return ""

    return ".dev+{}".format(commit[:7])


setup(
    name='s2geometry',
    version="0.0.0" + get_git_version(),
    packages=['s2geometry'],
    package_dir={'s2geometry': 'build/python'},
    package_data={'': ['_pywraps2.so']},
    distclass=BinaryDistribution
)

import os

from setuptools import find_packages, setup

import versioneer

setup(
    name='lsforce',
    version=versioneer.get_version(),
    packages=find_packages(),
    scripts=[os.path.join('bin', 'axisem2cps')],
    install_requires=[],
)

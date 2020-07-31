import os

from setuptools import find_packages, setup

setup(
    name='lsforce',
    packages=find_packages(),
    scripts=[os.path.join('bin', 'axisem2cps')],
    install_requires=[],
)

#!/usr/bin/env python3
# encoding: utf-8

try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup

import re
with open('mass_and_me.py') as file:
    version_pattern = re.compile("__version__ = '(.*)'")
    version = version_pattern.search(file.read()).group(1)

with open('README.rst') as file:
    readme = file.read()

setup(
    name='mass_and_me',
    version=version,
    author='Kale Kundert',
    author_email='kale.kundert@ucsf.edu',
    description='',
    long_description=readme,
    url='https://github.com/kalekundert/mass_and_me',
    py_modules=[
        'mass_and_me.py',
    ],
    keywords=[
        'mass_and_me',
    ],
    install_requires=[
        'docopt',
        'nonstdlib',
        'numpy',
        'pandas',
        'pylab',
        'scipy',
    ],
    entry_points = {
        'console_scripts': ['mass_and_me=mass_and_me.main'],
    },
    include_package_data=True,
    license='GPLv3',
    zip_safe=False,
)

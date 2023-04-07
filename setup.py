#!/usr/bin/env python
""" Setuptools setup file for FCpy. """
from setuptools import setup, find_packages
from FCpy.__init__ import (__version__, __pkgname__, __author__,
                           __url__, __license__, __description__)

with open('README.md', mode='rt', encoding='utf-8') as fh:
    __long_description__ = fh.read()

setup(
    name = __pkgname__,
    version = __version__,
    description = __description__,
    long_description = __long_description__,
    long_description_content_type='text/markdown',
    url = __url__,
    author = __author__,
    author_email = 'evan.groopman@nist.gov',
    license = __license__,
    classifiers = [
        'License :: OSI Approved :: NIST License',
        'Programming Language :: Python :: 3',
        'Development Status :: 4 - Beta',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering'
    ],
    keywords = 'python Feldman Cousins FC confidence interval CI',

    install_requires = [
        'numpy',
        'scipy',
    ],
    
    extras_require = {'PyQt5': ['PyQt5']},

    packages = find_packages(include=['FCpy', 'FCpy.*']),
    package_data = {},

    python_requires = '>=3.8',
    tests_require = ['pytest'],

    zip_safe = False
)
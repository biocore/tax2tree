#!/usr/bin/env python
from __future__ import division
from setuptools import setup
import os

__author__ = "Daniel McDonald"
__copyright__ = "Copyright 2011, The tax2tree project"
__credits__ = ["Daniel McDonald"]
__license__ = "BSD"
__version__ = "1.0-dev"
__maintainer__ = "Daniel McDonald"
__email__ = "mcdonadt@colorado.edu"
__status__ = "Development"

classes = """
    Development Status :: 4 - Beta
    License :: OSI Approved :: BSD License
    Topic :: Scientific/Engineering
    Topic :: Scientific/Engineering :: Bio-Informatics
    Programming Language :: Python
    Programming Language :: Python :: 2.7
    Operating System :: Unix
    Operating System :: POSIX
    Operating System :: MacOS :: MacOS X
"""
classifiers = [s.strip() for s in classes.split('\n') if s]

long_description = """The tax2tree project"""

# If readthedocs.org is building the project, we're not able to build the
# required numpy/scipy versions on their machines (nor do we want to, as that
# would take a long time). To build the docs, we don't need the latest versions
# of these dependencies anyways, so we use whatever is in their system's
# site-packages to make scikit-bio importable. See doc/rtd-requirements.txt for
# dependencies that RTD must install in order to build our docs.
#
# Code to check whether RTD is building our project is taken from
# http://read-the-docs.readthedocs.org/en/latest/faq.html
on_rtd = os.environ.get('READTHEDOCS', None) == 'True'
if on_rtd:
    install_requires = []
else:
    install_requires = ['numpy >= 1.7', 'future==0.13.1', 'scikit-bio',
                        'Click']

setup(name='tax2tree',
      version=__version__,
      description='Taxonomy to tree decoration tools',
      author=__maintainer__,
      author_email=__email__,
      maintainer=__maintainer__,
      maintainer_email=__email__,
      url='https://github.com/biocore/tax2tree',
      packages=['t2t'],
      scripts=['scripts/t2t'],
      install_requires=install_requires,
      extras_require={'test': ['nose >= 0.10.1', 'pep8'],
                      'doc': ['Sphinx >= 1.2.2']},
      long_description=long_description)

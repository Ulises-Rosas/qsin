#!/usr/bin/env python

import sys
import platform


import setuptools
from distutils.core import setup


dependencies = [
    'numpy',
    'numba',
    'scikit-learn',
]

with open('README.md') as readme_file:
    readme = readme_file.read()

setup(name = "qsin",
      version = '0.1',
      maintainer = 'Ulises Rosas',
      packages = ['qsin'],
      package_dir = {'qsin': 'src'},
      package_data = {'qsin': ['data/*']} ,
      include_package_data=True,
      install_requires = dependencies,
      zip_safe = False,
    #   entry_points = {
    #     'console_scripts': [
    #         'ggpy   = ggpy.cli:main'
    #         ]
    #   },
    #   scripts=[
    #       './scripts/root_groups.py',
    #   ],
      classifiers = [
          'Programming Language :: Python :: 3'
      ]
    )
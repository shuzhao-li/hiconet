#!/usr/bin/env python

from setuptools import setup

setup(
  name='hiconet',
  version='0.1.1',

  author='Shuzhao Li',
  author_email='shuzhao.li@gmail.com',
  description='Hierachical Community Network, data driven omics integration',
  long_description=open('README.md').read(),
  url='',
  license='BSD',

  keywords='bioinformatics systems biology immunology',

  classifiers=[
    'Development Status :: 4 - Beta',
    'Intended Audience :: Developers',
    'Intended Audience :: Science/Research',
    'License :: OSI Approved :: BSD License',
    'Natural Language :: English',
    'Operating System :: OS Independent',
    'Programming Language :: Python :: 2.7',
    'Topic :: Scientific/Engineering :: Bio-Informatics',
    'Topic :: Software Development :: Libraries :: Python Modules',
  ],

  packages=['hiconet'],
  include_package_data=True,
  zip_safe=True,

  install_requires=[
    'numpy',
    'scipy',
    'pandas',
    'sklearn',
    'leidenalg',
    'igraph',
    'fuzzywuzzy',

  ],

)

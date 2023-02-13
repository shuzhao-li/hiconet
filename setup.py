#!/usr/bin/env python

from setuptools import setup, find_packages

setup(
  name='hiconet',
  version='0.5.4',

  author='Shuzhao Li',
  author_email='shuzhao.li@gmail.com',
  description='Hierachical Community Network, data driven omics integration',
  long_description_content_type="text/markdown",
  long_description=open('README.md').read(),
  url='https://github.com/shuzhao-li/hiconet',
  license='BSD',

  keywords='bioinformatics systems biology immunology',

  classifiers=[
    'Development Status :: 4 - Beta',
    'Intended Audience :: Developers',
    'Intended Audience :: Science/Research',
    'License :: OSI Approved :: BSD License',
    'Natural Language :: English',
    'Operating System :: OS Independent',
    'Programming Language :: Python :: 3.7',
    'Topic :: Scientific/Engineering :: Bio-Informatics',
    'Topic :: Software Development :: Libraries :: Python Modules',
  ],

  # changed from earlier setuptools
  packages=find_packages(exclude=['datasets']),
  
  install_requires=[
    'numpy',
    'scipy',
    'pandas',
    'sklearn',
    'leidenalg',
    'python-igraph',
    'fuzzywuzzy',
    'pyyaml',
    'scanpy',
  ],

  python_requires='>=3',

  data_files=[
    ('datasets/SDY80', ['hiconet/datasets/SDY80/project.yaml', 
                        'hiconet/datasets/SDY80/BTMactivity_SDY80.txt',
                        'hiconet/datasets/SDY80/neut_ab_titer__data_matrix.txt',
                        'hiconet/datasets/SDY80/neut_ab_titer__observation_annotation.txt',
                        'hiconet/datasets/SDY80/SDY80_fcs__data_matrix.txt',
                        'hiconet/datasets/SDY80/SDY80_fcs__observation_annotation.txt',
                        'hiconet/datasets/SDY80/biosample.txt',
                        ]),
    ('datasets/Colombia', ['hiconet/datasets/Colombia/project.yaml', 
                        'hiconet/datasets/Colombia/formatted_BTMs_diag_basel.txt',
                        'hiconet/datasets/Colombia/HILIC_pos_diag_basel.txt',
                        'hiconet/datasets/Colombia/HILIC_pos_diag_basel_feature_annotation.txt',
                        'hiconet/datasets/Colombia/samples_meta.txt',
                        ])
  ],

)

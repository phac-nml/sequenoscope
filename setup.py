#!/usr/bin/env python3
import os
from setuptools import find_packages, setup

author = 'Abdallah Meknas'

classifiers = """
Development Status :: 3 - Alpha
Environment :: Console
License :: OSI Approved :: Apache Software License
Intended Audience :: Science/Research
Topic :: Scientific/Engineering
Topic :: Scientific/Engineering :: Bio-Informatics
Programming Language :: Python
Programming Language :: Python :: 3.7
Programming Language :: Python :: Implementation :: CPython
Operating System :: POSIX :: Linux
""".strip().split('\n')

# Load README.md as long description
with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

def read(fname):
    return open(os.path.join(os.path.dirname(__file__), fname)).read()

# Load package version
exec(open('sequenoscope/version.py').read())

setup(
    name='sequenoscope',
    version='1.0.0',
    python_requires='>=3.7.0',
    setup_requires=['pytest-runner'],
    tests_require=['pytest'],
    packages=find_packages(exclude=['tests']),
    include_package_data=True,
    url='https://github.com/phac-nml/sequenoscope',
    license='Apache-2.0',
    author='Abdallah Meknas',
    author_email='abdallahmeknas@gmail.com',
    description=(
        'Description'),
    long_description=long_description,
    long_description_content_type="text/markdown",   
    keywords='nanopore, ONT, adaptive sampling, sequencing, microbial genomics, metagenomics, visualization, bioinformatics, pipeline, read analysis',
    classifiers=classifiers,
    package_dir={'sequenoscope': 'sequenoscope'},
    package_data={
        "": ["*.csv", "*txt"],
    },

    install_requires=[
        'biopython',
        'numpy',
        'pandas>=1.5',
        'plotly',
        'psutil',
        'pysam',
        'scipy',
        'simplejson',
        'six'
    ],

    entry_points={
        'console_scripts': [
            'sequenoscope=sequenoscope.main:main',
        ],
    },
)
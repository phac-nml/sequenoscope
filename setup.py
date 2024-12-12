#!/usr/bin/env python3
import os
from distutils.core import setup
from setuptools import find_packages

author = 'Abdallah Meknas'

classifiers = """
Development Status :: 3 - Alpha
Environment :: Console
License :: OSI Approved :: GNU General Public License v3 or later (GPLv3+)
Intended Audience :: Science/Research
Topic :: Scientific/Engineering
Topic :: Scientific/Engineering :: Bio-Informatics
Programming Language :: Python
Programming Language :: Python :: 3.7
Programming Language :: Python :: Implementation :: CPython
Operating System :: POSIX :: Linux
""".strip().split('\n')


def read(fname):
    return open(os.path.join(os.path.dirname(__file__), fname)).read()


exec(open('sequenoscope/version.py').read())

setup(
    name='sequenoscope',
    include_package_data=True,
    version='1.0.0',
    python_requires='>=3.7.0',
    setup_requires=['pytest-runner'],
    tests_require=['pytest'],
    packages=find_packages(exclude=['tests']),
    url='https://github.com/phac-nml/sequenoscope',
    license='GPLv3',
    author='Abdallah Meknas',
    author_email='abdallah.meknas@phac-aspc.gc.ca',
    description=(
        'Description'),
    keywords='Keywords',
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
        'pytest',
        'scipy',
        'seqtk',
        'simplejson',
        'six'
    ],

    entry_points={
        'console_scripts': [
            'sequenoscope=sequenoscope.main:main',
        ],
    },
)
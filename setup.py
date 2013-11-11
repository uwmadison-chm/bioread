#!/usr/bin/env python
# coding: utf8
import ez_setup
ez_setup.use_setuptools()

from setuptools import setup, find_packages
setup(
    name = "bioread",
    version = "0.9.3",
    package_dir = {'':'src'},
    packages = find_packages('src'),
    install_requires = [
        "numpy",
    ],
    
    entry_points = {
        'console_scripts': [
            'acq2mat = bioread.runners.acq2mat:main',
            'acq2txt = bioread.runners.acq2txt:main',
            'acq_info = bioread.runners.acq_info:main'
        ]
    },

    package_data = {
    },

    # metadata for upload to PyPI
    author = "Nate Vack",
    author_email = "njvack@wisc.edu",
    description = ("Utilities to read BIOPAC AcqKnowledge files"),
    license = "GPL 2.0",
    keywords = "science research physiological biopac convert library",
    url = "http://github.com/njvack/bioread",
    classifiers = (
        "Development Status :: 4 - Beta",
        "Environment :: Console",
        "Intended Audience :: Developers",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: GNU General Public License (GPL)",
        "Natural Language :: English",
        "Operating System :: OS Independent",
        "Programming Language :: Python",
        "Topic :: Scientific/Engineering :: Information Analysis"
    ),
)

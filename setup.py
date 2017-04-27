#!/usr/bin/env python
# coding: utf8

from setuptools import setup, find_packages
import os


def get_locals(filename):
    l = {}
    exec(open(filename, 'r').read(), {}, l)
    return l

metadata = get_locals(os.path.join('bioread', '_metadata.py'))

setup(
    name="bioread",
    version=metadata['version'],
    packages=find_packages(),
    install_requires=["numpy"],
    entry_points={
        'console_scripts': [
            'acq2mat = bioread.runners.acq2mat:main',
            'acq2txt = bioread.runners.acq2txt:main',
            'acq_info = bioread.runners.acq_info:main',
            'acq_markers = bioread.runners.acq_markers:main',
            'acq2hdf5 = bioread.runners.acq2hdf5:main'
        ]
    },

    # metadata for upload to PyPI
    author=metadata['author'],
    author_email=metadata['author_email'],
    description=("Utilities to read BIOPAC AcqKnowledge files"),
    license="GPL 2.0",
    keywords="science research physiological biopac convert library",
    url="http://github.com/njvack/bioread",
    classifiers=[
        "Development Status :: 5 - Production/Stable",
        "Environment :: Console",
        "Intended Audience :: Developers",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: GNU General Public License (GPL)",
        "Natural Language :: English",
        "Operating System :: OS Independent",
        "Programming Language :: Python",
        "Topic :: Scientific/Engineering :: Information Analysis"
    ],
)

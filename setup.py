#!/usr/bin/env python
# coding: utf8

from setuptools import setup, find_packages
import sys
import os

SETUP_REQUIRES = ['setuptools >= 30.3.0']
SETUP_REQUIRES += ['wheel'] if 'bdist_wheel' in sys.argv else []


def get_locals(filename):
    l = {}
    exec(open(filename, 'r').read(), {}, l)
    return l


metadata = get_locals(os.path.join('bioread', '_metadata.py'))

setup(
    name="bioread",
    setup_requires=SETUP_REQUIRES,
    version=metadata['version'],
    packages=find_packages(),
    keywords="science research physiological biopac convert library",
)

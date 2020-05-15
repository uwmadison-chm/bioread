# coding: utf8
# Part of the bioread package for reading BIOPAC data.
#
# Copyright (c) 2020 Board of Regents of the University of Wisconsin System
#
# Written Nate Vack <njvack@wisc.edu> with research from John Ollinger
# at the Waisman Laboratory for Brain Imaging and Behavior, University of
# Wisconsin-Madison
# Project home: http://github.com/njvack/bioread

# I'm just going to test one file, write it to a temp dir
# Really I just want to make sure this thing basically runs so I don't
# push a totally broken executable in a release

from __future__ import absolute_import
import os
from os import path

from bioread.runners import acq2hdf5

DATA_PATH = path.join(path.dirname(path.abspath(__file__)), "data")

DATA_FILE = path.join(DATA_PATH, 'unicode', 'small-unicode-4.4.0.acq')


def test_acq2hdf5_runs(tmpdir):
    out_file = tmpdir.join("test.hdf5")
    fname = str(out_file)
    acq2hdf5.main([DATA_FILE, fname])
    assert os.stat(fname).st_size > 0

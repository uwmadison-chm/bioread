# coding: utf8
# Part of the bioread package for reading BIOPAC data.
#
# Copyright (c) 2020 Board of Regents of the University of Wisconsin System
#
# Written Nate Vack <njvack@wisc.edu> with research from John Ollinger
# at the Waisman Laboratory for Brain Imaging and Behavior, University of
# Wisconsin-Madison
# Project home: http://github.com/njvack/bioread

# I'm just going to test one file, capture stdout
# Really I just want to make sure this thing basically runs so I don't
# push a totally broken executable in a release

from __future__ import absolute_import
from os import path

from bioread.runners import acq_info

DATA_PATH = path.join(path.dirname(path.abspath(__file__)), "data")

DATA_FILE = path.join(DATA_PATH, 'unicode', 'small-unicode-4.4.0.acq')


def test_acq_info_runs(capsys):
    acq_info.main([DATA_FILE])
    out, err = capsys.readouterr()
    assert len(out) > 0
    assert '2016-02-02' in out

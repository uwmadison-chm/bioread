# coding: utf8
# Part of the bioread package for reading BIOPAC data.
#
# Copyright (c) 2025 Board of Regents of the University of Wisconsin System
#
# Written Nate Vack <njvack@wisc.edu> with research from John Ollinger
# at the Waisman Laboratory for Brain Imaging and Behavior, University of
# Wisconsin-Madison

from bioread.runners import acq_layout

def test_acq_layout_runs_all_files(any_acq_file, capsys):
    acq_layout.main([any_acq_file])
    out, err = capsys.readouterr()
    assert len(out) > 0

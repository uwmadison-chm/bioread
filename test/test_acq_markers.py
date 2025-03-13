# coding: utf8
# Part of the bioread package for reading BIOPAC data.
#
# Copyright (c) 2025 Board of Regents of the University of Wisconsin System
#
# Written Nate Vack <njvack@wisc.edu> with research from John Ollinger
# at the Waisman Laboratory for Brain Imaging and Behavior, University of
# Wisconsin-Madison

# I'm just going to test one file, capture stdout
# Really I just want to make sure this thing basically runs so I don't
# push a totally broken executable in a release

from bioread.runners import acq_markers


def test_acq_markers_runs_all_files(any_acq_file, capsys):
    """Test that acq_markers runs on all data files."""
    acq_markers.main([any_acq_file])
    out, err = capsys.readouterr()
    assert len(out) > 0

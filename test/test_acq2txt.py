# coding: utf8
# Part of the bioread package for reading BIOPAC data.
#
# Copyright (c) 2025 Board of Regents of the University of Wisconsin System
#
# Written Nate Vack <njvack@wisc.edu> with research from John Ollinger
# at the Waisman Laboratory for Brain Imaging and Behavior, University of
# Wisconsin-Madison

from bioread.runners import acq2txt


def test_acq2txt_runs_all_files(any_acq_file, capsys):
    """Test that acq2txt runs on all data files."""
    acq2txt.main([any_acq_file])
    out, err = capsys.readouterr()
    assert len(out) > 0  # Should always produce some output

# coding: utf8
# Part of the bioread package for reading BIOPAC data.
#
# Copyright (c) 2025 Board of Regents of the University of Wisconsin System
#
# Written Nate Vack <njvack@wisc.edu> with research from John Ollinger
# at the Waisman Laboratory for Brain Imaging and Behavior, University of
# Wisconsin-Madison


from bioread.runners import acq_info


def test_acq_info_runs_all_files(any_acq_file, capsys):
    """Test that acq_info runs on all data files."""
    acq_info.main([any_acq_file])
    out, err = capsys.readouterr()

    # All of these should produce output
    assert len(out) > 0

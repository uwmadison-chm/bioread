# coding: utf8
# Part of the bioread package for reading BIOPAC data.
#
# Copyright (c) 2025 Board of Regents of the University of Wisconsin System
#
# Written Nate Vack <njvack@wisc.edu> with research from John Ollinger
# at the Waisman Laboratory for Brain Imaging and Behavior, University of
# Wisconsin-Madison

# I'm just going to test one file, write it to a temp dir
# Really I just want to make sure this thing basically runs so I don't
# push a totally broken executable in a release

from pathlib import Path

from bioread.runners import acq2hdf5

def test_acq2hdf5_runs_all_files(any_acq_file, tmpdir):
    """Test that acq2hdf5 runs on all data files."""

    acq_path = Path(any_acq_file)
    base_name = acq_path.stem
    out_file = str(tmpdir / f"{base_name}.hdf5")
    acq2hdf5.main([any_acq_file, out_file])
    assert Path(out_file).stat().st_size > 0  # Should create a non-empty file

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

from bioread.runners import acq2mat

import pytest
from scipy.io import loadmat


def test_acq2mat_runs_all_files(any_acq_file, tmpdir):
    """Test that acq2mat runs on all data files."""
    acq_path = Path(any_acq_file)
    if acq_path.name == "small-unicode-4.4.0.acq":
        pytest.xfail("This file does not round-trip and I don't know why")

    base_name = acq_path.stem
    out_file = str(tmpdir / f"{base_name}.mat")
    acq2mat.main([any_acq_file, out_file])
    assert Path(out_file).stat().st_size > 0  # Should create a non-empty file
    data = loadmat(out_file)
    assert "channels" in data
    assert "event_markers" in data


def test_acq2mat_with_data_only(any_acq_file, tmpdir):
    """
    Test that acq2mat runs on all data files with data_only and does not write
    headers.
    """
    acq_path = Path(any_acq_file)
    if acq_path.name == "small-unicode-4.4.0.acq":
        pytest.xfail("This file does not round-trip and I don't know why")

    base_name = acq_path.stem
    out_file = str(tmpdir / f"{base_name}.mat")
    acq2mat.main([any_acq_file, out_file, "--data-only"])
    assert Path(out_file).stat().st_size > 0  # Should create a non-empty file
    data = loadmat(out_file)
    assert "channels" in data
    assert "event_markers" not in data

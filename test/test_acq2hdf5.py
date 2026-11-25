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

import bioread
from bioread.runners import acq2hdf5

import numpy as np
import pytest

h5py = pytest.importorskip("h5py")


def test_acq2hdf5_runs_all_files(any_acq_file, tmpdir):
    """Test that acq2hdf5 runs on all data files."""

    acq_path = Path(any_acq_file)
    base_name = acq_path.stem
    out_file = str(tmpdir / f"{base_name}.hdf5")
    acq2hdf5.main([any_acq_file, out_file])
    assert Path(out_file).stat().st_size > 0  # Should create a non-empty file
    bio_data = bioread.read(any_acq_file)
    with h5py.File(out_file, "r") as h5_file:
        assert "channels" in h5_file
        for i, channel in enumerate(bio_data.channels):
            h5_channel_name = f"channel_{i}"
            assert np.array_equal(channel.data, h5_file["channels"][h5_channel_name])


def test_acq2hdf5_rejects_unknown_values_as(
    original_uncompressed_acq_file, tmpdir, capsys
):
    out_file = tmpdir / "out.hdf5"
    with pytest.raises(SystemExit) as excinfo:
        acq2hdf5.main(
            [str(original_uncompressed_acq_file), str(out_file), "--values-as=bogus"]
        )
    # Usage errors exit with status 2, as with argparse.
    assert excinfo.value.code == 2
    assert "values-as" in capsys.readouterr().err
    assert not out_file.exists()


def test_acq2hdf5_exits_cleanly_without_h5py(
    original_uncompressed_acq_file, tmpdir, monkeypatch
):
    monkeypatch.setattr(acq2hdf5, "h5py", None)
    with pytest.raises(SystemExit) as excinfo:
        acq2hdf5.main([str(original_uncompressed_acq_file), str(tmpdir / "out.hdf5")])
    assert excinfo.value.code == 1

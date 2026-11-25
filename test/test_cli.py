# Part of the bioread package for reading BIOPAC data.
#
# Copyright (c) 2025 Board of Regents of the University of Wisconsin System

"""Tests for the bioread command and the options shared by every command."""

import importlib
import io
import sys
from pathlib import Path

import pytest

from bioread import _metadata as meta
from bioread import cli

SUBCOMMANDS = ["info", "markers", "layout", "to-txt", "to-hdf5", "to-mat"]


def run(argv, capsys):
    cli.main(argv)
    return capsys.readouterr()


def test_version(capsys):
    out = run(["--version"], capsys).out
    assert meta.version_description in out


@pytest.mark.parametrize("subcommand", SUBCOMMANDS)
def test_subcommand_help(subcommand, capsys):
    out = run([subcommand, "--help"], capsys).out
    assert f"bioread {subcommand}" in out
    assert "--verbose" in out


@pytest.mark.parametrize("subcommand", ["info", "markers", "layout", "to-txt"])
def test_printing_subcommands(subcommand, original_uncompressed_acq_file, capsys):
    out = run([subcommand, str(original_uncompressed_acq_file)], capsys).out
    assert len(out) > 0


@pytest.mark.parametrize("flag", ["-v", "--verbose", "-d", "--debug"])
def test_verbose_aliases(flag, original_uncompressed_acq_file, capsys):
    run(["info", flag, str(original_uncompressed_acq_file)], capsys)


@pytest.mark.parametrize(
    "indexes",
    [
        ["--channel-indexes", "0,1"],
        ["--channel-indexes=0,1"],
        ["--channel-indexes", "0", "--channel-indexes", "1"],
    ],
)
def test_to_txt_channel_indexes(indexes, original_uncompressed_acq_file, capsys):
    # Options before the file name, so the indexes can't swallow it.
    out = run(["to-txt", *indexes, str(original_uncompressed_acq_file)], capsys).out
    header = out.splitlines()[0].split("\t")
    assert len(header) == 3  # time plus two channels


@pytest.mark.parametrize("subcommand", ["markers", "to-txt"])
def test_outfile(subcommand, original_uncompressed_acq_file, tmp_path, capsys):
    out_file = tmp_path / "out.tsv"
    out = run(
        [subcommand, str(original_uncompressed_acq_file), "-o", str(out_file)], capsys
    ).out
    assert out == ""
    assert out_file.read_text()


@pytest.mark.parametrize("subcommand", ["info", "markers", "layout", "to-txt"])
def test_stdin(subcommand, original_uncompressed_acq_file, capsys, monkeypatch):
    data = Path(original_uncompressed_acq_file).read_bytes()
    monkeypatch.setattr(sys, "stdin", io.TextIOWrapper(io.BytesIO(data)))
    out = run([subcommand, "-"], capsys).out
    assert len(out) > 0


def test_to_hdf5(original_uncompressed_acq_file, tmp_path, capsys):
    pytest.importorskip("h5py")
    out_file = tmp_path / "out.hdf5"
    run(["to-hdf5", str(original_uncompressed_acq_file), str(out_file)], capsys)
    assert out_file.stat().st_size > 0


def test_to_mat(original_uncompressed_acq_file, tmp_path, capsys):
    pytest.importorskip("scipy")
    out_file = tmp_path / "out.mat"
    run(["to-mat", "-c", str(original_uncompressed_acq_file), str(out_file)], capsys)
    assert out_file.stat().st_size > 0


def test_unknown_command_is_a_usage_error(capsys):
    with pytest.raises(SystemExit) as excinfo:
        cli.main(["nope"])
    assert excinfo.value.code == 2


# The legacy command names keep accepting exactly the options they always
# did, print a deprecation notice to stderr, and run the new subcommand.
LEGACY_INVOCATIONS = [
    ("acq_info", ["-d", "{acq}"]),
    ("acq_info", ["--debug", "{acq}"]),
    ("acq_markers", ["-o", "{out}", "{acq}"]),
    ("acq_layout", ["-t", "4", "--max-data=8", "-x", "-d", "{acq}"]),
    ("acq_layout", ["--truncate=0", "--hex", "--debug", "{acq}"]),
    ("acq2txt", ["--channel-indexes=0,1", "-o", "{out}", "--missing-as=NA", "{acq}"]),
    ("acq2txt", ["--channel-indexes", "0", "{acq}"]),
    ("acq2txt", ["--outfile={out}", "{acq}"]),
    (
        "acq2hdf5",
        ["--values-as=raw", "--compress=lzf", "--data-only", "-v", "{acq}", "{out}"],
    ),
    ("acq2hdf5", ["--verbose", "{acq}", "{out}"]),
    ("acq2mat", ["-c", "-s", "--data-only", "{acq}", "{out}"]),
    ("acq2mat", ["--compress", "--single", "{acq}", "{out}"]),
]


@pytest.mark.parametrize("old_name, args", LEGACY_INVOCATIONS)
def test_legacy_commands(
    old_name, args, original_uncompressed_acq_file, tmp_path, capsys
):
    if old_name == "acq2hdf5":
        pytest.importorskip("h5py")
    if old_name == "acq2mat":
        pytest.importorskip("scipy")
    module = importlib.import_module(f"bioread.runners.{old_name}")
    out_file = tmp_path / "out"
    argv = [a.format(acq=original_uncompressed_acq_file, out=out_file) for a in args]
    module.main(argv)
    err = capsys.readouterr().err
    assert f"{old_name} is deprecated" in err
    assert "Use: bioread " in err
    if "{out}" in " ".join(args):
        assert out_file.stat().st_size > 0

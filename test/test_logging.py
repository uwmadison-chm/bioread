# Part of the bioread package for reading BIOPAC data.
#
# Copyright (c) 2025 Board of Regents of the University of Wisconsin System

"""bioread should log like a library: quietly, unless the application asks."""

import logging
import subprocess
import sys
from pathlib import Path

import bioread
from bioread import reader


def test_import_and_read_print_nothing(original_uncompressed_acq_file):
    # Run in a fresh interpreter so no test or plugin has configured logging.
    code = (
        "import bioread, logging, sys; "
        "logging.getLogger('bioread').setLevel(logging.DEBUG); "
        f"bioread.read({original_uncompressed_acq_file!r})"
    )
    result = subprocess.run(
        [sys.executable, "-c", code], capture_output=True, text=True, check=True
    )
    assert result.stdout == ""
    assert result.stderr == ""


def test_library_only_installs_a_null_handler():
    handlers = logging.getLogger("bioread").handlers
    assert any(isinstance(h, logging.NullHandler) for h in handlers)
    assert not any(type(h) is logging.StreamHandler for h in handlers)


def test_modules_use_their_own_loggers():
    assert reader.logger.name == "bioread.reader"


def test_debug_messages_reach_application_handlers(
    original_uncompressed_acq_file, caplog
):
    with caplog.at_level(logging.DEBUG, logger="bioread"):
        bioread.read_headers(original_uncompressed_acq_file)
    assert any(r.name.startswith("bioread.") for r in caplog.records)


def test_recovered_header_errors_are_recorded(
    original_uncompressed_acq_file, tmp_path, caplog
):
    # A truncated file: the graph header reads, the rest doesn't.
    data = Path(original_uncompressed_acq_file).read_bytes()[:3000]
    truncated = tmp_path / "truncated.acq"
    truncated.write_bytes(data)
    with caplog.at_level(logging.WARNING, logger="bioread"):
        r = reader.Reader.read_headers(str(truncated))
    assert r.read_errors
    assert "Error reading headers" in r.read_errors[0]
    assert "While reading" in r.read_errors[0]  # the note added at the source
    assert any(rec.levelno == logging.WARNING for rec in caplog.records)

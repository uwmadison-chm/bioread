# Part of the bioread package for reading BIOPAC data.
#
# Copyright (c) 2025 Board of Regents of the University of Wisconsin System
#
# Written Nate Vack <njvack@wisc.edu> with research from John Ollinger
# at the Waisman Laboratory for Brain Imaging and Behavior, University of
# Wisconsin-Madison

"""Pieces shared by all of bioread's command-line tools."""

import logging
import os
import shlex
import sys
from collections.abc import Iterator
from contextlib import contextmanager
from io import BytesIO
from pathlib import Path
from typing import Annotated, BinaryIO

from cyclopts import App, Parameter

from bioread import _metadata as meta

STDIN_NAME = "-"

AcqFile = Annotated[
    str,
    Parameter(
        allow_leading_hyphen=True,
        help="AcqKnowledge file to read. Use - to read from standard input.",
    ),
]

Verbose = Annotated[
    bool,
    Parameter(
        name=["--verbose", "-v"],
        alias=["--debug", "-d"],
        negative="",
        help="Print debugging messages.",
    ),
]

Outfile = Annotated[
    Path | None,
    Parameter(
        name=["--outfile", "-o"],
        help="Write to this file instead of standard output.",
    ),
]


def make_app(name: str, help: str | None = None) -> App:
    """Build an App with the options every bioread command shares."""
    return App(
        name=name,
        help=help,
        version=meta.version_description,
        result_action="return_value",
    )


def run_app(app: App, argv=None):
    """Run app, exiting quietly if stdout is closed early (as by | head)."""
    try:
        return app(argv)
    except BrokenPipeError:
        # Point stdout at devnull so the interpreter's final flush can't fail.
        devnull = os.open(os.devnull, os.O_WRONLY)
        os.dup2(devnull, sys.stdout.fileno())
        sys.exit(1)


def run_legacy(old_name: str, subcommand: str, argv=None):
    """Run an old command name (acq2txt, ...) as its bioread subcommand.

    The old names take exactly the same arguments as the subcommands, so the
    arguments pass through unchanged. A deprecation notice goes to stderr so
    output piped from stdout is unaffected.
    """
    from bioread import cli

    if argv is None:
        argv = sys.argv[1:]
    elif isinstance(argv, str):
        argv = shlex.split(argv)
    replacement = shlex.join(["bioread", subcommand, *argv])
    print(
        f"{old_name} is deprecated and will be removed in a future release. "
        f"Use: {replacement}",
        file=sys.stderr,
    )
    return cli.main([subcommand, *argv])


def setup_logging(verbose: bool) -> None:
    logger = logging.getLogger("bioread")
    logger.setLevel(logging.DEBUG if verbose else logging.INFO)
    if not logging.getLogger().handlers and not logger.handlers:
        handler = logging.StreamHandler()
        handler.setFormatter(logging.Formatter("%(message)s"))
        logger.addHandler(handler)


def fail(message: str) -> None:
    """Print an error to stderr and exit with status 1."""
    print(f"Error: {message}", file=sys.stderr)
    sys.exit(1)


@contextmanager
def open_acq(path: str) -> Iterator[BinaryIO]:
    """Open an AcqKnowledge file for reading, treating - as standard input.

    Standard input is read fully into memory, since reading needs to seek.
    """
    if path == STDIN_NAME:
        yield BytesIO(sys.stdin.buffer.read())
    else:
        with open(path, "rb") as f:
            yield f

#!/usr/bin/env python
# coding: utf8
# Part of the bioread package for reading BIOPAC data.
#
# Copyright (c) 2025 Board of Regents of the University of Wisconsin System

"""Check the generated Kaitai parser against bioread's own header reading.

Parses every .acq file in test/data with both, and compares the byte order, the
file revision, and every graph header field that the file's revision includes.

Usage:
    python notes/kaitai/gen_ksy.py > notes/kaitai/acq.ksy
    node notes/kaitai/compile_ksy.js notes/kaitai/acq.ksy notes/kaitai
    python notes/kaitai/validate.py
"""

import ctypes
import re
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO_ROOT = HERE.parents[1]
sys.path.insert(0, str(REPO_ROOT))
sys.path.insert(0, str(HERE))

from kaitaistruct import KaitaiStream  # noqa: E402

from bioread.reader import Reader  # noqa: E402

try:
    from acq import Acq
except ImportError:
    sys.exit(
        "No generated parser found. Run gen_ksy.py and compile_ksy.js first "
        "-- see the usage note at the top of this file."
    )


def snake(name):
    s = re.sub(r"(.)([A-Z][a-z]+)", r"\1_\2", name)
    s = re.sub(r"([a-z0-9])([A-Z])", r"\1_\2", s)
    return s.lower()


def as_bytes(value):
    """Byte arrays reach us as bytes from Kaitai and as ints from ctypes.

    gen_ksy.py emits `size: N` for byte arrays, so Kaitai hands back bytes,
    while ctypes gives a list of (possibly signed) ints. Reduce both to bytes
    so we compare values rather than representations.
    """
    if isinstance(value, (bytes, bytearray)):
        return bytes(value)
    return bytes(b & 0xFF for b in value)


def normalize(value):
    if isinstance(value, (bytes, bytearray)):
        return bytes(value)
    if hasattr(value, "__len__") and not isinstance(value, (str, bytes)):
        return list(value)
    return value


def main():
    data_dir = REPO_ROOT / "test" / "data"
    files = sorted(data_dir.glob("**/*.acq"))
    if not files:
        sys.exit(f"No .acq files under {data_dir}")

    compared = 0
    mismatches = []

    for path in files:
        reader = Reader.read_headers(str(path))
        graph_header = reader.graph_header

        with open(path, "rb") as f:
            parsed = Acq(KaitaiStream(f))

            if parsed.is_le != (graph_header.byte_order_char == "<"):
                mismatches.append((path.name, "byte order", None, None))
            if parsed.file_revision != graph_header.file_revision:
                mismatches.append(
                    (
                        path.name,
                        "file revision",
                        graph_header.file_revision,
                        parsed.file_revision,
                    )
                )

            kaitai_header = parsed.graph_header
            versioned = graph_header.__class__._versioned_fields
            for name, _ctype, min_revision in versioned:
                if min_revision > graph_header.file_revision:
                    continue
                compared += 1
                ours = getattr(graph_header._struct, name)
                theirs = getattr(kaitai_header, snake(name))

                if isinstance(ours, bytes):
                    # A c_char array. ctypes stops at the first NUL; Kaitai
                    # hands back the whole fixed-size field.
                    matched = theirs.split(b"\0")[0] == ours
                elif isinstance(ours, ctypes.Array):
                    if isinstance(theirs, bytes):
                        matched = as_bytes(ours) == theirs
                    else:
                        matched = list(ours) == list(theirs)
                else:
                    matched = ours == theirs

                if not matched:
                    mismatches.append(
                        (path.name, name, normalize(ours), normalize(theirs))
                    )

    print(f"files: {len(files)}   graph header fields compared: {compared}")
    print(f"mismatches: {len(mismatches)}")
    for entry in mismatches[:20]:
        print(f"   {entry}")
    return 1 if mismatches else 0


if __name__ == "__main__":
    sys.exit(main())

#!/usr/bin/env python
# coding: utf8
# Part of the bioread package for reading BIOPAC data.
#
# Copyright (c) 2025 Board of Regents of the University of Wisconsin System

"""Transpile bioread's header definitions into a Kaitai Struct spec.

The header layouts in bioread/headers.py are already declarative -- each class
carries a _versioned_fields table of (name, ctype, minimum_file_revision). That
is very nearly a Kaitai Struct spec already, so rather than hand-writing a .ksy
and letting it drift, we generate it from the same tables the reader uses.

Usage:
    python notes/kaitai/gen_ksy.py > notes/kaitai/acq.ksy

See README.md in this directory for what the generated spec can and cannot do.
"""

import ctypes
import re
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT))

from bioread import file_revisions as rev  # noqa: E402
from bioread import headers as bh  # noqa: E402

# Kaitai reads a `repeat: expr` array one element at a time, which is very slow
# for the large opaque padding blobs in these headers. Emitting `size: N`
# instead makes them a single read; the value comes back as bytes rather than a
# list of ints.
BYTES_FOR_BYTE_ARRAYS = True

SCALARS = {
    ctypes.c_int16: "s2",
    ctypes.c_uint16: "u2",
    ctypes.c_int32: "s4",
    ctypes.c_uint32: "u4",
    ctypes.c_int64: "s8",
    ctypes.c_uint64: "u8",
    ctypes.c_double: "f8",
    ctypes.c_float: "f4",
    ctypes.c_byte: "s1",
    ctypes.c_ubyte: "u1",
}

# The headers we know how to describe. Abstract bases and the dispatching
# classes (ChannelHeader, MarkerHeader, ...) are deliberately absent -- in
# Kaitai the revision dispatch happens in a `switch-on`, not a class hierarchy.
HEADER_CLASSES = [
    bh.GraphHeaderPre4,
    bh.GraphHeaderPost4,
    bh.UnknownPaddingHeader,
    bh.ChannelHeaderPre4,
    bh.ChannelHeaderPost4,
    bh.ForeignHeaderPre4,
    bh.ForeignHeaderPost4,
    bh.ChannelDTypeHeader,
    bh.JournalHeaderPre4,
    bh.JournalHeaderPost4,
    bh.JournalLengthHeader,
    bh.MainCompressionHeaderPre4,
    bh.MainCompressionHeaderPost4,
    bh.ChannelCompressionHeader,
    bh.MarkerHeaderPre4,
    bh.MarkerHeaderPost4,
    bh.MarkerPreItemMetadataHeaderPre4,
    bh.MarkerItemMetadataHeader,
    bh.MarkerItemHeaderPre4,
    bh.MarkerItemHeaderPost4,
]

# Every type resolves its endianness from the root, which works out which way
# round the file is by reading the version field both ways.
ENDIAN_META = [
    "    meta:",
    "      endian:",
    "        switch-on: _root.is_le",
    "        cases:",
    "          'true': le",
    "          'false': be",
]

PREAMBLE = """\
# GENERATED FILE -- do not edit.
# Regenerate with: python notes/kaitai/gen_ksy.py > notes/kaitai/acq.ksy
#
# A Kaitai Struct description of the AcqKnowledge header layouts, generated
# from the _versioned_fields tables in bioread/headers.py.
#
# This describes the headers only. Locating them in the file -- summing channel
# data sizes, scanning for the channel dtype headers -- is not expressible in
# Kaitai's expression language and stays in Python. See README.md.
"""


def snake(name):
    """CamelCase / hungarianNotation -> snake_case, as Kaitai wants."""
    s = re.sub(r"(.)([A-Z][a-z]+)", r"\1_\2", name)
    s = re.sub(r"([a-z0-9])([A-Z])", r"\1_\2", s)
    return s.lower()


def field_lines(name, ctype, min_revision):
    """The .ksy attribute lines for one _versioned_fields entry."""
    out = [f"      - id: {snake(name)}"]

    if ctype in SCALARS:
        out.append(f"        type: {SCALARS[ctype]}")
    elif hasattr(ctype, "_type_") and hasattr(ctype, "_length_"):
        element = ctype._type_
        count = ctype._length_
        byte_like = element in (ctypes.c_byte, ctypes.c_ubyte)
        if element is ctypes.c_char or (BYTES_FOR_BYTE_ARRAYS and byte_like):
            # Deliberately no `encoding:` on character arrays. Kaitai requires
            # it to be a compile-time literal, but bioread picks latin1 or
            # utf-8 based on the file revision it just read. Decoding stays the
            # caller's job.
            out.append(f"        size: {count}")
        else:
            out.append(f"        type: {SCALARS[element]}")
            out.append("        repeat: expr")
            out.append(f"        repeat-expr: {count}")
    else:
        raise TypeError(f"don't know how to describe {ctype} (field {name})")

    if min_revision > 0:
        out.append(f"        if: _root.file_revision >= {min_revision}")
    return out


def generate():
    lines = PREAMBLE.splitlines()
    lines += [
        "",
        "meta:",
        "  id: acq",
        "  file-extension: acq",
        "",
        "instances:",
        "  # bioread works out the byte order by reading the version field both",
        "  # ways and keeping whichever gives the smaller number. That heuristic",
        "  # is expressible here, so the spec is self-orienting.",
        "  version_le:",
        "    pos: 2",
        "    type: u4le",
        "  version_be:",
        "    pos: 2",
        "    type: u4be",
        "  is_le:",
        "    value: version_le <= version_be",
        "  file_revision:",
        "    value: 'is_le ? version_le : version_be'",
        "  graph_header:",
        "    pos: 0",
        "    type:",
        f"      switch-on: file_revision >= {rev.V_400B}",
        "      cases:",
        "        'true': graph_header_post4",
        "        'false': graph_header_pre4",
        "",
        "types:",
    ]

    for cls in HEADER_CLASSES:
        lines.append(f"  {snake(cls.__name__)}:")
        lines.extend(ENDIAN_META)
        lines.append("    seq:")
        for name, ctype, min_revision in cls._versioned_fields:
            lines.extend(field_lines(name, ctype, min_revision))
        lines.append("")

    return "\n".join(lines)


if __name__ == "__main__":
    print(generate())

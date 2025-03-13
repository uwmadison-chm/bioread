#!/usr/bin/env python
# coding: utf8
# Part of the bioread package for reading BIOPAC data.
#
# Copyright (c) 2025 Board of Regents of the University of Wisconsin System
#
# Written Nate Vack <njvack@wisc.edu> with research from John Ollinger
# at the Waisman Laboratory for Brain Imaging and Behavior, University of
# Wisconsin-Madison

# This contains the entry point for an executable to print the binary layout
# of an AcqKnowledge file.

"""Print the binary layout of an AcqKnowledge file.

Usage:
    acq_layout [options] <acq_file>
    acq_layout -h | --help
    acq_layout --version

Options:
  -x, --hex         print offsets and lengths in hex
  -d, --debug       print lots of debugging data
  
Note: Using - for <acq_file> reads from stdin.

"""

import csv
import ctypes
from itertools import groupby
import logging
import sys

from docopt import docopt

from bioread.reader import Reader
from bioread.file_revisions import version_string_guess
from bioread import _metadata as meta

logger = logging.getLogger("bioread")
logger.setLevel(logging.INFO)


def main(argv=None):
    if argv is None:
        argv = sys.argv[1:]
    parsed = docopt(__doc__, argv=argv, version=meta.version_description)
    if parsed["--debug"]:
        logger.setLevel(logging.DEBUG)
    logger.debug(parsed)

    alr = AcqLayoutRunner(parsed["<acq_file>"], parsed["--hex"])
    alr.run()


class AcqLayoutRunner:
    def __init__(self, acq_file, hex_output=False):
        self.acq_file = acq_file
        self.hex_output = hex_output

    def run(self):
        reader = Reader.read_headers(self.acq_file)
        self.version_string = reader.version_string
        writer = csv.writer(sys.stdout, delimiter="\t")
        tsv_columns = [
            "Header",
            "Header instance",
            "Byte order",
            "File revision",
            "Version",
            "Header offset",
            "Header struct length",
            "Header effective length",
            "Field name",
            "Field type",
            "Field offset",
            "Field length",
            "Field value",
            "Raw data",
        ]
        writer.writerow(tsv_columns)
        writer.writerows(self._header_rows(reader))
        writer.writerows(self._data_rows(reader))

    def _dec_or_hex(self, value):
        if self.hex_output:
            return hex(value)
        else:
            return value

    def _format_value(self, value):
        if isinstance(value, ctypes.Array):
            return list(value)

        return value

    def _header_rows(self, reader):
        rows = []
        for _, headers in groupby(reader.headers, key=lambda h: h.__class__):
            for header_instance, header in enumerate(headers):
                header_columns = self._format_header(header, header_instance)
                for field in header._struct_class._fields_:
                    field_columns = self._format_field(header, field)
                    rows.append(header_columns + field_columns)
        return rows

    def _format_header(self, header, i):
        header_class = header.__class__

        header_columns = [
            header_class.__name__,
            i,
            header.byte_order_char,
            header.file_revision,
            self.version_string,
            self._dec_or_hex(header.offset),
            self._dec_or_hex(header.struct_length),
            self._dec_or_hex(header.effective_len_bytes),
        ]
        return header_columns

    def _format_field(self, header, field):
        struct_class = header._struct_class
        field_name = field[0]
        field_class = getattr(struct_class, field_name)
        field_type = field[1].__name__
        field_value = getattr(header._struct, field_name)
        field_offset = (
            field_class.offset + header.offset
        )  # This is from the start of the file
        field_size = field_class.size
        data_bytes = header.raw_data[
            field_class.offset : field_class.offset + field_class.size
        ]

        return [
            field_name,
            field_type,
            self._dec_or_hex(field_offset),
            self._dec_or_hex(field_size),
            self._format_value(field_value),
            data_bytes.hex(),
        ]

    def _data_rows(self, reader):
        if reader.is_compressed:
            return self._data_rows_compressed(reader)
        return self._data_rows_uncompressed(reader)

    def _data_rows_compressed(self, reader):
        rows = []
        for i, cch in enumerate(reader.channel_compression_headers):
            row = [
                "Compressed data",
                i,
                cch.byte_order_char,
                cch.file_revision,
                self.version_string,
                self._dec_or_hex(cch.compressed_data_offset),
                self._dec_or_hex(cch.compressed_data_len),
                self._dec_or_hex(cch.uncompressed_data_len),
            ] + ([""] * 6)
            rows.append(row)
        return rows

    def _data_rows_uncompressed(self, reader):
        return [
            [
                "Uncompressed data",
                0,
                reader.byte_order_char,
                reader.file_revision,
                self.version_string,
                self._dec_or_hex(reader.data_start_offset),
                self._dec_or_hex(reader.data_length),
                self._dec_or_hex(reader.data_length),
            ]
            + ([""] * 6)
        ]


if __name__ == "__main__":
    main()

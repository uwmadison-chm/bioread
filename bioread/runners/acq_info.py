#!/usr/bin/env python
# coding: utf8
# Part of the bioread package for reading BIOPAC data.
#
# Copyright (c) 2025 Board of Regents of the University of Wisconsin System
#
# Written Nate Vack <njvack@wisc.edu> with research from John Ollinger
# at the Waisman Laboratory for Brain Imaging and Behavior, University of
# Wisconsin-Madison

# This contains the entry point for an executable to print basic information
# about an AcqKnowledge file.

"""Print some information about an AcqKnowledge file.

Usage:
    acq_info [options] <acq_file>
    acq_info -h | --help
    acq_info --version

Options:
  -d, --debug  print lots of debugging data

Note: Using - for <acq_file> reads from stdin.

"""

import sys
import logging

from io import BytesIO
from docopt import docopt

from bioread.reader import Reader
from bioread import _metadata as meta

logger = logging.getLogger("bioread")
logger.setLevel(logging.INFO)

def main(argv=None):
    if argv is None:
        argv = sys.argv[1:]

    parsed = docopt(
        __doc__,
        argv,
        version=meta.version_description)
    if parsed['--debug']:
        logger.setLevel(logging.DEBUG)
    logger.debug(parsed)
    acq_file = parsed['<acq_file>']
    if acq_file == '-':
        acq_file = BytesIO(sys.stdin.read())
    air = AcqInfoRunner(acq_file)
    air.run()


class AcqInfoRunner:

    def __init__(self, acq_file):
        self.acq_file = acq_file

    def run(self):
        reader = Reader.read_headers(self.acq_file)
        datafile = reader.datafile

        gh = reader.graph_header
        rev = gh.file_revision
        print(f"File revision: {rev} ({reader.version_string}), byte order: {gh.byte_order_char}")
        print(f"Sample time: {gh.sample_time}")
        print(f"Compressed: {gh.compressed}")
        print(f"Number of channels: {gh.channel_count}")
        for channel in datafile.channels:
            print(f"{channel.name}:")
            print(f"\tUnits: {channel.units}")
            print(f"\tNumber of samples: {channel.point_count}")
            print(f"\tFrequency divider: {channel.frequency_divider}")
            print(f"\tData type: {channel.dtype}")


if __name__ == '__main__':
    main()

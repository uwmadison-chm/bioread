#!/usr/bin/env python
# coding: utf8
# Part of the bioread package for reading BIOPAC data.
#
# Copyright (c) 2016 Board of Regents of the University of Wisconsin System
#
# Written Nate Vack <njvack@wisc.edu> with research from John Ollinger
# at the Waisman Laboratory for Brain Imaging and Behavior, University of
# Wisconsin-Madison
# Project home: http://github.com/njvack/bioread

# This contains the entry point for an executable to convert BIOPAC
# AcqKnowledge files into text files
"""Write the data from an AcqKnowledge file channel to a text file.

Usage:
  acq2txt [options] <acq_file>
  acq2txt -h | --help
  acq2txt --version

Options:
  --version                    Show program's version number and exit.
  -h, --help                   Show this help message and exit.
  --channel-indexes=<indexes>  The indexes of the channels to extract.
                               Separate numbers with commas. Default is to
                               extract all channels.
  -o, --outfile=<file>         Write to a file instead of standard out.
  --missing-as=<val>           What value to write where a channel is not
                               sampled. [default: ]

The first column will always be time in seconds. Channel raw values are
converted with scale and offset into native units.
"""

import sys

from bioread.vendor.docopt import docopt

import bioread
from bioread.writers import txtwriter
from bioread import version


def main(argv=None):
    if argv is None:
        argv = sys.argv[1:]

    amr = AcqToTxtRunner(argv)
    amr.run()


class AcqToTxtRunner(object):
    """The little wrapper class that converts acq files to text files"""

    def __init__(self, argv):
        self.argv = argv

    def run(self):
        pargs = docopt(
            __doc__,
            self.argv,
            version=version.description)
        infile = pargs['<acq_file>']
        channel_indexes = None
        if pargs['--channel-indexes']:
            channel_indexes = [
                int(i) for i in pargs['--channel-indexes'].split(',')]
        data = bioread.read(infile, channel_indexes=channel_indexes)
        mval = pargs['--missing-as']
        if pargs['--outfile']:
            with open(pargs['--outfile'], 'w') as f:
                txtwriter.write_text(data, f, channel_indexes, mval)
        else:
            txtwriter.write_text(data, sys.stdout, channel_indexes, mval)


if __name__ == '__main__':
    main()

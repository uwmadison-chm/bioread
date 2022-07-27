#!/usr/bin/env python
# coding: utf8
# Part of the bioread package for reading BIOPAC data.
#
# Copyright (c) 2022 Board of Regents of the University of Wisconsin System
#
# Written Nate Vack <njvack@wisc.edu> with research from John Ollinger
# at the Waisman Laboratory for Brain Imaging and Behavior, University of
# Wisconsin-Madison
# Project home: http://github.com/njvack/bioread

# This contains the entry point for an executable to convert BIOPAC
# AcqKnowledge files into Matlab files.

"""Convert an AcqKnowledge file to a MATLAB file.

Usage:
  acq2mat [options] <acq_file> <mat_file>
  acq2mat -h | --help
  acq2mat --version

Options:
  -c, --compress  Save compressed Matlab file
  --data-only     Only save data and required header information -- do not
                  save event markers.

Note: scipy is required for this program.
"""

from __future__ import absolute_import

import sys
from docopt import docopt

from bioread.reader import Reader
from bioread.writers.matlabwriter import MatlabWriter
from bioread import _metadata as meta


def main(argv=None):
    if argv is None:
        argv = sys.argv[1:]

    amr = AcqToMatRunner(argv)
    amr.run()


class AcqToMatRunner(object):
    """The little wrapper class that converts acq files to mat files"""

    def __init__(self, argv, err=None):
        self.argv = argv
        if err is None:
            err = sys.stderr
        self.err = err

    def run(self):
        old_err = sys.stderr
        sys.stderr = self.err
        pargs = docopt(
            __doc__,
            self.argv,
            version=meta.version_description)
        try:
            import scipy  # noqa -- catch this error before matlabwriter
        except Exception:
            sys.stderr.write("scipy is required for writing matlab files\n")
            sys.exit(1)
        infile = pargs['<acq_file>']
        matfile = pargs['<mat_file>']
        compress = pargs['--compress']

        data = Reader.read(infile).datafile
        MatlabWriter.write_file(
            data,
            matfile,
            compress=compress,
            data_only=pargs['--data-only'])

        sys.stderr = old_err


if __name__ == '__main__':
    main()

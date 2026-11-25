# Part of the bioread package for reading BIOPAC data.
#
# Copyright (c) 2025 Board of Regents of the University of Wisconsin System
#
# Written Nate Vack <njvack@wisc.edu> with research from John Ollinger
# at the Waisman Laboratory for Brain Imaging and Behavior, University of
# Wisconsin-Madison

"""bioread markers (also installed as acq_markers).

Writes the event markers from one or more AcqKnowledge files as tab-separated
values.
"""

import csv
import sys
from typing import Annotated

from cyclopts import Parameter

from bioread.runners._common import (
    Outfile,
    Verbose,
    open_acq,
    read_headers,
    run_legacy,
    setup_logging,
)

FIELDS = [
    "filename",
    "time (s)",
    "label",
    "channel",
    "date_created",
    "type_code",
    "type",
]


def markers(
    acq_files: Annotated[
        list[str],
        Parameter(
            allow_leading_hyphen=True,
            help="AcqKnowledge files to read. Use - to read from standard input.",
        ),
    ],
    /,
    *,
    outfile: Outfile = None,
    verbose: Verbose = False,
) -> None:
    """Print the event markers from AcqKnowledge files."""
    setup_logging(verbose)
    if outfile:
        acq_markers_output_file(acq_files, outfile)
    else:
        acq_markers(acq_files, sys.stdout)


def marker_formatter(acq_filename, graph_sample_msec):
    """Return a function that turns a marker into a dict."""

    def f(marker):
        return {
            "filename": acq_filename,
            "time (s)": (marker.sample_index * graph_sample_msec) / 1000,
            "label": marker.text,
            "channel": marker.channel_name or "Global",
            "date_created": marker.date_created_str,
            "type_code": marker.type_code or "None",
            "type": marker.type,
        }

    return f


def acq_markers_output_file(input_filenames, output_filename):
    with open(output_filename, "w", newline="") as f:
        return acq_markers(input_filenames, f)


def acq_markers(input_filenames, output_stream):
    csv_out = csv.DictWriter(output_stream, FIELDS, delimiter="\t")
    csv_out.writeheader()
    for fname in input_filenames:
        with open_acq(fname) as infile:
            _, df = read_headers(infile)
            mf = marker_formatter(fname, df.graph_header.sample_time)
            for m in df.event_markers:
                csv_out.writerow(mf(m))


def main(argv=None):
    return run_legacy("acq_markers", "markers", argv)


if __name__ == "__main__":
    main()

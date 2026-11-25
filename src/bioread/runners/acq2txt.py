# Part of the bioread package for reading BIOPAC data.
#
# Copyright (c) 2025 Board of Regents of the University of Wisconsin System
#
# Written Nate Vack <njvack@wisc.edu> with research from John Ollinger
# at the Waisman Laboratory for Brain Imaging and Behavior, University of
# Wisconsin-Madison

"""bioread to-txt (also installed as acq2txt)."""

import sys
from typing import Annotated

from cyclopts import Parameter

import bioread
from bioread.runners._common import (
    AcqFile,
    Outfile,
    Verbose,
    open_acq,
    run_legacy,
    setup_logging,
)
from bioread.writers import txtwriter


def _comma_separated_ints(type_, tokens):
    """Accept --channel-indexes 0,2 as well as a repeated --channel-indexes."""
    return [int(part) for t in tokens for part in t.value.split(",") if part.strip()]


ChannelIndexes = Annotated[
    list[int] | None,
    Parameter(
        converter=_comma_separated_ints,
        negative="",
        help="Indexes of the channels to write, separated by commas. "
        "Default is all channels.",
    ),
]


def to_txt(
    acq_file: AcqFile,
    /,
    *,
    channel_indexes: ChannelIndexes = None,
    outfile: Outfile = None,
    missing_as: Annotated[
        str,
        Parameter(help="Value to write where a channel is not sampled."),
    ] = "",
    verbose: Verbose = False,
) -> None:
    """Write the data from an AcqKnowledge file to tab-separated text.

    The first column is always time in seconds. Raw values are converted with
    each channel's scale and offset into native units.
    """
    setup_logging(verbose)
    with open_acq(acq_file) as f:
        data = bioread.read(f, channel_indexes=channel_indexes)
    if outfile:
        with open(outfile, "w") as out:
            txtwriter.write_text(data, out, channel_indexes, missing_as)
    else:
        txtwriter.write_text(data, sys.stdout, channel_indexes, missing_as)


def main(argv=None):
    return run_legacy("acq2txt", "to-txt", argv)


if __name__ == "__main__":
    main()

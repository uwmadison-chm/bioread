# Part of the bioread package for reading BIOPAC data.
#
# Copyright (c) 2025 Board of Regents of the University of Wisconsin System
#
# Written Nate Vack <njvack@wisc.edu> with research from John Ollinger
# at the Waisman Laboratory for Brain Imaging and Behavior, University of
# Wisconsin-Madison

"""bioread to-mat (also installed as acq2mat)."""

from pathlib import Path
from typing import Annotated

from cyclopts import Parameter

from bioread.reader import Reader
from bioread.runners._common import (
    AcqFile,
    Verbose,
    fail,
    open_acq,
    run_legacy,
    setup_logging,
)
from bioread.writers.matlabwriter import MatlabWriter


def to_mat(
    acq_file: AcqFile,
    mat_file: Annotated[Path, Parameter(help="MATLAB file to write.")],
    /,
    *,
    compress: Annotated[
        bool,
        Parameter(name=["--compress", "-c"], negative="", help="Compress the file."),
    ] = False,
    single: Annotated[
        bool,
        Parameter(
            name=["--single", "-s"],
            negative="",
            help="Save data in single precision.",
        ),
    ] = False,
    data_only: Annotated[
        bool,
        Parameter(
            negative="",
            help="Only save data and required metadata; skip event markers, "
            "the journal and most header information.",
        ),
    ] = False,
    verbose: Verbose = False,
) -> None:
    """Convert an AcqKnowledge file to a MATLAB (version 5) file.

    Deprecated: MATLAB reads HDF5 natively, so to-hdf5 is recommended instead.
    Requires scipy (pip install bioread[mat]).
    """
    setup_logging(verbose)
    try:
        import scipy  # noqa: F401 -- fail before reading the whole file
    except ImportError:
        fail("to-mat requires scipy; install it with: pip install bioread[mat]")
    with open_acq(acq_file) as f:
        data = Reader.read(f).datafile
    MatlabWriter.write_file(
        data, str(mat_file), compress=compress, data_only=data_only, single=single
    )


def main(argv=None):
    return run_legacy("acq2mat", "to-mat", argv)


if __name__ == "__main__":
    main()

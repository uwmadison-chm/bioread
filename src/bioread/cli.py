# Part of the bioread package for reading BIOPAC data.
#
# Copyright (c) 2025 Board of Regents of the University of Wisconsin System
#
# Written Nate Vack <njvack@wisc.edu> with research from John Ollinger
# at the Waisman Laboratory for Brain Imaging and Behavior, University of
# Wisconsin-Madison

"""The bioread command and its subcommands.

The older single-purpose commands (acq_info, acq2txt, ...) are still installed
and run the same functions.
"""

from bioread.runners import (
    acq2hdf5,
    acq2mat,
    acq2txt,
    acq_info,
    acq_layout,
    acq_markers,
)
from bioread.runners._common import make_app, run_app

app = make_app("bioread", help="Read and convert BIOPAC AcqKnowledge files.")
app.command(acq_info.info, name="info")
app.command(acq_markers.markers, name="markers")
app.command(acq_layout.layout, name="layout")
app.command(acq2txt.to_txt, name="to-txt")
app.command(acq2hdf5.to_hdf5, name="to-hdf5")
app.command(acq2mat.to_mat, name="to-mat")


def main(argv=None):
    return run_app(app, argv)


if __name__ == "__main__":
    main()

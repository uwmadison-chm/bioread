# Part of the bioread package for reading BIOPAC data.
#
# Copyright (c) 2025 Board of Regents of the University of Wisconsin System
#
# Written Nate Vack <njvack@wisc.edu> with research from John Ollinger
# at the Waisman Laboratory for Brain Imaging and Behavior, University of
# Wisconsin-Madison

"""bioread info (also installed as acq_info)."""

from bioread.runners._common import (
    AcqFile,
    Verbose,
    open_acq,
    read_headers,
    run_legacy,
    setup_logging,
)


def info(acq_file: AcqFile, /, *, verbose: Verbose = False) -> None:
    """Print some information about an AcqKnowledge file."""
    setup_logging(verbose)
    with open_acq(acq_file) as f:
        AcqInfoRunner(f).run()


class AcqInfoRunner:
    def __init__(self, acq_file):
        self.acq_file = acq_file

    def run(self):
        reader, datafile = read_headers(self.acq_file)

        gh = datafile.graph_header
        rev = gh.file_revision
        print(
            f"File revision: {rev} ({reader.version_string}), "
            f"byte order: {gh.byte_order_char}"
        )
        print(f"Sample time: {gh.sample_time}")
        print(f"Compressed: {gh.compressed}")
        print(f"Number of channels: {gh.channel_count}")
        for channel in datafile.channels:
            print(f"{channel.name}:")
            print(f"\tUnits: {channel.units}")
            print(f"\tNumber of samples: {channel.point_count}")
            print(f"\tFrequency divider: {channel.frequency_divider}")
            print(f"\tData type: {channel.dtype}")


def main(argv=None):
    return run_legacy("acq_info", "info", argv)


if __name__ == "__main__":
    main()

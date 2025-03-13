"""
This module contains the data groups for the tests.
"""

import itertools
import pytest

from pathlib import Path

TEST_DIR = Path(__file__).parent
DATA_DIR = TEST_DIR / "data"

PHYSIO = ["physio"]
NOJOURNAL = ["nojournal"]

ALL_DATASETS = PHYSIO + NOJOURNAL
# 3.8.1 is special, the compressed / uncompressed data don't agree in the last
# chunk. I think it's a bug in acqknowledge.
BADEND_VERSIONS = [
    "3.8.1",
]

ORIG_VERSION = ["4.4.0"]

# If you have more versions, send me updates!
NORMAL_VERSIONS = [
    "3.9.1",
    "4.1.0",
    "4.2.0",
    "4.3.0",
    "4.4.0",
    "5.0.1",
]

ALL_VERSIONS = BADEND_VERSIONS + NORMAL_VERSIONS

PRE4_VERSIONS = [v for v in ALL_VERSIONS if v < "4.0.0"]
POST4_VERSIONS = [v for v in ALL_VERSIONS if v >= "4.0.0"]

# Versions with text journals
TEXT_JOURNAL_VERSIONS = [
    "3.8.1",
    "3.9.1",
    "4.1.0",
]

# Versions with HTML journals
HTML_JOURNAL_VERSIONS = [
    "4.2.0",
    "4.3.0",
    "4.4.0",
    "5.0.1",
]

ALL_COMP = ["", "-c"]


def data_file_names(datasets, versions, compressions):
    for dset, ver, comp in itertools.product(datasets, versions, compressions):
        yield data_file_name(dset, ver, comp)


def data_file_name(dataset, version, compression):
    return str(DATA_DIR / dataset / f"{dataset}-{version}{compression}.acq")


def text_journal_files():
    return data_file_names(PHYSIO, TEXT_JOURNAL_VERSIONS, ALL_COMP)


def html_journal_files():
    return data_file_names(PHYSIO, HTML_JOURNAL_VERSIONS, ALL_COMP)


def all_files():
    return [str(p) for p in DATA_DIR.glob("**/*.acq")]


def uncompressed_files():
    return data_file_names(ALL_DATASETS, ALL_VERSIONS, [""])


def compressed_files():
    return data_file_names(ALL_DATASETS, ALL_VERSIONS, ["-c"])


def compressed_uncompressed_pairs():
    return zip(compressed_files(), uncompressed_files())


def good_data_compressed_pairs():
    return zip(
        data_file_names(ALL_DATASETS, NORMAL_VERSIONS, ["-c"]),
        data_file_names(ALL_DATASETS, NORMAL_VERSIONS, [""]),
    )


def buggy_data_compressed_pairs():
    return zip(
        data_file_names(ALL_DATASETS, BADEND_VERSIONS, ["-c"]),
        data_file_names(ALL_DATASETS, BADEND_VERSIONS, [""]),
    )

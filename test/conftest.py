"""
Configuration for pytest.
"""

import pytest

from pathlib import Path

import file_groups

# Base directory for test data
TEST_DIR = Path(__file__).parent
DATA_DIR = TEST_DIR / "data"
FIXTURES_DIR = TEST_DIR / "fixtures"

"""
This file contains the fixtures for the tests.
In addition to some "correct" comparison data, it has various categories of
files we want to test. For example, we have compressed and uncompressed files,
and files from pre-4.0 versions acqknowledge which have text instead of HTML
journals.

Not all the things here are fixtures, per se -- some of them are just iterables
of file paths.
"""


@pytest.fixture(scope="module")
def journal_text():
    with open(FIXTURES_DIR / "journal.txt") as f:
        return f.read()


@pytest.fixture(scope="module")
def journal_html():
    with open(FIXTURES_DIR / "journal.html") as f:
        return f.read()


"""
The "main" test data files are in the data/physio and data/nojournal
directories. The files are named like <dataset>-<version><compression>.acq
where <dataset> is physio or nojournal, <version> is a version number, and
<compression> is empty for uncompressed and '-c' for compressed. We should have
one file for each version and compression status.

The original file was collected in AcqKnowledge 4.4.0, saved in
3.8.1 format, and then converted up to other versions as we and other volunteers
have found other versions of AcqKnowledge.

3.8.1 in particular has a bug where the last chunk of data is incorrect in the
uncompressed files, and that error is propogated to all the other versions
below 4.4.0.

The files in the unicode directory explicitly has some non-ASCII characters.

The files in the misc directory are other files which have unusual revision
numbers and have been the subject of bug reports.
"""


def create_fixture(file_group):
    @pytest.fixture(params=file_group)
    def fx(request):
        return request.param

    return fx


text_journal_file = create_fixture(file_groups.text_journal_files())
html_journal_file = create_fixture(file_groups.html_journal_files())

any_acq_file = create_fixture(file_groups.all_files())


@pytest.fixture
def original_uncompressed_acq_file():
    return str(DATA_DIR / "physio-4.4.0.acq")


compressed_uncompressed_pair = create_fixture(
    file_groups.compressed_uncompressed_pairs()
)

good_data_compressed_pair = create_fixture(file_groups.good_data_compressed_pairs())

buggy_data_compressed_pair = create_fixture(file_groups.buggy_data_compressed_pairs())

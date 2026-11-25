# Part of the bioread package for reading BIOPAC data.
#
# Copyright (c) 2025 Board of Regents of the University of Wisconsin System
#
# Written Nate Vack <njvack@wisc.edu> with research from John Ollinger
# at the Waisman Laboratory for Brain Imaging and Behavior, University of
# Wisconsin-Madison

# The version lives in pyproject.toml; read it from the installed metadata.
from importlib.metadata import PackageNotFoundError, version

try:
    __version__ = version("bioread")
except PackageNotFoundError:  # running from a source tree without installing
    __version__ = "0+unknown"

version_description = f"bioread {__version__}"

author = "Nate Vack"
author_email = "njvack@wisc.edu"
license = "MIT"
copyright = (
    "Copyright (c) 2025 Boards of Regent of the University of Wisconsin System"  # noqa
)
url = "https://github.com/uwmadison-chm/bioread"

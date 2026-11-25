# Part of the bioread package for reading BIOPAC data.
#
# Copyright (c) 2025 Board of Regents of the University of Wisconsin System
#
# Written Nate Vack <njvack@wisc.edu> with research from John Ollinger
# at the Waisman Laboratory for Brain Imaging and Behavior, University of
# Wisconsin-Madison
#
# This script writes an AcqKnowledge file as HDF5, using h5py:
# http://www.h5py.org/

"""bioread to-hdf5 (also installed as acq2hdf5)."""

import logging
from pathlib import Path
from typing import Annotated, Literal

from cyclopts import Parameter

from bioread.runners._common import (
    AcqFile,
    Verbose,
    fail,
    open_acq,
    read_headers,
    run_legacy,
    setup_logging,
)

try:
    import h5py
except ImportError:
    # Checked when the command runs, so that importing this module never exits.
    h5py = None

logger = logging.getLogger(__name__)

COMPRESSION_OPTS = {
    "none": {},
    "gzip": {"compression": "gzip"},
    "lzf": {"compression": "lzf"},
}

SCALE_DATA = {"scaled": True, "raw": False}


def to_hdf5(
    acq_file: AcqFile,
    hdf5_file: Annotated[Path, Parameter(help="HDF5 file to write.")],
    /,
    *,
    values_as: Annotated[
        Literal["scaled", "raw"],
        Parameter(
            help="Save values stored as integers in the file either raw or "
            "scaled. Raw values can be converted using the scale and offset "
            "attributes on each channel; for scaled values, scale and offset "
            "are 1 and 0."
        ),
    ] = "scaled",
    compress: Annotated[
        Literal["gzip", "lzf", "none"],
        Parameter(help="How to compress the data."),
    ] = "gzip",
    data_only: Annotated[
        bool,
        Parameter(
            negative="",
            help="Only save data and required headers; skip the journal and "
            "event markers.",
        ),
    ] = False,
    verbose: Verbose = False,
) -> None:
    """Convert an AcqKnowledge file to an HDF5 file.

    Requires h5py (pip install bioread[hdf5]).
    """
    if h5py is None:
        fail("to-hdf5 requires h5py; install it with: pip install bioread[hdf5]")
    setup_logging(verbose)
    with open_acq(acq_file) as f:
        make_hdf5(
            f, hdf5_file, COMPRESSION_OPTS[compress], SCALE_DATA[values_as], data_only
        )


def make_hdf5(acq_file, hdf5_filename, compression_opts, scale, data_only):
    if h5py is None:
        fail("to-hdf5 requires h5py; install it with: pip install bioread[hdf5]")
    with h5py.File(hdf5_filename, "w") as hdf5_file:
        r, df = read_headers(acq_file)
        hdf5_file.attrs["acq_revision"] = r.file_revision
        hdf5_file.attrs["samples_per_second"] = r.samples_per_second
        hdf5_file.attrs["byte_order"] = r.byte_order_char
        if not data_only:
            hdf5_file.attrs["journal"] = df.journal or ""
        channel_datasets = None
        if r.is_compressed:
            channel_datasets = save_channels_compressed(
                acq_file, r, hdf5_file, compression_opts, scale
            )
        else:
            channel_datasets = save_channels_uncompressed(
                acq_file, r, hdf5_file, compression_opts, scale
            )
        channel_map = dict(
            [[ch.order_num, channel_datasets[i]] for i, ch in enumerate(df.channels)]
        )
        if not data_only:
            save_markers(hdf5_file, df, channel_map)


def cnum_formatter(channel_count):
    return f"channel_{{:0{len(str(channel_count))}d}}"


def create_channel_datasets(grp, channels, compression_opts, scale):
    cfstr = cnum_formatter(len(channels))

    channel_dsets = []
    for i, c in enumerate(channels):
        dset_kwargs = {"dtype": "f8" if scale else c.dtype, **compression_opts}
        cnum_name = cfstr.format(i)
        logger.debug("Creating dataset %s", cnum_name)
        dset = grp.create_dataset(cnum_name, (c.point_count,), **dset_kwargs)
        if c.dtype.kind == "i" and not scale:
            dset.attrs["scale"] = c.raw_scale_factor
            dset.attrs["offset"] = c.raw_offset
        else:
            dset.attrs["scale"] = 1.0
            dset.attrs["offset"] = 0.0
        dset.attrs["name"] = c.name
        dset.attrs["frequency_divider"] = c.frequency_divider
        dset.attrs["units"] = c.units
        dset.attrs["samples_per_second"] = c.samples_per_second
        dset.attrs["channel_number"] = c.order_num
        channel_dsets.append(dset)
    return channel_dsets


def save_channels_uncompressed(acq_file, reader, hdf5_file, compression_opts, scale):
    logger.debug("Saving uncompressed data to hdf5")
    cg = hdf5_file.create_group("/channels")
    df = reader.datafile
    chunker = reader.stream()
    channel_dsets = create_channel_datasets(cg, df.channels, compression_opts, scale)

    for chunk_num, chunk_buffers in enumerate(chunker):
        logger.debug("Got chunk %s", chunk_num)
        for buf, dset in zip(chunk_buffers, channel_dsets, strict=True):
            if scale:
                chan = buf.channel
                dset[buf.channel_slice] = (
                    buf.buffer[:] * chan.raw_scale_factor
                ) + chan.raw_offset
            else:
                dset[buf.channel_slice] = buf.buffer[:]
    return channel_dsets


def save_channels_compressed(acq_file, reader, hdf5_file, compression_opts, scale):
    logger.debug("Saving compressed data to hdf5")
    cg = hdf5_file.create_group("/channels")
    df = reader.datafile
    channel_dsets = create_channel_datasets(cg, df.channels, compression_opts, scale)
    for i, (c, dset) in enumerate(zip(df.channels, channel_dsets, strict=True)):
        logger.debug("Saving channel %s", i)
        # This magically populates c.raw_data; I know this is kind of bad
        reader._read_data(channel_indexes=[i])
        if scale:
            dset[:] = (c.raw_data[:] * c.raw_scale_factor) + c.raw_offset
        else:
            dset[:] = c.raw_data[:]
        # Release the channel's memory
        c.free_data()
    return channel_dsets


def save_markers(hdf5_file, datafile, dset_map):
    markers = datafile.event_markers
    mcount = len(markers)
    marker_name_formatter = f"marker_{{:0{len(str(mcount))}d}}"
    for i, m in enumerate(markers):
        mname = marker_name_formatter.format(i)
        mg = hdf5_file.create_group(f"/event_markers/{mname}")
        mg.attrs["label"] = m.text
        mg.attrs["global_sample_index"] = m.sample_index
        if m.type_code:
            mg.attrs["type_code"] = m.type_code
            mg.attrs["type"] = m.type
        mg.attrs["date_created"] = m.date_created_str
        if m.channel:
            mg.attrs["channel_number"] = m.channel_number
            cdset = dset_map[m.channel_number]
            mg["channel"] = cdset
            mg.attrs["channel_sample_index"] = (
                m.sample_index // cdset.attrs["frequency_divider"]
            )


def main(argv=None):
    return run_legacy("acq2hdf5", "to-hdf5", argv)


if __name__ == "__main__":
    main()

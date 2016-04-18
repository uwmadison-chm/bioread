#!/usr/bin/env python
# coding: utf8
# Part of the bioread package for reading BIOPAC data.
#
# Copyright (c) 2016 Board of Regents of the University of Wisconsin System
#
# Written Nate Vack <njvack@wisc.edu> with research from John Ollinger
# at the Waisman Laboratory for Brain Imaging and Behavior, University of
# Wisconsin-Madison
# Project home: http://github.com/njvack/bioread
#
# This script writes an AcqKnowledge file as HDF5, using h5py:
# http://www.h5py.org/

"""Convert an AcqKnowledge file to an HDF5 file.

Usage:
  acq2hdf5 [options] <acq_file> <hdf5_file>
  acq2hdf5 -h | --help
  acq2hdf5 --version

Options:
  --compress=<method>   How to compress data. Options are gzip, lzf, none.
                        [default: gzip]
  -v, --verbose         Print extra messages for debugging.

Note that channels are saved in raw format, with scale and offset attrs.
"""

import sys
import logging
# Re-adding the handler on reload causes duplicate log messages.
logger = logging.getLogger("bioread.runners.acq2hdf5")
logger.setLevel(logging.INFO)
log_handler = logging.StreamHandler()
log_handler.setLevel(logging.DEBUG)
log_handler.setFormatter(logging.Formatter("%(message)s"))

try:
    import h5py
except ImportError:
    logger.error("acq2hdf5 requires h5py")
    sys.exit(1)


import bioread
from bioread import version
from bioread import reader as br
from bioread.vendor.docopt import docopt

COMPRESSION_OPTS = {
    'none': {},
    'gzip': {'compression': 'gzip'},
    'lzf': {'compression': 'lzf'}
}


def main(argv=None):
    if argv is None:
        argv = sys.argv[1:]

    pargs = docopt(__doc__, argv, version=version.description)
    if pargs['--verbose']:
        logger.setLevel(logging.DEBUG)
    logger.debug(pargs)
    try:
        comp_opts = COMPRESSION_OPTS[pargs['--compress']]
    except KeyError:
        logger.error("Unknown compression: {0}".format(pargs['--compress']))
        sys.exit(1)
    make_hdf5(pargs['<acq_file>'], pargs['<hdf5_file>'], comp_opts)


def make_hdf5(acq_filename, hdf5_filename, compression_opts):
    acq_file = open(acq_filename, 'rb')
    hdf5_file = h5py.File(hdf5_filename, 'w')
    r = br.Reader.read_headers(acq_file)
    df = r.datafile
    hdf5_file.attrs['acq_revision'] = r.file_revision
    hdf5_file.attrs['samples_per_second'] = r.samples_per_second
    hdf5_file.attrs['byte_order'] = r.byte_order_char
    hdf5_file.attrs['journal'] = df.journal or ''
    channel_datasets = None
    if r.is_compressed:
        channel_datasets = save_channels_compressed(
            acq_file, r, hdf5_file, compression_opts)
    else:
        channel_datasets = save_channels_uncompressed(
            acq_file, r, hdf5_file, compression_opts)
    channel_map = dict(
        [[ch.order_num, channel_datasets[i]]
            for i, ch in enumerate(df.channels)])
    save_markers(hdf5_file, df, channel_map)


def cnum_formatter(channel_count):
    return "channel_{{:0{0}d}}".format(len(str(channel_count)))


def create_channel_datasets(grp, channels, compression_opts):
    cfstr = cnum_formatter(len(channels))

    channel_dsets = []
    for i, c in enumerate(channels):
        dset_kwargs = {'dtype': c.dtype}
        dset_kwargs.update(compression_opts)
        cnum_name = cfstr.format(i)
        logger.debug("Creating dataset {0}".format(cnum_name))
        dset = grp.create_dataset(
            cnum_name,
            (c.point_count,),
            **dset_kwargs)
        if c.dtype.kind == 'i':
            dset.attrs['scale'] = c.raw_scale_factor
            dset.attrs['offset'] = c.raw_offset
        else:
            dset.attrs['scale'] = 1.0
            dset.attrs['offset'] = 0.0
        dset.attrs['name'] = c.name
        dset.attrs['frequency_divider'] = c.frequency_divider
        dset.attrs['units'] = c.units
        dset.attrs['samples_per_second'] = c.samples_per_second
        channel_dsets.append(dset)
    return channel_dsets


def save_channels_uncompressed(acq_file, reader, hdf5_file, compression_opts):
    logger.debug("Saving uncompressed data to hdf5")
    cg = hdf5_file.create_group('/channels')
    df = reader.datafile
    channel_dsets = create_channel_datasets(cg, df.channels, compression_opts)
    acq_file.seek(reader.data_start_offset)

    chunker = br.make_chunk_reader(acq_file, df.channels)
    for chunk_num, chunk_buffers in enumerate(chunker):
        logger.debug("Got chunk {0}".format(chunk_num))
        for buf, dset in zip(chunk_buffers, channel_dsets):
            dset[buf.channel_slice] = buf.buffer[:]
    return channel_dsets


def save_channels_compressed(acq_file, reader, hdf5_file, compression_opts):
    logger.debug("Saving compressed data to hdf5")
    cg = hdf5_file.create_group('/channels')
    df = reader.datafile
    channel_dsets = create_channel_datasets(cg, df.channels, compression_opts)
    for i, (c, dset) in enumerate(zip(df.channels, channel_dsets)):
        logger.debug("Saving channel {0}".format(i))
        # This magically populates c.raw_data; I know this is kind of bad
        reader._read_data(channel_indexes=[i])
        dset[:] = c.raw_data[:]
        # Release the channel's memory
        c.free_data()
    return channel_dsets


def save_markers(hdf5_file, datafile, dset_map):
    markers = datafile.markers
    mcount = len(markers)
    marker_name_formatter = "marker_{{:0{0}d}}".format(len(str(mcount)))
    for i, m in enumerate(markers):
        mname = marker_name_formatter.format(i)
        mg = hdf5_file.create_group("/markers/{0}".format(mname))
        mg.attrs['label'] = m.text
        mg.attrs['global_sample_index'] = m.sample_index
        if m.style:
            mg.attrs['style'] = m.style
        if m.channel:
            mg.attrs['channel_index'] = m.channel
            cdset = dset_map[m.channel]
            mg['channel'] = cdset
            mg.attrs['channel_sample_index'] = (
                m.sample_index // cdset.attrs['frequency_divider'])


if __name__ == '__main__':
    argv = sys.argv[1:]
    main(argv)

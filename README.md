# Libraries for reading BIOPAC files

[![DOI](https://zenodo.org/badge/970625.svg)](https://zenodo.org/badge/latestdoi/970625)

These utilities are for reading the files produced by BIOPAC's AcqKnowledge software. Much of the information is based on [Application Note 156](http://www.biopac.com/Manuals/app_pdf/app156.pdf) from BIOPAC; however, newer file formats were decoded through the tireless efforts of John Ollinger, Nate Vack, and Mike Davison. [Updated file structure documentation from Mike](https://github.com/uwmadison-chm/bioread/raw/2e9a312342210c394b9bd3a6dbc220e1e9e89851/notes/acqknowledge_file_structure.pdf) is in the notes/ directory.

This library is mostly concerned with getting you the data, and less so with interpreting UI-related header values.

## Status

As far as I know, this should read any AcqKnowledge file you throw at it. Windows, Mac, uncompressed, compressed, old, new... it should happily read 'em all. If you have trouble with a file, I'd love to get a copy and make bioread work with it.

## Installation

We're up in [pypi](http://pypi.python.org/pypi), so installing should be as simple as:

```
pip install bioread
```

Some of the optional parts of bioread depend on external libraries. `bioread to-hdf5` depends on [h5py](https://www.h5py.org/) and `bioread to-mat` depends on [scipy](https://www.scipy.org/), but as neither of those are core parts of bioread (and can be hairy to get working on some systems), they aren't installed by default. To get them, do:

```
# Just h5py
pip install bioread[hdf5]
# Just scipy
pip install bioread[mat]
# The whole shebang
pip install bioread[all]
```

As of November 2026, we're testing with Python 3.11 — 3.15. bioread 1.0.4 and below should work with Python 2.7 and up.

## API Usage:

* [jupyter notebook: Quick Demo](http://uwmadison-chm.github.io/bioread/bioread_quick_demo.html)

## Command-line usage

Everything is under one `bioread` command:

```
bioread info       Print some information about an AcqKnowledge file
bioread markers    Print the event markers from AcqKnowledge files
bioread to-hdf5    Convert an AcqKnowledge file to an HDF5 file
bioread to-mat     Convert an AcqKnowledge file to a MATLAB file
bioread to-txt     Write the data from an AcqKnowledge file to tab-separated text
bioread layout     Print the binary layout of an AcqKnowledge file
```

`python -m bioread` works too. Run `bioread <command> --help` for the full options of any command.

The older command names (`acq_info`, `acq_markers`, `acq2hdf5`, `acq2mat`, `acq2txt`, `acq_layout`) are still installed and take the same arguments as the matching subcommand, so `acq2txt myfile.acq` and `bioread to-txt myfile.acq` do the same thing.

Every command takes:

* `-h, --help`: show the command's options
* `--version`: print bioread's version
* `-v, --verbose` (or `-d, --debug`): print debugging messages
* `-` as the input file name, to read from standard input. For example: `zcat myfile.acq.gz | bioread info -`

### bioread to-hdf5

If you want to convert files out of AcqKnowledge, this is probably what you want to use -- Matlab can read these out of the box and there are libraries for R and such. This converts the file, storing channels as datasets with names like `/channels/channel_0` and metadata in attributes. Event markers are stored in `/event_markers/marker_X`. Requires h5py (`pip install bioread[hdf5]`).

```
bioread to-hdf5 [OPTIONS] ACQ_FILE HDF5_FILE

--values-as {scaled,raw}   Save values stored as integers either scaled (the default)
                           or raw. Raw values can be converted using the scale and
                           offset attributes on each channel.
--compress {gzip,lzf,none} How to compress the data. [default: gzip]
--data-only                Only save data and required headers; skip the journal and
                           event markers.
```

Note this does *not* need to read the entire dataset into memory, so if you have a 2G dataset, this will work great.

To get the values you see in AcqKnowledge, leave the `--values-as` option to its default ('scaled'). For faster performance, less memory usage, and smaller files, you can use 'raw' and convert the channel later (if you care) with the scale and offset attributes.

Generally, gzip compression seems to work very well, but if you're making something really big you might want to use lzf (worse compression, much faster).

What you'll find in the file:

#### Root-level attributes:

* `file_revision` The internal AckKnowledge file version number
* `samples_per_second` The base sampling rate of the file
* `byte_order` The original file's byte ordering
* `journal` The file's journal data.

#### Channel-level attributes:

* `scale` The scale factor of raw data (for float-type data, will be 1)
* `offset` The offset of raw data (for float-type data, will be 0)
* `frequency_divider` The sampling rate divider for this channel
* `samples_per_second` The channel's sampling rate
* `name` The name of the channel
* `units` The units for the channel
* `channel_number` The display number for the channel (used in markers)

#### Markers

* `label` A text label for the channel
* `type` A description of this marker's type
* `type_code` A short, 4-character code for type
* `global_sample_index` The index, in units of the main sampling rate, of this marker
* `channel` A hard link to the referred channel (only for non-global events)
* `channel_number` The display number for the channel (only for non-global events)
* `channel_sample_index` The in the channel's data where this marker belongs (only for non-global events)


### bioread to-mat

**Note:** I recommend `bioread to-hdf5` for exporting to Matlab. This command is still around because hey: It works.

This creates a Matlab (version 5) file from an AcqKnowledge file. On the back-end, it uses [scipy.io.savemat](http://docs.scipy.org/doc/scipy/reference/generated/scipy.io.savemat.html), so it requires scipy (`pip install bioread[mat]`). Channels are stored in a cell array named 'channels'.

```
bioread to-mat [OPTIONS] ACQ_FILE MAT_FILE

-c, --compress   Save a compressed Matlab file
-s, --single     Save data in single precision
--data-only      Only save data and required metadata; skip event markers, the
                 journal and most header information
```

If you've saved a file as `myfile.mat`, you can, in Matlab:

```
>> data = load('myfile.mat')

data =

              channels: {1x2 cell}
               markers: {1x3 cell}
               headers: [1x1 struct]
    samples_per_second: 1000

>> data.channels{1}

ans =

                 units: 'Percent'
     frequency_divider: 1
    samples_per_second: 1000
                  data: [1x10002 double]
                  name: 'CO2'

>> plot(data.channels{1}.data)

(Plots the data)

>> data.markers{1}

ans =

           style: 'apnd'
    sample_index: 0
           label: 'Segment 1'
         channel: Global
```

### bioread to-txt

Writes the data in an AcqKnowledge file as tab-separated text. By default, all channels (plus a time index) are written. The first column is always time in seconds; raw values are converted with scale and offset into native units.

```
bioread to-txt [OPTIONS] ACQ_FILE

--channel-indexes INDEXES  Indexes of the channels to write, separated by commas
                           (0,2). Default is all channels.
-o, --outfile FILE         Write to a file instead of standard output.
--missing-as VALUE         Value to write where a channel is not sampled. [default: ""]
```

### bioread info

Prints some simple debugging information about an AcqKnowledge file. It'll do its best to print something out even for damaged files.

```
bioread info [OPTIONS] ACQ_FILE
```

### bioread markers

Prints all of the markers in one or more AcqKnowledge files in a tab-separated format, either to standard output or to a file. Fields are:

`filename    time (s)    label    channel    date_created    type_code    type`

```
bioread markers [OPTIONS] ACQ_FILE...

-o, --outfile FILE  Write to a file instead of standard output.
```

### bioread layout

Useful for debugging and/or adding support for new AcqKnowledge versions, `bioread layout` prints the data layout of an AcqKnowledge file, and what header classes are used for all the parts.

```
bioread layout [OPTIONS] ACQ_FILE

-t, --truncate N  Truncate arrays and byte strings to this length. 0 means no
                  truncation. [default: 16]
--max-data N      Maximum number of data bytes to print. [default: 16]
-x, --hex         Print offsets and lengths in hex.
```

## Development

bioread uses [uv](https://docs.astral.sh/uv/) for development and [just](https://just.systems/) as a task runner.

```
uv sync --all-extras   # create .venv with bioread, its extras, and dev tools
just                   # run all checks
just test -k markers   # run a subset of tests; arguments go to pytest
just build             # build the sdist and wheel into dist/
```

The version number lives in `pyproject.toml`; `bioread.__version__` reads it from the installed package metadata.

## Notes

I've tested all the various vintages of files I can think of and find, except very old (AcqKnowledge 2.x) files.

Also, the channel order I read is not the one displayed in the AcqKnowledge interface. Neither the order of the data nor any channel header value I can find seems to entirely control that. I'm gonna just assume it's not a very big deal.

## File Format Documentation

While there's no substite for code diving to see how things really work, I've written some [quick documentation of the file format.](https://github.com/njvack/bioread/blob/master/notes/file_format.md)

In addition, developer Mike Davison did a great job figuring out additional .acq file format information (far more than is implemented in bioread!); his contribuions are in notes/acqknowledge_file_structure.pdf

## Credits

This code was pretty much all written by Nate Vack <njvack@wisc.edu>, with a lot of initial research done by John Ollinger.

## Copyright & Disclaimers

bioread is distributed under the MIT license. For more details, see LICENSE.

BIOPAC and AcqKnowledge are trademarks of BIOPAC Systems, Inc. The authors of this software have no affiliation with BIOPAC Systems, Inc, and that company neither supports nor endorses this software package.

# Libraries for reading BIOPAC files

[![DOI](https://zenodo.org/badge/970625.svg)](https://zenodo.org/badge/latestdoi/970625)

These utilities are for reading the files produced by BIOPAC's AcqKnowledge software. Much of the information is based on [Application Note 156](http://www.biopac.com/Manuals/app_pdf/app156.pdf) from BIOPAC; however, newer file formats were decoded through the tireless efforts of John Ollinger and Nate Vack.

This library is mostly concerned with getting you the data, and less so with interpreting UI-related header values.

## Status

As far as I know, this should read any AcqKnowledge file you throw at it. Windows, Mac, uncompressed, compressed, old, new... it should happily read 'em all. If you have trouble with a file, I'd love to get a copy and make bioread work with it.

## Installation

We're up in [pypi](http://pypi.python.org/pypi), so installing should be as simple as:

```
pip install bioread
```

Some of the optional parts of bioread depend on external libraries. `acq2hdf5` depends on [h5py](https://www.h5py.org/) and `acq2mat` depends on [scipy](https://www.scipy.org/), but as neither of those are core parts of bioread (and can be hairy to get working on some systems), they aren't installed by default. To get them, do:

```
# Just h5py
pip install bioread[hdf5]
# Just scipy
pip install bioread[mat]
# The whole shebang
pip install bioread[all]
```

As of May 2020 (version 2), we now require Python 3.6 or later. Versions 1.0.4 and below should work with Python 2.7 and up.

## API Usage:

* [jupyter notebook: Quick Demo](http://uwmadison-chm.github.io/bioread/bioread_quick_demo.html)

## Command-line usage:

### acq2hdf5

If you want to convert files out of AcqKnowledge, this is probably what you want to use -- Matlab can read these out of the box and there are libraries for R and such. This converts the file, storing channels as datasets with names like `/channels/channel_0` and metadata in attributes. Event markers are stored in `/event_markers/marker_X`

```
Convert an AcqKnowledge file to an HDF5 file.

Usage:
  acq2hdf5 [options] <acq_file> <hdf5_file>
  acq2hdf5 -h | --help
  acq2hdf5 --version

Options:
  --values-as=<type>    Save raw measurement values, stored as integers in the
                        base file, as either 'raw' or 'scaled'. If stored as
                        raw, you can convert to scaled using the scale and
                        offset attributes on the channel. If storing scaled
                        values, scale and offset will be 1 and 0.
                        [default: scaled]
  --compress=<method>   How to compress data. Options are gzip, lzf, none.
                        [default: gzip]
  --data-only           Only save data and required headers -- do not save
                        journal or marker information.
  -v, --verbose         Print extra messages for debugging.
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


### acq2mat

**Note:** I recommend `acq2hdf5` for exporting to Matlab. This program is still around because hey: It works.

This program creates a Matlab (version 5) file from an AcqKnowledge file. On the back-end, it uses [scipy.io.savemat](http://docs.scipy.org/doc/scipy/reference/generated/scipy.io.savemat.html). Channels are stored in a cell array named 'channels'.

```
Convert an AcqKnowledge file to a MATLAB file.

Usage:
  acq2mat [options] <acq_file> <mat_file>
  acq2mat -h | --help
  acq2mat --version

Options:
  -c, --compress  Save compressed Matlab file
  --data-only     Only save data and required header information -- do not
                  save event markers.

Note: scipy is required for this program.
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

### acq2txt

acq2txt will take the data in an AcqKnowledge file and write it to a tab-delimited text file. By default, all channels (plus a time index) will be
written.

```
Write the data from an AcqKnowledge file channel to a text file.

Usage:
  acq2txt [options] <acq_file>
  acq2txt -h | --help
  acq2txt --version

Options:
  --version                    Show program's version number and exit.
  -h, --help                   Show this help message and exit.
  --channel-indexes=<indexes>  The indexes of the channels to extract.
                               Separate numbers with commas. Default is to
                               extract all channels.
  -o, --outfile=<file>         Write to a file instead of standard out.
  --missing-as=<val>           What value to write where a channel is not
                               sampled. [default: ]

The first column will always be time in seconds. Channel raw values are
converted with scale and offset into native units.
```

### acq_info

acq_info prints out some simple debugging information about an AcqKnowledge file. It'll do its best to print something out even for damaged files.

```
Print some information about an AcqKnowledge file.

Usage:
    acq_info [options] <acq_file>
    acq_info -h | --help
    acq_info --version

Options:
  -d, --debug  print lots of debugging data

Note: Using - for <acq_file> reads from stdin.
```

As noted in the usage instructions, acq_info will read from stdin, so if your files are gzipped, you can say:

```
zcat myfile.acq.gz | acq_info -
```

### acq_markers

Prints all of the markers in an AcqKnowlege file to a tab-delimited format, either to stdout or to a specified file. Fields are:

`filename    time (s)    label    channel    style`


```
Print the event markers from an AcqKnowledge file.

Usage:
  acq_markers [options] <file>...
  acq_markers -h | --help
  acq_markers --version

Options:
  -o <file>     Write to a file instead of standard output.
```

Note that this one does not read from stdin; in this case, printing the markers from a large number of files was more important than feeding from `zcat` or something.

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

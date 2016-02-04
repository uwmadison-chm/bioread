# Libraries for reading BIOPAC files

These utilities are for reading the files produced by BIOPAC's AcqKnowledge software. Much of the information is based on  [Application Note 156](http://www.biopac.com/Manuals/app_pdf/app156.pdf) from BIOPAC; however, newer file formats were decoded through the tireless efforts of John Ollinger and Nate Vack.

This library is mostly concerned with getting you the data, and less so with interpreting UI-related header values.

## Status

As far as I know, this should read any AcqKnowledge file you throw at it. Windows, Mac, uncompressed, compressed, old, new... it should happily read 'em all. If you have trouble with a file, I'd love to get a copy and make bioread work with it.

## Installation

We're up in [pypi](http://pypi.python.org/pypi), so installing should be as simple as:

```
pip install bioread
```

Note that bioread requires the excellent [NumPy][http://numpy.scipy.org/] package, and writing Matlab files requires [SciPy](http://scipy.org/).

## Command-line usage:

### acq2mat

This program creates a Matlab (version 5) file from an AcqKnowledge file. On the back-end, it uses [scipy.io.savemat](http://docs.scipy.org/doc/scipy/reference/generated/scipy.io.savemat.html). Channels are stored in a cell array named 'channels'.

```
Convert an AcqKnowledge file to a MATLAB file.

Usage:
  acq2mat [options] <acq_file> <mat_file>
  acq2mat -h | --help
  acq2mat --version

Options:
  -c, --compress  save compressed Matlab file

Note: scipy is required for this program.
```

Then, in Matlab:

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

As noted in the usage instructions, acq2mat will read from stdin, so if your files are gzipped, you can say:

```
zcat myfile.acq.gz | acq2mat - myfile.mat
```

### acq2txt

acq2mat will take the data in an AcqKnowledge file and write it to a text file.

```
Write the data from an AcqKnowledge file channel to a text file.

Usage:
  acq2txt [options] <acq_file>
  acq2txt -h | --help
  acq2txt --version

Options:
  --version          show program's version number and exit
  -h, --help         show this help message and exit
  --channel=CHANNEL  channel number to extract [default: 0]

Writing more than one channel is not supported at the current time, because
different channels can have different sampling rates, and it's hard to know
what to do in those cases.
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
Print the markers from an AcqKnowledge file.

Usage:
  acq_markers [options] <file>...
  acq_markers -h | --help
  acq_markers --version

Options:
  -o <file>     Write to a file instead of standard output.
```

Note that this one does not read from stdin; in this case, printing the markers from a large number of files was more important than feeding from `zcat` or something.

## API usage:

If you want to process the data as NumPy arrays instead, there's an easy API to work with it:

```
>>> import bioread
>>> data = bioread.read_file("myfile.acq")
>>> data.graph_header.file_revision
84
>>> len(data.channels)
2
>>> data.channels[0].samples_per_second
1000.0
>>> len(data.channels[0].data)
10002
>>> data.channels[1].samples_per_second
500.0
>>> len(data.channels[1].data) 
5001
>>> len(data.channels[1].upsampled_data)
10002
>>> data.channels[0].data[0]
1.23
>>> data.channels[0].raw_data[0] # ints are not scaled
13
>>> data.channels[0].name
'CO2'
>>> data.named_channels['CO2'].data[0]
1.23
>>> from bioread.writers import MatlabWriter
>>> MatlabWriter.write_file(data, "myfile.mat") # Creates a matlab file.
```

## Notes

I've tested all the various vintages of files I can think of and find, except very old (AcqKnowledge 2.x) files.

Also, the channel order I read is not the one displayed in the AcqKnowledge interface. Neither the order of the data nor any channel header value I can find seems to entirely control that. I'm gonna just assume it's not a very big deal.

## Credits

This code was pretty much all written by Nate Vack <njvack@wisc.edu>, with a lot of initial research done by John Ollinger, and some test code by Dan Fitch.

Bioread packages a couple great libraries:

[six](http://pythonhosted.org/six/) is copyright (c) 2010-2015 Benjamin Peterson.

[docopt](https://github.com/docopt/docopt) is copyright (c) 2012 Vladimir Keleshev, <vladimir@keleshev.com>.

## Copyright & Disclaimers

bioread is distributed under Version 2 of the GNU Public License. For more details, see LICENSE.

BIOPAC is a trademark of BIOPAC Systems, Inc. The authors of this software have no affiliation with BIOPAC Systems, Inc, and that company neither supports nor endorses this software package.

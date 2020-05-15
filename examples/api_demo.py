#!/usr/bin/env python

# Example taken from the README file.

import bioread

data = bioread.read('myfile.acq')

data.graph_header.file_revision
data.earliest_marker_created_at

len(data.channels)

data.channels[1].samples_per_second
len(data.channels[1].data)
len(data.channels[1].upsampled_data)

data.channels[0].samples_per_second
len(data.channels[0].data)

data.channels[0].data[0]
data.channels[0].raw_data[0]

data.channels[0].name
# let's assume the output is "CO2"
data.named_channels['CO2'].data[0]

from bioread.writers import MatlabWriter
MatlabWriter.write_file(data, "myfile.mat")

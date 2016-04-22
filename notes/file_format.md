# The AcqKnowledge File Format

BIOPAC's AcqKnowledge software saves its files in a fairly complex, badly-documented binary file format. For old versions of the format, there's [Application Note 156](http://www.biopac.com/Manuals/app_pdf/app156.pdf); however, this doesn't cover modern versions of the files and doesn't talk about things file compression and journal data. For those things, the fix was hours of screwing around in hex editors and guessing what values might mean. In the hope that no one else needs to undertake the same task, here are some notes on how all this works.

## File layout

There are two main styles of AcqKnowledge file: compressed and uncompressed. Many of the headers are the same in the two styles; however, the overall file layout is quite different. Whether or not the file is compressed is stored in a byte in the graph header.

### Uncompressed

Layout:

```
Graph Header
Channel Header 1
...
Channel Header N
Foreign Data Header
Channel Datatype Header 1
...
Channel Datatype Header N
Interleaved Channel Data
Marker Header
Marker Item 1
...
Marker Item M
Journal Header
Journal Contents
```

#### Reading uncompressed data

In uncompressed files, the data is *interleaved* — if you have two channels of data sampled at the same rate, the data will be ordered like `0 1 0 1 0 1 ...`.  Each channel, however, can be sampled at a fraction of the data file's main sampling rate; this fractional rate is the *frequency divider* for the channel. In AcqKnowledge, this is limited to powers of two. For a data file where channel 0 has a frequency divider of 1 and channel 1 has a frequency divider of 4, the data stream will look like `0 1 0 0 0 0 1 0 0 0 0 1 ...`. The data can be grouped into a repeating pattern as so: `01000 01000 01000...`

The algorithm to read the data looks like:

1. From each channel's `frequency_divider`, determine the sample pattern as follows:
   1. Compute the least common multiple of the frequency dividers; this yields the `base_length`, or the most times a sample could appear in the pattern.
   2. Make a `base_length x channel_count` dimension matrix. Count from 0 to `base_length` in each row; this yields `pattern_slots`.
   3. For every channel, take the `pattern_slots % frequency_divider` — this will be 0 where the channel is actually sampled. This boolean matrix is the `pattern_mask`.
   4. Make a `channel_slots` matrix in the same shape as `pattern_mask`, containing channel numbers. Select the elements of `channel_slots` where `pattern_mask` is true. This yields `sample_pattern`.
2. Store the total total number of samples in each channel in `channel_samples_remaining`.
3. While the sum of `channel_samples_remaining` > 0:
   1. Count the uses of each channel in `sample_pattern`. If we have fewer available samples in `channel_samples_remaining` in any channels, remove elements from `sample_pattern` for the empty channels, starting from the end of `sample_pattern`.
   2. Read `sample_pattern` samples.
   3. Using `sample_pattern`, select the samples in the data for each channel.
   4: Store those in a permenent channel data array.

Things are slightly more complicated because samples can be either 16-bit ints or 64-bit floats (this is indicated in the channel datatype header) but the basic algorithm is unchanged.

It's important not to take the more obvious route and just trim the end off `sample_pattern` when reading the end of the data, because in some cases, there can be "extra" samples in some channels, yielding a *different* pattern rather than simply a *shortened* pattern at the end of the data.

No, I am not making this up. See the physio test data included in bioread for an example.

### Compressed

In compressed files, the data is segregated by channel, and the data channels come *after* the marker and journal data. Data is compressed with zlib compression, and is little-endian regardless of the endianness of the rest of the file.

Layout:

```
Graph Header
Channel Header 1
...
Channel Header N
Foreign Data Header
Channel Datatype Header 1
...
Channel Datatype Header N
Marker Header
Marker Item 1
...
Marker Item M
Journal Header
Journal Contents
Channel Compression Header 1
...
Channel Compression Header N
```

In general, reading the data from compressed files is much easier than uncompressed files. Seek to the end of a compression header, read *compressed_data_len* bytes, and decompress it into a numpy array.

## Markers

Markers can contain a *channel number* to indicate that it's marking a particular channel of data (or -1 to indicate it's a global marker). This number does *not* correspond to its index in the interleaved or compressed data, but rather to the *order_num* field in the channel header, which controls channel display order in the AcqKnowledge UI.

## Journal data

Stored after the event markers, the journal is fundamentally just a big block of text. In more recent versions (4.2 +) of AcqKnowledge, it's stored as HTML instead of plain text.

## Header index

By and large, the best place to find these headers and their contents defined is in `bioread/headers.py` — here's a quick overview of what the headers are for, though:

* **Graph Header** The settings that apply to an entire file show up here. Both "functional" (eg, msec per sample, total number of channels) and "display" (eg, are we showing a toolbar?) data is stored in this section.
* **Channel Headers** There will be one of these per channel (either raw data or computed) in your data file. Similar to the Graph Header, you'll find both functionally-important and display-only things here.
* **Foreign Data Header** I don't know what this is for, even now.
* **Channel Datatype Headers** One per channel. Defines the datatype of the channel (int or float) and the size (in bytes) of each sample.
* **Channel Compression Headers** One per channel, in compressed files only.
* **Marker Header** Tells us how many event markers we should expect in the file.
* **Marker Item Headers** One of these per event marker in the file.
* **"Post-Marker" Header** I don't know what's in this, but it has a length and it's before the journal data.
* **Journal Header** The header telling us about the journal data.


#!/usr/bin/env python

ID = "$Id:$"[1:-1]

# Written by John Ollinger
#
# University of Wisconsin, 8/16/09

#Copyright (c) 2006-2007, John Ollinger, University of Wisconsin
#All rights reserved.
#
#Redistribution and use in source and binary forms, with or without
#modification, are permitted provided that the following conditions are
#met:
#
#    * Redistributions of source code must retain the above copyright
#       notice, this list of conditions and the following disclaimer.
#
#    * Redistributions in binary form must reproduce the above
#       copyright notice, this list of conditions and the following
#       disclaimer in the documentation and/or other materials provided
#       with the distribution.
#
#THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
#"AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
#LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
#A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
#OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
#SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
#LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
#DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
#THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
#(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
#OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#
# ** This software was designed to be used only for research purposes. **
# ** Clinical uses are not recommended, and have never been evaluated. **
# ** This software comes with no warranties of any kind whatsoever,    **
# ** and may not be useful for anything.  Use it at your own risk!     **
# ** If these terms are not acceptable, you aren't allowed to use the code.**

import os
import sys
import struct
from optparse import OptionParser

from numpy import zeros, fromstring, int16, arange, unique, prod, ones, \
                  reshape, nonzero, float, uint16, float64, int8, float32, \
                  int32
from wbl_util import except_msg


fmt_ctl_windows = '<'
fmt_ctl_mac = '>'
fmt_size = 'h2i'
fmt_main = '3h5d3hdh40s10s19hd3i4h2i4dihb7i' #260cid2i10240c4i'
fmt_rev42 = '260sidi'
fmt_rev44 = 'i10240s4i'
fmt_chan = 'ih40sih2d20si2d3hd128s2h2i'
fmt_chan_rev68 = 'ih40sih2d20si2d3hd128s2h2i'

class BiopacData():

    def __init__(self, biopac):
        self.contents = biopac.contents
        self.hdr = self.contents['hdr']
        self.channels = self.contents['channels']
        self.mhdr =  self.hdr['main']
        self.chdrs = self.hdr['channels']
        self.fhdrs = self.hdr['foreign']
        self.ctype = self.hdr['chan_type']

        self.nchan =    self.mhdr['num_channels']
        self.duration = self.mhdr['sample_time']

        self.chan_descriptions = []
        self.stepsizes = []
        self.units = []
        self.data = []
        self.npts = []
        for ch in xrange(self.nchan):
            self.chan_descriptions.append(self.channels[ch]['description'])
            self.stepsizes.append(self.channels[ch]['stepsize'])
            self.units.append(self.channels[ch]['units'])
            self.data.append(self.channels[ch]['data'])
            self.npts.append(self.channels[ch]['data'].shape[0])
            

class ReadBiopac():
    def __init__(self):
        pass

    def __call__(self):
        self.ParseOptions()
        self.ReadHeader()
        self.ReadData()
        return self

    def ParseOptions(self):
        usage = 'read_biopac <options> <filename>\n' + \
                '\tEnter --help for more information.\n\n'

        optparser = OptionParser(usage)

        optparser.add_option( "-v", "--verbose", action="store_true", \
            dest="verbose",default=False, help='Print stuff to screen.')
        optparser.add_option( "", "--print_roadmap", action="store_true", \
            dest="print_roadmap",default=False, \
            help='Print text file describing the structure of the output ' + \
                 'matlab file.')
        optparser.add_option( "", "--plot", action="store_true", \
            dest="plot",default=False, \
            help='Plot 2 seconds of data as a check.')
        optparser.add_option( "", "--format", action="store", type=str, \
            dest="format",default=None, \
            help='If present, write contents of .acq file to a file of the ' + \
                 'specified by this argument.  Allowable types are "matlab"' + \
                 ' "pickle", and "text". Defaults to "matlab".')
        optparser.add_option( "", "--prefix", action="store", type=str, \
            dest="prefix",default=None, \
            help='Ouput file name.  The appropriate extension will be appended.')
        opts, args = optparser.parse_args()

        if len(args) != 1:
            errstr = "\nExpecting 1 argument:\n" + usage
            sys.stderr.write(errstr)
            sys.exit(1)

        self.filename = args[0]
        self.verbose = opts.verbose
        if opts.prefix is None:
            self.prefix = self.filename.replace('.acq','')
        else:
            self.prefix = opts.prefix
        if opts.format is None:
            self.format = 'matlab'
        else:
            self.format = opts.format
        self.plot = opts.plot
        self.print_roadmap = opts.print_roadmap
        if self.verbose:
            self.print_roadmap = self.verbose

    def ReadHeader(self, filename=None):
        if filename is not None:
            self.filename = filename
        self.f = open(self.filename, 'r')
        self.hdr = {}

#       Read main header.
        self.hdr['main'] = self.ReadMainHeader()
        self.nchan = self.hdr['main']['num_channels']
        if self.verbose:
            print '\nMain header: '
            keys = self.hdr['main'].keys()
            keys.sort()
            for key in keys:
                print '%s: %s' % (key, self.hdr['main'][key])

#       Read Channel headers.
        posn = self.ext_item_header_len
        channel_headers = []
        self.max_buflen = 0
        for channel in xrange(self.hdr['main']['num_channels']):
            self.f.seek(posn)
            chan_hdr = self.ReadChanHeader()
            channel_headers.append(chan_hdr)
            posn += chan_hdr['chan_header_len']
            if chan_hdr['buf_length'] > self.max_buflen:
                self.max_buflen = chan_hdr['buf_length']

        for chan_hdr in channel_headers:
            if not chan_hdr.has_key('var_sample_divider'):
#               Must infer sample spacing.
                chan_hdr['var_sample_divider'] = \
                    int(round(float(self.max_buflen)/float(chan_hdr['buf_length'])))
        self.hdr['channels'] = channel_headers
        if self.verbose:
            ch = 0
            for chan_hdr in channel_headers:
                print '\nchannel: %d, chan_hdr: ' % ch, chan_hdr
                ch += 1

#       Read foreign header.
        self.hdr['foreign'] = self.ReadForeignHeader(posn)
        if self.verbose:
            print '\nforeign: ',self.hdr['foreign']
        posn += self.hdr['foreign']['length']

#        posn += self.HackHeader(posn)

#       Read per_chan_type header.
        self.f.seek(posn)
        self.per_chan_type = []
        for channel in xrange(self.hdr['main']['num_channels']):
            pcthdr, lgth = self.ReadPerChanType()
            self.per_chan_type.append(pcthdr)
            posn += lgth
            if self.verbose:
                print 'Channel: %d, dtype: %s' % (ch, pct['dtype'])
        self.hdr['chan_type'] = self.per_chan_type

        self.chan_hdrs = self.hdr['channels']
        self.duration = self.hdr['main']['sample_time']
        self.stepsizes = []
        for ich in xrange(self.nchan):
            if chan_hdr['buf_length'] > 0:
                self.stepsizes.append(self.duration/chan_hdr['var_sample_divider'])
#                self.stepsizes.append(float(self.duration) /\
#                                        self.chan_hdrs[ich]['buf_length'])
        else:
            self.stepsize = 0.

#        posn += self.nchan*self.len_chantype_hdr

#       Read Markers
#        self.f.seek(posn)
#        self.ReadMarkers()
        self.start_data = posn

    def ReadMainHeader(self):

        hdr = {}
        self.fmt_ctl = fmt_ctl_windows
        fmt = self.fmt_ctl + fmt_size
        data = struct.unpack(fmt, self.f.read(struct.calcsize(fmt)))
        if data[1] > 1000 or data[1] < 0:
#           Invalid file version. This must be a mac. Don't byteswap.
            self.fmt_ctl = fmt_ctl_mac
            fmt = self.fmt_ctl + fmt_size
            self.f.seek(0)
            data = struct.unpack(fmt, self.f.read(struct.calcsize(fmt)))
        self.file_version = data[1]
        if self.verbose:
            print 'File version : %s' % self.file_version

        self.ext_item_header_len = data[2]

        fmt += fmt_main
        if self.file_version > 42:
            fmt += fmt_rev42
        if self.file_version > 44:
            fmt += fmt_rev44

        self.f.seek(0)
        rawdata = self.f.read(struct.calcsize(fmt))
        data = struct.unpack(fmt, rawdata)

        hdr['unused'] = data[0]
        hdr['file_version'] = data[1]
        hdr['ext_item_header_len'] = data[2]
        hdr['num_channels'] = data[3]
        hdr['horiz_axis_type'] = data[4]
        hdr['curr_channel'] = data[5]
        hdr['sample_time'] = data[6]
        hdr['time_offset'] = data[7]
        hdr['time_scale'] = data[8]
        hdr['time_cursor1'] = data[9]
        hdr['time_cursor2'] = data[10]
        hdr['chart_window'] = data[11]
        hdr['measurement'] = data[12]
        hdr['hilite'] = data[13]
        hdr['first_time_offset'] = data[14]
        hdr['rescale'] = data[15]
        hdr['horiz_units1'] = data[16][:data[16].find('\x00')]
        hdr['horiz_units2'] = data[17][:data[17].find('\x00')]
        hdr['in_memory'] = data[18]
        hdr['grid'] = data[19]
        hdr['markers'] = data[20]
        hdr['plot_draft'] = data[21]
        hdr['display_mode'] = data[22]
        hdr['reserved'] = data[23]
        hdr['show_toolbar'] = data[24]
        hdr['show_chan_butt'] = data[25]
        hdr['show_measurement'] = data[26]
        hdr['show_marker'] = data[27]
        hdr['show_journal'] = data[28]
        hdr['cur_x_channel'] = data[29]
        hdr['mmt_precision'] = data[30]
        hdr['measurement_row'] = data[31]
        hdr['mmt'] = data[32]
        hdr['mmt_chan'] = data[33]
        hdr['mmt_calc_opnd1'] = data[34]
        hdr['mmt_calc_opnd2'] = data[35]
        hdr['mmt_calc_op'] = data[36]
        hdr['mmt_calc_constant'] = data[37]
        hdr['new_grid_minor'] = data[38]
        hdr['color_major_grid'] = data[39]
        hdr['color_minor_grid'] = data[40]
        hdr['major_grid_style'] = data[41]
        hdr['minor_grid_style'] = data[42]
        hdr['major_grid_width'] = data[43]
        hdr['minor_grid_width'] = data[44]
        hdr['fixed_units_div'] = data[45]
        hdr['mid_range_show'] = data[46]
        hdr['start_middle_point'] = data[47]
        hdr['offset_point'] = data[48]
        hdr['h_grid'] = data[49]
        hdr['v_grid'] = data[50]
        hdr['enable_wave_tool'] = data[51]
        hdr['horiz_precision'] = data[52]
        hdr['reserved2'] = data[53]
        hdr['overlap_mode'] = data[54]
        hdr['show_hardware'] = data[55]
        hdr['x_auto_plot'] = data[56]
        hdr['x_auto_scroll'] = data[57]
        hdr['start_butt_visible'] = data[58]
        hdr['compressed'] = data[59]
        hdr['always_start_butt_visible'] = data[60]

        if self.file_version > 44:
            hdr['path_video'] = data[61][:data[61].find('\x00')]
            hdr['opt_sync_delay'] = data[62]
            hdr['sync_delay'] = data[63]
            hdr['hrp_paste_measurement'] = data[64]

        if self.file_version > 44:
            hdr['graph_type'] = data[65]
            hdr['mmt_calc_expr'] = data[66][:data[66].find('\x00')]
            hdr['mmt_moment_order'] = data[67]
            hdr['mmt_time_delay'] = data[68]
            hdr['mmt_embed_dim'] = data[69]
            hdr['mmt_mi_delay'] = data[70]

        return hdr

    def ReadChanHeader(self):

        fmt_chan_rev68 = 'ih40sih2d20si2d3hd' # Up to rev 3.7
        fmt_chan_rev76 = 'ih40sih2d20si2d3hd128s2h2ii'
        fmt_chan_rev76 = 'ih40sih2d20si2d3hd68s2h2i'
        fmt_chan_rev76 = fmt_chan_rev68
#        fmt_chan_rev76 = 'ih40sih2d20si2d3hd42s2h3i' # guess 
        if self.file_version > 75:
            self.fmt_chan = self.fmt_ctl + fmt_chan_rev68
        elif self.file_version > 67:
            self.fmt_chan = self.fmt_ctl + fmt_chan_rev68
        else:
            self.fmt_chan = self.fmt_ctl + fmt_chan
        self.lgth_chan_hdr = struct.calcsize(self.fmt_chan)

        rawdata = self.f.read(self.lgth_chan_hdr)
        data = struct.unpack(self.fmt_chan, rawdata)

        chan_hdr = {}
        chan_hdr['chan_header_len'] = data[0]

        chan_hdr['num'] = data[1]
        chan_hdr['comment_text'] = data[2][:data[2].find('\x00')]
        chan_hdr['rgb_color'] = data[3]
        chan_hdr['disp_chan'] = data[4]
        chan_hdr['volt_offset'] = data[5]
        chan_hdr['volt_scale'] = data[6]
        chan_hdr['units_text'] = data[7][:data[7].find('\x00')]
        chan_hdr['buf_length'] = data[8]
        chan_hdr['ampl_scale'] = data[9]
        chan_hdr['ampl_offset'] = data[10]
        chan_hdr['chan_order'] = data[11]
        chan_hdr['disp_size'] = data[12]

        chan_hdr['plot_mode'] = data[13]
        chan_hdr['mid'] = data[14]
        return chan_hdr

    def HackHeader(self, posn):
        self.f.seek(posn)
        N = 10000
        Nc = 0
        fmt = self.fmt_ctl + Nc*'c' + N*'h' 
        from numpy import array, short, int, abs
        data = struct.unpack(fmt, self.f.read(struct.calcsize(fmt)))
        x = array(data[Nc:]).astype(int) #.reshape(N/2, 2)
        Nm = 2*self.nchan
        test = 2*ones(Nm).astype(short)
        for i in xrange(N-Nm):
            if (abs(x[i:i+Nm] - test)).sum() == 0:
                print 'Pattern matched at offset = %d bytes' % (2*i)
#                print 444,x[i:i+Nm], posn+i, self.f.tell()
#                self.DumpBinary(posn + 2*i - 200)
                return 2*i
        print 'Pattern not matched'
        sys.exit()

    def DumpBinary(self, posn):
        self.f.seek(posn + 0)
#        N = (6223 - 4)/2
#        fmt = self.fmt_ctl + N*'i'
        N = 100
        fmt1 = ''
        fmt = self.fmt_ctl + fmt1 + N*'i'
        data = struct.unpack(fmt, self.f.read(struct.calcsize(fmt)))
        M = 20
        print data[:len(fmt1)]
        for i in xrange(0, N, M):
            print '%4d:' % i,
            for j in xrange(M):
                if i+j < N:
                    print ' %6d' %  data[len(fmt1) + i + j],
            print ''

    def ReadForeignHeader(self, posn):
        self.f.seek(posn + 0)
        if self.file_version > 67:
            fmt = self.fmt_ctl + 'ii'
            info = struct.unpack(fmt, self.f.read(struct.calcsize(fmt)))[-2:]
#            print 'Foreign header length: %d' % info[0]
            if info[0] > struct.calcsize(fmt):
                data = self.f.read(info[0]-struct.calcsize(fmt))
            else:
                data = 'None'
#            self.DumpBinary(posn-8)# + info[-2])
#            self.DumpBinary(posn + info[-2])
            return {'length': info[0], 'id': info[1], 'data':data}
        else:
            fmt = self.fmt_ctl + '2h'
            data = struct.unpack(fmt, self.f.read(struct.calcsize(fmt)))
            return {'length': data[0], 'id': data[1]}

    def ReadPerChanType(self):
        type_cvt = {4:'float32', 8:'float64', 18:'int16'}
        fmt = self.fmt_ctl + 'hh'
        data = struct.unpack(fmt, self.f.read(struct.calcsize(fmt)))
        key = 16*(data[1] - 1) + data[0]
        chan_type =  {'size': data[0], \
                      'type':data[1], \
                      'dtype':type_cvt[16*(data[1]-1)+data[0]]}
        return chan_type, struct.calcsize(fmt)

    def ReadMarkers(self):
        fmt = self.fmt_ctl + 'ii'
        data = struct.unpack(fmt, self.f.read(struct.calcsize(fmt)))
        markers = {}
        return markers

    def ReadData(self):
#       Find start of data.
        self.chan_buflen = zeros(self.nchan, int)
        self.sample_spacing = zeros(self.nchan, int)
        data_lgth = 0
        for ich in xrange(self.nchan):
            self.chan_buflen[ich] = self.chan_hdrs[ich]['buf_length']
            self.sample_spacing[ich] = self.chan_hdrs[ich]['var_sample_divider']
            data_lgth += self.chan_buflen[ich]*self.per_chan_type[ich]['size']

        if self.file_version > 67:
            offset = -1
        else:
            offset = 0
        self.f.seek(self.start_data + offset)

#       Create a variable for number of points in each channel. This is not 
#       necessarily the number of samples in the data.
        self.npts = []
        for ich in xrange(self.nchan):
            self.npts.append(self.chan_buflen[ich])

#       Create list of output data arrays, one per channel and data
#       structures defining the type and length of data for each channel.
        self.dtypes = []
        self.dlgths = []
        for ich in xrange(self.nchan):
            self.dtypes.append(self.per_chan_type[ich]['dtype'])
            self.dlgths.append(self.per_chan_type[ich]['size'])

        for ich in xrange(self.nchan):
            scl = self.chan_hdrs[ich]['ampl_scale']
            offset = self.chan_hdrs[ich]['ampl_offset']
            voffset = self.chan_hdrs[ich]['volt_offset']
            xdtype = self.per_chan_type[ich]['dtype']
            if xdtype == 'float64':
                xmin = 0.
                xmax = 0.
            else:
                xmin = -32768*scl+offset
                xmax =  32767*scl + offset

#       Variable sampling intervals.
        self.SetupVarsampBlock()
        datalen = 0
        for ich in xrange(self.nchan):
            datalen += self.dlgths[ich]*self.chan_buflen[ich]

        nblock = datalen/self.raw_blocklen/2
        srawdata = self.f.read(datalen)
        pad_lgth = len(srawdata) % (2*nblock*self.raw_blocklen)
        xx = 2*nblock*self.raw_blocklen

        if pad_lgth:
#           Pad with zeros to ensure last block is full-length.
            pad = zeros(pad_lgth, int8)
            srawdata += pad.tostring()
        rawdata = fromstring(srawdata, uint16)
        blkoff = zeros(self.nchan, int)
        datoff = 0
        rawlen = self.raw_blocklen
        data = []
        for ich in xrange(self.nchan):
            data.append(zeros(self.chan_buflen[ich], float))
        for blk in xrange(nblock):
            for ich in xrange(self.nchan):
                Nsamp = self.Nsamps[ich]
                dat = rawdata[datoff:datoff+self.raw_blocklen].take(self.idcs_in[ich])
#               Convert to double (or keep as int16)
                dat = fromstring(dat.tostring(), self.dtypes[ich])
                data[ich][blkoff[ich]:blkoff[ich]+Nsamp] =  \
                                dat.astype(float)* \
                                self.chan_hdrs[ich]['ampl_scale'] + \
                                self.chan_hdrs[ich]['ampl_offset']
                blkoff[ich] += Nsamp
            datoff += self.raw_blocklen
#        sys.stdout.write('\n')

        self.contents = {'hdr':self.hdr, 'nchan':self.nchan, 'channels':[]}
        for ich in xrange(self.nchan):
            channel = { \
                     'data':data[ich], \
                     'stepsize':self.stepsizes[ich], \
                     'units':self.chan_hdrs[ich]['units_text'], \
                     'chan_num':self.chan_hdrs[ich]['chan_order'], \
                     'description':self.chan_hdrs[ich]['comment_text'], \
                     }
            self.contents['channels'].append(channel)
        self.data = data

    def DumpTabDelimited(self, nmax=None, channels=None):
        fname = '%s.txt' % self.prefix
        fd = open(fname, 'w')
        if channels is None:
            channels = range(self.nchan)
        fd.write('%s' % self.contents['channels'][channels[0]]['description'])
        for ich in channels[1:]:
            fd.write('\t%s' % self.contents['channels'][ich]['description'])
        fd.write('\n')
        for n in xrange(self.max_buflen):
            fd.write('%g' % self.contents['channels'][channels[0]]['data'][n])
#            for ich in xrange(1,self.nchan):
            for ich in channels[1:]:
                if n < self.chan_hdrs[ich]['buf_length']:
                    fd.write('\t%g' % self.contents['channels'][ich]['data'][n])
                else:
                    fd.write('\t')
            fd.write('\n')
            if nmax and n > nmax:
                break
        fd.close()

    def DumpPickle(self):
        import cPickle
        fname = '%s.pickle' % self.prefix
        f = open(fname, 'w')
        pickler = cPickle.Pickler(f)
        pickler.dump(self.contents)
        f.close()
        print 'File written to %s' % fname

    def DumpMatlab(self):
        from pickle_to_mat import PickleMat
        pm = PickleMat(self.contents, prefix=self.prefix)
        pm.Process()
        if self.print_roadmap:
            roadmap_prefix = self.prefix
        else:
            roadmap_prefix = None
        print 'Data written to %s.mat' % self.prefix

        if self.print_roadmap:
            pm.PrintRoadmap(prefix=roadmap_prefix)

    def DumpRaw(self, filename=None, fd=None):
        """
        Write results to a binary file.
        if "filename" is supplied, data will be written to it. Otherwise,
        it will be written to the stream specified by "fd". A IOError
        exception is raised if neither fd or filename is supplied.
        """
        header = '%d\n' % self.nchan

#       First dump the strings, one per line so Matlab can easily figure it out.
        for ich in xrange(self.nchan):
            header += '%s\n' % self.contents['channels'][ich]['description']
        for ich in xrange(self.nchan):
            header += '%s\n' % self.contents['channels'][ich]['units']
        ascii_lgth = len(header) + 6
        if not fd:
            if filename:
                fd = open(filename, 'w')
            else:
                raise IOError('Either file descriptor or filename must ' + \
                              'be specifed in call to DumpResults')

#       Save the header length as the first 5 bytes.
        fd.write('%05d\n' % ascii_lgth)
        fd.write(header)

#       Write the  numerical data, one array at a time.
        stepsizes = zeros(self.nchan, float32)
        buflgths = zeros(self.nchan, int32)
        chan_nums = zeros(self.nchan, int32)
        for ich in xrange(self.nchan):
            buflgths[ich] = self.contents['channels'][ich]['data'].shape[0]
            stepsizes[ich] = self.contents['channels'][ich]['stepsize']
            chan_nums[ich] = self.contents['channels'][ich]['chan_num']
        fd.write(buflgths.tostring())
        fd.write(stepsizes.tostring())
        fd.write(chan_nums.tostring())

#       Now write the actual data.
        for ich in xrange(self.nchan):
            fd.write(self.contents['channels'][ich]['data'].\
                                        astype(float32).tostring())
        if filename:
            fd.close()

    def FormatLine(self, fields, fieldwidths=None, labels=False, gap=3):
        out = [gap*' ']
        for i in xrange(len(fields)):
            field = fields[i]
            pad = (fieldwidths[i] - len(field))*' '
            out.append('%s%s%s' % \
                (field, pad, gap*' '))
        out.append('\n')
        out1 = []
        if labels:
            for field in out:
                sout = ''
                for c in field:
                    if c.isalnum():
                        sout += '-'
                    else:
                        sout += c
                out1.append(sout)
        return ''.join(out + out1)
        
    def DumpSummary(self, fd=sys.stdout):
        """
        Dump a summary to the file object specified by fd.
        """

#       Create a pad so comments are equal length.
        maxlab = 0
        maxunits = 0
        for ich in xrange(self.nchan):
            if len(self.chan_hdrs[ich]['comment_text']) > maxlab:
                maxlab = len(self.chan_hdrs[ich]['comment_text'])
            if len(self.chan_hdrs[ich]['units_text']) > maxunits:
                maxunits = len(self.chan_hdrs[ich]['units_text'])
        labpad = []
        unitspad = []
        for ich in xrange(self.nchan):
            labpad.append((maxlab - len(self.chan_hdrs[ich]['comment_text'])+1)*' ')
            unitspad.append((maxunits - len(self.chan_hdrs[ich]['units_text'])+1)*' ')
        pad = ((maxlab-11)/2)*' '
        title = ['Channel', 'Description', 'Minimum' ,'Maximum' ,'Units' ,'Stepsize', 'Length', 'Duration']
        maxfws = []
        for field in title:
            maxfws.append(len(field))
        lines = []
        for ich in xrange(self.nchan):
            ms = int(self.chan_buflen[ich]*self.stepsizes[ich])
            sec = int(ms/1000.)
            min = sec/60
            sec -= min*60
            ms -= 1000*(min*60 + sec)
            line = ['%d' % self.chan_hdrs[ich]['chan_order'], \
                    '%s' % self.chan_hdrs[ich]['comment_text'], \
                    '%7.2f' % self.data[ich].min(), \
                    '%7.2f' % self.data[ich].max(),  \
                    '%s' % self.chan_hdrs[ich]['units_text'],
                    '%5.1f' % self.stepsizes[ich]*self.sample_spacing[ich], \
                    '%d' % self.chan_buflen[ich], \
                    '%d:%02d:%03d' % (min, sec, ms)]
            lines.append(line)
            for i in xrange(len(line)):
                field = line[i]
                if len(field) > maxfws[i]:
                    maxfws[i] = len(field)
        line = self.FormatLine(title,maxfws, labels=True)
        fd.write(line)
        for line in lines:
            fd.write(self.FormatLine(line, maxfws))
        buf_size = 0

    def SetupVarsampBlock(self):
        """
        Setup parameters required to read data with variable sampling
        across channels. This is done by determining the minimum length block
        that defines the sampling pattern, generating adresses from it.
        """
        maxpts = self.chan_buflen.max()
        vsds = unique(self.sample_spacing.tolist())
        if self.verbose:
            print 'Relative sampling intervals: ',self.sample_spacing.tolist()
        self.blocklen = prod(vsds)
        if self.blocklen < 4096:
            self.blocklen = self.blocklen*int(4096/self.blocklen)
        self.chan_width = []
        self.act_chan = []
        self.npseudo_chan = 0
        for ich in xrange(self.nchan):
            self.npseudo_chan += self.per_chan_type[ich]['size']/2
            self.chan_width.append(self.per_chan_type[ich]['size']/2)
            for i in xrange(self.per_chan_type[ich]['size']/2):
                self.act_chan.append(ich)
#       Create a mask where each chan is represented by one to four columns (
#       inti6 channels use one columns, float64 channels use 4)
        mask = zeros([self.blocklen, self.npseudo_chan], int)
        self.Nsamps = []
        for ich in xrange(self.npseudo_chan):
            act_chan = self.act_chan[ich]
            Nsamp = self.blocklen/self.sample_spacing[act_chan]
            self.Nsamps.append(Nsamp)
            idcs = self.sample_spacing[act_chan]*arange(Nsamp)
            mask[:,ich].put(idcs, ones(Nsamp))
        self.raw_blocklen = mask.sum()

#       Create an array of indices mapping the samples for each channel to
#       specific indices in the input block.
        mask = mask.reshape(prod(mask.shape))
        idcs_in =  reshape(mask*(mask.cumsum()), \
                                [self.blocklen, self.npseudo_chan])
        self.idcs_in = []
        jch = 0
        mask = mask.reshape([self.blocklen, self.npseudo_chan])
        for ich in xrange(self.nchan):
            chdat = idcs_in[:,jch:jch+self.chan_width[ich]].\
                                reshape(self.blocklen*self.chan_width[ich])
            self.idcs_in.append((chdat-1).take(nonzero(chdat)))
            jch += self.chan_width[ich]

#       Create corresonding output indices.
        self.idcs_out = []
        for ich in xrange(self.nchan):
            self.idcs_out.append(arange(self.Nsamps[ich]))

    def DownsampleData(self, final_stepsize):
        if final_stepsize % self.stepsize > .001:
            raise OSError(\
            'read_biopac::DownSampleData: Final stepsize ' + \
            '(%f) must be a multiple of native stepsize (%f)' % \
            (final_stepsize, self.stepsize))
        else:
            data_dsamp = []
            npts = []
            for ich in xrange(self.nchan):
                stepsize = self.stepsize*self.sample_spacing[ich]
                dsamp_factor = int(final_stepsize/stepsize)
                npts = self.npts/dsamp_factor
                idcs = 2*arange(npts).astype(int)
                chdata = zeros([self.npts/dsamp_factor], float)
                for i in xrange(dsamp_factor):
                    chdata[:] += self.data[ich][:].take(idcs+i)
                data_dsamp.append(chdata)
                npts.append(npts)
#            self.data = data/float(dsamp_factor)

        return npts, data_dsamp

    def PlotData(self, xmin=None, npts=None, channels=None, width=8, height=10):
        from plotlib import ScatterPlot
        hdr = self.contents['hdr']
        chan_hdrs = hdr['channels']
        if npts is None:
            N = 2000
        else:
            N = npts
        if channels is None:
            chans = range(self.nchan)
        else:
            chans = channels
        sp = ScatterPlot(len(chans), 1, self.filename, width=width, height=height)
        for ch in chans:
#            sample_spacing = chan_hdrs[ch]['var_sample_divider']
            sample_spacing = self.contents['channels'][ch]['stepsize']
            buf_len = chan_hdrs[ch]['buf_length']
            if xmin is None:
                M = buf_len/2
            else:
                M = xmin
            x = .001*(M + sample_spacing*arange(N))
            y = self.contents['channels'][ch]['data'][M:M+N]
            if ch == self.nchan-1:
                xlab = 'Time in sec'
            else:
                xlab = ''
            sp.AddPlot( \
                    x, \
                    y, \
                    x_label=xlab, \
                    y_label='%s\n%s' % (chan_hdrs[ch]['comment_text'], \
                                        chan_hdrs[ch]['units_text']), \
                    overplot=False, \
                    markersize=4, \
                    xgrid_lines = True, \
                    lineonly=True)
#                    subtitle=chan_hdrs[ch]['comment_text'], \
        sp.Show()

def read_biopac():
    rh = ReadBiopac()
    rh.ParseOptions()
    rh.ReadHeader()
    rawdata = rh.ReadData()
#    if rawdata is not None:
#        data = rh.DownsampleData(.025)

#   Write a summary to the screen.
    rh.DumpSummary()

    if rh.format == 'matlab':
        rh.DumpMatlab()
    elif rh.format == 'pickle':
        rh.DumpPickle()
    elif rh.format == 'text':
        rh.DumpTabDelimited()
    else:
        raise RuntimeError( \
        'Invalid output format: %s' % rh.format)

    if rh.plot:
        rh.PlotData()

if __name__ == '__main__':
    try:
        read_biopac()
    except (RuntimeError, IOError, ValueError, AttributeError), errmsg:
        sys.stderr.write('\n%s\n%s\n' % (errmsg, except_msg()))

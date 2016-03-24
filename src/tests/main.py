import os
import unittest
from bioread import *
import numpy

SAMPLE_PATH = os.path.join(os.path.dirname(os.path.abspath(__file__)), "samples")
def sample(name):
    return os.path.join(SAMPLE_PATH, name)


class Mock(object):
    def __init__(self, **kwargs):
        self.__dict__.update(kwargs)


class AcqReaderTestCase(unittest.TestCase):
    """
    These are unit-style tests that check various bits and bobs
    inside the AcqReader class
    """

    def producesSampleIndexes(self, frequencies, expected):
        a = AcqReader("no_file_needed.acq")
        channels = [ Mock(frequency_divider=f) for f in frequencies ]
        indexes = a._AcqReader__stream_sample_indexes(channels)
        self.assertEqual(list(indexes), expected)

    def testSimpleSampleIndexes(self):
        self.producesSampleIndexes([1,2,4], [0,1,2,0,0,1,0])
        self.producesSampleIndexes([1,4,2], [0,1,2,0,0,2,0])
        self.producesSampleIndexes([1,4,4], [0,1,2,0,0,0])
        self.producesSampleIndexes([2,4],   [0,1,0])
        self.producesSampleIndexes([1,3],   [0,1,0,0])
        self.producesSampleIndexes([1,3,4], [0,1,2,0,0,0,1,0,2,0,0,1,0,0,2,0,1,0,0])


class IntegrationTestCase(unittest.TestCase):
    """
    These are integration-style tests that pass .acq files
    to bioread's top level `read` function and test
    that the results pass certain standards.
    """

    def sanityCheck(self):
        """
        Run general sanity checks to make sure numbers line up right.

        These are not actually true in many cases right now.
        """
        self.sanityCheckHeaderAndChannelCounts()

        # NOTE: These fail a lot right now! You can comment them out for easier fun-having.
        self.sanityCheckTotalSampleCount()
        self.sanityCheckPointCounts()

    def sanityCheckTotalSampleCount(self):
        sample_sum = sum([len(ch.raw_data) for ch in self.data.channels])
        self.assertEqual(self.reader.total_samples, sample_sum)

    def sanityCheckPointCounts(self):
        headers = [ch.point_count for ch in self.data.channel_headers]
        channels = [ch.point_count for ch in self.data.channels]
        self.assertEqual(headers, channels)

    def sanityCheckHeaderAndChannelCounts(self):
        self.assertEqual(len(self.data.channel_headers), len(self.data.channels))


    def channelLength(self, channel_index, length):
        """
        Helper to test if a given channel has a specific length
        in both the headers and the channel data.
        """
        self.assertEqual(self.data.channel_headers[channel_index].point_count, length)
        self.assertEqual(len(self.data.channels[channel_index].raw_data), length)


    def load(self, name):
        path = sample(name + ".acq")
        with open(path, 'rb') as f:
            self.reader = AcqReader(f)
            self.data = self.reader.read()
            self.sanityCheck()
            return self.data



    # Sample files were generated with AcqKnowledge 4.2.0
    # c1 prefix: a single normal "eyeblink" channel
    # c2 prefix: 2 channels; "eyeblink" channel, smoothed calculated channel
    # c3 prefix: 3 channels; "eyeblink" channel, "zygo" channel, "corrug" channel, rates are 500 250 100 (?) unless specified

    def check_c1(self):
        """
        Shared checks for the c1* files
        """
        d = self.data
        self.assertEqual(d.channel_headers[0].frequency_divider, 1)
        self.assertEqual(d.channel_headers[0].name, "Eyeblink - ERS100C")
        self.assertEqual(d.channel_headers[0].units, "mV")

        self.channelLength(0, 2001)

    def test_c1(self):
        """
        Simple case with a single channel
        """
        self.load("c1")
        self.check_c1()

    def test_c1_w3(self):
        """
        Simple case with a single channel, version 3
        """
        self.load("c1_w3")
        self.check_c1()


    def check_c2(self):
        """
        Shared checks for the c2* files
        """
        d = self.data
        # Channel 1 is the "eyeblink" data
        c1h = d.channel_headers[0]
        # Channel 2 is the calculated channel based on channel 1
        c2h = d.channel_headers[1]

        self.assertEqual(c1h.name, "Eyeblink - ERS100C")
        self.assertEqual(c1h.frequency_divider, 1)

        self.assertEqual(c2h.name, "C0 - Smoothing")
        self.assertEqual(c1h.frequency_divider, 1)

        self.channelLength(0, 2001)
        self.channelLength(1, 2001)

    def test_c2(self):
        """
        Simple case with smoothed calculated channel
        """
        d = self.load("c2")
        self.check_c2()

    def test_c2_w3(self):
        """
        Simple case with smoothed calculated channel, version 3
        """
        d = self.load("c2_w3")
        self.check_c2()


    def check_c3_rates1(self):
        """
        Shared checks for the c3_rates1* files
        """
        d = self.data
        c1h = d.channel_headers[0]
        c2h = d.channel_headers[1]
        c3h = d.channel_headers[2]

        self.assertEqual(c1h.name, "Eyeblink - ERS100C")
        self.assertEqual(c1h.frequency_divider, 1)

        self.assertEqual(c2h.name, "Zygo - ERS100C")
        self.assertEqual(c2h.frequency_divider, 16)

        self.assertEqual(c3h.name, "Corrug - ERS100C")
        self.assertEqual(c3h.frequency_divider, 2)

        self.channelLength(0, 2001)
        self.channelLength(1, 2001)

    def test_c3_rates1(self):
        """
        Different sample rates on each channel,
        including one at the base rate:
        1000, 62.5, 500
        """
        self.load("c3_rates1")
        self.check_c3_rates1()

    def test_c3_rates1_w3(self):
        """
        Different sample rates on each channel,
        including one at the base rate:
        1000, 62.5, 500, version 3
        """
        self.load("c3_rates1_w3")
        self.check_c3_rates1()



    def check_c3_rates2(self):
        """
        Shared checks for the c3_rates2* files
        """
        d = self.data
        c1h = d.channel_headers[0]
        c2h = d.channel_headers[1]
        c3h = d.channel_headers[2]

        self.assertEqual(c1h.name, "Eyeblink - ERS100C")
        self.assertEqual(c1h.frequency_divider, 2)

        self.assertEqual(c2h.name, "Zygo - ERS100C")
        self.assertEqual(c2h.frequency_divider, 4)

        self.assertEqual(c3h.name, "Corrug - ERS100C")
        self.assertEqual(c3h.frequency_divider, 8)

    def check_c3_rates2_cut(self):
        """
        Check the lengths for the simple cut version
        where a few frames are chopped off the first channel
        """
        self.channelLength(0, 998)
        self.channelLength(1, 500)
        self.channelLength(2, 250)

    def check_c3_rates2_bigcut(self):
        """
        Check the lengths for the big cut version
        where a ton of frames are chopped off the first channel
        """
        self.channelLength(0, 958)
        self.channelLength(1, 500)
        self.channelLength(2, 250)


    def test_c3_rates2(self):
        """
        None of the rates of any of the channels exactly
        matches the base rate
        """
        self.load("c3_rates2")
        self.check_c3_rates2()

    def test_c3_rates2_w3(self):
        """
        None of the rates of any of the channels exactly
        matches the base rate, version 3
        """
        self.load("c3_rates2_w3")
        self.check_c3_rates2()

    def test_c3_rates2_cut(self):
        """
        None of the rates of any of the channels exactly matches
        the base rate and the tail end of channel 1 is cut
        """
        self.load("c3_rates2_cut")
        self.check_c3_rates2()
        self.check_c3_rates2_cut()

    def test_c3_rates2_cut_w3(self):
        """
        None of the rates of any of the channels exactly matches
        the base rate and the tail end of channel 1 is cut, version 3
        """
        self.load("c3_rates2_cut_w3")
        self.check_c3_rates2()
        self.check_c3_rates2_cut()

    def test_c3_rates2_bigcut(self):
        """
        None of the rates of any of the channels exactly matches
        the base rate and the tail end of channel 1 is REALLY cut
        """
        self.load("c3_rates2_bigcut")
        self.check_c3_rates2()
        self.check_c3_rates2_bigcut()

    def test_c3_rates2_bigcut_w3(self):
        """
        None of the rates of any of the channels exactly matches
        the base rate and the tail end of channel 1 is REALLY cut, version 3
        """
        self.load("c3_rates2_bigcut_w3")
        self.check_c3_rates2()
        self.check_c3_rates2_bigcut()

    
    def check_c3_rates3_clearall(self):
        """
        Shared checks for the c3_rates3_clearall* files,
        where all the channels
        """
        d = self.data
        c1h = d.channel_headers[0]
        c2h = d.channel_headers[1]
        c3h = d.channel_headers[2]

        self.assertEqual(c1h.name, "Eyeblink - ERS100C")
        self.assertEqual(c1h.frequency_divider, 1)

        self.assertEqual(c2h.name, "Zygo - ERS100C")
        self.assertEqual(c2h.frequency_divider, 1)

        self.assertEqual(c3h.name, "Corrug - ERS100C")
        self.assertEqual(c3h.frequency_divider, 64)

        self.channelLength(0, 1789)
        self.channelLength(1, 1789)
        self.channelLength(2, 28)

    def test_c3_rates3_clearall(self):
        """
        Different sample rates on last channel and 
        'clear all' is used to break off end of file:
        1000, 1000, 15
        """
        self.load("c3_rates3_clearall")
        self.check_c3_rates3_clearall()

    def test_c3_rates3_clearall_w3(self):
        """
        Different sample rates on last channel and
        'clear all' is used to break off end of file:
        1000, 1000, 15, version 3
        """
        self.load("c3_rates3_clearall_w3")
        self.check_c3_rates3_clearall()


    def test_c1_compressed(self):
        """
        Confirm we can read compressed channel data fine,
        single channel
        """
        self.load("c1_compressed")

    def test_c2_compressed(self):
        """
        Confirm we can read compressed channel data fine,
        multiple channels
        """
        self.load("c2_compressed")

    def test_c3_rates3_clearall_compressed(self):
        """
        Confirm we can read crazy 'clear all' and compressed 
        channel data fine AT THE SAME TIME
        """
        self.load("c3_rates3_clearall_compressed")
        self.check_c3_rates3()


    def check_c1_marker(self):
        """
        Shared checks for the c1_marker* files
        """
        # TODO: Actual test of markers
        self.assertEqual(len(self.data.markers), 3)

    @unittest.skip("Not implemented")
    def test_c1_marker(self):
        """
        Confirm we can read markers
        """
        self.load("c1_marker")
        self.check_c1_marker()

    @unittest.skip("Not implemented")
    def test_c1_marker_w3(self):
        """
        Confirm we can read markers from version 3
        
        Note that markers are much more limited in v3; may not make sense
        to have a shared checking function.
        """
        self.load("c1_marker_w3")
        self.check_c1_marker()


    def test_c1_fast(self):
        """
        Confirm we can read extremely fast sample rate
        """
        self.load("c1_fast")
        self.channelLength(0, 2002)

    def test_c1_slow(self):
        """
        Confirm we can read extremely slow sample rate
        """
        self.load("c1_slow")
        self.channelLength(0, 101)



    def testBrokenFile(self):
        """
        Test loading strange recording of Dan during
        MIDUSREF audio testing that causes numpy reshape error

        (I broke it with my MIND)
        """
        self.load("broken")

    def testHuge(self):
        """
        Giant file that exposes issues with wacky channel lengths
        """
        self.load("huge")


unitSuite = unittest.makeSuite(AcqReaderTestCase,'test')
integrationSuite = unittest.makeSuite(IntegrationTestCase,'test')

if __name__ == "__main__":
    unittest.main()


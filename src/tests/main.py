
import os
import unittest
from bioread import *

SAMPLE_PATH = os.path.join(os.path.dirname(os.path.abspath(__file__)), "samples")
def sample(name):
    return os.path.join(SAMPLE_PATH, name)




class IntegrationTestCase(unittest.TestCase):
    """
    These are integration-style tests that pass .acq files
    to bioread's top level `read_file` function and test
    that the results pass certain standards.
    """

    def sanityCheck(self, data):
        self.sanityCheckHeaderAndChannelLengths(data)

    def sanityCheckHeaderAndChannelLengths(self, data):
        headers = [ch.point_count for ch in data.channel_headers]
        channels = [ch.point_count for ch in data.channels]

        self.assertEqual(len(headers), len(channels))

        # This is not actually true in the general case
        #self.assertEqual(headers, channels)

    def load(self, name):
        self.data = read_file(sample(name + ".acq"))
        self.sanityCheck(self.data)

    # Cases were generated with AcqKnowledge 4.2.0
    # c1 prefix: a single normal "eyeblink" channel
    # c2 prefix: 2 channels; "eyeblink" channel, smoothed calculated channel
    # c3 prefix: 3 channels; "eyeblink" channel, "zygo" channel, "corrug" channel, rates are 500 250 100 (?) unless specified

    def test_c1(self):
        """Simple case with a single channel"""
        self.load("c1")

    def test_c1_w3(self):
        """Simple case with a single channel, version 3"""
        self.load("c1_w3")

    def test_c2(self):
        """Simple case with smoothed calculated channel"""
        self.load("c2")

    def test_c2_w3(self):
        """Simple case with smoothed calculated channel, version 3"""
        self.load("c2_w3")


    def test_c3_rates1(self):
        """Different sample rates on each channel, including one at the base rate: 1000, 62.5, 500"""
        self.load("c3_rates1")

    def test_c3_rates1_w3(self):
        """Different sample rates on each channel, including one at the base rate: 1000, 62.5, 500, version 3"""
        self.load("c3_rates1_w3")

    def test_c3_rates2(self):
        """None of the rates of any of the channels exactly matches the base rate"""
        self.load("c3_rates2")

    def test_c3_rates2_w3(self):
        """None of the rates of any of the channels exactly matches the base rate, version 3"""
        self.load("c3_rates2_w3")

    def test_c3_rates2_cut(self):
        """None of the rates of any of the channels exactly matches the base rate and the tail end of channel 1 is cut"""
        self.load("c3_rates2_cut")

    def test_c3_rates2_cut_w3(self):
        """None of the rates of any of the channels exactly matches the base rate and the tail end of channel 1 is cut, version 3"""
        self.load("c3_rates2_cut_w3")

    def test_c3_rates2_bigcut(self):
        """None of the rates of any of the channels exactly matches the base rate and the tail end of channel 1 is REALLY cut"""
        self.load("c3_rates2_bigcut")

    def test_c3_rates2_bigcut_w3(self):
        """None of the rates of any of the channels exactly matches the base rate and the tail end of channel 1 is REALLY cut, version 3"""
        self.load("c3_rates2_bigcut_w3")

    
    def test_c3_rates3_clearall(self):
        """Different sample rates on last channel and 'clear all' is used to break off end of file: 1000, 1000, 15"""
        self.load("c3_rates3_clearall")

    def test_c3_rates3_clearall_w3(self):
        """Different sample rates on last channel and 'clear all' is used to break off end of file: 1000, 1000, 15, version 3"""
        self.load("c3_rates3_clearall_w3")

    def test_c1_compressed(self):
        """Confirm we can read compressed channel data fine, single channel"""
        self.load("c1_compressed")

    def test_c2_compressed(self):
        """Confirm we can read compressed channel data fine, multiple channels"""
        self.load("c2_compressed")

    def test_c3_rates3_clearall_compressed(self):
        """Confirm we can read crazy 'clear all' and compressed channel data fine"""
        self.load("c3_rates3_clearall_compressed")


    def test_c1_marker(self):
        """Confirm we can read markers"""
        self.load("c1_marker")

    def test_c1_marker_w3(self):
        """Confirm we can read markers from version 3"""
        self.load("c1_marker_w3")


    def test_c1_fast(self):
        """Confirm we can read extremely high sample rate"""
        self.load("c1_fast")

    def test_c1_slow(self):
        """Confirm we can read extremely slow sample rate"""
        self.load("c1_slow")



    def testBroken(self):
        """Strange file that causes numpy error"""
        self.load("broken")

    #@unittest.skip("Currently broken")
    def testHuge(self):
        """Giant file that exposes issues with wacky channel lengths"""
        self.load("huge")


simpleSuite = unittest.makeSuite(IntegrationTestCase,'test')

if __name__ == "__main__":
    unittest.main()


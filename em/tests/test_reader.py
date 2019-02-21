import unittest
import reader as rd

class TestReader(unittest.TestCase):


    def test_init(self):
        with self.assertRaises(IOError):
            myreader = rd.Reader('somefile.map')
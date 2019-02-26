import unittest
import reader as rd


class TestReader(unittest.TestCase):


    def test_init(self):
        with self.assertRaises(IOError):
            myreader = rd.Reader('somefile.map')


    def test_read(self):
        myreader = rd.Reader('tests/EMD-2677.map')
        otherreader = rd.Reader('tests/EMD-2677.map')
        mymolecule = myreader.read()
        othermolecule = otherreader.read()
        self.assertEqual(mymolecule, othermolecule)
 
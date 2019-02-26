import unittest
import reader as rd
import molecule as m


class TestReader(unittest.TestCase):


    def test_init(self):
        with self.assertRaises(IOError):
            myreader = rd.Reader('somefile.map')


    def test_read(self):
        myreader = rd.Reader('tests/EMD-2677.map')
        mym = myreader.read()
        
        othermolecule = m.Molecule(mym.rawHeader, mym.data(), mym.shape(), mym.start_point(), mym.grid_size(), mym.cell_dim(), mym.density_range(), mym.origin())
        self.assertEqual(mym, othermolecule)
 
    def test_bigEndianness(self):
        readerBig = rd.Reader('tests/EMD-2627_big.map')
        self.assertTrue(readerBig.is_endianness_reversed)

    def test_littleEndienness(self):
        readerLittle = rd.Reader('tests/EMD-2627_little.map')
        self.assertTrue(~readerLittle.is_endianness_reversed)
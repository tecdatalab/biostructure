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
 
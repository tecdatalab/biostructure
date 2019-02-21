import unittest
import reader as rd

class TestReader(unittest.TestCase):


    def test_init(self):
        with self.assertRaises(IOError):
            myreader = rd.Reader('somefile.map')


    def test_read(self):
        myreader = rd.Reader('tests/EMD-sample.map')
        mymolecule = myreader.read()
        
        self.assertEqual(mymolecule.shape,(140,140,140))
        self.assertEqual(mymolecule.start_point, (0,0,0))
        self.assertEqual(mymolecule.grid_size, (140,140,140))
        self.assertEqual(mymolecule.cell_dim, (246.4, 246.4, 246.4))
        self.assertEqual(mymolecule.density_range, (-0.24971596896648407,0.4126806855201721,5.458949090098031e-05))
        self.assertEqual(mymolecule.origin, (0., 0., 0.))
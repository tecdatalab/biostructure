
import reader
import visualizer
import processing
import molecule

import configparser
import argparse

import numpy as np



def main():
    #config = configparser.ConfigParser()
    #parser = argparse.ArgumentParser(description='Molecule Visualization and segmentation tool')

    #Open file
    #mapReader.open("../maps/1010/EMD-1010.map") #level 2.2479
    # #mapReader.open("../maps/EMD-2596.map") #level 0.1462
    #mapReader.open("../maps/1010/EMD-1010.map") #level 23
    #mapReader.open("../maps/1364/EMD-1364.map") #level 39.0086
    #mapReader.open("../maps/5017/EMD-5017.map") #level 17.3347


    myMolecule = molecule.Molecule("../maps/1010/EMD-1010.map", 7, [1])
    myMolecule.generateSegments(3,1)

    # Create visualizer with a map surface threshold level
    # Otherwise use otsu threshold

    #v1= visualizer.Visualizer(myMolecule.getEmMap().data(), 7, myMolecule.getSegments()[0])
    #
    # Watershed 


    


    #precomputed_descriptors = [6.00485,9.20974e-16,8.0491,1.96306,2.0434,1.52559,5.27135,4.58634,1.6972,3.12813,1.8177,1.5765,1.47143,6.2866,2.85751,1.24786,3.25437,2.00659,2.18337,1.35591,0.369011,5.01783,3.68149,1.2648,1.21123,2.71391,1.6935,2.89683,1.71589,1.11724,0.154151,2.49272,3.25498,1.80037,1.1011,1.08982,1.79274,1.74014,3.02078,2.16811,1.11224,0.989201,0.800181,2.05593,2.20992,2.07055,1.13324,1.01053,0.934,1.01314,1.89837,2.64202,1.93436,1.58645,0.767435,0.874263,0.939151,2.07539,2.01024,1.94173,1.05544,1.3689,0.698616,0.809511,0.870718,1.45429,2.21419,1.85029,1.4274,1.04587,0.557376,0.760452,0.946489,1.055,1.60476,1.50551,1.43065,1.2094,1.10033,0.501099,0.704349,0.77951,1.07687,1.76099,1.56281,1.43962,0.984525,0.773406,0.415912,0.654815,0.259077,0.811758,1.48604,1.14667,1.42665,1.6946,0.9882,0.831599,0.384764,0.610353,0.910067,0.884391,1.418,1.53506,0.853752,1.12984,0.71894,0.743243,0.343151,0.565689,0.551233,0.463538,0.908726,1.17644,1.31264,0.904315,1.30646,0.870354,0.619409,0.324997,0.526729]
    #import utils._zernike as z 
    #zd = z.computeDescriptors(myModel.getData()[0])

    #print(np.linalg.norm(precomputed_descriptors - zd/10))
 


    #v1.add_structure("../maps/pdb6acf.ent")
    #v2.add_structure("../maps/pdb6acf.ent")
    #v3.add_structure("../maps/pdb6acf.ent")

    #v1.show()
    #v2.show()
    #v3.show()


    # add corresponding atomic structure
    #v.add_structure("../maps/1010/pdb1mi6.ent")
    #v.map_structure_to_domain("../maps/1010/A-1GQE.aligned.pdb")
    #v.map_structure_to_domain("../maps/1010/B-1GQE.aligned.pdb")
    #v.map_structure_to_domain("../maps/1010/C-1GQE.aligned.pdb")
    #v.map_structure_to_domain("../maps/1010/D-1GQE.aligned.pdb")

    #v.add_structure("../maps/1364/pdb1pn6.ent")
    #v.map_structure_to_domain("../maps/1364/A-1FNM.aligned.pdb")
    #v.map_structure_to_domain("../maps/1364/B-1FNM.aligned.pdb")
    #v.map_structure_to_domain("../maps/1364/C-1FNM.aligned.pdb")
    #v.map_structure_to_domain("../maps/1364/D-1FNM.aligned.pdb")
    #v.map_structure_to_domain("../maps/1364/E-1FNM.aligned.pdb")

    #v.add_structure("../maps/5017/pdb3dny.ent")
    #v.map_structure_to_domain("../maps/5017/A-1N0U.aligned.pdb")
    #v.map_structure_to_domain("../maps/5017/B-1N0U.aligned.pdb")
    #v.map_structure_to_domain("../maps/5017/C-1N0U.aligned.pdb")

    #v.show_atom_matching()
    #v.show()



if __name__ == "__main__":
    main()
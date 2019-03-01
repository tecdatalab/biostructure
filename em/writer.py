import numpy as np 
import struct


class Writer():

    #Header size in bytes
    HEADER_SIZE = 1024

    # Write density map as ".map" format. Default write mode is set to mode 2. 
    def write(self,filename, molecule):
        #First generate buffer from array data in mode 2
        data = bytearray()
        for voxel in np.nditer(molecule.data()):
            data+=struct.pack('f',voxel)        
        try:
            with open(filename, 'wb') as output:
                header = bytearray(molecule.rawHeader)
                header[76:88] = struct.pack('<fff', *molecule.density_range())
                header[196:208] = struct.pack('<fff', *molecule.origin())

                output.write(header)
                output.write(data)
        except Exception as err:
            print("Could not write file, error ", err)



'''
import reader
filename = "tests/EMD-2677.map"
myreader = reader.Reader(filename)
D = myreader.read()
w = Writer()
w.write("emd_2677_2.map", D)
'''
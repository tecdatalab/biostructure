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
                header = bytearray(HEADER_SIZE)
                header[0:12] = struct.pack('<III', *molecule.shape())
                header[12:16] = struct.pack('<I', 2)
                header[16:28] = struct.pack('<III', *molecule.start_point())
                header[28:40] = struct.pack('<III' *molecule.grid_size())
                header[40:52] = struct.pack('<fff' *molecule.cell_dim())
                header[76:88] = struct.pack('<fff' *molecule.density_range())
                header[196:208] = struct.pack('<fff' *molecule.origin())

                output.write(molecule.rawHeader)
                output.write(data)
        except IOError as err:
            print("Could not write file, error ", err)




import reader
filename = "tests/EMD-2677.map"
myreader = reader.Reader(filename)
D = myreader.read()
w = Writer()
w.write("emd_2677_2.map", D)

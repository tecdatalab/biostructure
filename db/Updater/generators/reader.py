'''
   
  Based on https://github.com/tecdatalab/legacy/blob/master/core/em/mrc.h

  MRC files are a standard format for electron density maps.
  The first 1024 bytes represent the header followed by the actual
  densities.
  
  56 4-byte words compose the standard fields followed by
  custom text labels that can be added (e.g. to express symmetry)
 
  This is a specification excerpt that describes those fields
  (from http://ami.scripps.edu/software/mrctools/mrc_specification.php)
  (Here's another description http://bio3d.colorado.edu/imod/doc/mrc_format.txt)

  1	NX       number of columns (fastest changing in map)
  2	NY       number of rows   
  3	NZ       number of sections (slowest changing in map)
  4	MODE     data type :
       0        image : signed 8-bit bytes range -128 to 127 (char)
       1        image : 16-bit halfwords (signed short)
       2        image : 32-bit reals (float)
       3        transform : complex 16-bit integers (2 * signed short)
       4        transform : complex 32-bit reals (2 * float)
       6        image : unsigned 16-bit range 0 to 65535 (unsigned short)
  5	NXSTART number of first column in map (Default = 0)
  6	NYSTART number of first row in map
  7	NZSTART number of first section in map
  8	MX       number of intervals along X
  9	MY       number of intervals along Y
  10	MZ       number of intervals along Z
  11-13	CELLA    cell dimensions in angstroms (xlen, ylen, zlen)
  14-16	CELLB    cell angles in degrees
  17	MAPC     axis corresp to cols (1,2,3 for X,Y,Z)
  18	MAPR     axis corresp to rows (1,2,3 for X,Y,Z)
  19	MAPS     axis corresp to sections (1,2,3 for X,Y,Z)
  20	DMIN     minimum density value
  21	DMAX     maximum density value
  22	DMEAN    mean density value
  23	ISPG     space group number 0 or 1 (default=0)
  24	NSYMBT   number of bytes used for symmetry data (0 or 80)
  25-48(49??)	EXTRA    extra space used for anything   - 0 by default
  Note: There appears to be a discrepancy, it looks like the origin is 49-51
  50-52	ORIGIN   origin in X,Y,Z used for transforms
  53	MAP      character string 'MAP ' to identify file type
  54	MACHST   machine stamp
  55	RMS      rms deviation of map from mean density
  56	NLABL    number of labels being used
  57-256	LABEL(20,10) 10 80-character text labels
 
 '''

import struct
import array
import numpy as np
from generators.molecule import Molecule


class Reader():
    #Header size in bytes
    HEADER_SIZE = 1024
    
    def open(self, filename):
        try:
            with open(filename, 'rb') as mrc_file:
                mrc_buffer = mrc_file.read()
        except IOError as err:
            print("Could not open MRC file", err)
            raise
        self.mrc_header = mrc_buffer[0:self.HEADER_SIZE]
        self.is_endianness_reversed = self.test_and_set_endianness()
        self.mrc_data = mrc_buffer[self.HEADER_SIZE:]
        
    def read(self):
        #Read number of collumns(x), rows(y) and sections(z)
        nx = self.read_header_word(self.mrc_header[0:4])
        ny = self.read_header_word(self.mrc_header[4:8])
        nz = self.read_header_word(self.mrc_header[8:12])
        self.mode = self.read_header_word(self.mrc_header[12:16])
        if (self.mode > 2):
            print("Only map modes 0, 1 and 2 are supported. Mode %d found" % self.mode)
            raise ValueError
        #Read start x, y, z positions
        nxstart = self.read_header_word(self.mrc_header[16:20])
        nystart = self.read_header_word(self.mrc_header[20:24])
        nzstart = self.read_header_word(self.mrc_header[24:28])
        #Read grid size
        mx = self.read_header_word(self.mrc_header[28:32])
        my = self.read_header_word(self.mrc_header[32:36])
        mz = self.read_header_word(self.mrc_header[36:40])
        #Read cell dimensions
        xlen = self.read_header_float(self.mrc_header[40:44])
        ylen = self.read_header_float(self.mrc_header[44:48])
        zlen = self.read_header_float(self.mrc_header[48:52])
        #Skip 3 cell angle words and the 3 words corresponding to axes
        #Read density ranges
        dmin = self.read_header_float(self.mrc_header[76:80])
        dmax = self.read_header_float(self.mrc_header[80:84])
        dmean = self.read_header_float(self.mrc_header[84:88])
        #Read origin
        xorigin = self.read_header_float(self.mrc_header[196:200])
        yorigin = self.read_header_float(self.mrc_header[200:204])
        zorigin = self.read_header_float(self.mrc_header[204:208])
        #Read densities
        densities = self.read_densities((nz, ny, nx))
        #Generate Molecule object with parameters
        return Molecule(self.mrc_header, densities,  (nz, ny, nx), (nzstart, nystart, nxstart), (mz, my, mx), (zlen, ylen, xlen), (dmin, dmax, dmean), (zorigin, yorigin, xorigin))
            
    def test_and_set_endianness(self):
        regular_nx = struct.unpack('<I', self.mrc_header[0:4])
        reversed_nx = struct.unpack('>I', self.mrc_header[0:4])
        return reversed_nx < regular_nx

    def read_header_word(self,buffer):
        if self.is_endianness_reversed:
            return struct.unpack('>I', buffer)[0]
        else: 
            return struct.unpack('<I', buffer)[0]

    def read_header_float(self,buffer):
        if self.is_endianness_reversed:
            return struct.unpack('>f', buffer)[0]
        else: 
            return struct.unpack('<f', buffer)[0]


    def read_densities(self, shape):
        if self.is_endianness_reversed:
            endianness = '>'
        else:
            endianness = '<'
        if self.mode == 0:
            dt = np.dtype(np.int8)
        if self.mode == 1:
            dt = np.dtype(np.int16)
        if self.mode == 2:
            dt = np.dtype(np.float32)
        dt = dt.newbyteorder(endianness)
        density_array =  np.frombuffer(self.mrc_data, dtype=dt)
        data_array = density_array.reshape(shape).transpose(2,1,0).astype(float)
        return data_array

'''
filename = "../../emd_2847.map"
myreader = Reader(filename)
D = myreader.read()
myreader.write("emd_2847_2.map", D)
'''
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
import os.path

from em import emMap
import mrcfile

class Reader():
    
    def open(self, filepath):
        try:
            map_object = mrcfile.open(filepath)
        except ValueError as err:
            print("[ERROR]: Could not open MRC file", err)
            raise err
        self.map = map_object
        self.header = self.map.header.flat[0]
        self.filename = os.path.splitext(os.path.basename(filepath))[0]
        
    def read(self):
        #Read number of collumns(x), rows(y) and sections(z)
        nx = self.header[0]
        ny = self.header[1]
        nz = self.header[2]
        self.mode = self.header[3]
        #Read start x, y, z positions
        nxstart = self.header[4]
        nystart = self.header[5]
        nzstart = self.header[6]
        #Read grid size
        mx = self.header[7]
        my = self.header[8]
        mz = self.header[9]
        #Read cell dimensions
        (xlen, ylen, zlen) = self.header[10]
        #Skip 3 cell angle words and the 3 words corresponding to axes
        #Read density ranges
        (dmin, dmax, dmean) = self.header[11]
        #axis order
        x = self.header[12]
        y = self.header[13]
        z = self.header[14]
        #Read origin
        (xorigin, yorigin, zorigin) = self.header[24]
        #Read densities
        densities = np.copy(self.map.data)
        #Generate Molecule object with parameters
        return emMap.EMMap(self.header, densities,  (nz, ny, nx), (nzstart, nystart, nxstart), (mz, my, mx), (zlen, ylen, xlen), (dmin, dmax, dmean), (zorigin, yorigin, xorigin), (x,y,z), self.filename)
            
    def save(self, filename, data):
        with mrcfile.new(filename, overwrite=True) as mrc:
            mrc.set_data(data)       

'''
filename = "../../emd_2847.map"
myreader = Reader(filename)
D = myreader.read()
myreader.write("emd_2847_2.map", D)
'''

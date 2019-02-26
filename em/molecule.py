import numpy as np

class Molecule():

    def __init__(self, rawHeader, data, size, start, grid_size, cell_dim, density_ranges, origin):
        self.rawHeader = rawHeader
        self.array = data
        (self.ny, self.nx, self.nz) = size
        (self.nystart, self.nxstart, self.nzstart) = start
        (self.my, self.mx, self.mz) = grid_size
        (self.ylen, self.xlen, self.zlen) = cell_dim
        (self.dmin, self.dmax, self.dmean) = density_ranges
        (self.yorigin, self.xorigin, self.zorigin) = origin


    def __eq__(self, other):
        if isinstance(other, self.__class__):
            hasSameData = ~(~np.isclose(other.array, self.array)).sum()
            return hasSameData
        return False
    
    def data(self):
        return self.array

    def shape(self):
        return (self.ny, self.nx, self.nz)

    def start_point(self):
        return (self.nystart, self.nxstart, self.nzstart)

    def grid_size(self):
        return (self.my, self.mx, self.mz)

    def cell_dim(self):
        return (self.ylen, self.xlen, self.zlen)

    def density_range(self):
        return (self.dmin, self.dmax, self.dmean)

    def origin(self):
        return (self.yorigin, self.xorigin, self.zorigin)



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
        (self.zorigin, self.yorigin, self.xorigin) = origin


    def __eq__(self, other):
        if isinstance(other, self.__class__):
            hasSameData = ~(~np.isclose(other.array, self.array)).sum()
            return hasSameData
        return False

    def set_origin(self, new_origin):
        (self.zorigin, self.yorigin, self.xorigin) = new_origin
    
    def data(self):
        return self.array

    def shape(self):
        return (self.nz, self.ny, self.nx)

    def start_point(self):
        return (self.nzstart, self.nystart, self.nxstart)

    def grid_size(self):
        return (self.mz, self.my, self.mx)

    def cell_dim(self):
        return (self.zlen, self.ylen, self.xlen)

    def density_range(self):
        return (self.dmin, self.dmax, self.dmean)

    def origin(self):
        return (self.zorigin, self.yorigin, self.xorigin)



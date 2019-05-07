import numpy as np

class Molecule():

    def __init__(self, rawHeader, data, size, start, grid_size, cell_dim, density_ranges, origin):
        self.rawHeader = rawHeader
        self.array = data
        (self.nx, self.ny, self.nz) = size
        (self.nxstart, self.nystart, self.nzstart) = start
        (self.mx, self.my, self.mz) = grid_size
        (self.xlen, self.ylen, self.zlen) = cell_dim
        (self.xorigin, self.yorigin, self.zorigin) = origin
        # Assert volume info match header values
        d_mean = np.mean(self.array)
        d_min = np.min(self.array)
        d_max = np.max(self.array)
        if density_ranges != (d_min, d_max, d_mean):
            print("File density range does not match with data")
            (self.dmin, self.dmax, self.dmean) = (d_min, d_max, d_mean)
        else:
            (self.dmin, self.dmax, self.dmean) = density_ranges
        


    def __eq__(self, other):
        if isinstance(other, self.__class__):
            hasSameData = ~(~np.isclose(other.array, self.array)).sum()
            return hasSameData
        return False

    def set_origin(self, new_origin):
        (self.xorigin, self.yorigin, self.zorigin) = new_origin
    
    def data(self):
        return self.array

    def set_data(self, newData):
        self.array = newData

    def shape(self):
        return (self.nx, self.ny, self.nz)

    def start_point(self):
        return (self.nxstart, self.nystart, self.nzstart)

    def grid_size(self):
        return (self.mx, self.my, self.mz)

    def cell_dim(self):
        return (self.xlen, self.ylen, self.zlen)

    def density_range(self):
        return (self.dmin, self.dmax, self.dmean)

    def origin(self):
        return (self.xorigin, self.yorigin, self.zorigin)




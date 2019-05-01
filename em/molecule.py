import numpy as np

class Molecule():

    def __init__(self, rawHeader, data, size, start, grid_size, cell_dim, density_ranges, origin):
        self.rawHeader = rawHeader
        self.array = data
        (self.ny, self.nx, self.nz) = size
        (self.nystart, self.nxstart, self.nzstart) = start
        (self.my, self.mx, self.mz) = grid_size
        (self.ylen, self.xlen, self.zlen) = cell_dim
        (self.yorigin, self.xorigin, self.zorigin) = origin
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
        (self.zorigin, self.yorigin, self.xorigin) = new_origin
    
    def data(self):
        return self.array

    def set_data(self, newData):
        self.array = newData

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




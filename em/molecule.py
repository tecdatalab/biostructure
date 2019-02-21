
class Molecule():


    def __init__(self, rawHeader, data, size, start, grid_size, cell_dim, density_ranges, origin):
        self.rawHeader = rawHeader
        self.data = data
        (self.ny, self.nx, self.nz) = size
        (self.nystart, self.nxstart, self.nzstart) = start
        (self.my, self.mx, self.mz) = grid_size
        (self.ylen, self.xlen, self.zlen) = cell_dim
        (self.dmin, self.dmax, self.dmean) = density_ranges
        (self.yorigin, self.xorigin, self.zorigin) = origin

    
    def data(self):
        return self.data


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




class Molecule():


    def __init__(self, rawHeader, data, size, start, grid_size, cell_dim, density_ranges, origin):
        self.rawHeader = rawHeader
        self.data = data
        self.size = size
        self.start = start
        self.grid_size = grid_size
        self.cell_dim = cell_dim
        self.density_ranges = density_ranges
        self.origin = origin

    
    def getDataArray(self):
        return self.data

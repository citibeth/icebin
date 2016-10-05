import rtree.index
import collections
import math
import scipy

class Smoother(object):
    def __init__(self, nI, inserts_iter, sigma):
        """inserts_iter: (cell-ID, (x,y,x,y), weight) for each cell in grid...
            Information about the grid, as an iterator
        sigma:
            Standard deviation (spatial) for Gaussian smoothing.
        """
        self.nI = nI
        self.cells = collections.OrderedDict()
        self.rtree = rtree.index.Index()
        for tuple in inserts_iter:
            self.cells[tuple[0]] = tuple
            self.rtree.insert(*tuple)

        self.sigma = sigma


    def gaussian_weights(self, x0,y0,cell_mass0):
        """Gives the weights for a Gaussian centered at center=(x,y)
        x0, y0:
            Position (in real coordinate space) to compute Gaussian weights
        cell_mass0:
            Mass of presumed cell at (x0,y0).  Used to normalize.
        returns: [(cell-index, factor)]
            Weights for other cells in the grid.
        """
        two_sigma_squared_inv = 1./(2.*self.sigma*self.sigma)

#        index0,pos0,cell_mass0 = self.cells[index0]
#        x0 = pos0[0]
#        y0 = pos0[1]

        out_raw = list()
        mass_sum = 0.        # Sum of all the weights
        radius = 3.*self.sigma
        radius_squared = radius*radius
        for index in self.rtree.intersection((x0-radius, y0-radius, x0+radius, y0+radius)):
            index,pos,cell_mass = self.cells[index]

            x=pos[0]
            y=pos[1]

            # Theoretical weight for this cell's position
            # (don't worry about normalizing in theory, we will normalize ourselves)
            dx = (x-x0)
            dy = (y-y0)
            distance_squared = dx*dx + dy*dy
            if distance_squared < radius_squared:
                w0 = math.exp(-two_sigma_squared_inv * distance_squared)

                dmass = w0*cell_mass
                mass_sum += dmass
                out_raw.append((index, w0))

        # Normalize sum(w0*cell_mass) = cell_mass0
        factor = cell_mass0 / mass_sum
        return [(index,w*factor) for index,w in out_raw]

    def matrix(self):
        """Computes the smoothing matrix"""

        # Construct sparse matrix...
        M_data = list()
        M_ii = list()
        M_jj = list()
        n=0
        for index0,pos0,mass0 in self.cells.values():
            gweights = self.gaussian_weights(pos0[0], pos0[1], mass0)
            for index1,weight in gweights:
                M_ii.append(index1)
                M_jj.append(index0)
                M_data.append(weight)
            n += 1
            if (n % 100) == 0:
                print('Smoother finished {} cells'.format(n))

        # This is our smoothing matrix!
        M = scipy.sparse.coo_matrix((M_data, (M_ii, M_jj)), shape=(self.nI,self.nI))
        return M


def ibgrid_smoothing_matrix(ibgrid, wI, sigma):
    """Computes a smoothing matrix based on an ibgrid data structure.  See class ibgrid.Grid
    wI:
        Weights from rm.regrid('IvE')
    """

    # Obtain a single point and weight for each basis functions.
    # This is a simplification and approximate.  But it should work
    # pretty well as long as grid cells are small.  And it will work
    # for any kind of grid/mesh.
    def inserts():
        for cell in ibgrid.cells.values():
            if wI[cell.index] > 0:
                centroid = cell.centroid()
                weight = wI[cell.index]
                yield  (cell.index, (centroid[0], centroid[1], centroid[0], centroid[1]), wI[cell.index])
    nI = len(wI)
    return Smoother(nI, inserts(), sigma).matrix()


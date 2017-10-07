# IceBin: A Coupling Library for Ice Models and GCMs
# Copyright (c) 2013-2016 by Elizabeth Fischer
# 
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published
# by the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
# 
# You should have received a copy of the GNU Lesser General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

import netCDF4
import giss.proj
import giss.basemap
import numpy as np
import pyproj
import functools
import operator

# -------------------------------------------------------
class Indexing(object):
    def __init__(self, nc, vname):
        ncvar = nc.variables[vname]
        self.base = ncvar.base
        self.extent = ncvar.extent
        self.indices = ncvar.indices

        # Turn it all into Numpy arrays
        if not isinstance(self.base, np.ndarray):
            self.base = np.array([self.base], dtype=type(self.base))
        if not isinstance(self.extent, np.ndarray):
            self.extent = np.array([self.extent], dtype=type(self.extent))
        if not isinstance(self.indices, np.ndarray):
            self.indices = np.array([self.indices], dtype=type(self.indices))

        self.size = 1
        for ex in self.extent:
            self.size *= ex

        # Shape of a row-major array in memory
        # Row-order ordered indices
        if (self.indices[0] == 0):
            self.shape = self.extent
        else:
            self.shape = self.extent[::-1]    # Reverse

        self.make_strides()

    def __len__(self):
        return functools.reduce(operator.mul, self.extent)


    @property
    def rank(self):
        return len(self.extent)

    def make_strides(self):
        rank = self.rank
        self.strides = np.zeros((rank,),dtype='i')
        self.strides[self.indices[rank-1]] = 1;
        for d in range(rank-2,-1,-1):
            self.strides[self.indices[d]] = self.strides[self.indices[d+1]] * self.extent[self.indices[d+1]];

    def tuple_to_index(self, tuple):
        ix = 0
        for k in range(0,self.rank):
            ix += (tuple[k]-self.base[k]) * self.strides[k];
        return ix;

    def index_to_tuple(self, ix):
        ix0 = ix
        tpl = np.zeros(len(self.shape),dtype='i')
        for d in range(0,self.rank-1):       # indices by descending stride
            k = self.indices[d]
            tpl_k = ix // self.strides[k]
            ix -= tpl_k*self.strides[k]
            tpl[k] = tpl_k + self.base[k]

        tpl[self.indices[self.rank-1]] = ix
        return tuple(tpl)


# -------------------------------------------------------
class Vertex(object):
    __slots__ = ['index', 'x', 'y']

    def __init__(self, *args):
        for attr,val in zip(self.__slots__, args):
            setattr(self, attr, val)

    def __repr__(self):
        return 'Vertex({}, {},{})'.format(self.index, self.x, self.y)


class Cell(object):
    __slots__ = ['index', 'vertices', 'i', 'j', 'k']

    def __init__(self, *args):
        for attr,val in zip(self.__slots__, args):
            setattr(self, attr, val)

    def __repr__(self):
        return 'Cell({}, {})'.format(self.index, self.vertices)

    def area(self):
        sum = 0.
        v0 = self.vertices[-1]
        for v1 in self.vertices:
            sum += v0.x*v1.y - v1.x*v0.y
            v0 = v1
        return .5 * sum

    # https://en.wikipedia.org/wiki/Centroid#Bounded_region
    def centroid(self):
        A2 = 0.        # Will be 2A
        Cx = 0.
        Cy = 0.
        v0 = self.vertices[-1]
        for v1 in self.vertices:
            dA = v0.x*v1.y - v1.x*v0.y
            A2 += dA
            Cx += (v0.x + v1.x) * dA
            Cy += (v0.y + v1.y) * dA
            v0 = v1

        w = 1./(3.*A2)    # 1/6A
        return (w*Cx, w*Cy)

class Grid(object):
    # Read Grid from a netCDF file
    def __init__(self, indexing, vertices, cells, **attrs):
        self.__dict__.update(attrs)

        self.indexing = indexing

        if not isinstance(vertices, dict):
            vertices = {x.index : x for x in vertices}
        if not isinstance(cells, dict):
            cells = {x.index : x for x in cells}
            
        self.vertices = vertices
        self.cells = cells

        # Intuit the number of cells and vertices in the full vector space
        # if not given to us.
        if not hasattr(self, 'cells_num_full'):
            self.cells_num_full = max(cell.index for cell in self.cells.values()) + 1
        if not hasattr(self, 'vertices_num_full'):
            self.vertices_num_full = max(vertex.index for vertex in self.vertices.values()) + 1


        if self.coordinates == 'XY' :
            (self.llproj, self.xyproj) = giss.proj.make_projs(self.projection)
        else :
            self.xyproj = None
            self.llproj = None

    def plot_alt(self, basemap, **kwargs) :
        """Plots the grid cell outlines
    
        Args:
            basemap: Map on which to plot
            **kwargs: Any options passed through to Matplotlib plotting
        """

        npoly = len(self.cells)
        npoints = len(self.vertices)
        xdata = []
        ydata = []

        for cell in self.cells.values():
            xdata.append(cell.vertices[-1].x)
            ydata.append(cell.vertices[-1].y)

            for vertex in cell.vertices:
                xdata.append(vertex.x)
                ydata.append(vertex.y)

            # Add a NaN (to lift the pen)
            xdata.append(np.nan)
            ydata.append(np.nan)

        xdata = np.array(xdata)
        ydata = np.array(ydata)

        if self.xyproj is not None :    # translate xy->ll
            londata, latdata = pyproj.transform(self.xyproj, self.llproj, xdata, ydata)
            londata[np.isnan(xdata)] = np.nan
            latdata[np.isnan(ydata)] = np.nan

        else :      # Already in lon/lat coordinates
            londata = xdata
            latdata = ydata

        giss.basemap.plot(basemap, londata, latdata, **kwargs)


    def plot(self, basemap, cells=None, **kwargs) :
        """Plots the grid cell outlines
    
        basemap: [REQUIRED]
            Map on which to plot
        cells: [OPTIONAL]
            The grid cells to plot.  If None, plots entire grid.
        ax: [OPTIONAL]
            Plot on which to plot (if not implied in basemap)
        **kwargs:
            Other options passed through to Matplotlib plotting
        """

        # Plot entire grid, if user didn't specify
        if cells is None:
            cells = self.cells.values()

        # Get unique line segments
        segments = set()
        for cell in cells:
            v0 = cell.vertices[-1]
            for v1 in cell.vertices:
                segments.add((v0,v1))
                v0 = v1

#       # Batch segments into lines
#       lines = dict()   # vertex0 -> list(All segments starting with vertex0)
#       for segment in segments:
#           v0 = segment[0]
#           if v0 in lines:
#               lines[v0].append(segment[1])
#           else:
#               lines[v0] = [segment[1]]

        # Plot each line segment
        xdata = []
        ydata = []
        for v0,v1 in segments:
            xdata.append(v0.x)
            ydata.append(v0.y)
            xdata.append(v1.x)
            ydata.append(v1.y)
            # Add a NaN (to lift the pen)
            xdata.append(np.nan)
            ydata.append(np.nan)

        xdata = np.array(xdata)
        ydata = np.array(ydata)

        if self.xyproj is not None :    # translate xy->ll
            londata, latdata = pyproj.transform(self.xyproj, self.llproj, xdata, ydata)
            londata[np.isnan(xdata)] = np.nan
            latdata[np.isnan(ydata)] = np.nan

        else :      # Already in lon/lat coordinates
            londata = xdata
            latdata = ydata

        giss.basemap.plot(basemap, londata, latdata, **kwargs)



class Grid_XY(Grid):
    def plotter(self):
        return giss.plot.ProjXYPlotter(
            self.x_boundaries, self.y_boundaries,
            self.projection, transpose=(indexing.indices[0] == 0))

class Grid_LonLat(Grid):
    def plotter(self):
        return giss.plot.LonLatPlotter(
            self.lon_boundaries, self.lat_boundaries,
            boundaries=True, transpose=(indexing.indices[0] == 0))

# -------------------------------------------------------
def read_nc(nc, vname, set_area=False):
    """Read the Grid from a netCDF file.

        nc : netCDF4
            Open netCDF file handle.
        vname : str
            Name of variable from which to read the Grid.
        set_area:
            If True, make sure cell.area is set"""

    attrs = dict()      # Used to initialize grid

    variables = nc.variables
    info = variables[vname + '.info']

    indexing = Indexing(nc, vname + '.indexing')

    # -------------------------------
    # Read all attributes under .info:
    # name, type, cells.num_full, vertices.num_full
    attrs.update(info.__dict__)
    for key in info.__dict__.keys():
        nkey = key.replace('.', '_')
        attrs[nkey] = info.__dict__[key]

    # -------------------------------
    # Read the main grid into temporary arrays
    vertices_index = variables[vname + '.vertices.index'][:]
    vertices_xy = variables[vname + '.vertices.xy'][:]
    cells_index = variables[vname + '.cells.index'][:]
    if vname + '.cells.ijk' in variables:
        cells_ijk = variables[vname + '.cells.ijk'][:]
    else:
        cells_ijk = None
    if vname + '.cells.area' in variables:
        cells_area = variables[vname + '.cells.area'][:]
    else:
        cells_area = None
    cells_vertex_refs = variables[vname + '.cells.vertex_refs'][:]
    cells_vertex_refs_start = variables[vname + '.cells.vertex_refs_start'][:]

    # ---------------------------------
    # Construct data structures
    vertices = dict()
    for index, xy in zip(vertices_index, vertices_xy):
        vertices[index] = Vertex(index, xy[0], xy[1])

    cells = dict()
    for nc_i in range(0, len(cells_index)):
        index = cells_index[nc_i]
        ijk = cells_ijk[nc_i] if cells_ijk is not None else (0,0,0)

        r0 = cells_vertex_refs_start[nc_i]
        r1 = cells_vertex_refs_start[nc_i+1]
        vertex_indices = cells_vertex_refs[r0:r1]   # Vertex indices

        # Convert vertices to a list in the correct order
        vertices_as_list = [vertices[vertex_index] for vertex_index in vertex_indices]
        cells[index] = Cell(index, vertices_as_list, ijk[0], ijk[1], ijk[2])

        # Set the area ourselves if it wasn't set before
        if set_area and cells_area is None:
            cells[index].area = cells[index]._area()

    if attrs['type'] == 'XY':
        attrs['x_boundaries'] = variables[vname + '.x_boundaries'][:]
        attrs['y_boundaries'] = variables[vname + '.y_boundaries'][:]
        return Grid_XY(indexing, vertices, cells, **attrs)

    elif attrs['type'] == 'LONLAT':
        attrs['lon_boundaries'] = variables[vname + '.lon_boundaries'][:]
        attrs['lat_boundaries'] = variables[vname + '.lat_boundaries'][:]
        return Grid_LonLat(indexing, vertices, cells, **attrs)

    elif attrs['type'] == 'MESH':
        return Grid(indexing, vertices, cells, **attrs)

    elif attrs['type'] == 'GENERIC':
        return Grid(indexing, vertices, cells, **attrs)

    raise ValueError('Unknown grid type {}'.format(attrs['type']))

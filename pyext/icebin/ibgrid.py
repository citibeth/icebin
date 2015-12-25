import netCDF4
import giss.proj
import giss.basemap
import numpy as np
import pyproj

# -------------------------------------------------------
class Vertex(object):
	__slots__ = ['index', 'x', 'y']

	def __init__(self, *args):
		for attr,val in zip(self.__slots__, args):
			setattr(self, attr, val)

	def __repr__(self):
		return 'Vertex({}, {},{})'.format(self.index, self.x, self.y)


class Cell(object):
	__slots__ = ['index', 'vertices', 'i', 'j', 'k', 'area']

	def __init__(self, *args):
		for attr,val in zip(self.__slots__, args):
			setattr(self, attr, val)


def read_nc(nc, vname):
	"""Read the Grid from a netCDF file.

	    nc : netCDF4
	        Open netCDF file handle.
		vname : str
	        Name of variable from which to read the Grid."""

	attrs = dict()		# Used to initialize grid

	variables = nc.variables
	info = variables[vname + '.info']

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
	cells_ijk = variables[vname + '.cells.ijk'][:]
	cells_area = variables[vname + '.cells.area'][:]
	cells_vertex_refs = variables[vname + '.cells.vertex_refs'][:]
	cells_vertex_refs_start = variables[vname + '.cells.vertex_refs_start'][:]

	# Correct ISSM file (made with MATLAB?)
	# self.cells_vertex_refs = self.cells_vertex_refs - 1

	# ---------------------------------
	# Construct data structures
	vertices = dict()
	for index, xy in zip(vertices_index, vertices_xy):
		vertices[index] = Vertex(index, xy[0], xy[1])

	cells = dict()
	for nc_i in range(0, len(cells_index)):
		index = cells_index[nc_i]
		ijk = cells_ijk[nc_i]
		area = cells_area[nc_i]

		r0 = cells_vertex_refs_start[nc_i]
		r1 = cells_vertex_refs_start[nc_i+1]
		vertex_refs = cells_vertex_refs[r0:r1]
		vertex_indices = vertices_index[vertex_refs]

		vertices = [self.vertices[vertex_index] for vertex_index in vertex_indices]
		cells[index] = Cell(index, vertices, ijk[0], ijk[1], ijk[2], area)

	if attrs['type'] == 'XY' :
		attrs['x_boundaries'] = variables[vname + '.x_boundaries'][:]
		attrs['y_boundaries'] = variables[vname + '.y_boundaries'][:]

	return Grid(attrs, vertices, cells)


class Grid(object):
	# Read Grid from a netCDF file
	def __init__(self, vertices, cells, **attrs):
		self.__dict__.update(attrs)

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

		if self.xyproj is not None :	# translate xy->ll
			londata, latdata = pyproj.transform(self.xyproj, self.llproj, xdata, ydata)
			londata[np.isnan(xdata)] = np.nan
			latdata[np.isnan(ydata)] = np.nan

		else :		# Already in lon/lat coordinates
			londata = xdata
			latdata = ydata

		giss.basemap.plot(basemap, londata, latdata, **kwargs)


	def plot(self, cells=None, basemap=None, **kwargs) :
		"""Plots the grid cell outlines
	
		cells: [OPTIONAL]
			The grid cells to plot.  If None, plots entire grid.
		basemap: [REQUIRED]
			Map on which to plot
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

#		# Batch segments into lines
#		lines = dict()   # vertex0 -> list(All segments starting with vertex0)
#		for segment in segments:
#			v0 = segment[0]
#			if v0 in lines:
#				lines[v0].append(segment[1])
#			else:
#				lines[v0] = [segment[1]]

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

		if self.xyproj is not None :	# translate xy->ll
			londata, latdata = pyproj.transform(self.xyproj, self.llproj, xdata, ydata)
			londata[np.isnan(xdata)] = np.nan
			latdata[np.isnan(ydata)] = np.nan

		else :		# Already in lon/lat coordinates
			londata = xdata
			latdata = ydata

		giss.basemap.plot(basemap, londata, latdata, **kwargs)




	def plotter(self) :
		if self.type == 'XY' :
			return giss.plot.ProjXYPlotter(self.x_boundaries, self.y_boundaries, self.projection)
		return None

# -------------------------------------------------------


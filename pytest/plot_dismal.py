# GLINT2: A Coupling Library for Ice Models and GCMs
# Copyright (c) 2013 by Robert Fischer
# 
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

# Plot stuff on the searise grid

import netCDF4
import giss.basemap
import giss.modele
import sys
import glint2
import numpy as np
import giss.plot

def _default_plot_boundaries(basemap) :
    """Default function for 'plot_boundaries' output of plot_params()"""

    # ---------- Draw other stuff on the map
    # draw line around map projection limb.
    # color background of map projection region.
    # missing values over land will show up this color.
    basemap.drawmapboundary(fill_color='0.5')
    basemap.drawcoastlines()

    # draw parallels and meridians, but don't bother labelling them.
    basemap.drawparallels(np.arange(-90.,120.,30.))
    basemap.drawmeridians(np.arange(0.,420.,60.))
# --------------------------------
_reverse_scale = {'mass'}


nc_name = sys.argv[1]
var_name = sys.argv[2]

# ------------- Read SeaRise Grid
#nc = netCDF4.Dataset('searise.nc')
#grid = glint2.pyGrid(nc, 'grid')
#nc.close()



basemap = giss.basemap.greenland_laea()

nc = netCDF4.Dataset(nc_name)

#pp = giss.modele.plot_params(var_name, nc=nc)
pp = {}
pp['var_name'] = var_name
pp['val'] = nc.variables[var_name][:]
pp['title'] = "%s:%s" % (nc_name, var_name)
pp['plot_boundaries'] = _default_plot_boundaries
pp['plotter'] = glint2.Grid_read_plotter('searise.nc', 'grid')

plot_args = {}
pp['plot_args'] = plot_args
plot_args['norm'] = giss.plot.AsymmetricNormalize()
reverse = (var_name in _reverse_scale)
plot_args['cmap'] = giss.plot.cpt('giss-cpt/BlRe.cpt', reverse=reverse).cmap

giss.plot.plot_var(basemap=basemap, **pp)       # Plot, and show on screen

# Slightly more complex alternatives:
# Save figure:
#   giss.plot.plot_var(fname='plottest1.png', **pp)
# Save figure and snow on screen
#   giss.plot.plot_var(fname='plottest1.png', show=True, **pp)

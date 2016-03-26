# IceBin: A Coupling Library for Ice Models and GCMs
# Copyright (c) 2013-2016 by Elizabeth Fischer
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

import unittest
import sys
import glint2
import netCDF4
import numpy as np
import matplotlib.pyplot
import giss.basemap
import numpy.ma as ma
import giss.util
import giss.modele
import giss.searise

# -----------------------------------------------------------
class Glint2Tests(unittest.TestCase) :

    def assertSigEqual(self, a, b, epsilon) :
        self.assertTrue(abs(a - b) < epsilon * abs(b))

    def assertNotSigEqual(self, a, b, epsilon) :
        self.assertTrue(abs(a - b) >= epsilon * abs(b))

    def setUp(self) :
        self.grid1_fname = 'greenland_2x2_5.nc'
        self.grid2_fname = 'searise.nc'
        self.exgrid_fname = 'greenland_2x2_5-searise.nc'
        self.data_fname = 'data/JUL1950.ijhchc1k225.nc'
        self.vname = 'impm_lndice'
        self.searise_fname = 'data/Greenland_5km_v1.1.nc'

        # ---------- Read from Exchange Grid
        nc = netCDF4.Dataset(self.exgrid_fname)
        self.overlap = glint2.read_overlap(nc, 'grid')
        self.sproj = nc.variables['grid.info'].projection
        self.exgrid_name = nc.variables['grid.info'].name
        nc.close()

        # ---------- Read from grid1
        self.grid1 = glint2.Grid(self.grid1_fname, 'grid')
        self.proj_area1 = self.grid1.get_proj_area(self.sproj)
        self.native_area1 = self.grid1.get_native_area()

        nc = netCDF4.Dataset(self.grid1_fname)
        info_var = nc.variables['grid.info']
        self.grid1_name = info_var.name
        nc.close()

        # --------- Read from grid2
        nc = netCDF4.Dataset(self.grid2_fname)
        info_var = nc.variables['grid.info']
        self.grid2_name = info_var.name
        nc.close()

        # ----------- Read from Searise (more grid2-related stuff)
        (self.elev2, self.mask2) = giss.searise.read_elevation2_mask2(self.searise_fname)

        # ---------- Read some data
        nc = netCDF4.Dataset(self.data_fname)
        self.val1h = giss.modele.read_ncvar(nc, self.vname).astype(np.double)
        if self.vname == 'impm_lndice' :
            self.val1h *= 86400.        # Convert to more sensible unit
        nc.close()
        self.nhp = self.val1h.shape[0]


    def mytest_remap(self, n2p) :
        # Make sure we have the right file
        self.assertEqual(self.exgrid_name, self.grid1_name + '-' + self.grid2_name)

        # ---------- Read some data
        val1 = self.val1h[0,:]

        # =========== Mask out self.overlap matrix (maske=True means ignore)
        # (to get correct conservation numbers)
        mask1 = (np.isnan(val1)).astype(np.int32)
        mask2 = None
        overlap_masked = glint2.mask_out(self.overlap, mask1, mask2)
        area1_masked = giss.util.sum_by_rows(overlap_masked)
        area2_masked = giss.util.sum_by_cols(overlap_masked)

        # ========== Now try it out

        # Matrix to go from PROJECTED grid1 to grid2
        mat = glint2.grid1_to_grid2(overlap_masked)

        # Diagonal matrix converts from a vector of values in the native
        # grid to a vector of values in the projected grid.
        if n2p :
            diag = glint2.proj_native_area_correct(self.grid1, self.sproj, 'n2p')
            glint2.multiply_bydiag(mat, diag)

        # Regrid
        val2 = glint2.coo_multiply(mat, val1, fill=np.nan)

        # Conservation sum with "native" grid1 on a sphere
        total1_native = np.nansum(area1_masked * (self.native_area1 / self.proj_area1) * val1.reshape(-1))
        # Conservation sum with "projected" grid1 on a plane
        total1_proj = np.nansum(area1_masked * val1.reshape(-1))
        total2 = np.nansum(area2_masked * val2)

        print 'n2p           = %d' % n2p
        print 'total1_native = %f' % total1_native
        print 'total1_proj   = %f' % total1_proj
        print 'total2        = %f' % total2


        # Make sure we got it right
        if n2p :
            total1 = total1_native
        else :
            total1 = total1_proj
#       self.assertTrue(abs(total1 - total2) < 1e-12 * abs(total1))
        self.assertSigEqual(total1, total2, epsilon=1e-12)

        # ================ Plot It!
        # Plot multiple plots on one page
        figure = matplotlib.pyplot.figure(figsize=(8.5,11))

        # ---------- Plot the GCM Grid
        ax = figure.add_subplot(121)
        plotter = glint2.Grid_read_plotter(self.grid1_fname, 'grid')
        mymap = giss.basemap.greenland_laea(ax)
        #mymap = giss.basemap.north_laea(ax)
        im = plotter.pcolormesh(mymap, val1)
        mymap.colorbar(im, "right", size='5%')
        mymap.drawcoastlines()

        # ---------- Plot the ice grid
        ax = figure.add_subplot(122)
        plotter = glint2.Grid_read_plotter(self.grid2_fname, 'grid')
        mymap = giss.basemap.greenland_laea(ax)
        #mymap = giss.basemap.north_laea(ax)
        im = plotter.pcolormesh(mymap, val2)
        mymap.colorbar(im, "right", size='5%')
        mymap.drawcoastlines()

        fname = 'test_remap-%s.png' % ('native' if n2p else 'proj')
        figure.savefig(fname)
        print '====> Wrote output to %s' % fname


    # --------------------------------------------------------
    def mytest_hc(self, n2p) :

        # ---------- Set up height classes & corresponding height points
        hcmax = [25.]
        hcmax.extend(np.array(range(1,40))*100.0)
        hcmax = np.array(hcmax)
        print hcmax

        hpdefs = np.array([
                   0.,    50.,   150.,   250.,   350.,   450.,   550.,   650.,
                 750.,   850.,   950.,  1050.,  1150.,  1250.,  1350.,  1450.,
                1550.,  1650.,  1750.,  1850.,  1950.,  2050.,  2150.,  2250.,
                2350.,  2450.,  2550.,  2650.,  2750.,  2850.,  2950.,  3050.,
                3150.,  3250.,  3350.,  3450.,  3550.,  3650.,  3750.,  3850.])

        # ========== Height-classify overlap matrix
        self.assertEqual(len(hcmax), self.nhp)
        overlaph = glint2.height_classify(self.overlap, self.elev2, hcmax)

        # =========== Mask out overlap matrix
        # (to get correct conservation numbers)
        mask1h = (np.isnan(self.val1h)).astype(np.int32)
        #mask2 = mask2
        overlaph_m = glint2.mask_out(overlaph, mask1h, self.mask2)

        area1h_masked = giss.util.sum_by_rows(overlaph_m)
        area2_masked = giss.util.sum_by_cols(overlaph_m)

        # Matrix to go from PROJECTED grid1 to grid2
        mat_1hc_to_2 = glint2.grid1_to_grid2(overlaph_m)

        # ...and back
        mat_2_to_1hc = glint2.grid2_to_grid1(overlaph_m)

        # Diagonal matrix converts from a vector of values in the native
        # grid to a vector of values in the projected grid.
        if n2p :
            diag = glint2.proj_native_area_correct(self.grid1, self.sproj, 'n2p')
            # Convert the native area correction to height-classified space
            diag = np.tile(diag, (1,self.nhp)).reshape(-1)
            glint2.multiply_bydiag(mat_1hc_to_2, diag)

        # -------- Area-weighted remapping using height classes
        val2_hc = glint2.coo_multiply(mat_1hc_to_2, self.val1h, fill=np.nan)

        # ========== Interpolate to ice grid, using height point smoothing
        # Collapse height-classified mask1 into non-height-classified version
        mask1 = np.array(mask1h.sum(axis=0)).reshape(mask1h.shape[1:]).astype(np.int32)
        # Mask overlap matrix using this (for correct conservation checking)
        overlap_m = glint2.mask_out(self.overlap, mask1, self.mask2)
        mat_1hp_to_2 = glint2.hp_interp(overlap_m, self.elev2, hpdefs)
        val2_hp = glint2.coo_multiply(mat_1hp_to_2, self.val1h, fill=np.nan)

        # This is what you get when you start with val1h
        # (interpreted as height points), regrid to grid2,
        # then back to grid1h (interpreted as height classes)
        val1hcx = glint2.coo_multiply(mat_2_to_1hc, val2_hp, fill=np.nan)

        # ===============================================
        # Evaluate conservation

        # A bit more complication extending native_area1 to full hc scheme
        factor1h = np.tile(self.native_area1 / self.proj_area1, (self.nhp,1)).reshape(-1)

        total1h_native = np.nansum(area1h_masked * factor1h * self.val1h.reshape(-1))
        total1h_proj = np.nansum(area1h_masked * self.val1h.reshape(-1))
        total2_hc = np.nansum(area2_masked * val2_hc)
        total2_hp = np.nansum(area2_masked * val2_hp)
        total1hcx = np.nansum(area1h_masked * val1hcx)

        print 'total1h_native = %f' % total1h_native
        print 'total1h_proj   = %f' % total1h_proj
        print 'total2_hc      = %f' % total2_hc
        print 'total2_hp      = %f' % total2_hp
        print 'total1hcx      = %f' % total1hcx

        if n2p :
            self.assertSigEqual(total1h_native, total2_hc, epsilon=1e-12)
        else :
            self.assertSigEqual(total1h_proj, total2_hc, epsilon=1e-12)

        self.assertNotSigEqual(total1h_native, total1h_proj, epsilon=1e-12)
        self.assertSigEqual(total2_hp, total1hcx, epsilon=1e-12)


        # Plot multiple plots on one page
        figure = matplotlib.pyplot.figure(figsize=(8.5,11))

        # ---------- Plot the ice grid (height point interpolation)
        ax = figure.add_subplot(121)
        plotter = glint2.Grid_read_plotter(self.grid2_fname, 'grid')
        mymap = giss.basemap.greenland_laea(ax)
        #mymap = giss.basemap.north_laea(ax)
        im = plotter.pcolormesh(mymap, val2_hp)
        mymap.colorbar(im, "right", size='5%')
        mymap.drawcoastlines()


        # ---------- Plot with height classes, after going through ice grid
        ax = figure.add_subplot(122)
        plotter = glint2.Plotter_HC(
            glint2.Grid_read_plotter(self.grid2_fname, 'grid'),
            overlap_m, self.elev2, hcmax)
        mymap = giss.basemap.greenland_laea(ax)
        #mymap = giss.basemap.north_laea(ax)
        im = plotter.pcolormesh(mymap, val1hcx)
        mymap.colorbar(im, "right", size='5%')
        mymap.drawcoastlines()


        fname = 'test_hc-%s.png' % ('native' if n2p else 'proj')
        figure.savefig(fname)
        print '====> Wrote output to %s' % fname



    def test_remap_proj(self) :
        self.mytest_remap(False)

    def test_remap_native(self) :
        self.mytest_remap(True)

    def test_hc_proj(self) :
        self.mytest_hc(False)

    def test_hc_native(self) :
        self.mytest_hc(True)




suite = unittest.TestLoader().loadTestsFromTestCase(Glint2Tests)
unittest.TextTestRunner(verbosity=2).run(suite)

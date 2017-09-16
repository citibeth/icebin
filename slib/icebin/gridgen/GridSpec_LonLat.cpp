/*
 * IceBin: A Coupling Library for Ice Models and GCMs
 * Copyright (c) 2013-2016 by Elizabeth Fischer
 * 
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published
 * by the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 * 
 * You should have received a copy of the GNU Lesser General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include <functional>
#include <cmath>

#include <ibmisc/geodesy.hpp>
#include <ibmisc/netcdf.hpp>

#include <icebin/Grid.hpp>
#include <icebin/error.hpp>

#include <icebin/gridgen/GridSpec_LonLat.hpp>
#include <icebin/gridgen/gridutil.hpp>

using namespace netCDF;
using namespace ibmisc;

namespace icebin {


// Radians <--> Degrees
static const double D2R = M_PI / 180.0;
static const double R2D = 180.0 / M_PI;

template <typename T>
inline int sgn(T val) {
    return (val > T(0)) - (val < T(0));
}

// ---------------------------------------------------------
std::vector<double> GridSpec_LonLat::latc() const
{
    std::vector<double> ret;
    ret.reserve(nlat()+1);

    if (south_pole) ret.push_back(-90);
    for (int j=0; j<latb.size()-1; ++j)
        ret.push_back(.5 * (latb[j] + latb[j+1]));
    if (north_pole) ret.push_back(90);

    return ret;
}

std::vector<double> GridSpec_LonLat::lonc() const
{
    std::vector<double> ret;
    ret.reserve(nlon()-1);

    for (int j=0; j<lonb.size()-1; ++j)
        ret.push_back(.5 * (lonb[j] + lonb[j+1]));

    return ret;
}


// ---------------------------------------------------------

/** Computes the exact surface area of a lat-lon grid box on the
surface of the sphere.
<b>NOTES:</b>
<ol><li>lon0 must be numerically less than lon1.</li>
<li>lat0 and lat1 cannot cross the equator.</li></ol> */
static double graticule_area_exact(double eq_rad,
    double lat0_deg, double lat1_deg,
    double lon0_deg, double lon1_deg)
{
    double delta_lon_deg = lon1_deg - lon0_deg;
    // Normalize within a 360-degree range (gridutil.hpp)
    delta_lon_deg = loncorrect(delta_lon_deg, 0);

    double lat0 = lat0_deg * D2R;
    double lat1 = lat1_deg * D2R;
    double delta_lon = delta_lon_deg * D2R;

    return delta_lon * (eq_rad*eq_rad) * (sin(lat1) - sin(lat0));
}

/** The polar graticule is a (latitudinal) circle centered on the pole.
This computes its area. */
double polar_graticule_area_exact(double eq_rad,
    double radius_deg)
{
    // See http://en.wikipedia.org/wiki/Spherical_cap
    double theta = radius_deg * D2R;
    return 2.0 * M_PI * (eq_rad * eq_rad) * (1.0 - cos(theta));
}

// ---------------------------------------------------------
/** Set up a GridSpec_LonLat with a given specification.
@param spherical_clip Only realize grid cells that pass this test (before projection).
@see EuclidianClip, SphericalClip
*/
void GridSpec_LonLat::make_grid(Grid_LonLat &grid)
{
printf("BEGIN GridSpec_LonLat::make_grid()\n");
    grid.clear();
    grid.type = Grid::Type::LONLAT;
    grid.coordinates = Grid::Coordinates::LONLAT;
    grid.parameterization = Grid::Parameterization::L0;
    grid.name = this->name;
    grid.points_in_side = points_in_side;

    // Error-check the input parameters
    if (this->south_pole && this->latb[0] == -90.0) {
        (*icebin_error)(-1,
            "latb[] cannot include -90.0 if you're including the south pole cap");
    }
    if (this->north_pole && this->latb.back() == 90.0) {
        (*icebin_error)(-1,
            "latb[] cannot include 90.0 if you're including the north pole cap");
    }


    // Set up to project lines on sphere (and eliminate duplicate vertices) 
    VertexCache vcache(&grid);

    // ------------------- Set up the GCM Grid
    const int south_pole_offset = (this->south_pole ? 1 : 0);
    const int north_pole_offset = (this->north_pole ? 1 : 0);

    grid.cells._nfull = this->nlon() * this->nlat();

    // Get a bunch of points.  (i,j) is gridcell's index in canonical grid
    for (int ilat=0; ilat < this->latb.size()-1; ++ilat) {
        double lat0 = this->latb[ilat];
        double lat1 = this->latb[ilat+1];

        for (int ilon=0; ilon< this->lonb.size()-1; ++ilon) {
            Cell cell;
            double lon0 = this->lonb[ilon];
            double lon1 = this->lonb[ilon+1];

            // Figure out how to number this grid cell
            cell.j = ilat + south_pole_offset;  // 0-based 2-D index
            cell.i = ilon;
            cell.index = indexing.tuple_to_index<int,2>({cell.i, cell.j});
            cell.native_area = graticule_area_exact(this->eq_rad, lat0,lat1,lon0,lon1);

            if (!this->spherical_clip(cell.index, lon0, lat0, lon1, lat1)) continue;

            // Project the grid cell boundary to a planar polygon
            int n = this->points_in_side;

            // Pre-compute our points so we use exact same ones each time.
            std::vector<double> lons;
            lons.reserve(n+1);
            for (int i=0; i<=n; ++i)
                lons.push_back(lon0 + (lon1-lon0) * ((double)i/(double)n));

            std::vector<double> lats;
            lats.reserve(n+1);
            for (int i=0; i<=n; ++i)
                lats.push_back(lat0 + (lat1-lat0) * ((double)i/(double)n));

            // Build a square out of them (in lon/lat space)
            for (int i=0; i<n; ++i)
                vcache.add_vertex(cell, lons[i], lat0);

            for (int i=0; i<n; ++i)
                vcache.add_vertex(cell, lon1, lats[i]);

            // Try to keep calculations EXACTLY the same for VertexCache
            for (int i=n; i>0; --i)
                vcache.add_vertex(cell, lons[i], lat1);

            for (int i=n; i>0; --i)
                vcache.add_vertex(cell, lon0, lats[i]);

//printf("Adding lon/lat cell %d (%d, %d) area=%f\n", cell.index, cell.i, cell.j, cell.area);
            grid.cells.add(std::move(cell));
        }
    }

    // Make the polar caps (if this grid specifies them)

    // North Pole cap
    double lat = this->latb.back();
    Cell pole;
    pole.i = this->nlon()-1;
    pole.j = this->nlat();
    long index = (pole.j * this->nlon() + pole.i);
    if (this->north_pole && this->spherical_clip(index, 0, lat, 360, 90)) {
        for (int ilon=0; ilon< this->lonb.size()-1; ++ilon) {
            double lon0 = this->lonb[ilon];
            double lon1 = this->lonb[ilon+1];

            int n = this->points_in_side;
            for (int i=0; i<n; ++i) {
                double lon = lon0 + (lon1-lon0) * ((double)i/(double)n);
                pole.add_vertex(vcache.add_vertex(lon, lat));
            }
        }

        pole.index = index;
        pole.native_area = polar_graticule_area_exact(this->eq_rad, 90.0 - lat);

        grid.cells.add(std::move(pole));
    }

    // South Pole cap
    index = 0;
    lat = this->latb[0];
    if (this->south_pole && this->spherical_clip(index, 0, -90, 360, lat)) {
        Cell pole;
        for (int ilon=this->lonb.size()-1; ilon >= 1; --ilon) {
            double lon0 = this->lonb[ilon];     // Make the circle counter-clockwise
            double lon1 = this->lonb[ilon-1];

            int n = this->points_in_side;
            for (int i=0; i<n; ++i) {
                double lon = lon0 + (lon1-lon0) * ((double)i/(double)n);
                pole.add_vertex(vcache.add_vertex(lon, lat));
            }
        }
        pole.i = 0;
        pole.j = 0;
        pole.index = index;
        pole.native_area = polar_graticule_area_exact(this->eq_rad, 90.0 + lat);

        grid.cells.add(std::move(pole));
    }

    grid.lonb = std::move(lonb);
    grid.latb = std::move(latb);
    grid.south_pole = south_pole;
    grid.north_pole = north_pole;
    grid.indexing = indexing;
printf("END GridSpec_LonLat::make_grid()\n");
}

// ---------------------------------------------------------


}

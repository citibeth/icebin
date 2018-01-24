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

#include <icebin/gridgen/GridGen_LonLat.hpp>
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
static double polar_graticule_area_exact(double eq_rad,
    double radius_deg)
{
    // See http://en.wikipedia.org/wiki/Spherical_cap
    double theta = radius_deg * D2R;
    return 2.0 * M_PI * (eq_rad * eq_rad) * (1.0 - cos(theta));
}

// ---------------------------------------------------------
/** Set up a GridGen_LonLat with a given specification.
@param spherical_clip Only realize grid cells that pass this test (before projection).
@see EuclidianClip, SphericalClip
*/
Grid make_grid(
    std::string const &name,
    GridSpec_LonLat const &spec,
    std::function<bool(long, double,double,double,double)> spherical_clip)
{
    Indexing indexing({"lon", "lat"}, {0,0}, {spec.nlon(), spec.nlat()}, spec.indices);

    GridMap<Vertex> vertices(-1);    // Unknown and don't care how many vertices in full grid
    GridMap<Cell> cells(spec.nlon() * spec.nlat());

    VertexCache vcache(&vertices);

    // ------------------- Set up the GCM Grid
    const int south_pole_offset = (spec.south_pole ? 1 : 0);
    const int north_pole_offset = (spec.north_pole ? 1 : 0);

    // Get a bunch of points.  (i,j) is gridcell's index in canonical grid
    for (int ilat=0; ilat < spec.latb.size()-1; ++ilat) {
        double lat0 = spec.latb[ilat];
        double lat1 = spec.latb[ilat+1];

        for (int ilon=0; ilon< spec.lonb.size()-1; ++ilon) {
            Cell cell;
            double lon0 = spec.lonb[ilon];
            double lon1 = spec.lonb[ilon+1];

            // Figure out how to number this grid cell
            cell.j = ilat + south_pole_offset;  // 0-based 2-D index
            cell.i = ilon;
            cell.index = indexing.tuple_to_index<int,2>({cell.i, cell.j});
            cell.native_area = graticule_area_exact(spec.eq_rad, lat0,lat1,lon0,lon1);

            if (!spherical_clip(cell.index, lon0, lat0, lon1, lat1)) continue;

            // Project the grid cell boundary to a planar polygon
            int n = spec.points_in_side;

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

            cells.add(std::move(cell));
        }
    }

    // Make the polar caps (if this grid specifies them)

    // North Pole cap
    double lat = spec.latb.back();
    Cell pole;
    pole.i = spec.nlon()-1;
    pole.j = spec.nlat();
    long index = (pole.j * spec.nlon() + pole.i);
    if (spec.north_pole && spherical_clip(index, 0, lat, 360, 90)) {
        for (int ilon=0; ilon< spec.lonb.size()-1; ++ilon) {
            double lon0 = spec.lonb[ilon];
            double lon1 = spec.lonb[ilon+1];

            int n = spec.points_in_side;
            for (int i=0; i<n; ++i) {
                double lon = lon0 + (lon1-lon0) * ((double)i/(double)n);
                pole.add_vertex(vcache.add_vertex(lon, lat));
            }
        }

        pole.index = index;
        pole.native_area = polar_graticule_area_exact(spec.eq_rad, 90.0 - lat);

        cells.add(std::move(pole));
    }

    // South Pole cap
    index = 0;
    lat = spec.latb[0];
    if (spec.south_pole && spherical_clip(index, 0, -90, 360, lat)) {
        Cell pole;
        for (int ilon=spec.lonb.size()-1; ilon >= 1; --ilon) {
            double lon0 = spec.lonb[ilon];     // Make the circle counter-clockwise
            double lon1 = spec.lonb[ilon-1];

            int n = spec.points_in_side;
            for (int i=0; i<n; ++i) {
                double lon = lon0 + (lon1-lon0) * ((double)i/(double)n);
                pole.add_vertex(vcache.add_vertex(lon, lat));
            }
        }
        pole.i = 0;
        pole.j = 0;
        pole.index = index;
        pole.native_area = polar_graticule_area_exact(spec.eq_rad, 90.0 + lat);

        cells.add(std::move(pole));
    }

    return Grid(name,
        std::unique_ptr<GridSpec>(new GridSpec_LonLat(spec)),
        GridCoordinates::XY, "",
        GridParameterization::L0,
        std::move(indexing),
        std::move(vertices), std::move(cells));


}

AbbrGrid make_abbr_grid(
    std::string const &name,
    GridSpec_LonLat const &spec,
    spsparse::SparseSet<long,int> &&dim)    // Indices of gridcells we want to keep

    // std::function<bool(Cell const &)> spherical_clip = &SphericalClip::keep_all)
{
    Indexing indexing({"lon", "lat"}, {0,0}, {spec.nlon(), spec.nlat()}, spec.indices);

    int const N = dim.dense_extent();
    blitz::Array<int,2> ijk(N,3);
    blitz::Array<double,1> native_area(N);

    // Pre-compute area of grid cells
    blitz::Array<double,1> dxyp(spec.nlat());
    for (int j=0; j<spec.nlat(); ++j) {
        dxyp(j) = sin(spec.latb[j+1]) - sin(spec.latb[j]);
    }

    double const D2R_R2 = D2R * spec.eq_rad * spec.eq_rad;
    for (int iId=0; iId<dim.dense_extent(); ++iId) {
        long iI = dim.to_sparse(iId);
        int * const _ij = &ijk(iId,0);
        auto &i(_ij[0]);
        auto &j(_ij[1]);
        indexing.index_to_tuple(_ij, iI);    // Fills in i&j of ijk
        native_area(iId) = dxyp(j) * (spec.lonb[i+1] - spec.lonb[i]) * D2R_R2;
    }


    return AbbrGrid(
        std::unique_ptr<GridSpec>(new GridSpec_LonLat(spec)),
        GridCoordinates::LONLAT, GridParameterization::L0,
        std::move(indexing), name, "",
        std::move(dim), ijk, native_area,
        blitz::Array<double,2>());
}

// ---------------------------------------------------------


}

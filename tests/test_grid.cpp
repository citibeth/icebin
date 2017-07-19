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

// https://github.com/google/googletest/blob/master/googletest/docs/Primer.md

#include <iostream>
#include <cstdio>
#include <netcdf>
#include <gtest/gtest.h>
#include <icebin/Grid.hpp>
#include <icebin/gridgen/clippers.hpp>
#include <icebin/gridgen/GridSpec_LonLat.hpp>
#ifdef BUILD_MODELE
#include <icebin/modele/GridSpec_Hntr.hpp>
#endif

using namespace std::placeholders;  // for _1, _2, _3...
using namespace ibmisc;
using namespace icebin;
using namespace netCDF;

#if 0
bool operator==(Vertex const &a, Vertex const &b)
{
    return  ((a.index == b.index) && (a.x == b.x) && (a.y == b.y));
}

bool operator==(Cell const &a, Cell const &b)
{
    if (a.size() != b.size()) return false;

    if ((a.index != b.index) || (a.native_area != b.native_area) || (a.i != b.i) || (a.j != b.j) || (a.k != b.k)) return false;

    for (auto iia(a.begin()), iib(b.begin()); iia != a.end(); ++iia, ++iib) {
        if (!(*iia == *iib)) return false;
    }
    return true;
}

bool operator==(Grid const &a, Grid const &b)
{

    std::vector<Cell const *> acells(a.cells.sorted());
    std::vector<Cell const *> bcells(b.cells.sorted());
    if (acells.size() != bcells.size()) return false;
    for (auto iia(acells.begin()), iib(bcells.begin());
        iia != acells.end(); ++iia, ++iib)
    {
std::cout << "Comparing: " << **iia << " --- " << **iib << std::endl;
std::cout << "        == " << (**iia == **iib) << std::endl;
        if (!(**iia == **iib)) return false;
    }

    std::vector<Vertex const *> avertices(a.vertices.sorted());
    std::vector<Vertex const *> bvertices(b.vertices.sorted());
    if (avertices.size() != bvertices.size()) return false;
    for (auto iia(avertices.begin()), iib(bvertices.begin());
        iia != avertices.end(); ++iia, ++iib)
    {
std::cout << "Comparing: " << **iia << " --- " << **iib << std::endl;
        if (!(**iia == **iib)) return false;
    }

    return true;
}
#endif


// The fixture for testing class Foo.
class GridTest : public ::testing::Test {
protected:

    std::vector<std::string> tmpfiles;

    // You can do set-up work for each test here.
    GridTest() {}

    // You can do clean-up work that doesn't throw exceptions here.
    virtual ~GridTest()
    {
        for (auto ii(tmpfiles.begin()); ii != tmpfiles.end(); ++ii) {
//          ::remove(ii->c_str());
        }
    }

    // If the constructor and destructor are not enough for setting up
    // and cleaning up each test, you can define the following methods:

    // Code here will be called immediately after the constructor (right
    // before each test).
    virtual void SetUp() {}

    // Code here will be called immediately after each test (right
    // before the destructor).
    virtual void TearDown() {}

//    // The mock bar library shaed by all tests
//    MockBar m_bar;


    void expect_eq(Vertex const &a, Vertex const &b)
    {
        EXPECT_EQ(a.index, b.index);
        EXPECT_EQ(a.x, b.x);
        EXPECT_EQ(a.y, b.y);
    }

    void expect_eq(Cell const &a, Cell const &b)
    {
        EXPECT_EQ(a.size(), b.size());

        EXPECT_EQ(a.index, b.index);
        EXPECT_EQ(a.native_area, b.native_area);
        EXPECT_EQ(a.i, b.i);
        EXPECT_EQ(a.j, b.j);
        EXPECT_EQ(a.k, b.k);

        for (auto iia(a.begin()), iib(b.begin()); iia != a.end(); ++iia, ++iib) {
            expect_eq(*iia, *iib);
        }
    }

    void expect_eq(Grid const &a, Grid const &b)
    {

        std::vector<Cell const *> acells(a.cells.sorted());
        std::vector<Cell const *> bcells(b.cells.sorted());
        EXPECT_EQ(acells.size(), bcells.size());
        for (auto iia(acells.begin()), iib(bcells.begin());
            iia != acells.end(); ++iia, ++iib)
        {
            expect_eq(**iia, **iib);
        }

        std::vector<Vertex const *> avertices(a.vertices.sorted());
        std::vector<Vertex const *> bvertices(b.vertices.sorted());
        EXPECT_EQ(avertices.size(), bvertices.size());
        for (auto iia(avertices.begin()), iib(bvertices.begin());
            iia != avertices.end(); ++iia, ++iib)
        {
            expect_eq(**iia, **iib);
        }
    }


};

TEST_F(GridTest, create_grid)
{
    Grid grid;
    grid.type = Grid::Type::XY;
    grid.name = "Test Grid";
    grid.coordinates = Grid::Coordinates::XY;
    grid.parameterization = Grid::Parameterization::L0;

    Vertex *vertex;
    auto &vertices(grid.vertices);
    vertices.add(Vertex(0,0));
    vertices.add(Vertex(1,0));
    vertices.add(Vertex(2,0));
    vertices.add(Vertex(0,1));
    vertices.add(Vertex(1,1));
    vertex = vertices.add(Vertex(2,1));

    expect_eq(*vertex, *vertex);

    auto &cells(grid.cells);
    Cell *cell;
    cell = cells.add(Cell({vertices.at(0), vertices.at(1), vertices.at(4), vertices.at(3)}));
        cell->i = 0;
        cell->j = 0;
        cell->native_area = 2.;

    cell = cells.add(Cell({vertices.at(1), vertices.at(2), vertices.at(5), vertices.at(4)}));
        cell->i = 1;
        cell->j = 0;
        cell->native_area = 3.;

    expect_eq(*cell, *cell);

    EXPECT_EQ(2., cells.at(0)->native_area);
    EXPECT_EQ(3., cells.at(1)->native_area);

    EXPECT_DOUBLE_EQ(1., cells.at(0)->proj_area(NULL));
    EXPECT_DOUBLE_EQ(1., cells.at(1)->proj_area(NULL));

    expect_eq(grid, grid);

    // ---------------- Write to NetCDF
    // If the constructor and destructor are not enough for setting up
    std::string fname("__netcdf_test.nc");
    tmpfiles.push_back(fname);
    ::remove(fname.c_str());
    {
        ibmisc::NcIO ncio(fname, NcFile::replace);
        grid.ncio(ncio, "grid");
        ncio.close();
    }

    // ---------------- Read from NetCDF
    // If the constructor and destructor are not enough for setting up
    {
        Grid grid2;
        ibmisc::NcIO ncio(fname, NcFile::read);
        grid2.ncio(ncio, "grid");
        ncio.close();

        expect_eq(grid2, grid);
    }

}

TEST_F(GridTest, centroid)
{
    std::vector<Vertex> vertices;
    std::vector<Vertex *> vertices_p;

    vertices.push_back(Vertex(0,0));
    vertices.push_back(Vertex(2,0));
    vertices.push_back(Vertex(2,2));
    vertices.push_back(Vertex(0,2));
    for (auto &v : vertices) vertices_p.push_back(&v);

    {
        Cell cell0(std::move(vertices_p));
        Point pt(cell0.centroid());
        EXPECT_EQ(pt.x, 1.);
        EXPECT_EQ(pt.y, 1.);
    }

    vertices.clear();
    vertices.push_back(Vertex(0,0));
    vertices.push_back(Vertex(3.,0));
    vertices.push_back(Vertex(0,3.));
    vertices_p.clear();
    for (auto &v : vertices) vertices_p.push_back(&v);

    {
        Cell cell0(std::move(vertices_p));
        Point pt(cell0.centroid());
        EXPECT_EQ(pt.x, 1.);
        EXPECT_EQ(pt.y, 1.);
    }

}
// ------------------------------------------------------------
#ifdef BUILD_MODELE


// Bits...
const int GREENLAND = 1;
const int ANTARCTICA = 2;

bool clip(int zone, double lon0, double lat0, double lon1, double lat1)
{
    // Is it in Greenland range?
    if (zone & GREENLAND)
        if (SphericalClip::lonlat(-74., 59., -10., 87.5,
            lon0, lat0, lon1, lat1)) return true;

    // Is it in Antarctica range?
    if (zone & ANTARCTICA)
        if (lat0 <= -60. || lat1 <= -60) return true;

    // Not in range of either ice sheet, discard
    return false;
}


TEST_F(GridTest, hntr)
{
    // Create a grid with Hntr-style parameterization
    modele::GridSpec_Hntr hntr(modele::HntrGrid(4, 4, 0., 60.*45.));
    hntr.name = "hntr";
    hntr.spherical_clip = std::bind(&clip, GREENLAND|ANTARCTICA, _1, _2, _3, _4);
    hntr.pole_caps = false;
    hntr.points_in_side = 1;
    hntr.eq_rad = 1.;

    // Create a grid with traditional lon/lat parameterization
    GridSpec_LonLat ll;
    ll.name = "ll";
    ll.spherical_clip = hntr.spherical_clip;
    ll.latb = std::vector<double>{-90., -45., 0., 45., 90.};
    ll.lonb = std::vector<double>{0., 90., 180., 270., 360.};
    ll.indexing = ibmisc::Indexing(
        {"lon", "lat"}, {0,0}, {ll.nlon(), ll.nlat()}, {1,0});  // col major
    ll.south_pole = hntr.pole_caps;
    ll.north_pole = hntr.pole_caps;
    ll.points_in_side = hntr.points_in_side;
    ll.eq_rad = hntr.eq_rad;

    // Check the two grids are equal.
    Grid_LonLat ll_grid, hntr_grid;
    ll.make_grid(ll_grid);
    hntr.make_grid(hntr_grid);
    expect_eq(ll_grid, hntr_grid);

}
#endif // BUILD_MODELE
// ------------------------------------------------------------
int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}

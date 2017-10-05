#include <limits>
#include <blitz/array.h>
#include <ibmisc/netcdf.hpp>
#include <ibmisc/blitz.hpp>
#include <icebin/GCMRegridder.hpp>
#include <icebin/modele/z1qx1n_bs1.hpp>
#include <icebin/modele/hntr.hpp>
#include <icebin/modele/grids.hpp>

/* Demonstration program that takes a PISM ice sheet and puts it on a
   grid useful to Gary Russell's TOPO-generating code. */

using namespace ibmisc;
using namespace blitz;

namespace icebin {
namespace modele {

static double const NaN = std::numeric_limits<double>::quiet_NaN();

std::unique_ptr<GCMRegridder_Standard> load_AI_regridder(
    std::string const &gridA_fname,
    std::string const &gridA_vname,
    std::vector<std::tuple<std::string,std::string>> const &overlap_fnames)    // sheet name, fname
{
    std::unique_ptr<GCMRegridder_Standard> gcm_regridder(new GCMRegridder_Standard());

    // Read gridA
    printf("Opening %s\n", gridA_fname.c_str());
    NcIO ncio(gridA_fname, netCDF::NcFile::read);
    std::unique_ptr<Grid> gridA(new Grid); //new_grid(ncio, gridA_vname));
    gridA->ncio(ncio, gridA_vname);
    ncio.close();
    printf("Closing %s\n", gridA_fname.c_str());

    // Initialize the gcm_regridder with gridA
    std::vector<double> hcdefs {0};    // Dummy
    int nhc = hcdefs.size();
    gcm_regridder->init(
        std::move(gridA),
        std::move(hcdefs),
        Indexing({"A", "HC"}, {0,0}, {gridA->ndata(), nhc}, {1,0}),
        true);    // correctA; not used directly

    // Add gridI and exgrid for each ice sheet
    for (auto const &tuple : overlap_fnames) {
        auto const &sheet_name(std::get<0>(tuple));
        auto const &overlap_fname(std::get<1>(tuple));

        NcIO ncio(overlap_fname, netCDF::NcFile::read);

        std::unique_ptr<Grid> gridI(new Grid);//new_grid(ncio, "gridI"));
        gridI->ncio(ncio, "gridI");

        std::unique_ptr<Grid> exgrid(new_grid(ncio, "exgrid"));
        exgrid->ncio(ncio, "exgrid");

        ncio.close();

        blitz::Array<double,1> elevI(ibmisc::const_array<double,1>(blitz::shape(gridI->ndata()), 0));

        auto sheet(new_ice_regridder(gridI->parameterization));
        sheet->init(sheet_name, std::move(gridI), std::move(exgrid),
            InterpStyle::Z_INTERP, elevI);
        gcm_regridder->add_sheet(std::move(sheet));

    }

    printf("END load_AI_regridder()\n");
    return gcm_regridder;
}


template<class TypeT>
blitz::Array<TypeT,1> load_pism(netCDF::NcGroup *nc, std::string const &vname)
{
    // Variables have NetCDF dimensions (time, y, x)
    int ntime = nc->getDim("time").getSize();
    int ny = nc->getDim("y").getSize();
    int nx = nc->getDim("x").getSize();

    std::vector<size_t> startp {ntime-1,0,0};
    std::vector<size_t> countp {1,ny,nx};
    auto ncvar(nc->getVar(vname));

    blitz::Array<TypeT, 1> val(ny*nx);
    ncvar.getVar(startp, countp, val.data());

    return val;
}


void greenland_2m(
    std::string const &pism_fname,
    GCMRegridder &gcm_regridder,
    ArrayBundle<double,2> const &obundle)
{


    // ---------------------- Read PISM values
    printf("BEGIN load PISM: %s\n", pism_fname.c_str());
    NcIO ncio(pism_fname, netCDF::NcFile::read);

    auto maskI1(load_pism<uint8_t>(ncio.nc, "mask"));
    auto thkI1(load_pism<double>(ncio.nc, "thk"));
    auto topgI1(load_pism<double>(ncio.nc, "topg"));
    auto usurfI1(load_pism<double>(ncio.nc, "usurf"));

    ncio.close();
    printf("END load PISM: %s\n");

    // Different values for the mask
    const int MASK_UNKNOWN          = -1;
    const int MASK_ICE_FREE_BEDROCK = 0;
    const int MASK_GROUNDED         = 2;
    const int MASK_FLOATING         = 3;
    const int MASK_ICE_FREE_OCEAN   = 4;


    // ----------------------- Regrid to ModelE grid

    // Obtain regridding matrix
    SparseSetT dimA,dimI;
    RegridMatrices::Params params(true, false, {0.,0.,0.});
    //    params.scale = true;
    //    params.correctA = true;
    //    params.sigma = {0,0,0};
    //    params.conserve = true;
    IceRegridder *ice_regridder = gcm_regridder.ice_regridder("greenland");
    RegridMatrices rm(gcm_regridder.regrid_matrices(ice_regridder->name()));
    auto AvI(rm.matrix("AvI", {&dimA, &dimI}, params));

printf("dimA: ndense=%d nsparse%ld\n", dimA.dense_extent(), dimA.sparse_extent());
printf("dimI: ndense=%d nsparse%ld\n", dimI.dense_extent(), dimI.sparse_extent());

    // Copy PISM values into single regriddable densified II_d
    const int FCONT = 0;
    const int FGICE = 1;
    const int THK = 2;
    const int TOPG = 3;    // Solid ground topograpy
    const int NVAR = 4;

    blitz::Array<double,2> II_d(NVAR, dimI.dense_extent());
    for (int i_d=0; i_d<dimI.dense_extent(); ++i_d) {
        auto i_s(dimI.to_sparse(i_d));
        auto mask(maskI1(i_s));
        II_d(FCONT, i_d) = (mask != MASK_UNKNOWN && mask != MASK_ICE_FREE_OCEAN ? 1.0 : 0.0);
        II_d(FGICE, i_d) = (mask == MASK_GROUNDED || mask == MASK_FLOATING ? 1.0 : 0.0);
        II_d(THK, i_d) = thkI1(i_s);
        II_d(TOPG,  i_d) = topgI1(i_s);
    }

    TmpAlloc tmp;    // Must live for scope of AA_d
    auto AA_d(AvI->apply(II_d, NaN, false, tmp));
    II_d.free();

    // ---------------- Store in output
    // Extract output variables
    auto WT2_1(reshape1(obundle.array("WT")));
    auto FCONT2_1(reshape1(obundle.array("FCONT")));
    auto FOCEN2_1(reshape1(obundle.array("FOCEN")));
    auto FGICE_SHEET2_1(reshape1(obundle.array("FGICE_SHEET")));
    auto FGICE2_1(reshape1(obundle.array("FGICE")));
    auto dZGIC2_1(reshape1(obundle.array("dZGIC")));
    auto ZSOLD2_1(reshape1(obundle.array("ZSOLD")));
    auto ZSOLG2_1(reshape1(obundle.array("ZSOLG")));

    ibmisc::Proj_LL2XY proj(ice_regridder->gridI->sproj);
    for (int i_d=0; i_d<dimA.dense_extent(); ++i_d) {
        auto i_s(dimA.to_sparse(i_d));

        double area = gcm_regridder.gridA->cells.at(i_s)->proj_area(&proj);
        //double area = gcm_regridder.gridA->cells.at(i_s)->native_area;
        WT2_1(i_s) += AvI->wM(i_d) / area;

        FCONT2_1(i_s) = AA_d(FCONT, i_d);
        FOCEN2_1(i_s)= 1. - FCONT2_1(i_s);
        FGICE_SHEET2_1(i_s) = AA_d(FGICE, i_d);
        FGICE2_1(i_s) += FGICE_SHEET2_1(i_s);
        dZGIC2_1(i_s) = AA_d(THK, i_d);
        ZSOLD2_1(i_s) = std::max(0., AA_d(THK, i_d) + AA_d(TOPG, i_d));
        ZSOLG2_1(i_s) = AA_d(TOPG, i_d);
    }
}


ArrayBundle<double,2> make_bundle()
{
    ArrayBundle<double, 2> bundle;

    bundle.add("WT", {});            // Fortran Arrays
    bundle.add("FCONT", {});
    bundle.add("FOCEN", {});
    bundle.add("FGICE_SHEET", {});
    bundle.add("FGICE", {});
    bundle.add("dZGIC", {});
    bundle.add("ZSOLD", {});
    bundle.add("ZSOLG", {});

    return bundle;
}

void do_main()
{
    auto bundle2(make_bundle());
    bundle2.allocate(blitz::shape(IM2, JM2), {"im2", "jm2"});

    for (size_t i=0; i<bundle2.data.size(); ++i) {
        bundle2.data[i].arr = NaN;
    }
    bundle2.at("WT").arr = 0;


    // Load the IceBin Regridder
    std::string const gridA_fname = "modele_ll_g2mx2m.nc";
    std::string const overlap_fname = "modele_ll_g2mx2m-sr_g20_searise.nc";
    std::string const pism_fname = "/home2/rpfische/f15/modelE/init_cond/pism/std-greenland/ex_g20km_10ka.nc";

    std::unique_ptr<GCMRegridder_Standard> gcm_regridder(load_AI_regridder(
        overlap_fname, "gridA",
        {std::make_tuple("greenland", overlap_fname)}));

    icebin::modele::greenland_2m(pism_fname, *gcm_regridder, bundle2);

    // Regrid to 1/2-degree --- because the 2-minute file is too slow for ncview
    auto bundleH(make_bundle());
    bundleH.allocate(blitz::shape(IMH, JMH), {"im", "jm"});
    Hntr hntr(g2mx2m, g1qx1, NaN);
    for (size_t i=0; i<bundleH.data.size(); ++i) {
        bundleH.data[i].arr.reference(hntr.regrid(bundle2.at("WT").arr, bundle2.data[i].arr));
    }

    printf("BEGIN writing output\n");
    {NcIO ncio("greenland2m.nc", 'w');
        bundleH.ncio(ncio, {}, false, "", "double");
    }
    printf("END writing output\n");

}

}}

int main(int argc, char **argv)
{
    icebin::modele::do_main();
}

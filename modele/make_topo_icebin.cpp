#include <limits>
#include <blitz/array.h>
#include <ibmisc/netcdf.hpp>
#include <ibmisc/blitz.hpp>
#include <icebin/GCMRegridder.hpp>
#include <icebin/GCMCoupler.hpp>
#include <icebin/modele/z1qx1n_bs1.hpp>
#include <icebin/modele/hntr.hpp>
#include <icebin/modele/grids.hpp>
#include <everytrace.h>

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

        auto sheet(new_ice_regridder(gridI->parameterization));
        sheet->init(sheet_name, std::move(gridI), std::move(exgrid),
            InterpStyle::Z_INTERP);
        gcm_regridder->add_sheet(std::move(sheet));

    }

    printf("END load_AI_regridder()\n");
    return gcm_regridder;
}

// -------------------------------------------------------------------

ArrayBundle<uint8_t,1> pism_outputs_bundle_uint8_t()
{
    ArrayBundle<uint8_t,1> bundle;
    bundle.add("mask", {
        "description", "ice_free_bedrock grounded_ice floating_ice ice_free_ocean",
        "units", "m",
        "source", "PISM",
    });
    bundle.add("topg", {
        "description", "bedrock surface elevation",
        "units", "m",
        "source", "PISM",
    });
    bundle.add("usurf", {
        "description", "ice upper surface elevation",
        "units", "m",
        "source", "PISM",
    });
    return bundle;
}

ArrayBundle<double,1> pism_outputs_bundle_double()
{
    ArrayBundle<double,1> bundle;
    bundle.add("thk", {
        "description", "land ice thickness",
        "units", "m",
        "source", "PISM",
    });
    bundle.add("topg", {
        "description", "bedrock surface elevation",
        "units", "m",
        "source", "PISM",
    });
    bundle.add("usurf", {
        "description", "ice upper surface elevation",
        "units", "m",
        "source", "PISM",
    });
    return bundle;
}

struct PISMOutputs {
    ibmisc::ArrayBundle<uint8_t,1> bundle_uint8_t;
    ibmisc::ArrayBundle<double,1> bundle_double;

    blitz::Array<uint8_t,1> &mask;
    blitz::Array<double,1> &thk;
    blitz::Array<double,1> &topg;
    blitz::Array<double,1> &usurf;

    PISMOutputs();

    // Read the latest */
    void ncio_read(NcIO &nc);
};

PISMOutputs::PISMOutputs() :
    bundle_uint8_t(pism_outputs_bundle_uint8_t()),
    bundle_double(pism_outputs_bundle_double()),
    mask(bundle_uint8_t.array("mask")),
    thk(bundle_double.array("thk")),
    topg(bundle_double.array("topg")),
    usurf(bundle_double.array("usurf"))
{}


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

void PISMOutputs::ncio_read(NcIO &ncio)
{
    mask.reference((load_pism<uint8_t>(ncio.nc, "mask")));
    thk.reference((load_pism<double>(ncio.nc, "thk")));
    topg.reference((load_pism<double>(ncio.nc, "topg")));
    usurf.reference((load_pism<double>(ncio.nc, "usurf")));
}
// ---------------------------------------------------------------------
#if 0
template<class TypeT, int RANK>
ArrayBundle<TypeT,1> reshape1(
    ArrayBundle<TypeT, RANK> &bundle,
    int lbound = 0,
    std::array<std::string,1> const &sdims = {""})
{
    ArrayBundle<TypeT,1> bundle1;
    for (size_t i=0; i<bundle.index.size(); ++i) {
        auto &meta(bundle.data[i]);
        typename ArrayBundle<TypeT,1>::Meta meta1(
            meta.name, ibmisc::reshape1(meta.arr, lbound),
            blitz::shape(meta1.arr.extent(0)), sdims, meta.attr);
        bundle1.index.insert(meta.name);
        bundle1.data.push_back(meta1);
    }
    return bundle1;
}
#endif

void write_regrid_nc(blitz::Array<double,2> &AA_d, SparseSetT &dimA)
{
    ArrayBundle<double,2> b2;
    b2.add("thk", {});
    b2.add("cont", {});
    b2.add("ice", {});
    b2.add("topg", {});
    b2.allocate(shape(IM,JM), {"im", "jm"}, true, fortranArray);

    ArrayBundle<double,1> b1(reshape1(b2,0));

    for (int j=0; j<4; ++j) b1.data[j].arr = NaN;

    for (int i_d=0; i_d<dimA.dense_extent(); ++i_d) {
        auto i_s(dimA.to_sparse(i_d));
        for (int j=0; j<4; ++j) {
            b1.data[j].arr(i_s) = AA_d(j,i_d);
        }
    }

    NcIO ncio("AA_d.nc", 'w');
    b2.ncio(ncio, {}, false, "", "double", fortranArray);
}


void pism_replace_greenland(
    TopoOutputs<2> &tout,
    PISMOutputs &pout,
    GCMRegridder_Standard &gcm_regridder)
{
    // ----------------------- Regrid to ModelE grid

    // Obtain regridding matrix
    SparseSetT dimA,dimI;
    RegridMatrices::Params params(true, false, {0.,0.,0.});
    //    params.scale = true;
    //    params.correctA = false;
    //    params.sigma = {0,0,0};
    std::string const sheet = "greenland";
    auto sheet_ix = gcm_regridder.ice_regridders().index.at(sheet);
    IceRegridder *ice_regridder = &*gcm_regridder.ice_regridders()[sheet_ix];
    blitz::Array<double,1> elevmaskI(ibmisc::const_array<double,1>(blitz::shape(ice_regridder->gridI->ndata()), 0));
    RegridMatrices rm(gcm_regridder.regrid_matrices(sheet_ix, elevmaskI));
    auto AvI(rm.matrix("AvI", {&dimA, &dimI}, params));

printf("dimA: ndense=%d nsparse%ld\n", dimA.dense_extent(), dimA.sparse_extent());
printf("dimI: ndense=%d nsparse%ld\n", dimI.dense_extent(), dimI.sparse_extent());

    // Copy PISM values into single regriddable densified II_d
    const int THK = 0;
    const int GRND = 1;
    const int ICE = 2;
    const int TOPG = 3;
    const int NVAR = 4;

    blitz::Array<double,2> II_d(NVAR, dimI.dense_extent());
    II_d = 0;

    for (int i_d=0; i_d<dimI.dense_extent(); ++i_d) {
        auto i_s(dimI.to_sparse(i_d));
        //double areaI = gcm_regridder.gridI->cells.at(i_s).native_area;
        auto mask(pout.mask(i_s));

        switch(mask) {
            case IceMask::GROUNDED_ICE:
            case IceMask::FLOATING_ICE:
                II_d(ICE,i_d) = 1.0;
            break;
            case IceMask::ICE_FREE_BEDROCK:
                II_d(GRND,i_d) = 1.0;
            break;
            case IceMask::ICE_FREE_OCEAN:
            break;
            default:
                (*icebin_error)(-1,
                    "Illegal mask[%ld] = %d", i_s, mask);
        }
        II_d(TOPG,i_d) = pout.topg(i_s);
        II_d(THK,i_d) = pout.thk(i_s);
    }

    TmpAlloc tmp;    // Must live for scope of AA_d
    auto AA_d(AvI->apply(II_d, NaN, false, tmp));
    II_d.free();

    //write_regrid_nc(AA_d, dimA);    // Testing / Debugging

    TopoOutputs<1> tout1(reshape1(tout.bundle));
    ibmisc::Proj_LL2XY proj(ice_regridder->gridI->sproj);
    for (int i_d=0; i_d<dimA.dense_extent(); ++i_d) {
        auto i_s(dimA.to_sparse(i_d));

        // Original Land surface types
        double const fgice0 = tout1.FGICE(i_s);
        double const fgrnd0 = tout1.FGRND(i_s);
        double const flake0 = tout1.FLAKE(i_s);

        // Land surface types before the ocean was roudned
        double const fcont_lg0 = 1.-tout1.FOCENF(i_s);
        double const fgice_lg0 = fcont_lg0 * fgice0;    // Land has already been rounded
        double const fgrnd_lg0 = fcont_lg0 * fgrnd0;
        double const flake_lg0 = fcont_lg0 * flake0;

        // Portion of cell covered by PISM continent
        double Aportion = AvI->wM(i_d) /
            gcm_regridder.gridA->cells.at(i_s)->proj_area(&proj);
        double const fgice_gr0 = Aportion * AA_d(ICE,i_d);
        double const fgrnd_gr0 = Aportion * AA_d(GRND,i_d);
        double const flake_gr0 = 0;
        double const fcont_gr0 = fgice_gr0 + fgrnd_gr0 + flake_gr0;

        double const fcont_all0 = fcont_lg0 + fcont_gr0;

        // Only layer ourselves in if we have new land to add...
        if (fcont_gr0 == 0) continue;

        // Remove newly discovered land area from ocean
        tout1.FOCENF(i_s) -= fcont_gr0;

        // Round the ocean...
        double const focean = std::round(tout1.FOCENF(i_s));

        if (focean == 0.0) {
            // Land grid cell
            double const by_fcont_all0 = 1. / fcont_all0;

            // Adjust FGICE & FGRND because the full cell is now ice or bare land
            tout1.FOCEAN(i_s) = focean;
            tout1.FGICE_greenland(i_s) = fgice_gr0 * by_fcont_all0;
            tout1.FGICE(i_s) = fgice_lg0 * by_fcont_all0;
            tout1.FGRND(i_s) = (fgrnd_lg0 + fgrnd_gr0) * by_fcont_all0;
            tout1.FLAKE(i_s) = (flake_lg0 + flake_gr0) * by_fcont_all0;
            tout1.dZOCEN(i_s) = 0.;

            {double const lg_portion = fcont_lg0 / (fcont_lg0 + fcont_gr0);
            double const gr_portion = 1. - lg_portion;
                tout1.ZATMO(i_s) = tout1.ZATMO(i_s) * lg_portion +
                    (AA_d(THK,i_d) + AA_d(TOPG,i_d)) * gr_portion;
                tout1.ZSOLDG(i_s) = tout1.ZSOLDG(i_s) * lg_portion +
                    AA_d(TOPG,i_d) * gr_portion;
            }

            if (fgice_lg0 + fgice_gr0 == 0) {
                tout1.dZGICE(i_s) = 0;
            } else {
                double const lg_portion = fgice_lg0 / (fgice_lg0 + fgice_gr0);
                double const gr_portion = 1. - lg_portion;
                tout1.dZGICE(i_s) = tout1.dZGICE(i_s) * lg_portion +
                    AA_d(THK,i_d) * gr_portion;
            }

            // Is this needed???
            tout1.ZSGLO(i_s) = tout1.ZSOLDG(i_s);
            tout1.ZSGHI(i_s) = tout1.ZSOLDG(i_s);
            tout1.ZGRND(i_s) = tout1.ZATMO(i_s);
        } else {
            if (tout1.FOCEAN(i_s) == 0) {
                fprintf(stderr, "Adding ocean at %ld, is this correct?", i_s);
                tout1.FOCEAN(i_s) = 2.0;
            }

            tout1.FOCEAN(i_s) = focean;
            tout1.FGICE_greenland(i_s) = 0;
            tout1.FGICE(i_s) = 0;
            tout1.FGRND(i_s) = 0;
            tout1.FLAKE(i_s) = 0;
        }


    }
}








void do_main()
{
    std::string const TOPO_nogr_fname = "z1qx1n_bs1-nogr.nc";
    std::string const overlap_fname = "modele_ll_g1qx1-sr_g20_searise.nc";
    std::string const pism_fname = "/home2/rpfische/f15/modelE/init_cond/pism/std-greenland/ex_g20km_10ka.nc";


    // Load outputs from Greenland-free TOPO generation
    TopoOutputs<2> tout(topo_outputs_bundle2(true));
    {NcIO ncio(TOPO_nogr_fname, 'r');
        tout.bundle.ncio(ncio,
            {"FOCEAN", "FLAKE", "FGRND", "FGICE",
            "ZATMO", "dZOCEN", "dZLAKE", "dZGICE", "ZSOLDG", "ZSGLO",
            "ZLAKE", "ZGRND", "ZSGHI", "FOCENF"},
            false, "", "double", blitz::FortranArray<2>());
    }

    // LoadPISM Values
    printf("BEGIN load PISM: %s\n", pism_fname.c_str());
    PISMOutputs pout;
    {NcIO ncio(pism_fname, netCDF::NcFile::read);
        pout.ncio_read(ncio);
    }

    // Load the IceBin Regridder
    std::unique_ptr<GCMRegridder_Standard> gcm_regridder(load_AI_regridder(
        overlap_fname, "gridA",
        {std::make_tuple("greenland", overlap_fname)}));

    pism_replace_greenland(tout, pout, *gcm_regridder);

    printf("BEGIN writing output\n");
    {NcIO ncio("out-pism.nc", 'w');
        tout.bundle.ncio(ncio, {}, false, "", "double");
    }
    printf("END writing output\n");

}

}}

int main(int argc, char **argv)
{
    everytrace_init();
    icebin::modele::do_main();
}

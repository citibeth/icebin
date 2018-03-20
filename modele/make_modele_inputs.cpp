#include <icebin/modele/grids.hpp>
#include <icebin/modele/hntr.hpp>


static std::vector<std::string> const itopo_vars
    {"FOCEAN", "FLAKE", "FGRND", "FGICE", "ZATMO", "ZLAKE"};
static std::vector<std::string> const otopo_vars
    {"focean", "flake", "fgrnd", "fgice", "zatmo", "hlake"};


/** Reads relevant variables from Ocean-grid TOPO file, regrids them,
    and writes to Atmosphere-grid TOPO file.  Also renames in the
    process, to conform to names expected by ModelE. */
void topoo_to_topoa(NcIO &ncout, std::string const &TOPOO_nc, std::string const &TOPOA_nc)
{

    // No new ice sheets to merge.  We just need to regrid TOPOO
    // and rename...
    NcIO ncin(TOPOO_nc, 'r');

    // Define output variables
    auto odims(get_or_add_dims(ncout, {"jm", "im"}, {hspecA.jm, hspecA.im}));
#if 0
    for (int i=0; i<otopo_vars.size(); ++i) {
        get_or_add_var(ncout, otopo_vars[i], "double", odims);
    }
#endif

    // Create the AvO matrix (Eigen format)
    Hntr hntr_AvO(17.17, args.hspecA, args.hspecO);
    TupleListT<2> AvO_tp;
    hntrAvO.scaled_regrid_matrix(spsparse::accum::ref(AvO_tp));
    EigenSparseMatrixT AvO_e(hntrA.size(), hntrO.size());
    AvO_e.setFromTriplets(AvO_tp.begin(), AvO_tp.end());

    for (int i=0; i<itopo_vars.size(); ++i) {
        // Read on O grid
        blitz::Array<double,2> valO2(hspecO.jm, hspecO.im);
        std::vector<std::pair<std::string, NcAttValue>> atts;
        get_or_put_all_atts(
            ncio_blitz(ncin, valO2, itopo_vars[i], ""),
            ncin.rw, atts);
        blitz::Array<double,1> valO(reshape1(valO2));

        // Regrid to A grid
        auto &valA2(ncout.tmp.make<blitz::Array<double,2>>(hspecA.jm, hspecA.im));
        blitz::Array<double,1> valA(reshape1(valA2));
        map_eigen_colvector(valA) = AvO_e * map_eigen_colvector(valO);

        // Write it out, transferring attributes from original variable
        get_or_put_all_atts(
            ncio_blitz(ncout, valA2, otopo_vars[i], "double", odims),
            ncout.rw, attss);
    }
}

void add_fhc_legacy(NcIO &ncout, std::string const &global_ec_nc, std::string const &etopo1_ice_nc)
{
    blitz::Array<double,1> zatmo(hspecA.jm, hspecA.im);
    {auto fname(files.locate(TOPOA));
        NcIO ncio(fname, 'r');
        ncio_blitz(ncio, zatmo, "ZATMO", "double", {});
    }

    auto all(blitz::Range::all());
    fhc(ec_base, all, all) = zatmo;
    elevE(ec_base, all, all) = elevA
    underice(ec_base, all, all) = UI_NOTHING;
}


void add_fhc_sealand()
{

TODO: Read hspecA, hspecO, hspecI from global_ec file

    blitz::Array<double,1> elevA(hspecA.jm * hspecA.im);
    {
        // Read in ice extent and elevation
        blitz::Array<int16_t,2> fgiceI(hspecI.jm, hspecI.im);    // 0 or 1
        blitz::Array<int16_t,2> elevI(hspecI.jm, hspecI.im);
        {auto fname(files.locate(args.nc_fname));
            NcIO ncio(fname, 'r');
            ncio_blitz(ncio, fgiceI, args.fgiceI_vname, "short", {});
            ncio_blitz(ncio, elevI, args.elevI_vname, "short", {});
        }

        // Read AvI
        blitz::unique_ptr<linear::Weighted> AvI;
        {NcIO ncio(files.locate("global_ec-standard.nc"), 'r');
            AvI = nc_read_weighted(ncio.nc, "AvI");
        }

        // Produce elevA: Average elevation of ICE-COVERED
        // region
        elevA = NaN;
        AvI->apply_M(reshape1(elevI), elevA);
    }


    auto all(blitz::Range::all());
    fhc(ec_base, all, all) = 1;
    elevE(ec_base, all, all) = elevA
    underice(ec_base, all, all) = UI_NOTHING;
}

void add_fhc_globalec()
{
}    







    {
        // Allocate arrays
        blitz::Array<int16_t,2> fgiceI(hspecI.jm, hspecI.im);    // 0 or 1
        blitz::Array<int16_t,2> elevI(hspecI.jm, hspecI.im);

        // Read in ice extent and elevation
        {auto fname(files.locate(args.nc_fname));
            NcIO ncio(fname, 'r');
            ncio_blitz(ncio, fgiceI, args.fgiceI_vname, "short", {});
            ncio_blitz(ncio, elevI, args.elevI_vname, "short", {});
        }
        // -----------------------------------------

        // Generate fgiceO
        auto wtI(const_array(fgiceI.shape(), 1.0));
        Hntr hntrOvI(17.17, args.hspecO, args.hspecI);
        hntrOvI.regrid(wtI, fgiceI, fgiceO);

#if 0
        // Generate elevA
        // (This is a shortcut, compared to creating a GCMRegridder for A).
        HntrSpec const hspecA(make_hntrA(args.hspecO));
        blitz::Array<double,2> elevA(
        Hntr hntrAvI(17.17, hspecA, args.hspecI);
        hntrAvI.regrid(fgiceI, elevI, elevA);
#endif

    }



}

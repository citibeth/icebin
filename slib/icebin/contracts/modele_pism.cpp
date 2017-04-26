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

#include <mpi.h>        // Must be first
#include <limits>
#include <pism/base/enthalpyConverter.hh>
#include <ibmisc/VarTransformer.hpp>
#include <icebin/contracts/contracts.hpp>
#include <icebin/modele/GCMCoupler_ModelE.hpp>
#include <icebin/pism/IceCoupler_PISM.hpp>

using namespace ibmisc;
using namespace icebin::modele;

// --------------------------------------------------------
namespace icebin {
namespace contracts {

// Aliases
static auto const &UNIT(VarTransformer::UNIT);
static auto const &INPUT(IceCoupler::INPUT);
static auto const &OUTPUT(IceCoupler::OUTPUT);
static double const nan = std::numeric_limits<double>::quiet_NaN();

static void _reconstruct_ice_ivalsI(
    GCMCoupler_ModelE const *gcm_coupler,
    gpism::IceCoupler_PISM const *ice_coupler,
    blitz::Array<double,2> &ice_ivalsI,
    double dt)
{
    // ------------------ Pick relevant variables out of ice input & output

    // --------- Inputs of this Computation
    // surface_senth: State of top of PISM ice sheet [J kg-1]
    blitz::Array<double,1> surface_senth(ice_coupler->ice_ovalsI(
        ice_coupler->contract[OUTPUT].index.at("surface_senth"),
        blitz::Range::all()));
    // deltah: Change in enthalpy to effect in PISM by B.C.
    blitz::Array<double,1> deltah(ice_ivalsI(
        ice_coupler->contract[INPUT].index.at("deltah"),
        blitz::Range::all()));

    blitz::Array<double,1> pism_surface_senth(ice_ivalsI(
        ice_coupler->contract[INPUT].index.at("pism_surface_senth"),
        blitz::Range::all()));


    // --------- Outputs of this Computation
    blitz::Array<double,1> bc_senth(ice_ivalsI(
        ice_coupler->contract[INPUT].index.at("bc_senth"),
        blitz::Range::all()));
    blitz::Array<double,1> bc_temp(ice_ivalsI(
        ice_coupler->contract[INPUT].index.at("bc_temp"),
        blitz::Range::all()));
    blitz::Array<double,1> bc_watercontent(ice_ivalsI(
        ice_coupler->contract[INPUT].index.at("bc_watercontent"),
        blitz::Range::all()));
    blitz::Array<double,1> pism_surface_temp(ice_ivalsI(
        ice_coupler->contract[INPUT].index.at("pism_surface_temp"),
        blitz::Range::all()));


    // kappa: see modelE/model/landice/lisnowsubs.F90
    // Layer 0 = ice borrowed from GCM
    // Layer 1 = top layer of PISM
    double const RHOI = gcm_coupler->gcm_constants.get_as("constant.rhoi", "kg m-3");
    double const SHI = gcm_coupler->gcm_constants.get_as("constant.shi", "J kg-1 K-1");
    double const TF = gcm_coupler->gcm_constants.get_as("constant.tf", "K");
    double const STIEGLITZ_8B = gcm_coupler->gcm_constants.get_as("constant.stieglitz_8b", "W m-1 K-1");

    // PISM will "see" a layer of PISM-style ice, with a certain enthalpy.
    // Make sure that this layer produces ~deltah of heat flux

    // Pretend both layers have density RHOI; because that is what PISM will see
    // [3.2217] [kg m-3] [kg m-3] = [W m-1 K-1]  (Stieglitz eq 8b)
    double const ksn = STIEGLITZ_8B * (RHOI*RHOI);    // [W m-1 K-1]
    // Pretend that both layers are 40m; the standard depth of a PISM layer
    double const dz = 40;    // [m]  TODO: Find where this is in PISM
    // Constant for Fourier's law:
    //     downward flux (layer 0->1) [W m-2] = kappa (H0 - H1)
    //     where H=specific enth
    // [J s-1 m-1 K-1] [J-1 kg K] [m-1] = [kg s-1 m-2]
    double const kappa = ksn / (SHI * dz);   // [kg m-2 s-1]

    // Get a PISM Enthalpy Converter
    pism::EnthalpyConverter enth(*ice_coupler->pism_config());

    for (int i=0; i<ice_coupler->nI(); ++i) {
        double const H1 = surface_senth(i);    // [J kg-1]
        double const Q = deltah(i) / dt;          // [W m-2]

        // Specific enthalpy for the boundary condition
        // Q/kappa: [W m-2 kg-1 s m^2] = [J kg-1]
        double const H0 = (Q/kappa) + H1;        // [J kg-1]
        bc_senth(i) = H0;

        // Pressure at top of PISM ice sheet; by convention, 0
        double const P = 0;

        // Convert specific enthalpy to T and water content
        bc_temp(i) = enth.temperature(H0, P);
        bc_watercontent(i) = enth.water_fraction(H0, P);
        pism_surface_senth(i) = surface_senth(i);

        pism_surface_temp(i) = enth.temperature(pism_surface_senth(i), P) - TF;

    }
}

/** GCM-specific contract */
void setup_modele_pism(GCMCoupler const *_gcm_coupler, IceCoupler *_ice_coupler)
{
    // Get arguments we need from coupler
    auto gcm_coupler(dynamic_cast<modele::GCMCoupler_ModelE const *>(_gcm_coupler));
    auto ice_coupler(dynamic_cast<icebin::gpism::IceCoupler_PISM *>(_ice_coupler));

    printf("BEGIN setup_modele_pism()\n");

    ice_coupler->reconstruct_ice_ivalsI = std::bind(
        &_reconstruct_ice_ivalsI,
        gcm_coupler, ice_coupler, std::placeholders::_1, std::placeholders::_2);

    // =========== Transfer  constants from Icebin's gcm_constants --> PISM
    // ice_coupler->transfer_constant(PISM_destionation, icebin_src, multiply_by)
    // (see IceCoupler_PISM.cpp)
    ice_coupler->transfer_constant("constants.standard_gravity", "constant.grav");
    ice_coupler->transfer_constant("constants.ice.beta_Clausius_Clapeyron", "seaice.dtdp", -1.0);
    ice_coupler->transfer_constant("constants.fresh_water.melting_point_temperature", "constant.tf");
    ice_coupler->transfer_constant("constants.fresh_water.latent_heat_of_fusion", "constant.lhm");
    ice_coupler->transfer_constant("constants.fresh_water.specific_heat_capacity", "constant.shw");

    ice_coupler->transfer_constant("constants.ice.density", "constant.rhoi");
    ice_coupler->transfer_constant("constants.ice.thermal_conductivity", "seaice.alami0");
    ice_coupler->transfer_constant("constants.ice.specific_heat_capacity", "constant.shi");
    ice_coupler->transfer_constant("constants.fresh_water.density", "constant.rhow");
    ice_coupler->transfer_constant("constants.sea_water.density", "constant.rhows");
    ice_coupler->transfer_constant("constants.ideal_gas_constant", "constant.gasc");

    // To set this, see (in ModelE): Function SHCGS in ocnfuntab.f is
    // used for the Russell ocean. I. The simple models use SHW=4185.
    // This probably doesn't matter much at this point (May 2014)
    //      ice_coupler->transfer_constant("sea_water_specific_heat_capacity", "");

    /* The following constants were not transferred 
    pism_config:fill_value = -2e9;
    pism_config:fill_value_doc = "_FillValue used when saving diagnostic quantities";
    */

    // In PISM and ModelE Clausius-Clapeyron equation, surfce_pressure is the DIFFERENCE
    // from 1atm.  Thus, surface_pressure=0 implies the ice sheet existing at 1atm
    ice_coupler->set_constant("surface.pressure", 0, "Pa");       // Match ModelE thermodynam

    // No need to set enthalpy_reference_temperature.  pism::EnthalpyConverter is used (below)
    // to convert enthalpy values between ModelE and PISM.
    // ice_coupler->transfer_constant("enthalpy_reference_temperature", "enthalpy_reference_temperature");


    // ============ GCM -> Ice
    VarSet &ice_input(ice_coupler->contract[IceCoupler::INPUT]);

    // ------ Decide on the coupling contract for this ice sheet
    ice_input.add("massxfer", 0., "kg m-2 s-1", 0,
        "Mass of ice being transferred Stieglitz --> Icebin");
    ice_input.add("enthxfer", 0., "W m-2", 0,
        "Enthalpy of ice being transferred Stieglitz --> Icebin");
    ice_input.add("gcm_bottom_temp", 0., "degC", contracts::PRIVATE,
        "Temperature at bottom of GCM snow/firn model");


    // Temporary variables of the ice input: not bound to PISM variables.
    // deltah: Regridded from GCM output's deltah
    ice_input.add("deltah", 0., "J m-2", contracts::PRIVATE,
        "Change of enthalpy to apply to top layer in PISM");
    // bc_senth: Boundary condition, computed based on deltah
    ice_input.add("bc_senth", 0., "J kg-1", contracts::PRIVATE,
        "Enthalpy of the Dirichlet B.C. being applied to the top of the ice sheet");
    ice_input.add("pism_surface_senth", 0., "J kg-1", contracts::PRIVATE,
        "surface_senth output by PISM on the last timestep");

    // These variables are the actual boundary condition provided to PISM.
    // They are computed based on bc_senth.
    ice_input.add("bc_temp", 0., "K", 0,
        "Temperature of the Dirichlet B.C.");
    ice_input.add("bc_watercontent", 0., "1", 0,
        "Water content of the Dirichlet B.C.");
    ice_input.add("pism_surface_temp", 0., "degC", contracts::PRIVATE,
        "Temperature of pism_surface_senth");


    // Figure out the conversion between GCM and PISM enthalpy
    // ModelE's reference state is 1atm, 0C, 100% liquid water.
    // The enthalpy for that reference state would be the top end
    // of PISM's EnthalpyInterval.
    // NOTE: Pressure in PISM is RELATIVE to atmospheric pressure.
    //       Thus, p=0 is the correct to use at the top surface of
    //       the ice sheet (where ModelE operates).
    pism::EnthalpyConverter enth(*ice_coupler->pism_config());
    double const pressure = 0;
    double E_l;
    // getEnthalpyInterval() replaced with enthalpy_liquid()
    // https://github.com/pism/pism/commit/c820dfd8
    E_l = enth.enthalpy_liquid(pressure);
    double const enth_modele_to_pism = E_l;     // [J/kg]: Add to convert ModelE specific enthalpies [J/kg] to PISM specific enthalpies [J/kg]
    // NOTE: enth_modele_to_pism == 437000 [J/kg]
    printf("enth_modele_to_pism = %g\n", enth_modele_to_pism);

    bool ok = true;

    double const RHOW = gcm_coupler->gcm_constants.get_as("constant.rhow", "kg m-3");
    double const byRHOW = 1.0 / RHOW;


    // ------------- Convert the contract to a var transformer
    // ------------- of I <- E   (Ice <- GCM)
    {VarTransformer &vt(ice_coupler->var_trans_inE);
    vt.set_dims(
        ice_input.keys(),             // outputs
        gcm_coupler->gcm_outputsE.keys(),  // inputs
        gcm_coupler->scalars.keys());      // scalars

    ok = ok && vt.set("massxfer", "massxfer", "by_dt", 1.0);
    ok = ok && vt.set("enthxfer", "enthxfer", "by_dt", 1.0);
    ok = ok && vt.set("enthxfer", "massxfer", "by_dt", enth_modele_to_pism);
    ok = ok && vt.set("gcm_bottom_temp", "gcm_bottom_temp", UNIT, 1.0);
    ok = ok && vt.set("deltah", "deltah", UNIT, 1.0);
    }

    // ============== Ice -> GCM
    VarSet &ice_output(ice_coupler->contract[IceCoupler::OUTPUT]);
    auto &standard_names(ice_coupler->standard_names[IceCoupler::OUTPUT]);

    // All these outputs are on the ICE grid.
    // Icebin requires that all ice models return elev2, so that it can regrid in the vertical.

    standard_names["elevI"] =
    ice_output.add("ice_surface_elevation", nan, "m", contracts::INITIAL, "ice upper surface elevation");
    ice_output.add("ice_thickness", nan, "m", contracts::INITIAL, "thickness of ice");
    ice_output.add("bed_topography", nan, "m", contracts::INITIAL, "topography of bedrock");


    ice_output.add("mask", nan, "", contracts::INITIAL, "PISM land surface type");

    ice_output.add("surface_senth", nan, "J kg-1", contracts::INITIAL, "");

    ice_output.add("basal_frictional_heating", nan, "W m-2", 0, "");
    ice_output.add("strain_heating", nan, "W m-2", 0, "");
    ice_output.add("geothermal_flux", nan, "W m-2", 0, "");
    ice_output.add("upward_geothermal_flux", nan, "W m-2", 0, "");

    ice_output.add("calving.mass", nan, "kg m-2 s-1", 0, "");
    ice_output.add("calving.enth", nan, "W m-2", 0, "");
    ice_output.add("icebin_smb.mass", nan, "kg m-2 s-1", 0, "");
    ice_output.add("icebin_smb.enth", nan, "W m-2", 0, "");
    ice_output.add("pism_smb.mass", nan, "kg m-2 s-1", 0, "");
    ice_output.add("pism_smb.enth", nan, "W m-2", 0, "");

    // basal_runoff (GCM input) = melt_grounded + melt_floatig (PISM outputs)
    ice_output.add("melt_grounded.mass", nan, "kg m-2 s-1", 0, "");
    ice_output.add("melt_grounded.enth", nan, "W m-2", 0, "");
    ice_output.add("melt_floating.mass", nan, "kg m-2 s-1", 0, "");
    ice_output.add("melt_floating.enth", nan, "W m-2", 0, "");

    ice_output.add("internal_advection.mass", nan, "kg m-2 s-1", 0, "");
    ice_output.add("internal_advection.enth", nan, "W m-2", 0, "");
    ice_output.add("epsilon.mass", nan, "kg m-2 s-1", 0, "");
    ice_output.add("epsilon.enth", nan, "W m-2", 0, "");

    std::cout << "========= Ice Model Outputs (" << ice_coupler->name() << ") modele_pism.cpp:" << std::endl;
    std::cout << ice_output << std::endl;

    // ------- Variable and unit conversions, GCM <- Ice
    {
    VarTransformer &vtA(ice_coupler->var_trans_outAE[GridAE::A]);
    VarTransformer &vtE(ice_coupler->var_trans_outAE[GridAE::E]);

    // NOTE: coupler->gcm_inputs is set up through calls to add_gcm_input_xx() in LANDICE_COM.f
    vtA.set_dims(
        gcm_coupler->gcm_inputsAE[GridAE::A].keys(),    // outputs
        ice_output.keys(),            // inputs
        gcm_coupler->scalars.keys());    // scalars

    vtE.set_dims(
        gcm_coupler->gcm_inputsAE[GridAE::E].keys(),    // outputs
        ice_output.keys(),            // inputs
        gcm_coupler->scalars.keys());    // scalars

    ok = ok && vtA.set("elevA", "ice_surface_elevation", UNIT, 1.0);
    ok = ok && vtE.set("elevE", "ice_surface_elevation", UNIT, 1.0);

    // Top layer state from ice model
    ok = ok && vtE.set("surface_senth", "surface_senth", UNIT, 1.0);
    ok = ok && vtE.set("surface_senth", UNIT, UNIT, -enth_modele_to_pism);

    ok = ok && vtA.set("basal_frictional_heating", "basal_frictional_heating", UNIT, 1.0);
    ok = ok && vtA.set("strain_heating", "strain_heating", UNIT, 1.0);

    ok = ok && vtA.set("geothermal_flux", "geothermal_flux", UNIT, 1.0);
    ok = ok && vtA.set("upward_geothermal_flux", "upward_geothermal_flux", UNIT, 1.0);

    ok = ok && vtA.set("basal_runoff.mass", "melt_grounded.mass", UNIT, 1.0);
    ok = ok && vtA.set("basal_runoff.enth", "melt_grounded.enth", UNIT, 1.0);
    ok = ok && vtA.set("basal_runoff.enth", "melt_grounded.mass", UNIT, -enth_modele_to_pism);

    ok = ok && vtA.set("basal_runoff.mass", "melt_floating.mass", UNIT, 1.0);
    ok = ok && vtA.set("basal_runoff.enth", "melt_floating.enth", UNIT, 1.0);
    ok = ok && vtA.set("basal_runoff.enth", "melt_floating.mass", UNIT, -enth_modele_to_pism);

    ok = ok && vtA.set("calving.mass", "calving.mass", UNIT, 1.0);
    ok = ok && vtA.set("calving.enth", "calving.enth", UNIT, 1.0);
    ok = ok && vtA.set("calving.enth", "calving.mass", UNIT, -enth_modele_to_pism);


    ok = ok && vtA.set("internal_advection.mass", "internal_advection.mass", UNIT, 1.0);
    ok = ok && vtA.set("internal_advection.enth", "internal_advection.enth", UNIT, 1.0);
    ok = ok && vtA.set("internal_advection.enth", "internal_advection.mass", UNIT, -enth_modele_to_pism);


    ok = ok && vtA.set("epsilon.mass", "epsilon.mass", UNIT, 1.0);
    ok = ok && vtA.set("epsilon.enth", "epsilon.enth", UNIT, 1.0);
    ok = ok && vtA.set("epsilon.enth", "epsilon.mass", UNIT, -enth_modele_to_pism);

    }

    // Catch all our errors at once
    if (!ok) (*ibmisc_error)(-1,
        "Error(s) setting up contract modele_pism.");
    printf("END setup_modele_pism()\n");
}


}}
// --------------------------------------------------------


// Contracts should also specify how constants are agreed upon between the two parties.
// 
// PISM needs at least the following constants:
// 
// EnthalpyConverter::EnthalpyConverter(const PISMConfig &config) {
//   beta  = config.get("beta_CC");                                 // K Pa-1
//   c_i   = config.get("ice_specific_heat_capacity");              // J kg-1 K-1
//   g     = config.get("standard_gravity");                        // m s-2
//   L     = config.get("water_latent_heat_fusion");                // J kg-1
//   p_air = config.get("surface_pressure");                        // Pa
//   rho_i = config.get("ice_density");                             // kg m-3
//   T_melting = config.get("water_melting_point_temperature");       // K  
//   T_tol = config.get("cold_mode_is_temperate_ice_tolerance");    // K 
//   T_0   = config.get("enthalpy_converter_reference_temperature");// K  
// 
//   do_cold_ice_methods  = config.get_flag("do_cold_ice_methods");
// }

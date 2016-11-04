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
#include <icebin/contracts/contracts.hpp>
#include <icebin/modele/GCMCoupler_ModelE.hpp>
#include <icebin/pism/IceModel_PISM.hpp>

// --------------------------------------------------------
namespace icebin {
namespace contracts {

static double const nan = std::numeric_limits<double>::quiet_NaN();

/** GCM-specific contract */
void setup_modele_pism(GCMCoupler const *_coupler, IceModel *_model)
{
    // Get arguments we need from coupler
    auto coupler(dynamic_cast<modele::GCMCoupler_ModelE const *>(_coupler));
    auto model(dynamic_cast<icebin::gpism::IceModel_PISM *>(_model));
    

    printf("BEGIN setup_modele_pism()\n");

    // =========== Transfer  constants from Icebin's gcm_constants --> PISM
    // model->transfer_constant(PISM_destionation, icebin_src, multiply_by)
    // (see IceModel_PISM.cpp)
    model->transfer_constant("standard_gravity", "constant::grav");
    model->transfer_constant("beta_CC", "seaice::dtdp", -1.0);
    model->transfer_constant("water_melting_point_temperature", "constant::tf");
    model->transfer_constant("water_latent_heat_fusion", "constant::lhm");
    model->transfer_constant("water_specific_heat_capacity", "constant::shw");

    model->transfer_constant("ice_density", "constant::rhoi");
    model->transfer_constant("ice_thermal_conductivity", "seaice::alami0");
    model->transfer_constant("ice_specific_heat_capacity", "constant::shi");
    model->transfer_constant("fresh_water_density", "constant::rhow");
    model->transfer_constant("sea_water_density", "constant::rhows");
    model->transfer_constant("standard_gravity", "constant::grav");
    model->transfer_constant("ideal_gas_constant", "constant::gasc");

    // To set this, see (in ModelE): Function SHCGS in ocnfuntab.f is
    // used for the Russell ocean. I. The simple models use SHW=4185.
    // This probably doesn't matter much at this point (May 2014)
    //      model->transfer_constant("sea_water_specific_heat_capacity", "");

    /* The following constants were not transferred 
    pism_config:fill_value = -2e9;
    pism_config:fill_value_doc = "_FillValue used when saving diagnostic quantities";
    */

    // In PISM and ModelE Clausius-Clapeyron equation, surfce_pressure is the DIFFERENCE
    // from 1atm.  Thus, surface_pressure=0 implies the ice sheet existing at 1atm
    model->set_constant("surface_pressure", 0, "Pa");       // Match ModelE thermodynam

    // No need to set enthalpy_reference_temperature.  pism::EnthalpyConverter is used (below)
    // to convert enthalpy values between ModelE and PISM.
    // model->transfer_constant("enthalpy_reference_temperature", "enthalpy_reference_temperature");


    // ============ GCM -> Ice
    VarSet &ice_input(contract[IceModel::INPUT]);

    // ------ Decide on the coupling contract for this ice sheet
    ice_input.add_field("massxfer", 0., "kg m-2 s-1", contracts::ELEVATION,
        "Mass of ice being transferred Stieglitz --> Icebin");
    ice_input.add_field("enthxfer", 0., "W m-2", contracts::ELEVATION,
        "Enthalpy of ice being transferred Stieglitz --> Icebin");
    ice_input.add_field("deltah", 0., "J m-2", contracts::ELEVATION,
        "Change of enthalpy of top layer in PISM");

    // Figure out the conversion between GCM and PISM enthalpy
    // ModelE's reference state is 1atm, 0C, 100% liquid water.
    // The enthalpy for that reference state would be the top end
    // of PISM's EnthalpyInterval.
    // NOTE: Pressure in PISM is RELATIVE to atmospheric pressure.
    //       Thus, p=0 is the correct to use at the top surface of
    //       the ice sheet (where ModelE operates).
    pism::EnthalpyConverter enth(*config);
    double const pressure = 0;
    double E_l;
    // getEnthalpyInterval() replaced with enthalpy_liquid()
    // https://github.com/pism/pism/commit/c820dfd8
    E_l = enth.enthalpy_liquid(pressure);
    double const enth_modele_to_pism = E_l;     // (J/kg): Add to convert ModelE specific enthalpies (J/kg) to PISM specific enthalpies (J/kg)
    // NOTE: enth_modele_to_pism == 437000 J/kg
    if (pism_rank == 0) printf("enth_modele_to_pism = %g\n", enth_modele_to_pism);

    bool ok = true;

    double const RHOW = coupler->gcm_constants.get_as("constant::rhow", "kg m-3");
    double const byRHOW = 1.0 / RHOW;

    // ------------- Convert the contract to a var transformer
    {VarTransformer &vt(var_transformer[IceModel::INPUT]);
    vt.set_names(VarTransformer::INPUTS, &coupler->gcm_outputs);
    vt.set_names(VarTransformer::OUTPUTS, &ice_input);
    vt.set_names(VarTransformer::SCALARS, &coupler->ice_input_scalars);
    vt.allocate();

    ok = ok && vt.set("massxfer", "massxfer", "by_dt", 1.0);
    ok = ok && vt.set("enthxfer", "enthxfer", "by_dt", 1.0);
    ok = ok && vt.set("enthxfer", "massxfer", "by_dt", enth_modele_to_pism);
    ok = ok && vt.set("deltah", "deltah", "unit", 1.0);
    }

    // ============== Ice -> GCM
    CouplingContract &ice_output(contract[IceModel::OUTPUT]);

    // Icebin requires that all ice models return elev2, so that it can regrid in the vertical.

    ice_output.add_field("ice_surface_elevation", nan, "m", contracts::ICE|contracts::INITIAL, "ice upper surface elevation");
    ice_output.add_field("ice_thickness", nan, "m", contracts::ICE|contracts::INITIAL, "thickness of ice");
    ice_output.add_field("bed_topography", nan, "m", contracts::ICE|contracts::INITIAL, "topography of bedrock");


    ice_output.add_field("mask", nan, "", contracts::ICE|contracts::INITIAL, "PISM land surface type");

    ice_output.add_field("M1", nan, "kg m-2", contracts::ICE|contracts::INITIAL, "");
    ice_output.add_field("M2", nan, "kg m-2", contracts::ICE|contracts::INITIAL, "");
    ice_output.add_field("H1", nan, "J m-2", contracts::ICE|contracts::INITIAL, "");
    ice_output.add_field("H2", nan, "J m-2", contracts::ICE|contracts::INITIAL, "");
    ice_output.add_field("V1", nan, "m^3 m-2", contracts::ICE|contracts::INITIAL, "");
    ice_output.add_field("V2", nan, "m^3 m-2", contracts::ICE|contracts::INITIAL, "");

    ice_output.add_field("basal_frictional_heating", nan, "W m-2", contracts::ICE, "");
    ice_output.add_field("strain_heating", nan, "W m-2", contracts::ICE, "");
    ice_output.add_field("geothermal_flux", nan, "W m-2", contracts::ICE, "");
    ice_output.add_field("upward_geothermal_flux", nan, "W m-2", contracts::ICE, "");

    ice_output.add_field("calving.mass", nan, "kg m-2 s-1", contracts::ICE, "");
    ice_output.add_field("calving.enth", nan, "W m-2", contracts::ICE, "");
    ice_output.add_field("icebin_smb.mass", nan, "kg m-2 s-1", contracts::ICE, "");
    ice_output.add_field("icebin_smb.enth", nan, "W m-2", contracts::ICE, "");
    ice_output.add_field("pism_smb.mass", nan, "kg m-2 s-1", contracts::ICE, "");
    ice_output.add_field("pism_smb.enth", nan, "W m-2", contracts::ICE, "");

    // basal_runoff (GCM input) = melt_grounded + melt_floatig (PISM outputs)
    ice_output.add_field("melt_grounded.mass", nan, "kg m-2 s-1", contracts::ICE, "");
    ice_output.add_field("melt_grounded.enth", nan, "W m-2", contracts::ICE, "");
    ice_output.add_field("melt_floating.mass", nan, "kg m-2 s-1", contracts::ICE, "");
    ice_output.add_field("melt_floating.enth", nan, "W m-2", contracts::ICE, "");

    ice_output.add_field("internal_advection.mass", nan, "kg m-2 s-1", contracts::ICE, "");
    ice_output.add_field("internal_advection.enth", nan, "W m-2", contracts::ICE, "");
    ice_output.add_field("epsilon.mass", nan, "kg m-2 s-1", contracts::ICE, "");
    ice_output.add_field("epsilon.enth", nan, "W m-2", contracts::ICE, "");

    ice_output.add_field("unit", nan, "", 0, "Dimensionless identity");

    std::cout << "========= Ice Model Outputs (" << sheet->name << ") modele_pism.cpp:" << std::endl;
    std::cout << ice_output << std::endl;

    // ------- Variable and unit conversions, Ice -> GCM
    {VarTransformer &vt(var_transformer[IceModel::OUTPUT]);

    // NOTE: coupler->gcm_inputs is set up through calls to add_gcm_input_xx() in LANDICE_COM.f
    vt.set_dim(VarTransformer::INPUTS, &ice_output);
    vt.set_dim(VarTransformer::OUTPUTS, &coupler->gcm_inputs);
    vt.set_dim(VarTransformer::SCALARS, &coupler->ice_input_scalars);
    vt.allocate();

//  ok = ok && vt.set("elev2", "ice_surface_elevation", "unit", 1.0);
    ok = ok && vt.set("elev1", "ice_surface_elevation", "unit", 1.0);

    // Top layer state from ice model
    ok = ok && vt.set("M1", "M1", "unit", 1.0); // Divide by RHOW to convert to m water equiv
    ok = ok && vt.set("H1", "H1", "unit", 1.0);
    ok = ok && vt.set("H1", "M1", "unit", -enth_modele_to_pism);
    ok = ok && vt.set("V1", "V1", "unit", 1.0);

    // Second-top layer state from ice model
    ok = ok && vt.set("M2", "M2", "unit", 1.0); // Divide by RHOW to convert to m water equiv
    ok = ok && vt.set("H2", "H2", "unit", 1.0);
    ok = ok && vt.set("H2", "M2", "unit", -enth_modele_to_pism);
    ok = ok && vt.set("V2", "V2", "unit", 1.0);


    ok = ok && vt.set("basal_frictional_heating", "basal_frictional_heating", "unit", 1.0);
    ok = ok && vt.set("strain_heating", "strain_heating", "unit", 1.0);

    ok = ok && vt.set("geothermal_flux", "geothermal_flux", "unit", 1.0);
    ok = ok && vt.set("upward_geothermal_flux", "upward_geothermal_flux", "unit", 1.0);

    ok = ok && vt.set("basal_runoff.mass", "melt_grounded.mass", "unit", 1.0);
    ok = ok && vt.set("basal_runoff.enth", "melt_grounded.enth", "unit", 1.0);
    ok = ok && vt.set("basal_runoff.enth", "melt_grounded.mass", "unit", -enth_modele_to_pism);

    ok = ok && vt.set("basal_runoff.mass", "melt_floating.mass", "unit", 1.0);
    ok = ok && vt.set("basal_runoff.enth", "melt_floating.enth", "unit", 1.0);
    ok = ok && vt.set("basal_runoff.enth", "melt_floating.mass", "unit", -enth_modele_to_pism);

    ok = ok && vt.set("calving.mass", "calving.mass", "unit", 1.0);
    ok = ok && vt.set("calving.enth", "calving.enth", "unit", 1.0);
    ok = ok && vt.set("calving.enth", "calving.mass", "unit", -enth_modele_to_pism);


    ok = ok && vt.set("internal_advection.mass", "internal_advection.mass", "unit", 1.0);
    ok = ok && vt.set("internal_advection.enth", "internal_advection.enth", "unit", 1.0);
    ok = ok && vt.set("internal_advection.enth", "internal_advection.mass", "unit", -enth_modele_to_pism);


    ok = ok && vt.set("epsilon.mass", "epsilon.mass", "unit", 1.0);
    ok = ok && vt.set("epsilon.enth", "epsilon.enth", "unit", 1.0);
    ok = ok && vt.set("epsilon.enth", "epsilon.mass", "unit", -enth_modele_to_pism);

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

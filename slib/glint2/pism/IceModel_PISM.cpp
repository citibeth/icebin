#include <mpi.h>	// For Intel MPI, mpi.h must be included before stdio.h
#include <PISMStressBalance.hh>
#include <glint2/pism/IceModel_PISM.hpp>
#include <giss/ncutil.hpp>
#include <boost/filesystem.hpp>
#include <boost/format.hpp>
#include <cmath>
#include <cstdio>
#include <sstream>
#include <string>
#include <glint2/GCMCoupler.hpp>
#include <glint2/contracts/contracts.hpp>

using namespace giss;
using namespace pism;

namespace glint2 {
namespace gpism {


IceModel_PISM::IceModel_PISM(std::string const &_name, GCMCoupler const *_coupler)
	: IceModel(IceModel::Type::PISM, _name, _coupler),
	write_pism_inputs(true)
{
	printf("BEGIN/END IceModel_PISM::IceModel_PISM\n");
}


IceModel_PISM::~IceModel_PISM()
{
	if (deallocate() != 0) {
		PetscPrintf(pism_comm, "IceModel_PISM::IceModel_PISM(...): allocate() failed\n");
		PISMEnd();
	}
}


void IceModel_PISM::update_ice_sheet(
	NcFile &nc,
	std::string const &vname,
	IceSheet *sheet)
{
	printf("BEGIN IceModel_PISM::update_ice_sheet(%s)\n", vname.c_str());

	auto pism_var = nc.get_var((vname + ".pism").c_str());	// PISM parameters
	auto pism_i_att(giss::get_att(pism_var, "i"));	// PISM -i argument (input file)
	std::string pism_i = boost::filesystem::absolute(boost::filesystem::path(
		pism_i_att->as_string(0)), coupler->gcm_params.config_dir).string();

	// Read variables from PISM input file
	// byte mask(time, x, y) ;
	// 		mask:units = "" ;
	// 		mask:coordinates = "lat lon" ;
	// 		mask:flag_meanings = "ice_free_bedrock grounded_ice floating_ice ice_free_ocean" ;
	// 		mask:grid_mapping = "mapping" ;
	// 		mask:long_name = "ice-type (ice-free/grounded/floating/ocean) integer mask" ;
	// 		mask:pism_intent = "diagnostic" ;
	// 		mask:flag_values = 0b, 2b, 3b, 4b ;
	// double thk(time, x, y) ;
	// 		thk:units = "m" ;
	// 		thk:valid_min = 0. ;
	// 		thk:coordinates = "lat lon" ;
	// 		thk:grid_mapping = "mapping" ;
	// 		thk:long_name = "land ice thickness" ;
	// 		thk:pism_intent = "model_state" ;
	// 		thk:standard_name = "land_ice_thickness" ;
	// double topg(time, x, y) ;
	// 		topg:units = "m" ;
	// 		topg:coordinates = "lat lon" ;
	// 		topg:grid_mapping = "mapping" ;
	// 		topg:long_name = "bedrock surface elevation" ;
	// 		topg:pism_intent = "model_state" ;
	// 		topg:standard_name = "bedrock_altitude" ;
printf("Opening PISM file for elev2 and mask2: %s\n", pism_i.c_str());
	NcFile ncin(pism_i.c_str());
	long ntime = ncin.get_dim("time")->size();
	long nx = ncin.get_dim("x")->size();
	long ny = ncin.get_dim("y")->size();

	blitz::Array<ncbyte,2> mask(nx, ny);
	NcVar *mask_var = ncin.get_var("mask");
	mask_var->set_cur(ntime-1, 0, 0);
	mask_var->get(mask.data(), 1, nx, ny);

	NcVar *thk_var = ncin.get_var("thk");
	blitz::Array<double,2> thk(nx, ny);
	thk_var->set_cur(ntime-1, 0, 0);
	thk_var->get(thk.data(), 1, nx, ny);

	NcVar *topg_var = ncin.get_var("topg");
	blitz::Array<double,2> topg(nx, ny);
	topg_var->set_cur(ntime-1, 0, 0);
	topg_var->get(topg.data(), 1, nx, ny);

	ncin.close();

	// Transpose and copy the data
	if (!sheet->mask2.get()) sheet->mask2.reset(
		new blitz::Array<int,1>(glint2_grid->ndata()));
	for (int i=0; i<nx; ++i) {
	for (int j=0; j<ny; ++j) {
		int ix2 = glint2_grid->ij_to_index(i, j);
		// Set update_elevation=false to temporarily "repair" fields that are
		// broken due to a difference in the elevation classes used to generate
		// the SMB in a one-way coupled run, and the elevation classes
		// as defined by the ice sheet.
		if (update_elevation) sheet->elev2(ix2) = topg(i,j) + thk(i,j);

		// Mask uses same convention as MATPLOTLIB: 1 = masked out
		// Even if not update_elevation, it is still important to update the
		// mask, or else PISM will just not work.
		(*sheet->mask2)(ix2) = (mask(i,j) == 2 ? 0 : 1);
	}}

	printf("END IceModel_PISM::update_ice_sheet()\n");
}



// See: http://stackoverflow.com/questions/3418231/replace-part-of-a-string-with-another-string
bool replace(std::string& str, const std::string& from, const std::string& to) {
    size_t start_pos = str.find(from);
    if(start_pos == std::string::npos)
        return false;
    str.replace(start_pos, from.length(), to);
    return true;
}

void IceModel_PISM::transfer_constant(std::string const &dest, std::string const &src, double multiply_by, bool set_new)
{

	// Make sure the PISM constant already exists
	if (!set_new && !config->is_set(dest)) {
		fprintf(stderr, "IceModel_PISM::transfer_constant: Trying to set '%s', which is not a PISM configuration parameter.  Is it misspelled?\n", dest.c_str());
		throw std::exception();
	}

	// Discover the units PISM requires.
	std::string doc = config->get_string(dest + "_doc");
	std::string units = doc.substr(0, doc.find(';'));
	double val = coupler->gcm_constants.get_as(src, units) * multiply_by;
	config->set_double(dest, val);
printf("IceModel_PISM::transfer_constant: %s = %g %s (from %s in GCM)\n", dest.c_str(), val, units.c_str(), src.c_str());
}

void IceModel_PISM::set_constant(std::string const &dest, double src_val, std::string const &src_units, bool set_new)
{
	// Make sure the PISM constant already exists
	if (!set_new && !config->is_set(dest)) {
		fprintf(stderr, "IceModel_PISM::set_constant: Trying to set '%s', which is not a PISM configuration parameter.  Is it misspelled?\n", dest.c_str());
		throw std::exception();
	}

	ConstantSet const &gcm_constants(coupler->gcm_constants);

	// Discover the units PISM requires.
	std::string doc = config->get_string(dest + "_doc");
	std::string dest_units = doc.substr(0, doc.find(';'));

	UTUnit usrc(gcm_constants.ut_system->parse(src_units));
	UTUnit udest(gcm_constants.ut_system->parse(dest_units));
	CVConverter cv(usrc, udest);
	double dest_val = cv.convert(src_val);

	config->set_double(dest, dest_val);
printf("IceModel_PISM::transfer_constant: %s = %g %s (from %s in GCM)\n", dest.c_str(), dest_val, dest_units.c_str(), usrc.c_str());
}



// Arguments that are paths, and thus need pathname resolution
// For stable0.5 branch
// static std::set<std::string> path_args = {"config_override", "i", "o", "surface_given_file", "extra_file", "ts_file"};
// For dev branch
static std::set<std::string> path_args = {"i", "o", "surface_given_file", "ocean_kill", "ocean_kill_file", "extra_file", "ts_file"};

void IceModel_PISM::init(
	std::shared_ptr<glint2::Grid> const &grid2,
	NcFile &nc,
	std::string const &vname_base)
{
	printf("BEGIN IceModel_PISM::init(%s)\n", vname_base.c_str());
	IceModel::init(grid2);
	GCMParams const &_gcm_params(coupler->gcm_params);

	this->glint2_grid = std::dynamic_pointer_cast<Grid_XY const>(grid2);

	// General args passed to the ice sheet, regardless of which ice model is being used
	auto info_var = nc.get_var((vname_base + ".info").c_str());
	// PISM parameters, passed to PISM via argv
	auto pism_var = nc.get_var((vname_base + ".pism").c_str());

	// Get simple arguments
	update_elevation = giss::get_att_as_bool(info_var, "update_elevation");

	// Create arguments from PISM configuration
	pism_args.push_back("glint2_pism");

	// Get arguments from GLINT2 configuration
	for (int i=0; i<pism_var->num_atts(); ++i) {
		auto att = giss::get_att(pism_var, i);
		std::string name = att->name();
		std::string val = att->as_string(0);

		if (path_args.find(name) == path_args.end()) {
			// Regular case, just use the value there
			val = std::string(att->as_string(0));
		} else {
			// Resolve path names according to the configuration directory
			val = boost::filesystem::absolute(
				boost::filesystem::path(att->as_string(0)),
				coupler->gcm_params.config_dir).string();
			printf("IceModel_PISM resolving %s: %s --> %s\n", name.c_str(), att->as_string(0), val.c_str());
		}

		pism_args.push_back("-" + name);
		pism_args.push_back(val);
	}

#if 0
	// Hard-code these variables because we will always need them
	pism_args.push_back("-extra_vars");
	pism_args.push_back("climatic_mass_balance_cumulative,nonneg_flux_cumulative,grounded_basal_flux_cumulative,floating_basal_flux_cumulative,flux_divergence");
#endif
}

void IceModel_PISM::start_time_set()
{
	if (allocate() != 0) {
		PetscPrintf(coupler->gcm_params.gcm_comm, "IceModel_PISM::start_time_set(...): allocate() failed\n");
		PISMEnd();
	}
}

PetscErrorCode IceModel_PISM::allocate()
{
	// Convert PISM arguments to old C style
	int argc = pism_args.size();
	char *argv_array[argc];
	std::vector<char> all_str;
	for (int i=0; i<argc; ++i) {
		std::string &arg = pism_args[i];
		for (unsigned int j=0; j<arg.size(); ++j) all_str.push_back(arg[j]);
		all_str.push_back('\0');
	}
	char *pos = &all_str[0];
	for (int i=0; i<argc; ++i) {
		std::string &arg = pism_args[i];
		argv_array[i] = pos;
		pos +=arg.size() + 1;
	}
	char **argv = argv_array;

printf("*** PISM Args:");
for (int i=0; i<argc; ++i) printf(" %s", argv[i]);
printf("\n");


	// Set up communicator for PISM to use
	// Use same group of processes.
	// No spawning or intercommunicators for now --- maybe not ever.
//	MPI_Comm_dup(gcm_params.gcm_comm, &pism_comm);
	pism_comm = coupler->gcm_params.gcm_comm;
	PetscErrorCode ierr;
	ierr = MPI_Comm_rank(pism_comm, &pism_rank); CHKERRQ(ierr);
	ierr = MPI_Comm_size(pism_comm, &pism_size); CHKERRQ(ierr);

printf("[%d] pism_size = %d\n", pism_rank, pism_size);

	// Initialize Petsc
printf("Initializing PETSc\n");
	petsc_context.reset(new PetscContext(pism_comm, argc, argv));

    unit_system.reset(new pism::UnitSystem(NULL));
    config.reset(new pism::Config(pism_comm, "pism_config", *unit_system));
	overrides.reset(new pism::Config(pism_comm, "pism_overrides", *unit_system));
	ierr = init_config(pism_comm, *config, *overrides, true); CHKERRQ(ierr);

	// Get arguments from the GCM
	giss::time::tm const &tb(coupler->gcm_params.time_base);
	std::string reference_date = (boost::format("%04d-%02d-%02d") % tb.year() % tb.month() % tb.mday()).str();
	config->set_string("reference_date", reference_date);

#if 0
	// We don't need to write this file, it is essentially the
	// same as one PISM already writes.
	if (pism_rank == 0) {
		// PISM isn't able to write a new file, so we must create it first.
		NcFile nc(full_fname.c_str(), NcFile::Replace);
		nc.close();
	}
	auto full_fname(coupler->gcm_params.config_dir / "glint2_pism_config.nc");
	printf("IceModel_PISM writing config (1) to: %s\n", full_fname.c_str());
	config->write(full_fname.c_str());
#endif

    pism_grid.reset(new ::IceGrid(pism_comm, *config));
printf("pism_grid=%p: (xs,xm,ys,ym,Mx,My) = %d %d %d %d %d %d %ld %ld\n", &*pism_grid, pism_grid->xs, pism_grid->xm, pism_grid->ys, pism_grid->ym, pism_grid->Mx, pism_grid->My, pism_grid->x.size(), pism_grid->y.size());
	PISMIceModel::Params params;
		params.time_start_s = coupler->gcm_params.time_start_s;
		params.output_dir = boost::filesystem::absolute(coupler->gcm_params.config_dir).string();
    ice_model.reset(new PISMIceModel(*pism_grid, *config, *overrides, params));

	// Transfer constants from GCM to PISM, and also set up coupling contracts.
	// This is the right place to do it, since the PISM systme is fully up and functional,
	// and all PISM config files have been read.
	// This call through the GCMCoupler will call back to setup_contracts_xxx().
	coupler->setup_contracts(*this);

printf("[%d] start = %f\n", pism_rank, pism_grid->time->start());
printf("[%d] end = %f\n", pism_rank, pism_grid->time->end());

    ierr = ice_model->setExecName("pismr"); CHKERRQ(ierr);
	// This has the following stack trace:
	// 	IceModel::init()					[iceModel.cc]
	// 	IceModel::model_state_setup()		[iMinit.cc]
	// 	IceModel::init_couplers()			[iMinit.cc]
	// 	surface->init()
    ierr = ice_model->init(); CHKERRQ(ierr);

	// ============== Set up variables for INPUT contract

	// During the ice_model->init() call above the PISMIceModel
	// class (derived from PISM's IceModel) allocated an instance
	// of PSConstantGLINT2. This instance is owned and will be
	// de-allocated by PISMIceModel ice_model.
	pism_surface_model = ice_model->ps_constant_glint2();

	// Set up corresponence between GLINT2 fields and variables
	// in the PISM data structures.
	int ix;
	pism_ivars.resize(contract[INPUT].size_nounit(), NULL);
	ix = contract[INPUT].index("smb_mass");
		pism_ivars[ix] = &pism_surface_model->glint2_smb_mass;

	// We don't really use this, but we do need to store and pass through for conservation computations
	ix = contract[INPUT].index("smb_enth");
		pism_ivars[ix] = &pism_surface_model->glint2_smb_enth;

	// Ignore surface_temp, it is not useful...
	ix = contract[INPUT].index("surface_temp");
		pism_ivars[ix] = &pism_surface_model->glint2_surface_temp;

	ix = contract[INPUT].index("heat_flux");	// Positive is down
		pism_ivars[ix] = &pism_surface_model->glint2_heat_flux;

	// Initialize scatter/gather stuff
	ierr = pism_grid->get_dm(1, pism_grid->max_stencil_width, da2); CHKERRQ(ierr);

	ierr = DMCreateGlobalVector(da2, &g2); CHKERRQ(ierr);


	// note we want a global Vec but reordered in the natural ordering
	// so when it is scattered to proc zero it is not all messed up;
	// see above
	ierr = DMDACreateNaturalVector(da2, &g2natural); CHKERRQ(ierr);

	// next get context *and* allocate samplep0 (on proc zero only, naturally)
	ierr = VecScatterCreateToZero(g2natural, &scatter, &Hp0); CHKERRQ(ierr);


	// ============== Set up variables for OUTPUT contract
	// -------------- Allocate blitz:Array<double,1> output variables,
	// which are passed back to Glint2.

	// -------------- Link to PISM-format output variables, used to fill ovars
	pism_ovars.resize(contract[OUTPUT].size_nounit(), NULL);
	ix = contract[OUTPUT].index("usurf");		// Elevation of top surface of ice sheet
		pism_ovars[ix] = &ice_model->ice_surface_elevation;	// see PISM's iceModel.hh
	ix = contract[OUTPUT].index("ice_surface_enth");		// Specific enthalpy of top surface
		pism_ovars[ix] = &ice_model->ice_surface_enth;
	ix = contract[OUTPUT].index("ice_surface_enth_depth");
		pism_ovars[ix] = &ice_model->ice_surface_enth_depth;

	ix = contract[OUTPUT].index("effective_surface_temp");
		pism_ovars[ix] = &ice_model->ps_constant_glint2()->glint2_surface_temp;

	// For MassEnergyBudget variables that have a contract name specified,
	// link them up into pism_ovars now.
	for (auto ii = ice_model->rate.all_vecs.begin(); ii != ice_model->rate.all_vecs.end(); ++ii) {
		if (ii->contract_name == "") continue;

		int ix = contract[OUTPUT].index(ii->contract_name);
		pism_ovars[ix] = &ii->vec;
	}

	// ============== Miscellaneous
	// Check that grid dimensions match
	if ((pism_grid->Mx != glint2_grid->nx()) || (pism_grid->My != glint2_grid->ny())) {
		fprintf(stderr, "Grid mismatch: pism=(%d, %d) glint2=(%d, %d)\n", pism_grid->Mx, pism_grid->My, glint2_grid->nx(), glint2_grid->ny());
		throw std::exception();
	}

	printf("END IceModel_PISM::allocate()\n");
	return 0;
}

PetscErrorCode IceModel_PISM::deallocate()
{
	PetscErrorCode ierr;

	ierr = VecDestroy(&g2); CHKERRQ(ierr);
	ierr = VecDestroy(&g2natural); CHKERRQ(ierr);
	// ierr = VecScatterDestroy(&scatter); CHKERRQ(ierr);
//	ierr = VecDestroy(&Hp0); CHKERRQ(ierr);

	return 0;
}


// --------------------------------------------------------

// --------------------------------------------------
/** glint2_var Variable, already allocated, to receive data
@param glint2_var_xy The array to write into (on the root node).
If this array is not yet allocated (ROOT NODE ONLY), it will be allocated.*/
PetscErrorCode IceModel_PISM::iceModelVec2S_to_blitz_xy(IceModelVec2S &pism_var, blitz::Array<double,2> &ret)
{
	PetscErrorCode ierr;

	if (am_i_root()) {
		auto xy_shape(blitz::shape(ny(), nx()));
		if (ret.size() == 0) {
			ret.reference(blitz::Array<double,2>(ny(), nx()));
		} else {
			if (ret.extent(0) != xy_shape[0] || ret.extent(1) != xy_shape[1]) {
				fprintf(stderr, "IceModel_PISM::iceModelVec2S_to_blitz_xy(): ret(%d, %d) should be (%d, %d)\n", ret.extent(0), ret.extent(1), xy_shape[0], xy_shape[1]);
				throw std::exception();
			}
		}
	}

	if (pism_var.get_dof() != 1)
		SETERRQ(pism_grid->com, 1, "This method only supports IceModelVecs with dof == 1");

	// Gather data to one processor
	PetscScalar **bHp0;
	ierr = pism_var.put_on_proc0(Hp0, scatter, g2, g2natural); CHKERRQ(ierr);

	if (am_i_root()) {
		// Copy it to blitz array (on the root node only)
		ierr = VecGetArray2d(Hp0, pism_grid->Mx, pism_grid->My, 0, 0, &bHp0);
		for (PetscInt i=0; i < pism_grid->Mx; i++) {
			for (PetscInt j=0; j < pism_grid->My; j++) {
				ret(j, i) = bHp0[i][j];
			}
		}
		ierr = VecRestoreArray2d(Hp0, pism_grid->Mx, pism_grid->My, 0, 0, &bHp0); CHKERRQ(ierr);
	}


	return 0;
}
// --------------------------------------------------
void IceModel_PISM::run_timestep(double time_s,
	blitz::Array<int,1> const &indices,
	std::vector<blitz::Array<double,1>> const &ivals2)
{
	printf("BEGIN IceModel_PISM::run_timestep(%f)\n", time_s);

	if (run_timestep_petsc(time_s, indices, ivals2) != 0) {
		PetscPrintf(pism_comm, "IceModel_PISM::runtimestep() failed\n");
		PISMEnd();
	}
	printf("END IceModel_PISM::run_timestep(%f)\n", time_s);
}

// --------------------------------------------------
//static double const nan = std::numeric_limits<double>::quiet_NaN();

PetscErrorCode IceModel_PISM::run_timestep_petsc(double time_s,
	blitz::Array<int,1> const &indices,
	std::vector<blitz::Array<double,1>> const &ivals2)
{
printf("[%d] BEGIN IceModel_PISM::run_timestep_petsc(%f)\n", pism_rank, time_s);
	PetscErrorCode ierr = 0;

	int nvals = indices.extent(0);

	// Check number of variables matches for input
	if (ivals2.size() != pism_ivars.size()) {
		fprintf(stderr, "[%d] IceModel_PISM::run_timestep_petsc: ivals2.size()=%ld does not match pism_ivars.size()=%ld\n", pism_rank, ivals2.size(), pism_ivars.size());
		ierr = -1;
	}
	CHKERRQ(ierr);

	// Check Petsc types
	if (sizeof(double) != sizeof(PetscScalar)) {
		fprintf(stderr, "PetscScalar must be same as double\n");
		ierr = -1;
	}
	CHKERRQ(ierr);

printf("AA1\n");
#if 0
	// Set to NaN, so we can tell where we have no data
	PSConstantGLINT2 *surface = ice_model->ps_constant_glint2();
	ierr = surface->glint2_surface_temp.set(nan); CHKERRQ(ierr);
#endif

	// Needed because indices is not contiguous in RAM (has non-unit stride).
	// Also because Glint2 and PISM index cells in different order
	std::unique_ptr<int[]> g2_ix(new int[nvals]);
	for (int ix=0; ix<nvals; ++ix) g2_ix[ix] = glint2_to_pism1d(indices(ix));
printf("AA1\n");

	std::unique_ptr<PetscScalar[]> g2_y(new PetscScalar[nvals]);
for (int ix=0; ix<nvals; ++ix) g2_y[ix] = 17.;
printf("AA1\n");
	for (unsigned int i=0; i<pism_ivars.size(); ++i) {

printf("BB1\n");
		giss::CoupledField const &cf(contract[IceModel::INPUT].field(i));

		// Get matching input (val) and output (pism_var) variables
		blitz::Array<double,1> const &val(ivals2[i]);
		IceModelVec2S *pism_var = pism_ivars[i];

printf("BB1\n");
		// Inputs specified in the contract are not (necessarily) attached
		// to any PISM var.  If they are not, just drop them on the ground.
		if (!pism_var) continue;

		// Copy sparse array
		for (int ix=0; ix<nvals; ++ix) g2_y[ix] = val(ix);
for (int ix=0; ix<nvals; ++ix) g2_y[ix] += 1.;

		// Put into a natural-ordering global distributed Petsc Vec
		ierr = VecSet(g2natural, cf.default_value); CHKERRQ(ierr);
		ierr = VecSetValues(g2natural, nvals, g2_ix.get(), g2_y.get(), ADD_VALUES); CHKERRQ(ierr);

		ierr = VecAssemblyBegin(g2natural); CHKERRQ(ierr);
		ierr = VecAssemblyEnd(g2natural); CHKERRQ(ierr);

		// Copy to Petsc-ordered global vec
		ierr = DMDANaturalToGlobalBegin(da2, g2natural, INSERT_VALUES, g2); CHKERRQ(ierr);
		ierr = DMDANaturalToGlobalEnd(da2, g2natural, INSERT_VALUES, g2); CHKERRQ(ierr);

		// Copy to the output variable
		// (Could we just do DMDANaturalToGlobal() directly to this?)
		ierr = pism_var->copy_from_vec(g2); CHKERRQ(ierr);

	}

	// At this point, glint2_surface_temp will be NaN in places off the ice sheet.
	// So replace those NaN values with temperatures from the Enth3 array.
	PSConstantGLINT2 * const surface = ice_model->ps_constant_glint2();
	ice_model->merge_surface_temp(surface->glint2_surface_temp,
		contract[INPUT].field("surface_temp").default_value);

	// Write out what wev'e now calculated (for debugging)
	if (write_pism_inputs) {
	for (unsigned int i=0; i<pism_ivars.size(); ++i) {

		// ================ BEGIN Write PISM Inputs
		long time_day = (int)(time_s / 86400. + .5);
		std::stringstream fname;
		std::string const &fnpart = contract[INPUT].name(i);

		fname << time_day << "-" << fnpart << ".nc";
		boost::filesystem::path pfname(coupler->gcm_params.config_dir / "pism_inputs" / fname.str());

//		printf("ICeModel_PISM writing (2) to: %s\n", pfname.c_str());
		IceModelVec2S *pism_var = pism_ivars[i];
		pism_var->dump(pfname.c_str());
		// ================ END Write PISM Inputs
	}}	// For each pism_var



printf("[%d] BEGIN ice_model->run_to(%f -> %f) %p\n", pism_rank, pism_grid->time->current(), time_s, ice_model.get());
	// =========== Run PISM for one coupling timestep
	// Time of last time we coupled
	auto old_pism_time(pism_grid->time->current());
	ierr = ice_model->run_to(time_s); CHKERRQ(ierr);	// See glint2::gpism::PISMIceModel::run_to()
	if ((ice_model->mass_t() != time_s) || (ice_model->enthalpy_t() != time_s)) {
		fprintf(stderr, "ERROR: PISM time (mass=%f, enthalpy=%f) doesn't match GLINT2 time %f\n", ice_model->mass_t(), ice_model->enthalpy_t(), time_s);
		throw std::exception();
	}

	// ice_model->enthalpy_t() == time_s here
	ierr = ice_model->prepare_outputs(old_pism_time); CHKERRQ(ierr);
	ierr = get_state_petsc(0); CHKERRQ(ierr);	// Copy PISM->Glint2 output vars
	ierr = ice_model->reset_rate(); CHKERRQ(ierr);

printf("Current time is pism: %f-%f, GLINT2: %f\n", old_pism_time, pism_grid->time->current(), time_s);
printf("[%d] END ice_model->run_to()\n", pism_rank);

	return 0;
}


/** Copies PISM->Glint2 output variables from PISM variables to
the Glint2-supplied variables (on the root node).
@param mask Only do it for variables where (flags & mask) == mask.  Set to 0 for "all." */
PetscErrorCode IceModel_PISM::get_state_petsc(unsigned int mask)
{
	giss::CouplingContract const &ocontract(contract[IceModel::OUTPUT]);

	// Copy the outputs to the blitz arrays
	if (am_i_root()) allocate_ice_ovals_I();		// ROOT in PISM communicator
	for (unsigned int i=0; i<pism_ovars.size(); ++i) {
		if (!pism_ovars[i]) {
			fprintf(stderr, "IceModel_PISM: Contract output %s is not linked up to a pism_ovar\n", ocontract.field(i).name.c_str());
			throw std::exception();
		}


		giss::CoupledField const &cf(ocontract.field(i));
		if ((cf.flags & mask) != mask) continue;

		printf("[%d] IceModel_PISM::get_state_petsc(mask=%d) copying field %s\n", pism_rank, mask, cf.name.c_str());

		if (am_i_root()) {		// ROOT in PISM communicator

			// Check number of variables matches for output
			if (ice_ovals_I.size() != pism_ovars.size()) {
				fprintf(stderr, "[%d] IceModel_PISM::run_timestep_petsc: ice_ovals_I.size()=%ld does not match pism_ovars.size()=%ld\n", pism_rank, ice_ovals_I.size(), pism_ovars.size());
				throw std::exception();
			}

			// Reshape 1D blitz variable to 2D for use with PISM
			blitz::Array<double,2> oval2_xy(
				giss::reshape<double,1,2>(ice_ovals_I[i], blitz::shape(ny(), nx())));
		
	
			// Get matching input (val) and output (pism_var) variables
			iceModelVec2S_to_blitz_xy(*pism_ovars[i], oval2_xy);	// Allocates oval2_xy if needed

		} else {
			blitz::Array<double,2> oval2_xy;	// dummy
			iceModelVec2S_to_blitz_xy(*pism_ovars[i], oval2_xy);
		}

		// Now send those data from the PISM root to the GCM root
		// (DUMMY for now, just make sure PISM and GCM have the same root)
		if (pism_root != coupler->gcm_params.gcm_root) {
			fprintf(stderr, "PISM and the GCM must share the same root!\n");
			throw std::exception();
		}
	}

	return 0;
}


PetscErrorCode IceModel_PISM::get_initial_state_petsc()
{
	PetscErrorCode ierr;

	// Only prepare PISMIceModel outputs for things we need at init time.
	ierr = ice_model->prepare_initial_outputs(); CHKERRQ(ierr);

	// Copy outputs to Glint2-supplied variables, only for INITIAL variables
printf("[%d] Calling get_state_petsc(%d)\n", pism_rank, contracts::INITIAL);
	ierr = get_state_petsc(contracts::INITIAL); CHKERRQ(ierr);


	return 0;
}

void IceModel_PISM::get_initial_state()
{
	printf("[%d] BEGIN IceModel_PISM::get_initial_state()\n", pism_rank);

	if (get_initial_state_petsc() != 0) {
		PetscPrintf(pism_comm, "IceModel_PISM::get_initial_state() failed\n");
		PISMEnd();
	}
	printf("[%d] END IceModel_PISM::get_initial_state()\n", pism_rank);
}











}}

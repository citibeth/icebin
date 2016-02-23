/*
 * GLINT2: A Coupling Library for Ice Models and GCMs
 * Copyright (c) 2013 by Robert Fischer
 * 
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include <mpi.h>		// Intel MPI wants to be first
#include <glint2/GCMCoupler.hpp>
#include <glint2/MatrixMaker.hpp>
#include <glint2/MultiMatrix.hpp>
#include <glint2/contracts/contracts.hpp>
#include <giss/exit.hpp>

#ifdef USE_PISM
#include <icebin/pism/IceModel_PISM.hpp>
#endif

namespace icebin {

std::unique_ptr<IceModel> new_ice_model(IceModel::Type type,
	GCMCoupler const *_coupler, IceRegridder const *_sheet)
{
	std::unique_ptr<IceModel> ice_model;

	switch(type.index()) {
		case IceModel::Type::DISMAL :
			ice_model.reset(new IceModel_DISMAL);
		break;
		case IceModel::Type::WRITER :
			ice_model.reset(new IceModel_Writer);
		break;
#ifdef USE_PISM
		case IceModel::Type::PSIM :
			ice_model.reset(new IceModel_PSIM);
		break;
#endif
		default :
			(*icebin_error)(-1,
				"Unknown IceModel::Type %s", type.str());
	}


	// Do basic initialization...
	ice_model->coupler = _coupler;
	ice_model->sheet = _sheet;
	ice_model->ice_constants.init(_coupler->ut_system);

	return ice_model;


}

std::unique_ptr<IceModel> new_ice_model(NcIO &ncio, std::string vname,
	GCMCoupler const *_coupler, IceRegridder const *_sheet)
{
	std::string vn(vname + ".info");
	auto info_v = get_or_add_var(ncio, vn, netCDF::ncInt64, {});

	IceModel::Type type;
	get_or_put_att_enum(info_v, ncio.rw, "ice_model", type);

	return new_ice_model(type, _coupler, _sheet);
}

IceModel::IceModel(IceModel::Type _type) : type(_type) {}

IceModel::ncread(
ibmisc::NcIO &ncio, std::string const &vname_sheet)
{
	gcm_per_ice_sheet_params = coupler->ncread_gcm_per_ice_sheet_params(ncio, vname_sheet);
}



IceModel::~IceModel() {}

giss::CouplingContract *IceModel::new_CouplingContract() {
	_extra_contracts.push_back(
		std::unique_ptr<giss::CouplingContract>(
		new giss::CouplingContract()));
	return _extra_contracts[_extra_contracts.size()-1].get();
}

// ==========================================================
bool IceModel::am_i_root() const
	{ return coupler->am_i_root(); }

/** Allocate vectors in preparation of calling an ice model. */
void IceModel::allocate_ice_ovals_I()
{
	// Check for program errors
	if (!coupler->am_i_root()) {
		fprintf(stderr, "IceModel::allocate_ice_ovals_I() should only be called from GCM root MPI node.  Fix the code.\n");
		giss::exit(1);
	}
	if (ice_ovals_I.size() != 0) {
		fprintf(stderr, "[%d] IceModel::allocate_ice_ovals_I(): called twice without a free() inbetween.  Fix the code. (old size is %ld)\n", coupler->gcm_params.gcm_rank, ice_ovals_I.size());
		giss::exit(1);
	}

	// Allocate for direct output from ice model
	giss::CouplingContract const &ocontract(contract[IceModel::OUTPUT]);
	int nfields = ocontract.size_nounit();
	for (int i=0; i < nfields; ++i) {
		giss::CoupledField const &cf(ocontract.field(i));
		long n2 = ndata();
		ice_ovals_I.push_back(blitz::Array<double,1>(n2));
	}
}


/** Allocate in preparation of var transformations (but not regridding yet) */
void IceModel::allocate_gcm_ivals_I()
{
	// Check for program errors
	if (!coupler->am_i_root()) {
		fprintf(stderr, "IceModel::allocate_ice_ivals_I() should only be called from GCM root MPI node.  Fix the code.\n");
		giss::exit(1);
	}
	if (gcm_ivals_I.size() != 0) {
		fprintf(stderr, "IceModel::allocate_gcm_ivals_I(): called twice without a free() inbetween.  Fix the code.\n");
		giss::exit(1);
	}

	giss::CouplingContract const &gcm_inputs(coupler->gcm_inputs);
	int nfields = gcm_inputs.size_nounit();
	for (int i=0; i < nfields; ++i) {
		giss::CoupledField const &cf(gcm_inputs.field(i));
		long n2 = ndata();
		gcm_ivals_I.push_back(blitz::Array<double,1>(n2));
	}
}

/** Free portions not needed after finished calling ice model and
applying variable transform.  This will be variables desired on
anything other than the ELEVATION grid. */
void IceModel::free_ice_ovals_I()
{
	// Check for program errors
	if (!coupler->am_i_root()) {
		fprintf(stderr, "IceModel::free_ice_ovals_I() should only be called from GCM root MPI node.  Fix the code.\n");
		giss::exit(1);
	}

	ice_ovals_I.clear();
}

/** Free all memory used by this.  Called when we're done with a coupling timestep. */
void IceModel::free_ovals_ivals_I()
{
	// Check for program errors
	if (!coupler->am_i_root()) {
		fprintf(stderr, "IceModel::free_ovals_ovals_I() should only be called from GCM root MPI node.  Fix the code.\n");
		giss::exit(1);
	}

	ice_ovals_I.clear();
	gcm_ivals_I.clear();
}

// -----------------------------------------------------------
/** Allocates and sets gcm_ivals_I variable
@param mask Control which GCM input variables to set (according to, etc, INITIAL flag in contract) */
void IceModel::set_gcm_inputs(unsigned int mask)
{
  printf("BEGIN IceModel::set_gcm_inputs()\n");
	allocate_gcm_ivals_I();

	// Compute the variable transformation
	giss::VarTransformer &vt(var_transformer[IceModel::OUTPUT]);
	giss::CSRAndUnits trans = vt.apply_scalars({
		std::make_pair("unit", 1.0)});

	giss::CouplingContract const &gcm_inputs(coupler->gcm_inputs);

	// Apply the variable transformation
	for (int xi=0; xi<vt.dimension(giss::VarTransformer::OUTPUTS).size_nounit(); ++xi) {	// xi is index of output variable
		gcm_ivals_I[xi] = 0;	// Vector operation: clear before sum
		giss::CoupledField const &cf(gcm_inputs.field(xi));

		if ((cf.flags & mask) != mask) continue;

		// Consider each output variable separately...
		std::vector<std::pair<int, double>> const &row(trans.mat[xi]);
		for (auto xjj=row.begin(); xjj != row.end(); ++xjj) {
			int xj = xjj->first;		// Index of input variable
			double io_val = xjj->second;	// Amount to multiply it by
			gcm_ivals_I[xi] += ice_ovals_I[xj] * io_val;		// blitz++ vector operation
		}
		gcm_ivals_I[xi] += trans.units[xi];
	}

	printf("END IceModel::set_gcm_inputs()\n");
}

// -----------------------------------------------------------









// ==========================================================
/** @param nc The Glint22 configuration file */
void GCMCoupler::ncread(
	std::string const &_fname,
	std::string const &_vname,
	ibmisc::Domain<int> &&_domainA)
{
	printf("BEGIN GCMCoupler::read_from_netcdf() %s\n", vname.c_str());

	fname = _fname;
	vname = _vname;
	domainA = std::move(_domainA);

	NcIO ncio(fname, NcFile::ReadOnly);

	// Load the MatrixMaker	(filtering by our domain, of course)
	// Also load the ice sheets
	regridder.ncio(ncio, vname);
	regridder.filter_cellsA(domainA);

	// The root node needs to make the full regridding matrices, not just those
	// for its own domain.  Therefore, create a second MatrixMaker for it.
	if (am_i_root()) {
		regridder_full.reset(new GCMRegridder);
		regridder_full->ncio(ncio, vname);
	}

	// Read gcm_out_file, an optional variable telling the GCM-specific
	// part of IceBin to write out exactly what it sees coming from the GCM
	// (so it can be replayed later with desm)
	auto info_v = get_or_add_var(ncio, vname + ".info", netCDF::ncInt64, {});

	gcm_out_file = "";
	get_or_put_att(info_v, ncio.rw, "gcm_out_file", gcm_out_file, false);
	if (gcm_out_file.length() > 0) {
		gcm_out_file = boost::filesystem::absolute(
			boost::filesystem::path(gcm_out_file),
			gcm_params.run_dir).string();
	}

	gcm_in_file = "";
	get_or_put_att(info_v, ncio.rw, "gcm_in_file", gcm_in_file, false);
	if (gcm_in_file.length() > 0) {
		gcm_in_file = boost::filesystem::absolute(
			boost::filesystem::path(gcm_in_file),
			gcm_params.run_dir).string();
	}



#if 1
	std::cout << "========= GCM Constants" << std::endl;
	std::cout << gcm_constants;
	std::cout << "========= GCM Outputs" << std::endl;
	std::cout << gcm_outputs;
	std::cin << "========= GCM Inputs" << std::endl;
	std::cout << gcm_inputs;


#endif

	ice_models.clear();
	for (size_t i=0; i < regridder.sheets.size(); ++i) {
		IceRegridder const *sheet = &*regridder.sheets[i];
		std::string vname_sheet = vname + "." + sheet->name;

		// Create an IceModel corresponding to this IceSheet.
		std::unique_ptr<IceModel> ice_model(new_ice_model(ncio, vname_sheet, this, sheet));
		ice_model->ncread(ncio, vname_sheet);

		setup_contracts(*this, *ice_model);	// where does this go w.r.t ncread() and upate?
		ice_model->update_ice_sheet(ncio, vname_sheet);

		// ---- Create the affiliated writers
		// Writer for ice model input
		auto iwriter(new_ice_model(IceModel::Type::WRITER, this, sheet));
		iwriter->init(IceModel::INPUT, ice_model);		// Writer-specific initialization
		ice_model->_iwriter = std::move(iwriter);

		// Writer for ice model output
		auto owriter(new_ice_model(IceModel::Type::WRITER, this, sheet));
		owriter->init(IceModel::OUTPUT, ice_model);		// Writer-specific initialization
		ice_model->_owriter = std::move(owriter);
	}

	nc.close();
	printf("END GCMCoupler::read_from_netcdf()\n");
}



void GCMCoupler::set_start_time(
	giss::time::tm const &time_base,
	double time_start_s)
{
	gcm_params.set_start_time(time_base, time_start_s);

	for (int sheetix=0, auto ice_model = ice_models.begin();
		ice_model != ice_models.end(); ++sheetix, ++ice_model)
	{

		model->start_time_set();

		// Handle writing stuff, if the
		IceModel_Writer *iwriter = model->iwriter();	// The affiliated input-writer (if it exists).
		if (iwriter) iwriter->start_time_set();

		IceModel_Writer *owriter = model->owriter();
		if (owriter) owriter->start_time_set();

		// This function is the last phase of initialization.
		// Only now can we assume that the contracts are fully set up.
#if 1
		// Print out the contract and var transformations
		std::cout << "========= Contract for " << model->name << std::endl;
		std::cout << "---- GCM->Ice     Output Variables:" << std::endl;
		std::cout << model->contract[IceModel::INPUT];
		std::cout << "TRANSFORMATIONS:" << std::endl;
		std::cout << model->var_transformer[IceModel::INPUT];
		std::cout << "---- Ice->GCM     Output Variables:" << std::endl;
		std::cout << model->contract[IceModel::OUTPUT];
		std::cout << "TRANSFORMATIONS:" << std::endl;
		std::cout << model->var_transformer[IceModel::OUTPUT];
#endif
	}

}



// ===================================================
// SMBMsg

MPI_Datatype SMBMsg::new_MPI_struct(int nfields)
{
	int nele = 2 + nfields;
	int blocklengths[] = {1, 1, nfields};
	MPI_Aint displacements[] = {offsetof(SMBMsg,sheetno), offsetof(SMBMsg,i2), offsetof(SMBMsg, vals)};
	MPI_Datatype types[] = {MPI_INT, MPI_INT, MPI_DOUBLE};
	MPI_Datatype ret;
	MPI_Type_create_struct(3, blocklengths, displacements, types, &ret);
	MPI_Type_commit(&ret);
	return ret;
}

/** for use with qsort */
int SMBMsg::compar(void const *a, void const *b)
{
	SMBMsg const *aa = reinterpret_cast<SMBMsg const *>(a);
	SMBMsg const *bb = reinterpret_cast<SMBMsg const *>(b);
	int cmp = aa->sheetno - bb->sheetno;
	if (cmp != 0) return cmp;
	return aa->i2 - bb->i2;
}

// ===================================================
// GCMCoupler

/** Parameters to the call_ice_model() method.  Calls to the
ice model are planned first, then executed separately. */
class CallIceModelParams {
public:
	int sheetno;
	SMBMsg *begin;
	SMBMsg *next;

	CallIceModelParams() {}

	CallIceModelParams(int _sheetno, SMBMsg *_begin, SMBMsg *_next) :
		sheetno(_sheetno), begin(_begin), next(_next) {}
};

/** PROTECTED method.  Called by all MPI nodes. */
void GCMCoupler::call_ice_model(
	IceModel *model,
	int sheetno,
	double time_s,
	giss::DynArray<SMBMsg> &rbuf,
	SMBMsg *begin, SMBMsg *end)		// Messages have the MAX number of fields for any ice model contract
{
	// ----------------- Construct input arrays
	// The number of fields for THIS ice sheet will be <= the number
	// of fields in SMBMsg
	int nfields = model->contract[IceModel::INPUT].size();

printf("BEGIN GCMCoupler::call_ice_model(%s, nfields=%ld)\n", model->name.c_str(), nfields);

	// Construct indices vector
	blitz::TinyVector<int,1> shape(rbuf.diff(end, begin));
	blitz::TinyVector<int,1> stride(rbuf.ele_size / sizeof(int));
	blitz::Array<int,1> indices(&begin->i2,
		shape, stride, blitz::neverDeleteData);

	// Construct input values vectors: sparse vectors pointing into the SMBMsgs
	std::vector<blitz::Array<double,1>> ivals2;
	stride[0] = rbuf.ele_size / sizeof(double);
	for (int i=0; i<nfields; ++i) {
		SMBMsg &rbegin(*begin);		// Reference to begin

		ivals2.push_back(blitz::Array<double,1>(
			&rbegin[i], shape, stride, blitz::neverDeleteData));
	}


	// -------------- Run the model
	// Record exactly the same inputs that this ice model is seeing.
	IceModel_Writer *iwriter = model->iwriter();
	if (iwriter) iwriter->run_timestep(time_s, indices, ivals2);

	// Now call to the ice model
//	printf("[%d] BEGIN GCMCoupler::call_ice_model() calling model->run_timestep()\n", rank);
	model->run_timestep(time_s, indices, ivals2);
//	printf("[%d]END GCMCoupler::call_ice_model() calling model->run_timestep()\n", rank);

	// Record what the ice model produced.
	// NOTE: This shouldn't change model->ice_ovals_I
	IceModel_Writer *owriter = model->owriter();
	if (owriter) {
		printf("BEGIN owriter->run_timestep()\n");
		owriter->run_decoded(time_s, model->ice_ovals_I);
		printf("END owriter->run_timestep()\n");
	}

printf("END GCMCoupler::call_ice_model(nfields=%ld)\n", nfields);
};

/** @param sbuf the (filled) array of ice grid values for this MPI node. */
void GCMCoupler::couple_to_ice(
double time_s,
int nfields,			// Number of fields in sbuf.  Not all will necessarily be filled, in the case of heterogeneous ice models.
giss::DynArray<SMBMsg> &sbuf,	// Values, already converted to ice model inputs (from gcm outputs)
std::vector<giss::VectorSparseVector<int,double>> &gcm_ivals)	// Root node only: Already-allocated space to put output values.  Members as defined by the CouplingContract GCMCoupler::gcm_inputs
{
	// TODO: Convert this to use giss::gather_msg_array() instead!!!

	// Gather buffers on root node
	int num_mpi_nodes;
	MPI_Comm_size(gcm_params.gcm_comm, &num_mpi_nodes); 

printf("BEGIN GCMCoupler::couple_to_ice() time_s=%f, sbuf.size=%d, sbuf.ele_size=%d\n", time_s, sbuf.size, sbuf.ele_size);

	// MPI_Gather the count
	std::unique_ptr<int[]> rcounts;
	rcounts.reset(new int[num_mpi_nodes]);
	for (int i=0; i<num_mpi_nodes; ++i) rcounts[i] = 0;

	int nele_l = sbuf.size;
	MPI_Gather(&nele_l, 1, MPI_INT, &rcounts[0], 1, MPI_INT, gcm_params.gcm_root, gcm_params.gcm_comm);

	// Compute displacements as prefix sum of rcounts
	std::unique_ptr<int[]> displs;
	std::unique_ptr<giss::DynArray<SMBMsg>> rbuf;

	displs.reset(new int[num_mpi_nodes+1]);
	displs[0] = 0;
	for (int i=0; i<num_mpi_nodes; ++i) displs[i+1] = displs[i] + rcounts[i];
	int nele_g = displs[num_mpi_nodes];

	// Create receive buffer, and gather into it
	// (There's an extra item in the array for a sentinel)
	rbuf.reset(new giss::DynArray<SMBMsg>(SMBMsg::size(nfields), nele_g+1));

	MPI_Datatype mpi_type(SMBMsg::new_MPI_struct(nfields));
	MPI_Gatherv(sbuf.begin().get(), sbuf.size, mpi_type,
		rbuf->begin().get(), &rcounts[0], &displs[0], mpi_type,
		gcm_params.gcm_root, gcm_params.gcm_comm);
	MPI_Type_free(&mpi_type);

	if (am_i_root()) {
		// (ONLY ON GCM ROOT)
		// Clear output arrays, which will be filled in additively
		// on each ice model
//		for (auto ov=gcm_ivals.begin(); ov != gcm_ivals.end(); ++ov) *ov = 0;

		// Add a sentinel
		(*rbuf)[rbuf->size-1].sheetno = 999999;

		// Sort the receive buffer so items in same ice sheet
		// are found together
		qsort(rbuf->begin().get(), rbuf->size, rbuf->ele_size, &SMBMsg::compar);

		// (ONLY ON GCM ROOT)
		// Figure out which ice sheets we have data for
		auto lscan(rbuf->begin());
		auto rscan(lscan);
		std::map<int, CallIceModelParams> im_params;
		while (rscan < rbuf->end()) {
			if (rscan->sheetno != lscan->sheetno) {
				int sheetno = lscan->sheetno;
				auto cimp(CallIceModelParams(sheetno, lscan.get(), rscan.get()));
				im_params[sheetno] = cimp;
				lscan = rscan;
			}

			++rscan;
		}

		// (ONLY ON GCM ROOT)
		// NOTE: call_ice_model() is called (below) even on NON-ROOT
		// Call all our ice models
		for (int sheetno=0, auto model = ice_models.data.begin();
			model != ice_models.data.end(); ++sheetno, ++model)
		{
			// Assume we have data for all ice models
			// (So we can easily maintain MPI SIMD operation)
			auto params(im_params.find(sheetno));
			call_ice_model(&*model, sheetno, time_s, *rbuf,
				params->second.begin, params->second.next);

			// Convert to variables the GCM wants (but still on the ice grid)
			model->set_gcm_inputs(0);	// Fills in gcm_ivals_I

			// Free ice_ovals_I
			model->free_ice_ovals_I();
		}

		regrid_gcm_inputs_onroot(time_s, gcm_ivals, 0);
	} else {
		// (ONLY ON NOT GCM ROOT)
		// We're not root --- we have no data to send to ice
		// models, we just call through anyway because we will
		// receive data in an upcomming MPI_Scatter
		// Call all our ice models
		for (int sheetno=0, auto model = ice_models.data.begin();
			model != ice_models.data.end(); ++sheetno, ++model)
		{
			// Assume we have data for all ice models
			// (So we can easily maintain MPI SIMD operation)
			call_ice_model(&*model, sheetno, time_s, *rbuf,
				NULL, NULL);

		}		// if (gcm_params.gcm_rank == gcm_params.gcm_root)

	}

	printf("END GCMCoupler::couple_to_ice()\n");
}

/** Implicit input: gcm_ivals_I, set by set_gcm_inputs() */
void GCMCoupler:: regrid_gcm_inputs_onroot(
double time_s,
std::vector<giss::VectorSparseVector<int,double>> &gcm_ivals,	// Root node only: Already-allocated space to put output values.  Members as defined by the CouplingContract GCMCoupler::gcm_inputs
unsigned int mask)
{
	// This method is meant to be run only on the GCM root.
	if (!am_i_root()) return;

	printf("BEGIN GCMCoupler::regrid_gcm_inputs_onroot\n");

	// (ONLY ON GCM ROOT)
	// =============== Regrid to the grid requested by the GCM

	// (ONLY ON GCM ROOT)
	// ----------- Create the MultiMatrix used to regrid to atmosphere (I -> A)
	MultiMatrix ice2atm;
	for (auto model = models.begin(); model != models.end(); ++model) {
		// Get matrix for this single ice model.
		giss::MapSparseVector<int,double> area1_m;
		int sheetno = model.key();		// IceSheet::index
		IceSheet *sheet = (*maker_full)[sheetno];
		std::unique_ptr<giss::VectorSparseMatrix> M(
			sheet->iceinterp_to_projatm(area1_m, IceInterp::ICE));

		// Add on correction for projection
		if (maker_full->correct_area1) {
			sheet->atm_proj_correct(area1_m, ProjCorrect::PROJ_TO_NATIVE);
		}

		// Store it away...
		ice2atm.add_matrix(std::move(M), area1_m);
	}

	// (ONLY ON GCM ROOT)
	// ------------ Regrid each GCM input from ice grid to whatever grid it needs.
	for (int var_ix=0; var_ix < gcm_inputs.size_nounit(); ++var_ix) {
		giss::CoupledField const &cf(gcm_inputs.field(var_ix));

		if ((cf.flags & mask) != mask) continue;

		if ((cf.flags & contracts::GRID_BITS) == contracts::ATMOSPHERE) {
			// --- Assemble all inputs, to multiply by ice_to_hp matrix

			std::vector<blitz::Array<double, 1>> ival_I;
			for (auto model = models.begin(); model != models.end(); ++model) {
				// Assemble vector of the same GCM input variable from each ice model.
				ival_I.push_back(model->gcm_ivals_I[var_ix]);
			}

			ice2atm.multiply(ival_I, gcm_ivals[var_ix], true);
			gcm_ivals[var_ix].consolidate();

		} else if ((cf.flags & contracts::GRID_BITS) == contracts::ELEVATION) {
			// --- Assemble all inputs, to send to Glint2 QP regridding

			std::map<int, blitz::Array<double,1>> f4s;
			for (auto model = models.begin(); model != models.end(); ++model) {
				int sheetno = model.key();
				f4s.insert(std::make_pair(sheetno, model->gcm_ivals_I[var_ix]));
			}

			// Use previous return as our initial guess
			blitz::Array<double,1> initial3(maker_full->n3());
			giss::to_blitz(gcm_ivals[var_ix], initial3);

			// --------- Devise to write out this QP problem before we solve.
			// (it will actually be written out in maker_full->iceinterp_to_hp())
			std::string qpt_fname;
			if (true) {
				long time_day = (int)(time_s / 86400. + .5);
				std::stringstream fname;
				fname << time_day << "-" << cf.name << ".nc";
				boost::filesystem::path output_dir(
					gcm_params.run_dir / "qp_problems");
				boost::filesystem::create_directory(output_dir);	// Make sure it exists
				boost::filesystem::path pfname(output_dir / fname.str());

				qpt_fname = pfname.string();
			} else {
				qpt_fname = "";
			}

			// Do the regridding (involves a QP problem)
			gcm_ivals[var_ix] = maker_full->iceinterp_to_hp(
				f4s, initial3, IceInterp::ICE, QPAlgorithm::SINGLE_QP,
				qpt_fname);
			gcm_ivals[var_ix].consolidate();
		}
	}

	// ----------------- Free Memory
	// (ONLY ON GCM ROOT)
	for (auto model = models.begin(); model != models.end(); ++model) {
		model->free_ovals_ivals_I();
	}

	printf("END GCMCoupler::regrid_gcm_inputs_onroot\n");
}


/** Follows the pattern of couple_to_ice()
@param sbuf the (filled) array of ice grid values for this MPI node. */
void GCMCoupler::get_initial_state(
double time_s,
std::vector<giss::VectorSparseVector<int,double>> &gcm_ivals)	// Root node only: Already-allocated space to put output values.  Members as defined by the CouplingContract GCMCoupler::gcm_inputs
{
	printf("BEGIN GCMCoupler::get_initial_state()\n");

	if (am_i_root()) {

		// (ONLY ON GCM ROOT)
		// NOTE: call_ice_model() is called (below) even on NON-ROOT
		// Call all our ice models
		for (auto model = models.begin(); model != models.end(); ++model) {
			int sheetno = model.key();
			model->get_initial_state(time_s);

			// Record what the ice model produced.
			// NOTE: This shouldn't change model->ice_ovals_I
			IceModel_Writer *owriter = model->owriter();
			const double time_s = gcm_params.time_start_s;
			if (owriter) {
				printf("BEGIN owriter->run_timestep()\n");
				owriter->run_decoded(time_s, model->ice_ovals_I);
				printf("END owriter->run_timestep()\n");
			}

			// Convert to variables the GCM wants (but still on the ice grid)
			model->set_gcm_inputs(contracts::INITIAL);	// Fills in gcm_ivals_I

			// Free ice_ovals_I
			model->free_ice_ovals_I();
		}

		// Fill in gcm_ivals
		regrid_gcm_inputs_onroot(time_s, gcm_ivals, contracts::INITIAL);

	} else {
		// (NOT GCM ROOT)
		// We're not root --- we have no data to send to ice
		// models, we just call through anyway because we will
		// receive data in an upcomming MPI_Scatter
		// Call all our ice models
		for (auto model = models.begin(); model != models.end(); ++model) {
			model->get_initial_state(time_s);
		}		// if (gcm_params.gcm_rank == gcm_params.gcm_root)

	}

	printf("END GCMCoupler::get_initial_state()\n");
}

int GCMCoupler::rank() const
	{ return gcm_params.gcm_rank; }

}		// namespace

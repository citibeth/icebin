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

namespace glint2 {

/** @param nc The Glint22 configuration file */
void GCMCoupler::read_from_netcdf(
	NcFile &nc,
	std::string const &vname,
	std::unique_ptr<GridDomain> &&mdomain)
{
	printf("BEGIN GCMCoupler::read_from_netcdf()\n");


	// Load the MatrixMaker	(filtering by our domain, of course)
	// Also load the ice sheets
	maker.reset(new MatrixMaker(true, std::move(mdomain)));
	maker->read_from_netcdf(nc, vname);


	std::vector<std::string> const &sheet_names(maker->get_sheet_names());
    giss::MapDict<std::string, IceSheet> &sheets(maker->sheets);

	// Read gcm_out_file, an optional variable telling the GCM-specific
	// part of GLINT2 to write out exactly what it sees coming from the GCM
	// (so it can be replayed later with desm)
	NcVar *info_var = giss::get_var_safe(nc, vname + ".info");
	auto attr(giss::get_att(info_var, "gcm_out_file"));
	if (!attr.get()) {
		gcm_out_file = "";
	} else {
		gcm_out_file = attr->as_string(0);
		if (gcm_out_file.length() > 0) {
		    gcm_out_file = boost::filesystem::absolute(
				boost::filesystem::path(gcm_out_file),
				gcm_params.config_dir).string();
		}
	}

#if 1
	std::cout << "========= GCM Constants" << std::endl;

	int ix=0;
	for (auto field = gcm_constants.begin(); field != gcm_constants.end(); ++field, ++ix) {
		std::cout << "    " << *field << " = " << gcm_constants[ix] << std::endl;
	}
	std::cout << "========= GCM Outputs" << std::endl;
	std::cout << gcm_outputs << std::endl;

#endif

	int i = 0;
	for (auto name = sheet_names.begin(); name != sheet_names.end(); ++name) {
		std::string vname_sheet = vname + "." + *name;

		// Create an IceModel corresponding to this IceSheet.
		IceSheet *sheet = sheets[*name];
		models.insert(sheet->index,
			read_icemodel(*name, this, nc, vname_sheet,
				read_gcm_per_ice_sheet_params(nc, vname_sheet),
				sheet));
		IceModel *ice_model = models[i];

		// Set up the contracts specifying the variables to be passed
		// between the GCM and the ice model.  This contract is specific
		// to both the GCM and the ice model.  Note that setup_contracts()
		// is a virtual method.
		
		// This code MUST call GCMCoupler::setup_contracts() somewhere inside.
		ice_model->init(sheet->grid2, nc, vname_sheet);
		ice_model->update_ice_sheet(nc, vname_sheet, sheet);

		// Check the contract for errors
		giss::CouplingContract const &ocontract(ice_model->contract[IceModel::OUTPUT]);
		int nfields = ocontract.size_nounit();
		for (int i=0; i<nfields; ++i) {
			giss::CoupledField const &cf(ocontract.field(i));

			if (cf.grid != "ICE") {
				fprintf(stderr, "ERROR: Ice model outputs must be all on the ice grid, field %s is not\n", cf.name.c_str());
				throw std::exception();
			}
		}


		// ---- Create the affiliated writers
		// Writer for ice model input
		std::unique_ptr<IceModel_Writer> iwriter(new IceModel_Writer(*name, IceModel::INPUT, this));
		iwriter->init(sheet->grid2, ice_model, *name);
		writers[IceModel::INPUT].insert(i, std::move(iwriter));

		// Writer for ice model output
		std::unique_ptr<IceModel_Writer> owriter(new IceModel_Writer(*name, IceModel::OUTPUT, this));
		owriter->init(sheet->grid2, ice_model, *name);
		writers[IceModel::OUTPUT].insert(i, std::move(owriter));

		++i;
	}
	printf("END GCMCoupler::read_from_netcdf()\n");
}


void GCMCoupler::set_start_time(
	giss::time::tm const &time_base,
	double time_start_s)
{
	gcm_params.set_start_time(time_base, time_start_s);

	for (auto model = models.begin(); model != models.end(); ++model) {
		int sheetno = model.key();

		model->start_time_set();

		IceModel_Writer *iwriter = writers[IceModel::INPUT][sheetno];		// The affiliated input-writer (if it exists).
		if (iwriter) iwriter->start_time_set();

		IceModel_Writer *owriter = writers[IceModel::OUTPUT][sheetno];		// The affiliated input-writer (if it exists).
		if (owriter) owriter->start_time_set();

		// This function is the last phase of initialization.
		// Only now can we assume that the contracts are fully set up.
#if 1
		// Print out the contract and var transformations
		std::cout << "========= Contract for " << model->name << std::endl;
		std::cout << "---- GCM->Ice     Output Variables:" << std::endl;
		std::cout << model->contract[IceModel::INPUT];
		std::cout << model->var_transformer[IceModel::INPUT];
		std::cout << "---- Ice->GCM     Output Variables:" << std::endl;
		std::cout << model->contract[IceModel::OUTPUT];
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

/** PROTECTED method */
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
	int nfields = model->contract[IceModel::INPUT].size_nounit();

printf("BEGIN GCMCoupler::call_ice_model(nfields=%ld)\n", nfields);

	// Construct indices vector
	blitz::TinyVector<int,1> shape(rbuf.diff(end, begin));
	blitz::TinyVector<int,1> stride(rbuf.ele_size / sizeof(int));
	blitz::Array<int,1> indices(&begin->i2,
		shape, stride, blitz::neverDeleteData);

	// Construct values vectors: sparse vectors pointing into the SMBMsgs
	std::vector<blitz::Array<double,1>> ivals2;
	stride[0] = rbuf.ele_size / sizeof(double);
	for (int i=0; i<nfields; ++i) {
		SMBMsg &rbegin(*begin);		// Reference to begin

		ivals2.push_back(blitz::Array<double,1>(
			&rbegin[i], shape, stride, blitz::neverDeleteData));
	}


	// -------------- Run the model
	model->allocate0();	// Allocates ovals_I

	// Record exactly the same inputs that this ice model is seeing.
	IceModel_Writer *iwriter = writers[IceModel::INPUT][sheetno];		// The affiliated input-writer (if it exists).
	if (iwriter) iwriter->run_timestep(time_s, indices, ivals2, model->ovals_I);

	// Now call to the ice model
	model->run_timestep(time_s, indices, ivals2, model->ovals_I);

	// Record what the ice model produced
	IceModel_Writer *owriter = writers[IceModel::OUTPUT][sheetno];		// The affiliated input-writer (if it exists).
	if (owriter) owriter->run_timestep(time_s, indices, ivals2, model->ovals_I);

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

printf("[%d] BEGIN GCMCoupler::couple_to_ice() time_s=%f, sbuf.size=%d, sbuf.ele_size=%d\n", gcm_params.gcm_rank, time_s, sbuf.size, sbuf.ele_size);

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
	if (gcm_params.gcm_rank == gcm_params.gcm_root) {
		// Clear output arrays, which will be filled in additively
		// on each ice model
//		for (auto ov=gcm_ivals.begin(); ov != gcm_ivals.end(); ++ov) *ov = 0;

		// Add a sentinel
		(*rbuf)[rbuf->size-1].sheetno = 999999;

		// Sort the receive buffer so items in same ice sheet
		// are found together
		qsort(rbuf->begin().get(), rbuf->size, rbuf->ele_size, &SMBMsg::compar);

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

		// Call all our ice models
		for (auto model = models.begin(); model != models.end(); ++model) {
			int sheetno = model.key();
			// Assume we have data for all ice models
			// (So we can easily maintain MPI SIMD operation)
			auto params(im_params.find(sheetno));
			call_ice_model(&*model, sheetno, time_s, *rbuf,
				params->second.begin, params->second.next);

			// Convert to variables the GCM wants (but still on the ice grid)
			model->set_gcm_inputs();
		}

		// =============== Regrid to the grid requested by the GCM

		// ----------- Create the MultiMatrix used to regrid to atmosphere
		MultiMatrix ice2atm;
		for (auto model = models.begin(); model != models.end(); ++model) {
			// Get matrix for this single ice model.
			giss::MapSparseVector<int,double> area1_m;
			int sheetno = model.key();		// IceSheet::index
			IceSheet *sheet = (*maker)[sheetno];
			std::unique_ptr<giss::VectorSparseMatrix> M(
				sheet->iceinterp_to_projatm(area1_m, IceInterp::ICE));

			// Add on correction for projection
			if (maker->correct_area1)
				sheet->atm_proj_correct(area1_m, ProjCorrect::PROJ_TO_NATIVE);

			// Store it away...
			ice2atm.add_matrix(std::move(M), area1_m);
		}

		// ------------ Regrid each GCM input from ice grid to whatever grid it needs.
		for (int var_ix=0; var_ix < gcm_inputs.size_nounit(); ++var_ix) {
			giss::CoupledField const &cf(gcm_inputs.field(var_ix));

			if (cf.grid == "ATMOSPHERE") {
				// --- Assemble all inputs, to multiply by ice_to_hp matrix

				std::vector<blitz::Array<double, 1>> ival_I;
				for (auto model = models.begin(); model != models.end(); ++model) {
					// Assemble vector of the same GCM input variable from each ice model.
					ival_I.push_back(model->ivals_I[var_ix]);
				}

				giss::VectorSparseVector<int,double> &gcm_ivals_var_ix(gcm_ivals[var_ix]);
				ice2atm.multiply(ival_I, gcm_ivals_var_ix, true);
				gcm_ivals[var_ix].consolidate();
			} else if (cf.grid == "ELEVATION") {
				// --- Assemble all inputs, to send to Glint2 QP regridding

				std::map<int, blitz::Array<double,1>> f4s;
				for (auto model = models.begin(); model != models.end(); ++model) {
					int sheetno = model.key();
					f4s.insert(std::make_pair(sheetno, model->ivals_I[var_ix]));
				}

				// Use previous return as our initial guess
				blitz::Array<double,1> initial3(maker->n3());
				giss::to_blitz(gcm_ivals[var_ix], initial3);
				gcm_ivals[var_ix] = maker->iceinterp_to_hp(
					f4s, initial3, IceInterp::ICE, QPAlgorithm::SINGLE_QP);
				gcm_ivals[var_ix].consolidate();
			}
		}
	} else {
		// We're not root --- we have no data to send to ice
		// models, we just call through anyway because we will
		// receive data in an upcomming MPI_Scatter
		// Call all our ice models
		for (auto model = models.begin(); model != models.end(); ++model) {
			int sheetno = model.key();
			// Assume we have data for all ice models
			// (So we can easily maintain MPI SIMD operation)
			call_ice_model(&*model, sheetno, time_s, *rbuf,
				NULL, NULL);

		}		// if (gcm_params.gcm_rank == gcm_params.gcm_root)

	}
}

int GCMCoupler::rank() const
	{ return gcm_params.gcm_rank; }





}

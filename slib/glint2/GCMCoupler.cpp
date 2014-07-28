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

namespace glint2 {

/** @param glint2_config_dir Director the GLINT2 config file is found in.  Used to interpret relative paths found in that file.
@param nc The GLINT2 configuration file */
void GCMCoupler::read_from_netcdf(
	NcFile &nc,
	std::string const &vname,
	std::vector<std::string> const &sheet_names,
    giss::MapDict<std::string, IceSheet> &sheets)
{
	printf("BEGIN GCMCoupler::read_from_netcdf()\n");

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
		models.insert(i,
			read_icemodel(this, nc, vname_sheet,
				read_gcm_per_ice_sheet_params(nc, vname_sheet),
				sheet));
		IceModel *ice_model = models[i];

		// Set up the contracts specifying the variables to be passed
		// between the GCM and the ice model.  This contract is specific
		// to both the GCM and the ice model.  Note that setup_contracts()
		// is a virtual method.
		
		// THis code MUST call GCMCoupler::setup_contracts() somewhere inside.
		ice_model->init(sheet->grid2, nc, vname_sheet);
		ice_model->update_ice_sheet(nc, vname_sheet, sheet);


		// Create the affiliated writer
		std::unique_ptr<IceModel_Writer> writer(new IceModel_Writer(this));
		writer->init_from_ice_model(ice_model, *name);
		writers.insert(i, std::move(writer));


#if 1
		// Print out the contract and var transformations
		std::cout << "========= Contract for " << *name << std::endl;
		std::cout << "---- GCM->Ice     Output Variables:" << std::endl;
		std::cout << ice_model->contract[IceModel::INPUT];
		std::cout << ice_model->var_transformer[IceModel::INPUT];
		std::cout << "---- Ice->GCM     Output Variables:" << std::endl;
		std::cout << ice_model->contract[IceModel::OUTPUT];
		std::cout << ice_model->var_transformer[IceModel::OUTPUT];
#endif



		++i;
	}
	printf("END GCMCoupler::read_from_netcdf()\n");
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
	IceModel_Writer *writer,		// The affiliated input-writer (if it exists).
	double time_s,
	giss::DynArray<SMBMsg> &rbuf,
	SMBMsg *begin, SMBMsg *end)		// Messages have the MAX number of fields for any ice model contract
{
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
	std::vector<blitz::Array<double,1>> vals2;
	stride[0] = rbuf.ele_size / sizeof(double);
	for (int i=0; i<nfields; ++i) {
		SMBMsg &rbegin(*begin);		// Reference to begin

		vals2.push_back(blitz::Array<double,1>(
			&rbegin[i], shape, stride, blitz::neverDeleteData));
	}

	// Record exactly the same inputs that this ice model is seeing.
	if (writer) writer->run_timestep(time_s, indices, vals2);

	// Now call to the ice model
	model->run_timestep(time_s, indices, vals2);

printf("[%d] END GCMCoupler::call_ice_model(nfields=%ld)\n", gcm_params.gcm_rank, nfields);
};



/** @param sbuf the (filled) array of ice grid values for this MPI node. */
void GCMCoupler::couple_to_ice(
double time_s,
int nfields,			// Number of fields in sbuf.  Not all will necessarily be filled, in the case of heterogeneous ice models.
giss::DynArray<SMBMsg> &sbuf)	// Values, already converted to ice model inputs (from gcm outputs)
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
			call_ice_model(&*model, writers[sheetno], time_s, *rbuf,
				params->second.begin, params->second.next);

		}		// if (gcm_params.gcm_rank == gcm_params.gcm_root)
	} else {
		// We're not root --- we have no data to send to ice
		// models, we just call through anyway because we will
		// receive data in an upcomming MPI_Scatter
		// Call all our ice models
		for (auto model = models.begin(); model != models.end(); ++model) {
			int sheetno = model.key();
			// Assume we have data for all ice models
			// (So we can easily maintain MPI SIMD operation)
			call_ice_model(&*model, writers[sheetno], time_s, *rbuf,
				NULL, NULL);

		}		// if (gcm_params.gcm_rank == gcm_params.gcm_root)

	}
}

int GCMCoupler::rank() const
	{ return gcm_params.gcm_rank; }





}

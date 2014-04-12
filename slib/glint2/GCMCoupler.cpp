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
	int i = 0;
	for (auto name = sheet_names.begin(); name != sheet_names.end(); ++name) {
		// Create an IceModel corresponding to this IceSheet.
		IceSheet *sheet = sheets[*name];
		models.insert(i, read_icemodel(gcm_params, nc, vname + "." + *name, sheet));
		IceModel *mod = models[i];

		// Set up the contracts specifying the variables to be passed
		// between the GCM and the ice model.  This contract is specific
		// to both the GCM and the ice model.  Note that setup_contracts()
		// is a virtual method.
		setup_contracts(*mod, nc, vname + "." + *name);

#if 1
		// Print out the contract and var transformations
		std::cout << "========= Contract for " << *name << std::endl;
		std::cout << "---- GCM->Ice     Output Variables:" << std::endl;
		std::cout << mod->contract[IceModel::INPUT];
		std::cout << mod->var_transformer[IceModel::INPUT];
		std::cout << "---- Ice->GCM     Output Variables:" << std::endl;
		std::cout << mod->contract[IceModel::OUTPUT];
		std::cout << mod->var_transformer[IceModel::OUTPUT];
#endif

		// Finish initializing the IceModel.
		// This code MUST come after setup_contracts() above.
		auto const_var = nc.get_var("const");	// Physical constants
		mod->init(gcm_params, sheet->grid2, nc, vname, const_var);
		mod->update_ice_sheet(nc, vname, sheet);

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
	double time_s,
	giss::DynArray<SMBMsg> &rbuf,
	SMBMsg *begin, SMBMsg *end)		// Messages have the MAX number of fields for any ice model contract
{
	// The number of fields for THIS ice sheet will be <= the number
	// of fields in SMBMsg
	int nfields = model->contract[IceModel::INPUT].size_nounit();

printf("BEGIN call_ice_model(nfields=%ld)\n", nfields);

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

	model->run_timestep(time_s, indices, vals2);
printf("END call_ice_model(nfields=%ld)\n", nfields);
};



/** @param sbuf the (filled) array of ice grid values for this MPI node. */
void GCMCoupler::couple_to_ice(
double time_s,
int nfields,			// Number of fields in sbuf.  Not all will necessarily be filled, in the case of heterogeneous ice models.
giss::DynArray<SMBMsg> &sbuf)	// Values, already converted to ice model inputs (from gcm outputs)
{
	// Gather buffers on root node
	int num_mpi_nodes;
	MPI_Comm_size(gcm_params.gcm_comm, &num_mpi_nodes); 

printf("[%d] BEGIN couple_to_ice() time_s=%f, sbuf.size=%d, sbuf.ele_size=%d\n", gcm_params.gcm_rank, time_s, sbuf.size, sbuf.ele_size);

	// MPI_Gather the count
	std::unique_ptr<int[]> rcounts;
	rcounts.reset(new int[num_mpi_nodes]);
	for (int i=0; i<num_mpi_nodes; ++i) rcounts[i] = 0;

#if 0
for (int i=0; i<sbuf.size; ++i) {
	printf("sscan1 = %d %d %f %f\n", sbuf[i].sheetno, sbuf[i].i2, (sbuf[i])[0], (sbuf[i])[1]);
}

for (SMBMsg *sscan = sbuf.begin(); sscan < sbuf.end(); sbuf.incr(sscan)) {
	printf("sscan2 = %d %d %f %f (%p %p)\n", sscan->sheetno, sscan->i2, (*sscan)[0], (*sscan)[1], sscan, sbuf.end());
}
#endif

	int nele_l = sbuf.size;
printf("[%d] MPI_Gather sbuf.size=%ld\n", gcm_params.gcm_rank, sbuf.size);
	MPI_Gather(&nele_l, 1, MPI_INT, &rcounts[0], 1, MPI_INT, gcm_params.gcm_root, gcm_params.gcm_comm);
printf("[%d] DONE MPI_Gather\n", gcm_params.gcm_rank);

	// Compute displacements as prefix sum of rcounts
	std::unique_ptr<int[]> displs;
	std::unique_ptr<giss::DynArray<SMBMsg>> rbuf;

printf("[%d] CC\n", gcm_params.gcm_rank);
	displs.reset(new int[num_mpi_nodes+1]);
	displs[0] = 0;
	for (int i=0; i<num_mpi_nodes; ++i) displs[i+1] = displs[i] + rcounts[i];
	int nele_g = displs[num_mpi_nodes];

	// Create receive buffer, and gather into it
	// (There's an extra item in the array for a sentinel)
	rbuf.reset(new giss::DynArray<SMBMsg>(SMBMsg::size(nfields), nele_g+1));

	MPI_Datatype mpi_type(SMBMsg::new_MPI_struct(nfields));
printf("[%d] MPI_Gatherv: %p, %ld, %p, %p, %p, %p\n", gcm_params.gcm_rank, sbuf.begin(), sbuf.size, mpi_type, rbuf->begin(), &rcounts[0], &displs[0]);
	MPI_Gatherv(sbuf.begin(), sbuf.size, mpi_type,
		rbuf->begin(), &rcounts[0], &displs[0], mpi_type,
		gcm_params.gcm_root, gcm_params.gcm_comm);
printf("[%d] MPI_Gatherv DONE\n", gcm_params.gcm_rank);
	MPI_Type_free(&mpi_type);
	if (gcm_params.gcm_rank == gcm_params.gcm_root) {
		// Add a sentinel
		(*rbuf)[rbuf->size-1].sheetno = 999999;

		// Sort the receive buffer so items in same ice sheet
		// are found together
		qsort(rbuf->begin(), rbuf->size, rbuf->ele_size, &SMBMsg::compar);

#if 1
printf("[%d] BB1\n", gcm_params.gcm_rank);
int nprt=0;
for (int i=0; i<rbuf->size; ++i) {
	SMBMsg &msg((*rbuf)[i]);
//	if (msg.vals[0] != 0.0 || msg.i2 > 168860) {
		printf("    msg %d: %d %d %f\n", i, msg.sheetno, msg.i2, msg.vals[0]);
		++nprt;
		if (nprt == 10) break;
//	}
}
#endif

		// Figure out which ice sheets we have data for
printf("[%d] rbuf->begin()=%p, rbuf->end()=%p\n", gcm_params.gcm_rank, rbuf->begin(), rbuf->end());
printf("[%d] rbuf->size() = %d\n", rbuf->size);
		SMBMsg *lscan = rbuf->begin();
		SMBMsg *rscan = lscan;
		std::map<int, CallIceModelParams> im_params;
		while (rscan < rbuf->end()) {
//printf("rscan = %d %d %f %f (%p %p)\n", rscan->sheetno, rscan->i2, (*rscan)[0], (*rscan)[1], rscan, rbuf->end());
			if (rscan->sheetno != lscan->sheetno) {
printf("[%d] rscan = =========================================\n", gcm_params.gcm_rank);
				int sheetno = lscan->sheetno;
				auto cimp(CallIceModelParams(sheetno, lscan, rscan));
				im_params[sheetno] = cimp;
				lscan = rscan;
			}

			rbuf->incr(rscan);
		}
printf("[%d] BB2\n", gcm_params.gcm_rank);

		// Call all our ice models
		for (auto model = models.begin(); model != models.end(); ++model) {
			int sheetno = model.key();
printf("[%d] Calling to model sheetno=%d\n", rank(), sheetno);
			// Assume we have data for all ice models
			// (So we can easily maintain MPI SIMD operation)
			auto params(im_params.find(sheetno));
			call_ice_model(&*model, time_s, *rbuf,
				params->second.begin, params->second.next);

printf("[%d] BB\n", gcm_params.gcm_rank);
		}		// if (gcm_params.gcm_rank == gcm_params.gcm_root)
	} else {
		// We're not root --- we have no data to send to ice
		// models, we just call through anyway because we will
		// receive data in an upcomming MPI_Scatter
		// Call all our ice models
		for (auto model = models.begin(); model != models.end(); ++model) {
			int sheetno = model.key();
printf("[%d] Calling to model sheetno=%d: NULL\n", rank(), sheetno);
			// Assume we have data for all ice models
			// (So we can easily maintain MPI SIMD operation)
			call_ice_model(&*model, time_s, *rbuf,
				NULL, NULL);

printf("[%d] BB\n", gcm_params.gcm_rank);
		}		// if (gcm_params.gcm_rank == gcm_params.gcm_root)

	}
}

int GCMCoupler::rank()
	{ return gcm_params.gcm_rank; }





}

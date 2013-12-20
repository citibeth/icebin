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

/** Query all the ice models to figure out what fields they need */
std::set<IceField> GCMCoupler::get_required_fields()
{
	std::set<IceField> ret;
	for (auto model = models.begin(); model != models.end(); ++model) {
		model->get_required_fields(ret);
	}
	return ret;
}


/** @param glint2_config_dir Director the GLINT2 config file is found in.  Used to interpret relative paths found in that file.
@param nc The GLINT2 configuration file */
void GCMCoupler::read_from_netcdf(
	boost::filesystem::path const &glint2_config_dir,
	NcFile &nc,
	std::string const &vname,
	std::vector<std::string> const &sheet_names,
    giss::MapDict<std::string, IceSheet> const &sheets)
{
	int rank;
	MPI_Comm_rank(comm, &rank);

	// Only load up ice model proxies on root node
	if (rank != root) return;

	printf("BEGIN GCMCoupler::read_from_netcdf()\n");
	int i = 0;
	for (auto name = sheet_names.begin(); name != sheet_names.end(); ++name) {
		models.insert(i, read_icemodel(comm, glint2_config_dir, nc, vname + "." + *name, sheets[*name]));
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

void GCMCoupler::call_ice_model(double time_s,
	giss::DynArray<SMBMsg> &rbuf,
	std::vector<IceField> const &fields,
	SMBMsg *begin, SMBMsg *end)
{
	int nfields = fields.size();
	int sheetno = begin->sheetno;

printf("BEGIN call_ice_model(sheetno=%d, nfields=%ld)\n", sheetno, fields.size());

	// Construct indices vector
	blitz::TinyVector<int,1> shape(rbuf.diff(end, begin));
	blitz::TinyVector<int,1> stride(rbuf.ele_size / sizeof(int));
	blitz::Array<int,1> indices(&begin->i2,
		shape, stride, blitz::neverDeleteData);

	// Construct values vectors
	std::map<IceField, blitz::Array<double,1>> vals2;
	stride[0] = rbuf.ele_size / sizeof(double);
	for (int i=0; i<nfields; ++i) {
		SMBMsg &rbegin(*begin);

		vals2.insert(std::make_pair(fields[i],
			blitz::Array<double,1>(&rbegin[i],
				shape, stride, blitz::neverDeleteData)));
	}

	models[sheetno]->run_timestep(time_s, indices, vals2);
printf("END call_ice_model(sheetno=%d, nfields=%ld)\n", sheetno, fields.size());
};



/** @param sbuf the (filled) array of ice grid values for this MPI node. */
void GCMCoupler::couple_to_ice(
double time_s,
std::vector<IceField> const &fields,
giss::DynArray<SMBMsg> &sbuf)
{
	int nfields = fields.size();

	// Gather buffers on root node
	int num_mpi_nodes, rank;
	MPI_Comm_size(comm, &num_mpi_nodes); 
	MPI_Comm_rank(comm, &rank);

printf("[%d] BEGIN couple_to_ice() time_s=%f, sbuf.size=%d, sbuf.ele_size=%d\n", rank, time_s, sbuf.size, sbuf.ele_size);

	// MPI_Gather the count
	std::unique_ptr<int[]> rcounts;
//	if (rank == root)
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
printf("[%d] MPI_Gather sbuf.size=%ld\n", rank, sbuf.size);
	MPI_Gather(&nele_l, 1, MPI_INT, &rcounts[0], 1, MPI_INT, root, comm);
printf("[%d] DONE MPI_Gather\n", rank);

//if (rank == root) for (int i=0; i<num_mpi_nodes; ++i) printf("[%d] rcounts[%d] = %d\n", rank, i, rcounts[i]);

	// Compute displacements as prefix sum of rcounts
	std::unique_ptr<int[]> displs;
	std::unique_ptr<giss::DynArray<SMBMsg>> rbuf;
//	if (rank == root) {
printf("[%d] CC\n", rank);
		displs.reset(new int[num_mpi_nodes+1]);
		displs[0] = 0;
		for (int i=0; i<num_mpi_nodes; ++i) displs[i+1] = displs[i] + rcounts[i];
#if 0
if (rank == 0) {
for (int i=0; i<num_mpi_nodes; ++i) printf("[%d] rcounts[%d] = %d\n", rank, i, rcounts[i]);
for (int i=0; i<=num_mpi_nodes; ++i) printf("[%d] displs[%d] = %d\n", rank, i, displs[i]);
}
#endif
		int nele_g = displs[num_mpi_nodes];

		// Create receive buffer, and gather into it
		// (There's an extra item in the array for a sentinel)
//printf("[%d] nfields=%d, SMBMsg::size(nfields)=%ld, nele_g=%d\n", rank, nfields, SMBMsg::size(nfields), nele_g);
		rbuf.reset(new giss::DynArray<SMBMsg>(SMBMsg::size(nfields), nele_g+1));
//printf("[%d] rbuf->size = %ld\n", rank, rbuf->size);
//	}

	MPI_Datatype mpi_type(SMBMsg::new_MPI_struct(nfields));
printf("[%d] MPI_Gatherv: %p, %ld, %p, %p, %p, %p\n", rank, sbuf.begin(), sbuf.size, mpi_type, rbuf->begin(), &rcounts[0], &displs[0]);
	MPI_Gatherv(sbuf.begin(), sbuf.size, mpi_type,
		rbuf->begin(), &rcounts[0], &displs[0], mpi_type,
		root, comm);
printf("[%d] MPI_Gatherv DONE\n", rank);
	MPI_Type_free(&mpi_type);
	if (rank == root) {
		// Add a sentinel
		(*rbuf)[rbuf->size-1].sheetno = 999999;

		// Sort the receive buffer so items in same ice sheet
		// are found together
		qsort(rbuf->begin(), rbuf->size, rbuf->ele_size, &SMBMsg::compar);

		// Make a set of all the ice sheets in this model run
		std::set<int> sheets_remain;
		for (auto sheet=models.begin(); sheet != models.end(); ++sheet)
			sheets_remain.insert(sheet.key());
#if 1
printf("[%d] BB1\n", rank);
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

		// Call each ice sheet (that we have data for)
printf("[%d] rbuf->begin()=%p, rbuf->end()=%p\n", rank, rbuf->begin(), rbuf->end());
printf("[%d] rbuf->size() = %d\n", rbuf->size);
		SMBMsg *lscan = rbuf->begin();
		SMBMsg *rscan = lscan;
		while (rscan < rbuf->end()) {
//printf("rscan = %d %d %f %f (%p %p)\n", rscan->sheetno, rscan->i2, (*rscan)[0], (*rscan)[1], rscan, rbuf->end());
			if (rscan->sheetno != lscan->sheetno) {
printf("[%d] rscan = =========================================\n", rank);
				int sheetno = lscan->sheetno;
				call_ice_model(time_s, *rbuf, fields, lscan, rscan);
				sheets_remain.erase(sheetno);
				lscan = rscan;
			}
			rbuf->incr(rscan);
		}
printf("[%d] BB2\n", rank);

		// Call the ice sheets we received no data for (should not usually happen)
		for (auto sheet=sheets_remain.begin(); sheet != sheets_remain.end(); ++sheet) {
			int sheetno = *sheet;
			std::vector<int> indices;
			std::map<IceField, blitz::Array<double,1>> vals2;
			// Run with null data
			models[sheetno]->run_timestep(time_s,
				giss::vector_to_blitz(indices), vals2);
		}
printf("[%d] BB\n", rank);
	}		// if (rank == root)
}

int GCMCoupler::rank()
{
	int ret;
	MPI_Comm_rank(comm, &ret);
	return ret;
}





}

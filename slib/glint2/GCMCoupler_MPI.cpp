#include <mpi.h>

#include <cstdlib>
#include <glint2/GCMCoupler_MPI.hpp>

namespace glint2 {

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
// GCMCoupler_MPI

void GCMCoupler_MPI::call_ice_model(
	giss::DynArray<SMBMsg> &rbuf,
	std::vector<IceField> const &fields,
	SMBMsg *begin, SMBMsg *end)
{
	int nfields = fields.size();
	int sheetno = begin->sheetno;

printf("call_ice_model(sheetno=%d, nfields=%ld)\n", sheetno, fields.size());

	// Construct indices vector
	blitz::TinyVector<int,1> shape(rbuf.size);
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

	models[sheetno]->run_timestep(indices, vals2);
};



/** @param sbuf the (filled) array of ice grid values for this MPI node. */
void GCMCoupler_MPI::couple_to_ice(
std::vector<IceField> const &fields,
giss::DynArray<SMBMsg> &sbuf)
{
	int nfields = fields.size();

	// Gather buffers on root node
	int num_mpi_nodes, rank;
	MPI_Comm_size(comm, &num_mpi_nodes); 
	MPI_Comm_rank(comm, &rank);

	// MPI_Gather the count
	std::unique_ptr<int[]> rcounts;
//	if (rank == root)
		rcounts.reset(new int[num_mpi_nodes]);
#if 0
for (int i=0; i<sbuf.size; ++i) {
	printf("sscan1 = %d %d %f %f\n", sbuf[i].sheetno, sbuf[i].i2, (sbuf[i])[0], (sbuf[i])[1]);
}

for (SMBMsg *sscan = sbuf.begin(); sscan < sbuf.end(); sbuf.incr(sscan)) {
	printf("sscan2 = %d %d %f %f (%p %p)\n", sscan->sheetno, sscan->i2, (*sscan)[0], (*sscan)[1], sscan, sbuf.end());
}
#endif

	int nele_l = sbuf.size;
printf("[%d] MPI_Gather\n", rank);
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
for (int i=0; i<=num_mpi_nodes; ++i) printf("[%d] displs[%d] = %d\n", rank, i, displs[i]);
		int nele_g = displs[num_mpi_nodes];

		// Create receive buffer, and gather into it
		// (There's an extra item in the array for a sentinel)
		rbuf.reset(new giss::DynArray<SMBMsg>(SMBMsg::size(nfields), nele_g+1));
printf("[%d] CC\n", rank);
//	}
printf("[%d] AA\n", rank);

	MPI_Datatype mpi_type(SMBMsg::new_MPI_struct(nfields));
printf("[%d] MPI_Gatherv: %p, %ld, %p, %p, %p, %p\n", rank, sbuf.begin(), sbuf.size, mpi_type, rbuf->begin(), &rcounts[0], &displs[0]);
	MPI_Gatherv(sbuf.begin(), sbuf.size, mpi_type,
		rbuf->begin(), &rcounts[0], &displs[0], mpi_type,
		root, comm);
printf("[%d] MPI_Gatherv DONE\n", rank);
	MPI_Type_free(&mpi_type);
	if (rank == root) {
		// Sort the receive buffer so items in same ice sheet
		// are found together
		qsort(rbuf->begin(), rbuf->size, rbuf->ele_size, &SMBMsg::compar);

		// Add a sentinel
		(*rbuf)[rbuf->size-1].sheetno = 999999;

		// Make a set of all the ice sheets in this model run
		std::set<int> sheets_remain;
		for (auto sheet=models.begin(); sheet != models.end(); ++sheet)
			sheets_remain.insert(sheet.key());
printf("BB\n");

		// Call each ice sheet (that we have data for)
		SMBMsg *lscan = rbuf->begin();
		SMBMsg *rscan = lscan;
		while (rscan < rbuf->end()) {
//printf("rscan = %d %d %f %f (%p %p)\n", rscan->sheetno, rscan->i2, (*rscan)[0], (*rscan)[1], rscan, rbuf->end());
			if (rscan->sheetno != lscan->sheetno) {
printf("rscan = =========================================\n");
				int sheetno = lscan->sheetno;
				call_ice_model(*rbuf, fields, lscan, rscan);
				sheets_remain.erase(sheetno);
				lscan = rscan;
			}
			rbuf->incr(rscan);
		}
printf("BB\n");

		// Call the ice sheets we don't have data for
		for (auto sheet=sheets_remain.begin(); sheet != sheets_remain.end(); ++sheet) {
			int sheetno = *sheet;
			std::vector<int> indices;
			std::map<IceField, blitz::Array<double,1>> vals2;
			// Run with null data
			models[sheetno]->run_timestep(
				giss::vector_to_blitz(indices), vals2);
		}
printf("BB\n");
	}		// if (rank == root)
}

void GCMCoupler_MPI::read_from_netcdf(NcFile &nc, std::string const &vname,
	std::vector<std::string> const &sheet_names,
    giss::MapDict<std::string, IceSheet> const &sheets)
{
	int rank;
	MPI_Comm_rank(comm, &rank);

	// Only load up ice model proxies on root node
	if (rank != root) return;

	GCMCoupler::read_from_netcdf(nc, vname, sheet_names, sheets);
}



} 	// namespace glint2

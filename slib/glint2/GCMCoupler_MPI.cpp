#include <mpi.h>

#include <cstdlib>

namespace glint2 {

static void call_ice_model(
	DynArray<SMBMsg> &rbuf,
	std::vector<IceField> const &fields,
	SMBMsg *begin, SMBMsg *end)
{
	int nfields = fields.size();
	int sheetno = begin->sheetno;

	// Construct indices vector
	blitz::TinyVector<int,1> shape(rbuf.size);
	blitz::TinyVector<int,1> stride(rbuf.ele_size / sizeof(int));
	blitz::Array<int,1> indices(&begin->i2,
		shape, stride, blitz::neverDeleteData);

	// Construct values vectors
	std::vector<blitz::Array<double,1>> vals2;
	stride[0] = rbuf.ele_size / sizeof(double);
	for (int i=0; i<nfields; ++i) {
		vals2.insert(std::make_pair(fields[i],
			blitz::Array<double,1>(&(*begin)[i],
				shape, stride, blitz::neverDeleteData)));
	}

	models[sheetno].run_timestep(indices, vals2);
};



/** @param sbuf the (filled) array of ice grid values for this MPI node. */
void GCMCoupler_MPI::couple_to_ice(
std::vector<IceField> const &fields,
DynArray<SMBMsg> &sbuf)
{
	int nfields = fields.size();

	// Gather buffers on root node
	int num_mpi_nodes, rank;
	MPI::Comm::size(comm, &num_mpi_nodes); 
	MPI::Comm::rank(comm, &rank);

	// MPI_Gather the count
	std::unique_ptr<int[]> rcounts;
	if (rank == root) rcounts.reset(new int[num_mpi_nodes]);
	MPI::Gather(&nele_l, 1, MPI_INT, &rcounts[0], 1, MPI_INT, root, comm);

	// Compute displacements as prefix sum of rcounts
	std::unique_ptr<int[]> displs;
	std::unique_ptr<DynArray<SMBMsg>> rbuf;
	if (rank == root) {
		displs.reset(new int[num_mpi_nodes+1]);
		displs[0] = 0;
		for (int i=0; i<num_mpi_nodes; ++i) displs[i+1] = displs[i] + rcounts[i];
		int nele_g = displs[num_mpi_nodes];

		// Create receive buffer, and gather into it
		// (There's an extra item in the array for a sentinel)
		rbuf.reset(new DynArray<SMBMsg>(SMBMsg::size(nfields), nele_g+1));
	}

	MPI::Datatype mpi_type(SMBMsg::new_MPI_struct(nfields));
	MPI::Gatherv(sbuf.begin(), sbuf.size, mpi_type,
		rbuf->begin(), rcounts, displs, mpi_type,
		root, comm);
	mpi_type::Free();

	if (rank == root) {
		// Sort the receive buffer so items in same ice sheet
		// are found together
		qsort(rbuf.begin(), rbuf->size, rbuf->ele_size, &SMBMsg::compar);

		// Add a sentinel
		(*rbuf)[rbuf->size-1].sheetno = 999999;

		// Make a set of all the ice sheets in this model run
		std::set<int> sheets_remain;
		for (auto sheet=sheets.begin(); sheet != sheets.end(); ++sheet)
			sheets_remain.add(sheet.key());

		// Call each ice sheet (that we have data for)
		SMBMsg *lscan = rbuf->begin();
		SMBMsg *rscan = lscan;
		while (rscan < rbuf->end()) {
			if (rscan->sheetno != lscan->sheetno) {
				int sheetno = lscan->sheetno;
				call_ice_model(rbuf, fields, lscan, rscan);
				sheets_remain.remove(sheetno);
				lscan = rscan;
			}
			rbuf->incr(rscan);
		}

		// Call the ice sheets we don't have data for
		for (auto sheet=sheets_remain.begin(); sheet != sheets_remain.end(); ++sheet) {
			int sheetno = *sheet;
			std::vector<int> indices;
			std::vector<blitz::Array<double,1>> vals2;
			models[sheetno].run_timestep(indices, vals2);
		}
	}		// if (rank == root)
}

void GCMCoupler_MPI::read_from_netcdf(NcFile &nc, std::string const &vname,
	std::vector<string> const &sheet_names)
{
	int rank;
	MPI::Comm::rank(comm, &rank);

	// Only load up ice model proxies on root node
	if (rank != root) return;

	GCMCoupler::read_from_netcdf(nc, vname, sheet_names);
}



} 	// namespace glint2
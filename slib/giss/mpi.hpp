#pragma once

#include <giss/DynArray.hpp>
#include <memory>

namespace giss {

template<class MPIMsgT>
static std::unique_ptr<DynArray<MPIMsgT>> gather_msg_array(
MPI_Comm comm,
int root,		// The root rank
DynArray<MPIMsgT> &sbuf,		// Source buffer (msg on this MPI node)
int nfields,		// # of fields in MPIMsgT
int nele_l,					// # slots in source buffer
int nele_r_extra = 0)		// # extra slots to allocate in receive buffer
{
	int rank;
	MPI_Comm_rank(comm, &rank);

	// Gather buffers on root node
	int num_mpi_nodes;
	MPI_Comm_size(comm, &num_mpi_nodes); 

	// MPI_Gather the count
	std::unique_ptr<int[]> rcounts;
	rcounts.reset(new int[num_mpi_nodes]);
	for (int i=0; i<num_mpi_nodes; ++i) rcounts[i] = 0;

	MPI_Gather(&nele_l, 1, MPI_INT, &rcounts[0], 1, MPI_INT, root, comm);

	// Compute displacements as prefix sum of rcounts
	std::unique_ptr<int[]> displs;
	std::unique_ptr<DynArray<MPIMsgT>> rbuf;

	displs.reset(new int[num_mpi_nodes+1]);
	displs[0] = 0;
	for (int i=0; i<num_mpi_nodes; ++i) displs[i+1] = displs[i] + rcounts[i];
	int nele_g = displs[num_mpi_nodes];

	// Create receive buffer, and gather into it
	// (There's an extra item in the array for a sentinel)
	rbuf.reset(new DynArray<MPIMsgT>(sbuf.ele_size, nele_g+nele_r_extra));

	MPI_Datatype mpi_type(MPIMsgT::new_MPI_struct(nfields));
	MPI_Gatherv(sbuf.begin().get(), sbuf.size, mpi_type,
		rbuf->begin().get(), &rcounts[0], &displs[0], mpi_type,
		root, comm);
	MPI_Type_free(&mpi_type);


	return rbuf;
}

}

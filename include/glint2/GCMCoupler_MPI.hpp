#pragma once

#include <cstddef>	// Get offsetof()
#include <new>		// Get placement new
#include <mpi.h>
#include <cstdlib>
#include <giss/DynArray.hpp>
#include <glint2/GCMCoupler.hpp>

namespace glint2 {

struct SMBMsg {
	int sheetno;
	int i2;			// Index into ice model
	double vals[1];		// Always at least one val; but this could be extended

	double &get(int i) { return *(vals + i); }

	/** @return size of the struct, given a certain number of values */
	static size_t size(int nfields)
		{ return sizeof(SMBMsg) + (nfields-1) * sizeof(double); }

	static MPI_Datatype new_MPI_struct(int nfields);

	/** for use with qsort */
	static int compar(void const * a, void const * b);

};

// =======================================================

class GCMCoupler_MPI : public GCMCoupler {
protected :

	MPI_Comm comm;

	// Rank of the root process, which will receive data on the ice grid
	int root;

	void call_ice_model(
		giss::DynArray<SMBMsg> &rbuf,
		std::vector<IceField> const &fields,
		SMBMsg *begin, SMBMsg *end);

public :
	GCMCoupler_MPI(MPI_Comm _comm, int _root) :
		comm(_comm), root(_root) {}

	/** @param sbuf the (filled) array of ice grid values for this MPI node. */
	void couple_to_ice(
		std::vector<IceField> const &fields,
		giss::DynArray<SMBMsg> &sbuf);

	virtual void read_from_netcdf(NcFile &nc, std::string const &vname,
		std::vector<std::string> const &sheet_names);

};

}


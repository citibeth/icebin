#include <new>		// Get placement new
#include <mpi.h>
#include <cstdlib>

namespace glint2 {


struct SMBMsg {
	int sheetno;
	int i2;			// Index into ice model
	double vals[1];		// Always at least one val; but this could be extended

	/** @return size of the struct, given a certain number of values */
	static size_t size(int nfields)
		{ return sizeof(SMBMsg) + (nfields-1) * sizeof(double); }

	static MPI::Datatype new_MPI_struct(int nfields)
	{
		int nele = 2 + nfields;
		int blocklengths[] = {1, 1, nfields};
		MPI::Aint displacements[] = {offset(SMBMsg,sheetno), offset(SMBMsg,i2), offset(SMBMsg, vals)};
		MPI::Datatype types[] = {MPI_INT, MPI_INT, MPI_DOUBLE};
		return MPI::Datatype MPI::Datatype::Create_struct(3,
			blocklengths, displacements, types);
	}

	/** for use with qsort */
	static int compar(const void *a, const void *b)
	{
		SMBMsg *aa = reinterpret_cast<SMBMsg *>(a);
		SMBMsg *bb = reinterpret_cast<SMBMsg *>(b);
		int cmp = aa->sheetno - bb->sheetno;
		if (cmp != 0) return cmp;
		return aa->i2 - bb->i2;
	}

};

// =======================================================

class GCMCoupler_MPI : public GCMCoupler {
protected :

	MPI_Comm comm;

	// Rank of the root process, which will receive data on the ice grid
	int root;

public :
	public GCMCoupler_MPI(MPI_Comm _comm, int _root) :
		comm(_comm), root(_root) {}

	/** @param sbuf the (filled) array of ice grid values for this MPI node. */
	void GCMCoupler_MPI::couple_to_ice(
		std::vector<IceField> const &fields,
		DynArray<SMBMsg> &sbuf);

	virtual void read_from_netcdf(NcFile &nc, std::string const &vname,
		std::vector<string> const &sheet_names);


}

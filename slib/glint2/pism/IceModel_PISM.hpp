#pragma once

// --------------------------------
// PISM Includes... want to be included first
#include <petsc.h>
#include "IceGrid.hh"
#include "iceModel.hh"

#include "pism_options.hh"
#include "PAFactory.hh"
#include "POFactory.hh"
#include "PSFactory.hh"

#include <PISMTime.hh>
// --------------------------------
#include <mpi.h>
#include <glint2/pism/PetscContext.hpp>
#include <glint2/IceModel.hpp>
#include <glint2/IceModel_Decode.hpp>
#include <memory>
#include <glint2/pism/PSConstantGLINT2.hpp>
#include <glint2/pism/PISMIceModel.hpp>		// Our specialized subclass of PISM's ::IceModel
#include <glint2/Grid_XY.hpp>
#include <glint2/GCMParams.hpp>

namespace glint2 {
namespace gpism {

// =================================================


class IceModel_PISM : public IceModel_Decode
{
	std::shared_ptr<glint2::Grid_XY const> glint2_grid;
	MPI_Comm pism_comm;			// Commnicator used by ice model
	PetscMPIInt pism_rank, pism_size;

	/** Arguments to be passed to PISM at PISM initialization */
	std::vector<std::string> pism_args;

	/** Should we write inputs to PISM via PISM's dump() functionality?
	This is (almost certainly) supercede by the IceModel_Writer class,
	except in cases of extreme debugging. */
	const bool write_pism_inputs = false;

	// These should be declared in the same order they're created,
	// so they get destroyed in the proper reverse order.
	// See: http://msdn.microsoft.com/en-us/library/8183zf3x%28v=vs.110%29.aspx
	std::unique_ptr<PetscContext> petsc_context;
	std::unique_ptr<pism::UnitSystem> unit_system;
	std::unique_ptr<pism::Config> config;
	std::unique_ptr<pism::Config> overrides;
	std::unique_ptr<pism::IceGrid> pism_grid;
	std::unique_ptr<glint2::gpism::PISMIceModel> ice_model;
	PSConstantGLINT2 *pism_surface_model;	// We don't own this.

	// Stuff used for Scatter/Gather
	DM da2;
	Vec g2, g2natural;  //!< global Vecs used to transfer data to/from processor 0.
	VecScatter scatter; //!< VecScatter used to transfer data to/from processor 0.
	Vec Hp0;			//!< Resulting vector on process 0

	// Corresponding PISM variable for each input field
	std::vector<pism::IceModelVec2S *> pism_ivars;

	// --- Corresponding PISM variable for each output field
	// The variable, as it exists in PISM
	std::vector<pism::IceModelVec2S *> pism_ovars;
//	// An MPI-collected version of the variable
//	std::vector<blitz::Array<double,2> glint2_ovars;


	double BY_ICE_DENSITY;		// CONSTANT Used to prepare input for PISM

	/** Should we upate the elevation field in update_ice_sheet()?  Normally, yes.
	But in some TEST CASES ONLY --- when the SMB field was created with a different
	set of elevations than the ice model is using --- then this can cause problems
	in the generated SMB fields. */
	bool update_elevation = true;

	// ------------------------
public:

	/** Initialize any grid information, etc. from the IceSheet struct.
	@param vname_base Construct variable name from this, out of which to pull parameters from netCDF */
	void init(
		std::shared_ptr<glint2::Grid> const &grid2,
		NcFile &nc,
		std::string const &vname_base);

	/** Event handler to let IceModels know the start time is (finally) set */
	void start_time_set();

	int nx() { return pism_grid->Mx; }
	int ny() { return pism_grid->My; }

	// ----------------------------------------------
	/** Converts ij indices to PISM
	@param i Index in x ("horizontal") direction
	@param j Index in y ("vertical") direction, toward the poles */
	long ij_to_pism1d(int i, int j) const
		{ return i * glint2_grid->ny() + j; }

	void pism1d_to_ij(long ix1, int &i, int &j)
	{
		i = ix1 / ny();
		j = ix1 - (i * ny());
	}

	// ----------------------------------------------
	long pism1d_to_glint2(long ix_pism)
	{
		int i, j;
		pism1d_to_ij(ix_pism, i, j);
		return glint2_grid->ij_to_index(i, j);
	}

	long glint2_to_pism1d(long ix_glint2)
	{
		int i, j;
		glint2_grid->index_to_ij(ix_glint2, i, j);
		return ij_to_pism1d(i, j);
	}

	// ----------------------------------------------
	PetscErrorCode iceModelVec2S_to_blitz_xy(pism::IceModelVec2S &pism_var, blitz::Array<double,2> &ret);
	// ----------------------------------------------

	void update_ice_sheet(
		NcFile &nc,
		std::string const &vname,
		IceSheet *sheet);

	IceModel_PISM(std::string const &_name, GCMCoupler const *_coupler);

	~IceModel_PISM();

protected:

	PetscErrorCode allocate();

	PetscErrorCode deallocate();

public:
	void setup_contracts_modele();

protected:
	/** Transfers a constant from GCMCoupler::gcm_constants to PISM's configuration variable.
	Runs from within transfer_constants_xxxx() */
	void transfer_constant(std::string const &dest, std::string const &src, double multiply_by=1.0, bool set_new = false);

	/** @param set_new If true, PISM constant will be set, even if it
	was not already set in the configuration.  This defaults to false,
	as a program error check against misspelled parameter names. */
	void set_constant(std::string const &dest, double src_val, std::string const &src_units, bool set_new = false);

public:
	void run_decoded(double time_s,
		std::vector<blitz::Array<double,1>> const &vals2,
		std::vector<blitz::Array<double,1>> &ovals2);			// Output values; must be pre-allocated.
private:
	PetscErrorCode run_decoded_petsc(double time_s,
		std::vector<blitz::Array<double,1>> const &vals2);

};

}}	// namespace glint2::pism

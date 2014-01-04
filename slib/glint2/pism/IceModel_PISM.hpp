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
#include <glint2/IceModel_DISMAL.hpp>
#include <memory>
#include <glint2/pism/PSConstantGLINT2.hpp>
#include <glint2/pism/PISMIceModel.hpp>		// Our specialized subclass of PISM's ::IceModel
#include <glint2/Grid_XY.hpp>

namespace glint2 {
namespace pism {

// =================================================


class IceModel_PISM : public IceModel_Decode
{
	std::shared_ptr<glint2::Grid_XY const> glint2_grid;
	MPI_Comm pism_comm;			// Commnicator used by ice model
	PetscMPIInt pism_rank, pism_size;

	// These should be declared in the same order they're created,
	// so they get destroyed in the proper reverse order.
	// See: http://msdn.microsoft.com/en-us/library/8183zf3x%28v=vs.110%29.aspx
	std::unique_ptr<PetscContext> petsc_context;
	std::unique_ptr<::PISMUnitSystem> unit_system;
	std::unique_ptr<::NCConfigVariable> config;
	std::unique_ptr<::NCConfigVariable> overrides;
	std::unique_ptr<::IceGrid> pism_grid;
	std::unique_ptr<PISMIceModel> ice_model;
	/** We retain ownership to this because IceModel::attach_surface_model() does not take ownership. */
//	std::unique_ptr<PSConstantGLINT2> pism_surface_model;
	PSConstantGLINT2 *pism_surface_model;	// We don't own this.

	// Stuff used for Scatter/Gather
	DM da2;
	Vec g2, g2natural;  //!< global Vecs used to transfer data to/from processor 0.
	// VecScatter scatter; //!< VecScatter used to transfer data to/from processor 0.
	// Vec Hp0;			//!< Resulting vector on process 0


	// Corresponding PISM variable for each field
	std::map<IceField, IceModelVec2S *> pism_vars;

	double BY_ICE_DENSITY;		// CONSTANT Used to prepare input for PISM

	/** Use a DISMAL ice model to save stuff easily (for debugging) */
	std::unique_ptr<IceModel_DISMAL> dismal;

	// ------------------------
public:
	/** Not sure if this is used... */
	int process_options();

	/** Initialize any grid information, etc. from the IceSheet struct.
	@param vname_base Construct variable name from this, out of which to pull parameters from netCDF */
	void init(
		IceModel::GCMParams const &gcm_params,
		std::shared_ptr<glint2::Grid> const &grid2,
		NcFile &nc,
		std::string const &vname_base,
		NcVar *const_var);

	~IceModel_PISM();

	/** Query all the ice models to figure out what fields they need */
	void get_required_fields(std::set<IceField> &fields);

	PetscErrorCode allocate(
		std::shared_ptr<const glint2::Grid_XY> &,
		NcVar *pism_var, NcVar *const_var);

	PetscErrorCode deallocate();

#if 0
public:
	/** Inherited from IceModel */
	void run_timestep(double time_s,
		blitz::Array<int,1> const &indices,
		std::map<IceField, blitz::Array<double,1>> const &vals2);
private:
	PetscErrorCode run_timestep_petsc(double time_s,
		blitz::Array<int,1> const &indices,
		std::map<IceField, blitz::Array<double,1>> const &vals2);

#else
public:
	void run_decoded(double time_s,
		std::map<IceField, blitz::Array<double,1>> const &vals2);
private:
	PetscErrorCode run_decoded_petsc(double time_s,
		std::map<IceField, blitz::Array<double,1>> const &vals2);

#endif

};

}}	// namespace glint2::pism

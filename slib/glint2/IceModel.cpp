#include <mpi.h>	// Must be first
#include <glint2/IceModel.hpp>
#include <glint2/GCMCoupler.hpp>

namespace glint2 {

IceModel::IceModel(IceModel::Type _type, std::string const &_name, GCMCoupler const *_coupler)
	: type(_type), name(_name), coupler(_coupler), ice_constants(&_coupler->ut_system)
//	contract({giss::CouplingContract(), giss::CouplingContract()})
{
}

IceModel::~IceModel() {}


giss::CouplingContract *IceModel::new_CouplingContract() {
	_extra_contracts.push_back(
		std::unique_ptr<giss::CouplingContract>(
		new giss::CouplingContract()));
	return _extra_contracts[_extra_contracts.size()-1].get();
}

// ==========================================================

/** Allocate vectors in preparation of calling an ice model. */
void IceModel::allocate0()
{
	// Allocate for direct output from ice model
	int nfields = ocontract->size_nounit();
	for (int i=0; i < nfields; ++i) {
		CoupledField &cf(ocontract->field(i));
		std::string const &grid(cf.get_grid());
		long n2 = model->ndata();
		ovals_I.push_back(blitz::Array<double,1>(n2));
	}


/** Allocate in preparation of var transformations (but not regridding yet) */
void IceModel::allocate1()
{

	CouplingContract const *icontract = all->icontract;
	int nfields = icontract->size_nounit();
	for (int i=0; i < nfields; ++i) {
		CoupledField &cf(icontract->field(i));
		std::string const &grid(cf.get_grid());
		long n2 = model->ndata();
		ivals_I.push_back(blitz::Array<double,1>(n2));
	}
}

/** Free portions not needed after finished calling ice model and
applying variable transform.  This will be variables desired on
anything other than the ELEVATION grid. */
void IceModel::free1()
{
	ovals_I.clear();

	int nfields = icontract->size_nounit();
	for (int i=0; i < nfields; ++i) {
		CoupledField &cf(icontract->field(i));
		std::string const &grid(cf.get_grid());

		if (grid != "ELEVATION") giss::free_array(ivals_I(i));
	}
}

/** Free all memory used by this.  Called when we're done with a coupling timestep. */
void IceModel::free0()
{
	ovals_I.clear();
	ivals_I.clear();
}

// -----------------------------------------------------------
/** Allocates and sets ivals_I variable */
void IceModel::set_gcm_inputs()
{
	allocate1();		// Allocate ivals_I

	// Compute the variable transformation
	giss::VarTransformer &vt(var_transformer[IceModel::OUTPUT]);
	giss::CSRAndUnits trans = vt.apply_scalars({
//		std::make_pair("by_dt", 1.0 / ((itime - api->itime_last) * api->dtsrc)),
		std::make_pair("unit", 1.0)});

	// Apply the variable transformation
	for (int xi=0; xi<vt.dimension(giss::VarTransformer::OUTPUTS).size_nounit(); ++xi) {	// xi is index of output variable

		// Consider each output variable separately...
		std::vector<std::pair<int, double>> const &row(trans.mat[xi]);
		for (auto xjj=row.begin(); xjj != row.end(); ++xjj) {
			int xj = xjj->first;		// Index of input variable
			double io_val = xjj->second;	// Amount to multiply it by
			
			ivals_ += ovals2[xj] * io_val;		// blitz++ vector operation
		}
	}

}

#if 0
		// -------- Regrid just the things going to ATMOSPHERE grid.
		CoupledField &cf(model->contract[IceModel::OUTPUT].field(xi));



		// Regrid to the request grid for the GCM
		CoupledField &cf(model->contract[IceModel::OUTPUT].field(xi));
		if (cf.grid == "ICE") {
			// PASS: grid=ICE variables are just for internal Glint2 consumption.
		} else if (cf.grid == "ELEVATION") {
			// PASS: Elevation stuff from all ice sheets must be regridded at once.
		} else if (cf.grid == "ATMOSPHERE") {
			auto mat(model->maker->iceinterp_to_projatm(

	}

#endif



}
// -----------------------------------------------------------

}

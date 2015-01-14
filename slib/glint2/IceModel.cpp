#include <mpi.h>	// Must be first
#include <glint2/IceModel.hpp>
#include <glint2/IceModel_Writer.hpp>
#include <glint2/GCMCoupler.hpp>

namespace glint2 {

IceModel::IceModel(IceModel::Type _type, std::string const &_name, GCMCoupler const *_coupler)
	: type(_type), name(_name), coupler(_coupler), ice_constants(&_coupler->ut_system)
//	contract({giss::CouplingContract(), giss::CouplingContract()})
{
}

IceModel::~IceModel() {}


void IceModel::set_writers(
	std::unique_ptr<IceModel_Writer> &&iwriter,
	std::unique_ptr<IceModel_Writer> &&owriter)
{
	_iwriter = std::move(iwriter);
	_owriter = std::move(owriter);
}


giss::CouplingContract *IceModel::new_CouplingContract() {
	_extra_contracts.push_back(
		std::unique_ptr<giss::CouplingContract>(
		new giss::CouplingContract()));
	return _extra_contracts[_extra_contracts.size()-1].get();
}

// ==========================================================
bool IceModel::am_i_root() const
	{ return coupler->am_i_root(); }

/** Allocate vectors in preparation of calling an ice model. */
void IceModel::allocate_ice_ovals_I()
{
	// Check for program errors
	if (!coupler->am_i_root()) {
		fprintf(stderr, "IceModel::allocate_ice_ovals_I() should only be called from GCM root MPI node.  Fix the code.\n");
		throw std::exception();
	}
	if (ice_ovals_I.size() != 0) {
		fprintf(stderr, "[%d] IceModel::allocate_ice_ovals_I(): called twice without a free() inbetween.  Fix the code. (old size is %ld)\n", coupler->gcm_params.gcm_rank, ice_ovals_I.size());
		throw std::exception();
	}

	// Allocate for direct output from ice model
	giss::CouplingContract const &ocontract(contract[IceModel::OUTPUT]);
	int nfields = ocontract.size_nounit();
	for (int i=0; i < nfields; ++i) {
		giss::CoupledField const &cf(ocontract.field(i));
		long n2 = ndata();
		ice_ovals_I.push_back(blitz::Array<double,1>(n2));
	}
}


/** Allocate in preparation of var transformations (but not regridding yet) */
void IceModel::allocate_gcm_ivals_I()
{
	// Check for program errors
	if (!coupler->am_i_root()) {
		fprintf(stderr, "IceModel::allocate_ice_ivals_I() should only be called from GCM root MPI node.  Fix the code.\n");
		throw std::exception();
	}
	if (gcm_ivals_I.size() != 0) {
		fprintf(stderr, "IceModel::allocate_gcm_ivals_I(): called twice without a free() inbetween.  Fix the code.\n");
		throw std::exception();
	}

	giss::CouplingContract const &gcm_inputs(coupler->gcm_inputs);
	int nfields = gcm_inputs.size_nounit();
	for (int i=0; i < nfields; ++i) {
		giss::CoupledField const &cf(gcm_inputs.field(i));
		long n2 = ndata();
		gcm_ivals_I.push_back(blitz::Array<double,1>(n2));
	}
}

/** Free portions not needed after finished calling ice model and
applying variable transform.  This will be variables desired on
anything other than the ELEVATION grid. */
void IceModel::free_ice_ovals_I()
{
	// Check for program errors
	if (!coupler->am_i_root()) {
		fprintf(stderr, "IceModel::free_ice_ovals_I() should only be called from GCM root MPI node.  Fix the code.\n");
		throw std::exception();
	}

	ice_ovals_I.clear();
}

/** Free all memory used by this.  Called when we're done with a coupling timestep. */
void IceModel::free_ovals_ivals_I()
{
	// Check for program errors
	if (!coupler->am_i_root()) {
		fprintf(stderr, "IceModel::free_ovals_ovals_I() should only be called from GCM root MPI node.  Fix the code.\n");
		throw std::exception();
	}

	ice_ovals_I.clear();
	gcm_ivals_I.clear();
}

// -----------------------------------------------------------
/** Allocates and sets gcm_ivals_I variable
@param mask Control which GCM input variables to set (according to, etc, INITIAL flag in contract) */
void IceModel::set_gcm_inputs(unsigned int mask)
{
	printf("[%d] BEGIN IceModel::set_gcm_inputs()\n", coupler->gcm_params.gcm_rank);
	allocate_gcm_ivals_I();

	// Compute the variable transformation
	giss::VarTransformer &vt(var_transformer[IceModel::OUTPUT]);
	giss::CSRAndUnits trans = vt.apply_scalars({
		std::make_pair("unit", 1.0)});

	giss::CouplingContract const &gcm_inputs(coupler->gcm_inputs);

	// Apply the variable transformation
	for (int xi=0; xi<vt.dimension(giss::VarTransformer::OUTPUTS).size_nounit(); ++xi) {	// xi is index of output variable
		giss::CoupledField const &cf(gcm_inputs.field(xi));

		if ((cf.flags & mask) != mask) continue;

		// Consider each output variable separately...
		std::vector<std::pair<int, double>> const &row(trans.mat[xi]);
		for (auto xjj=row.begin(); xjj != row.end(); ++xjj) {
			int xj = xjj->first;		// Index of input variable
			double io_val = xjj->second;	// Amount to multiply it by
			gcm_ivals_I[xi] += ice_ovals_I[xj] * io_val;		// blitz++ vector operation
		}
		gcm_ivals_I[xi] += trans.units[xi];
	}

	printf("END IceModel::set_gcm_inputs()\n");
}

// -----------------------------------------------------------

}	// namespace glint2

#include <cstdio>
#include <unordered_set>

#include <spsparse/netcdf.hpp>
#include <spsparse/multiply_sparse.hpp>

#include <icebin/GCMRegridder.hpp>
#include <icebin/IceRegridder_L0.hpp>

using namespace netCDF;
using namespace ibmisc;
using namespace std::placeholders;  // for _1, _2, _3...
using namespace spsparse;

namespace icebin {

// -----------------------------------------------------
IceRegridder::IceRegridder() : interp_style(InterpStyle::Z_INTERP), name("icesheet") {}

IceRegridder::~IceRegridder() {}


std::unordered_map<long,double> IceRegridder::elevI_hash()
{
	// Convert elevI to a hash table, so we can look up in it easily
	std::unordered_map<long,double> elevIh;
	for (auto ii=elevI.begin(); ii != elevI.end(); ++ii)
		elevIh.insert(std::make_pair(ii.index(0), ii.val()));
	return elevIh;
}
// --------------------------------------------------------

// ==============================================================
// Write out the parts that this class computed --- so we can test/check them
// -------------------------------------------------------------
/** Produces the diagonal matrix [Atmosphere projected] <-- [Atmosphere]
NOTE: wAvAp == sApvA */
void IceRegridder::wAvAp(SparseVector &w)
{
	ibmisc::Proj_LL2XY proj(gridI->sproj);
	for (auto cell=gcm->gridA->cells.begin(); cell != gcm->gridA->cells.end(); ++cell) {
		w.add({cell->index}, cell->native_area / cell->proj_area(&proj));
	}
}

/** Produces the diagonal matrix [Atmosphere projected] <-- [Atmosphere]
NOTE: wAvAp == sApvA */
void IceRegridder::wEvEp(SparseVector &w)
{
	ibmisc::Proj_LL2XY proj(gridI->sproj);
	for (auto cell=gcm->gridA->cells.begin(); cell != gcm->gridA->cells.end(); ++cell) {
		long nhp = gcm->nhp(cell->index);
		long tuple[2] = {cell->index, 0};
		long &ihp(tuple[1]);
		for (ihp=0; ihp<nhp; ++ihp) {
			long indexE = gcm->indexingHP.tuple_to_index(tuple);
			w.add({indexE}, cell->native_area / cell->proj_area(&proj));
		}
	}
}
// -------------------------------------------------------------
void IceRegridder::clear()
{
	gridI.reset();
	exgrid.reset();
	elevI.clear();
}
// -------------------------------------------------------------
void IceRegridder::ncio(NcIO &ncio, std::string const &vname)
{
	if (ncio.rw == 'r') clear();

	auto info_v = get_or_add_var(ncio, vname + ".info", netCDF::ncInt64, {});
	get_or_put_att(info_v, ncio.rw, "name", name);
	get_or_put_att_enum(info_v, ncio.rw, "interp_style", interp_style);

	gridI->ncio(ncio, vname + ".gridI");
	exgrid->ncio(ncio, vname + ".exgrid");
	ncio_spsparse(ncio, elevI, true, vname + ".elevI");
}
// ========================================================
void GCMRegridder::init(
	std::unique_ptr<Grid> &&_gridA,
	ibmisc::Domain<int> &&_domainA,		// Tells us which cells in gridA to keep...
	std::vector<double> &&_hpdefs,	// [nhp]
	ibmisc::Indexing<long,long> &&_indexingHP,
	bool _correctA)
{
	gridA = std::move(_gridA);
	domainA = std::move(_domainA);
	hpdefs = std::move(_hpdefs);
	indexingHP = std::move(_indexingHP);
	correctA = _correctA;

	if (indexingHP.rank() != 2) (*icebin_error)(-1,
		"indexingHP has rank %d, it must have rank=2", indexingHP.rank());
}
// -------------------------------------------------------------
std::unique_ptr<IceRegridder> new_ice_regridder(IceRegridder::Type type)
{
	switch(type.index()) {
		case IceRegridder::Type::L0 :
			return std::unique_ptr<IceRegridder>(new IceRegridder_L0);
		break;
		default :
			(*icebin_error)(-1,
				"Unknown IceRegridder::Type %s", type.str());
		break;
	}
}
// -------------------------------------------------------------
// ==============================================================

void GCMRegridder::clear()
{
	gridA.reset();
	hpdefs.clear();
	sheets.clear();
}

void GCMRegridder::ncio(NcIO &ncio, std::string const &vname)
{
	if (ncio.rw == 'r') clear();
	auto info_v = get_or_add_var(ncio, vname + ".info", netCDF::ncInt64, {});
	// get_or_put_att_enum(info_v, ncio.rw, "type", type);

	gridA->ncio(ncio, vname + ".gridA");
	indexingHP.ncio(ncio, ncInt, vname + ".indexingHP");
	ncio_vector(ncio, hpdefs, true, vname + ".hpdefs", ncDouble,
		get_or_add_dims(ncio, {vname + ".nhp"}, {hpdefs.size()} ));
	domainA.ncio(ncio, ncInt, vname + ".domainA");
	get_or_put_att(info_v, ncio.rw, vname + ".correctA", &correctA, 1);


	// Read/Write the Main Sheets
	if (ncio.rw == 'r') {
		std::vector<std::string> sheet_names;
		get_or_put_att(info_v, ncio.rw, vname + ".sheets", ncString, sheet_names);

		for (auto sheet_name : sheet_names) {
			std::string sheet_vname = vname + "." + sheet_name;
			IceRegridder::Type sheet_type;
			auto sheet_info_v = get_or_add_var(ncio, sheet_vname + ".info", netCDF::ncInt64, {});

			get_or_put_att_enum(sheet_info_v, ncio.rw, "type", sheet_type);
			std::unique_ptr<IceRegridder> sheet(new_ice_regridder(sheet_type));
			sheet->ncio(ncio, vname + "." + sheet_name);
			add_sheet(std::move(sheet));
		}
	} else {	// Write
		for (auto ii=sheets.begin(); ii != sheets.end(); ++ii) {
			IceRegridder *sheet = &*(ii->second);
			sheet->ncio(ncio, vname + "." + sheet->name);
		}
	}
}
// -------------------------------------------------------------

// ---------------------------------------------------------------------
/** Made for binding... */
static bool in_good(std::unordered_set<long> const *set, long index_c)
{
	return (set->find(index_c) != set->end());
}

void IceRegridder::filter_cellsA(std::function<bool (long)> const &useA)
{

  // Figure out which cells to keep

	// List of cells in gridI / exgrid that overlap a cell we want to keep
	std::unordered_set<long> good_index_gridI;
	std::unordered_set<long> good_index_exgrid;


	std::unordered_set<int> good_j;
	for (auto excell = exgrid->cells.begin(); excell != exgrid->cells.end(); ++excell) {
		int index1 = excell->i;
		if (useA(index1)) {
			good_index_gridI.insert(excell->j);
			good_index_exgrid.insert(excell->index);
		}
	}

	// Remove unneeded cells from gridI
	gridI->filter_cells(std::bind(&in_good, &good_index_gridI, _1));
	exgrid->filter_cells(std::bind(&in_good, &good_index_exgrid, _1));
}

void GCMRegridder::filter_cellsA(std::function<bool(long)> const &keepA)
{

	// Now remove cells from the exgrids and gridIs that
	// do not interact with the cells we've kept in grid1.
	for (auto sheet=sheets.begin(); sheet != sheets.end(); ++sheet) {
		sheet->second->filter_cellsA(keepA);
	}

	gridA->filter_cells(keepA);
}

void GCMRegridder::wA(SparseVector &w) const
{
	for (auto cell=gridA->cells.begin(); cell != gridA->cells.end(); ++cell)
		w.add({cell->index}, cell->native_area);
}
// ---------------------------------------------------------------------
void IceRegridder::init(
	std::string const &_name,
	std::unique_ptr<Grid> &&_gridI,
	std::unique_ptr<Grid> &&_exgrid,
	InterpStyle _interp_style,
	SparseVector &&elevI)
{
	name = (_name != "" ? _name : gridI->name);
	gridI = std::move(_gridI);
	exgrid = std::move(_exgrid);
	interp_style = _interp_style;
	elevI = std::move(elevI);
}
// ==============================================================
/** Gives weights for linear interpolation with a bunch of points.
If our point is off the end of the range, just continue the slope in extrapolation.
@param xpoints This is not blitz::Array<double,1> because Blitz++ does not (yet) implement STL-compatibl
e iterators. */
extern void linterp_1d(
	std::vector<double> const &xpoints,
	double xx,
	int *indices, double *weights)  // Size-2 arrays
{
	int n = xpoints.size();

	// This is the point ABOVE our value.
	// (i0 = i1 - 1, xpoints[i0] < xx <= xpoints[i1])
	// See: http://www.cplusplus.com/reference/algorithm/lower_bound/
	int i1 = lower_bound(xpoints.begin(), xpoints.end(), xx) - xpoints.begin();

	if (i1 <= 0) i1 = 1;
	if (i1 >= n) i1 = n-1;

	int i0 = i1-1;
	indices[0] = i0;
	indices[1] = i1;
	double ratio = (xx - xpoints[i0]) / (xpoints[i1] - xpoints[i0]);
	weights[0] = (1.0 - ratio);
	weights[1] = ratio;
}

// =================================================================





// ----------------------------------------------------------------
std::unique_ptr<SparseMatrix> new_regrid(
	std::array<int, SparseMatrix::rank> shape,
	std::function<void(SparseMatrix &)> const &fill_regrid)
{
	std::unique_ptr<SparseMatrix> M(new SparseMatrix);
	M->set_shape(shape);
	fill_regrid(*M);

	/* Since all raw matrices going into this have a destination G
	(exchange grid), they must all be sorted row-major.  The
	transpose matrices all go FROM G, so they must be sorted
	column-major.  But sorting the transpose of a matrix
	column-major is the same as sorting the original matrix row
	major.  So... all matrices must be sorted the same way
	(row-major) */
	M->consolidate({0,1});

	return M;
}

std::unique_ptr<SparseVector> new_diag(
	std::array<int, SparseVector::rank> shape,
	std::function<void(SparseVector &)> const &fill_diag)
{
	std::unique_ptr<SparseVector> M(new SparseVector);
	M->set_shape(shape);
	fill_diag(*M);
	M->consolidate({0});
	return M;
}


std::unique_ptr<SparseMatrix> new_regrid_with_elevI(
	IceRegridder *sheet,
	std::array<int, SparseMatrix::rank> shape,
	std::function<void(SparseMatrix &, std::unordered_map<long,double> const &elevIh)> const &fill_regrid)
{
	std::unordered_map<long,double> elevIh(sheet->elevI_hash());

	std::unique_ptr<SparseMatrix> M(new SparseMatrix);
	M->set_shape(shape);
	fill_regrid(*M, elevIh);

	/* Since all raw matrices going into this have a destination G
	(exchange grid), they must all be sorted row-major.  The
	transpose matrices all go FROM G, so they must be sorted
	column-major.  But sorting the transpose of a matrix
	column-major is the same as sorting the original matrix row
	major.  So... all matrices must be sorted the same way
	(row-major) */
	mGvZ->consolidate({0,1});

	return M;
}



SparseMatrix *get_regrid(
	RegridMatrices *rm,
	std::string const &urname)
{ return &*(rm->regrids.at(urname).M); }

SparseVector *get_diag(
	RegridMatrices *rm,
	std::string const &urname)
{ return &*(rm->diags.at(urname)); }


std::unique_ptr<SparseVector> new_weight(
	RegridMatrices *rm,
	std::string const &urname,
	int dim)	// Dimension to sum over
{
	SparseMatrix *M = get_regrid(rm, urname);

	std::unique_ptr<SparseVector> weight(new SparseVector);
	weight.set_shape({M->shape[dim]});

	for (auto ii=M->begin(); ii != M->end(); ++ii) {
		weight->add({ii.index(dim)}, ii.val());
	}
	weight->consolidate({0});
}

std::unique_ptr<SparseVector> new_invert(
	RegridMatrices *rm,
	std::string const &urname)
{
	SparseVector *V = get_diag(rm, urname);

	std::unique_ptr<SparseVector> scale(new SparseVector);
	scale.set_shape({V->shape[0]});

	for (auto ii=M->begin(); ii != M->end(); ++ii) {
		scale->add({ii.index(dim)}, 1./ii.val());
	}
	return scale;
}

void RegridMatrices::add_regrid(
	std::string const &G,
	std::string const &Z,
	std::function<std::unique_ptr<SparseMatrix>()> const &new_GvZ)
{
	// Make sure destination space is G
	if (G != "G") (*icebin_error)(-1,
		"Destination vector space for add_GvZ must be \"G\" (exchnage grid)!");

	std::string GvZ = G+"v"+Z;
	std::string ZvG = Z+"v"+G;
	regrids[GvZ] = make_transpose(LazyPtr<SparseMatrix>(new_GvZ), false);
	regrids[ZvG] = make_transpose(
		LazyPtr<SparseMatrix>(std::bind(&get_regrid, this, GvZ)), true);

	diags["w"+GvZ] = LazyPtr<SparseVector>(std::bind(&new_weight, this, GvZ, 0));
	diags["w"+ZvG] = LazyPtr<SparseVector>(std::bind(&new_weight, this, GvZ, 1));
	diags["s"+GvZ] = LazyPtr<SparseVector>(std::bind(&new_invert, this, "w"+GvZ));
	diags["s"+ZvG] = LazyPtr<SparseVector>(std::bind(&new_invert, this, "w"+ZvG));
}
// ----------------------------------------------------------------
void RegridMatrices::add_scale(
	std::string const &B,
	std::string const &A,
	std::function<std::unique_ptr<SparseVector>()> const &scale_fn)
{
	std::string BvA = B+"v"+A;
	std::string AvB = A+"v"+B;
	diags["s"+BvA] = LazyPtr<SparseVector>(scale_fn);
	diags["w"+BvA] = LazyPtr<SparseVector>(std::bind(&new_invert, "s"+BvA));

	diags["w"+AvB] = LazyPtr<SparseVector>(std::bind(&get_diag, "s"+BvA));
	diags["s"+AvB] = LazyPtr<SparseVector>(std::bind(&get_diag, "w"+BvA));
}
// ----------------------------------------------------------------
std::unique_ptr<SparseMatrix> compose_regrid(
	RegridMatrices *rm,
	std::array<std::string, 5> const &spec)
{
	Transpose<LazyPtr<SparseMatrix>> const &A = rm->regrids.at(spec[1]);
	Transpose<LazyPtr<SparseMatrix>> const &B = rm->regrids.at(spec[3]);

	SparseVector *scalei = spec[0] == "" ? NULL : &*rm->diags.at(spec[0]);
	SparseVector *scalej = spec[2] == "" ? NULL : &*rm->diags.at(spec[2]);
	SparseVector *scalek = spec[4] == "" ? NULL : &*rm->diags.at(spec[4]);

	std::unique_ptr<SparseMatrix> M(new SparseMatrix);
	M->set_shape({A.M->shape[0], B.M->shape[1]});
	multiply(*M, 1.0,
		scalei,
		*A.M, A.transpose,
		scalej,
		*B.M, B.transpose,
		scakek);
	return M;
}

Transpose<LazyPtr<SparseMatrix>> RegridMatrix::add_compose(
	std::string const &spec_basename,
	std::vector<std::array<std::string, 2>> const &scale_variants,
	std::array<std::string, 5> spec)
{
	for (auto &variant : variants) {
		std::string spec_name = spec_basename + "(" + variant[0] + ")";
		spec[0] = variant[1];

		regrids[spec_name] = make_transpose(
			LazyPtr<SparseMatrix>(std::bind(&compose_regrid, this, spec)),
			false);
	}
}

RegridMatrices::add_weight(
	std::string const &wName,
	std::array<int, 1> shape,
	std::function<void(SparseVector &)> const &fill_diag)
{
	std::string sName = wName;
	sName[0] = "s";

	diags[wName] = LazyPtr<SparseVector>(std::bind(&new_diag, shape, fill_diag));
	diags[sName] = LazyPtr<SparseVector>(std::bind(&new_invert, this, wName));
}

// ----------------------------------------------------------------
RegridMatrices::RegridMatrices(IceRegridder *_sheet)
{
	sheet = _sheet;

	std::unordered_map<long,double> elevIh(sheet->elevI_hash());

	// ----- Set up Ur matrices
	add_regrid("G", "Ep", 
		std::bind(&new_regrid_with_elevI, sheet, {sheet->nG(), sheet->gcm->nE()},
			std::bind(&IceRegridder::GvEp_noweight, sheet, _1, _2)));
		
	add_regrid("G", "I", 
		std::bind(&new_regrid_with_elevI, sheet, {sheet->nG(), sheet->nI()},
			std::bind(&IceRegridder::GvI_noweight, sheet, _1, _2)));

	add_regrid("G", "Ap",
		std::bind(&new_regrid, {sheet->nG(), sheet->gcm->nA()},
			std::bind(&IceRegridder::GvAp_noweight, sheet, _1)));

	// ----- Set up computed scaling matrices
	add_weight("wAvG", std::bind(&IceRegridder::wAvG, sheet, _1));
	add_weight("wEvG", std::bind(&IceRegridder::wAvG, sheet, _1));
	add_weight("wA", std::bind(GCMRegridder::wA, sheet->gcm->wA, _1));

	// ----- Set up computed regrid matrices
	std::vector<std::array<std::string,2>> A_variants =
		{{"NONE", ""}, {"PARTIAL_CELL", "sAvG"}, {"FULL_CELL", "sA"}};
	std::vector<std::array<std::string,2>> E_variants =
		{{"NONE", ""}, {"PARTIAL_CELL", "sEvG"}};
	std::vector<std::array<std::string,2>> I_variants =
		{{"NONE", ""}, {"PARTIAL_CELL", "sIvG"}};

	add_compose("AvI", A_variants, {"", "ApvG", "sGvI",  "GvI",  ""});
	add_compose("EvI", E_variants, {"", "EpvG", "sGvI",  "GvI",  ""});
	add_compose("IvA", I_variants, {"", "IvG",  "sGvAp", "GvAp", "sApvA"});
	add_compose("IvE", I_variants, {"", "IvG",  "sGvEp", "GvEp", "sEpvE"});

	add_compose("AvE", A_variants, {"", "EpvG", "sGvAp", "GvAp", "sApvA"});
	add_compose("EvA", E_variatns, {"",     "ApvG", "sGvEp", "GvEp", "sEpvE"});



	// ----- Show what we have!
	printf("Regrids:");
	for (auto ii = regrids.begin(); ii != regrids.end(); ++ii) {
		printf(" %s", ii->first.c_str());
	}
	printf("\n");

	printf("Diags:");
	for (auto ii = diags.begin(); ii != diags.end(); ++ii) {
		printf(" %s", ii->first.c_str());
	}
	printf("\n");

}

void RegridMatrices::wAvG(SparseVector &ret) const {
	SparseVector wAvAp;
	wAvAp.set_shape({sheet->gcm->nA()});
	sheet->wAvAp(wAvAp);
	multiply_ele(ret, wAvAp, *diags.at("wApvG"));
}
void RegridMatrices::wEvG(SparseVector &ret) const {
	SparseVector wEvEp;
	wEvEp.set_shape({sheet->gcm->nE()});
	sheet->wEvEp(wEvEp);
	multiply_ele(ret, wEvEp, *diags.at("wEpvG"));
}



}	// namespace icebin

#include <algorithm>
#include <ibmisc/netcdf.hpp>
#include <ibmisc/cython.hpp>
#include "icebin_cython.hpp"

using namespace ibmisc;
using namespace ibmisc::cython;

namespace icebin {
namespace cython {

void GCMRegridder_init(GCMRegridder *cself,
	std::string const &gridA_fname,
	std::string const &gridA_vname,
	std::vector<double> &hpdefs,
	bool _correctA)
{
	// Read gridA
	NcIO ncio(gridA_fname, netCDF::NcFile::read);
	std::unique_ptr<Grid> gridA = new_grid(ncio, gridA_vname);
	gridA->ncio(ncio, gridA_vname);
	ncio.close();

	// Construct a domain to be the full extent of indexing.
	int rank = gridA->indexing.base.size();
	std::vector<int> high(rank);
	for (int k=0; k<rank; ++k)
		high[k] = gridA->indexing.base[k] + gridA->indexing.extent[k];


	// Put it together
	long nhp = hpdefs.size();
	cself->init(
		std::move(gridA),
		Domain<int>(
			std::vector<int>(gridA->indexing.base), std::move(high)),
		std::move(hpdefs),
		Indexing<long,long>(
			{0,0}, {gridA->ndata(), nhp}, {1,0}),
		_correctA);

}


#if 0
template<class T>
static blitz::Array<T,1> np_to_blitz_1d(
PyObject *ovec,
std::string const &vname,
std::array<int,1> dims)
{
	// Check that it's type PyArrayObject
	if (!PyArray_Check(ovec)) {
		(*icebin_error)(-1,
			"check_dimensions: Object %s is not a Numpy array", vname.c_str());
	}
	PyArrayObject *vec = (PyArrayObject *)ovec;


	// Set up shape and strides
    int const T_size = sizeof(T);
	size_t len = 1;
	auto min_stride = PyArray_STRIDE(vec, 0);
    for (int i=0;i<PyArray_NDIM(vec); ++i) {
		len *=  PyArray_DIM(vec, i);

		// Python/Numpy strides are in bytes, Blitz++ in sizeof(T) units.
		min_stride = std::min(min_stride, PyArray_STRIDE(vec, i) / T_size);
    }


	assert(T_size == PyArray_ITEMSIZE(vec));

    blitz::TinyVector<int,1> shape(0);
	shape[0] = len;
    blitz::TinyVector<int,1> strides(0);
	strides[0] = min_stride;

    return blitz::Array<T,1>((T*) PyArray_DATA(vec),shape,strides,
		blitz::neverDeleteData);
}
#endif

void GCMRegridder_add_sheet(GCMRegridder *cself,
	std::string name,
	std::string const &gridI_fname, std::string const &gridI_vname,
	std::string const &exgrid_fname, std::string const &exgrid_vname,
	std::string const &sinterp_style,
	PyObject *elevI_py, PyObject *maskI_py)
{
	NcIO ncio_I(gridI_fname, netCDF::NcFile::read);
	std::unique_ptr<Grid> gridI(new_grid(ncio_I, "grid"));
	gridI->ncio(ncio_I, "grid");
	ncio_I.close();

	NcIO ncio_exgrid(exgrid_fname, netCDF::NcFile::read);
	std::unique_ptr<Grid> exgrid(new_grid(ncio_exgrid, exgrid_vname));
	exgrid->ncio(ncio_exgrid, exgrid_vname);
	ncio_exgrid.close();

	auto interp_style(parse_enum<InterpStyle>(sinterp_style));
	auto elevI(np_to_blitz<double,1>(elevI_py, "elevI", {gridI->ndata()}));
	auto maskI(np_to_blitz<int,1>(maskI_py, "maskI", {gridI->ndata()}));

	SparseVector elevI_sp({elevI.extent(0)});
	for (int i=0; i<elevI.extent(0); ++i)
		if (!maskI(i)) elevI_sp.add({i}, elevI(i));

	auto sheet(new_ice_regridder(gridI->parameterization));
	sheet->init(name, std::move(gridI), std::move(exgrid),
		interp_style, std::move(elevI_sp));
	cself->add_sheet(std::move(sheet));
}

/** Computes yy = M xx.  yy is allocated, not necessarily set. */
void coo_matvec(PyObject *yy_py, PyObject *xx_py, bool ignore_nan,
	size_t M_nrow, size_t M_ncol, PyObject *M_row_py, PyObject *M_col_py, PyObject *M_data_py)
{
	auto xx(np_to_blitz<double,1>(xx_py, "xx", {M_ncol}));
	auto yy(np_to_blitz<double,1>(yy_py, "yy", {M_nrow}));
	auto M_row(np_to_blitz<int,1>(M_row_py, "M_row_py", {-1}));
	auto M_col(np_to_blitz<int,1>(M_col_py, "M_col_py", {-1}));
	auto M_data(np_to_blitz<double,1>(M_data_py, "M_data_py", {-1}));

	// Keep track of which items we've written to.
	std::vector<bool> written(M_nrow, false);

	// Do the multiplication, and we're done!
	int nnz = M_data.size();
	for (int n=0; n<nnz; ++n) {
		size_t row = M_row(n);
		size_t col = M_col(n);
		double data = M_data(n);

		// Ignore NaN in input vector
		if (ignore_nan && std::isnan(xx(col))) continue;

		// Just do Snowdrift-style "REPLACE".  "MERGE" was never used.
		double old_yy;
		if (written[row]) {
			old_yy = yy(row);
		} else {
			old_yy = 0;
			written[row] = true;
		}
		yy(row) = old_yy + data * xx(col);
	}

}


}}


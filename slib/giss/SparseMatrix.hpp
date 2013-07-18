#pragma once 

#include <map>
#include <boost/function.hpp>
#include <boost/bind.hpp>
#include <galahad/zd11_c.hpp>
#include "ncutil.hpp"
#include "IndexTranslator.hpp"
#include <giss/blitz.hpp>

class NcFile;

namespace giss {

// ------------------------------------------------------------

// ------------------------------------------------------------


// Matches up with Fortran's sparsecoord_t

/** Common header for sparse matrix types.
Describes the properties and storage format of the matrix.  Derived from the DESCRA array of the NIST SparseBlas proposal:
<pre>int descra[] Descriptor argument.  Nine element integer array
             descra[0] matrix structure
                     0 : general
                     1 : symmetric
                     2 : Hermitian
                     3 : Triangular
                     4 : Skew(Anti)-Symmetric
                     5 : Diagonal
             descra[1] upper/lower triangular indicator
                     1 : lower
                     2 : upper
             descra[2] main diagonal type
                     0 : non-unit
                     1 : unit
             descra[3] Array base
                     0 : C/C++ compatible
                     1 : Fortran compatible
             descra[4] repeated indices?
                     0 : unknown
                     1 : no repeated indices</pre>

@see http://www.nektar.info/browser/trunk/library/LibUtilities/LinearAlgebra/SparseBlas.hpp */
class SparseDescr {
public:

enum class MatrixStructure {GENERAL, SYMMETRIC, HERMETIAN, TRIANGULAR, ANTI_SYMMETRIC, DIAGONAL};
enum class TriangularType {GENERAL, LOWER, UPPER};
enum class MainDiagonalType {NON_UNIT, UNIT};
enum class DuplicatePolicy {REPLACE, ADD};
enum class SortOrder {ROW_MAJOR, COLUMN_MAJOR};

	MatrixStructure matrix_structure;
	TriangularType triangular_type;	
	MainDiagonalType main_diagonal_type;
	int index_base;

	// DuplicatePolicy const duplicate_policy;

	/** Number of rows in the matrix. */
	int const nrow;
	/** Number of columns in the matrix. */
	int const ncol;
	// int const nnz;

	SparseDescr &descr() { return *this; }

	/** Construct a sparse matrix header. */
	SparseDescr(
		int const _nrow, int const _ncol, int _index_base = 0,
		MatrixStructure _matrix_structure = MatrixStructure::GENERAL,
		TriangularType _triangular_type = TriangularType::GENERAL,
		MainDiagonalType _main_diagonal_type = MainDiagonalType::NON_UNIT) :

		nrow(_nrow),
		ncol(_ncol),
		index_base(_index_base),
		matrix_structure(_matrix_structure),
		triangular_type(_triangular_type),
		main_diagonal_type(_main_diagonal_type)
	{}
};

/** Base class for all sparse matrix types. */
class SparseMatrix : public SparseDescr {
public:
	SparseMatrix(SparseDescr const &descr) : SparseDescr(descr) {}

	virtual ~SparseMatrix() {}

	/** Returns the number of non-zero elements in the SparseMatrix. */
	virtual size_t size() const = 0;

	/** Removes all elements from the matrix. */
	virtual void clear() = 0;

	/** Set the value of an element in the matrix.
	Checks that row and column indices are within range.
	If the matrix is triangular, swaps row and col if needed, to get the element in the correct upper/lower triangle.
	@param row Row of element to set (base=0).
	@param col Column of element to set (base=0).
	@param val Value to set at (row, col).
	@param dups What should be done if there was an existing element in (row,col)?
	<dl><dt>dups=DuplicatePolicy::REPLACE</dt><dd> Replace with new value.</dd>
	<dt>dups=DuplicatePolicy::ADD</dt><dd> Add new value to existing value.</dd></dl>
	Note that only MapSparseMatrix supports these modes.  VectorSparseMatrix and ZD11SparseMatrix ignore the values of dups. */
	virtual void set(int row, int col, double const val, DuplicatePolicy dups = DuplicatePolicy::REPLACE) = 0;

	/** Convenience function, equal to set(row, col, val, DuplicatePolicy::ADD) */
	void add(int row, int col, double const val)
		{ set(row, col, val, DuplicatePolicy::ADD); }

	/** Used to write this data structure to a netCDF file.
	Defines the required variables.  Call the returned boost::function
	later to write the data.
	@param nc NetCDF file to write
	@param vname Variable name to use in writing this sparse matrix.
	@return Function to call later to write the data. */
	virtual boost::function<void ()> netcdf_define(NcFile &nc, std::string const &vname) const = 0;

	/** Multiply this matrix A by a vector.
	Computes y = A * x
	@param x IN: A vector of length ncol
	@param y OUT: A vector of length nrow
	@param clear_y If true, y = Ax.  Otherwise, y += Ax. */
	virtual void multiply(double const * x, double *y, bool clear_y = true) const = 0;

	/** Multiply transpose of this matrix A by a vector.
	Computes y = A^T * x
	@param x IN: A vector of length nrow
	@param y OUT: A vector of length ncol
	@param clear_y If true, y = A^T x.  Otherwise, y += A^T x. */
	virtual void multiplyT(double const * x, double *y, bool clear_y = true) const = 0;

	/** Computes the sum of each row of this matrix.
	@return Vector[nrow], each element containing the sum of the respective row from the matrix. */
	virtual std::vector<double> sum_per_row() const = 0;

	/** Computes the sum of each column of this matrix.
	@return Vector[ncol], each element containing the sum of the respective column from the matrix. */
	virtual std::vector<double> sum_per_col() const = 0;

	/** Computes the sum of each row of this matrix.
	@return A map, providing a value only for rows with at least one element. */
	virtual std::map<int,double> sum_per_row_map() const = 0;

	/** Computes the sum of each column of this matrix.
	@return A map, providing a value only for columns with at least one element. */
	virtual std::map<int,double> sum_per_col_map() const = 0;

};

// =================================================================
// Mix-ins common to all sparse matrix types

/** Mix-in common to all sparse matrix types.  Not part of the API. */
template<class SparseMatrix0T>
class SparseMatrix1 : public SparseMatrix0T
{
protected:
	SparseMatrix1(SparseDescr const &descr) :
		SparseMatrix0T(descr) {}
public:
	void set(int row, int col, double const val, SparseMatrix::DuplicatePolicy dups = SparseMatrix::DuplicatePolicy::REPLACE);
	boost::function<void ()> netcdf_define(NcFile &nc, std::string const &vname) const;

	void multiply(double const * x, double *y, bool clear_y = true) const;
	void multiplyT(double const * x, double *y, bool clear_y = true) const;
	std::vector<double> sum_per_row() const;
	std::vector<double> sum_per_col() const;
	std::map<int,double> sum_per_row_map() const;
	std::map<int,double> sum_per_col_map() const;

	template<class SparseMatrixT>
	void append(SparseMatrixT const &mat) {
		if (mat.nrow != this->nrow || mat.ncol != this->ncol) {
			fprintf(stderr, "SparseMatrix::append() has wrong size argument (%d, %d) vs. (%d, %d) expected\n", mat.nrow, mat.ncol, this->nrow, this->ncol);
			throw std::exception();
		}
		for (auto ii=mat.begin(); ii != mat.end(); ++ii)
			this->add(ii.row(), ii.col(), ii.val());
	}

private:
	void netcdf_write(NcFile *nc, std::string const &vname) const;
};


template<class SparseMatrix0T>
void SparseMatrix1<SparseMatrix0T>::set(int row, int col, double const val, SparseMatrix::DuplicatePolicy dups)
{
	// Fix row and col for possible triangular type
	switch(this->triangular_type) {
		case SparseMatrix::TriangularType::UPPER :
			if (row > col) std::swap(row, col);
		break;
		case SparseMatrix::TriangularType::LOWER :
			if (row < col) std::swap(row, col);
		break;
	}

	// Check range
	if (row >= this->nrow || row < 0) {
		fprintf(stderr, "SparseMatrix1<>::set(), row=%d >= nrow=%d or <0\n", row, this->nrow);
		throw std::exception();
	}
	if (col >= this->ncol || col < 0) {
		fprintf(stderr, "SparseMatrix1<>::set(), col=%d >= ncol=%d or <0\n", col, this->ncol);
		throw std::exception();
	}

	// Adjust for index_base
	row += this->index_base;
	col += this->index_base;

	this->_set(row, col, val, dups);
}

template<class SparseMatrix0T>
void SparseMatrix1<SparseMatrix0T>::netcdf_write(NcFile *nc, std::string const &vname) const
{
	NcVar *grid_indexVar = nc->get_var((vname + ".index").c_str());
	NcVar *areaVar = nc->get_var((vname + ".val").c_str());

	int i=0;
	for (auto ov = this->begin(); ov != this->end(); ++ov) {
		grid_indexVar->set_cur(i,0);
		int index[2] = {ov.row() + this->index_base, ov.col() + this->index_base};
		grid_indexVar->put(index, 1,2);

		areaVar->set_cur(i);
		areaVar->put(&ov.val(), 1);

		++i;
	}
}


template<class SparseMatrix0T>
boost::function<void ()> SparseMatrix1<SparseMatrix0T>::netcdf_define(
	NcFile &nc, std::string const &vname) const
{
	auto lenDim = nc.add_dim((vname + ".num_elements").c_str(), this->size());
	auto num_gridsDim = nc.add_dim((vname + ".rank").c_str(), 2);
	auto grid_indexVar = nc.add_var((vname + ".index").c_str(), ncInt, lenDim, num_gridsDim);
	auto areaVar = nc.add_var((vname + ".val").c_str(), ncDouble, lenDim);

	auto oneDim = get_or_add_dim(nc, "one", 1);
	auto descrVar = nc.add_var((vname + ".descr").c_str(), ncInt, oneDim);	// TODO: This should be ".info"
	descrVar->add_att("nrow", this->nrow);
	descrVar->add_att("ncol", this->ncol);
	descrVar->add_att("index_base", this->index_base);
	descrVar->add_att("matrix_structure", (int)this->matrix_structure);
	descrVar->add_att("triangular_type", (int)this->triangular_type);
	descrVar->add_att("main_diagonal_type", (int)this->main_diagonal_type);

	return boost::bind(&SparseMatrix1<SparseMatrix0T>::netcdf_write,
		this, &nc, vname);
}

/// Computes y = A * x
template<class SparseMatrix0T>
void SparseMatrix1<SparseMatrix0T>::multiply(double const * x, double *y, bool clear_y) const
{
	int nx = this->ncol;
	int ny = this->nrow;
	if (clear_y) for (int iy = 0; iy < ny; ++iy) y[iy] = 0;
	for (auto ii = this->begin(); ii != this->end(); ++ii) {
		int ix = ii.col();
		int iy = ii.row();
		y[iy] += ii.val() * x[ix];
	}
}

/// Computes y = A^T * x
template<class SparseMatrix0T>
void SparseMatrix1<SparseMatrix0T>::multiplyT(double const * x, double *y, bool clear_y) const
{
	int nx = this->nrow;
	int ny = this->ncol;
	if (clear_y) for (int iy = 0; iy < ny; ++iy) y[iy] = 0;
	for (auto ii = this->begin(); ii != this->end(); ++ii) {
		int iy = ii.col();
		int ix = ii.row();
		y[iy] += ii.val() * x[ix];
	}
}

// ------------------------------------------------------------
template<class SparseMatrix0T>
std::vector<double> SparseMatrix1<SparseMatrix0T>::sum_per_row() const {
	std::vector<double> ret(this->nrow);
	for (auto ii = this->begin(); ii != this->end(); ++ii) {
		ret[ii.row()] += ii.val();
	}
	return ret;
}

template<class SparseMatrix0T>
std::vector<double> SparseMatrix1<SparseMatrix0T>::sum_per_col() const {
	std::vector<double> ret(this->ncol);
	for (auto ii = this->begin(); ii != this->end(); ++ii) {
		ret[ii.col()] += ii.val();
	}
	return ret;
}


template<class SparseMatrix0T>
std::map<int,double> SparseMatrix1<SparseMatrix0T>::sum_per_row_map() const {
	std::map<int,double> ret;
	for (auto ii = this->begin(); ii != this->end(); ++ii) {
		auto f = ret.find(ii.row());
		if (f == ret.end()) {
			ret.insert(std::make_pair(ii.row(), ii.val()));
		} else {
			f->second += ii.val();
		}
	}
	return ret;
}

template<class SparseMatrix0T>
std::map<int,double> SparseMatrix1<SparseMatrix0T>::sum_per_col_map() const {
	std::map<int,double> ret;
	for (auto ii = this->begin(); ii != this->end(); ++ii) {
		auto f = ret.find(ii.col());
		if (f == ret.end()) {
			ret.insert(std::make_pair(ii.col(), ii.val()));
		} else {
			f->second += ii.val();
		}
	}
	return ret;
}

// ------------------------------------------------------------




// =================================================================
// Three different kinds of Sparse matrices

// ------------------------------------------------------------
/** Mix-in, not part of the API. */
class ZD11SparseMatrix0 : public SparseMatrix
{
protected:
	// Current number of elements in matrix.  Matrix is not
	// valid until this equals zd11.ne
	int _nnz_cur;

	// Pointers/references to main storage
	galahad::zd11_c *_zd11;

	ZD11SparseMatrix0(SparseDescr const &descr) : SparseMatrix(descr) {}
public:

	galahad::zd11_c &zd11() const { return *_zd11; }

	// --------------------------------------------------
	/** Standard STL-type iterator for iterating through a ZD11SparseMatrix. */
	class iterator {
	protected:
		ZD11SparseMatrix0 *parent;
		int i;
	public:
		iterator(ZD11SparseMatrix0 *z, int _i) : parent(z), i(_i) {}
		bool operator==(iterator const &rhs) { return i == rhs.i; }
		bool operator!=(iterator const &rhs) { return i != rhs.i; }
		void operator++() { ++i; }
		int row() { return parent->zd11().row[i] - parent->index_base; }
		int col() { return parent->zd11().col[i] - parent->index_base; }
		double &val() { return parent->zd11().val[i]; }
		double &value() { return val(); }
	};
	iterator begin() { return iterator(this, 0); }
	iterator end() { return iterator(this, _nnz_cur); }
	// --------------------------------------------------
	// --------------------------------------------------
	/** Standard STL-type const_iterator for iterating through a ZD11SparseMatrix. */
	class const_iterator {
	protected:
		ZD11SparseMatrix0 const *parent;
		int i;
	public:
		const_iterator(ZD11SparseMatrix0 const *z, int _i) : parent(z), i(_i) {}
		bool operator==(const_iterator const &rhs) { return i == rhs.i; }
		bool operator!=(const_iterator const &rhs) { return i != rhs.i; }
		void operator++() { ++i; }
		int row() { return parent->zd11().row[i] - parent->index_base; }
		int col() { return parent->zd11().col[i] - parent->index_base; }
		double const &val() { return parent->zd11().val[i]; }
		double const &value() { return val(); }
	};
	const_iterator begin() const { return const_iterator(this, 0); }
	const_iterator end() const { return const_iterator(this, _nnz_cur); }
	// --------------------------------------------------

	void clear() { _nnz_cur = 0; }

	bool is_complete() { return _nnz_cur == zd11().ne; }

	size_t size() const { return _nnz_cur; }

protected:
	/** Internal set() function without any error checking or adjustments for index_base. */
	void _set(int row, int col, double const val, DuplicatePolicy dups)
	{
		if (_nnz_cur >= zd11().ne) {
			fprintf(stderr, "ZD11SparseMatrix is full with %d elements\n", zd11().ne);
			throw std::exception();
		}
		zd11().row[_nnz_cur] = row;
		zd11().col[_nnz_cur] = col;
		zd11().val[_nnz_cur] = val;
		++_nnz_cur;
	}
};
// -----------------------------------------------------------------
/** A SparseMatrix based on external Fortran HSL_ZD11 storage.
<p>A standardized sparse matrix interface to the Fortran ZD11 data structure.
ZD11SparseMatrix wraps an existing ZD11 data structure --- which itself
is a peer to an underlying hsl_zd11_double::zd11_type Fortran derived type.</p>
<p>Use this class when you need to construct a ZD11-type matrix and pass
it to/from a Fortran subroutine.  <b>NOTE:</b> Not only the dimensions,
but also the number of non-zero elements must be pre-determined.  This
data structure is also used to copy from an existing SparseMatrix of another type.</p>
@see hsl_zd11_double::zd11_type, ZD11, hsl_zd11_double_x::hsl_zd11_c_init, copy */
class ZD11SparseMatrix : public SparseMatrix1<ZD11SparseMatrix0>
{
public:
	/** Call this after ZD11 has been initialized.
	@param zd11 C++ peer of Fortran hsl_zd11d::zd11_f sparse matrix structure.
	@param nnz_cur The number of non-zero elements currently held in _zd11.  Set to zero to clear the matrix.  <b>NOTE:</b> The sparse matrix must be filled with exactly ZD11::ne elements before it may be considered valid and passed to a Fortran subroutine.
 */
	ZD11SparseMatrix(galahad::zd11_c &zd11, int nnz_cur,
		MatrixStructure matrix_structure = MatrixStructure::GENERAL,
		TriangularType triangular_type = TriangularType::GENERAL,
		MainDiagonalType main_diagonal_type = MainDiagonalType::NON_UNIT)
	: SparseMatrix1<ZD11SparseMatrix0>(SparseDescr(zd11.m, zd11.n, 1,
	matrix_structure, triangular_type, main_diagonal_type))
	{
		_nnz_cur = nnz_cur;
		_zd11 = &zd11;
	}
};
// ==================================================================
// ---------------------------------------------------------
/** Mix-in, not part of the API. */
class VectorSparseMatrix0 : public SparseMatrix
{
friend class BlitzSparseMatrix;
protected:
	std::vector<int> indx;
	std::vector<int> jndx;
	std::vector<double> val;

	VectorSparseMatrix0(SparseDescr const &descr) : SparseMatrix(descr) {}

public:

	std::vector<int> const &rows() const { return indx; }
	std::vector<int> const &cols() const { return jndx; }
	std::vector<int> const &rowcols(int i) const
		{ return (i == 0 ? rows() : cols()); }
	std::vector<double> const &vals() const { return val; }

	// --------------------------------------------------
	/** Standard STL-type iterator for iterating through a ZD11SparseMatrix. */
	class iterator {
	protected:
		VectorSparseMatrix0 *parent;
		int i;
	public:
		int position() const { return i; }
		iterator(VectorSparseMatrix0 *p, int _i) : parent(p), i(_i) {}
		bool operator==(iterator const &rhs) const { return i == rhs.i; }
		bool operator!=(iterator const &rhs) const { return i != rhs.i; }
		void operator++() { ++i; }
		int row() const { return parent->indx[i] - parent->index_base; }
		int col() const { return parent->jndx[i] - parent->index_base; }
		int rowcol(int i) const
			{ return (i == 0 ? row() : col()); }
		double &val() { return parent->val[i]; }
		double &value() { return val(); }
	};
	iterator begin() { return iterator(this, 0); }
	iterator end() { return iterator(this, val.size()); }
	// --------------------------------------------------
	class const_iterator {
	protected:
		VectorSparseMatrix0 const *parent;
		int i;
	public:
		int position() const { return i; }
		const_iterator(VectorSparseMatrix0 const *p, int _i) : parent(p), i(_i) {}
		bool operator==(const_iterator const &rhs) const { return i == rhs.i; }
		bool operator!=(const_iterator const &rhs) const { return i != rhs.i; }
		void operator++() { ++i; }
		int row() const { return parent->indx[i] - parent->index_base; }
		int col() const { return parent->jndx[i] - parent->index_base; }
		int rowcol(int i) const
			{ return (i == 0 ? row() : col()); }
		double const &val() { return parent->val[i]; }
		double const &value() { return val(); }
	};
	const_iterator begin() const { return const_iterator(this, 0); }
	const_iterator end() const { return const_iterator(this, val.size()); }

	// --------------------------------------------------

	void clear() {
		indx.clear();
		jndx.clear();
		val.clear();
	}
	void reserve(size_t n) {
		indx.reserve(n);
		jndx.reserve(n);
		val.reserve(n);
	}
	size_t size() const { return val.size(); }

protected :
	/** Internal set() function without any error checking or adjustments for index_base. */
	void _set(int row, int col, double _val, DuplicatePolicy dups)
	{
		indx.push_back(row);
		jndx.push_back(col);
		val.push_back(_val);
	}
};

/** A SparseMatrix type made out of std::vector.  This is good when you don't
know how many elements you will have in the matrix in the end. */
class VectorSparseMatrix : public SparseMatrix1<VectorSparseMatrix0>
{
public:
	/** Construct a new sparse matrix to the given specifications.
	Memory will be allocated later as elements are added. */
	explicit VectorSparseMatrix(SparseDescr const &descr) :
	SparseMatrix1<VectorSparseMatrix0>(descr)
	{}


	/** Construct from existing vectors.
	(For example, when read from a netCDF file).
	@param descr Specifications for the matrix.
	@param _indx The row of each non-zero element (base via index_base).
	@param _jndx The column of each non-zero element (base via index_base).
	@param _val The value of each non-zero element. */
	VectorSparseMatrix(SparseDescr const &descr,
		std::vector<int> &&_indx,
		std::vector<int> &&_jndx,
		std::vector<double> &&_val) :
	SparseMatrix1<VectorSparseMatrix0>(descr)
	{
		indx = std::move(_indx);
		jndx = std::move(_jndx);
		val = std::move(_val);
	}

	/** Sorts the elements in the sparse matrix by index.
	@param sort_order Sort to ROW_MAJOR or COLUMN_MAJOR ordering. */
	void sort(SparseMatrix::SortOrder sort_order = SortOrder::ROW_MAJOR);

	/** Sums together items in the matrix with duplicate (row, col) */
	void sum_duplicates(
		SparseMatrix::SortOrder sort_order = SparseMatrix::SortOrder::ROW_MAJOR);

	/** Construct a VectorSparseMatrix based on arrays in a netCDF file.
	@param nc The netCDF file
	@param vname Name of the variable in the netCDF file.
	(will be base name of the arrays needed to store the sparse matrix). */
	static std::unique_ptr<VectorSparseMatrix> netcdf_read(NcFile &nc, std::string const &vname);

};
// ====================================================================
class BlitzSparseMatrix0 : public SparseMatrix
{
protected:
//	// Hold a pointer to previous SparseMatrix, in case we're wrapping one.
//	std::unique_ptr<VectorSparseMatrix> wrapped;

	// Current number of elements in matrix.  Matrix is not
	// valid until this equals zd11.ne
	int _nnz_cur;

	blitz::Array<int,1> indx;
	blitz::Array<int,1> jndx;
	blitz::Array<double,1> val;

	BlitzSparseMatrix0(SparseDescr const &descr) : SparseMatrix(descr) {}

public:
#if 0
	/** NOTE: The "copy" constructor below actually creates references
	to existing Blitz data. */
	BlitzSparseMatrix0(SparseDescr const &descr,
		blitz::Array<int,1> &_rows,
		blitz::Array<int,1> &_cols,
		blitz::Array<double,1> &_vals,
		int nnz_cur = -1)
	: SparseMatrix(descr), indx(_rows), jndx(_cols), val(_vals)
	{
		_nnz_cur = (nnz_cur < 0 ? vals().size() : nnz_cur);
	}
#endif

//	size_t size() const { return vals.extent(0); }

	blitz::Array<int,1> const &rows() const { return indx; }
	blitz::Array<int,1> const &cols() const { return jndx; }
	blitz::Array<double,1> const &vals() const { return val; }


	// --------------------------------------------------
	/** Standard STL-type iterator for iterating through a ZD11SparseMatrix. */
	class iterator {
	protected:
		BlitzSparseMatrix0 *parent;
		int i;
	public:
		int position() const { return i; }
		iterator(BlitzSparseMatrix0 *p, int _i) : parent(p), i(_i) {}
		bool operator==(iterator const &rhs) const { return i == rhs.i; }
		bool operator!=(iterator const &rhs) const { return i != rhs.i; }
		void operator++() { ++i; }
		int row() const { return parent->indx(i) - parent->index_base; }
		int col() const { return parent->jndx(i) - parent->index_base; }
		double &val() { return parent->val(i); }
		double &value() { return val(); }
	};
	iterator begin() { return iterator(this, 0); }
	iterator end() { return iterator(this, _nnz_cur); }
	// --------------------------------------------------
	class const_iterator {
	protected:
		BlitzSparseMatrix0 const *parent;
		int i;
	public:
		int position() const { return i; }
		const_iterator(BlitzSparseMatrix0 const *p, int _i) : parent(p), i(_i) {}
		bool operator==(const_iterator const &rhs) const { return i == rhs.i; }
		bool operator!=(const_iterator const &rhs) const { return i != rhs.i; }
		void operator++() { ++i; }
		int row() const { return parent->indx(i) - parent->index_base; }
		int col() const { return parent->jndx(i) - parent->index_base; }
		double const &val() { return parent->val(i); }
		double const &value() { return val(); }
	};
	const_iterator begin() const { return const_iterator(this, 0); }
	const_iterator end() const { return const_iterator(this, _nnz_cur); }

	// --------------------------------------------------

	void clear() { _nnz_cur = 0; }

	bool is_complete() { return _nnz_cur == vals().size(); }

	size_t size() const { return _nnz_cur; }


protected :
	/** Internal set() function without any error checking or adjustments for index_base. */
	void _set(int row, int col, double const val, DuplicatePolicy dups)
	{
		if (_nnz_cur >= vals().size()) {
			fprintf(stderr, "ZD11SparseMatrix is full with %d elements\n", vals().size());
			throw std::exception();
		}
		this->indx(_nnz_cur) = row;
		this->jndx(_nnz_cur) = col;
		this->val(_nnz_cur) = val;
		++_nnz_cur;
	}
};

/** A SparseMatrix type made out of std::vector.  This is good when you don't
know how many elements you will have in the matrix in the end. */
class BlitzSparseMatrix : public SparseMatrix1<BlitzSparseMatrix0>
{
public:
	/** Construct a new sparse matrix to the given specifications.
	Memory will be allocated later as elements are added. */
	explicit BlitzSparseMatrix(SparseDescr const &descr, int nnz_max) :
	SparseMatrix1<BlitzSparseMatrix0>(descr)
	{
		this->_nnz_cur = 0;
		this->indx.reference(blitz::Array<int,1>(nnz_max));
		this->jndx.reference(blitz::Array<int,1>(nnz_max));
		this->val.reference(blitz::Array<double,1>(nnz_max));
	}

	explicit BlitzSparseMatrix(VectorSparseMatrix &mat) :
	SparseMatrix1<BlitzSparseMatrix0>(mat.descr())
	{
		this->_nnz_cur = mat.size();
		this->indx.reference(vector_to_blitz(mat.indx));
		this->jndx.reference(vector_to_blitz(mat.jndx));
		this->val.reference(vector_to_blitz(mat.val));
	}

#if 0
	BlitzSparseMatrix(BlitzSparseMatrix &mat) :
	SparseMatrix1<BlitzSparseMatrix0>(mat.descr())
	{
		this->_nnz_cur = mat.size();
		this->indx.reference(mat.indx);
		this->jndx.reference(mat.jndx);
		this->val.reference(mat.val);
	}

	explicit BlitzSparseMatrix(std::unique_ptr<VectorSparseMatrix> &&mat) :
	SparseMatrix1<BlitzSparseMatrix0>(*mat)
	{
		this->wrapped = std::move(mat);
		this->_nnz_cur = mat->size();
		this->indx.reference(vector_to_blitz(mat->indx));
		this->jndx.reference(vector_to_blitz(mat->jndx));
		this->val.reference(vector_to_blitz(mat->val));
	}
#endif


	/** Construct from existing vectors.
	(For example, when read from a netCDF file).
	@param descr Specifications for the matrix.
	@param _indx The row of each non-zero element (base via index_base).
	@param _jndx The column of each non-zero element (base via index_base).
	@param _val The value of each non-zero element. */
	BlitzSparseMatrix(SparseDescr const &descr,
		blitz::Array<int,1> &_rows,
		blitz::Array<int,1> &_cols,
		blitz::Array<double,1> &_vals,
		int nnz_cur = -1)
	: SparseMatrix1<BlitzSparseMatrix0>(descr)
	{
		indx.reference(_rows);
		jndx.reference(_cols);
		val.reference(_vals);
		_nnz_cur = (nnz_cur < 0 ? val.extent(0) : nnz_cur);
//printf("Constructed nnz = %d (%d %d)\n", _nnz_cur, val.extent(0), _vals.extent(0));
	}


#if 0
//	/** Sorts the elements in the sparse matrix by index.
//	@param sort_order Sort to ROW_MAJOR or COLUMN_MAJOR ordering. */
//	void sort(SparseMatrix::SortOrder sort_order = SortOrder::ROW_MAJOR);

//	/** Sums together items in the matrix with duplicate (row, col) */
//	void sum_duplicates();

	/** Construct a BlitzSparseMatrix based on arrays in a netCDF file.
	@param nc The netCDF file
	@param vname Name of the variable in the netCDF file.
	(will be base name of the arrays needed to store the sparse matrix). */
	static std::unique_ptr<BlitzSparseMatrix> netcdf_read(NcFile &nc, std::string const &vname);
#endif

};
// ====================================================================



// ----------------------------------------------------------
/** Mix-in, not part of the API. */
class MapSparseMatrix0 : public SparseMatrix {
protected :
	std::map<std::pair<int,int>, double> _cells;
	typedef std::map<std::pair<int,int>, double>::iterator ParentIterator;
	typedef std::map<std::pair<int,int>, double>::const_iterator const_ParentIterator;
	MapSparseMatrix0(SparseDescr const &descr) : SparseMatrix(descr) {}

public:

	// --------------------------------------------------
	/** Standard STL-type iterator for iterating through a MapSparseMatrix. */
	class iterator {
		ParentIterator ii;
		MapSparseMatrix0 *parent;
	public:

		iterator(ParentIterator const &_i, MapSparseMatrix0 *p) : ii(_i), parent(p) {}
		bool operator==(iterator const &rhs) { return ii == rhs.ii; }
		bool operator!=(iterator const &rhs) { return ii != rhs.ii; }
		void operator++() { ++ii; }
		int row() { return ii->first.first - parent->index_base; }
		int col() { return ii->first.second - parent->index_base; }
		double &val() { return ii->second; }
		double &value() { return val(); }
	};
	iterator begin() { return iterator(_cells.begin(), this); }
	iterator end() { return iterator(_cells.end(), this); }
	// --------------------------------------------------
	/** Standard STL-type const_iterator for iterating through a MapSparseMatrix. */
	class const_iterator {
		const_ParentIterator ii;
		MapSparseMatrix0 const *parent;
	public:

		const_iterator(const_ParentIterator const &_i, MapSparseMatrix0 const *p) : ii(_i), parent(p) {}
		bool operator==(const_iterator const &rhs) { return ii == rhs.ii; }
		bool operator!=(const_iterator const &rhs) { return ii != rhs.ii; }
		void operator++() { ++ii; }
		int row() { return ii->first.first - parent->index_base; }
		int col() { return ii->first.second - parent->index_base; }
		double const &val() { return ii->second; }
		double const &value() { return val(); }
	};
	const_iterator begin() const { return const_iterator(_cells.begin(), this); }
	const_iterator end() const { return const_iterator(_cells.end(), this); }
	// --------------------------------------------------

	void clear() { _cells.clear(); }

	size_t size() const { return _cells.size(); }

protected :
	/** Internal set() function without any error checking or adjustments for index_base. */
	void _set(int row, int col, double const val, DuplicatePolicy dups)
	{
		// Could make this find-insert operation a bit more efficient
		// by only indexing into the std::map once...
		auto ii = _cells.find(std::make_pair(row, col));
		if (ii != _cells.end()) {
			if (dups == DuplicatePolicy::ADD) ii->second += val;
			else ii->second = val;
		} else {
			_cells.insert(std::make_pair(std::make_pair(row, col), val));
		}
	}
};
// ---------------------------------------------------------------
/** A SparseMatrix type made out of std::map.
This is good if you need random access to existing elements.
DuplicatePolicy::ADD and DuplicatePolicy::REPLACE are both respected. */
class MapSparseMatrix : public SparseMatrix1<MapSparseMatrix0>
{
public:
	/** Construct a new sparse matrix to the given specifications.
	Memory will be allocated later as elements are added. */
	MapSparseMatrix(SparseDescr const &descr) :
		SparseMatrix1<MapSparseMatrix0>(descr) {}
};
// ============================================================


// ------------------------------------------------------------
/** Copy a to b.  Does not clear b */
template<class SparseMatrixT1, class SparseMatrixT2>
void copy(SparseMatrixT1 &a, SparseMatrixT2 &b, SparseMatrix::DuplicatePolicy dups = SparseMatrix::DuplicatePolicy::REPLACE)
{
	b.clear();
	for (typename SparseMatrixT1::iterator ii = a.begin(); ii != a.end(); ++ii) {
		b.set(ii.row(), ii.col(), ii.val(), dups);
	}
}
// ------------------------------------------------------------
/** Copy a to b while transposing.  Does not clear b. */
template<class SparseMatrixT1, class SparseMatrixT2>
void transpose(SparseMatrixT1 &a, SparseMatrixT2 &b, SparseMatrix::DuplicatePolicy dups = SparseMatrix::DuplicatePolicy::REPLACE)
{
	b.clear();
	for (typename SparseMatrixT1::iterator ii = a.begin(); ii != a.end(); ++ii) {
		b.set(ii.col(), ii.row(), ii.val(), dups);
	}
}
// ------------------------------------------------------------
/** Converts from a used-set to a translator */
template<class SparseMatrixT1>
inline void make_used_translators(SparseMatrixT1 &a,
IndexTranslator &trans_row,
IndexTranslator &trans_col,
std::set<int> *_used_row = NULL,	// If not NULL, output to here.
std::set<int> *_used_col = NULL)	// If not NULL, output to here.

{
	// Figure out what is used
	std::set<int> used_row;
	std::set<int> used_col;
	for (typename SparseMatrixT1::iterator ii = a.begin(); ii != a.end(); ++ii) {
		used_row.insert(ii.row());
		used_col.insert(ii.col());
	}

	// Convert used sets to translators
	trans_row.init(a.nrow, used_row);
	trans_col.init(a.ncol, used_col);

	if (_used_row) *_used_row = std::move(used_row);
	if (_used_col) *_used_col = std::move(used_col);
}


template<class SparseMatrixT1>
inline void make_used_row_translator(SparseMatrixT1 &a,
IndexTranslator &trans_row,
std::set<int> *_used_row = NULL)	// If not NULL, output to here.
{
	// Figure out what is used
	std::set<int> used_row;
	for (typename SparseMatrixT1::iterator ii = a.begin(); ii != a.end(); ++ii) {
		used_row.insert(ii.row());
	}

	// Convert used sets to translators
	trans_row.init(a.nrow, used_row);
	if (_used_row) *_used_row = std::move(used_row);
}

template<class SparseMatrixT1>
inline void make_used_col_translator(SparseMatrixT1 &a,
IndexTranslator &trans_col,
std::set<int> *_used_col = NULL)	// If not NULL, output to here.
{
	// Figure out what is used
	std::set<int> used_col;
	for (typename SparseMatrixT1::iterator ii = a.begin(); ii != a.end(); ++ii) {
		used_col.insert(ii.col());
	}

	// Convert used sets to translators
	trans_col.init(a.ncol, used_col);
	if (_used_col) *_used_col = std::move(used_col);
}



template<class SparseMatrixT1, class SparseMatrixT2>
inline void translate_indices(SparseMatrixT1 &a, SparseMatrixT2 &b,
IndexTranslator const &trans_row,
IndexTranslator const &trans_col,
double *row_sum = NULL,		// Place to sum row and column values, if it exists
double *col_sum = NULL,
SparseMatrix::DuplicatePolicy dups = SparseMatrix::DuplicatePolicy::REPLACE)
{
	b.clear();
	if (row_sum) for (int i=0; i<b.nrow; ++i) row_sum[i] = 0;
	if (col_sum) for (int i=0; i<b.ncol; ++i) col_sum[i] = 0;
	for (typename SparseMatrixT1::iterator ii = a.begin(); ii != a.end(); ++ii) {
		int arow = ii.row();
		int brow = trans_row.a2b(arow);
		int acol = ii.col();
		int bcol = trans_col.a2b(acol);

		b.set(brow, bcol, ii.val(), dups);
		if (row_sum) row_sum[brow] += ii.val();
		if (col_sum) col_sum[bcol] += ii.val();
	}
}
// ----------------------------------------------------

// =============================================================

/** Computes M * diag and stores back in M
@param diag [ncol] Elements of the diagonal of the matrix
*/
template<class SparseMatrixT>
void multiply_bydiag(SparseMatrixT &mat,
blitz::Array<double,1> const &diag)
{
	int ndiag = diag.extent(0);
	if (ndiag != mat.ncol) {
		fprintf(stderr, "Matrix-diagonal multiply with mismatched dimensions %d vs %d", ndiag, mat.ncol);
		throw std::exception();
	}

	// Multiply by it
	for (auto ii = mat.begin(); ii != mat.end(); ++ii) {
		ii.val() *= diag(ii.col());
	}
}

/** Computes diag * M, stores back in M
@param diag [nrow] Elements of the diagonal of the matrix
*/
template<class SparseMatrixT>
void multiply_bydiag(
blitz::Array<double,1> const &diag,
SparseMatrixT &mat)
{
	int ndiag = diag.extent(0);
	if (ndiag != mat.nrow) {
		fprintf(stderr, "Matrix-diagonal multiply with mismatched dimensions %d vs %d", ndiag, mat.nrow);
		throw std::exception();
	}

	// Multiply by it
	for (auto ii = mat.begin(); ii != mat.end(); ++ii)
		ii.val() *= diag(ii.row());
}

// ===============================================================
// ======== Extra Functions

extern std::unique_ptr<VectorSparseMatrix> multiply_eigen_algorithm(VectorSparseMatrix &a, VectorSparseMatrix &b);
extern std::unique_ptr<VectorSparseMatrix> multiply_giss_algorithm(VectorSparseMatrix &a, VectorSparseMatrix &b);

inline std::unique_ptr<VectorSparseMatrix> multiply(VectorSparseMatrix &a, VectorSparseMatrix &b)
	{ return multiply_eigen_algorithm(a, b); }

extern std::vector<int> get_rowcol_beginnings(
	VectorSparseMatrix const &a,
	int const rowcol);

/** Only use on a matrix that's been sorted row-major */
inline std::vector<int> get_row_beginnings(VectorSparseMatrix const &a)
	{ return get_rowcol_beginnings(a, 0); }

/** Only use on a matrix that's been sorted column-major */
inline std::vector<int> get_col_beginnings(VectorSparseMatrix const &a)
	{ return get_rowcol_beginnings(a, 1); }


}	// namespace giss

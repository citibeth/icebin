#pragma once

#include <vector>
#include <algorithm>


namespace giss{

class SparseVector {
public:

	enum class DuplicatePolicy {REPLACE, ADD};

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
	virtual void set(int index, double const val) = 0;
	virtual void add(int index, double const val) { set(index, val); }

	/** Sums or otherwise consolidates the sparse representation so
	there is no more than one entry per index. */
	virtual void consolidate(DuplicatePolicy dups = DuplicatePolicy::ADD);

#if 0
	/** Used to write this data structure to a netCDF file.
	Defines the required variables.  Call the returned boost::function
	later to write the data.
	@param nc NetCDF file to write
	@param vname Variable name to use in writing this sparse matrix.
	@return Function to call later to write the data. */
	virtual boost::function<void ()> netcdf_define(NcFile &nc, std::string const &vname) const = 0;
#endif

};

// ================================================================
/**
Implements a "sparse vector", stored as an array of <index, value>
pairs.  This is analogous to a coordinate-style sparse matrix
representation (VectorSparseMatrix).  The SparseAccumulator class does
a similar job, but it is analogous to a lookup-based sparse matrix
representation (MapSparseMatrix).  This class is templated in order to
support multiple <index, value> types.
@see SparseAccumulator, VectorSparseMatrix, MapSparseMatrix */
template<class IndexT, class ValT>
class VectorSparseVector : public SparseVector
{
	std::vector<std::pair<IndexT, ValT>> vals;

public:

	auto begin() -> decltype(vals.begin()) { return vals.begin(); }
	auto end() -> decltype(vals.end()) { return vals.end(); }

	auto begin() const -> decltype(vals.begin()) { return vals.begin(); }
	auto end() const -> decltype(vals.end()) { return vals.end(); }

	/** Insert a new <index, value> pair to the end of the sparse
	vector list.  Does not detect or eliminate duplicates. */
	void set(IndexT const &index, ValT const &val)
		{ vals.push_back(std::make_pair(index, val)); }

	/** Sort the <index, value> pairs by index. */
	void sort_stable()
		{ std::stable_sort(vals.begin(), vals.end()); }

	/** Sums duplicate indices, resulting in a vector with no
	duplicates.  Also sorts in order to do this. */
	void sum_duplicates();

	/** Uses only the last-set value for a given index*/
	void remove_duplicates();

	void consolidate(DuplicatePolicy dups = DuplicatePolicy::ADD)
	{
		if (dups == DuplicatePolicy::ADD) sum_duplicates();
		else remove_duplicates();
	}

};


template<class IndexT, class ValT>
void VectorSparseVector<IndexT,ValT>::sum_duplicates()
{
	sort_stable();

	// Scan through, overwriting our array
	// New array will never be bigger than original
	int j=-1;		// Last-written item in output array
	double last_index = -1;	// Last index we saw in input array
	for (int i=0; i<this->size(); ++i) {
		int index = (*this)[i].first;
		ValT val = (*this)[i].second;

		if (index == last_index) {
			(*this)[j].second += val;
		} else {
			++j;
			(*this)[j].second = val;
			last_index = index;
		}
	}
	int n = j+1;	// Size of output array

	this->resize(n);
	this->shrink_to_fit();
}
// ---------------------------------------------------
template<class IndexT, class ValT>
void VectorSparseVector<IndexT,ValT>::remove_duplicates()
{
	sort_stable();

	// Scan through, overwriting our array
	// New array will never be bigger than original
	int j=-1;		// Last-written item in output array
	int last_index = -1;	// Last index we saw in input array
	for (int i=0; i<this->size(); ++i) {
		int index = (*this)[i].first;
		ValT val = (*this)[i].second;

		if (index == last_index) {
			(*this)[j].second == val;
		} else {
			++j;
			(*this)[j].second = val;
			last_index = index;
		}
	}
	int n = j+1;	// Size of output array

	this->resize(n);
	this->shrink_to_fit();
}


// ================================================================
/** Implements a "sparse vector", stored as an
std::unordered_map<index, value>.  This is analogous to a lookup-style
sparse matrix representation (MapSparseMatrix).  The CooVecotr class
does a similar job, but it is analogous to a coordinate-style sparse
matrix representation (VectorSparseMatrix).  This class is templated
in order to support multiple <index, value> types.
@see CooVector, VectorSparseMatrix, MapSparseMatrix */
template<class IndexT, class ValT>
class MapSparseVector : public SparseVector
{
	std::unordered_map<IndexT, ValT> vals;
public :
	auto begin() -> decltype(vals.begin()) { return vals.begin(); }
	auto end() -> decltype(vals.end()) { return vals.end(); }

	auto begin() const -> decltype(vals.begin()) { return vals.begin(); }
	auto end() const -> decltype(vals.end()) { return vals.end(); }

	/** Inserts a new <index, value> pair to sparse vector.  If that
	element was already non-zero, adds to it. */
	void add(IndexT index, ValT val) {
		auto ii(vals.find(index));
		if (ii == vals.end()) vals.insert(std::make_pair(index, val));
		else ii->second += val;
	}

	/** Computes the sum of two SparseAccumulator vectors.  Puts result in this. */
	void add(MapSparseVector const &b) {
		for (auto ii = b.begin(); ii != b.end(); ++ii) {
			add(ii->first, ii->second);
		}
	}


	/** Inserts a new <index, value> pair to sparse vector.  If that
	element was already non-zero, replaces it. */
	void set(IndexT index, ValT val) {
		auto ii(vals.find(index));
		if (ii == vals.end()) vals.insert(std::make_pair(index, val));
		else ii->second = val;
	}

#if 0
	/** Looks up a value in the sparse index by index.  Throws an
	exception if the index does not exist. */
	ValT &operator[](IndexT index) {
		auto ii(vals.find(index));
		if (ii == vals.end()) {
			std::cout << "SparseAccumulator[" << index << "] doesn't exist" << std::endl;
			throw std::exception();
		}
		return ii->second;
	}

	/** Looks up a value in the sparse index by index.  Throws an
	exception if the index does not exist. */
	ValT const &operator[](IndexT const &index) const {
 		auto ii(vals.find(index));
 		if (ii == vals.end()) throw std::exception();
		return ii->second;
	}
#endif

	void divide_by(MapSparseVector const &b) {
		for (auto ii=begin(); ii != end(); ++ii) {
			IndexT ix = ii->first;
			auto jj(b.find(ix));
			if (jj == b.end()) {
				// Not found, so there's a zero value there.
				// Go ahead, divide by zero!  It will create the appropriate
				// NaN value as if we were dividing by a dense vector.
				// This is probably an error, and it will give the user an indication
 				// of where to fix things.
				ii->second /= 0.;
			} else {
				ii->second /= jj->second;
			}
		}
	}

	void multiply_by(MapSparseVector const &b) {
		for (auto ii=begin(); ii != end(); ++ii) {
			IndexT ix = ii->first;
			auto jj(b.find(ix));
			if (jj == b.end()) {
				// Not found, so there's a zero value there.
				ii->second *= 0.;
			} else {
				ii->second *= jj->second;
			}
		}
	}

};

template<class SparseMatrixT>
inline void accum_per_row(SparseMatrixT const &mat,
	MapSparseVector<int,double> &accum)
{
	for (auto ii = mat.begin(); ii != mat.end(); ++ii)
		accum.add(ii.row(), ii.val());
}

template<class SparseMatrixT>
inline void accum_per_col(SparseMatrixT const &mat,
	MapSparseVector<int,double> &accum)
{
	for (auto ii = mat.begin(); ii != mat.end(); ++ii)
		accum.add(ii.col(), ii.val());
}




}

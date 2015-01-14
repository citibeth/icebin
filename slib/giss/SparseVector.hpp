#pragma once

#include <vector>
#include <algorithm>
#include <unordered_map>
#include <blitz/array.h>
#include <giss/sort.hpp>

namespace giss{

template<class IndexT, class ValT>
class SparseVector {
public:

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
	virtual void set(IndexT const &index, ValT const &val) = 0;
	virtual void add(IndexT const &index, ValT const &val) { set(index, val); }

	/** Sums or otherwise consolidates the sparse representation so
	there is no more than one entry per index. */
	virtual void consolidate(giss::DuplicatePolicy dups = giss::DuplicatePolicy::ADD) = 0;

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
class VectorSparseVector : public SparseVector<IndexT, ValT>
{
	std::vector<std::pair<IndexT, ValT>> vals;

public:
	VectorSparseVector() {}
	VectorSparseVector(VectorSparseVector &&b) : vals(std::move(b.vals)) {}
	void operator=(VectorSparseVector &&b) { vals = std::move(b.vals); }

#if 0
	// These don't make sense for VectorSparseVector!!!
	auto operator[](IndexT const &ix) -> decltype(vals[ix])
		{ return vals[ix]; }
	auto operator[](IndexT const &ix) const -> decltype(vals[ix])
		{ return vals[ix]; }
#endif


	void sort() { std::sort(vals.begin(), vals.end()); }
	size_t size() const { return vals.size(); }
	void clear() { vals.clear(); }

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

	void consolidate(DuplicatePolicy dups = DuplicatePolicy::ADD);

};

template<class SparseVectorT, class ValT>
void to_blitz(SparseVectorT const &mat, blitz::Array<ValT,1> &out)
{
	out = 0;
	for (auto ii = mat.begin(); ii != mat.end(); ++ii)
		out(ii->first) = ii->second;
}


template<class SparseVectorT, class IndexT, class ValT>
void to_parallel_arrays(SparseVectorT const &mat,
	std::vector<IndexT> &indices, std::vector<ValT> &vals)
{
	for (auto ii = mat.begin(); ii != mat.end(); ++ii) {
		indices.push_back(ii->first);
		vals.push_back(ii->second);
	}
}


// ---------------------------------------------------
template<class IndexT, class ValT>
void VectorSparseVector<IndexT,ValT>::consolidate(DuplicatePolicy dups)
{
	if (size() == 0) return;

	sort_stable();

	// Scan through, overwriting our array
	// New array will never be bigger than original
	int j=0;		// Last-written item in output array
	IndexT last_index = vals[0].first;	// Last index we saw in input array
	for (int i=1; i<this->size(); ++i) {
		IndexT &index = vals[i].first;
		ValT &val = vals[i].second;

		if (index == last_index) {
			if (dups == DuplicatePolicy::ADD) {
				vals[j].second += val;
			} else {
				vals[j].second = val;
			}
		} else {
			++j;
			vals[j].first = index;
			vals[j].second = val;
			last_index = index;
		}
	}
	int n = j+1;	// Size of output array

	this->vals.resize(n);
	this->vals.shrink_to_fit();

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
class MapSparseVector : public SparseVector<IndexT, ValT>
{
	std::unordered_map<IndexT, ValT> vals;
public :
	auto find(IndexT const &ix) -> decltype(vals.find(ix))
		{ return vals.find(ix); }
	auto find(IndexT const &ix) const -> decltype(vals.find(ix))
		{ return vals.find(ix); }
	auto operator[](IndexT const &ix) -> decltype(vals[ix])
		{ return vals[ix]; }

	size_t size() const { return vals.size(); }
	void clear() { vals.clear(); }

	auto begin() -> decltype(vals.begin()) { return vals.begin(); }
	auto end() -> decltype(vals.end()) { return vals.end(); }

	auto begin() const -> decltype(vals.begin()) { return vals.begin(); }
	auto end() const -> decltype(vals.end()) { return vals.end(); }

	/** Inserts a new <index, value> pair to sparse vector.  If that
	element was already non-zero, adds to it. */
	void add(IndexT const &index, ValT const &val) {
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

	/** Nothing to do for consolidate with MapSparseVector. */
	void consolidate(DuplicatePolicy dups = DuplicatePolicy::ADD) {}

	/** Inserts a new <index, value> pair to sparse vector.  If that
	element was already non-zero, replaces it. */
	void set(IndexT const &index, ValT const &val) {
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


#if 0
	/** Element-by-element division. */
	void divide_by(MapSparseVector const &b) {
		for (auto ii=begin(); ii != end(); ++ii) {
			IndexT ix = ii->first;
			auto jj(b.vals.find(ix));
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
#endif

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

/** Element-by-element division. */
template<class SparseVectorT, class IndexT, class ValT>
void divide_by(SparseVectorT &vec, MapSparseVector<IndexT, ValT> const &b)
{
	for (auto ii=vec.begin(); ii != vec.end(); ++ii) {
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

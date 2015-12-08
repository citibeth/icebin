template<int RANK>
class ArrayDescr
{
public:
	enum class DuplicatePolicy {LEAVE_ALONE, ADD, REPLACE};


	std::array<size_t, RANK> const shape;		// Extent of each dimension
	std::array<std::string, RANK> const dims;				// Name of each dimension

	ArrayDescr(
		std::array<size_t, RANK> const _shape,
		std::array<std::string, RANK> const _dims)
	: shape(_shape), dims(_dims) {}

	ArrayDescr(std::array<size_t, RANK> const _shape)
	: shape(_shape) {}
}


// Values for sort_order formal parameter below
const std::array<int,2> ROW_MAJOR = {{0,1}};
const std::array<int,2> COLUMN_MAJOR = {{1,0}};

template<class IndexT, int RANK>
struct CmpIndex {
	std::array<std::vector<IndexT> *, RANK> const indices;

	CmpIndex(
		std::array<std::vector<IndexT>, RANK> const *_indices,
		std::vector<int,RANK> const &sort_order)
	{
		for (int i=0; i<RANK; ++i)
			indices[i] = &_indices[sort_order[i]];
	}

	bool operator()(int i, int j)
	{
		for (int k=0; k<RANK-1; ++k) {
			if (indices[k][i] < indices[k][j]) return true;
			if (indices[k][i] > indices[k][j]) return false;
		}
		return (indices[RANK-1][i] < indices[RANK-1][j]);
	}


template<class IndexT, class ValueT, int RANK>
class CooArray : public ArrayDescr
{
	typedef CooArray<class IndexT, class ValueT, int RANK> ThisCooArray;
public:
	std::array<std::vector<IndexT>, RANK> indices;
	std::vector<ValueT> vals;

private:
	bool in_edit;		// Are we in edit mode?
	const std::array<int,RANK> sort_order;	// Non-negative elements if this is sorted

public:
	CooArray() : in_edit(true) {}
	int rank() { return RANK; }

	// Move semantics
	CooArray(CooArray &&other) :
		indices(std::move(other.indices)),
		vals(std::move(other.vals)),
		in_edit(other.in_edit),
		sort_order(other.sort_order) {}

	void operator=(ThisCooArray &&other) {
		indices = std::move(other.indices);
		vals = std::move(other.vals);
		in_edit = other.in_edit;
		sort_order = other.sort_order;
	}

	// Copy semantics
	CooArray(CooArray &other) :
		indices(other.indices),
		vals(other.vals),
		in_edit(other.in_edit),
		sort_order(other.sort_order) {}

	void operator=(ThisCooArray &other) {
		indices = other.indices;
		vals = other.vals;
		in_edit = other.in_edit;
		sort_order = other.sort_order;
	}


	// -------------------------------------------------
	size_t size()
		{ return vals.size(); }

	void clear() {
		for (int k=0; k<RANK; ++k) indicies[k].clear();
		vals.clear();
	}

	void reserve(size_t size) {
		vals.reserve(size);
		for (int k=0; k<RANK; ++k) indicies[k].reserve(size);
	}

	/** Adds an element to the sparse array. */
	void add(std::array<IndexT,RANK> const index, ValueT const &val)
	{
		for (int k=0; k<RANK; ++k) indicies[k].push_back(index[k]);
		vals.push_back(val);
	}
	void add(std::vector<IndexT> const &index, ValueT const &val)
	{
		for (int k=0; k<RANK; ++k) indicies[k].push_back(index[k]);
		vals.push_back(val);
	}

	// -------------------------------------------------
	// --------------------------------------------------
	/** Standard STL-type iterator for iterating through a VectorSparseMatrix. */
	class iterator {
	protected:
		ThisCooArray *parent;
		size_t i;
	public:
		iterator(ThisCooArray *p, size_t _i) : parent(p), i(_i) {}

		size_t position() const { return i; }
		bool operator==(iterator const &rhs) const { return i == rhs.i; }
		bool operator!=(iterator const &rhs) const { return i != rhs.i; }
		bool operator<(iterator const &rhs) const { return i < rhs.i; }
		void operator++() { ++i; }
		ValueT &val() { return parent->vals[i]; }

		IndexT &index(int k) const { return parent->indices[k][i]; }
		std::array<IndexT, RANK> index() const {
			std::array<size_t, RANK> ret;
			for (int k=0; k<RANK; ++k) ret[k] = index(k);
			return ret;
		}
		std::vector<IndexT> index_vec() const {
			std::vector<size_t> ret;
			for (int k=0; k<RANK; ++k) ret.push_back(index(k));
			return ret;
		}

		// Convenience methods
		IndexT &row() const { return index(0); }
		IndexT &col() const { return index(1); }
		ValueT &value() { return val(); }
	};
	iterator begin() { return iterator(this, 0); }
	iterator end() { return iterator(this, vals.size()); }
	// --------------------------------------------------
	class const_iterator {
	protected:
		ThisCooArray const *parent;
		size_t i;
	public:
		const_iterator(ThisCooArray *p, size_t _i) : parent(p), i(_i) {}
		size_t position() const { return i; }
		bool operator==(iterator const &rhs) const { return i == rhs.i; }
		bool operator!=(iterator const &rhs) const { return i != rhs.i; }
		bool operator<(iterator const &rhs) const { return i < rhs.i; }
		void operator++() { ++i; }
		ValueT &val() { return parent->vals[i]; }

		IndexT const &index(int k) const { return parent->indices[k][i]; }
		std::array<IndexT, RANK> index() const {
			std::array<size_t, RANK> ret;
			for (int k=0; k<RANK; ++k) ret[k] = index(k);
			return ret;
		}
		std::vector<IndexT> index_vec() const {
			std::vector<size_t> ret;
			for (int k=0; k<RANK; ++k) ret.push_back(index(k));
			return ret;
		}


		// Convenience methods
		IndexT const &row() const { return index(0); }
		IndexT const &col() const { return index(1); }
		ValueT const &value() { return val(); }
	};
	iterator begin() const { return const_iterator(this, 0); }
	iterator end() const { return const_iterator(this, vals.size()); }
	// --------------------------------------------------

	// -------------------------------------------------
	/** Goes in to add mode: legal to add more things to the vector. */
	void edit()
	{
		in_edit = true;
		sort_order = {{-1, -1}};
	}

	void add(std::array<IndexT, RANK> const index, ValueT const val)
	{
		if (in_edit) {
			std::cerr << "Must be in edit mode to use CooArray::add()" << std::endl;
			giss::exit(-1);
		}
		for (int i=0; i<RANK; ++i) indices[i].push_back(index[i]);
		vals.push_back(val);
	}

	void sort(std::array<int, 2> const _sort_order)
	{
		// Decide on how we'll sort
		CmpIndex cmp(&indices, _sort_order);

		// Generate a permuatation
		int n = size();
		std::vector<int> perm; perm.reserve(n);
		for (int i=0; i<n; ++i) perm.push_back(i);
		std::sort(perm.begin(), perm.end(), cmp);

		// Apply permutation to vals
		std::vector<ValueT> dtmp; dtmp.reserve(n);
		for (int i=0; i<n; ++i) dtmp.push_back(vals[perm[i]]);
		vals = std::move(dtmp);

		// Apply permutation to indices
		std::vector<IndexT> itmp; itmp.reserve(n);
		for (int k=0; k<RANK; ++k) {
			itmp.clear();
			for (int i=0; i<n; ++i) itmp.push_back(indices[k][perm[i]]);
			std::swap(itmp, indices[k]);
		}

		in_edit = false;
		sort_order = _sort_order;
	}

	void consolidate(std::array<int, RANK> sort_order, DuplicatePolicy duplicate_policy, handle_nan = false)
	{
		// Nothing to do for zero-size matrices (or ones with just one element)
		if (size() <= 1) return;

		// Decide on how we'll sort
		CmpIndex cmp(&indices, _sort_order);

		// Generate a permuatation
		int n = size();
		std::vector<int> perm; perm.reserve(n);
		for (int i=0; i<n; ++i) perm.push_back(i);
		std::stable_sort(perm.begin(), perm.end(), cmp);

		// Output arrays
		std::array<std::vector<IndexT>, RANK> nindices;
		std::vector<ValueT> nval;

		// Identify duplicates
		i = 0;
		while (vals[perm[i]] == 0 || std::isnan(vals[perm[i]])) ++i;		// Remove nans
		for (int k=0; k<RANK; ++k) nindices[k].push_back(indices[k][perm[i]]);
		nvals.push_back(vals[perm[i]]);
		for (int i=1; i<indicies.size(); ++i) {
			while (vals[perm[i]] == 0 || std::isnan(vals[perm[i]])) ++i;	// Remove nans
			std::array<IndexT, RANK> index;
			bool match_back = true;
			for (int k=0; k<RANK; ++k) {
				int ix = indicies[k][perm[i]];
				if (ix != nindices[k].back()) {
					match_back = false;
					break;
				}
			}
			if (match_back) {
				ValueT val = vals[perm[i]];
				if ((!handle_nan || !std::isnan(val)) && val != 0) {
					if (duplicate_policy == DuplicatePolicy::ADD)
						nvals.back() += val;
					else if (duplicate_policy == DuplicatePolicy::REPLACE)
						nvals.back() = vals[perm[i]];
				}
			} else {
				for (int k=0; k<RANK; ++k) nindices[k].push_back(indices[k][perm[i]]);
				ValueT val = vals[perm[i]];
				if ((!handle_nan || !std::isnan(val)) && val != 0) {
					nvals.push_back(vals[perm[i]]);
				}
			}

			// Change over to our new copy (double buffering here)
			for (int k=0; k<RANK; ++k) indicies[k] = std::move(nindicies[k])
			vals = std::move(nindices[vals]);
		}

		in_edit = false;
		this->sort_order = sort_order;
	}


	 void netcdf_define(
		netCDF::NcFile &nc, std::string const &vname,
		std::vector<std::function<void ()>> &writes) const
	{
		NcDim size_d = nc.addDim(vname + ".size", this->size());
		NcDim rank_d = nc.addDim(vname + ".rank", RANK);
		nc.add_var(vname + ".indices", ncInt, {size_d, rank_d});
		nc.add_var(vname + ".vals", ncDouble, {size_d});

		one_d = getOrAddDim(nc, "one", 1);
		auto descrVar = nc.add_var(vname + ".descr", ncInt, {one_d});	// TODO: This should be ".info"
		descrVar.putAtt("shape", ncLong, RANK, &shape);

		writes.push_back(&CooArray<IndexT,ValueT,RANK>::netcdf_write,
			this, &nc, vname);
	}


	void netcdf_write(NcFile &nc, std::string const &vname) const
	{
		NcVar indices_v = nc.getVar(vname + ".index");
		NcVar vals_v = nc.getVar(vname + ".val");

		vals_v.putVar(vals, size());

		std::vector<size_t> startp(2) = {0, 0};		// SIZE, RANK
		std::vector<size_t> countp(2) = {1, RANK};	// Write RANK elements at a time
		for (auto ov = this->begin(); ov != this->end(); ++ov) {
			std::array<size_t> index = index();

			indices_v.putVar(startp, countp, &index[0]);	// 2-D
			vals_v.putVar(startp, countp, &ov.val());	// 1-D

			++startp[0];
		}
	}
};


template<class IndexT, class ValueT>
using CooMatrix = CooArray<IndexT, ValueT, 2>

template<class IndexT, class ValueT>
using CooVector = CooArray<IndexT, ValueT, 1>

// ========================================================

/** (Sparse Matrix) * (Dense Vector) */
template<class CooMatrixT, class AccumulatorT>
void multiply(
	CooMatrixT const &M,
	blitz::Array<double,1> const &x,
	AccumulatorT &y,
	bool handle_nan = false,
	bool transpose = false)
{
	for (auto ii = M.begin(); ii != M.end(); ++ii) {
		std::array<IndexT,2> index(ii.index());
		if (transpose) std::swap(index[0], index[1]);
		double val = ii.val() * x(index[1]);
		if (!handle_nan || !(std::isnan(val) || std::isinf(val))) {
			y.add(  {{index[0]}},  val);
		}
	}
}

template<class CooMatrixT, class AccumulatorT>
void multiplyT(
	CooMatrixT const &M,
	blitz::Array<double,1> const &x,
	AccumulatorT &y,
	bool handle_nan = false)
{ return multiply(M, x, y, handle_nan, true); }



/** Computes [row] diag(scale) [[col]].
Or: row_j scale_j col_j   (out of overall matrix computation row_ij scale_j col_jk

NOTE: row, col and scale must NOT have duplicate elements.  Consolidate() must have been run before --- which gets ride of nans too, via handle_nan.
*/
static void multiply_row_col(
	AccumulatorT &ret,

	IndexT[] row_cols,		// Column for each element in the row
	ValueT[] row_vals,		// Value of each element in the row
	size_t row_size,		// Number of non-zero elements in the row

	bool has_scale,
	IndexT[] scale_indices,
	ValueT[] scale_vals,
	size_t scale_size,

	IndexT[] col_rows,		// Row for each element in the column
	ValueT[] col_vals,		// Value of each element in the column
	size_t col_size)		// Number of non-zero elements in the column
{
	if (row_size == 0 || col_size == 0) return 0.0;

	size_t ri = 0;
	size_t ci = 0;
	size_t si = 0;
	// Invariant: next_match = max(row_cols[ri], col_rows[ci], scale_indices[si])
	IndexT next_match = std::max(row_cols[ri], col_rows[ci])
	if (has_scale) next_match = std::max(next_match, scale_indices[si]);
	for (;; ++ri, ++ci, ++si) {

		// Scan forward in the row and column, looking for a match.
		for (;;++ri) {
			if (ri >= row_size) return ret;
			if (row_cols[ri] >= next_match) break;
		}
		next_match = row_cols[ri];	// >= the old next_match

		for (;;++ci) {
			if (ci >= col_size) return ret;
			if (col_rows[ci] >= next_match) break;
		}
		next_match = col_rows[ci];	// >= the old next_match

		if (has_scale) for (;;++si) {
			if (si >= scale_size) return ret;
			if (scale_indices[si] >= next_match) break;
		}
		next_match = scale_indices[si];	// >= the old next_match

		IndexT const rcol = row_cols[ri];
		IndexT const crow = col_rows[ci];
		if (rcol == crow &&
			((!scale_vals) || scale_indices[si] == rcol))
		{
			// Get current value of row, column and scale
			ValueT rval = row_vals[ri];
			ValueT cval = col_vals[ci];
			ValueT sval = (has_scale ? scale_vals[si] : 1.0);

			ValueT term = rval * sval * cval;
			ret.add({{rcol}}, term);
		}
	}
}

/** A must be sorted properly. */
std::vector<int> get_rowcol_beginnings(
	CooMatrix const &A,
	int const select_index)	// 0 for beginning of rows, 1 for beginning of cols
{
	std::vector<int> abegin;

	// Get beginning of each row in a (including sentinel at end)
	int last_row = -1;
	for (auto ai(A.begin()); ; ++ai) {
		if (ai == A.end()) {
			abegin.push_back(ai.position());
			break;
		}
		if (ai.index(select_index) != last_row) {
			abegin.push_back(ai.position());
			last_row = ai.index(select_index);
		}
	}

	return abegin;
}


/** NOTES:
A must be consolidated ROW_MAJOR, B must be consolidated COLUMN_MAJOR
*/
template<class IndexT, class ValueT, class AccumulatorT>
static double multiply(
	AccumulatorT &ret,
	CooVector<IndexT, ValueT> const *scalei,
	CooMatrix<IndexT, ValueT> const &A,
	CooVector<IndexT, ValueT> const *scalej,
	CooMatrix<IndexT, ValueT> const &B,
	CooVector<IndexT, ValueT> const *scalek)
{
	std::vector<size_t> abegin(get_rowcol_beginnings(A, 0));
	std::vector<size_t> bbegin(get_rowcol_beginnings(B, 1));

	// Multiply each row by each column

	// ---------- Outer loop: rows in A
	size_t ai = 0;
	size_t sii = 0;		// Index into sparse scalei representation
	IndexT next_match_a = A.indices[0][ai];
	if (scalei) next_match_a = std::max(next_match_a, scalei->indices[0][sii]);
	for (;; ++ai, ++sii) {
		// ---- Increment forward to the next row in A w/ matching scalei
		if (!scalei) {
			if (ai >= abegin.size()-1) goto break_outer;
		} else {
			for (;;++ai) {
				if (ai >= abegin.size()-1) goto break_outer;
				if (A.indices[0][ai] >= netxt_match_a) break;
			}
			next_match_a = A.indices[0][ai];

			for (;;++sii) {
				if (si >= scalei->size()) goto break_outer;
				if (scalei->indices[0][sii] >= next_match_a) break;
			}
			next_match_a = scalei->indices[0][sii];
		}

		ValueT sival;
		if (!scalei) {
			sival = 1.0;
		} else {
			sival = scalei->indices[0][sii];
			if (sival == 0) continue;
		}

		// -----------------------------------------------------
		// ------------ BEGIN Inner loop: columns in B
		size_t bi = 0;
		size_t ski = 0;		// Index into sparse scalei representation
		IndexT next_match_b = B.indices[1][bi];
		if (scalek) next_match_b = std::max(next_match_b, scalek->indices[0][ski]);
		for (;;++bi, ++ski) {
			// ---- Increment forward to the next row in A w/ matching scalek
			if (!scalek) {
				if (bi >= bbegin.size()-1) goto break_inner;
			} else {
				for (;;++bi) {
					if (bi >= bbegin.size()-1) goto break_inner;
					if (B.indices[1][bi] >= next_match_b) break;
				}
				next_match_b = B.indices[1][bi];

				for (;;++ski) {
					if (ski >= scalek->size()) goto break_inner;
					if (scalek->indices[0][ski] >= next_match_b) break;
				}
				next_match_b = scalek->indices[0][ski];
			}

			ValueT skval;
			if (!scalei) {
				skval = 1.0;
			} else {
				skval = scalek->indices[0][ski];
				if (skval == 0) continue;
			}

			// Multiply a row by a column
			IndexT const a0 = abegin[ai];
			IndexT const b0 = bbegin[bi];
			ScalarAccumulator<IndexT, ValueT, RANK> rcval;
			multiply_row_col(rcval,
				&A.indices[0][a0], &A.vals[a0], abegin[ai+1]-a0,
				scalej,
					scalej ? &scalej->indices[0][0] : 0,
					scalej ? &scalej->vals[0] : 0,
					scalej ? scalej->size() : 0,
				&B.indices[1][b0], &B.vals[b0], bbegin[bi+1]-b0)

			double val = sival * rcval.value() * skval;
			if (val == 0) continue;

			IndexT aRow = a.indices[0][a0];
			IndexT bCol = b.indices[1][b0];
			ret.add({{aRow, bCol}}, val);
		}
		break_inner: ;
		// ------------ END Inner loop: columns in B
	}
	break_outer: ;

}





/** NOTES:
A must be consolidated ROW_MAJOR, B must be consolidated COLUMN_MAJOR
*/
template<class IndexT, class ValueT, class AccumulatorT>
static double multiply(
	AccumulatorT &ret,
	CooVector<IndexT, ValueT> const *scalei,
	CooMatrix<IndexT, ValueT> const &A,
	CooVector<IndexT, ValueT> const &b)
{
	std::vector<size_t> abegin(get_rowcol_beginnings(A, 0));
	std::vector<size_t> bbegin(get_rowcol_beginnings(B, 1));

	// Multiply each row by each column

	// ---------- Outer loop: rows in A
	size_t ai = 0;
	size_t sii = 0;		// Index into sparse scalei representation
	IndexT next_match_a = A.indices[0][ai];
	if (scalei) next_match_a = std::max(next_match_a, scalei->indices[0][sii]);
	for (;; ++ai, ++sii) {
		// ---- Increment forward to the next row in A w/ matching scalei
		if (!scalei) {
			if (ai >= abegin.size()-1) goto break_outer;
		} else {
			for (;;++ai) {
				if (ai >= abegin.size()-1) goto break_outer;
				if (A.indices[0][ai] >= netxt_match_a) break;
			}
			next_match_a = A.indices[0][ai];

			for (;;++sii) {
				if (si >= scalei->size()) goto break_outer;
				if (scalei->indices[0][sii] >= next_match_a) break;
			}
			next_match_a = scalei->indices[0][sii];
		}

		ValueT sival;
		if (!scalei) {
			sival = 1.0;
		} else {
			sival = scalei->indices[0][sii];
			if (sival == 0) continue;
		}

		// ------------ Just one column of b

		// Multiply a row by a column
		IndexT const a0 = abegin[ai];
		ScalarAccumulator<IndexT, ValueT, RANK> rcval;
		multiply_row_col(rcval,
			&A.indices[0][a0], &A.vals[a0], abegin[ai+1]-a0,
			0, 0, 0, 0,
			&b.indices[0][0], &b.vals[0], b.size())
			handle_nan);

		double val = sival * rcval.value();
		if (val == 0) continue;

		IndexT aRow = a.indices[0][a0];
		ret.add({{aRow}}, val);

	}
	break_outer: ;

}






template<class IndexT, class ValueT, class AccumulatorT>
static double multiply(
	AccumulatorT &ret,
	CooVector<IndexT, ValueT> const &A,
	CooVector<IndexT, ValueT> const &B,
	bool handle_nan = false)
{
	multiply_row_col(ret,
		&a.indices[0][0], &a.vals[0], a.size(),
		0, 0, 0, 0,
		&b.indices[0][0], &b.vals[0], b.size())
		handle_nan);
}



/** Copy a to b while transposing.  Does not clear b. */
template<class CooMatrixT>
CooMatrixT transpose(
	AccumulatorT &ret,
	CooMatrixT &A)
{
	CooMatrixT ret;
	for (auto ii = A.begin(); ii != A.end(); ++ii) {
		ret.add(ii.col(), ii.row(), ii.val(), dups);
	}
	return ret;
}

template<class AccumulatorT, class IndexT, class ValueT, int RANK>
void invert(
	AccumulatorT &ret,
	CooArray<IndexT, ValueT, RANK> &A)
{
	size_t size = A.size();
	ret.reserve(size);
	for (auto ii = A.begin(); ii != A.end(); ++ii) {
		std::array<ValueT, RANK> index;
		for (int k=0; i<RANK; ++k) index[k] = A.indices[k][i];
		ret.add(index, A.vals[i]);
	}
}


template<class AccumulatorT, class IndexT, class ValueT, int RANK>
void invert(
	AccumulatorT &ret,
	CooArray<IndexT, ValueT, RANK> &A)


// --------------------------------------------------------------
template<int IndexT, class ValueT, int RANK>
class BlankAccumulator {
public:
	int rank() { return 0; }
	void reserve(size_t size) {}
	void add(std::array<IndexT, RANK> const index, ValueT const &val) {}
	void add(std::vector<IndexT> const &index, ValueT const &val) {}
};


/** Accumulates into a dense array */
template<class IndexT, class ValueT, int RANK, int DUPLICATE_POLICY>
class DenseAccumulator : public BlankAccumulator {
	blitz::Array<ValueT, RANK> &dense;
	DuplicatePolicy duplicate_policy;
public:
	int rank() { return RANK; }
	DenseAccumulator(blitz::Array<ValueT, RANK> &_dense, DuplicatePolicy _duplicate_policy)
		: dense(_dense), duplicate_policy(_duplicate_policy) {}


	void add(TinyVector<int,RANK> &index, ValueT const &val)
	{
		ValueT &oval(dense(bindex));

		switch(DUPLICATE_POLICY) {
			case DuplicatePolicy::LEAVE_ALONE :
				if (!std::isnan(oval)) oval = val;
			break;
			case DuplicatePolicy::ADD :
				oval += val;
			break;
			case DuplicatePolicy::REPLACE :
				oval = val;
			break;
		}
	}

	void add(std::array<IndexT, RANK> const index, ValueT const &val)
	{
		TinyVector<int,RANK> bindex;
		for (int k=0; k<RANK; ++k) bindex(k) = index[k];
		add(bindex, val);
	}

	void add(std::vector<IndexT> const index, ValueT const &val)
	{
		TinyVector<int,RANK> bindex;
		for (int k=0; k<RANK; ++k) bindex(k) = index[k];
		add(bindex, val);
	}


	blitz::Array<ValueT, RANK> &value()
		{ return dense; }
};

template<class IndexT, class ValueT, class LesserAccumulatorT>
class CollapseAccumulator : public BlankAccumulator
{
	LesserAccumulatorT sub;
	std::vector<int> keep_dims;		// Keep these dimensions; sum over the rest

	std::vector<IndexT> sub_index;	// Temporary
public:
	int rank() { return sub.rank(); }
	CollapseAccumulator(LesserAccumulatorT &&_sub, std::vector<int> &&_keep_dims)
		: sub(std::move(_sub)), keep_dims(std::move(_keep_dims))
	{
		if (keep_dims.size() != sub.rank()) {
			fprintf(stderr, "Rank mismatch: %d vs %d\n", keep_dims.size(), sub.rank());
			giss::exit(-1);
		}
		index.reserve(sub.rank());
	}

	void reserve(size_t size) { sub.reserve(size); }
	void add(std::array<IndexT, RANK> const index, ValueT const &val)
	{
		for (int k=0; k<sub.rank(); ++k) sub_index[k] = index[keep_dims[i]];
		sub.add(sub_index, val);
	}

	LesserAccumulatorT &value()
		{ return sub; }
}

template<int IndexT, class ValueT, int RANK>
class ScalarAccumulator : public BlankAccumulator {
	ValueT val = 0;

public:

	void add(std::array<IndexT, RANK> const index, ValueT const &val) { this->val += val; }
	void add(std::vector<IndexT> const &index, ValueT const &val) { this->val += val; }

	ValueT &value()
		{ return val; }
};

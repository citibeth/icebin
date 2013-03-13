#include <vector>
#include <giss/SparseMatrix.hpp>

namespace giss {

static std::vector<int> get_rowcol_beginnings(
	VectorSparseMatrix &a,
	int const rowcol,
	std::vector<int> &abegin)
{
	// Get beginning of each row in a (including sentinel at end)
	a.sort(SparseMatrix::SortOrder::ROW_MAJOR);
	int last_row = -1;
	for (auto ai(a.begin()); ; ++ai) {
		if (ai == a.end()) {
			abegin.push_back(ai.position());
			break;
		}
		if (ai.rowcol(rowcol) != last_row) {
			// arow.push_back(ai.rowcol(rowcol));
			abegin.push_back(ai.position());
			last_row = ai.rowcol(rowcol);
		}
	}
}

static double multiply_row_col(
	VectorSparseMatrix const &r,
	int const r0,	// Position within the sparsematrix vectors
	int const r1,
	VectorSparseMatrix const &c,
	int const c0,
	int const c1)
{
	double ret = 0;
	int ri = r0;
	int ci = c0;
	for (;;) {
		if (ri == r1) break;
		if (ci == c1) break;

		while (r.cols()[ri] < c.rows()[ci]) ++ri;
		while (c.rows()[ci] < r.cols()[ri]) ++ci;

		if (r.cols()[ri] == c.rows()[ci]) {
			// Duplicate entries ==> add together
			int rcol = r.cols()[ri];
			double rval = 0;
			do {
				rval += r.vals()[ri++];
			} while (r.cols()[ri] == rcol);

			// Duplicate entries ==> add together
			int crow = c.rows()[ci];
			double cval = 0;
			do { cval += c.vals()[ci++]; } while (c.rows()[ci] == crow);

			ret += rval * cval;
		}
	}
	return ret;
}


std::unique_ptr<VectorSparseMatrix> multiply(VectorSparseMatrix &a, VectorSparseMatrix &b)
{
	// Get beginning of each row in a (including sentinel at end)
	a.sort(SparseMatrix::SortOrder::ROW_MAJOR);
	std::vector<int> abegin(get_rowcol_beginnings(a, 0, abegin));

	// Get beginning of each col in b (including sentinel at end)
	b.sort(SparseMatrix::SortOrder::COLUMN_MAJOR);
	std::vector<int> bbegin(get_rowcol_beginnings(b, 1, bbegin));

	// Multiply each row by each column
	std::unique_ptr<VectorSparseMatrix> ret(new VectorSparseMatrix(
		SparseDescr(a.nrow, b.ncol)));;

	for (int ai = 0; ai < abegin.size()-1; ++ai) {
	for (int bi = 0; bi < bbegin.size()-1; ++bi) {
		// Multiply a row by a column
		double val = multiply_row_col(
			a, abegin[ai], abegin[ai+1],
			b, bbegin[bi], bbegin[bi+1]);
		ret->add(abegin[ai], bbegin[bi], val);
	}}

	return ret;
}

} // namespace giss

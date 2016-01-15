/*
 * GLINT2: A Coupling Library for Ice Models and GCMs
 * Copyright (c) 2013 by Robert Fischer
 * 
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include <vector>
#include <giss/SparseMatrix.hpp>
//#include <giss/eigen.hpp>

namespace giss {

// =======================================================================

#if 0
std::unique_ptr<VectorSparseMatrix> multiply_eigen_algorithm(VectorSparseMatrix &a, VectorSparseMatrix &b)
{
        // Get two matrices and convert to Eigen format.
        auto a_e(giss_to_Eigen(a));
		auto b_e(giss_to_Eigen(b));   // TODO: Consider making this one column-major to ease multiplication below.

        // Multiply the matices in Eigen format
        auto e_ret((*a_e) * (*b_e));

        // Convert back to GISS format sparse matrices
        return giss::Eigen_to_giss(e_ret);
}
#endif

// =======================================================================

std::vector<int> get_rowcol_beginnings(
	VectorSparseMatrix const &a,
	int const rowcol)
{
	std::vector<int> abegin;

	// Get beginning of each row in a (including sentinel at end)
//	a.sort(SparseMatrix::SortOrder::ROW_MAJOR);
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

	return abegin;
}

/** Computes [row] diag(scale) [[col]].
Or: row_j scale_j col_j   (out of overall matrix computation row_ij scale_j col_jk */

static double multiply_row_col(
	VectorSparseMatrix const &r,
	int const r0,	// Position within the sparsematrix vectors
	int const r1,
	blitz::Array<double,1> const * const scale,
	VectorSparseMatrix const &c,
	int const c0,
	int const c1)
{

	if (r0 == r1 || c0 == c1) return 0.0;

#if 0
int orow = r.rows()[r0];
int ocol = c.cols()[c0];

if (orow == 1 && ocol == 0) {
	printf("\n");
	printf("row: ");
	for (int i=r0; i<r1; ++i) printf("(%d,%d,%g) ", r.rows()[i], r.cols()[i], r.vals()[i]);
	printf("<<\n");

	printf("col: ");
	for (int i=c0; i<c1; ++i) printf("(%d,%d,%g) ", c.rows()[i], c.cols()[i], c.vals()[i]);
	printf("<<\n");
}
#endif


	double ret = 0;
	int ri = r0;
	int ci = c0;
	for (;;) {
		if (ri == r1) break;
		if (ci == c1) break;

		for (;;++ri) {
			if (ri >= r1) return ret;
			if (r.cols()[ri] >= c.rows()[ci]) break;
		}

		for (;;++ci) {
			if (ci >= c1) return ret;
			if (c.rows()[ci] >= r.cols()[ri]) break;
		}

//if (orow == 1 && ocol == 0) printf("ri=%d ci=%d\n", ri-r0, ci-c0);


		int rcol = r.cols()[ri];
		int crow = c.rows()[ci];
		if (rcol == crow) {			// Duplicate entries ==> add together
			double rval = 0;
			do {
				rval += r.vals()[ri++];
			} while (ri < r1 && r.cols()[ri] == rcol);

			// Duplicate entries ==> add together
			double cval = 0;
			do {
				cval += c.vals()[ci++];
			} while (ci < c1 && c.rows()[ci] == crow);

//if (orow == 1 && ocol == 0)  printf("      + %g * %g = %g\n", rval, cval, rval*cval);
			ret += rval * (scale ? (*scale)(crow) : 1.0) * cval;
		}
	}
	return ret;
}


/** Computes the sparse matrix:
      A diag(scaleB) B
*/
std::unique_ptr<VectorSparseMatrix> multiply(
//	blitz::Array<double,1> const *scale_i,			// scaleA_i
	VectorSparseMatrix &a,							// A_ij; sorts row major
	blitz::Array<double,1> const *scale_j			// scaleB_j
	VectorSparseMatrix &b,							// B_jk; sorts col major
	blitz::Array<double,1> const *scale_k			// scaleB_j
	VectorSparseMatrix &out)
{
	// Get beginning of each row in a (including sentinel at end)
	a.sort(SparseMatrix::SortOrder::ROW_MAJOR);
	std::vector<int> abegin(get_rowcol_beginnings(a, 0));

#if 0
printf("----- A\n");
{ int i=0; for (auto ii=a.begin(); ii != a.end(); ++ii) printf("%d: (%d, %d) --> %g\n", i++, ii.row(), ii.col(), ii.val()); }
for (auto ii=abegin.begin(); ii != abegin.end(); ++ii) printf("%d, ", *ii);
printf("\n");
#endif

	// Get beginning of each col in b (including sentinel at end)
	b.sort(SparseMatrix::SortOrder::COLUMN_MAJOR);
	std::vector<int> bbegin(get_rowcol_beginnings(b, 1));

#if 0
printf("----- B\n");
{ int i=0; for (auto ii=b.begin(); ii != b.end(); ++ii) printf("%d: (%d, %d) --> %g\n", i++, ii.row(), ii.col(), ii.val()); }
for (auto ii=bbegin.begin(); ii != bbegin.end(); ++ii) printf("%d, ", *ii);
printf("\n");
#endif

	// Multiply each row by each column
	std::unique_ptr<VectorSparseMatrix> ret(new VectorSparseMatrix(
		SparseDescr(a.nrow, b.ncol)));;

	for (int ai = 0; ai < abegin.size()-1; ++ai) {
		// Make sure this row is non-empty
		if (abegin[ai] == abegin[ai+1]) continue;
		int aRow = a.rows()[ai];	// Index of this row
		double iScale = scale_i ? (*scale_i)(aRow) : 1.0;

		for (int bi = 0; bi < bbegin.size()-1; ++bi) {
			// Make sure this column is non-empty
			if (bbegin[bi] == bbegin[bi+1]) continue;
			int bCol = b.rows()[bi];	// Index of this column
			double kScale = scale_k ? scale_k[bRow] : 1.0;

			// Multiply a row by a column
			double val =
				iScale *
				multiply_row_col(
					a, abegin[ai], abegin[ai+1],
					scale_j,
					b, bbegin[bi], bbegin[bi+1])
				* kScale;

			if (val == 0.0) continue;

		int row = a.rows()[abegin[ai]];
		int col = b.cols()[bbegin[bi]];
		ret->add(row, col, val);
	}}

	return ret;
}

} // namespace giss
